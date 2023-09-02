#include "macroLegalizer.h"
#include "placerBase.h"
#include "nesterovBase.h"
#include <algorithm>
#include <random>
#include "point.h"
#include <cmath>
#include <ctime>
#include "log.h"
#include <stdlib.h>

namespace replace
{
    ///////////////////////////////////////
    // MacroLegalizerVars

    MacroLegalizerVars::MacroLegalizerVars()
        : maxPostLegalizeIter(2000)
    {
    }

    /////////////////////////////////////////
    // MacroLegalizer

    MacroLegalizer::MacroLegalizer(MacroLegalizerVars lgVars, std::shared_ptr<PlacerBase> pb)
        : pb_(pb), lgVars_(lgVars)
    {
    }

    void MacroLegalizer::doLegalization()
    {
        for (Die *die : pb_->dies())
        {
            macros_.clear();
            die_ = die;
            for (Instance *inst : die->insts())
            {
                // regardless of fixed instance
                if (inst->isMacro())
                    macros_.push_back(inst);
            }
            bool isLegal = saLegalize2();
            if(!isLegal)
              postLegalize(macros_, die_);
        }
    }

    ////////////////////////////////
    // saLegalize

    static std::mt19937 rng;

    double MacroLegalizer::getOverlapArea(const Instance *inst1, const Instance *inst2)
    {
        int lx = std::max(inst1->lx(), inst2->lx());
        int ux = std::min(inst1->ux(), inst2->ux());
        int ly = std::max(inst1->ly(), inst2->ly());
        int uy = std::min(inst1->uy(), inst2->uy());
        if (lx < ux && ly < uy)
            return static_cast<double>(ux - lx) * (uy - ly);
        else
            return 0;
    }

    double MacroLegalizer::getMacroOverlap(const Instance *macro, const std::vector<Instance *> &macros)
    {
        double totalOvarlap = 0;
        for (Instance *mac : macros)
        {
            if (mac == macro)
                continue;
            totalOvarlap += getOverlapArea(mac, macro);
        }
        return totalOvarlap;
    }

    double MacroLegalizer::getAllMacroOverlap(const std::vector<Instance *> &macros)
    {
        double totalOverlap = 0;
        for (int i = 0; i < macros.size(); i++)
        {
            for (int j = i + 1; j < macros.size(); j++)
                totalOverlap += getOverlapArea(macros[i], macros[j]);
        }
        return totalOverlap;
    }

    double MacroLegalizer::getMacroHpwl(const Instance *macro)
    {
        int64_t macroHpwl = 0;
        for (Pin *pin : macro->pins())
        {
            Net *net = pin->net();
            macroHpwl += net->hpwl();
        }
        return static_cast<double>(macroHpwl);
    }

    void MacroLegalizer::updateMacroNetBox(Instance *macro)
    {
        for (Pin *pin : macro->pins())
        {
            Net *net = pin->net();
            net->updateBox();
        }
    }

    bool MacroLegalizer::saLegalize2()
    {
        if(macros_.size() == 0)
            return true;
        sa_init_top();
        sa_mac_leg_top();
        return ovlp_free_flg;
    }

    void MacroLegalizer::sa_init_top()
    {
        for (Instance* macro : macros_)
        {
            int lx = macro->lx();
            int ly = macro->ly();
            // 对齐到row
            int rowHeight = die_->isSetRow() ? die_->rowHeight() : 1;
            ly = die_->coreLy() + static_cast<int>((ly - die_->coreLy()) / (double)rowHeight + 0.5) * rowHeight;
            // 边界处理
            lx = std::max(std::min(lx, die_->coreUx() - macro->dx()), die_->coreLx());
            ly = std::max(std::min(ly, die_->coreUy() - macro->dy()), die_->coreLy());
            macro->setLocation(lx, ly);
            updateMacroNetBox(macro);
        }

        tot_mac_hpwl = pb_->hpwl();
        tot_mac_den = 0;
        init_bins();
        for(const Instance* mac : macros_)
        {
            tot_mac_den += get_mac_den(mac);
        }
        tot_mac_ovlp = getAllMacroOverlap(macros_);
        ovlp_free_flg = (tot_mac_ovlp == 0);

        sa_param_init_top();
    }

    void MacroLegalizer::sa_param_init_top()
    {
        sa_hpwl_wgt = 1.0;
        // sa_den_wgt = 0.0; // we never use den
        sa_den_wgt = sa_hpwl_wgt * tot_mac_hpwl / tot_mac_den;
        sa_ovlp_wgt = (tot_mac_hpwl * sa_hpwl_wgt + tot_mac_den * sa_den_wgt) / tot_mac_ovlp;
        sa_hpwl_cof = 1.0;
        sa_den_cof = 1.0;
        sa_ovlp_cof = 1.5;
        sa_max_iter0 = 20;
    }

    void MacroLegalizer::sa_mac_leg_top()
    {
        for(int iter = 0; iter < sa_max_iter0 && !ovlp_free_flg; iter++)
        {
            LOG_DEBUG("iter {}, hpwl: {}, den: {}, ovlp: {}", iter+1, tot_mac_hpwl, tot_mac_den, tot_mac_ovlp);
            sa_param_init(iter);
            sa_mac_leg(iter);
            sa_param_update();
        }
    }

    void MacroLegalizer::sa_param_init(int iter)
    {
        double sa_coef = 0.03;
        sa_iter_cof = 1;
        sa_max_iter = 1000;
        sa_max_iter2 = sa_iter_cof * macros_.size();

        sa_init_neg_rate = sa_coef * std::pow(1.5, static_cast<double>(iter));
        sa_last_neg_rate = 0.0001 * std::pow(1.5, static_cast<double>(iter));

        sa_init_t = sa_init_neg_rate / std::log(2.0);
        sa_t = sa_init_t;
        sa_t_cof = std::pow(sa_last_neg_rate / sa_init_neg_rate,
                            1.0 / static_cast<double>(sa_max_iter));

        sa_n_x = sa_n_y = std::sqrt(static_cast<double>(macros_.size()));
        sa_ncof_x = sa_ncof_y = 0.05 * std::pow(1.5, static_cast<double>(iter));
        min_sa_r_x = 1.0;
        min_sa_r_y = die_->isSetRow() ? die_->rowHeight() : 1;
        max_sa_r_x = static_cast<double>(die_->coreDx()) / sa_n_x * sa_ncof_x;
        max_sa_r_y = static_cast<double>(die_->coreDy()) / sa_n_y * sa_ncof_y;
        max_sa_r_x = std::max(max_sa_r_x, min_sa_r_x);
        max_sa_r_y = std::max(max_sa_r_y, min_sa_r_y);
        sa_r_x = max_sa_r_x;
        sa_r_y = max_sa_r_y;
        sa_r_stp_x = (max_sa_r_x - min_sa_r_x) / static_cast<double>(sa_max_iter);
        sa_r_stp_y = (max_sa_r_y - min_sa_r_y) / static_cast<double>(sa_max_iter);
    }

    void MacroLegalizer::sa_mac_leg(int iter)
    {
        for(int i = 0; i < sa_max_iter && !ovlp_free_flg; i++)
        {
            for(int j = 0; j < sa_max_iter2 && !ovlp_free_flg; j++)
                sa_mac_leg_sub();
            sa_param_update_sub();
        }
    }

    void MacroLegalizer::sa_param_update()
    {
        sa_hpwl_wgt *= sa_hpwl_cof;
        sa_den_wgt *= sa_den_cof;
        sa_ovlp_wgt *= sa_ovlp_cof;
    }

    void MacroLegalizer::sa_mac_leg_sub()
    {
        //  选择一个宏单元
        int mac_idx = sa_get_mac_idx();
        Instance* mac = macros_[mac_idx];
        // 计算宏单元的成本函数
        double old_mac_hpwl, old_mac_den, old_mac_ovlp;
        double mac_c0 = sa_get_mac_cost(mac, old_mac_hpwl, old_mac_den, old_mac_ovlp);
        // 随机移动
        std::pair<int, int> mov = sa_get_mac_mov(mac);
        sa_do_mac_mov(mac, mov);
        // 计算宏单元的成本函数
        double new_mac_hpwl, new_mac_den, new_mac_ovlp;
        double mac_c1 = sa_get_mac_cost(mac, new_mac_hpwl, new_mac_den, new_mac_ovlp);
        // 判断是否接受新的解
        bool accept = sa_mov_accept(mac_c0, mac_c1);

        if(accept)
        {
            tot_mac_den += (new_mac_den - old_mac_den);
            tot_mac_ovlp = getAllMacroOverlap(macros_);
            tot_mac_hpwl += (new_mac_hpwl - old_mac_hpwl);
            ovlp_free_flg = (tot_mac_ovlp == 0);
        }
        else
        {
            mov.first = -mov.first;
            mov.second = -mov.second;
            sa_do_mac_mov(mac, mov);
        }
    }

    void MacroLegalizer::sa_param_update_sub()
    {
        sa_t *= sa_t_cof;  // 0.999999 ;
        sa_r_x -= sa_r_stp_x;
        sa_r_y -= sa_r_stp_y;
    }

    int MacroLegalizer::sa_get_mac_idx()
    {
        std::uniform_int_distribution<int> uni(0, macros_.size() - 1);
        return uni(rng);
    }

    std::pair<int, int> MacroLegalizer::sa_get_mac_mov(const Instance* mac)
    {
        std::uniform_real_distribution<float> unf(-0.5f, 0.5f);
        float rndx = unf(rng);
        float rndy = unf(rng);
        int rowHeight = die_->isSetRow() ? die_->rowHeight() : 1;
        int mov_x = static_cast<int>(rndx * sa_r_x);
        int mov_y = static_cast<int>(rndy * sa_r_y / rowHeight) * rowHeight;

        int target_lx = mac->lx() + mov_x;
        int target_ly = mac->ly() + mov_y;
        target_lx = std::max(die_->coreLx(), target_lx);
        target_lx = std::min(die_->coreUx() - mac->dx(), target_lx);
        target_ly = std::max(die_->coreLx(), target_ly);
        target_ly = std::min(die_->coreUy() - mac->dy(), target_ly);

        mov_x = target_lx - mac->lx();
        mov_y = target_ly - mac->ly();

        return std::make_pair(mov_x, mov_y);
    }

    void MacroLegalizer::sa_do_mac_mov(Instance* mac, const std::pair<int, int>& mov)
    {
        mac->setLocation(mac->lx() + mov.first, mac->ly() + mov.second);
        updateMacroNetBox(mac);
    }

    double MacroLegalizer::sa_get_mac_cost(const Instance* mac, double& hpwl, double& den, double& ovlp)
    {
        hpwl = getMacroHpwl(mac);
        // double den = 0;
        den = get_mac_den(mac);
        ovlp = getMacroOverlap(mac, macros_);
        double tot_cost = sa_hpwl_wgt * hpwl 
                         + sa_den_wgt * den
                         + sa_ovlp_wgt * ovlp;
        return tot_cost;
    }

    bool MacroLegalizer::sa_mov_accept(double old_cost, double new_cost)
    {
        std::uniform_real_distribution<double> und(0.0, 1.0);
        double dc = (new_cost - old_cost) / old_cost;
        double drnd = und(rng);
        double exp_val = std::exp(-1.0 * dc / sa_t);
        if(drnd < exp_val)
            return true;
        else
            return false;
    }

 void MacroLegalizer::init_bins()
    {
        // 初始化bin的大小为最小的macro的1/16
        bin_size_x = bin_size_y = std::numeric_limits<double>::max();
        for(const Instance* mac : macros_)
        {
            bin_size_x = std::min(static_cast<double>(mac->dx()), bin_size_x);
            bin_size_y = std::min(static_cast<double>(mac->dy()), bin_size_y);
        }
        bin_size_x /= 4.0;
        bin_size_y /= 4.0;
        bin_count_x = static_cast<int>(std::round(die_->coreDx() / bin_size_x));
        bin_count_y = static_cast<int>(std::round(die_->coreDy() / bin_size_y));
        bin_size_x = die_->coreDx() / static_cast<double>(bin_count_x);
        bin_size_y = die_->coreDy() / static_cast<double>(bin_count_y);

        bin_density.resize(bin_count_x * bin_count_y, 0.0);
        for(const Instance* cell : die_->insts())
        {
            if(cell->isMacro())
                continue;
            int idx_x_low = static_cast<int>((cell->lx() - die_->coreLx()) / bin_size_x);
            int idx_x_high = static_cast<int>((cell->ux() - die_->coreLx()) / bin_size_x);
            int idx_y_low = static_cast<int>((cell->ly() - die_->coreLy()) / bin_size_y);
            int idx_y_high = static_cast<int>((cell->uy() - die_->coreLy()) / bin_size_y);
            idx_x_low = std::min(std::max(idx_x_low, 0), bin_count_x-1);
            idx_x_high = std::min(std::max(idx_x_high, 0), bin_count_x-1);
            idx_y_low = std::min(std::max(idx_y_low, 0), bin_count_y-1);
            idx_y_high = std::min(std::max(idx_y_high, 0), bin_count_y-1);

            for(int x = idx_x_low; x <= idx_x_high; x++)
            {
                double bin_lx = x * bin_size_x;
                double bin_ux = (x + 1) * bin_size_x;
                double ovlp_lx = std::max(static_cast<double>(cell->lx()), bin_lx);
                double ovlp_ux = std::min(static_cast<double>(cell->ux()), bin_ux);
                for(int y = idx_y_low; y <= idx_y_high; y++)
                {
                    double bin_ly = y * bin_size_y;
                    double bin_uy = (y + 1) * bin_size_y;
                    double ovlp_ly = std::max(static_cast<double>(cell->ly()), bin_ly);
                    double ovlp_uy = std::min(static_cast<double>(cell->uy()), bin_uy);

                    if(ovlp_lx < ovlp_ux && ovlp_ly < ovlp_uy)
                    {
                        double ovlp_area = (ovlp_ux - ovlp_lx) * (ovlp_uy - ovlp_ly);
                        bin_density[x * bin_count_y + y] += ovlp_area;
                    }
                }
            }
        }

        double bin_area = bin_size_x * bin_size_y;
        for(int x = 0; x < bin_count_x; x++)
        {
            for(int y = 0; y < bin_count_y; y++)
                bin_density[x * bin_count_y + y] /= bin_area;
        }
    }

    double MacroLegalizer::get_mac_den(const Instance* mac)
    {
        int idx_x_low = static_cast<int>((mac->lx() - die_->coreLx()) / bin_size_x);
        int idx_x_high = static_cast<int>((mac->ux() - die_->coreLx()) / bin_size_x);
        int idx_y_low = static_cast<int>((mac->ly() - die_->coreLy()) / bin_size_y);
        int idx_y_high = static_cast<int>((mac->uy() - die_->coreLy()) / bin_size_y);
        idx_x_low = std::min(std::max(idx_x_low, 0), bin_count_x-1);
        idx_x_high = std::min(std::max(idx_x_high, 0), bin_count_x-1);
        idx_y_low = std::min(std::max(idx_y_low, 0), bin_count_y-1);
        idx_y_high = std::min(std::max(idx_y_high, 0), bin_count_y-1);

        double den = 0.0;
        for(int x = idx_x_low; x <= idx_x_high; x++)
        {
            double bin_lx = x * bin_size_x;
            double bin_ux = (x + 1) * bin_size_x;
            double ovlp_lx = std::max(static_cast<double>(mac->lx()), bin_lx);
            double ovlp_ux = std::min(static_cast<double>(mac->ux()), bin_ux);
            for(int y = idx_y_low; y <= idx_y_high; y++)
            {
                double bin_ly = y * bin_size_y;
                double bin_uy = (y + 1) * bin_size_y;
                double ovlp_ly = std::max(static_cast<double>(mac->ly()), bin_ly);
                double ovlp_uy = std::min(static_cast<double>(mac->uy()), bin_uy);

                if(ovlp_lx < ovlp_ux && ovlp_ly < ovlp_uy)
                {
                    double ovlp_area = (ovlp_ux - ovlp_lx) * (ovlp_uy - ovlp_ly);
                    den += ovlp_area * bin_density[x * bin_count_y + y];
                }
            }
        }
        return den;
    }

    //////////////////////////////////
    // Post Legalization

    static bool repel(int l1, int u1, int l2, int u2, int ll, int uu, int *pd1, int *pd2)
    {
        if (l1 <= l2 && l2 < u1)
        {
            // |--------|===========|-------|
            // ^        ^           ^       ^
            // l1       l2          u1      u2
            if (u1 < u2)
            {
                int d = u1 - l2;
                // try: inst1 <--- d/2
                int d1 = std::min(d / 2, l1 - ll);
                // try: inst2 --> d - d1
                int d2 = std::min(d - d1, uu - u2);
                // if inst2 can't --> d - d1, then inst1 <-- d - d2
                d1 = std::min(d - d2, l1 - ll);

                *pd1 = -d1;
                *pd2 = d2;
                return d1 + d2 == d;
            }

            // |--------|===========|-------|
            // ^        ^           ^       ^
            // l1       l2          u2      u1
            else
            {
                int d = std::min(u1 - l2, u2 - l1);
                // try: inst1 <-- d / 2
                int d1 = std::min(d / 2, l1 - ll);
                // try: inst2 --> d - d1
                int d2 = std::min(d - d1, uu - u2);
                // if inst2 can't --> d - d1, then inst1 <-- d - d2
                d1 = std::min(d - d2, l1 - ll);

                *pd1 = -d1;
                *pd2 = d2;
                return d1 + d2 == d;
            }
        }
        if (l2 <= l1 && l1 < u2)
        {
            // |--------|===========|-------|
            // ^        ^           ^       ^
            // l2       l1          u2      u1
            if (u1 < u2)
            {
                int d = u2 - l1;
                int d2 = std::min(d / 2, l2 - ll);
                int d1 = std::min(d - d2, uu - u1);
                d2 = std::min(d - d1, l2 - ll);

                *pd1 = d1;
                *pd2 = -d2;
                return d1 + d2 == d;
            }

            // |--------|===========|-------|
            // ^        ^           ^       ^
            // l2       l1          u1      u2
            else
            {
                int d = std::min(u1 - l2, u2 - l1);
                int d2 = std::min(d / 2, l2 - ll);
                int d1 = std::min(d - d2, uu - u1);
                d2 = std::min(d - d1, l2 - ll);

                *pd1 = -d1;
                *pd2 = d2;
                return d1 + d2 == d;
            }
        }

        *pd1 = 0;
        *pd2 = 0;
        return true;
    }

    bool MacroLegalizer::postLegalize(const std::vector<Instance *> &insts, const Die* die)
    {
        rng.seed(514);
        std::uniform_real_distribution<float> unf(0.f, 1.f);
        for (int iter = 0; iter < lgVars_.maxPostLegalizeIter; iter++)
        {
            bool isLegal = true;
            for (int i = 0; i < insts.size(); i++)
            {
                // check boundary
                Instance *inst1 = insts[i];
                if (inst1->lx() < die->coreLx())
                {
                    inst1->setLocation(die->coreLx(), inst1->ly());
                    isLegal = false;
                }
                if (inst1->ly() < die->coreLy())
                {
                    inst1->setLocation(inst1->lx(), die->coreLy());
                    isLegal = false;
                }
                if (inst1->ux() > die->coreUx())
                {
                    inst1->setLocation(die->coreUx() - inst1->dx(), inst1->ly());
                    isLegal = false;
                }
                if (inst1->uy() > die->coreUy())
                {
                    inst1->setLocation(inst1->lx(), die->coreUy() - inst1->dy());
                    isLegal = false;
                }

                for (int j = 0; j < insts.size(); j++)
                {
                    Instance *inst2 = insts[j];
                    if (i == j)
                        continue;

                    // check overlap
                    int lx = std::max(inst1->lx(), inst2->lx());
                    int ly = std::max(inst1->ly(), inst2->ly());
                    int ux = std::min(inst1->ux(), inst2->ux());
                    int uy = std::min(inst1->uy(), inst2->uy());

                    if (lx < ux && ly < uy)
                    {
                        isLegal = false;

                        int dx1, dx2, dy1, dy2;
                        // x overlap
                        bool xOk = repel(inst1->lx(), inst1->ux(), inst2->lx(), inst2->ux(),
                                         die->coreLx(), die->coreUx(), &dx1, &dx2);
                        // y overlap
                        auto yOk = repel(inst1->ly(), inst1->uy(), inst2->ly(), inst2->uy(),
                                         die->coreLy(), die->coreUy(), &dy1, &dy2);
                        double thres = (double)(std::abs(dx1) + std::abs(dx2)) / (std::abs(dx1) + std::abs(dx2) + std::abs(dy1) + std::abs(dy2));
                        if (unf(rng) < thres)
                        {
                            inst1->setLocation(inst1->lx() + dx1, inst1->ly());
                            inst2->setLocation(inst2->lx() + dx2, inst2->ly());
                        }
                        else
                        {
                            inst1->setLocation(inst1->lx(), inst1->ly() + dy1);
                            inst2->setLocation(inst2->lx(), inst2->ly() + dy2);
                        }
                    }
                }
            }

            if (isLegal)
            {
                return true;
            }
        }
        return false;
    }

    bool MacroLegalizer::checkLegal(const std::vector<Instance *> &macros, const Die* die)
    {
        bool isLegal = true;
        for (int i = 0; i < macros.size() && isLegal; i++)
        {
            // check boundary
            Instance *m1 = macros[i];
            if (m1->lx() < die->coreLx())
                isLegal = false;
            if (m1->ly() < die->coreLy())
                isLegal = false;
            if (m1->ux() > die->coreUx())
                isLegal = false;
            if (m1->uy() > die->coreUy())
                isLegal = false;
            for (int j = 0; j < i && isLegal; j++)
            {
                Instance *m2 = macros[j];
                if(getOverlapArea(m1, m2) > 0)
                    isLegal = false;
            }
        }
        return isLegal;
    }
}