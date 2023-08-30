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
        : maxPostLegalizeIter(100000)
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
        std::vector<Instance *> macros;
        for (Die *die : pb_->dies())
        {
            macros.clear();
            for (Instance *inst : die->insts())
            {
                // regardless of fixed instance
                if (inst->isMacro())
                {
                    macros.push_back(inst);
                }
            }
            if(!checkLegal(macros, die))
              saLegalize(macros, die);
            if(!checkLegal(macros, die))
              postLegalize(macros, die);
        }
    }

    ////////////////////////////////
    // saLegalize

    static std::mt19937 rng;

    bool MacroLegalizer::saLegalize(const std::vector<Instance *> &macros, Die *die)
    {
        if (macros.size() == 0)
            return true;
        // 初始化参数
        double tot_mac_hpwl = 0;
        for(Instance* macro : macros)
            tot_mac_hpwl += getMacroHpwl(macro);
        double tot_mac_ovlp = getAllMacroOverlap(macros);
        lgVars_.sa_hpwl_wgt = 1.0;
        lgVars_.sa_den_wgt = 0.0; // ignore
        lgVars_.sa_ovlp_wgt = tot_mac_hpwl / tot_mac_ovlp;
        lgVars_.sa_hpwl_cof = 1.0;
        lgVars_.sa_den_cof = 1.0;
        lgVars_.sa_ovlp_cof = 1.5;
        lgVars_.sa_max_iter0 = 1000;

        double isLegal = false;
        for (int iter = 0; iter < lgVars_.sa_max_iter && !isLegal; iter++)
        {
            // 初始化参数
            double cooling_rate;
            double temp_0 = 0.03 * pow(1.5, (double)iter) / log(2.0);
            double temp_max = 0.0001 * pow(1.5, (double)iter) / log(2.0);
            double rcof = 0.05 * pow(1.5, (double)iter);

            int64_t tot_place_region = die->coreDx() * (int64_t)die->coreDy();
            int max_sa_r_x = static_cast<int>(tot_place_region / sqrt(macros.size()) * rcof);
            int max_sa_r_y = static_cast<int>(tot_place_region / sqrt(macros.size()) * rcof);

            // 为使结果可复现，应手动设置seed
            rng.seed(iter * 114 + 514);
            std::uniform_int_distribution<int> uni(0, macros.size() - 1);
            std::uniform_real_distribution<double> unf(0.0, 1.0);
            for (int i = 0; i < lgVars_.sa_max_iter0 && !isLegal; i++)
            {
                double temp = temp_0 + (temp_max - temp_0) * (double)(i / lgVars_.sa_max_iter0);
                //  选择一个宏单元
                int index = uni(rng);
                Instance *macro = macros[index];
                // 计算宏单元的成本函数
                double old_cost = getMacroCost(macro, macros);
                // 随机移动
                std::pair<int, int> move = getRandomMove(macro, iter, die, max_sa_r_x, max_sa_r_y);
                // 尝试将宏单元移动到新的位置
                macro->setLocation(macro->lx() + move.first, macro->ly() + move.second);
                updateMacroNetBox(macro);
                // 计算宏单元的成本函数
                double new_cost = getMacroCost(macro, macros);
                // 判断是否接受新的解
                double delta = (new_cost - old_cost) / old_cost;
                double tau = unf(rng);
                double p = std::exp(-1.0 * delta / temp);
                int isAccept = 0;
                LOG_DEBUG("[MacroLegalization] old_cost: {} new_cost: {} move: {}-{} delta{}", old_cost, new_cost, move.first, move.second, delta);
                if (p > tau)
                {
                    // 接受新的解
                    // 判断是否满足要求
                    isAccept = 1;
                    int om = getAllMacroOverlap(macros);
                    isLegal = (om == 0);
                }
                else
                {
                    // 拒绝新的解
                    isAccept = 0;
                    macro->setLocation(macro->lx() - move.first, macro->ly() - move.second);
                    updateMacroNetBox(macro);
                }
                LOG_DEBUG("[MacroLegalization] Iter {} new_cost: {} accept: {}", i, new_cost, isAccept);
            }

            // update param
            lgVars_.sa_hpwl_wgt *= lgVars_.sa_hpwl_cof;
            lgVars_.sa_den_wgt *= lgVars_.sa_den_cof;
            lgVars_.sa_ovlp_wgt *= lgVars_.sa_ovlp_cof;
        }

        return isLegal;
    }

    std::pair<int, int> MacroLegalizer::getRandomMove(Instance *cell, int iter, Die *die, int max_sa_r_x, int max_sa_r_y)
    {
        int clx = cell->lx();
        int cly = cell->ly();
        int cux = cell->ux();
        int cuy = cell->uy();
        // 获得die的大小
        int dlx = die->dieLx();
        int dly = die->dieLy();
        int dux = die->dieUx();
        int duy = die->dieUy();
        int diewidth = dux - dlx;
        int dieheight = duy - dly;
        std::uniform_real_distribution<float> unf(-0.5f, 0.5f);
        float rndx = unf(rng);
        float rndy = unf(rng);
        int rh = die->rowHeight();
        double rcof = 0.05 * pow(1.5, (double)iter);
        double rx = (max_sa_r_x - rh) / (double)(iter + 1); // 1 = rowheight
        double ry = (max_sa_r_y - rh) / (double)(iter + 1);
        int moveX = static_cast<int>(rndx * rx);
        int moveY = static_cast<int>(rndy * ry);
        LOG_DEBUG("[MacroLegalization] moveX: {} moveY: {} ", moveX, moveY);
        // 计算移动后的 Macro 的四个角坐标
        int newcux = cux + moveX;
        int newcuy = cuy + moveY;
        int newclx = clx + moveX;
        int newcly = cly + moveY;
        // 检查移动后的坐标是否超过 Die 的边界范围，若超过则修正坐标
        if (newclx < dlx)
        {
            int diff = dlx - newclx;
            moveX += diff;
        }
        if (newcuy > duy)
        {
            int diff = newcuy - duy;
            moveY -= diff;
        }
        if (newcux > dux)
        {
            int diff = newcux - dux;
            moveX -= diff;
        }
        if (newcly < dly)
        {
            int diff = newcly - dly;
            moveY += diff;
        }
        return std::make_pair(moveX, moveY);
    }

    double MacroLegalizer::getMacroCost(Instance *macro, const std::vector<Instance *> &macros)
    {
        double hpwlCost = getMacroHpwl(macro);
        double denCost = 0;
        double overlapCost = getMacroOverlap(macro, macros);
        double totalCost = lgVars_.sa_hpwl_wgt * hpwlCost + lgVars_.sa_den_wgt * denCost + lgVars_.sa_ovlp_wgt * overlapCost;
        return totalCost;
    }

    double MacroLegalizer::getOverlapArea(Instance *inst1, Instance *inst2)
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

    double MacroLegalizer::getMacroOverlap(Instance *macro, const std::vector<Instance *> &macros)
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

    double MacroLegalizer::getMacroHpwl(Instance *macro)
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

    bool MacroLegalizer::postLegalize(const std::vector<Instance *> insts, Die *die)
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

    bool MacroLegalizer::checkLegal(const std::vector<Instance *> macros, Die *die)
    {
        bool isLegal = true;
        for (int i = 0; i < macros.size(); i++)
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