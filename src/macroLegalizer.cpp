#include "macroLegalizer.h"
#include "placerBase.h"
#include "nesterovBase.h"
#include <algorithm>
#include <random>
#include "point.h"
#include <cmath>
#include <ctime>
#include "log.h"

namespace replace
{
    ///////////////////////////////////////
    // MacroLegalizerVars

    MacroLegalizerVars::MacroLegalizerVars()
        : maxPostLegalizeIter(1000)
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
            saLegalize(macros, die);
            postLegalize(macros, die);
            // LOG_INFO("[MacroLegalization] Die {} ", die->name());
        }
    }

    void MacroLegalizer::saLegalize(const std::vector<Instance *> macros, Die *die)
    {
        if (macros.size() == 0)
            return;
        for (int iter = 0; iter < vars_.sa_max_iter; iter++)
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
            srand(iter * 114 + 514);
            // double cooling_rate_ = cooling_rate;
            //  选择一个宏单元
            int index = rand() % macros.size();

            for (int i = 0; i < vars_.sa_max_iter0; i++)
            {
                double temp = temp_0 + (temp_max - temp_0) * (double)(i / vars_.sa_max_iter0);
                Instance *cell = macros[index];
                // 计算宏单元的成本函数
                double old_cost = calc_cost(macros, die);
                // 随机移动
                std::pair<int, int> move = getRandomMove(cell, iter, die, max_sa_r_x, max_sa_r_y);
                // 尝试将宏单元移动到新的位置
                cell->setLocation(cell->lx() + move.first, cell->ly() + move.second);
                // 计算宏单元的成本函数
                double new_cost = calc_cost(macros, die);
                // 判断是否接受新的解
                double delta = (new_cost - old_cost) / old_cost;
                double tau = (double)rand() / RAND_MAX;
                double p = exp(-1.0 * delta / temp);
                int isAccept=0;
                LOG_INFO("[MacroLegalization] old_cost: {} new_cost: {} move: {}-{} delta{}", old_cost, new_cost, move.first, move.second, delta);
                if (p > tau)
                {
                    // 接受新的解
                    // 判断是否满足要求
                    isAccept=1;
                    int om = getMacrosOverlap(macros, die);
                    if (om == 0)
                    {
                        break;
                    }
                }
                else
                {
                    // 拒绝新的解
                    isAccept=0;
                    cell->setLocation(cell->lx() - move.first, cell->ly() - move.second);
                }

                LOG_INFO("[MacroLegalization] Iter {} new_cost: {} accept: {}", i, new_cost, isAccept);
            }
        }
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
        int rndx = rand();
        int rndy = rand();
        int rh = die->rowHeight();
        double rcof = 0.05 * pow(1.5, (double)iter);
        double rx = (max_sa_r_x - rh) / (double)(iter + 1); // 1 = rowheight
        double ry = (max_sa_r_y - rh) / (double)(iter + 1);
        int moveX = int((rndx / RAND_MAX - 0.5) * rx);
        int moveY = int((rndy / RAND_MAX - 0.5) * ry);
        LOG_INFO("[MacroLegalization] moveX: {} moveY: {} ", moveX, moveY);
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

    // 计算两个矩形的重叠面积
    int MacroLegalizer::overlapArea(Instance *cell1, Instance *cell2)
    {
        int x_overlap = std::max(0, std::min(cell1->ux(), cell2->ux()) - std::max(cell2->lx(), cell2->lx()));
        int y_overlap = std::max(0, std::min(cell1->uy(), cell2->uy()) - std::max(cell2->ly(), cell2->ly()));
        return x_overlap * y_overlap;
    }

    int MacroLegalizer::get_hpwl(const std::vector<Instance *> &macros, Die *die)
    {
        return pb_->hpwl();
        /*
        int64_t getHpwl();

        nt64_t NesterovBase::getHpwl()
        {
        int64_t hpwl = 0;
        for(auto& gNet : gNets_)
        {
            gNet->updateBox();
            hpwl += gNet->hpwl();
        }
        return hpwl;
        }

        int64_t
        GNet::hpwl() {
        return static_cast<int64_t>((ux_ - lx_) + (uy_ - ly_));
        }

        void Net::updateBox()
        {
            lx_ = INT_MAX;
            ly_ = INT_MAX;
            ux_ = INT_MIN;
            uy_ = INT_MIN;

            for (const Pin *p : pins_)
            {
            lx_ = std::min(p->cx(), lx_);
            ux_ = std::max(p->cx(), ux_);
            ly_ = std::min(p->cy(), ly_);
            uy_ = std::max(p->cy(), uy_);
            }
        }
        */
    }

    int MacroLegalizer::getCellMacroOverlap(const std::vector<Instance *> &macros, Die *die)
    {
        // TODO: use quadtree

        int totalOverlap = 0;
        // 外层macro里层stdcell
        for (int i = 0; i < macros.size(); i++)
        {
            for (int j = 0; j < die->insts().size(); j++)
            {
                if (!die->insts()[j]->isMacro())
                {
                    totalOverlap += overlapArea(macros[i], die->insts()[j]);
                }
            }
        }
        return totalOverlap;
    }

    int MacroLegalizer::getMacrosOverlap(const std::vector<Instance *> &macros, Die *die)
    {
        int totalOverlap = 0;
        for (int i = 0; i < macros.size(); i++)
        {
            for (int j = i + 1; j < macros.size(); j++)
            {
                totalOverlap += overlapArea(macros[i], macros[j]);
            }
        }
        return totalOverlap;
    }

    // 计算宏单元的成本函数
    double MacroLegalizer::calc_cost(const std::vector<Instance *> &macros, Die *die)
    {
        // HPWL
        double hpwl = get_hpwl(macros, die);
        // 宏单元与标准单元重叠
        double den = getCellMacroOverlap(macros, die);
        // 宏单元与宏单元重叠
        double ov = getMacrosOverlap(macros, die);
        return vars_.sa_hpwl_wgt * hpwl + vars_.sa_den_wgt * den + vars_.sa_ovlp_wgt * ov;
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

    void MacroLegalizer::postLegalize(const std::vector<Instance *> insts, Die *die)
    {
        srand(114);
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

                        if (((float)rand() / RAND_MAX) < 0.5)
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
                break;
            }
        }
    }
}