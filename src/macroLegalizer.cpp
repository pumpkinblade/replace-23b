#include "macroLegalizer.h"
#include "placerBase.h"
#include "nesterovBase.h"
#include <algorithm>
#include <random>
#include "point.h"
#include <cmath>
#include <ctime>

namespace replace
{
    int max_sa_r_x, max_sa_r_y; 

    bool repel(int l1, int u1, int l2, int u2, int ll, int uu, int* pd1, int* pd2)
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
    
    MacroLegalizer::MacroLegalizer()
    {
    }

    MacroLegalizer::MacroLegalizer(NesterovBase *nb_)
    {
        nb_ = nb_;
        // pb_ = nb_->pb_;
        macros = std::vector<GCell*>();
        for(GCell* gcell : nb_->gCells()){
            if(gcell->isMacroInstance()){
                macros.emplace_back(gcell);
            }
        }
        vars_ = MacroLegalizerVars();

    }
    
    MacroLegalizer::~MacroLegalizer()
    {
    }

    void MacroLegalizer::doSimulatedAnnealing(double temp_0, double temp_max, int iter)
    {
        
        //double cooling_rate_ = cooling_rate;
        // 选择一个宏单元
        srand(static_cast<unsigned int>(time(0)));
        int index = rand() % macros.size();
        for (int i = 0; i < vars_.sa_max_iter0; i++) {
            srand(static_cast<unsigned int>(time(0)));
            double temp = temp_0 + (temp_max - temp_0) * (double)(i / vars_.sa_max_iter0);
            GCell* cell= macros[index];
            // 计算宏单元的成本函数
            double old_cost = calc_cost();
            // 随机移动
            std::pair<int, int> move = getRandomMove(cell , iter);
            // 尝试将宏单元移动到新的位置
            cell->setLocation(cell->lx() + move.first, cell->ly() + move.second);
            // 计算宏单元的成本函数
            double new_cost = calc_cost();
            // 判断是否接受新的解
            double delta = (new_cost - old_cost) / old_cost;
            double tau = (double)rand() / RAND_MAX;
            double p = exp(-1.0 * delta / temp);
            if (p > tau) {
                // 接受新的解
                // 判断是否满足要求
                int om = getMacrosOverlap();
                if (om == 0){
                    break;
                }
            } else {
                // 拒绝新的解
                cell->setLocation(cell->lx() - move.first, cell->ly() - move.second);
            }

            // The temperature tj,k at each iteration (j, k) is determined based on the maximum cost 
            // increase fmax(j, k) that will be accepted by more than 50% probabil-
            // -ity, thus we set tj,k = (fmax(j, k)/ln 2). We set fmax(j, 0) (fmax(j, kmax))as0.03×βj (0.0001×βj), 
            // denoting that cost increase by less than 3% (0.01%) at the first (last) SA iteration will be accepted by more than 50% probability. 
            // We initialize fmax(j, k) by fmax(j, 0) and linearly decrease it toward fmax(j, kmax). 
            // The radius rj,k of macro motion range is dependent on both the penalty factor and the amount of macros.
            /////////////////////////////////////////////////////////////////////////
            /*  function accept (programmed by Replace_master)
                int rnd = 0;
                int flg = 0;
                prec dc = (c1 - c0) / c0;
                prec drnd = 0;
                prec exp_val = 0;
                sa_t = 0.0001 * pow(1.5, (prec)iter) / log(2.0);
                rnd = rand();
                drnd = (prec)rnd / RAND_MAX;
                exp_val = exp(-1.0 * dc / sa_t);
                if(drnd < exp_val)
                    flg = 1;
                else
                    flg = 0;
                return flg;
            */

            // 退火
            // temp_ *= cooling_rate_;
        }
    }

    

    std::pair<int, int> MacroLegalizer::getRandomMove(GCell* cell,int iter)
    {
        // 随机种子，使用当前时间
        srand(static_cast<unsigned int>(time(0)));
        int clx = cell->lx();
        int cly = cell->ly();
        int cux = cell->lx();
        int cuy = cell->ly();
        // 获得die的大小
        int dlx = pb_->dies()[0]->dieLx();
        int dly = pb_->dies()[0]->dieLy();
        int dux = pb_->dies()[0]->dieUx();
        int duy = pb_->dies()[0]->dieUy();
        int diewidth = dux - dlx;
        int dieheight = duy - dly;
        // maxMoveX，maxMoveY为die的长宽的1/10 
        int maxMoveX = diewidth / 10;
        int maxMoveY = dieheight / 10;
        /*int moveX = (-maxMoveX) + rand() % (maxMoveX - (-maxMoveX) + 1);
        int moveY = (-maxMoveY) + rand() % (maxMoveY - (-maxMoveY) + 1);*/
        int rndx = rand();
        int rndy = rand();
        double rcof = 0.05 * pow(1.5, (double)iter);
        double rx = (max_sa_r_x - 1)/iter; // 1 = rowheight
        double ry = (max_sa_r_y - 1)/iter;
        int moveX = (int) (rndx / RAND_MAX - 0.5) * rx;
        int moveY = (int) (rndy / RAND_MAX - 0.5) * ry;
        // 计算移动后的 Macro 的四个角坐标
        int newcux = cux + moveX; 
        int newcuy = cuy + moveY;
        int newclx = clx + moveX;
        int newcly = cly + moveY;
        // 检查移动后的坐标是否超过 Die 的边界范围，若超过则修正坐标
        if (newclx < dlx) {
            int diff = dlx - newclx;
            moveX += diff;
        }
        if (newcuy > duy) {
            int diff = newcuy - duy;
            moveY -= diff;
        }
        if (newcux > dux) {
            int diff = newcux - dux;
            moveX -= diff;
        }
        if (newcly < dly) {
            int diff = newcly - dly;
            moveY += diff;
        }
        return std::make_pair(moveX, moveY);
    }

    void MacroLegalizer::doMacroLegalization()
    {
        
        for (int iter = 0; iter < vars_.sa_max_iter; iter++) {
        //     // 初始化参数
            double cooling_rate;
            double temp_0 = 0.03 * pow(1.5, (double)iter) / log(2.0);
            double temp_max = 0.0001 * pow(1.5, (double)iter) / log(2.0);
            double rcof = 0.05 * pow(1.5, (double) iter);
            max_sa_r_x = 1 / sqrt(macros.size()) * rcof; // 1 = total place region
            max_sa_r_y = 1 / sqrt(macros.size()) * rcof;
        
            doSimulatedAnnealing(temp_0, temp_max, iter);
            
        //     // 选择一个宏单元
        //     for (int j = 0; j < num_macros; j++) {
        //         // 计算宏单元的成本函数
        //         double old_cost = calc_cost(j);

        //         // 尝试将宏单元移动到新的位置
        //         move_macro(j);

        //         // 计算宏单元的成本函数
        //         double new_cost = calc_cost(j);

        //         // 判断是否接受新的解
        //         if (accept(new_cost, old_cost, temp)) {
        //             // 接受新的解
        //             update_cost(j, new_cost);
        //         } else {
        //             // 拒绝新的解
        //             undo_move(j);
        //         }
        //     }

        //     // 退火
        //     temp *= cooling_rate;
        // }
    // return std::make_pair(0, 0);
    }
  }

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
    std::vector<Instance *> insts;
    for (Die *die : pb_->dies())
    {
      insts.clear();
      for (Instance *inst : die->insts())
      {
        // regardless of fixed instance
        if (inst->isMacro())
        {
          insts.push_back(inst);
        }
      }
      postLegalize(insts, die);
    }
  }

  // 计算两个矩形的重叠面积
  int MacroLegalizer::overlapArea(GCell *cell1, GCell *cell2)
  {
    int x_overlap = std::max(0, std::min(cell1->ux(), cell2->ux()) - std::max(cell2->lx(), cell2->lx()));
    int y_overlap = std::max(0, std::min(cell1->uy(), cell2->uy()) - std::max(cell2->ly(), cell2->ly()));
    return x_overlap * y_overlap;
  }

  int MacroLegalizer::getCellMacroOverlap()
  {
    int totalOverlap = 0;
    for (auto &cell : nb_->gCells())
    {
        int totalOverlap = 0;
        // 外层macro里层stdcell
        for(int i=0; i<nb_->gCells().size();i++){
            if(!nb_->gCells()[i]->isMacroInstance()){
                continue;
            }else{
                for(int j=i+1; j<nb_->gCells().size();j++){
                    if(nb_->gCells()[j]->isStdInstance()){
                        totalOverlap += overlapArea(nb_->gCells()[i], nb_->gCells()[j]);
                    }else{
                        continue;
                    }
                }
            }
        }
        return totalOverlap;
    }
  }

    int MacroLegalizer::getMacrosOverlap()
    {
        int totalOverlap = 0;
        for(int i=0; i<nb_->gCells().size();i++){
            if(!nb_->gCells()[i]->isMacroInstance()){
                continue;
            }else{
                for(int j=i+1; j<nb_->gCells().size();j++){
                    if(nb_->gCells()[j]->isMacroInstance()){
                        totalOverlap += overlapArea(nb_->gCells()[i], nb_->gCells()[j]);
                    }else{
                        continue;
                    }
                }
            }
        }
        return totalOverlap;
    }

    // 计算宏单元的成本函数
    double MacroLegalizer::calc_cost()
    {
        // HPWL
        double hpwl = nb_->getHpwl();
        // 宏单元与标准单元重叠
        double den= getCellMacroOverlap();
        // 宏单元与宏单元重叠
        double ov = getMacrosOverlap();
        return vars_.sa_hpwl_wgt * hpwl + vars_.sa_den_wgt * den + vars_.sa_ovlp_wgt * ov;
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
                    if(i == j)
                        continue;

                    // check overlap
                    int lx = std::max(inst1->lx(), inst2->lx());
                    int ly = std::max(inst1->ly(), inst2->ly());
                    int ux = std::min(inst1->ux(), inst2->ux());
                    int uy = std::min(inst1->uy(), inst2->uy());

                    if(lx < ux && ly < uy)
                    {
                        isLegal = false;

                        int dx1, dx2, dy1, dy2;
                        // x overlap
                        bool xOk = repel(inst1->lx(), inst1->ux(), inst2->lx(), inst2->ux(),
                                         die->coreLx(), die->coreUx(), &dx1, &dx2);
                        // y overlap
                        auto yOk = repel(inst1->ly(), inst1->uy(), inst2->ly(), inst2->uy(),
                                         die->coreLy(), die->coreUy(), &dy1, &dy2);

                        float thres = (float)(dx1 + dx2) / (float)(dy1 + dy2 + dx1 + dx2);

                        if(((float)rand() / RAND_MAX) < thres)
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