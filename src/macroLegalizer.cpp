#include "macroLegalizer.h"
#include "placerBase.h"
#include "nesterovBase.h"
#include <algorithm>
#include <random>
#include "point.h"

namespace replace
{
    std::pair<int, int> adjust(int l1, int u1, int l2, int u2, int ll, int uu)
    {
        if (l1 <= l2 && l2 < u1)
        {
        // |--------|===========|-------|
        // ^        ^           ^       ^
        // l1       l2          u1      u2
        if (u1 < u2)
        {
            int d = u1 - l2;
            int d1 = std::min(d / 2, l1 - ll);
            int d2 = std::min(d - d1, uu - u2);
            return std::make_pair(-d1, d2);
        }

        // |--------|===========|-------|
        // ^        ^           ^       ^
        // l1       l2          u2      u1
        else
        {
            int d = std::min(u1 - l2, u2 - l1);
            int d1 = std::min(d / 2, l1 - ll);
            int d2 = std::min(d - d1, uu - u2);
            return std::make_pair(-d1, d2);
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
            return std::make_pair(d1, -d2);
        }

        // |--------|===========|-------|
        // ^        ^           ^       ^
        // l2       l1          u1      u2
        else
        {
            int d = std::min(u1 - l2, u2 - l1);
            int d2 = std::min(d / 2, l2 - ll);
            int d1 = std::min(d - d2, uu - u1);
            return std::make_pair(d1, -d2);
        }
        }

        return std::make_pair(0, 0);
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

    void MacroLegalizer::doSimulatedAnnealing(double temp, double cooling_rate)
    {
        double temp_ = temp;
        double cooling_rate_ = cooling_rate;
        // 选择一个宏单元
        random_shuffle(macros.begin(), macros.end());
        for (int i = 0; i < vars_.sa_max_iter0; i++) {
            GCell* cell= macros[0];
            // 计算宏单元的成本函数
            double old_cost = calc_cost();
            // 随机移动

            // 计算宏单元的成本函数

            // 退火
            temp *= cooling_rate;
        }
    }

    IntPoint MacroLegalizer::getRandomMove(GCell* cell)
    {
        // 暂未完成
        int lx = cell->lx();
        int ly = cell->ly();
        // int x_max = nb_->pb()->die()[0].width() - w;
        // int y_max = nb_->die().height() - h;
        // std::random_device rd;
        // std::mt19937 gen(rd());
        // std::uniform_int_distribution<> dis_x(x, x_max);
        // std::uniform_int_distribution<> dis_y(y, y_max);
        // return IntPoint(dis_x(gen), dis_y(gen));
        return IntPoint(lx, ly);
    }

    void MacroLegalizer::doMacroLegalization()
    {
        // for (int i = 0; i < vars_.max_iter; i++) {
        //     // 初始化参数
            
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

        for (int j = 0; j < i; j++)
        {
          // check overlap
          Instance *inst2 = insts[j];

          // x overlap
          auto xMove = adjust(inst1->lx(), inst1->ux(), inst2->lx(), inst2->ux(),
                              die->coreLx(), die->coreUx());
          inst1->setLocation(inst1->lx() + xMove.first, inst1->ly());
          inst2->setLocation(inst2->lx() + xMove.second, inst2->ly());
          isLegal &= (xMove.first == 0) && (xMove.second == 0);

          // y overlap
          auto yMove = adjust(inst1->ly(), inst1->uy(), inst2->ly(), inst2->uy(),
                              die->coreLy(), die->coreUy());
          inst1->setLocation(inst1->lx(), inst1->ly() + yMove.first);
          inst2->setLocation(inst2->lx(), inst2->ly() + yMove.second);
          isLegal &= (yMove.first == 0) && (yMove.second == 0);
        }
      }

      if (isLegal)
      {
        break;
      }
    }
  }
}