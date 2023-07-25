#include "macroLegalizer.h"
#include "placerBase.h"
#include "nesterovBase.h"
#include <algorithm>
#include <random>
#include "point.h"

namespace replace{
    // using namespace std;


    MacroLegalizerVars::MacroLegalizerVars()
    {

    }

    macroLegalizer::macroLegalizer()
    {
    }

    macroLegalizer::macroLegalizer(NesterovBase *nb_)
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
    
    macroLegalizer::~macroLegalizer()
    {
    }

    void macroLegalizer::doSimulatedAnnealing(double temp, double cooling_rate)
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

    IntPoint macroLegalizer::getRandomMove(GCell* cell)
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

    void macroLegalizer::doMacroLegalization()
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
    }

    // 计算两个矩形的重叠面积
    int macroLegalizer::overlapArea(GCell* cell1, GCell* cell2) {
        int x_overlap = std::max(0, std::min(cell1->ux(), cell2->ux()) - std::max(cell2->lx(), cell2->lx()));
        int y_overlap = std::max(0, std::min(cell1->uy(), cell2->uy()) - std::max(cell2->ly(), cell2->ly()));
        return x_overlap * y_overlap;
    }

    int macroLegalizer::getCellMacroOverlap()
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

    int macroLegalizer::getMacrosOverlap()
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
    double macroLegalizer::calc_cost()
    {
        // HPWL
        double hpwl = nb_->getHpwl();
        // 宏单元与标准单元重叠
        double den= getCellMacroOverlap();
        // 宏单元与宏单元重叠
        double ov = getMacrosOverlap();
        return vars_.sa_hpwl_wgt * hpwl + vars_.sa_den_wgt * den + vars_.sa_ovlp_wgt * ov;
    }
}