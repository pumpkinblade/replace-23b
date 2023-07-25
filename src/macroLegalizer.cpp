#include "macroLegalizer.h"
#include "placerBase.h"
#include "nesterovBase.h"
#include <algorithm>

namespace replace{
    // using namespace std;


    MacroLegalizerVars::MacroLegalizerVars()
    {
    }

    macroLegalizer::macroLegalizer(/* args */)
    {
    }
    
    macroLegalizer::~macroLegalizer()
    {
    }

    void macroLegalizer::doMacroLegalization()
    {
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
        for(auto& cell : nb_->gCells()){
        }
        return 0;
    }
    int macroLegalizer::getMacrosOverlap()
    {
        return 0;
    }
}