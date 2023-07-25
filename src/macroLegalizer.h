#ifndef __REPLACE_MACRO_LEGALIZER__
#define __REPLACE_MACRO_LEGALIZER__

#include <memory>
#include <vector>
#include <cassert>

namespace replace{
    class Instance;
    class Die;
    class PlacerBase;
    class NesterovBase;
    class GCell;

    class MacroLegalizerVars
    {
        public:
            double sa_hpwl_wgt = 1.0;
            double sa_den_wgt;
            double sa_ovlp_wgt;
            double sa_hpwl_cof = 1.0;
            double sa_den_cof = 1.0;
            double sa_ovlp_cof = 1.5;    /// need tuning
            double sa_max_iter0 = 1000;  /// need tuning


            MacroLegalizerVars();
    };

    class macroLegalizer
    {
        public:
            macroLegalizer(/* args */);
            macroLegalizer(NesterovBase *nb_);
            ~macroLegalizer();
            void doSimulatedAnnealing(double temp, double cooling_rate);
            IntPoint getRandomMove(GCell *cell);
            void doMacroLegalization();
            int overlapArea(GCell* , GCell*);
            int getCellMacroOverlap();
            int getMacrosOverlap();
            double calc_cost();

        private:
        /* data */
            std::shared_ptr<PlacerBase> pb_;
            std::shared_ptr<NesterovBase> nb_;
            std::vector<GCell*> macros;
            MacroLegalizerVars vars_;
    };
    
    

}

#endif