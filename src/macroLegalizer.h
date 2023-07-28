#ifndef __REPLACE_MACRO_LEGALIZER__
#define __REPLACE_MACRO_LEGALIZER__

#include <vector>
#include <cassert>
#include <memory>
#include "point.h"

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
            int maxPostLegalizeIter;

            MacroLegalizerVars();
    };

    class MacroLegalizer
    {
        public:
            MacroLegalizer(/* args */);
            MacroLegalizer(NesterovBase *nb_);
            MacroLegalizer(MacroLegalizerVars lgVars, std::shared_ptr<PlacerBase> pb);
            void doLegalization();
            ~MacroLegalizer();
            void doSimulatedAnnealing(double temp, double cooling_rate);
            std::pair<int, int> getRandomMove(GCell *cell);
            void doMacroLegalization();
            int overlapArea(GCell* , GCell*);
            int getCellMacroOverlap();
            int getMacrosOverlap();
            double calc_cost();

        private:
        /* data */
        void postLegalize(const std::vector<Instance*> insts, Die* die);
            std::shared_ptr<PlacerBase> pb_;
        std::shared_ptr<NesterovBase> nb_;
        std::vector<GCell*> macros;
        MacroLegalizerVars vars_;
        MacroLegalizerVars lgVars_;
    };
    
    





}

#endif