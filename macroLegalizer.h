#ifndef __REPLACE_MACRO_LEGALIZER__
#define __REPLACE_MACRO_LEGALIZER__

#include <vector>
#include <cassert>
#include <memory>
#include "point.h"

namespace replace
{
    class Instance;
    class Die;
    class PlacerBase;
    class NesterovBase;
    class GCell;

    class MacroLegalizerVars
    {
    public:
        double sa_hpwl_wgt = 1.0;
        double sa_den_wgt = 1.0;
        double sa_ovlp_wgt = 1.0;
        double sa_hpwl_cof = 1.0;
        double sa_den_cof = 1.0;
        double sa_ovlp_cof = 1.5;   /// need tuning
        double sa_max_iter = 500;  /// need tuning
        double sa_max_iter0 = 20; /// need tuning
        int maxPostLegalizeIter;

        MacroLegalizerVars();
    };

    class MacroLegalizer
    {
    public:
        MacroLegalizer() = default;
        MacroLegalizer(MacroLegalizerVars lgVars, std::shared_ptr<PlacerBase> pb);
        ~MacroLegalizer() = default;

        void doLegalization();

        void saLegalize(const std::vector<Instance *> macros, Die *die);
        std::pair<int, int> getRandomMove(Instance *cell, int iter, Die *die, int max_sa_r_x, int max_sa_r_y);
        int overlapArea(Instance *inst1, Instance *inst2);
        int getCellMacroOverlap(const std::vector<Instance *> &macros, Die *die);
        int getMacrosOverlap(const std::vector<Instance *> &macros, Die *die);
        int get_hpwl(const std::vector<Instance *> &macros, Die *die);
        double calc_cost(const std::vector<Instance *> &macros, Die *die);

        void postLegalize(const std::vector<Instance *> macros, Die *die);

        ///OVLP///
        void build_seg_tree(int i, int left, int right) ;
        int get_mac_ovlp(const std::vector<Instance *> macros,int index, Die *die);
        int get_all_macro_ovlp(const std::vector<Instance *> &macros ,Die *die);
        void updata(int i, struct NODE b);
        void callen(int i);
    private:
        /* data */
        std::shared_ptr<PlacerBase> pb_;
        std::shared_ptr<NesterovBase> nb_;

        MacroLegalizerVars vars_;
        MacroLegalizerVars lgVars_;
    };

}

#endif