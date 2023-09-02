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
        double sa_ovlp_cof = 1.05; /// need tuning
        double sa_max_iter = 500; /// need tuning
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

        // 计算两个instance之间重叠面积
        double getOverlapArea(const Instance *inst1, const Instance *inst2);
        // 计算一个macro与其余macro的重叠面积
        double getMacroOverlap(const Instance *macro, const std::vector<Instance *> &macros);
        // 计算所有macro之间的重叠面积
        double getAllMacroOverlap(const std::vector<Instance *> &macros);
        // 计算一个macro所影响的hpwl
        double getMacroHpwl(const Instance *macro);
        // 若一个macro的位置更新了，需要对其所在的所有net更新
        void updateMacroNetBox(Instance *macro);

        bool saLegalize2();
        void sa_init_top();
        void sa_param_init_top();
        void sa_mac_leg_top();
        void sa_param_init(int iter);
        void sa_mac_leg(int iter);
        void sa_param_update();
        void sa_mac_leg_sub();
        void sa_param_update_sub();
        int sa_get_mac_idx();
        std::pair<int, int> sa_get_mac_mov(const Instance* mac);
        void sa_do_mac_mov(Instance* mac, const std::pair<int, int>& mov);
        double sa_get_mac_cost(const Instance* mac, double& hpwl, double& den, double& ovlp);
        bool sa_mov_accept(double old_cost, double new_cost);

        void init_bins();
        double get_mac_den(const Instance* mac);

        bool postLegalize(const std::vector<Instance *> &macros, const Die* die);

        bool checkLegal(const std::vector<Instance *> &macros, const Die* die);

    private:
        /* data */
        std::shared_ptr<PlacerBase> pb_;
        std::vector<Instance*> macros_;
        Die* die_;
        MacroLegalizerVars lgVars_;

        // 与搜索半径有关的参数
        double max_sa_r_x;
        double max_sa_r_y;
        double min_sa_r_x;
        double min_sa_r_y;
        double sa_r_x;
        double sa_r_y;
        double sa_r_stp_x;
        double sa_r_stp_y;
        double sa_ncof_x;
        double sa_ncof_y;
        double sa_n_x;
        double sa_n_y;

        // 与接受率有关的参数
        double sa_t;
        double sa_init_t;
        double sa_t_cof;
        double sa_init_neg_rate;
        double sa_last_neg_rate;

        // 与迭代次数有关的参数
        int sa_max_iter;
        int sa_max_iter2;
        int sa_max_iter0;
        int sa_iter_cof;

        // 权重
        double sa_hpwl_wgt;
        double sa_hpwl_cof;
        double sa_den_wgt;
        double sa_den_cof;
        double sa_ovlp_wgt;
        double sa_ovlp_cof;

        double tot_mac_hpwl;
        double tot_mac_den;
        double tot_mac_ovlp;

        bool ovlp_free_flg;

        double bin_size_x;
        double bin_size_y;
        int bin_count_x;
        int bin_count_y;
        std::vector<double> bin_density;
    };

}

#endif