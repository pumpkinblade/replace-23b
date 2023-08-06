#include "dreamLegalizer.h"
#include "placerBase.h"
#include "plot.h"
#include "log.h"
#include <algorithm>
#include <string>
#include "dreamplace/function_cpu.h"
#include "dreamplace/abacus_legalize_cpu.h"

namespace replace
{
  ////////////////////////////////////////////
  // DreamLegalizerVars

  DreamLegalizerVars::DreamLegalizerVars()
      : weightOpt(One) {}

  ////////////////////////////////////////////
  // DreamLegalizer

  DreamLegalizer::DreamLegalizer(DreamLegalizerVars vars, std::shared_ptr<PlacerBase> pb)
  {
    vars_ = vars;
    pb_ = pb;
  }

  void DreamLegalizer::doLegalization()
  {
    int64_t hpwlBeforeLG = pb_->hpwl();
    LOG_DEBUG("hpwl Before DreamLegalization: {}", hpwlBeforeLG);

    for(Die* die : pb_->dies())
    {
      if(!die->isSetRow())
        continue;

      std::vector<float> init_x(die->insts().size());
      std::vector<float> init_y(die->insts().size());
      std::vector<float> node_size_x(die->insts().size());
      std::vector<float> node_size_y(die->insts().size());
      std::vector<float> node_weights(die->insts().size());
      std::vector<float> legal_x(die->insts().size());
      std::vector<float> legal_y(die->insts().size());
      float xl = die->coreLx();
      float yl = die->coreLy();
      float xh = die->coreUx();
      float yh = die->coreUy();
      float site_width = std::numeric_limits<float>::infinity();
      float row_height = die->rowHeight();
      int num_bins_x = 1;
      int num_bins_y = 64;
      int num_moveable_nodes = 0;

      int idx = 0;
      for(Instance* inst : die->insts())
      {
        if(inst->isFixed() || inst->isMacro())
          continue;

        inst->setExtId(idx);
        init_x[idx] = inst->lx();
        init_y[idx] = inst->ly();
        node_size_x[idx] = inst->dx();
        node_size_y[idx] = inst->dy();
        node_weights[idx] = 1;
        site_width = std::min(site_width, (float)inst->dx());
        idx++;
      }
      num_moveable_nodes = idx;

      for(Instance* inst : die->insts())
      {
        if(!inst->isFixed() && !inst->isMacro())
          continue;

        inst->setExtId(idx);
        init_x[idx] = inst->lx();
        init_y[idx] = inst->ly();
        node_size_x[idx] = inst->dx();
        node_size_y[idx] = inst->dy();
        node_weights[idx] = 1;
        idx++;
      }

      DreamPlace::LegalizationDB<float> db;
      db.init_x = init_x.data();
      db.init_y = init_y.data();
      db.node_size_x = node_size_x.data();
      db.node_size_y = node_size_y.data();
      db.node_weights = node_weights.data();
      db.flat_region_boxes = nullptr;
      db.flat_region_boxes_start = nullptr;
      db.node2fence_region_map = nullptr;
      db.x = legal_x.data();
      db.y = legal_y.data();
      db.xl = xl;
      db.yl = yl;
      db.xh = xh;
      db.yh = yh;
      db.site_width = site_width;
      db.row_height = row_height;
      db.bin_size_x = static_cast<float>((xh - xl) / num_bins_x);
      db.bin_size_y = static_cast<float>((yh - yl) / num_bins_y);
      db.num_bins_x = num_bins_x;
      db.num_bins_y = num_bins_y;
      db.num_nodes = static_cast<int>(die->insts().size());
      db.num_movable_nodes = num_moveable_nodes;
      db.num_regions = 1;

      DreamPlace::greedyLegalizationCPU(
        db,
        init_x.data(), init_y.data(),
        node_size_x.data(), node_size_y.data(),
        legal_x.data(), legal_y.data(),
        xl, yl, xh, yh,
        site_width, row_height,
        num_bins_x, num_bins_y, static_cast<int>(die->insts().size()),
        num_moveable_nodes
      );

      for (Instance* inst : die->insts())
      {
        if (inst->isFixed() || inst->isMacro())
          continue;
        int lgLx = legal_x[inst->extId()];
        int lgLy = legal_y[inst->extId()];
        inst->setLocation(lgLx, lgLy);
      }
      LOG_DEBUG("hpwl After greedyLegaization: {}", pb_->hpwl());
      Plot::plot(pb_.get(), "./plot/cell", "greedy");

      std::copy(legal_x.begin(), legal_x.end(), init_x.begin());
      std::copy(legal_y.begin(), legal_y.end(), init_y.begin());

      DreamPlace::abacusLegalizationCPU(
        init_x.data(), init_y.data(),
        node_size_x.data(), node_size_y.data(), node_weights.data(),
        legal_x.data(), legal_y.data(),
        xl, yl, xh, yh,
        site_width, row_height,
        num_bins_x, 64, static_cast<int>(die->insts().size()),
        num_moveable_nodes
      );

      for(Instance* inst : die->insts())
      {
        if(inst->isFixed() || inst->isMacro())
          continue;
        int lgLx = legal_x[inst->extId()];
        int lgLy = legal_y[inst->extId()];
        inst->setLocation(lgLx, lgLy);
      }
      LOG_DEBUG("hpwl After abacusLegaization: {}", pb_->hpwl());
      Plot::plot(pb_.get(), "./plot/cell", "abacus");
    }
  }
}