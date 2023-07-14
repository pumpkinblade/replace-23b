#include "abacusLegalizer.h"
#include "placerBase.h"
#include "plot.h"
#include "log.h"
#include <algorithm>
#include <string>

namespace replace
{
  ////////////////////////////////////////////
  // AbacuseLegalizerVars

  AbacusLegalizerVars::AbacusLegalizerVars()
      : weightOpt(One) {}

  ////////////////////////////////////////////
  // AbacusCell

  AbacusCell::AbacusCell(Instance* inst)
  {
    inst_ = inst;
    gpLx_ = static_cast<float>(inst->lx());
    gpLy_ = static_cast<float>(inst->ly());
    width_ = static_cast<float>(inst->dx());
    height_ = static_cast<float>(inst->dy());
    weight_ = 1.0f;
  }

  /////////////////////////////////////////////
  // AbacusCluster
  
  AbacusCluster::AbacusCluster()
      : ec_(0.f), qc_(0.f), wc_(0.f)
  {
  }

  void AbacusCluster::addCell(AbacusCell* cell)
  {
    cells_.push_back(cell);
    ec_ += cell->weight();
    qc_ += cell->weight() * (cell->lgLx() - wc_);
    wc_ += cell->width();
  }

  void AbacusCluster::addCluster(AbacusCluster* cluster)
  {
    cells_.insert(cells_.end(), cluster->cells_.begin(), cluster->cells_.end());
    ec_ += cluster->ec_;
    qc_ += cluster->qc_ - cluster->ec_ * wc_;
    wc_ += cluster->wc_;
    cluster->reset();
  }

  void AbacusCluster::place()
  {
    float x = xc_;
    for(AbacusCell* cell : cells_)
    {
      cell->setLgLx(x);
      x += cell->width();
    }
  }

  void AbacusCluster::reset()
  {
    ec_ = qc_ = wc_ = 0.0f;
    cells_.clear();
  }

  /////////////////////////////////////////////
  // AbacusRow

  AbacusRow::AbacusRow()
    : lx_(0.0f), ly_(0.0f), width_(0.0f), height_(0.0f), usedWidth_(0.0f)
  {
  }

  AbacusRow::AbacusRow(float lx, float ly, float w, float h)
    : lx_(lx), ly_(ly), width_(w), height_(h), usedWidth_(0.0f)
  {
  }

  void AbacusRow::pushCell(AbacusCell* cell)
  {
    cells_.push_back(cell);
    cell->setLgLy(ly_);
    usedWidth_ += cell->width();
  }

  void AbacusRow::popCell()
  {
    usedWidth_ -= cells_.back()->width();
    cells_.pop_back();
  }

  void AbacusRow::placeRow()
  {
    clusters_.clear();
    for(AbacusCell* cell : cells_)
    {
      AbacusCluster* cluster = clusters_.size() == 0 ? nullptr : &clusters_.back();
      if(cluster == nullptr || cluster->xc() + cluster->wc() <= cell->gpLx())
      {
        clusters_.emplace_back();
        clusters_.back().addCell(cell);
        clusters_.back().setXc(cell->gpLx());
      }
      else
      {
        cluster->addCell(cell);
        collapse();
      }
    }

    for (AbacusCluster& cluster : clusters_)
    {
        cluster.place();
    }
  }

  void AbacusRow::collapse()
  {
    bool needCollapse = clusters_.size() > 0;
    while(needCollapse)
    {
      AbacusCluster* curr = &clusters_[clusters_.size() - 1];
      float xc = curr->qc() / curr->ec();
      xc = std::max(xc, lx_);
      xc = std::min(xc, lx_ + width_ - curr->wc());
      curr->setXc(xc);

      AbacusCluster* prev = clusters_.size() <= 1 ? nullptr : &clusters_[clusters_.size() - 2];
      if(prev != nullptr && prev->xc() + prev->wc() > xc)
      {
        prev->addCluster(curr);
        clusters_.pop_back();
        needCollapse = true;
      }
      else
      {
        needCollapse = false;
      }
    }
  }

  ////////////////////////////////////////////
  // AbacusLegalizer

  AbacusLegalizer::AbacusLegalizer(AbacusLegalizerVars lgVars, std::shared_ptr<PlacerBase> pb)
  {
    lgVars_ = lgVars;
    pb_ = pb;

    generateCells();
    generateRows();
  }

  void AbacusLegalizer::doLegalization()
  {
    std::sort(cells_.begin(), cells_.end(), [](const AbacusCell *left, const AbacusCell *right)
              { return left->gpLx() < right->gpLx(); });
    for(AbacusCell* cell : cells_)
    {
      double cbest = std::numeric_limits<double>::infinity();
      AbacusRow* rbest = nullptr;
      for(AbacusRow* row : rows_)
      {
        double clowerBound = cell->gpLy() - row->ly();
        clowerBound *= clowerBound;
        if(clowerBound > cbest)
          continue;
        if(row->usedWidth() + cell->width() > row->width())
          continue;

        row->pushCell(cell);
        row->placeRow();
        double dx = cell->gpLx() - cell->lgLx();
        double dy = cell->gpLy() - cell->lgLy();
        double c = dx * dx + dy * dy;
        if(c < cbest)
        {
          cbest = c;
          rbest = row;
        }
        row->popCell();
      }
      rbest->pushCell(cell);
      rbest->placeRow();
    }

#ifdef ENABLE_CIMG_LIB
    PlotEnv pe;
    pe.setPlacerBase(pb_);
    pe.Init();
    pe.SaveCellPlotAsJPEG(std::string("BeforeLegalization"),
        false, std::string("./plot/cell/before_lg"));
#endif

    int64_t hpwlBeforeLG = pb_->hpwl();
    LOG_INFO("hpwl Before AbacusLegalization: {}", hpwlBeforeLG);
    for(AbacusCell* cell : cells_)
    {
      int lx = static_cast<int>(cell->lgLx() + 0.5f);
      int ly = static_cast<int>(cell->lgLy() + 0.5f);
      cell->instance()->setLocation(lx, ly);
    }
    int64_t hpwlAfterLG = pb_->hpwl();
    LOG_INFO("hpwl After AbacusLegalization: {}", hpwlAfterLG);

#ifdef ENABLE_CIMG_LIB
    pe.SaveCellPlotAsJPEG(std::string("AfterLegalization"),
                          false, std::string("./plot/cell/after_lg"));
#endif
  }

  void AbacusLegalizer::generateCells()
  {
    // TODO: Macros and fixed instances should not be treated as AbacusCell

    cellStor_.reserve(pb_->insts().size());
    for(Instance* inst : pb_->insts())
    {
      cellStor_.emplace_back(inst);
      switch (lgVars_.weightOpt)
      {
      case AbacusLegalizerVars::One:
        cellStor_.back().setWeight(1.0f);
        break;
      case AbacusLegalizerVars::Area:
        cellStor_.back().setWeight(static_cast<float>(inst->dx()) * inst->dy());
        break;
      case AbacusLegalizerVars::NumPins:
        cellStor_.back().setWeight(static_cast<float>(inst->pins().size()));
      default:
        break;
      }
    }

    cells_.resize(cellStor_.size());
    for(size_t i = 0; i < cells_.size(); i++)
      cells_[i] = &cellStor_[i];
  }

  void AbacusLegalizer::generateRows()
  {
    // TODO: If a row is blocked by macros or fixed instances, we should split it

    rowStor_.reserve(pb_->die()->rowRepeatCount());
    for(int i = 0; i < pb_->die()->rowRepeatCount(); i++)
    {
      rowStor_.emplace_back(static_cast<float>(pb_->die()->rowStartX()),
                            static_cast<float>(pb_->die()->rowStartY() + i * pb_->die()->rowHeight()),
                            static_cast<float>(pb_->die()->rowWidth()),
                            static_cast<float>(pb_->die()->rowHeight()));
      rows_.push_back(&rowStor_.back());
    }
  }
}