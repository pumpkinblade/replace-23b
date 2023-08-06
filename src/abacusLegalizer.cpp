#include "abacusLegalizer.h"
#include "placerBase.h"
#include "plot.h"
#include "log.h"
#include <algorithm>
#include <string>
#include <cassert>

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
    lgLx_ = gpLx_; lgLy_ = gpLy_;
    weight_ = 1.0f;
  }

  /////////////////////////////////////////////
  // AbacusCluster

  AbacusCluster::AbacusCluster()
      : startIdx_(0), endIdx_(0), xc_(0.f), ec_(0.f), qc_(0.f), wc_(0.f)
  {
  }

  AbacusCluster::AbacusCluster(int startIdx)
      : AbacusCluster()
  {
    startIdx_ = startIdx;
    endIdx_ = startIdx;
  }

  void AbacusCluster::addCell(AbacusCell* cell)
  {
    endIdx_++;
    ec_ += cell->weight();
    qc_ += cell->weight() * (cell->lgLx() - wc_);
    wc_ += cell->width();
  }

  void AbacusCluster::addCluster(AbacusCluster* other)
  {
    // Note: two clusters can add togather only if they are adjacent
    endIdx_ = other->endIdx_;
    ec_ += other->ec_;
    qc_ += other->qc_ - other->ec_ * wc_;
    wc_ += other->wc_;
    other->reset();
  }

  void AbacusCluster::reset()
  {
    xc_ = ec_ = qc_ = wc_ = 0.f;
    startIdx_ = endIdx_ = 0;
  }

  /////////////////////////////////////////////
  // AbacusSubrow

  AbacusSubrow::AbacusSubrow()
    : lx_(0.0f), ux_(0.0f), usedWidth_(0.0f)
  {
  }

  AbacusSubrow::AbacusSubrow(float lx, float ux)
    : lx_(lx), ux_(ux), usedWidth_(0.0f)
  {
  }

  void AbacusSubrow::tryAddCell(AbacusCell* cell)
  {
    // push the cell
    cells_.push_back(cell);

    // make a copy
    std::vector<AbacusCluster> copy;
    copy.assign(clusters_.begin(), clusters_.end());

    appendCell(copy, cell);

    // find the lgLx of the cell
    float x = copy.back().xc();
    for(int i = copy.back().startIdx(); i < copy.back().endIdx(); i++)
    {
      x += cells_[i]->width();
    }
    cell->setLgLx(x - cell->width());

    // pop the cell
    cells_.pop_back();
  }

  void AbacusSubrow::addCell(AbacusCell* cell)
  {
    // push the cell
    cells_.push_back(cell);
    usedWidth_ += cell->width();

    appendCell(clusters_, cell);

    // place the last cluster
    AbacusCluster* last = &clusters_.back();
    float x = last->xc();
    for(int i = last->startIdx(); i < last->endIdx(); i++)
    {
      cells_[i]->setLgLx(x);
      x += cells_[i]->width();
    }
  }

  void AbacusSubrow::appendCell(std::vector<AbacusCluster>& clusters, AbacusCell* cell)
  {
    // find last cluster
    AbacusCluster* last = clusters.size() == 0 ? nullptr : &clusters.back();
    float xc = cell->gpLx();
    xc = std::max(xc, lx_);
    xc = std::min(xc, ux_ - cell->width());
    if(last == nullptr || last->xc() + last->wc() <= xc)
    {
      // create a new cluster
      clusters.emplace_back(static_cast<int>(cells_.size() - 1));
      clusters.back().addCell(cell);
      clusters.back().setXc(xc);
    }
    else
    {
      last->addCell(cell);
      collapse(clusters);
    }
  }

  void AbacusSubrow::collapse(std::vector<AbacusCluster>& clusters)
  {
    bool needCollapse = clusters.size() > 0;
    while(needCollapse)
    {
      AbacusCluster* curr = &clusters[clusters.size() - 1];
      float xc = curr->qc() / curr->ec();
      xc = std::max(xc, lx_);
      xc = std::min(xc, ux_- curr->wc());
      curr->setXc(xc);

      AbacusCluster* prev = clusters.size() <= 1 ? nullptr : &clusters[clusters.size() - 2];
      if(prev != nullptr && prev->xc() + prev->wc() > xc)
      {
        prev->addCluster(curr);
        clusters.pop_back();
        needCollapse = true;
      }
      else
      {
        needCollapse = false;
      }
    }
  }

  /////////////////////////////////////////////
  // AbacusRow

  AbacusRow::AbacusRow()
      : lx_(0.f), ly_(0.f), ux_(0.f), uy_(0.f)
  {
  }

  AbacusRow::AbacusRow(float lx, float ly, float ux, float uy)
      : lx_(lx), ly_(ly), ux_(ux), uy_(uy)
  {
  }

  void AbacusRow::generateSubrows(const std::vector<Instance*>& obstacles)
  {
    // assume obs has been sorted
    int currLx = static_cast<int>(lx_);
    for(Instance* obs : obstacles)
    {
      // find overlap
      int lx = std::max(currLx, obs->lx());
      int ly = std::max(static_cast<int>(ly_), obs->ly());
      int ux = std::min(static_cast<int>(ux_), obs->ux());
      int uy = std::min(static_cast<int>(uy_), obs->uy());

      // overlap exists
      if(lx < ux && ly < uy)
      {
        // |        |        |******|     |
        // ^        ^        ^      ^     ^
        // rowLx    currLx   lx     ux    rowUx
        if(currLx < lx)
        {
          // create subrow
          subrows_.emplace_back(static_cast<float>(currLx),
                                static_cast<float>(lx));
        }
        // update currLx
        currLx = ux;
      }
    }
    
    if(currLx < ux_)
    {
      subrows_.emplace_back(static_cast<float>(currLx),
                            static_cast<float>(ux_));
    }
  }

  ////////////////////////////////////////////
  // AbacusLegalizer

  AbacusLegalizer::AbacusLegalizer(AbacusLegalizerVars lgVars, std::shared_ptr<PlacerBase> pb)
  {
    lgVars_ = lgVars;
    pb_ = pb;
  }

  void AbacusLegalizer::doLegalization()
  {
    int64_t hpwlBeforeLG = pb_->hpwl();
    LOG_DEBUG("hpwl Before AbacusLegalization: {}", hpwlBeforeLG);
    Plot::plot(pb_.get(), "./plot/cell", "before_lg");

    for(Die* die : pb_->dies())
    {
      if(!die->isSetRow())
        continue;

      reset(die);
      std::sort(cellStor_.begin(), cellStor_.end(), [](const AbacusCell& left, const AbacusCell& right)
                { return left.gpLx() < right.gpLx(); });
      LOG_TRACE("finish sorting in AbacusLegalizer::doLegalization");
      for(AbacusCell& cell : cellStor_)
      {
        float cbest = std::numeric_limits<float>::infinity();
        AbacusSubrow* srbest = nullptr;
        AbacusRow* rbest = nullptr;

        // find nearest row
        AbacusRow tmp(cell.gpLx(), cell.gpLy(), 0.f, 0.f);
        auto it = std::upper_bound(rows_.begin(), rows_.end(), tmp, 
                                   [](const AbacusRow& left, const AbacusRow& right)
                                   { return left.lx() < right.lx(); });
        int idx = static_cast<int>(it - rows_.begin());
        bool ascend = true;
        int step = 0;
        int touchBound = 0;
        while(true)
        {
          idx = idx + (ascend ? step : -step);
          step += 1;
          ascend = !ascend;
          if(idx < 0 || idx >= rows_.size())
          {
            int nextIdx = idx + (ascend ? step : -step);
            if(nextIdx < 0 || nextIdx >= rows_.size())
              break;
            else
              continue;
          }
          
          AbacusRow& row = rows_[idx];
          if(std::abs(cell.gpLy() - row.ly()) > cbest)
          {
            touchBound++;
            if(touchBound > 2)
              break;
            else
              continue;
          }

          for(AbacusSubrow& subrow : row.subrows())
          {
            if(subrow.usedWidth() + cell.width() > subrow.width())
              continue;

            cell.setLgLy(row.ly());
            subrow.tryAddCell(&cell);
            float c = std::abs(cell.gpLx() - cell.lgLx()) + std::abs(cell.gpLy() - cell.lgLy());
            if(c < cbest)
            {
              cbest = c;
              rbest = &row;
              srbest = &subrow;
            }
          }
        }
        if(rbest == nullptr)
        {
          LOG_ERROR("Lack of area. Unable to do stdcell legalization on die `{}`", die->name());
          break;
        }
        cell.setLgLy(rbest->ly());
        srbest->addCell(&cell);
      }
      LOG_TRACE("finish row assignment");

      LOG_TRACE("replace instance's location with its legalized location");
      for(AbacusCell& cell : cellStor_)
      {
        int lx = static_cast<int>(cell.lgLx() + 0.5f);
        int ly = static_cast<int>(cell.lgLy() + 0.5f);
        cell.instance()->setLocation(lx, ly);
      }
    }

    int64_t hpwlAfterLG = pb_->hpwl();
    LOG_DEBUG("hpwl After AbacusLegalization: {}", hpwlAfterLG);
    Plot::plot(pb_.get(), "./plot/cell", "after_lg");
  }

  void AbacusLegalizer::generateCells()
  {
    // Macros and fixed instances should not be treated as AbacusCell
    cellStor_.reserve(die_->insts().size());
    for(Instance* inst : die_->insts())
    {
      if(inst->isMacro() || inst->isFixed())
        continue;
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
  }

  void AbacusLegalizer::generateRows()
  {
    // TODO: If a row is blocked by macros or fixed instances, we should split it
    std::vector<Instance*> obstacles;
    for(Instance* inst : die_->insts())
    {
      if(inst->isFixed() || inst->isMacro())
      {
        obstacles.push_back(inst);
      }
    }
    // sort obstacles by lx
    std::sort(obstacles.begin(), obstacles.end(), [](const Instance* left, const Instance* right)
              { return left->lx() < right->lx(); });

    rows_.reserve(die_->rowRepeatCount());
    for(int i = 0; i < die_->rowRepeatCount(); i++)
    {
      int rowLx = die_->rowStartX();
      int rowLy = die_->rowStartY() + i * die_->rowHeight();
      int rowUx = die_->rowStartX() + die_->rowWidth();
      int rowUy = die_->rowStartY() + (i + 1) * die_->rowHeight();
      rows_.emplace_back(static_cast<float>(rowLx),
                         static_cast<float>(rowLy),
                         static_cast<float>(rowUx),
                         static_cast<float>(rowUy));
      rows_.back().generateSubrows(obstacles);
    }
  }

  void AbacusLegalizer::reset(Die* die)
  {
    rows_.clear();
    cellStor_.clear();
    die_ = die;

    generateCells();
    generateRows();
  }
}