#include "abaxLegalizer.h"
#include "placerBase.h"
#include "plot.h"
#include "log.h"
#include <algorithm>
#include <string>
#include <cassert>

namespace replace
{
  /////////////////////////////////////////////
  // AbaxCluster

  AbaxCluster::AbaxCluster()
      : startIdx_(0), endIdx_(0), xc_(0.f), ec_(0.f), qc_(0.f), wc_(0.f)
  {
  }

  AbaxCluster::AbaxCluster(int startIdx)
      : AbaxCluster()
  {
    startIdx_ = startIdx;
    endIdx_ = startIdx;
  }

  void AbaxCluster::addCell(const AbaxCell& cell)
  {
    endIdx_++;
    const float weight = 1.0f;
    ec_ += weight;
    qc_ += weight * (cell.gpLx() - wc_);
    wc_ += cell.width();
  }

  void AbaxCluster::addCluster(AbaxCluster& other)
  {
    // Note: two clusters can add togather only if they are adjacent
    endIdx_ = other.endIdx_;
    ec_ += other.ec_;
    qc_ += other.qc_ - other.ec_ * wc_;
    wc_ += other.wc_;
    other.reset();
  }

  void AbaxCluster::reset()
  {
    xc_ = ec_ = qc_ = wc_ = 0.f;
    startIdx_ = endIdx_ = 0;
  }

  /////////////////////////////////////////////
  // AbaxSubrow

  AbaxSubrow::AbaxSubrow()
    : lx_(0.0f), ux_(0.0f), usedWidth_(0.0f)
  {
    reset();
  }

  AbaxSubrow::AbaxSubrow(float lx, float ux)
    : lx_(lx), ux_(ux), usedWidth_(0.0f)
  {
    reset();
  }

  void AbaxSubrow::reset()
  {
    cellStartIdx_ = cellEndIdx_ = 0;
    clusters_.clear();
    usedWidth_ = 0.f;
  }

  void AbaxSubrow::addCell(int cellIdx, std::vector<AbaxCell>& cells)
  {
    // push the cell
    if (cellStartIdx_ >= cellEndIdx_)
      cellStartIdx_ = cellEndIdx_ = cellIdx;
    cellEndIdx_++;
    usedWidth_ += cells[cellIdx].width();

    appendCell(clusters_, cellIdx, cells);

    // place the last cluster
    AbaxCluster* last = &clusters_.back();
    float x = last->xc();
    for(int i = last->startIdx(); i < last->endIdx(); i++)
    {
      cells[i].setLgLx(x);
      x += cells[i].width();
    }
  }

  void AbaxSubrow::appendCell(std::vector<AbaxCluster>& clusters, int cellIdx, std::vector<AbaxCell>& cells)
  {
    // find last cluster
    AbaxCluster* last = clusters.size() == 0 ? nullptr : &clusters.back();
    float xc = cells[cellIdx].gpLx();
    xc = std::max(xc, lx_);
    xc = std::min(xc, ux_ - cells[cellIdx].width());
    if(last == nullptr || last->xc() + last->wc() <= xc)
    {
      // create a new cluster
      clusters.emplace_back(static_cast<int>(cellIdx));
      clusters.back().addCell(cells[cellIdx]);
      clusters.back().setXc(xc);
    }
    else
    {
      last->addCell(cells[cellIdx]);
      collapse(clusters);
    }
  }

  void AbaxSubrow::collapse(std::vector<AbaxCluster>& clusters)
  {
    bool needCollapse = clusters.size() > 0;
    while(needCollapse)
    {
      AbaxCluster& curr = clusters[clusters.size() - 1];
      float xc = curr.qc() / curr.ec();
      xc = std::max(xc, lx_);
      xc = std::min(xc, ux_- curr.wc());
      curr.setXc(xc);

      AbaxCluster* prev = clusters.size() <= 1 ? nullptr : &clusters[clusters.size() - 2];
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
  // AbaxRow

  AbaxRow::AbaxRow()
      : lx_(0.f), ly_(0.f), ux_(0.f), uy_(0.f)
  {
  }

  AbaxRow::AbaxRow(float lx, float ly, float ux, float uy)
      : lx_(lx), ly_(ly), ux_(ux), uy_(uy)
  {
  }

  void AbaxRow::addInstance(Instance* inst, int subrowIdx)
  {
    cells_.emplace_back(inst);
    int cellIdx = static_cast<int>(cells_.size() - 1);
    recursiveMoves(cellIdx, cellIdx + 1, cells_, subrowIdx, subrows_);
  }

  void AbaxRow::tryAddInstance(Instance* inst, int& subrowIdx, float& cost)
  {
    cost = std::numeric_limits<float>::infinity();

    // make a copy
    trySubrows_.assign(subrows_.begin(), subrows_.end());
    tryCells_.assign(cells_.begin(), cells_.end());

    // push the cell
    tryCells_.emplace_back(inst);
    int cellIdx = static_cast<int>(tryCells_.size() - 1);

    int idx = -1;
    // find the last nonempty subrow
    for(idx = trySubrows_.size() - 1; idx >= 0; idx--)
    {
      if(trySubrows_[idx].cellStartIdx() < trySubrows_[idx].cellEndIdx())
        break;
    }

    if(idx >= 0)
    {
      // try to place cell in nonempty subrows_[idx]
      if(recursiveMoves(cellIdx, cellIdx + 1, tryCells_, idx, trySubrows_))
      {
        cost = std::abs(tryCells_.back().gpLx() - tryCells_.back().lgLx())
               + std::abs(tryCells_.back().gpLy() - ly_);
        subrowIdx = idx;
      }
    }

    // find next subrow which can contain the cell
    for(idx += 1; idx < trySubrows_.size(); idx++)
    {
      // subrow must be empty
      if(trySubrows_[idx].width() < tryCells_.back().width())
        continue;
      // try to place cell in subrows_[idx]
      trySubrows_[idx].addCell(cellIdx, tryCells_);
      // calculate the cost
      float c = std::abs(tryCells_.back().gpLx() - tryCells_.back().lgLx()) 
                + std::abs(tryCells_.back().gpLy() - ly_);
      if(c < cost)
      {
        cost = c;
        subrowIdx = idx;
      }
    }
  }

  float AbaxRow::getCost(const std::vector<AbaxCell> cells) const
  {
    float cost = 0.0f;
    for(const AbaxCell& cell : cells)
      cost += std::abs(cell.gpLx() - cell.lgLx()) + std::abs(cell.gpLy() - ly_);
    return cost;
  }

  void AbaxRow::applyLegalization()
  {
    for(AbaxCell& cell :  cells_)
    {
      int lx = static_cast<int>(cell.lgLx() + 0.5f);
      int ly = static_cast<int>(ly_ + 0.5f);
      cell.instance()->setLocation(lx, ly);
    }
  }

  bool AbaxRow::recursiveMoves(int groupStartIdx, int groupEndIdx, std::vector<AbaxCell>& cells,
                               int subrowIdx, std::vector<AbaxSubrow>& subrows)
  {
    // empty group
    if(groupStartIdx >= groupEndIdx)
      return true;

    // calculate the required space
    float reqSpace = 0.0f;
    for(int i = groupStartIdx; i < groupEndIdx; i++)
      reqSpace += cells[i].width();

    // find the front cells that should move to the previous subrow
    AbaxSubrow& subrow = subrows[subrowIdx];
    if(!subrow.isEmpty())
      assert(subrow.cellEndIdx() == groupStartIdx);
    float currSpace = subrow.width() - subrow.usedWidth();
    int moveGroupStartIdx = subrow.cellStartIdx();
    int moveGroupEndIdx = subrow.cellStartIdx();
    while(currSpace < reqSpace && moveGroupEndIdx < subrow.cellEndIdx())
    {
      currSpace += cells[moveGroupEndIdx].width();
      moveGroupEndIdx++;
    }
    
    // the whole subrow couldn't contain this cell
    if(currSpace < reqSpace)
      return false;

    // some cells should be removed from this subrow
    if (moveGroupStartIdx < moveGroupEndIdx)
    {
      // there is no pervious subrow
      if (subrowIdx == 0)
        return false;
      // move cells to pervious subrow
      // try to place cell in nonempty subrows_[idx]
      else if (!recursiveMoves(moveGroupStartIdx, moveGroupEndIdx, cells, subrowIdx - 1, subrows))
        return false;
    }

    // some cells have been removed from this subrow
    int newStartIdx, newEndIdx;
    if (moveGroupStartIdx < moveGroupEndIdx)
    {
      subrow.reset();
      newStartIdx = moveGroupEndIdx;
      newEndIdx = groupEndIdx;
    }
    else
    {
      newStartIdx = groupStartIdx;
      newEndIdx = groupEndIdx;
    }

    // redo abacus legalization in this subrow
    for(int i = newStartIdx; i < newEndIdx; i++)
      subrow.addCell(i, cells);
    return true;
  }

  void AbaxRow::generateSubrows(const std::vector<Instance*>& obstacles)
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
  // AbaxLegalizer

  AbaxLegalizer::AbaxLegalizer(std::shared_ptr<PlacerBase> pb)
  {
    pb_ = pb;
  }

  void AbaxLegalizer::doLegalization()
  {
    int64_t hpwlBeforeLG = pb_->hpwl();
    LOG_DEBUG("hpwl Before AbaxLegalization: {}", hpwlBeforeLG);
    Plot::plot(pb_.get(), "./plot/cell", "before_lg");

    for(Die* die : pb_->dies())
    {
      if(!die->isSetRow())
        continue;

      reset(die);
      std::sort(die->insts().begin(), die->insts().end(), [](const Instance* left, const Instance* right)
                { return left->lx() < right->lx(); });
      LOG_TRACE("finish sorting in AbaxLegalizer::doLegalization");
      for(Instance* inst : die->insts())
      {
        if(inst->isFixed() || inst->isMacro())
          continue;

        float cbest = std::numeric_limits<float>::infinity();
        AbaxRow* rbest = nullptr;
        int subrowIdxBest = -1;

        // find nearest row
        int idx = (inst->ly() - die->coreLy()) / die->rowHeight();

        bool ascend = false;
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
          
          AbaxRow& row = rows_[idx];
          if(std::abs(static_cast<float>(inst->ly()) - row.ly()) > cbest)
          {
            touchBound++;
            if(touchBound > 2)
              break;
            else
              continue;
          }

          int subrowIdx;
          float cost;
          row.tryAddInstance(inst, subrowIdx, cost);
          if(cost < cbest)
          {
            cbest = cost;
            rbest = &row;
            subrowIdxBest = subrowIdx;
          }
        }
        if(rbest == nullptr)
          LOG_ERROR("Lack of area. Unable to do place inst `{}` on die `{}`", inst->name(), die->name());
        else
          rbest->addInstance(inst, subrowIdxBest);
      }
      LOG_TRACE("finish row assignment");

      LOG_TRACE("replace instance's location with its legalized location");
      for(AbaxRow& row : rows_)
        row.applyLegalization();
    }

    int64_t hpwlAfterLG = pb_->hpwl();
    LOG_DEBUG("hpwl After AbaxLegalization: {}", hpwlAfterLG);
    Plot::plot(pb_.get(), "./plot/cell", "after_lg");
  }

  void AbaxLegalizer::generateRows()
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

  void AbaxLegalizer::reset(Die* die)
  {
    rows_.clear();
    die_ = die;

    generateRows();
  }
}