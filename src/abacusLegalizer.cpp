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
  
  // AbacusCluster::AbacusCluster()
  //     : ec_(0.f), qc_(0.f), wc_(0.f)
  // {
  // }

  // void AbacusCluster::addCell(AbacusCell* cell)
  // {
  //   cells_.push_back(cell);
  //   ec_ += cell->weight();
  //   qc_ += cell->weight() * (cell->lgLx() - wc_);
  //   wc_ += cell->width();
  // }

  // void AbacusCluster::addCluster(AbacusCluster* cluster)
  // {
  //   cells_.insert(cells_.end(), cluster->cells_.begin(), cluster->cells_.end());
  //   ec_ += cluster->ec_;
  //   qc_ += cluster->qc_ - cluster->ec_ * wc_;
  //   wc_ += cluster->wc_;
  //   cluster->reset();
  // }

  // void AbacusCluster::place()
  // {
  //   float x = xc_;
  //   for(AbacusCell* cell : cells_)
  //   {
  //     cell->setLgLx(x);
  //     x += cell->width();
  //   }
  // }

  // void AbacusCluster::reset()
  // {
  //   ec_ = qc_ = wc_ = 0.0f;
  //   cells_.clear();
  // }

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
  // AbacusRow

  // AbacusRow::AbacusRow()
  //   : lx_(0.0f), ly_(0.0f), width_(0.0f), height_(0.0f), usedWidth_(0.0f)
  // {
  // }

  // AbacusRow::AbacusRow(float lx, float ly, float w, float h)
  //   : lx_(lx), ly_(ly), width_(w), height_(h), usedWidth_(0.0f)
  // {
  // }

  // void AbacusRow::pushCell(AbacusCell* cell)
  // {
  //   cells_.push_back(cell);
  //   cell->setLgLy(ly_);
  //   usedWidth_ += cell->width();
  // }

  // void AbacusRow::popCell()
  // {
  //   usedWidth_ -= cells_.back()->width();
  //   cells_.pop_back();
  // }

  // void AbacusRow::placeRow()
  // {
  //   clusterStor_.clear();
  //   for(AbacusCell* cell : cells_)
  //   {
  //     AbacusCluster* cluster = clusterStor_.size() == 0 ? nullptr : &clusterStor_.back();
  //     if(cluster == nullptr || cluster->xc() + cluster->wc() <= cell->gpLx())
  //     {
  //       // create a new cluster
  //       clusterStor_.emplace_back();
  //       clusterStor_.back().addCell(cell);
  //       clusterStor_.back().setXc(cell->gpLx());
  //     }
  //     else
  //     {
  //       cluster->addCell(cell);
  //       collapse();
  //     }
  //   }

  //   // place cluster
  //   for(AbacusCluster& cluster : clusterStor_)
  //   {
  //     cluster.place();
  //   }
  // }

  // void AbacusRow::collapse()
  // {
  //   bool needCollapse = clusterStor_.size() > 0;
  //   while(needCollapse)
  //   {
  //     AbacusCluster* curr = &clusterStor_[clusterStor_.size() - 1];
  //     float xc = curr->qc() / curr->ec();
  //     xc = std::max(xc, lx_);
  //     xc = std::min(xc, lx_ + width_ - curr->wc());
  //     curr->setXc(xc);

  //     AbacusCluster* prev = clusterStor_.size() <= 1 ? nullptr : &clusterStor_[clusterStor_.size() - 2];
  //     if(prev != nullptr && prev->xc() + prev->wc() > xc)
  //     {
  //       prev->addCluster(curr);
  //       clusterStor_.pop_back();
  //       needCollapse = true;
  //     }
  //     else
  //     {
  //       needCollapse = false;
  //     }
  //   }
  // }

  AbacusRow::AbacusRow()
    : lx_(0.0f), ly_(0.0f), width_(0.0f), height_(0.0f), usedWidth_(0.0f)
  {
  }

  AbacusRow::AbacusRow(float lx, float ly, float w, float h)
    : lx_(lx), ly_(ly), width_(w), height_(h), usedWidth_(0.0f)
  {
  }

  void AbacusRow::tryAddCell(AbacusCell* cell)
  {
    // push the cell
    cell->setLgLy(ly_);
    cells_.push_back(cell);

    // make a copy
    std::vector<AbacusCluster> copy;
    copy.assign(clusterStor_.begin(), clusterStor_.end());

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

  void AbacusRow::addCell(AbacusCell* cell)
  {
    // push the cell
    cell->setLgLy(ly_);
    cells_.push_back(cell);
    usedWidth_ += cell->width();

    appendCell(clusterStor_, cell);

    // place the last cluster
    AbacusCluster* last = &clusterStor_.back();
    float x = last->xc();
    for(int i = last->startIdx(); i < last->endIdx(); i++)
    {
      cells_[i]->setLgLx(x);
      x += cells_[i]->width();
    }
  }

  void AbacusRow::appendCell(std::vector<AbacusCluster>& clusters, AbacusCell* cell)
  {
    // find last cluster
    AbacusCluster* last = clusters.size() == 0 ? nullptr : &clusters.back();
    if(last == nullptr || last->xc() + last->wc() <= cell->gpLx())
    {
      // create a new cluster
      clusters.emplace_back(static_cast<int>(cells_.size() - 1));
      clusters.back().addCell(cell);
      clusters.back().setXc(cell->gpLx());
    }
    else
    {
      last->addCell(cell);
      collapse(clusters);
    }
  }

  void AbacusRow::collapse(std::vector<AbacusCluster>& clusters)
  {
    bool needCollapse = clusters.size() > 0;
    while(needCollapse)
    {
      AbacusCluster* curr = &clusters[clusters.size() - 1];
      float xc = curr->qc() / curr->ec();
      xc = std::max(xc, lx_);
      xc = std::min(xc, lx_ + width_ - curr->wc());
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

  void AbacusRow::findOverlap() const
  {
    for(AbacusCell* cell1 : cells_)
    {
      float lx1 = cell1->lgLx();
      float ux1 = cell1->lgLx() + cell1->width();
      for(AbacusCell* cell2 : cells_)
      {
        if(cell2 == cell1)
          break;
        float lx2 = cell2->lgLx();
        float ux2 = cell2->lgLx() + cell2->width();

        if(std::max(lx1, lx2) < std::min(ux1, ux2))
        {
          LOG_ERROR("Overlap found in abacus legalizer: {}: [{}, {}] -- {}: [{}, {}]", 
                    cell1->instance()->name(), cell1->lgLx(), cell1->lgLx() + cell1->width(),
                    cell2->instance()->name(), cell2->lgLx(), cell2->lgLx() + cell2->width());
        }
      }
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
    LOG_INFO("hpwl Before AbacusLegalization: {}", hpwlBeforeLG);
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
        double cbest = std::numeric_limits<double>::infinity();
        AbacusRow* rbest = nullptr;
        for(AbacusRow& row : rowStor_)
        {
          double clowerBound = cell.gpLy() - row.ly();
          clowerBound *= clowerBound;
          if(clowerBound > cbest)
            continue;
          if(row.usedWidth() + cell.width() > row.width())
            continue;

          // row.pushCell(&cell);
          // row.placeRow();
          row.tryAddCell(&cell);
          double dx = cell.gpLx() - cell.lgLx();
          double dy = cell.gpLy() - cell.lgLy();
          double c = dx * dx + dy * dy;
          if(c < cbest)
          {
            cbest = c;
            rbest = &row;
          }
          // assert(rbest != nullptr);
          // row.popCell();
        }
        if(rbest == nullptr)
        {
          LOG_ERROR("Lack of area. Unable to do stdcell legalization on die `{}`", die->name());
          break;
        }
        // rbest->pushCell(&cell);
        // rbest->placeRow();
        rbest->addCell(&cell);
      }
      LOG_TRACE("finish row assignment");

      // for(AbacusRow& row : rowStor_)
      // {
      //   row.findOverlap();
      // }

      LOG_TRACE("replace instance's location with its legalized location");
      for(AbacusCell& cell : cellStor_)
      {
        int lx = static_cast<int>(cell.lgLx() + 0.5f);
        int ly = static_cast<int>(cell.lgLy() + 0.5f);
        cell.instance()->setLocation(lx, ly);
      }
    }

    int64_t hpwlAfterLG = pb_->hpwl();
    LOG_INFO("hpwl After AbacusLegalization: {}", hpwlAfterLG);
    Plot::plot(pb_.get(), "./plot/cell", "after_lg");
  }

  void AbacusLegalizer::generateCells()
  {
    // Macros and fixed instances should not be treated as AbacusCell
    cellStor_.reserve(die_->placeInsts().size());
    for(Instance* inst : die_->placeInsts())
    {
      if(inst->isMacro())
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

    rowStor_.reserve(die_->rowRepeatCount());
    for(int i = 0; i < die_->rowRepeatCount(); i++)
    {
      // rowStor_.emplace_back(static_cast<float>(die_->rowStartX()),
      //                       static_cast<float>(die_->rowStartY() + i * die_->rowHeight()),
      //                       static_cast<float>(die_->rowWidth()),
      //                       static_cast<float>(die_->rowHeight()));
      // rows_.push_back(&rowStor_.back());

      int rowLx = die_->rowStartX();
      int rowLy = die_->rowStartY() + i * die_->rowHeight();
      int rowUx = die_->rowStartX() + die_->rowWidth();
      int rowUy = die_->rowStartY() + (i + 1) * die_->rowHeight();

      int currLx = rowLx;
      for(Instance* obs : obstacles)
      {
        // find overlap
        int lx = std::max(currLx, obs->lx());
        int ly = std::max(rowLy, obs->ly());
        int ux = std::min(rowUx, obs->ux());
        int uy = std::min(rowUy, obs->uy());

        // overlap exists
        if(lx < ux && ly < uy)
        {
          // |        |        |******|     |
          // ^        ^        ^      ^     ^
          // rowLx    currLx   lx     ux    rowUx
          if(currLx < lx)
          {
            // create subrow
            rowStor_.emplace_back(static_cast<float>(currLx),
                                  static_cast<float>(rowLy),
                                  static_cast<float>(lx - currLx),
                                  static_cast<float>(rowUy - rowLy));
          }
          // update currLx
          currLx = ux;
        }
      }
      
      if(currLx < rowUx)
      {
        rowStor_.emplace_back(static_cast<float>(currLx),
                              static_cast<float>(rowLy),
                              static_cast<float>(rowUx - currLx),
                              static_cast<float>(rowUy - rowLy));
      }
    }
  }

  void AbacusLegalizer::reset(Die* die)
  {
    rowStor_.clear();
    cellStor_.clear();
    die_ = die;

    generateCells();
    generateRows();
  }
}