#ifndef __REPLACE_ABAX_LEGALIZER__
#define __REPLACE_ABAX_LEGALIZER__

#include <memory>
#include <vector>
#include <cassert>
#include "placerBase.h"

namespace replace
{
  class AbaxCell
  {
  public:
    AbaxCell() : inst_(nullptr),  lgLx_(0.f) {}
    AbaxCell(Instance* inst) : inst_(inst), lgLx_(0.f) {}
    ~AbaxCell() = default;

    Instance *instance() const { return inst_; }

    float gpLx() const { return static_cast<float>(inst_->lx()); }
    float gpLy() const { return static_cast<float>(inst_->ly()); }
    float width() const { return static_cast<float>(inst_->dx()); }
    float height() const { return static_cast<float>(inst_->dy()); }

    float lgLx() const { return lgLx_; }
    void setLgLx(float x) { lgLx_ = x; } 

  private:
    Instance* inst_;
    float lgLx_;
  };

  class AbaxCluster
  {
  public:
    AbaxCluster();
    AbaxCluster(int startIdx);
    ~AbaxCluster() = default;

    void addCell(const AbaxCell& cell);
    void addCluster(AbaxCluster& other);
    void reset();

    int startIdx() const { return startIdx_; }
    int endIdx() const { return endIdx_; }
    float xc() const { return xc_; }
    float wc() const { return wc_; }
    float qc() const { return qc_; }
    float ec() const { return ec_; }

    void setXc(float x) { xc_ = x; }

  private:
    int startIdx_;
    int endIdx_;
    float xc_;
    float ec_;
    float qc_;
    float wc_;
  };

  class AbaxSubrow
  {
  public:
    AbaxSubrow();
    AbaxSubrow(float lx, float ux);
    ~AbaxSubrow() = default;

    void reset();
    void addCell(int cellIdx, std::vector<AbaxCell>& cells);

    float lx() const { return lx_; }
    float ux() const { return ux_; }
    float width() const { return ux_ - lx_; }
    float usedWidth() const { return usedWidth_; }

    int cellStartIdx() const { return cellStartIdx_; }
    int cellEndIdx() const { return cellEndIdx_; }
    bool isEmpty() const { return cellStartIdx_ >= cellEndIdx_; }

  private:
    void appendCell(std::vector<AbaxCluster>& clusters, int cellIdx, std::vector<AbaxCell>& cells);
    void collapse(std::vector<AbaxCluster>& clusters);

  private:
    float lx_, ux_;
    std::vector<AbaxCluster> clusters_;
    int cellStartIdx_, cellEndIdx_;
    float usedWidth_;
  };

  class AbaxRow
  {
  public:
    AbaxRow();
    AbaxRow(float lx, float ly, float ux, float uy);
    ~AbaxRow() = default;

    float lx() const { return lx_; }
    float ly() const { return ly_; }
    float ux() const { return ux_; }
    float uy() const { return uy_; }
    float width() const { return ux_ - lx_; }
    float height() const { return uy_ - ly_; }

    void addInstance(Instance* inst, int subrowIdx);
    void tryAddInstance(Instance* inst, int& subrowIdx, float& cost);
    void generateSubrows(const std::vector<Instance*>& obstacles);
    std::vector<AbaxSubrow>& subrows() { return subrows_; }
    const std::vector<AbaxSubrow>& subrows() const { return subrows_; }

    float getCost(const std::vector<AbaxCell> cells) const;
    void applyLegalization();

  private:
    bool recursiveMoves(int groupStartIdx, int groupEndIdx, std::vector<AbaxCell>& cells,
                        int subrowIdx, std::vector<AbaxSubrow>& subrows);

  private:
    float lx_, ly_, ux_, uy_;

    std::vector<AbaxSubrow> subrows_;
    std::vector<AbaxCell> cells_;

    std::vector<AbaxSubrow> trySubrows_;
    std::vector<AbaxCell> tryCells_;
  };

  class AbaxLegalizer
  {
  public:
    AbaxLegalizer() = default;
    AbaxLegalizer(std::shared_ptr<PlacerBase> pb);
    ~AbaxLegalizer() = default;

    void doLegalization();

  private:
    void generateRows();

    void reset(Die* die);

  private:
    std::shared_ptr<PlacerBase> pb_;

    Die* die_;

    std::vector<AbaxRow> rows_;
  };

}

#endif