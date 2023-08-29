#ifndef __REPLACE_ABACUS_LEGALIZER__
#define __REPLACE_ABACUS_LEGALIZER__

#include <memory>
#include <vector>
#include <cassert>

namespace replace
{
  class Instance;
  class Die;
  class PlacerBase;

  class AbacusLegalizerVars
  {
  public:
    enum { One, Area, NumPins } weightOpt;

    AbacusLegalizerVars();
  };

  class AbacusCell
  {
  public:
    AbacusCell() = default;
    AbacusCell(Instance* gcell);
    ~AbacusCell() = default;

    Instance *instance() const { return inst_; }

    float gpLx() const { return gpLx_; }
    float gpLy() const { return gpLy_; }
    float width() const { return width_; }
    float height() const { return height_; }

    float lgLx() const { return lgLx_; }
    void setLgLx(float x) { lgLx_ = x; } 
    float lgLy() const { return lgLy_; }
    void setLgLy(float y) { lgLy_ = y; }

    float weight() const { return weight_; }
    void setWeight(float w) { weight_ = w; }

  private:
    Instance* inst_;

    float gpLx_; // lx after global placement
    float lgLx_; // lx after legalization
    float width_;

    float lgLy_;
    float gpLy_;
    float height_;

    float weight_;
  };

  class AbacusCluster
  {
  public:
    AbacusCluster();
    AbacusCluster(int startIdx);
    ~AbacusCluster() = default;

    void addCell(AbacusCell* cell);
    void addCluster(AbacusCluster* other);
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

  class AbacusSubrow
  {
  public:
    AbacusSubrow();
    AbacusSubrow(float lx, float ux);
    ~AbacusSubrow() = default;

    void tryAddCell(AbacusCell* cell);
    void addCell(AbacusCell* cell);

    float lx() const { return lx_; }
    float ux() const { return ux_; }
    float width() const { return ux_ - lx_; }
    float usedWidth() const { return usedWidth_; }

    std::vector<AbacusCell*>& cells() { return cells_; }
    const std::vector<AbacusCell*>& cells() const { return cells_; }

  private:
    void appendCell(std::vector<AbacusCluster>& clusters, AbacusCell* cell);
    void collapse(std::vector<AbacusCluster>& clusters);

  private:
    std::vector<AbacusCluster> clusters_;
    std::vector<AbacusCell*> cells_;

    float lx_, ux_;
    float usedWidth_;
  };

  class AbacusRow
  {
  public:
    AbacusRow();
    AbacusRow(float lx, float ly, float ux, float uy);
    ~AbacusRow() = default;

    float lx() const { return lx_; }
    float ly() const { return ly_; }
    float ux() const { return ux_; }
    float uy() const { return uy_; }
    float width() const { return ux_ - lx_; }
    float height() const { return uy_ - ly_; }

    void generateSubrows(const std::vector<Instance*>& obstacles);
    std::vector<AbacusSubrow>& subrows() { return subrows_; }
    const std::vector<AbacusSubrow>& subrows() const { return subrows_; }

  private:
    float lx_, ly_, ux_, uy_;
    std::vector<AbacusSubrow> subrows_;
  };

  class AbacusLegalizer
  {
  public:
    AbacusLegalizer() = default;
    AbacusLegalizer(AbacusLegalizerVars lgVars, std::shared_ptr<PlacerBase> pb);
    ~AbacusLegalizer() = default;

    void doLegalization();

  private:
    void generateCells();
    void generateRows();

    void reset(Die* die);

  private:
    std::shared_ptr<PlacerBase> pb_;
    AbacusLegalizerVars lgVars_;

    Die* die_;

    std::vector<AbacusRow> rows_;
    std::vector<AbacusCell> cellStor_;
  };

}

#endif