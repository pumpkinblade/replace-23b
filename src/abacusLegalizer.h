#ifndef __REPLACE_ABACUS_LEGALIZER__
#define __REPLACE_ABACUS_LEGALIZER__

#include <memory>
#include <vector>

namespace replace
{
  class Instance;
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
    ~AbacusCluster() = default;

    void addCell(AbacusCell* cell);
    void addCluster(AbacusCluster* cluster);
    void place();
    void reset();

    float xc() const { return xc_; }
    float wc() const { return wc_; }
    float qc() const { return qc_; }
    float ec() const { return ec_; }

    void setXc(float x) { xc_ = x; }
  private:
    std::vector<AbacusCell *> cells_;
    float xc_;
    float ec_;
    float qc_;
    float wc_;
  };

  class AbacusRow
  {
  public:
    AbacusRow();
    AbacusRow(float lx, float ly, float w, float h);
    ~AbacusRow() = default;

    void pushCell(AbacusCell* cell);
    void popCell();
    void placeRow();

    float lx() const { return lx_; }
    float ly() const { return ly_; }
    float width() const { return width_; }
    float height() const { return height_; }
    float usedWidth() const { return usedWidth_; }

    const std::vector<AbacusCell*>& cells() const { return cells_; }

  private:
    void collapse();

  private:
    std::vector<AbacusCluster> clusters_;
    std::vector<AbacusCell*> cells_;

    float lx_, ly_;
    float width_, height_;
    float usedWidth_;
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

  private:
    std::shared_ptr<PlacerBase> pb_;
    AbacusLegalizerVars lgVars_;

    std::vector<AbacusRow> rowStor_;
    std::vector<AbacusCell> cellStor_;

    std::vector<AbacusRow*> rows_;
    std::vector<AbacusCell*> cells_;
  };

}

#endif