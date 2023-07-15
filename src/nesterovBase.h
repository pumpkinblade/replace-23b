#ifndef __NESTEROV_BASE__
#define __NESTEROV_BASE__

#include <vector>
#include <memory>
#include <unordered_map>

#include "point.h"

namespace replace
{
  class Instance;
  class Die;
  class PlacerBase;

  class Instance;
  class Pin;
  class Net;

  class GPin;
  class FFT;
  class BinGrid;

  class GCell
  {
  public:
    GCell();

    // instance cells
    GCell(Instance *inst);
    GCell(std::vector<Instance *> &insts);

    // filler cells
    GCell(int cx, int cy, int dx, int dy);
    ~GCell();

    Instance *instance() const;
    const std::vector<Instance *> &insts() const { return insts_; }
    const std::vector<GPin *> &gPins() const { return gPins_; }

    void addGPin(GPin *gPin);

    void setClusteredInstance(std::vector<Instance *> &insts);
    void setInstance(Instance *inst);
    void setFiller();
    void setMacroInstance();
    void setStdInstance();

    // normal coordinates
    int lx() const { return lx_; }
    int ly() const { return ly_; }
    int ux() const { return ux_; }
    int uy() const { return uy_; }
    int cx() const { return (lx_ + ux_) / 2; }
    int cy() const { return (ly_ + uy_) / 2; }
    int dx() const { return (ux_ - lx_); }
    int dy() const { return (uy_ - ly_); }

    // virtual density coordinates
    int dLx() const { return dLx_; }
    int dLy() const { return dLy_; }
    int dUx() const { return dUx_; }
    int dUy() const { return dUy_; }
    int dCx() const { return (dLx_ + dUx_) / 2; }
    int dCy() const { return (dLy_ + dUy_) / 2; }
    int dDx() const { return (dUx_ - dLx_); }
    int dDy() const { return (dUy_ - dLy_); }

    void setLocation(int lx, int ly);
    void setCenterLocation(int cx, int cy);
    void setSize(int dx, int dy);

    void setDensityLocation(int dLx, int dLy);
    void setDensityCenterLocation(int dCx, int dCy);
    void setDensitySize(int dDx, int dDy);

    void setDensityScale(float densityScale);
    void setGradientX(float gradX);
    void setGradientY(float gradY);

    float gradientX() const { return gradientX_; }
    float gradientY() const { return gradientY_; }
    float densityScale() const { return densityScale_; }

    bool isInstance() const;
    bool isClusteredInstance() const;
    bool isFiller() const;
    bool isMacroInstance() const;
    bool isStdInstance() const;

    void setBinGrid(BinGrid *bg) { bg_ = bg; }
    BinGrid *binGrid() const { return bg_; }

  private:
    std::vector<Instance *> insts_;
    std::vector<GPin *> gPins_;
    int lx_;
    int ly_;
    int ux_;
    int uy_;

    int dLx_;
    int dLy_;
    int dUx_;
    int dUy_;

    float densityScale_;
    float gradientX_;
    float gradientY_;

    // need to be stored for
    // MS replace
    bool isMacroInstance_;

    BinGrid *bg_;
  };

  class GNet
  {
  public:
    GNet();
    GNet(Net *net);
    GNet(std::vector<Net *> &nets);
    ~GNet();

    Net *net() const;
    const std::vector<Net *> &nets() const { return nets_; }
    const std::vector<GPin *> &gPins() const { return gPins_; }

    int lx() const { return lx_; }
    int ly() const { return ly_; }
    int ux() const { return ux_; }
    int uy() const { return uy_; }

    void setCustomWeight(float customWeight);
    float customWeight() const { return customWeight_; }
    float netWeight() const { return weight_; }

    void addGPin(GPin *gPin);
    void updateBox();
    int64_t hpwl();

    void setDontCare();
    bool isDontCare();

    // clear WA(Weighted Average) variables.
    void clearWaVars();

    void addWaExpMinSumX(float waExpMinX) { waExpMinSumX_ += waExpMinX; }
    void addWaXExpMinSumX(float waXExpMinX) { waXExpMinSumX_ += waXExpMinX; }

    void addWaExpMinSumY(float waExpMinY) { waExpMinSumY_ += waExpMinY; }
    void addWaYExpMinSumY(float waYExpMinY) { waYExpMinSumY_ += waYExpMinY; }

    void addWaExpMaxSumX(float waExpMaxX) { waExpMaxSumX_ += waExpMaxX; }
    void addWaXExpMaxSumX(float waXExpMaxX) { waXExpMaxSumX_ += waXExpMaxX; }

    void addWaExpMaxSumY(float waExpMaxY) { waExpMaxSumY_ += waExpMaxY; }
    void addWaYExpMaxSumY(float waYExpMaxY) { waYExpMaxSumY_ += waYExpMaxY; }

    float waExpMinSumX() const { return waExpMinSumX_; }
    float waXExpMinSumX() const { return waXExpMinSumX_; }

    float waExpMinSumY() const { return waExpMinSumY_; }
    float waYExpMinSumY() const { return waYExpMinSumY_; }

    float waExpMaxSumX() const { return waExpMaxSumX_; }
    float waXExpMaxSumX() const { return waXExpMaxSumX_; }

    float waExpMaxSumY() const { return waExpMaxSumY_; }
    float waYExpMaxSumY() const { return waYExpMaxSumY_; }

  private:
    std::vector<GPin *> gPins_;
    std::vector<Net *> nets_;
    int lx_;
    int ly_;
    int ux_;
    int uy_;

    float customWeight_;
    float weight_;

    //
    // weighted average WL model stor for better indexing
    // Please check the equation (4) in the ePlace-MS paper.
    //
    // WA: weighted Average
    // saving four variable will be helpful for
    // calculating the WA gradients/wirelengths.
    //
    // gamma: modeling accuracy.
    //
    // X forces.
    //
    // waExpMinSumX_: store sigma {exp(x_i/gamma)}
    // waXExpMinSumX_: store sigma {x_i*exp(x_i/gamma)}
    // waExpMaxSumX_ : store sigma {exp(-x_i/gamma)}
    // waXExpMaxSumX_: store sigma {x_i*exp(-x_i/gamma)}
    //
    float waExpMinSumX_;
    float waXExpMinSumX_;

    float waExpMaxSumX_;
    float waXExpMaxSumX_;

    //
    // Y forces.
    //
    // waExpMinSumY_: store sigma {exp(y_i/gamma)}
    // waYExpMinSumY_: store signa {y_i*exp(e_i/gamma)}
    // waExpMaxSumY_ : store sigma {exp(-y_i/gamma)}
    // waYExpMaxSumY_: store sigma {y_i*exp(-y_i/gamma)}
    //
    float waExpMinSumY_;
    float waYExpMinSumY_;

    float waExpMaxSumY_;
    float waYExpMaxSumY_;

    unsigned char isDontCare_ : 1;
  };

  class GPin
  {
  public:
    GPin();
    GPin(Pin *pin);
    GPin(std::vector<Pin *> &pins);
    ~GPin();

    Pin *pin() const;
    const std::vector<Pin *> &pins() const { return pins_; }

    GCell *gCell() const { return gCell_; }
    GNet *gNet() const { return gNet_; }

    void setGCell(GCell *gCell);
    void setGNet(GNet *gNet);

    int cx() const { return cx_; }
    int cy() const { return cy_; }

    // clear WA(Weighted Average) variables.
    void clearWaVars();

    void setMaxExpSumX(float maxExpSumX);
    void setMaxExpSumY(float maxExpSumY);
    void setMinExpSumX(float minExpSumX);
    void setMinExpSumY(float minExpSumY);

    float maxExpSumX() const { return maxExpSumX_; }
    float maxExpSumY() const { return maxExpSumY_; }
    float minExpSumX() const { return minExpSumX_; }
    float minExpSumY() const { return minExpSumY_; }

    bool hasMaxExpSumX() const { return (hasMaxExpSumX_ == 1); }
    bool hasMaxExpSumY() const { return (hasMaxExpSumY_ == 1); }
    bool hasMinExpSumX() const { return (hasMinExpSumX_ == 1); }
    bool hasMinExpSumY() const { return (hasMinExpSumY_ == 1); }

    void setCenterLocation(int cx, int cy);
    void updateLocation(const GCell *gCell);
    void updateDensityLocation(const GCell *gCell);

  private:
    GCell *gCell_;
    GNet *gNet_;
    std::vector<Pin *> pins_;

    int offsetCx_;
    int offsetCy_;
    int cx_;
    int cy_;

    // weighted average WL vals stor for better indexing
    // Please check the equation (4) in the ePlace-MS paper.
    //
    // maxExpSum_: holds exp(x_i/gamma)
    // minExpSum_: holds exp(-x_i/gamma)
    // the x_i is equal to cx_ variable.
    //
    float maxExpSumX_;
    float maxExpSumY_;

    float minExpSumX_;
    float minExpSumY_;

    // flag variables
    //
    // check whether
    // this pin is considered in a WA models.
    unsigned char hasMaxExpSumX_ : 1;
    unsigned char hasMaxExpSumY_ : 1;

    unsigned char hasMinExpSumX_ : 1;
    unsigned char hasMinExpSumY_ : 1;
  };

  class Bin
  {
  public:
    Bin();
    Bin(int x, int y, int lx, int ly, int ux, int uy, float targetDensity);

    ~Bin();

    int x() const { return x_; }
    int y() const { return y_; }

    int lx() const { return lx_; }
    int ly() const { return ly_; }
    int ux() const { return ux_; }
    int uy() const { return uy_; }
    int cx() const { return (lx_ + ux_) / 2; }
    int cy() const { return (ly_ + uy_) / 2; }
    int dx() const { return (ux_ - lx_); }
    int dy() const { return (uy_ - ly_); }

    float electroPhi() const;
    float electroForceX() const;
    float electroForceY() const;
    float targetDensity() const;
    float density() const;

    void setDensity(float density);
    void setTargetDensity(float density);
    void setElectroForce(float electroForceX, float electroForceY);
    void setElectroPhi(float phi);

    void setNonPlaceArea(int64_t area) { nonPlaceArea_ = area; }
    void setInstPlacedArea(int64_t area) { instPlacedArea_ = area; }
    void setFillerArea(int64_t area) { fillerArea_ = area; }

    void addNonPlaceArea(int64_t area) { nonPlaceArea_ += area; }
    void addInstPlacedArea(int64_t area) { instPlacedArea_ += area; }
    void addFillerArea(int64_t area) { fillerArea_ += area; }

    const int64_t binArea() const;
    const int64_t nonPlaceArea() const { return nonPlaceArea_; }
    const int64_t instPlacedArea() const { return instPlacedArea_; }
    const int64_t fillerArea() const { return fillerArea_; }

  private:
    // index
    int x_;
    int y_;

    // coordinate
    int lx_;
    int ly_;
    int ux_;
    int uy_;

    int64_t nonPlaceArea_;
    int64_t instPlacedArea_;
    int64_t fillerArea_;

    float density_;
    float targetDensity_; // will enable bin-wise density screening
    float electroPhi_;
    float electroForceX_;
    float electroForceY_;
  };

  //
  // The bin can be non-uniform because of
  // "integer" coordinates
  //
  class BinGrid
  {
  public:
    BinGrid();
    BinGrid(Die *die);
    BinGrid(BinGrid&& other);
    BinGrid& operator=(BinGrid&& other);
    BinGrid(const BinGrid& other) = delete;
    BinGrid& operator=(const BinGrid& other) = delete;
    ~BinGrid() = default;

    void setDie(Die *die);
    void setBinCntX(int binCntX);
    void setBinCntY(int binCntY);
    void setTargetDensity(float density);
    void updateBinsGCellDensityArea();

    void initBins();

    // lx, ly, ux, uy will hold coreArea
    int lx() const { return lx_; }
    int ly() const { return ly_; }
    int ux() const { return ux_; }
    int uy() const { return uy_; }
    int cx() const { return (lx_ + ux_) / 2; }
    int cy() const { return (ly_ + uy_) / 2; }
    int dx() const { return (ux_ - lx_); }
    int dy() const { return (uy_ - ly_); }

    int binCntX() const { return binCntX_; }
    int binCntY() const { return binCntY_; }
    int binSizeX() const { return binSizeX_; }
    int binSizeY() const { return binSizeY_; }

    int64_t overflowArea() const { return overflowArea_; }

    // return bins_ index with given gcell
    std::pair<int, int> getDensityMinMaxIdxX(GCell *gcell);
    std::pair<int, int> getDensityMinMaxIdxY(GCell *gcell);

    std::pair<int, int> getMinMaxIdxX(Instance *inst);
    std::pair<int, int> getMinMaxIdxY(Instance *inst);

    const std::vector<Bin *> &bins() const { return bins_; }
    Die *die() const { return die_; }
    const std::vector<GCell *> gCells() const { return gCells_; }
    float sumPhi() const { return sumPhi_; }

    void addGCell(GCell *gc);
    void updateDensityForceBin();
    void updateGCellDensityScaleAndSize();

  private:
    std::vector<Bin> binStor_;
    std::vector<Bin *> bins_;
    int lx_;
    int ly_;
    int ux_;
    int uy_;
    int binCntX_;
    int binCntY_;
    int binSizeX_;
    int binSizeY_;
    float targetDensity_;
    int64_t overflowArea_;
    bool isSetBinCntX_;
    bool isSetBinCntY_;
    float sumPhi_;

    Die *die_;
    std::vector<GCell *> gCells_;
    std::unique_ptr<FFT> fft_;

    void updateBinsNonPlaceArea();
  };

  class NesterovBaseVars
  {
  public:
    float targetDensity;
    float minAvgCut;
    float maxAvgCut;
    int binCntX;
    int binCntY;
    float minWireLengthForceBar;
    unsigned char isSetBinCntX : 1;
    unsigned char isSetBinCntY : 1;

    NesterovBaseVars();
    void reset();
  };

  class NesterovBase
  {
  public:
    NesterovBase();
    NesterovBase(NesterovBaseVars nbVars, std::shared_ptr<PlacerBase> pb);
    ~NesterovBase();

    const std::vector<GCell *> &gCells() const { return gCells_; }
    // const std::vector<GCell *> &gCellInsts() const { return gCellInsts_; }
    // const std::vector<GCell *> &gCellFillers() const { return gCellFillers_; }

    const std::vector<GNet *> &gNets() const { return gNets_; }
    const std::vector<GPin *> &gPins() const { return gPins_; }
    const std::vector<BinGrid *> &binGrids() const { return binGrids_; }

    //
    // placerBase To NesterovBase functions
    //
    GCell *placerToNesterov(Instance *inst);
    GPin *placerToNesterov(Pin *pin);
    GNet *placerToNesterov(Net *net);
    BinGrid *placerToNesterov(Die* die);

    // update gCells with lx, ly
    void updateGCellLocation(
        std::vector<FloatPoint> &points);

    // update gCells with cx, cy
    void updateGCellCenterLocation(
        std::vector<FloatPoint> &points);

    void updateGCellDensityCenterLocation(
        std::vector<FloatPoint> &points);

    int64_t overflowArea() const;
    float sumPhi() const;
    float targetDensity() const;

    void updateDensityCoordiLayoutInside(GCell *gcell);

    float getDensityCoordiLayoutInsideX(GCell *gCell, float cx);
    float getDensityCoordiLayoutInsideY(GCell *gCell, float cy);

    // WL force update based on WeightedAverage model
    // wlCoeffX : WireLengthCoefficient for X.
    //            equal to 1 / gamma_x
    // wlCoeffY : WireLengthCoefficient for Y.
    //            equal to 1 / gamma_y
    //
    // Gamma is described in the ePlaceMS paper.
    //
    void updateWireLengthForceWA(
        float wlCoeffX,
        float wlCoeffY);

    FloatPoint getWireLengthGradientPinWA(
        GPin *gPin, float wlCoeffX, float wlCoeffY);

    FloatPoint getWireLengthGradientWA(
        GCell *gCell, float wlCoeffX, float wlCoeffY);

    // for preconditioner
    FloatPoint getWireLengthPreconditioner(GCell *gCell);

    FloatPoint getDensityPreconditioner(GCell *gCell);

    FloatPoint getDensityGradient(GCell *gCell);

    int64_t getHpwl();

    // update electrostatic forces within Bin
    void updateDensityForceBin();

  private:
    NesterovBaseVars nbVars_;
    std::shared_ptr<PlacerBase> pb_;

    std::vector<BinGrid> binGridStor_;
    std::vector<BinGrid *> binGrids_;

    std::vector<GCell> gCellStor_;
    std::vector<GNet> gNetStor_;
    std::vector<GPin> gPinStor_;

    std::vector<GCell *> gCells_;
    // std::vector<GCell *> gCellInsts_;
    // std::vector<GCell *> gCellFillers_;

    std::vector<GNet *> gNets_;
    std::vector<GPin *> gPins_;

    std::unordered_map<Instance *, GCell *> gCellMap_;
    std::unordered_map<Pin *, GPin *> gPinMap_;
    std::unordered_map<Net *, GNet *> gNetMap_;
    std::unordered_map<Die *, BinGrid *> binGridMap_;

    float sumPhi_;

    void init();
    void initFillerGCells(Die* die);

    void reset();
  };

}

#endif
