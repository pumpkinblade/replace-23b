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

    // filler cells
    GCell(int cx, int cy, int dx, int dy);
    ~GCell() = default;

    Instance *instance() const { return inst_; }
    const std::vector<GPin *> &gPins() const { return gPins_; }

    void addGPin(GPin *gPin);

    void setInstance(Instance *inst);
    void setMacro(bool on) { isMacro_ = on; }

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

    void setDensityScale(prec densityScale) { densityScale_ = densityScale; }
    void setGradientX(prec gradX) { gradientX_ = gradX; }
    void setGradientY(prec gradY) { gradientY_ = gradY; }

    prec gradientX() const { return gradientX_; }
    prec gradientY() const { return gradientY_; }
    prec densityScale() const { return densityScale_; }

    bool isInstance() const { return inst_ != nullptr; }
    bool isFiller() const { return inst_ == nullptr; }
    bool isMacro() const { return isMacro_; }

    void setBinGrid(BinGrid *bg) { bg_ = bg; }
    BinGrid *binGrid() const { return bg_; }

  private:
    Instance *inst_;
    std::vector<GPin *> gPins_;
    int lx_;
    int ly_;
    int ux_;
    int uy_;

    int dLx_;
    int dLy_;
    int dUx_;
    int dUy_;

    prec densityScale_;
    prec gradientX_;
    prec gradientY_;

    // need to be stored for
    // MS replace
    bool isMacro_;

    BinGrid *bg_;
  };

  class GNet
  {
  public:
    GNet();
    GNet(Net *net);
    ~GNet() = default;

    Net *net() const { return net_; }
    const std::vector<GPin *> &gPins() const { return gPins_; }

    int lx() const { return lx_; }
    int ly() const { return ly_; }
    int ux() const { return ux_; }
    int uy() const { return uy_; }

    void addGPin(GPin *gPin);
    void updateBox();
    int64_t hpwl() const;

    // clear WA(Weighted Average) variables.
    void clearWaVars();

    void addWaExpMinSumX(prec waExpMinX) { waExpMinSumX_ += waExpMinX; }
    void addWaXExpMinSumX(prec waXExpMinX) { waXExpMinSumX_ += waXExpMinX; }

    void addWaExpMinSumY(prec waExpMinY) { waExpMinSumY_ += waExpMinY; }
    void addWaYExpMinSumY(prec waYExpMinY) { waYExpMinSumY_ += waYExpMinY; }

    void addWaExpMaxSumX(prec waExpMaxX) { waExpMaxSumX_ += waExpMaxX; }
    void addWaXExpMaxSumX(prec waXExpMaxX) { waXExpMaxSumX_ += waXExpMaxX; }

    void addWaExpMaxSumY(prec waExpMaxY) { waExpMaxSumY_ += waExpMaxY; }
    void addWaYExpMaxSumY(prec waYExpMaxY) { waYExpMaxSumY_ += waYExpMaxY; }

    prec waExpMinSumX() const { return waExpMinSumX_; }
    prec waXExpMinSumX() const { return waXExpMinSumX_; }

    prec waExpMinSumY() const { return waExpMinSumY_; }
    prec waYExpMinSumY() const { return waYExpMinSumY_; }

    prec waExpMaxSumX() const { return waExpMaxSumX_; }
    prec waXExpMaxSumX() const { return waXExpMaxSumX_; }

    prec waExpMaxSumY() const { return waExpMaxSumY_; }
    prec waYExpMaxSumY() const { return waYExpMaxSumY_; }

  private:
    Net* net_;
    std::vector<GPin *> gPins_;
    int lx_;
    int ly_;
    int ux_;
    int uy_;

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
    prec waExpMinSumX_;
    prec waXExpMinSumX_;

    prec waExpMaxSumX_;
    prec waXExpMaxSumX_;

    //
    // Y forces.
    //
    // waExpMinSumY_: store sigma {exp(y_i/gamma)}
    // waYExpMinSumY_: store signa {y_i*exp(e_i/gamma)}
    // waExpMaxSumY_ : store sigma {exp(-y_i/gamma)}
    // waYExpMaxSumY_: store sigma {y_i*exp(-y_i/gamma)}
    //
    prec waExpMinSumY_;
    prec waYExpMinSumY_;

    prec waExpMaxSumY_;
    prec waYExpMaxSumY_;
  };

  class GPin
  {
  public:
    GPin();
    GPin(Pin *pin);
    ~GPin() = default;

    Pin *pin() const { return pin_; }

    GCell *gCell() const { return gCell_; }
    GNet *gNet() const { return gNet_; }

    void setGCell(GCell *gCell) { gCell_ = gCell; }
    void setGNet(GNet *gNet) { gNet_ = gNet; }

    int cx() const { return cx_; }
    int cy() const { return cy_; }

    // clear WA(Weighted Average) variables.
    void clearWaVars();

    void setMaxExpSumX(prec maxExpSumX);
    void setMaxExpSumY(prec maxExpSumY);
    void setMinExpSumX(prec minExpSumX);
    void setMinExpSumY(prec minExpSumY);

    prec maxExpSumX() const { return maxExpSumX_; }
    prec maxExpSumY() const { return maxExpSumY_; }
    prec minExpSumX() const { return minExpSumX_; }
    prec minExpSumY() const { return minExpSumY_; }

    bool hasMaxExpSumX() const { return (hasMaxExpSumX_ == 1); }
    bool hasMaxExpSumY() const { return (hasMaxExpSumY_ == 1); }
    bool hasMinExpSumX() const { return (hasMinExpSumX_ == 1); }
    bool hasMinExpSumY() const { return (hasMinExpSumY_ == 1); }

    void setCenterLocation(int cx, int cy);
    void updateLocation(const GCell *gCell);
    void updateDensityLocation(const GCell *gCell);

  private:
    Pin* pin_;
    GCell *gCell_;
    GNet *gNet_;

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
    prec maxExpSumX_;
    prec maxExpSumY_;

    prec minExpSumX_;
    prec minExpSumY_;

    // flag variables
    //
    // check whether
    // this pin is considered in a WA models.
    bool hasMaxExpSumX_;
    bool hasMaxExpSumY_;
    bool hasMinExpSumX_;
    bool hasMinExpSumY_;
  };

  class Bin
  {
  public:
    Bin();
    Bin(int x, int y, int lx, int ly, int ux, int uy, prec targetDensity);
    ~Bin() = default;

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

    prec electroPhi() const { return electroPhi_; }
    prec electroForceX() const { return electroForceX_; }
    prec electroForceY() const { return electroForceY_; }
    prec targetDensity() const { return targetDensity_; }
    prec density() const { return density_; }

    void setDensity(prec density) { density_ = density; }
    void setTargetDensity(prec density) { targetDensity_ = density; }
    void setElectroForceX(prec force) { electroForceX_ = force; }
    void setElectroForceY(prec force) { electroForceY_ = force; }
    void setElectroPhi(prec phi) { electroPhi_ = phi; }

    void setNonPlaceArea(int64_t area) { nonPlaceArea_ = area; }
    void setInstPlacedArea(int64_t area) { instPlacedArea_ = area; }
    void setFillerArea(int64_t area) { fillerArea_ = area; }

    void addNonPlaceArea(int64_t area) { nonPlaceArea_ += area; }
    void addInstPlacedArea(int64_t area) { instPlacedArea_ += area; }
    void addFillerArea(int64_t area) { fillerArea_ += area; }

    int64_t binArea() const { return (int64_t)(ux_ - lx_) * (uy_ - ly_); }
    int64_t nonPlaceArea() const { return nonPlaceArea_; }
    int64_t instPlacedArea() const { return instPlacedArea_; }
    int64_t fillerArea() const { return fillerArea_; }

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

    prec density_;
    prec targetDensity_; // will enable bin-wise density screening
    prec electroPhi_;
    prec electroForceX_;
    prec electroForceY_;
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
    BinGrid(BinGrid &&other);
    BinGrid &operator=(BinGrid &&other);
    BinGrid(const BinGrid &other) = delete;
    BinGrid &operator=(const BinGrid &other) = delete;
    ~BinGrid() = default;

    void setDie(Die *die);
    void setBinCntX(int binCntX);
    void setBinCntY(int binCntY);
    void setTargetDensity(prec density);
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
    prec sumPhi() const { return sumPhi_; }

    void addGCell(GCell *gc);
    void updateDensityForceBin();
    void updateGCellDensityScaleAndSize();

    int placeInstanceCount() const { return placeInstCnt_; }
    int64_t placeStdCellArea() const { return placeStdCellArea_; }
    int64_t placeMacroArea() const { return placeMacroArea_; }
    int64_t fixedInstanceArea() const { return fixedInstArea_; }

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
    prec targetDensity_;
    int64_t overflowArea_;
    bool isSetBinCntX_;
    bool isSetBinCntY_;
    prec sumPhi_;

    Die *die_;
    std::vector<GCell *> gCells_;
    std::unique_ptr<FFT> fft_;
    int placeInstCnt_;
    int64_t placeStdCellArea_;
    int64_t placeMacroArea_;
    int64_t fixedInstArea_;

    void updateBinsNonPlaceArea();
  };

  class NesterovBaseVars
  {
  public:
    prec targetDensity;
    prec minAvgCut;
    prec maxAvgCut;
    int binCntX;
    int binCntY;
    prec minWireLengthForceBar;
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
    ~NesterovBase() = default;

    const std::vector<GCell *> &gCells() const { return gCells_; }
    // const std::vector<GCell *> &gCellInsts() const { return gCellInsts_; }
    // const std::vector<GCell *> &gCellFillers() const { return gCellFillers_; }

    const std::vector<GNet *> &gNets() const { return gNets_; }
    const std::vector<GPin *> &gPins() const { return gPins_; }
    const std::vector<BinGrid *> &binGrids() const { return binGrids_; }
    const std::shared_ptr<PlacerBase> &pb() const { return pb_; }

    ////
    //// placerBase To NesterovBase functions
    ////
    //GCell *placerToNesterov(Instance *inst);
    //GPin *placerToNesterov(Pin *pin);
    //GNet *placerToNesterov(Net *net);
    //BinGrid *placerToNesterov(Die *die);

    // update gCells with lx, ly
    void updateGCellLocation(
        std::vector<Point> &points);

    // update gCells with cx, cy
    void updateGCellCenterLocation(
        std::vector<Point> &points);

    void updateGCellDensityCenterLocation(
        std::vector<Point> &points);

    int64_t overflowArea() const;
    prec sumPhi() const;
    prec targetDensity() const;

    void updateDensityCoordiLayoutInside(GCell *gcell);

    prec getDensityCoordiLayoutInsideX(GCell *gCell, prec cx);
    prec getDensityCoordiLayoutInsideY(GCell *gCell, prec cy);

    // WL force update based on WeightedAverage model
    // wlCoeffX : WireLengthCoefficient for X.
    //            equal to 1 / gamma_x
    // wlCoeffY : WireLengthCoefficient for Y.
    //            equal to 1 / gamma_y
    //
    // Gamma is described in the ePlaceMS paper.
    //
    void updateWireLengthForceWA(
        prec wlCoeffX,
        prec wlCoeffY);

    Point getWireLengthGradientPinWA(
        GPin *gPin, prec wlCoeffX, prec wlCoeffY);

    Point getWireLengthGradientWA(
        GCell *gCell, prec wlCoeffX, prec wlCoeffY);

    // for preconditioner
    Point getWireLengthPreconditioner(GCell *gCell);

    Point getDensityPreconditioner(GCell *gCell);

    Point getDensityGradient(GCell *gCell);

    int64_t hpwl();
    prec overflow() const;

    // update electrostatic forces within Bin
    void updateDensityForceBin();
  
  private:
    void init();
    void initFillerGCells(BinGrid* bg);

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

    //std::unordered_map<Instance *, GCell *> gCellMap_;
    //std::unordered_map<Pin *, GPin *> gPinMap_;
    //std::unordered_map<Net *, GNet *> gNetMap_;
    //std::unordered_map<Die *, BinGrid *> binGridMap_;

    prec sumPhi_;
  };

}

#endif
