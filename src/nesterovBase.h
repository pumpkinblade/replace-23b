#ifndef __NESTEROV_BASE__
#define __NESTEROV_BASE__

#include <vector>
#include <memory>
#include <unordered_map>

#include "fft.h"
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
    GCell(prec cx, prec cy, prec dx, prec dy);
    ~GCell() = default;

    Instance *instance() const { return inst_; }
    const std::vector<GPin *> &gPins() const { return gPins_; }

    void addGPin(GPin *gPin);
    void setInstance(Instance *inst);
    void setMacro(bool on) { isMacro_ = on; }

    bool isInstance() const { return inst_ != nullptr; }
    bool isFiller() const { return inst_ == nullptr; }
    bool isMacro() const { return isMacro_; }

    // normal coordinates
    prec lx() const { return lx_; }
    prec ly() const { return ly_; }
    prec ux() const { return ux_; }
    prec uy() const { return uy_; }
    prec cx() const { return (lx_ + ux_) / 2; }
    prec cy() const { return (ly_ + uy_) / 2; }
    prec dx() const { return (ux_ - lx_); }
    prec dy() const { return (uy_ - ly_); }
    prec theta() const { return theta_; }

    void setLocation(prec lx, prec ly);
    void setCenterLocation(prec cx, prec cy);
    void setCenterLocationTheta(prec cx, prec cy, prec theta);
    void setSize(prec dx, prec dy);
    void setThetaNoUpdatePin(prec theta) { theta_ = theta; }

    prec densityScale() const { return densityScale_; }
    void setDensityScale(prec densityScale) { densityScale_ = densityScale; }

    BinGrid* binGrid() const { return bg_; }
    void setBinGrid(BinGrid* bg) { bg_ = bg; }

  private:
    Instance *inst_;
    std::vector<GPin *> gPins_;
    // need to be stored for MS replace
    bool isMacro_;

    prec lx_;
    prec ly_;
    prec ux_;
    prec uy_;
    prec theta_;

    prec densityScale_;
    BinGrid* bg_;
  };

  class GNet
  {
  public:
    GNet();
    GNet(Net *net);
    ~GNet() = default;

    Net *net() const { return net_; }
    const std::vector<GPin *> &gPins() const { return gPins_; }

    prec lx() const { return lx_; }
    prec ly() const { return ly_; }
    prec ux() const { return ux_; }
    prec uy() const { return uy_; }

    void addGPin(GPin *gPin);
    void updateBox();
    prec hpwl() const { return (ux_ - lx_) + (uy_ - ly_); }

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
    prec lx_;
    prec ly_;
    prec ux_;
    prec uy_;

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

    prec cx() const { return cx_; }
    prec cy() const { return cy_; }
    prec offsetCx() const { return offsetCx_; }
    prec offsetCy() const { return offsetCy_; }

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

    void setCenterLocation(prec cx, prec cy);
    void updateLocation(const GCell *gCell);
    void updateLocationWithTheta(const GCell *gCell);

  private:
    Pin* pin_;
    GCell *gCell_;
    GNet *gNet_;

    prec offsetCx_;
    prec offsetCy_;
    prec cx_;
    prec cy_;

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
    Bin(int x, int y, prec lx, prec ly, prec ux, prec uy, prec targetDensity);
    ~Bin() = default;

    int x() const { return x_; }
    int y() const { return y_; }

    prec lx() const { return lx_; }
    prec ly() const { return ly_; }
    prec ux() const { return ux_; }
    prec uy() const { return uy_; }
    prec cx() const { return (lx_ + ux_) / 2; }
    prec cy() const { return (ly_ + uy_) / 2; }
    prec dx() const { return (ux_ - lx_); }
    prec dy() const { return (uy_ - ly_); }

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

    void setNonPlaceArea(prec area) { nonPlaceArea_ = area; }
    void setInstPlacedArea(prec area) { instPlacedArea_ = area; }
    void setFillerArea(prec area) { fillerArea_ = area; }

    void addNonPlaceArea(prec area) { nonPlaceArea_ += area; }
    void addInstPlacedArea(prec area) { instPlacedArea_ += area; }
    void addFillerArea(prec area) { fillerArea_ += area; }

    prec binArea() const { return (ux_ - lx_) * (uy_ - ly_); }
    prec nonPlaceArea() const { return nonPlaceArea_; }
    prec instPlacedArea() const { return instPlacedArea_; }
    prec fillerArea() const { return fillerArea_; }

  private:
    // index
    int x_;
    int y_;

    // coordinate
    prec lx_;
    prec ly_;
    prec ux_;
    prec uy_;

    prec nonPlaceArea_;
    prec instPlacedArea_;
    prec fillerArea_;

    prec density_;
    prec targetDensity_; // will enable bin-wise density screening
    prec electroPhi_;
    prec electroForceX_;
    prec electroForceY_;
  };

  class BinGrid
  {
  public:
    BinGrid();
    BinGrid(Die *die);
    ~BinGrid() = default;

    void setDie(Die *die);
    void setBinCntX(int binCntX);
    void setBinCntY(int binCntY);
    void setTargetDensity(prec density);

    void updateBinsGCellDensityArea();
    void addFillerBinArea(const GCell* gcell);
    void addStdCellBinArea(const GCell* gcell);
    void addMacroBinArea(const GCell* gcell);
    void addMacroBinAreaWithTheta(const GCell* gcell);

    void initBins();

    // lx, ly, ux, uy will hold coreArea
    prec lx() const { return lx_; }
    prec ly() const { return ly_; }
    prec ux() const { return ux_; }
    prec uy() const { return uy_; }
    prec cx() const { return (lx_ + ux_) / 2; }
    prec cy() const { return (ly_ + uy_) / 2; }
    prec dx() const { return (ux_ - lx_); }
    prec dy() const { return (uy_ - ly_); }

    int binCntX() const { return binCntX_; }
    int binCntY() const { return binCntY_; }
    prec binSizeX() const { return binSizeX_; }
    prec binSizeY() const { return binSizeY_; }

    double overflowArea() const { return overflowArea_; }

    // return overlap bins_ index
    std::pair<int, int> getMinMaxIdxX(prec lx1, prec ux1);
    std::pair<int, int> getMinMaxIdxY(prec ly1, prec uy1);

    const std::vector<Bin *> &bins() const { return bins_; }
    Die *die() const { return die_; }
    const std::vector<GCell *> gCells() const { return gCells_; }
    double sumPhi() const { return sumPhi_; }

    void addGCell(GCell *gc);
    void updateDensityForceBin();
    void updateGCellDensityScaleAndSize();

    int placeInstanceCount() const { return placeInstCnt_; }
    double placeStdCellArea() const { return placeStdCellArea_; }
    double placeMacroArea() const { return placeMacroArea_; }
    double fixedInstanceArea() const { return fixedInstArea_; }
    double totalCellArea() const { return totalCellArea_; }

  private:
    void updateBinsNonPlaceArea();

  private:
    Die *die_;
    std::vector<Bin> binStor_;
    std::vector<Bin *> bins_;
    std::vector<GCell *> gCells_;
    FFT fft_;

    prec lx_;
    prec ly_;
    prec ux_;
    prec uy_;

    int binCntX_;
    int binCntY_;
    bool isSetBinCntX_;
    bool isSetBinCntY_;
    prec binSizeX_;
    prec binSizeY_;

    prec targetDensity_;
    double sumPhi_;
    double overflowArea_;

    int placeInstCnt_;
    double placeStdCellArea_;
    double placeMacroArea_;
    double fixedInstArea_;
    double totalCellArea_;
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
    bool isSetBinCntX;
    bool isSetBinCntY;

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
    const std::vector<GNet *> &gNets() const { return gNets_; }
    const std::vector<GPin *> &gPins() const { return gPins_; }
    const std::vector<BinGrid *> &binGrids() const { return binGrids_; }
    const std::shared_ptr<PlacerBase> &pb() const { return pb_; }

    // update gCells with lx, ly
    void updateGCellLocation(std::vector<Point> &points);
    // update gCells with cx, cy
    void updateGCellCenterLocation(std::vector<Point> &points);

    double overflowArea() const;
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
    void updateWireLengthForceWA(prec wlCoeffX, prec wlCoeffY);

    Point getWireLengthGradientPinWA(GPin *gPin, prec wlCoeffX, prec wlCoeffY);

    Point getWireLengthGradientWA(GCell *gCell, prec wlCoeffX, prec wlCoeffY);
    Point getWireLengthGradientWAWithTheta(GCell *gCell, prec wlCoeffX, prec wlCoeffY, prec& gradTheta);

    // for preconditioner
    Point getWireLengthPreconditioner(GCell *gCell);

    Point getDensityPreconditioner(GCell *gCell);

    Point getDensityGradient(GCell *gCell);
    Point getDensityGradientWithTheta(GCell *gCell, prec& gradTheta);

    Point getDensityGradientLocal(GCell *gCell, prec alpha, prec beta, prec& cellDelta);
    Point getDensityGradientLocalWithTheta(GCell *gCell, prec alpha, prec beta, prec& cellDelta, prec& gradTheta);

    double hpwl();
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
    std::vector<GCell> gCellStor_;
    std::vector<GNet> gNetStor_;
    std::vector<GPin> gPinStor_;

    std::vector<BinGrid *> binGrids_;
    std::vector<GCell *> gCells_;
    std::vector<GNet *> gNets_;
    std::vector<GPin *> gPins_;

    prec sumPhi_;
  };

}

#endif
