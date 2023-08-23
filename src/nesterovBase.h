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
    GCell(double cx, double cy, double dx, double dy);
    ~GCell() = default;

    Instance *instance() const { return inst_; }
    const std::vector<GPin *> &gPins() const { return gPins_; }

    void setInstance(Instance *inst);
    void addGPin(GPin *gPin);
    void updatePins(bool useTheta);

    void setMacro(bool on) { isMacro_ = on; }
    bool isInstance() const { return inst_ != nullptr; }
    bool isMacro() const { return isMacro_; }
    bool isFiller() const { return inst_ == nullptr; }

    // normal coordinates
    double lx() const { return cx_ - 0.5f * dx_; }
    double ly() const { return cy_ - 0.5f * dy_; }
    double ux() const { return cx_ + 0.5f * dx_; }
    double uy() const { return cy_ + 0.5f * dy_; }
    double cx() const { return cx_; }
    double cy() const { return cy_; }
    double dx() const { return dx_; }
    double dy() const { return dy_; }
    double theta() const { return theta_; }

    void setCenterLocation(double cx, double cy)
    {
      cx_ = cx;
      cy_ = cy;
    }
    void setSize(double dx, double dy)
    {
      dx_ = dx;
      dy_ = dy;
    }
    void setTheta(double theta) { theta_ = theta; }

    double densityScale() const { return densityScale_; }
    void setDensityScale(double densityScale) { densityScale_ = densityScale; }

    BinGrid *binGrid() const { return bg_; }
    void setBinGrid(BinGrid *bg) { bg_ = bg; }

  private:
    Instance *inst_;
    std::vector<GPin *> gPins_;
    // need to be stored for MS replace
    bool isMacro_;

    double cx_;
    double cy_;
    double dx_;
    double dy_;
    double theta_;

    double densityScale_;
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

    double lx() const { return lx_; }
    double ly() const { return ly_; }
    double ux() const { return ux_; }
    double uy() const { return uy_; }

    void addGPin(GPin *gPin);
    void updateBox();
    double hpwl() const { return (ux_ - lx_) + (uy_ - ly_); }

    // clear WA(Weighted Average) variables.
    void clearWaVars();

    void addWaExpMinSumX(double waExpMinX) { waExpMinSumX_ += waExpMinX; }
    void addWaXExpMinSumX(double waXExpMinX) { waXExpMinSumX_ += waXExpMinX; }

    void addWaExpMinSumY(double waExpMinY) { waExpMinSumY_ += waExpMinY; }
    void addWaYExpMinSumY(double waYExpMinY) { waYExpMinSumY_ += waYExpMinY; }

    void addWaExpMaxSumX(double waExpMaxX) { waExpMaxSumX_ += waExpMaxX; }
    void addWaXExpMaxSumX(double waXExpMaxX) { waXExpMaxSumX_ += waXExpMaxX; }

    void addWaExpMaxSumY(double waExpMaxY) { waExpMaxSumY_ += waExpMaxY; }
    void addWaYExpMaxSumY(double waYExpMaxY) { waYExpMaxSumY_ += waYExpMaxY; }

    double waExpMinSumX() const { return waExpMinSumX_; }
    double waXExpMinSumX() const { return waXExpMinSumX_; }

    double waExpMinSumY() const { return waExpMinSumY_; }
    double waYExpMinSumY() const { return waYExpMinSumY_; }

    double waExpMaxSumX() const { return waExpMaxSumX_; }
    double waXExpMaxSumX() const { return waXExpMaxSumX_; }

    double waExpMaxSumY() const { return waExpMaxSumY_; }
    double waYExpMaxSumY() const { return waYExpMaxSumY_; }

  private:
    Net *net_;
    std::vector<GPin *> gPins_;
    double lx_;
    double ly_;
    double ux_;
    double uy_;

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
    double waExpMinSumX_;
    double waXExpMinSumX_;

    double waExpMaxSumX_;
    double waXExpMaxSumX_;

    //
    // Y forces.
    //
    // waExpMinSumY_: store sigma {exp(y_i/gamma)}
    // waYExpMinSumY_: store signa {y_i*exp(e_i/gamma)}
    // waExpMaxSumY_ : store sigma {exp(-y_i/gamma)}
    // waYExpMaxSumY_: store sigma {y_i*exp(-y_i/gamma)}
    //
    double waExpMinSumY_;
    double waYExpMinSumY_;

    double waExpMaxSumY_;
    double waYExpMaxSumY_;
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

    double cx() const { return cx_; }
    double cy() const { return cy_; }
    double offsetCx() const { return offsetCx_; }
    double offsetCy() const { return offsetCy_; }

    // clear WA(Weighted Average) variables.
    void clearWaVars();

    void setMaxExpSumX(double maxExpSumX)
    {
      hasMaxExpSumX_ = true;
      maxExpSumX_ = maxExpSumX;
    }
    void setMaxExpSumY(double maxExpSumY)
    {
      hasMaxExpSumY_ = true;
      maxExpSumY_ = maxExpSumY;
    }
    void setMinExpSumX(double minExpSumX)
    {
      hasMinExpSumX_ = true;
      minExpSumX_ = minExpSumX;
    }
    void setMinExpSumY(double minExpSumY)
    {
      hasMinExpSumY_ = true;
      minExpSumY_ = minExpSumY;
    }

    double maxExpSumX() const { return maxExpSumX_; }
    double maxExpSumY() const { return maxExpSumY_; }
    double minExpSumX() const { return minExpSumX_; }
    double minExpSumY() const { return minExpSumY_; }

    bool hasMaxExpSumX() const { return hasMaxExpSumX_; }
    bool hasMaxExpSumY() const { return hasMaxExpSumY_; }
    bool hasMinExpSumX() const { return hasMinExpSumX_; }
    bool hasMinExpSumY() const { return hasMinExpSumY_; }

    void setCenterLocation(double cx, double cy);
    void updateLocation(const GCell *gCell);
    void updateLocationWithTheta(const GCell *gCell);

  private:
    Pin *pin_;
    GCell *gCell_;
    GNet *gNet_;

    double offsetCx_;
    double offsetCy_;
    double cx_;
    double cy_;

    // weighted average WL vals stor for better indexing
    // Please check the equation (4) in the ePlace-MS paper.
    //
    // maxExpSum_: holds exp(x_i/gamma)
    // minExpSum_: holds exp(-x_i/gamma)
    // the x_i is equal to cx_ variable.
    //
    double maxExpSumX_;
    double maxExpSumY_;

    double minExpSumX_;
    double minExpSumY_;

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
    Bin(double lx, double ly, double ux, double uy);
    ~Bin() = default;

    double lx() const { return lx_; }
    double ly() const { return ly_; }
    double ux() const { return ux_; }
    double uy() const { return uy_; }
    double cx() const { return (lx_ + ux_) / 2; }
    double cy() const { return (ly_ + uy_) / 2; }
    double dx() const { return (ux_ - lx_); }
    double dy() const { return (uy_ - ly_); }

    double electroPhi() const { return electroPhi_; }
    double electroForceX() const { return electroForceX_; }
    double electroForceY() const { return electroForceY_; }
    double density() const { return density_; }

    void setDensity(double density) { density_ = density; }
    void setElectroForceX(double force) { electroForceX_ = force; }
    void setElectroForceY(double force) { electroForceY_ = force; }
    void setElectroPhi(double phi) { electroPhi_ = phi; }

    void setNonPlaceArea(double area) { nonPlaceArea_ = area; }
    void setInstPlacedArea(double area) { instPlacedArea_ = area; }
    void setFillerArea(double area) { fillerArea_ = area; }

    void addNonPlaceArea(double area) { nonPlaceArea_ += area; }
    void addInstPlacedArea(double area) { instPlacedArea_ += area; }
    void addFillerArea(double area) { fillerArea_ += area; }

    double binArea() const { return (ux_ - lx_) * (uy_ - ly_); }
    double nonPlaceArea() const { return nonPlaceArea_; }
    double instPlacedArea() const { return instPlacedArea_; }
    double fillerArea() const { return fillerArea_; }

  private:
    // coordinate
    double lx_;
    double ly_;
    double ux_;
    double uy_;

    double nonPlaceArea_;
    double instPlacedArea_;
    double fillerArea_;
    double density_;

    double electroPhi_;
    double electroForceX_;
    double electroForceY_;
  };

  class BinGrid
  {
  public:
    BinGrid();
    BinGrid(Die *die);
    ~BinGrid() = default;

    void setDie(Die *die);
    void setTargetDensity(double density) { targetDensity_ = density; }

    void initBins();
    void initBins(int cntX, int cntY);

    void updateBinsNonPlaceArea();
    void updateBinsGCellDensityArea(bool useTheta);
    void addFillerBinArea(const GCell *gcell);
    void addStdCellBinArea(const GCell *gcell);
    void addMacroBinArea(const GCell *gcell);
    void addMacroBinAreaWithTheta(const GCell *gcell);

    // lx, ly, ux, uy will hold coreArea
    double lx() const { return lx_; }
    double ly() const { return ly_; }
    double ux() const { return ux_; }
    double uy() const { return uy_; }
    double cx() const { return (lx_ + ux_) / 2; }
    double cy() const { return (ly_ + uy_) / 2; }
    double dx() const { return (ux_ - lx_); }
    double dy() const { return (uy_ - ly_); }

    int binCntX() const { return binCntX_; }
    int binCntY() const { return binCntY_; }
    double binSizeX() const { return binSizeX_; }
    double binSizeY() const { return binSizeY_; }

    double overflowArea() const { return overflowArea_; }

    // return overlap bins_ index
    std::pair<int, int> getMinMaxIdxX(double lx1, double ux1);
    std::pair<int, int> getMinMaxIdxY(double ly1, double uy1);

    const std::vector<Bin *> &bins() const { return bins_; }
    Die *die() const { return die_; }
    const std::vector<GCell *> gCells() const { return gCells_; }
    double sumPhi() const { return sumPhi_; }

    void addGCell(GCell *gCell);
    void updateDensityForceBin();
    void updateGCellDensityScaleAndSize();

    int placeInstanceCount() const { return placeInstCnt_; }
    double placeStdCellArea() const { return placeStdCellArea_; }
    double placeMacroArea() const { return placeMacroArea_; }
    double fixedInstanceArea() const { return fixedInstArea_; }
    double totalCellArea() const { return totalCellArea_; }

  private:
    Die *die_;
    std::vector<Bin> binStor_;
    std::vector<Bin *> bins_;
    std::vector<GCell *> gCells_;
    FFT fft_;

    double lx_;
    double ly_;
    double ux_;
    double uy_;

    int binCntX_;
    int binCntY_;
    double binSizeX_;
    double binSizeY_;

    double targetDensity_;
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
    double targetDensity;
    double minAvgCut;
    double maxAvgCut;
    int binCntX;
    int binCntY;
    double minWireLengthForceBar;
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

    void init();
    void initFillerGCells(const Die *die);

    const std::vector<GCell *> &gCells() const { return gCells_; }
    const std::vector<GNet *> &gNets() const { return gNets_; }
    const std::vector<GPin *> &gPins() const { return gPins_; }
    const std::vector<BinGrid *> &binGrids() const { return binGrids_; }
    const std::shared_ptr<PlacerBase> &pb() const { return pb_; }

    // update gCells with cx, cy
    void updateGCellCenterLocation(const std::vector<Point> &coordis);
    void updateGCellCenterLocationWithTheta(const std::vector<Point> &coordis, const std::vector<double> &macrosTheta);

    double overflowArea() const;
    double sumPhi() const;
    double targetDensity() const;

    void updateDensityCoordiLayoutInside(GCell *gcell, bool useTheta);
    double getDensityCoordiLayoutInsideX(const GCell *gcell, double newCx, bool useTheta);
    double getDensityCoordiLayoutInsideY(const GCell *gcell, double newCy, bool useTheta);

    // WL force update based on WeightedAverage model
    // wlCoeffX : WireLengthCoefficient for X.
    //            equal to 1 / gamma_x
    // wlCoeffY : WireLengthCoefficient for Y.
    //            equal to 1 / gamma_y
    //
    // Gamma is described in the ePlaceMS paper.
    //
    void updateWireLengthForceWA(double wlCoeffX, double wlCoeffY);
    Point getWireLengthGradientPinWA(const GPin *gPin, double wlCoeffX, double wlCoeffY);

    // gradient of wire length
    Point getWireLengthGradientWA(const GCell *gCell, double wlCoeffX, double wlCoeffY);
    Point getWireLengthGradientWAWithTheta(const GCell *gCell, double wlCoeffX, double wlCoeffY, double &gradTheta);

    // gradient of density
    Point getDensityGradient(const GCell *gCell);
    Point getDensityGradientWithTheta(const GCell *gCell, double &gradTheta);

    // gradient of local density
    Point getDensityGradientLocal(const GCell *gCell, double alpha, double beta, double &cellDelta);
    Point getDensityGradientLocalWithTheta(const GCell *gCell, double alpha, double beta, double &cellDelta, double &gradTheta);

    // for preconditioner
    Point getWireLengthPreconditioner(const GCell *gCell);
    Point getDensityPreconditioner(const GCell *gCell);
    double getWireLengthPreconditionerTheta(const GCell *gCell);
    double getDensityPreconditionerTheta(const GCell *gCell);

    // update electrostatic forces within Bin
    void updateDensityForceBin();

    void updateNetsBox();
    double getHpwl() const;
    double getOverflow() const;

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

    double sumPhi_;
  };

}

#endif
