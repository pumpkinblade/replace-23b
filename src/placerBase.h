#ifndef __PLACER_BASE__
#define __PLACER_BASE__

#include <vector>
#include <unordered_map>
#include <memory>

namespace replace
{
  class Pin;
  class Net;

  class Instance
  {
  public:
    Instance();
    ~Instance();

    // a cell that no need to be moved.
    bool isFixed() const { return isFixed_; }
    void setFixed(bool on) { isFixed_ = on; }

    bool isMacro() const { return isMacro_; }
    void setMacro(bool on) { isMacro_ = on; }

    void setLocation(int x, int y);
    void setCenterLocation(int x, int y);
    void setSize(int w, int h);
    void setBox(int lx, int ly, int ux, int uy);

    int lx() const { return lx_; }
    int ly() const { return ly_; }
    int ux() const { return ux_; }
    int uy() const { return uy_; }
    int cx() const { return (lx_ + ux_) / 2; }
    int cy() const { return (ly_ + uy_) / 2; }
    int dx() const { return ux_ - lx_; }
    int dy() const { return uy_ - ly_; }

    void setExtId(int extId) { extId_ = extId; }
    int extId() const { return extId_; }

    void addPin(Pin *pin);
    const std::vector<Pin *> &pins() const { return pins_; }

  private:
    std::vector<Pin *> pins_;
    int lx_;
    int ly_;
    int ux_;
    int uy_;
    bool isFixed_;
    bool isMacro_;

    // For initialPlace
    int extId_;
  };

  class Pin
  {
  public:
    Pin();
    ~Pin();

    bool isMinPinX() const { return minPinX_; }
    bool isMaxPinX() const { return maxPinX_; }
    bool isMinPinY() const { return minPinY_; }
    bool isMaxPinY() const { return maxPinY_; }

    void setMinPinX(bool on) { minPinX_ = on; }
    void setMaxPinX(bool on) { maxPinX_ = on; }
    void setMinPinY(bool on) { minPinY_ = on; }
    void setMaxPinY(bool on) { maxPinY_ = on; }

    int cx() const { return cx_; }
    int cy() const { return cy_; }

    int offsetCx() const { return offsetCx_; }
    int offsetCy() const { return offsetCy_; }

    void updateLocation(const Instance *inst);
    void updateLocation(const Instance *inst, int offsetX, int offsetY);

    void setInstance(Instance *inst);
    void setNet(Net *net);

    Instance *instance() const { return inst_; }
    Net *net() const { return net_; }

  private:
    Instance *inst_;
    Net *net_;

    // pin center coordinate is enough
    // Pins' placed location.
    int cx_;
    int cy_;

    // offset coordinates inside instance.
    // origin point is center point of instance.
    // (e.g. (DX/2,DY/2) )
    // This will increase efficiency for bloating
    int offsetCx_;
    int offsetCy_;

    // For initialPlace
    bool minPinX_;
    bool minPinY_;
    bool maxPinX_;
    bool maxPinY_;
  };

  class Net
  {
  public:
    Net();
    ~Net();

    int lx() const { return lx_; }
    int ly() const { return ly_; }
    int ux() const { return ux_; }
    int uy() const { return uy_; }
    int cx() const { return (lx_ + ux_) / 2; }
    int cy() const { return (ly_ + uy_) / 2; }

    // HPWL: half-parameter-wire-length
    int64_t hpwl() const { return static_cast<int64_t>((ux_ - lx_) + (uy_ - ly_)); }

    void updateBox();

    const std::vector<Pin *> &pins() const { return pins_; }

    void addPin(Pin *pin);

  private:
    std::vector<Pin *> pins_;
    int lx_;
    int ly_;
    int ux_;
    int uy_;
  };

  class Die
  {
  public:
    Die();
    ~Die() = default;

    void setDieBox(int lx, int ly, int ux, int uy);

    int dieLx() const { return dieLx_; }
    int dieLy() const { return dieLy_; }
    int dieUx() const { return dieUx_; }
    int dieUy() const { return dieUy_; }

    int coreLx() const { return coreLx_; }
    int coreLy() const { return coreLy_; }
    int coreUx() const { return coreUx_; }
    int coreUy() const { return coreUy_; }

    int dieCx() const { return (dieLx_ + dieUx_) / 2; }
    int dieCy() const { return (dieLy_ + dieUy_) / 2; }
    int dieDx() const { return dieUx_ - dieLy_; }
    int dieDy() const { return dieUy_ - dieLy_; }
    int coreCx() const { return (coreLx_ + coreUx_) / 2; }
    int coreCy() const { return (coreLy_ + coreUy_) / 2; }
    int coreDx() const { return coreUx_ - coreLx_; }
    int coreDy() const { return coreUy_ - coreLy_; }

    void setRowParams(int startX, int startY, int width, int height, int repeatCount);
    int rowStartX() const { return rowStartX_; }
    int rowStartY() const { return rowStartY_; }
    int rowWidth() const { return rowWidth_; }
    int rowHeight() const { return rowHeight_; }
    int rowRepeatCount() const { return rowRepeatCount_; }

    const std::vector<Instance*>& insts() const { return insts_; }
    const std::vector<Instance*>& placeInsts() const { return placeInsts_; }
    const std::vector<Instance*>& fixedInsts() const { return fixedInsts_; }
    void addInstance(Instance* inst);

    int64_t placeInstsArea() const { return placeStdcellsArea_ + placeMacrosArea_; }
    int64_t placeStdcellsArea() const { return placeStdcellsArea_; }
    int64_t placeMacrosArea() const { return placeMacrosArea_; }
    int64_t fixedInstsArea() const { return fixedInstsArea_; }

  private:
    int dieLx_;
    int dieLy_;
    int dieUx_;
    int dieUy_;

    int coreLx_;
    int coreLy_;
    int coreUx_;
    int coreUy_;

    int rowStartX_;
    int rowStartY_;
    int rowWidth_;
    int rowHeight_;
    int rowRepeatCount_;

    std::vector<Instance*> insts_;
    std::vector<Instance*> placeInsts_;
    std::vector<Instance*> fixedInsts_;

    // placeInstsArea_ = placeStdcellArea_ + placeMacroArea_
    int64_t placeStdcellsArea_;
    int64_t placeMacrosArea_;
    // fixedInstsArea_ = fixedStdcellArea_ + fixedMacroArea_
    int64_t fixedInstsArea_;
  };

  class PlacerBase
  {
    friend class Parser;

  public:
    PlacerBase();
    ~PlacerBase();

    const std::vector<Instance *> &insts() const { return insts_; }
    const std::vector<Pin *> &pins() const { return pins_; }
    const std::vector<Net *> &nets() const { return nets_; }

    //
    // placeInsts : a real instance that need to be placed
    // fixedInsts : a real instance that is fixed (e.g. macros, tapcells)
    //
    const std::vector<Instance *> &placeInsts() const { return placeInsts_; }
    const std::vector<Instance *> &fixedInsts() const { return fixedInsts_; }

    // Note: each die should have the same size
    Die* die() { return dies_.front(); }
    const std::vector<Die*>& dies() const { return dies_; }

    int64_t hpwl() const;
    void printInfo() const;

    int64_t placeInstsArea() const;
    int64_t placeStdcellsArea() const;
    int64_t placeMacrosArea() const;
    int64_t fixedInstsArea() const;

  private:
    std::vector<Die> dieStor_;
    std::vector<Instance> instStor_;
    std::vector<Pin> pinStor_;
    std::vector<Net> netStor_;

    std::vector<Die *> dies_;
    std::vector<Instance *> insts_;
    std::vector<Pin *> pins_;
    std::vector<Net *> nets_;

    std::vector<Instance *> placeInsts_;
    std::vector<Instance *> fixedInsts_;

    void reset();
  };

}

#endif
