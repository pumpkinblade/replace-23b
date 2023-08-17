#ifndef __PLACER_BASE__
#define __PLACER_BASE__

#include <vector>
#include <unordered_map>
#include <map>
#include <memory>
#include <cassert>
#include "technology.h"

namespace replace
{
  class Pin;
  class Net;

  // Instance is either a cell or a macro related to a given technology node,
  // Instance can be fixed or not, has shape and location.
  // It also contains a pin list
  class Instance
  {
  public:
    Instance();
    ~Instance() = default;

    // a cell that no need to be moved.
    bool isFixed() const { return isFixed_; }
    void setFixed(bool on) { isFixed_ = on; }

    bool isMacro() const { return isMacro_; }
    void setMacro(bool on) { isMacro_ = on; }

    void setLocation(int x, int y);
    void setCenterLocation(int x, int y);
    // w for x dimension, h for y dimension
    void setSize(int w, int h);
    // use libCellName_ to get LibCell data under tech. and set
    // instance size to the LibCell size
    void setSize(Technology &tech);
    void setBox(int lx, int ly, int ux, int uy);

    int lx() const { return lx_; }
    int ly() const { return ly_; }
    int ux() const { return ux_; }
    int uy() const { return uy_; }
    int cx() const { return (lx_ + ux_) / 2; }
    int cy() const { return (ly_ + uy_) / 2; }
    // instance size in x dimension
    int dx() const { return ux_ - lx_; }
    // instance size in y dimension
    int dy() const { return uy_ - ly_; }
    int size() const { return dx() * dy(); }

    void setExtId(int extId) { extId_ = extId; }
    int extId() const { return extId_; }

    void addPin(Pin *pin);
    const std::vector<Pin *> &pins() const { return pins_; }

    const std::string &name() const { return name_; }
    void setName(const std::string &name) { name_ = name; }
    int libCellId() const { return libCellId_; }
    void setLibCellId(int id) { libCellId_ = id; }

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

    // For 23b
    std::string name_;
    int libCellId_;
  };

  // A Pin belongs to one instance, can be connected to no more than one net,
  // and has location info
  class Pin
  {
  public:
    Pin();
    ~Pin() = default;

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

    void setOffset(int offsetX, int offsetY);
    void setLocation(int cx, int cy);

    void updateLocation(const Instance *inst);
    // this is the only medifiter of offsetX nad offsetY
    void updateLocation(const Instance *inst, int offsetX, int offsetY);

    void setInstance(Instance *inst);
    void setNet(Net *net);

    Instance *instance() const { return inst_; }
    Net *net() const { return net_; }

    int libPinId() const { return libPinId_; }
    void setLibPinId(int id) { libPinId_ = id; }

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

    // For 23b
    int libPinId_;
  };

  // Net can have name, a vector of pins,
  // and hold a boundary to calculate HPWL
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

    // NOTE: this should be called sometime.
    void updateBox();

    const std::vector<Pin *> &pins() const { return pins_; }

    void addPin(Pin *pin);
    void removePin(Pin *pin);

  private:
    std::vector<Pin *> pins_;
    int lx_;
    int ly_;
    int ux_;
    int uy_;
  };

  // Die class containes shape info, row params, instances placed on it
  // and area statistics
  class Die
  {
  public:
    Die();
    ~Die() = default;
    void copyDieBoxAndRowParamsFrom(const Die &ano);

    void setDieBox(int lx, int ly, int ux, int uy);
    void setCoreBox(int lx, int ly, int ux, int uy);
    void setDieBox(const Die &ano);

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
    void setRowParams(const Die &ano);
    int rowStartX() const { return rowStartX_; }
    int rowStartY() const { return rowStartY_; }
    int rowWidth() const { return rowWidth_; }
    int rowHeight() const { return rowHeight_; }
    int rowRepeatCount() const { return rowRepeatCount_; }

    const std::vector<Instance *> &insts() const { return insts_; }
    std::vector<Instance *> &insts() { return insts_; }
    void addInstance(Instance *inst);
    void removeInstance(Instance *inst);

    const std::string &name() const { return name_; }
    void setName(const std::string &name) { name_ = name; }
    Technology *tech() const { return tech_; }
    void setTech(Technology *tech) { tech_ = tech; }
    float maxUtil() const { return maxUtil_; }
    void setMaxUtil(float util) { maxUtil_ = util; }

    bool isSetRow() { return isSetRow_; }

  private:
    int dieLx_;
    int dieLy_;
    int dieUx_;
    int dieUy_;

    int coreLx_;
    int coreLy_;
    int coreUx_;
    int coreUy_;

    bool isSetRow_;
    int rowStartX_;
    int rowStartY_;
    int rowWidth_;
    int rowHeight_;
    int rowRepeatCount_;

    std::vector<Instance *> insts_;

    // For 23b
    std::string name_;
    Technology *tech_;
    float maxUtil_;
  };

  // PlacerBase is a class containing all data, including instance set,
  // pin set, net set, die set. Also hold statistics about area and hpwl.
  // It lead out data by providing vector<T *> so T is writeable by
  // other classes. This class itself does not perform any placement.
  class PlacerBase
  {
    friend class Parser;
    friend class Partitioner;

  public:
    PlacerBase();
    ~PlacerBase() = default;

    const std::vector<Instance*>& insts() const { return insts_; }
    const std::vector<Die*>& dies() const { return dies_; }
    const std::vector<Technology*>& techs() const { return techs_; }
    const std::vector<Net*>& nets() const { return nets_; }
    const std::vector<Pin*>& pins() const { return pins_; }

    int64_t hpwl();
    void printDebugInfo() const;

    // query object by name
    Die *die(const std::string &name) const;
    Technology *tech(const std::string &name) const;

    // terminal params
    int terminalSizeX() const { return termSizeX_; }
    int terminalSizeY() const { return termSizeY_; }
    int terminalSpacing() const { return termSpace_; }
    int terminalCost() const { return termCost_; }

    // some modifier methods
    // create a new instance in instStor_, return reference
    Instance &emplaceInstance(bool isFixed, bool isMacro);
    Pin &emplacePin();
    // NOTE: this will move net into netStor_. Address will change
    void addNet(const Net &net);

  private:
    // clean vector of pointers
    void cleanIPNs();
    // derive insts_, pins_ and nets_ from XXXStor_. Assuming vector is empty
    // before call this method
    void deriveIPNs();
  private:
    std::vector<Die> dieStor_;
    std::vector<Instance> instStor_;
    std::vector<Pin> pinStor_;
    std::vector<Net> netStor_;
    std::vector<Technology> techStor_;

    // For partition
    std::vector<Instance> extraInstStor_;
    std::vector<Pin> extraPinStor_;
    std::vector<Net> extraNetStor_;

    // cheat: Save memory. After partition, it will be clear.
    std::vector<std::string> netNameStor_;

    std::vector<Die *> dies_;
    std::vector<Instance *> insts_;
    std::vector<Pin *> pins_;
    std::vector<Net *> nets_;
    std::vector<Technology*> techs_;

    std::unordered_map<std::string, Die *> dieNameMap_;
    std::unordered_map<std::string, Technology *> techNameMap_;

    int termSizeX_;
    int termSizeY_;
    int termSpace_;
    int termCost_;
  };
}

#endif