#include "placerBase.h"
#include "log.h"

#include <iostream>

namespace replace
{

  using namespace std;

  static int fastModulo(const int input, const int ceil);

  static std::pair<int, int> getMinMaxIdx(int ll, int uu, int coreLL, int siteSize, int minIdx, int maxIdx);

  static bool isCoreAreaOverlap(Die &die, Instance &inst);

  static int64_t getOverlapWithCoreArea(Die &die, Instance &inst);

  ////////////////////////////////////////////////////////
  // Instance

  Instance::Instance()
      : lx_(0), ly_(0), ux_(0), uy_(0), extId_(INT_MIN)
  {
  }

  Instance::~Instance()
  {
    lx_ = ly_ = 0;
    ux_ = uy_ = 0;
    pins_.clear();
  }

  void Instance::setLocation(int x, int y)
  {
    setBox(x, y, x + (ux_ - lx_), y + (uy_ - ly_));
  }

  void Instance::setCenterLocation(int x, int y)
  {
    int halfX = (ux_ - lx_) / 2;
    int halfY = (uy_ - ly_) / 2;
    setBox(x - halfX, y - halfY, x + halfX, y + halfY);
  }

  void Instance::setSize(int w, int h)
  {
    setBox(lx_, ly_, lx_ + w, ly_ + h);
  }

  void Instance::setBox(int lx, int ly, int ux, int uy)
  {
    lx_ = lx;
    ly_ = ly;
    ux_ = ux;
    uy_ = uy;

    for (auto &pin : pins_)
    {
      pin->updateLocation(this);
    }
  }

  void Instance::addPin(Pin *pin)
  {
    pins_.push_back(pin);
  }

  ////////////////////////////////////////////////////////
  // Pin

  Pin::Pin()
      : inst_(nullptr), net_(nullptr),
        cx_(0), cy_(0),
        offsetCx_(0), offsetCy_(0),
        minPinX_(false), minPinY_(false),
        maxPinX_(false), maxPinY_(false)
  {
  }

  void Pin::updateLocation(const Instance *inst)
  {
    cx_ = inst->cx() + offsetCx_;
    cy_ = inst->cy() + offsetCy_;
  }

  void Pin::updateLocation(const Instance *inst, int offsetX, int offsetY)
  {
    offsetCx_ = offsetX;
    offsetCy_ = offsetY;
    cx_ = inst->cx() + offsetCx_;
    cy_ = inst->cy() + offsetCy_;
  }

  void Pin::setInstance(Instance *inst)
  {
    inst_ = inst;
  }

  void Pin::setNet(Net *net)
  {
    net_ = net;
  }

  Pin::~Pin()
  {
    inst_ = nullptr;
    net_ = nullptr;
  }

  ////////////////////////////////////////////////////////
  // Net

  Net::Net() : lx_(0), ly_(0), ux_(0), uy_(0)
  {
  }

  Net::~Net()
  {
    lx_ = ly_ = ux_ = uy_ = 0;
  }

  void Net::updateBox()
  {
    lx_ = INT_MAX;
    ly_ = INT_MAX;
    ux_ = INT_MIN;
    uy_ = INT_MIN;

    for (const Pin *p : pins_)
    {
      lx_ = std::min(p->cx(), lx_);
      ux_ = std::max(p->cx(), ux_);
      ly_ = std::min(p->cy(), ly_);
      uy_ = std::max(p->cy(), uy_);
    }
  }

  void Net::addPin(Pin *pin)
  {
    pins_.push_back(pin);
  }

  ////////////////////////////////////////////////////////
  // Die

  Die::Die() : dieLx_(0), dieLy_(0), dieUx_(0), dieUy_(0),
               coreLx_(0), coreLy_(0), coreUx_(0), coreUy_(0)
  {
  }

  Die::~Die()
  {
    dieLx_ = dieLy_ = dieUx_ = dieUy_ = 0;
    coreLx_ = coreLy_ = coreUx_ = coreUy_ = 0;
  }

  void Die::setDieBox(int lx, int ly, int ux, int uy)
  {
    dieLx_ = lx;
    dieLy_ = ly;
    dieUx_ = ux;
    dieUy_ = uy;
  }

  void Die::setRowParams(int startX, int startY, int width, int height, int repeatCount)
  {
    rowStartX_ = startX;
    rowStartY_ = startY;
    rowWidth_ = width;
    rowHeight_ = height;
    rowRepeatCount_ = repeatCount;

    coreLx_ = rowStartX_;
    coreLy_ = rowStartY_;
    coreUx_ = rowStartX_ + rowWidth_;
    coreUy_ = rowStartY_ + rowHeight_ * rowRepeatCount_;
  }

  ////////////////////////////////////////////////////////
  // PlacerBase

  PlacerBase::PlacerBase()
  // : placeInstsArea_(0), nonPlaceInstsArea_(0),
  //   macroInstsArea_(0), stdInstsArea_(0)
  {
  }

  PlacerBase::~PlacerBase()
  {
    reset();
  }

  void PlacerBase::reset()
  {
    instStor_.clear();
    pinStor_.clear();
    netStor_.clear();

    pins_.clear();
    nets_.clear();
    insts_.clear();

    placeInsts_.clear();
    fixedInsts_.clear();
    // dummyInsts_.clear();
    // nonPlaceInsts_.clear();
  }

  int64_t
  PlacerBase::hpwl() const
  {
    int64_t hpwl = 0;
    for (auto &net : nets_)
    {
      net->updateBox();
      hpwl += net->hpwl();
    }
    return hpwl;
  }

  void PlacerBase::printInfo() const
  {
    LOG_INFO("NumInstances: {}", instStor_.size());
    LOG_INFO("NumPlaceInstances: {}", placeInsts_.size());
    LOG_INFO("NumFixedInstances: {}", fixedInsts_.size());
    LOG_INFO("NumNets: {}", nets_.size());
    LOG_INFO("NumPins: {}", pins_.size());

    int maxFanout = INT_MIN;
    int sumFanout = 0;
    for (const Net *net : nets_)
    {
      maxFanout = std::max((int)net->pins().size(), maxFanout);
      sumFanout += (int)net->pins().size();
    }

    LOG_INFO("MaxFanout: {}", maxFanout);
    LOG_INFO("AvgFanout: {}", sumFanout / (float)nets_.size());

    LOG_INFO("DieAreaLxLy: ({}, {})", die_.dieLx(), die_.dieLy());
    LOG_INFO("DieAreaUxUy: ({}, {})", die_.dieUx(), die_.dieUy());
    LOG_INFO("CoreAreaLxLy: ({}, {})", die_.coreLx(), die_.coreLy());
    LOG_INFO("CoreAreaUxUy: ({}, {})", die_.coreUx(), die_.coreUy());

    int64_t coreArea =
        static_cast<int64_t>(die_.coreUx() - die_.coreLx()) *
        static_cast<int64_t>(die_.coreUy() - die_.coreLy());

    LOG_INFO("CoreArea: {}", coreArea);
    LOG_INFO("PlaceInstArea: {}", placeStdcellsArea_ + placeMacrosArea_);
    LOG_INFO("PlaceStdcellsArea: {}", placeStdcellsArea_);
    LOG_INFO("PlaceMacrosArea: {}", placeMacrosArea_);
    LOG_INFO("FixedInstArea: {}", fixedStdcellsArea_ + fixedMacrosArea_);
    LOG_INFO("FixedStdcellsArea: {}", fixedStdcellsArea_);
    LOG_INFO("FixedMacrosArea: {}", fixedMacrosArea_);
  }
}
