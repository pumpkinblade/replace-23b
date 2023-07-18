#include "placerBase.h"
#include "log.h"

#include <iostream>

namespace replace
{

  using namespace std;

  ////////////////////////////////////////////////////////
  // Instance

  Instance::Instance()
      : lx_(0), ly_(0), ux_(0), uy_(0), extId_(INT_MIN)
  {
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
    pinNameMap_.emplace(pin->name(), pin);
  }

  Pin* Instance::pin(const std::string& name) const
  {
    auto it = pinNameMap_.find(name);
    return it == pinNameMap_.end() ? nullptr : it->second;
  }

  ////////////////////////////////////////////////////////
  // Pin

  Pin::Pin()
      : inst_(nullptr), net_(nullptr),
        cx_(0), cy_(0),
        offsetCx_(0), offsetCy_(0)
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

  Die::Die() 
      : dieLx_(0), dieLy_(0), dieUx_(0), dieUy_(0),
        coreLx_(0), coreLy_(0), coreUx_(0), coreUy_(0),
        fixedInstsArea_(0), placeStdcellsArea_(0), placeMacrosArea_(0)
  {
  }

  void Die::setDieBox(int lx, int ly, int ux, int uy)
  {
    dieLx_ = lx;
    dieLy_ = ly;
    dieUx_ = ux;
    dieUy_ = uy;
  }

  void Die::setCoreBox(int lx, int ly, int ux, int uy)
  {
    coreLx_ = lx;
    coreLy_ = ly;
    coreUx_ = ux;
    coreUy_ = uy;
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

  void Die::addInstance(Instance *inst)
  {
    insts_.push_back(inst);
    if (inst->isFixed())
    {
      fixedInsts_.push_back(inst);
    }
    else
    {
      placeInsts_.push_back(inst);
    }

    int64_t area = static_cast<int64_t>(inst->dx()) *
                   static_cast<int64_t>(inst->dy());
    if (inst->isFixed())
    {
      fixedInstsArea_ += area;
    }
    else if (inst->isMacro())
    {
      placeMacrosArea_ += area;
    }
    else
    {
      placeStdcellsArea_ += area;
    }
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
    for (Technology* tech : techs_)
    {
      LOG_INFO("Technology {}", tech->name());
      tech->printInfo();
    }

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

    Die *die = dies_.front();
    LOG_INFO("CoreAreaLxLy: ({}, {})", die->coreLx(), die->coreLy());
    LOG_INFO("CoreAreaUxUy: ({}, {})", die->coreUx(), die->coreUy());

    int64_t coreArea =
        static_cast<int64_t>(die->coreUx() - die->coreLx()) *
        static_cast<int64_t>(die->coreUy() - die->coreLy());

    LOG_INFO("CoreArea: {}", coreArea);
    int dieIdx = 0;
    for (const Die *die : dies_)
    {
      LOG_INFO("Die {}", dieIdx);
      LOG_INFO("PlaceInstArea: {}", die->placeInstsArea());
      LOG_INFO("PlaceStdcellsArea: {}", die->placeStdcellsArea());
      LOG_INFO("PlaceMacrosArea: {}", die->placeMacrosArea());
      LOG_INFO("FixedInstArea: {}", die->fixedInstsArea());

      int numStdCell = 0, numMacro = 0;
      for (Instance* inst : die->insts())
      {
        if (inst->isMacro())
          numMacro++;
        else
          numStdCell++;
      }
      LOG_INFO("NumStdCells: {}", numStdCell);
      LOG_INFO("NumMacros: {}", numMacro);

      dieIdx++;
    }
  }

  int64_t PlacerBase::placeInstsArea() const
  {
    int64_t area = 0;
    for (const Die *die : dies_)
    {
      area += die->placeInstsArea();
    }
    return area;
  }

  int64_t PlacerBase::placeStdcellsArea() const
  {
    int64_t area = 0;
    for (const Die *die : dies_)
    {
      area += die->placeStdcellsArea();
    }
    return area;
  }

  int64_t PlacerBase::placeMacrosArea() const
  {
    int64_t area = 0;
    for (const Die *die : dies_)
    {
      area += die->placeMacrosArea();
    }
    return area;
  }

  int64_t PlacerBase::fixedInstsArea() const
  {
    int64_t area = 0;
    for (const Die *die : dies_)
    {
      area += die->fixedInstsArea();
    }
    return area;
  }

  Instance *PlacerBase::inst(const std::string &name) const
  {
    auto it = instNameMap_.find(name);
    return it == instNameMap_.end() ? nullptr : it->second;
  }

  Die *PlacerBase::die(const std::string &name) const
  {
    auto it = dieNameMap_.find(name);
    return it == dieNameMap_.end() ? nullptr : it->second;
  }

  Net *PlacerBase::net(const std::string &name) const
  {
    auto it = netNameMap_.find(name);
    return it == netNameMap_.end() ? nullptr : it->second;
  }

  Technology *PlacerBase::tech(const std::string &name) const
  {
    auto it = techNameMap_.find(name);
    return it == techNameMap_.end() ? nullptr : it->second;
  }
}
