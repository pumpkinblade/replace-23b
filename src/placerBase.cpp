#include "placerBase.h"
#include "log.h"

#include <iostream>
#include <algorithm>

namespace replace
{

  using namespace std;

  ////////////////////////////////////////////////////////
  // Instance

  Instance::Instance()
      : lx_(0), ly_(0), ux_(0), uy_(0), 
        isMacro_(false), isFixed_(false),
        libCellId_(-1), extId_(-1), orient_(Orientation::R0)
  {
  }

  void Instance::setLocation(int x, int y)
  {
    setBox(x, y, x + (ux_ - lx_), y + (uy_ - ly_));
  }

  void Instance::setCenterLocation(int x, int y)
  {
    int halfXLeft = (ux_ - lx_) / 2;
    int halfXRight = ux_ - lx_ - halfXLeft;
    int halfYTop = (uy_ - ly_) / 2;
    int halfYBottom = uy_ - ly_ - halfYTop;
    setBox(x - halfXLeft, y - halfYTop, x + halfXRight, y + halfYBottom);
  }

  void Instance::setSize(const Technology& tech)
  {
    const LibCell* libcell = tech.libCell(libCellId());
    orient_ = Orientation::R0;
    setSize(libcell->sizeX(), libcell->sizeY());
    for(Pin* pin : pins_)
    {
      LibPin* libpin = libcell->libPin(pin->libPinId());
      pin->updateLocation(this, libpin->x(), libpin->y());
    }
  }

  void Instance::setOrientSize(const Technology& tech, Orientation ori)
  {
    const LibCell* libcell = tech.libCell(libCellId());

    orient_ = ori;
    switch(ori)
    {
    case Orientation::R0:
    case Orientation::R180:
    case Orientation::R360:
      setSize(libcell->sizeX(), libcell->sizeY());
      break;
    case Orientation::R90:
    case Orientation::R270:
      setSize(libcell->sizeY(), libcell->sizeX());
      break;
    };

    switch(ori)
    {
    case Orientation::R0:
    case Orientation::R360:
      for(Pin* pin : pins_)
      {
        LibPin* libpin = libcell->libPin(pin->libPinId());
        pin->updateLocation(this, libpin->x(), libpin->y());
      }
      break;
    case Orientation::R90:
      for(Pin* pin : pins_)
      {
        LibPin* libpin = libcell->libPin(pin->libPinId());
        pin->updateLocation(this, -libpin->y(), libpin->x());
      }
      break;
    case Orientation::R180:
      for(Pin* pin : pins_)
      {
        LibPin* libpin = libcell->libPin(pin->libPinId());
        pin->updateLocation(this, -libpin->x(), -libpin->y());
      }
      break;
    case Orientation::R270:
      for(Pin* pin : pins_)
      {
        LibPin* libpin = libcell->libPin(pin->libPinId());
        pin->updateLocation(this, libpin->y(), -libpin->x());
      }
      break;
    };
  }

  void Instance::setSize(int w, int h)
  {
    int x = cx(), y = cy();
    int hw1 = w / 2;
    int hw2 = w - hw1;
    int hh1 = h / 2;
    int hh2 = h - hh1;
    setBox(x - hw1, y - hh1, x + hw2, y + hh2);
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
    pin->setInstance(this);
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

  void Pin::setOffset(int offsetX, int offsetY)
  {
    offsetCx_ = offsetX;
    offsetCy_ = offsetY;
    cx_ = offsetCx_ + (inst_ == nullptr ? 0 : inst_->cx());
    cy_ = offsetCy_ + (inst_ == nullptr ? 0 : inst_->cy());
  }

  void Pin::setLocation(int cx, int cy)
  {
    cx_ = cx;
    cy_ = cy;
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
    pin->setNet(this);
  }

  void Net::removePin(Pin *pin){
    pins_.erase(
      std::remove(pins_.begin(), pins_.end(), pin),
      pins_.end()
    );
  }

  ////////////////////////////////////////////////////////
  // Die

  Die::Die() 
      : dieLx_(0), dieLy_(0), dieUx_(0), dieUy_(0),
        coreLx_(0), coreLy_(0), coreUx_(0), coreUy_(0),
        rowStartX_(0), rowStartY_(0), rowWidth_(0), rowHeight_(0),
        rowRepeatCount_(0), isSetRow_(false)
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

  void Die::copyDieBoxAndRowParamsFrom(const Die & ano)
  {
    setDieBox(ano);
    setRowParams(ano);
  }

  void replace::Die::setDieBox(const Die & ano)
  {
    dieLx_ = ano.dieLx();
    dieLy_ = ano.dieLy();
    dieUx_ = ano.dieUx();
    dieUy_ = ano.dieUy();
  }

  void Die::setRowParams(int startX, int startY, int width, int height, int repeatCount)
  {
    rowStartX_ = startX;
    rowStartY_ = startY;
    rowWidth_ = width;
    rowHeight_ = height;
    rowRepeatCount_ = repeatCount;
    isSetRow_ = true;

    coreLx_ = rowStartX_;
    coreLy_ = rowStartY_;
    coreUx_ = rowStartX_ + rowWidth_;
    coreUy_ = rowStartY_ + rowHeight_ * rowRepeatCount_;
  }

  void Die::setRowParams(const Die &ano)
  {
    rowStartX_      = ano.rowStartX()      ;
    rowStartY_      = ano.rowStartY()      ;
    rowWidth_       = ano.rowWidth()       ; 
    rowHeight_      = ano.rowHeight()      ;
    rowRepeatCount_ = ano.rowRepeatCount() ;
    isSetRow_       = ano.isSetRow_        ;

    coreLx_         = ano.coreLx()         ;
    coreLy_         = ano.coreLy()         ;
    coreUy_         = ano.coreUy()         ;
    coreUx_         = ano.coreUx()         ;
  }

  void Die::addInstance(Instance *inst)
  {
    insts_.push_back(inst);
  }

  void Die::removeInstance(Instance *inst)
  {
    insts_.erase(std::remove(insts_.begin(), insts_.end(), inst), insts_.end());
  }

  ////////////////////////////////////////////////////////
  // PlacerBase

  PlacerBase::PlacerBase()
      : termSizeX_(0), termSizeY_(0),
        termSpace_(0), termCost_(0)
  {
  }

  Instance& PlacerBase::emplaceInstance(bool isFixed, bool isMacro){
    // TODO: instStor_ maybe moves to get more space
    instStor_.emplace_back();
    Instance& inst = instStor_.back();
    inst.setFixed(isFixed); inst.setMacro(isMacro);
    insts_.push_back(&inst);
    return inst;
  }

  Pin& PlacerBase::emplacePin(){
    pinStor_.emplace_back();
    pins_.push_back(&pinStor_.back());
    return pinStor_.back();
  }

  void PlacerBase::addNet(const Net& net){
    netStor_.push_back(std::move(net));
    // address of net will change after the above move 
    nets_.push_back(&netStor_.back());
    // point to the new address of net
    for(Pin* pin : net.pins()){
      pin->setNet(nets_.back());
    }
  }

  int64_t PlacerBase::hpwl()
  {
    int64_t hpwl = 0;
    for (auto &net : nets_)
    {
      net->updateBox();
      hpwl += net->hpwl();
    }
    return hpwl;
  }

  void PlacerBase::printDebugInfo() const
  {
    for (Technology* tech : techs_)
    {
      LOG_DEBUG("Technology {}", tech->name());
      tech->printDebugInfo();
    }

    LOG_DEBUG("NumInstances: {}", instStor_.size());
    // LOG_DEBUG("NumPlaceInstances: {}", placeInsts_.size());
    // LOG_DEBUG("NumFixedInstances: {}", fixedInsts_.size());
    LOG_DEBUG("NumNets: {}", nets_.size());
    LOG_DEBUG("NumPins: {}", pins_.size());

    int maxFanout = INT_MIN;
    int sumFanout = 0;
    for (const Net *net : nets_)
    {
      maxFanout = std::max((int)net->pins().size(), maxFanout);
      sumFanout += (int)net->pins().size();
    }

    LOG_DEBUG("MaxFanout: {}", maxFanout);
    LOG_DEBUG("AvgFanout: {}", sumFanout / (float)nets_.size());

    Die *die = dies_.front();
    LOG_DEBUG("CoreAreaLxLy: ({}, {})", die->coreLx(), die->coreLy());
    LOG_DEBUG("CoreAreaUxUy: ({}, {})", die->coreUx(), die->coreUy());

    int64_t coreArea =
        static_cast<int64_t>(die->coreUx() - die->coreLx()) *
        static_cast<int64_t>(die->coreUy() - die->coreLy());

    LOG_DEBUG("CoreArea: {}", coreArea);
    int dieIdx = 0;
    for (const Die *die : dies_)
    {
      LOG_DEBUG("Die {}", dieIdx);

      int numStdCell = 0, numMacro = 0;
      for (Instance* inst : die->insts())
      {
        if (inst->isMacro())
          numMacro++;
        else
          numStdCell++;
      }
      LOG_DEBUG("NumStdCells: {}", numStdCell);
      LOG_DEBUG("NumMacros: {}", numMacro);

      dieIdx++;
    }
  }

  Die *PlacerBase::die(const std::string &name) const
  {
    auto it = dieNameMap_.find(name);
    return it == dieNameMap_.end() ? nullptr : it->second;
  }

  Technology *PlacerBase::tech(const std::string &name) const
  {
    auto it = techNameMap_.find(name);
    return it == techNameMap_.end() ? nullptr : it->second;
  }

  void PlacerBase::cleanIPNs()
  {
    insts_.clear();
    pins_.clear();
    nets_.clear();
    dies_.clear();
    techs_.clear();
  }

  void PlacerBase::deriveIPNs()
  {
    insts_.reserve(instStor_.size());
    pins_.reserve(pinStor_.size());
    nets_.reserve(netStor_.size());
    dies_.reserve(dieStor_.size());
    techs_.reserve(techStor_.size());

    for (auto &inst : instStor_)
    {
      insts_.push_back(&inst);
    }

    for (auto &pin : pinStor_)
    {
      pins_.push_back(&pin);
    }

    for (auto &net : netStor_)
    {
      nets_.push_back(&net);
    }

    for (auto& die : dieStor_)
    {
      dies_.push_back(&die);
    }

    for (auto& tech : techStor_)
    {
      techs_.push_back(&tech);
    }
  }
}
