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
      : lx_(0), ly_(0), ux_(0), uy_(0), extId_(INT_MIN)
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

  void Instance::setSize(int w, int h)
  {
    setBox(lx_, ly_, lx_ + w, ly_ + h);
  }

  void Instance::setSize(Technology &tech)
  {
    auto& libcell = *tech.libCell(libCellName());
    setSize(libcell.sizeX(), libcell.sizeY());
    for(Pin* pin : pins_)
    {
      LibPin* libpin = libcell.libPin(pin->name());
      if(libpin != nullptr)
        pin->setOffset(libpin->x(), libpin->y());
    }
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
    pin->setInstance(this);
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

  void Pin::setOffset(int offsetX, int offsetY)
  {
    offsetCx_ = offsetX;
    offsetCy_ = offsetY;
    cx_ = offsetCx_ + (inst_ == nullptr ? 0 : inst_->cx());
    cy_ = offsetCx_ + (inst_ == nullptr ? 0 : inst_->cy());
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
        fixedInstsArea_(0), placeStdcellsArea_(0), placeMacrosArea_(0),
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
    isSetRow_       = ano.isSetRow_;

    coreLx_         = ano.coreLx()         ;
    coreLy_         = ano.coreLy()         ;
    coreUy_         = ano.coreUy()         ;
    coreUx_         = ano.coreUx()         ;
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

  void Die::removeInstance(Instance *inst)
  {
    insts_.erase(std::remove(insts_.begin(), insts_.end(), inst), insts_.end());
    if (inst->isFixed())
    {
      fixedInsts_.erase(std::remove(fixedInsts_.begin(), fixedInsts_.end(), inst), fixedInsts_.end());
    }
    else
    {
      placeInsts_.erase(std::remove(placeInsts_.begin(), placeInsts_.end(), inst), placeInsts_.end());
    }

    int64_t area = static_cast<int64_t>(inst->dx()) *
                   static_cast<int64_t>(inst->dy());
    if (inst->isFixed())
    {
      fixedInstsArea_ -= area;
    }
    else if (inst->isMacro())
    {
      placeMacrosArea_ -= area;
    }
    else
    {
      placeStdcellsArea_ -= area;
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

  Instance& PlacerBase::emplaceInstance(bool isFixed, bool isMacro){
    // TODO: instStor_ maybe moves to get more space
    instStor_.emplace_back();
    Instance& inst = instStor_.back();
    inst.setFixed(isFixed); inst.setMacro(isMacro);
    insts_.push_back(&inst);
    pushToInstsByType(isFixed, inst);
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

  void PlacerBase::printDebugInfo() const
  {
    for (Technology* tech : techs_)
    {
      LOG_DEBUG("Technology {}", tech->name());
      tech->printDebugInfo();
    }

    LOG_DEBUG("NumInstances: {}", instStor_.size());
    LOG_DEBUG("NumPlaceInstances: {}", placeInsts_.size());
    LOG_DEBUG("NumFixedInstances: {}", fixedInsts_.size());
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
      LOG_DEBUG("PlaceInstArea: {}", die->placeInstsArea());
      LOG_DEBUG("PlaceStdcellsArea: {}", die->placeStdcellsArea());
      LOG_DEBUG("PlaceMacrosArea: {}", die->placeMacrosArea());
      LOG_DEBUG("FixedInstArea: {}", die->fixedInstsArea());

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

  void PlacerBase::cleanIPNs(){
    insts_.clear(); pins_.clear(); nets_.clear();
  }

  void PlacerBase::deriveIPNs()
  {
    insts_.reserve(instStor_.size());
    pins_.reserve(pinStor_.size());
    nets_.reserve(netStor_.size());

    for (auto &inst : instStor_)
    {
      insts_.push_back(&inst);
      if (inst.isFixed())
      {
        fixedInsts_.push_back(&inst);
      }
      else
      {
        placeInsts_.push_back(&inst);
      }
    }

    for (auto &pin : pinStor_)
    {
      pins_.push_back(&pin);
    }

    for (auto &net : netStor_)
    {
      nets_.push_back(&net);
    }
  }

}
