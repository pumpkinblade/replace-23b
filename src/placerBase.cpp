#include "placerBase.h"
#include "log.h"

#include <iostream>

namespace replace
{

  using namespace std;

  // static odb::Rect getCoreRectFromDb(dbSet<odb::dbRow> &rows);

  // static int fastModulo(const int input, const int ceil);

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

  // void Pin::updateCoordi(odb::dbITerm *iTerm)
  // {
  //   int offsetLx = INT_MAX;
  //   int offsetLy = INT_MAX;
  //   int offsetUx = INT_MIN;
  //   int offsetUy = INT_MIN;

  //   for (dbMPin *mPin : iTerm->getMTerm()->getMPins())
  //   {
  //     for (dbBox *box : mPin->getGeometry())
  //     {
  //       offsetLx = std::min(box->xMin(), offsetLx);
  //       offsetLy = std::min(box->yMin(), offsetLy);
  //       offsetUx = std::max(box->xMax(), offsetUx);
  //       offsetUy = std::max(box->yMax(), offsetUy);
  //     }
  //   }

  //   int lx = iTerm->getInst()->getBBox()->xMin();
  //   int ly = iTerm->getInst()->getBBox()->yMin();

  //   int instCenterX = iTerm->getInst()->getMaster()->getWidth() / 2;
  //   int instCenterY = iTerm->getInst()->getMaster()->getHeight() / 2;

  //   // Pin SHAPE is NOT FOUND;
  //   // (may happen on OpenDB bug case)
  //   if (offsetLx == INT_MAX || offsetLy == INT_MAX ||
  //       offsetUx == INT_MIN || offsetUy == INT_MIN)
  //   {

  //     // offset is center of instances
  //     offsetCx_ = offsetCy_ = 0;
  //   }
  //   // usual case
  //   else
  //   {

  //     // offset is Pin BBoxs' center, so
  //     // subtract the Origin coordinates (e.g. instCenterX, instCenterY)
  //     //
  //     // Transform coordinates
  //     // from (origin: 0,0)
  //     // to (origin: instCenterX, instCenterY)
  //     //
  //     offsetCx_ = (offsetLx + offsetUx) / 2 - instCenterX;
  //     offsetCy_ = (offsetLy + offsetUy) / 2 - instCenterY;
  //   }

  //   cx_ = lx + instCenterX + offsetCx_;
  //   cy_ = ly + instCenterY + offsetCy_;
  // }

  //
  // for BTerm, offset* will hold bbox info.
  //
  // void Pin::updateCoordi(odb::dbBTerm *bTerm)
  // {
  //   int lx = INT_MAX;
  //   int ly = INT_MAX;
  //   int ux = INT_MIN;
  //   int uy = INT_MIN;

  //   for (dbBPin *bPin : bTerm->getBPins())
  //   {
  //     lx = std::min(bPin->getBox()->xMin(), lx);
  //     ly = std::min(bPin->getBox()->yMin(), ly);
  //     ux = std::max(bPin->getBox()->xMax(), ux);
  //     uy = std::max(bPin->getBox()->yMax(), uy);
  //   }

  //   if (lx == INT_MAX || ly == INT_MAX ||
  //       ux == INT_MIN || uy == INT_MIN)
  //   {
  //     string msg = string(bTerm->getConstName()) + " toplevel port is not placed!\n";
  //     msg += "       Replace will regard " + string(bTerm->getConstName()) + " is placed in (0, 0)";
  //     slog_->warn(msg, 1);
  //   }

  //   // Just center
  //   offsetCx_ = offsetCy_ = 0;

  //   cx_ = (lx + ux) / 2;
  //   cy_ = (ly + uy) / 2;
  // }

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
  // Row

  Row::Row() : siteWidth_(0), lx_(0), ly_(0), ux_(0), uy_(0)
  {
  }

  Row::Row(int siteWidth, int lx, int ly, int ux, int uy)
      : siteWidth_(siteWidth), lx_(lx), ly_(ly), ux_(ux), uy_(uy)
  {
  }

  Row::~Row()
  {
    siteWidth_ = lx_ = ly_ = ux_ = uy_ = 0;
  }

  ////////////////////////////////////////////////////////
  // Die

  Die::Die() : dieLx_(0), dieLy_(0), dieUx_(0), dieUy_(0),
               coreLx_(0), coreLy_(0), coreUx_(0), coreUy_(0) {}

  Die::~Die()
  {
    dieLx_ = dieLy_ = dieUx_ = dieUy_ = 0;
    coreLx_ = coreLy_ = coreUx_ = coreUy_ = 0;
    rows_.clear();
  }

  void Die::setDieBox(int lx, int ly, int ux, int uy)
  {
    dieLx_ = lx;
    dieLy_ = ly;
    dieUx_ = ux;
    dieUy_ = uy;
  }

  void Die::addRow(const Row &row)
  {
    rows_.push_back(row);
  }

  void Die::updateCoreBox()
  {
    coreLx_ = INT_MAX;
    coreLy_ = INT_MAX;
    coreUx_ = INT_MIN;
    coreUy_ = INT_MIN;

    for (const Row &row : rows_)
    {
      coreLx_ = std::min(row.lx(), coreLx_);
      coreLy_ = std::min(row.ly(), coreLy_);
      coreUx_ = std::max(row.ux(), coreUx_);
      coreUy_ = std::max(row.uy(), coreUy_);
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

  // void PlacerBase::init()
  // {
  //   slog_ = log_;

  //   // log_->infoInt("DBU", db_->getTech()->getDbUnitsPerMicron());

  //   // dbBlock *block = db_->getChip()->getBlock();
  //   // dbSet<dbInst> insts = block->getInsts();

  //   // // die-core area update
  //   // dbSet<dbRow> rows = block->getRows();
  //   // odb::Rect coreRect = getCoreRectFromDb(rows);
  //   // die_ = Die(block->getBBox(), &coreRect);

  //   // // siteSize update
  //   // dbRow *firstRow = *(rows.begin());
  //   // siteSizeX_ = firstRow->getSite()->getWidth();
  //   // siteSizeY_ = firstRow->getSite()->getHeight();

  //   // log_->infoIntPair("SiteSize", siteSizeX_, siteSizeY_);
  //   // log_->infoIntPair("CoreAreaLxLy", die_.coreLx(), die_.coreLy());
  //   // log_->infoIntPair("CoreAreaUxUy", die_.coreUx(), die_.coreUy());

  //   // // insts fill with real instances
  //   // instStor_.reserve(insts.size());
  //   // for (dbInst *inst : insts)
  //   // {
  //   //   Instance myInst(inst);
  //   //   instStor_.push_back(myInst);
  //   // }

  //   // insts fill with fake instances (fragmented row)
  //   initInstsForFragmentedRow();

  //   // init inst ptrs and areas
  //   insts_.reserve(instStor_.size());
  //   for (auto &inst : instStor_)
  //   {
  //     if (inst.isInstance())
  //     {
  //       if (inst.isFixed())
  //       {
  //         // Check whether fixed instance is
  //         // within the corearea
  //         //
  //         // outside of corearea is none of RePlAce's business
  //         if (isCoreAreaOverlap(die_, inst))
  //         {
  //           fixedInsts_.push_back(&inst);
  //           nonPlaceInsts_.push_back(&inst);
  //           nonPlaceInstsArea_ +=
  //               getOverlapWithCoreArea(die_, inst);
  //         }
  //       }
  //       else
  //       {
  //         placeInsts_.push_back(&inst);
  //         int64_t instArea = static_cast<int64_t>(inst.dx()) * static_cast<int64_t>(inst.dy());
  //         placeInstsArea_ += instArea;
  //         // macro cells should be
  //         // macroInstsArea_
  //         if (inst.dy() > siteSizeY_ * 6)
  //         {
  //           macroInstsArea_ += instArea;
  //         }
  //         // smaller or equal height cells should be
  //         // stdInstArea_
  //         else
  //         {
  //           stdInstsArea_ += instArea;
  //         }
  //       }
  //       instMap_[inst.dbInst()] = &inst;
  //     }
  //     else if (inst.isDummy())
  //     {
  //       dummyInsts_.push_back(&inst);
  //       nonPlaceInsts_.push_back(&inst);
  //       nonPlaceInstsArea_ += static_cast<int64_t>(inst.dx()) * static_cast<int64_t>(inst.dy());
  //     }
  //     insts_.push_back(&inst);
  //   }

  //   // nets fill
  //   dbSet<dbNet> nets = block->getNets();
  //   netStor_.reserve(nets.size());
  //   for (dbNet *net : nets)
  //   {
  //     dbSigType netType = net->getSigType();

  //     // escape nets with VDD/VSS/reset nets
  //     if (netType == dbSigType::GROUND ||
  //         netType == dbSigType::POWER ||
  //         netType == dbSigType::RESET)
  //     {
  //       continue;
  //     }

  //     Net myNet(net);
  //     netStor_.push_back(myNet);

  //     // this is safe because of "reserve"
  //     Net *myNetPtr = &netStor_[netStor_.size() - 1];
  //     netMap_[net] = myNetPtr;

  //     for (dbITerm *iTerm : net->getITerms())
  //     {
  //       Pin myPin(iTerm);
  //       myPin.setNet(myNetPtr);
  //       myPin.setInstance(dbToPlace(iTerm->getInst()));
  //       pinStor_.push_back(myPin);
  //     }

  //     for (dbBTerm *bTerm : net->getBTerms())
  //     {
  //       Pin myPin(bTerm);
  //       myPin.setNet(myNetPtr);
  //       pinStor_.push_back(myPin);
  //     }
  //   }

  //   // pinMap_ and pins_ update
  //   pins_.reserve(pinStor_.size());
  //   for (auto &pin : pinStor_)
  //   {
  //     if (pin.isITerm())
  //     {
  //       pinMap_[(void *)pin.dbITerm()] = &pin;
  //     }
  //     else if (pin.isBTerm())
  //     {
  //       pinMap_[(void *)pin.dbBTerm()] = &pin;
  //     }
  //     pins_.push_back(&pin);
  //   }

  //   // instStor_'s pins_ fill
  //   for (auto &inst : instStor_)
  //   {
  //     if (!inst.isInstance())
  //     {
  //       continue;
  //     }
  //     for (dbITerm *iTerm : inst.dbInst()->getITerms())
  //     {
  //       // note that, DB's ITerm can have
  //       // VDD/VSS pins.
  //       //
  //       // Escape those pins
  //       Pin *curPin = dbToPlace(iTerm);
  //       if (curPin)
  //       {
  //         inst.addPin(curPin);
  //       }
  //     }
  //   }

  //   // nets' pin update
  //   nets_.reserve(netStor_.size());
  //   for (auto &net : netStor_)
  //   {
  //     for (dbITerm *iTerm : net.dbNet()->getITerms())
  //     {
  //       net.addPin(dbToPlace(iTerm));
  //     }
  //     for (dbBTerm *bTerm : net.dbNet()->getBTerms())
  //     {
  //       net.addPin(dbToPlace(bTerm));
  //     }
  //     nets_.push_back(&net);
  //   }

  //   printInfo();
  // }

  // void PlacerBase::initInstsForFragmentedRow()
  // {
  //   dbSet<dbRow> rows = db_->getChip()->getBlock()->getRows();

  //   // dummy cell update to understand fragmented-row
  //   //

  //   int siteCountX = (die_.coreUx() - die_.coreLx()) / siteSizeX_;
  //   int siteCountY = (die_.coreUy() - die_.coreLy()) / siteSizeY_;

  //   enum PlaceInfo
  //   {
  //     Empty,
  //     Row,
  //     FixedInst
  //   };

  //   //
  //   // Initialize siteGrid as empty
  //   //
  //   std::vector<PlaceInfo>
  //       siteGrid(
  //           siteCountX * siteCountY,
  //           PlaceInfo::Empty);

  //   // fill in rows' bbox
  //   for (dbRow *row : rows)
  //   {
  //     Rect rect;
  //     row->getBBox(rect);

  //     std::pair<int, int> pairX = getMinMaxIdx(rect.xMin(), rect.xMax(),
  //                                              die_.coreLx(), siteSizeX_, 0, siteCountX);

  //     std::pair<int, int> pairY = getMinMaxIdx(rect.yMin(), rect.yMax(),
  //                                              die_.coreLy(), siteSizeY_, 0, siteCountY);

  //     for (int i = pairX.first; i < pairX.second; i++)
  //     {
  //       for (int j = pairY.first; j < pairY.second; j++)
  //       {
  //         siteGrid[j * siteCountX + i] = Row;
  //       }
  //     }
  //   }

  //   // fill fixed instances' bbox
  //   for (auto &inst : instStor_)
  //   {
  //     if (!inst.isFixed())
  //     {
  //       continue;
  //     }
  //     std::pair<int, int> pairX = getMinMaxIdx(inst.lx(), inst.ux(),
  //                                              die_.coreLx(), siteSizeX_, 0, siteCountX);
  //     std::pair<int, int> pairY = getMinMaxIdx(inst.ly(), inst.uy(),
  //                                              die_.coreLy(), siteSizeY_, 0, siteCountY);

  //     for (int i = pairX.first; i < pairX.second; i++)
  //     {
  //       for (int j = pairY.first; j < pairY.second; j++)
  //       {
  //         siteGrid[j * siteCountX + i] = FixedInst;
  //       }
  //     }
  //   }

  //   //
  //   // Search the "Empty" coordinates on site-grid
  //   // --> These sites need to be dummyInstance
  //   //
  //   for (int j = 0; j < siteCountY; j++)
  //   {
  //     for (int i = 0; i < siteCountX; i++)
  //     {
  //       // if empty spot found
  //       if (siteGrid[j * siteCountX + i] == Empty)
  //       {
  //         int startX = i;
  //         // find end points
  //         while (i < siteCountX &&
  //                siteGrid[j * siteCountX + i] == Empty)
  //         {
  //           i++;
  //         }
  //         int endX = i;
  //         Instance myInst(
  //             die_.coreLx() + siteSizeX_ * startX,
  //             die_.coreLy() + siteSizeY_ * j,
  //             die_.coreLx() + siteSizeX_ * endX,
  //             die_.coreLy() + siteSizeY_ * (j + 1));
  //         instStor_.push_back(myInst);
  //       }
  //     }
  //   }
  // }

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
    float util = 
      static_cast<float>(placeInstsArea_) 
      / (coreArea - nonPlaceInstsArea_) * 100;

    LOG_INFO("CoreArea: {}", coreArea);
    LOG_INFO("NonPlaceInstArea: {}", nonPlaceInstsArea_);
    LOG_INFO("PlaceInstsArea: {}", placeInstsArea_);
    LOG_INFO("Util(%): {}", util);
    LOG_INFO("StdInstArea: {}", stdInstsArea_);
    LOG_INFO("MacroInstArea: {}", macroInstsArea_);
  }

  // static odb::Rect getCoreRectFromDb(dbSet<odb::dbRow> &rows)
  // {
  //   int minX = INT_MAX, minY = INT_MAX;
  //   int maxX = INT_MIN, maxY = INT_MIN;

  //   for (dbRow *row : rows)
  //   {
  //     Rect rowRect;
  //     row->getBBox(rowRect);

  //     minX = std::min(rowRect.xMin(), minX);
  //     minY = std::min(rowRect.yMin(), minY);
  //     maxX = std::max(rowRect.xMax(), maxX);
  //     maxY = std::max(rowRect.yMax(), maxY);
  //   }
  //   return odb::Rect(minX, minY, maxX, maxY);
  // }

  // https://stackoverflow.com/questions/33333363/built-in-mod-vs-custom-mod-function-improve-the-performance-of-modulus-op
  static int fastModulo(const int input, const int ceil)
  {
    return input >= ceil ? input % ceil : input;
  }

  static std::pair<int, int> getMinMaxIdx(int ll, int uu, int coreLL, int siteSize, int minIdx, int maxIdx)
  {
    int lowerIdx = (ll - coreLL) / siteSize;
    int upperIdx =
        (fastModulo((uu - coreLL), siteSize) == 0) ? (uu - coreLL) / siteSize : (uu - coreLL) / siteSize + 1;
    return std::make_pair(
        std::max(minIdx, lowerIdx),
        std::min(maxIdx, upperIdx));
  }

  static bool isCoreAreaOverlap(Die &die, Instance &inst)
  {
    int rectLx = max(die.coreLx(), inst.lx()),
        rectLy = max(die.coreLy(), inst.ly()),
        rectUx = min(die.coreUx(), inst.ux()),
        rectUy = min(die.coreUy(), inst.uy());
    return !(rectLx >= rectUx || rectLy >= rectUy);
  }

  static int64_t getOverlapWithCoreArea(Die &die, Instance &inst)
  {
    int rectLx = max(die.coreLx(), inst.lx()),
        rectLy = max(die.coreLy(), inst.ly()),
        rectUx = min(die.coreUx(), inst.ux()),
        rectUy = min(die.coreUy(), inst.uy());
    return static_cast<int64_t>(rectUx - rectLx) * static_cast<int64_t>(rectUy - rectLy);
  }

}
