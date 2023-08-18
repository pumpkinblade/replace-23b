#include "nesterovBase.h"
#include "placerBase.h"
#include "fft.h"
#include "log.h"

#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>


#define REPLACE_SQRT2 1.414213562373095048801L

namespace replace {

using namespace std;

static prec getOverlapArea(const Bin* bin, const Instance* inst);
static prec getOverlapDensityArea(const Bin* bin, const GCell* cell);
static prec fastExp(prec exp);
static prec bellShape(prec theta_i, prec theta);
static prec bellShape1(prec theta_i);
static prec bellShape2(prec theta_i);
static prec bellShapeDerive(prec theta_i, prec theta);
static prec bellShapeDerive1(prec theta_i);
static prec bellShapeDerive2(prec theta_i);


////////////////////////////////////////////////
// GCell 

GCell::GCell() 
    : inst_(nullptr), isMacro_(false),
      lx_(0), ly_(0), ux_(0), uy_(0), theta_(0),
      densityScale_(0)
{
}

GCell::GCell(Instance* inst) 
    : GCell()
{
  setInstance(inst);
}

GCell::GCell(prec cx, prec cy, prec dx, prec dy) 
    : GCell() 
{
  lx_ = cx - dx / 2;
  ly_ = cy - dy / 2;
  ux_ = cx + dx / 2;
  uy_ = cy + dy / 2;
}

void GCell::setInstance(Instance* inst)
{
  inst_ = inst;
  // density coordi has the same center points.
  lx_ = static_cast<prec>(inst->lx());
  ly_ = static_cast<prec>(inst->ly());
  ux_ = static_cast<prec>(inst->ux());
  uy_ = static_cast<prec>(inst->uy());
  isMacro_ = inst->isMacro();
}

void GCell::addGPin(GPin* gPin)
{
  gPins_.push_back(gPin);
}

void GCell::setLocation(prec lx, prec ly)
{
  ux_ = lx + (ux_ - lx_);
  uy_ = ly + (uy_ - ly_);
  lx = lx_;
  ly = ly_;

  for(GPin* pin : gPins_)
    pin->updateLocation(this);
}

void GCell::setCenterLocation(prec cx, prec cy)
{
  prec halfDx = dx() / 2;
  prec halfDy = dy() / 2;

  lx_ = cx - halfDx;
  ly_ = cy - halfDy;
  ux_ = cx + halfDx;
  uy_ = cy + halfDy;

  for(GPin* pin : gPins_)
    pin->updateLocation(this);
}

// changing size and preserve center coordinates
void GCell::setSize(prec dx, prec dy)
{
  prec centerX = cx();
  prec centerY = cy();

  lx_ = centerX - dx / 2;
  ly_ = centerY - dy / 2;
  ux_ = centerX + dx / 2;
  uy_ = centerY + dy / 2;
}

void GCell::setCenterLocationTheta(prec cx, prec cy, prec theta)
{
  prec halfDx = dx() / 2;
  prec halfDy = dy() / 2;

  lx_ = cx - halfDx;
  ly_ = cy - halfDy;
  ux_ = cx + halfDx;
  uy_ = cy + halfDy;
  theta_ = theta;

  for(GPin* pin : gPins_)
    pin->updateLocationWithTheta(this);
}

////////////////////////////////////////////////
// GNet

GNet::GNet()
    : net_(nullptr),
      lx_(0), ly_(0), ux_(0), uy_(0),
      waExpMinSumX_(0), waXExpMinSumX_(0),
      waExpMaxSumX_(0), waXExpMaxSumX_(0),
      waExpMinSumY_(0), waYExpMinSumY_(0),
      waExpMaxSumY_(0), waYExpMaxSumY_(0)
{
}

GNet::GNet(Net* net)
{
  net_ = net;
}


void GNet::addGPin(GPin* gPin)
{
  gPins_.push_back(gPin);
}

void GNet::updateBox()
{
  lx_ = ly_ = std::numeric_limits<prec>::max();
  ux_ = uy_ = std::numeric_limits<prec>::lowest();

  for(auto& gPin : gPins_)
  {
    lx_ = std::min(gPin->cx(), lx_);
    ly_ = std::min(gPin->cy(), ly_);
    ux_ = std::max(gPin->cx(), ux_);
    uy_ = std::max(gPin->cy(), uy_);
  } 
}

void GNet::clearWaVars()
{
  waExpMinSumX_ = 0;
  waXExpMinSumX_ = 0;

  waExpMaxSumX_ = 0;
  waXExpMaxSumX_ = 0;
    
  waExpMinSumY_ = 0;
  waYExpMinSumY_ = 0;

  waExpMaxSumY_ = 0;
  waYExpMaxSumY_ = 0;
}

////////////////////////////////////////////////
// GPin 

GPin::GPin()
  : pin_(nullptr), gCell_(nullptr), gNet_(nullptr),
    offsetCx_(0), offsetCy_(0),
    cx_(0), cy_(0),
    maxExpSumX_(0), maxExpSumY_(0),
    minExpSumX_(0), minExpSumY_(0),
    hasMaxExpSumX_(0), hasMaxExpSumY_(0), 
    hasMinExpSumX_(0), hasMinExpSumY_(0)
{
}

GPin::GPin(Pin* pin)
    : GPin()
{
  pin_ = pin;
  cx_ = static_cast<prec>(pin->cx());
  cy_ = static_cast<prec>(pin->cy());
  offsetCx_ = static_cast<prec>(pin->offsetCx());
  offsetCy_ = static_cast<prec>(pin->offsetCy());
}

void GPin::setCenterLocation(prec cx, prec cy)
{
  cx_ = cx;
  cy_ = cy;
}

void GPin::clearWaVars()
{
  hasMaxExpSumX_ = false;
  hasMaxExpSumY_ = false;
  hasMinExpSumX_ = false;
  hasMinExpSumY_ = false;
    
  maxExpSumX_ = maxExpSumY_ = 0;
  minExpSumX_ = minExpSumY_ = 0;
}

void GPin::setMaxExpSumX(prec maxExpSumX)
{
  hasMaxExpSumX_ = true;
  maxExpSumX_ = maxExpSumX;
}

void GPin::setMaxExpSumY(prec maxExpSumY)
{
  hasMaxExpSumY_ = true;
  maxExpSumY_ = maxExpSumY;
}

void GPin::setMinExpSumX(prec minExpSumX)
{
  hasMinExpSumX_ = true;
  minExpSumX_ = minExpSumX;
}

void GPin::setMinExpSumY(prec minExpSumY)
{
  hasMinExpSumY_ = true;
  minExpSumY_ = minExpSumY;
}

void GPin::updateLocation(const GCell* gCell)
{
  cx_ = gCell->cx() + offsetCx_;
  cy_ = gCell->cy() + offsetCy_;
}

void GPin::updateLocationWithTheta(const GCell* gCell)
{
  prec cosTheta = std::cos(gCell->theta());
  prec sinTheta = std::sin(gCell->theta());
  cx_ = gCell->cx() + offsetCx_ * cosTheta - offsetCy_ * sinTheta;
  cy_ = gCell->cy() + offsetCx_ * sinTheta + offsetCy_ * cosTheta;
}

////////////////////////////////////////////////////////
// Bin

Bin::Bin() 
    : x_(0), y_(0), lx_(0), ly_(0),
      ux_(0), uy_(0), 
      nonPlaceArea_(0), instPlacedArea_(0),
      fillerArea_(0),
      density_ (0),
      targetDensity_(0),
      electroPhi_(0), 
      electroForceX_(0), electroForceY_(0)
{
}

Bin::Bin(int x, int y, prec lx, prec ly, prec ux, prec uy, prec targetDensity) 
    : Bin()
{
  x_ = x;
  y_ = y;
  lx_ = lx; 
  ly_ = ly;
  ux_ = ux;
  uy_ = uy;
  targetDensity_ = targetDensity;
}

////////////////////////////////////////////////
// BinGrid

BinGrid::BinGrid()
    : die_(nullptr),
      lx_(0), ly_(0), ux_(0), uy_(0),
      binCntX_(0), binCntY_(0),
      binSizeX_(0), binSizeY_(0),
      targetDensity_(0), overflowArea_(0),
      isSetBinCntX_(0), isSetBinCntY_(0),
      totalCellArea_(0)
{
}

BinGrid::BinGrid(Die* die)
   : BinGrid()
{
  setDie(die);
}

void BinGrid::setDie(Die* die)
{
  die_ = die;
  lx_ = static_cast<prec>(die->coreLx());
  ly_ = static_cast<prec>(die->coreLy());
  ux_ = static_cast<prec>(die->coreUx());
  uy_ = static_cast<prec>(die->coreUy());

  placeInstCnt_ = 0;
  placeStdCellArea_ = 0;
  placeMacroArea_ = 0;
  fixedInstArea_ = 0;
  for(Instance* inst : die_->insts())
  {
    double area = (double)(inst->dx()) * inst->dy();
    if(inst->isFixed())
      fixedInstArea_ += area;
    else
    {
      placeInstCnt_++;
      if(inst->isMacro())
        placeMacroArea_ += area;
      else
        placeStdCellArea_ += area;
    }
  }
}

void BinGrid::setTargetDensity(prec density)
{
  targetDensity_ = density;
}

void BinGrid::setBinCntX(int binCntX)
{
  isSetBinCntX_ = true;
  binCntX_ = binCntX;
}

void BinGrid::setBinCntY(int binCntY)
{
  isSetBinCntY_ = true;
  binCntY_ = binCntY;
}

void BinGrid::initBins()
{
  double totalBinArea = (double)(ux_ - lx_) * (uy_ - ly_);

  int foundBinCnt = 2;
  if (placeInstCnt_ != 0)
  {
    // assume no fixed instance
    double avgPlaceInstArea = (placeMacroArea_ + placeStdCellArea_) / placeInstCnt_;
    double idealBinArea = avgPlaceInstArea / targetDensity_;
    int idealBinCnt = static_cast<int>(totalBinArea / idealBinArea);
    
    LOG_DEBUG("TargetDensity: {}", targetDensity_);
    LOG_DEBUG("AveragePlaceInstArea: {}", avgPlaceInstArea);
    LOG_DEBUG("IdealBinArea: {}", idealBinArea);
    LOG_DEBUG("IdealBinCnt: {}", idealBinCnt);
    LOG_DEBUG("TotalBinArea: {}", totalBinArea);

    // find binCnt: 2, 4, 8, 16, 32, 64, ...
    // s.t. binCnt^2 <= idealBinCnt <= (binCnt*2)^2.
    for(foundBinCnt = 2; foundBinCnt <= 1024; foundBinCnt *= 2)
    {
      if(foundBinCnt * foundBinCnt <= idealBinCnt && 4 * foundBinCnt * foundBinCnt > idealBinCnt)
        break;
    }
  }
  // setBinCntX_;
  if(!isSetBinCntX_)
  {
    binCntX_ = foundBinCnt;
  }
  // setBinCntY_;
  if(!isSetBinCntY_)
  {
    binCntY_ = foundBinCnt;
  }

  LOG_DEBUG("BinCnt: {}, {}", binCntX_, binCntY_);
  
  binSizeX_ = static_cast<prec>((ux_ - lx_) / binCntX_);
  binSizeY_ = static_cast<prec>((uy_ - ly_) / binCntY_);
  
  LOG_DEBUG("BinSize: {}, {}", binSizeX_, binSizeY_);

  // initialize binStor_, bins_ vector
  binStor_.reserve(binCntX_ * binCntY_);
  bins_.reserve(binCntX_ * binCntY_);
  for (int idxY = 0; idxY < binCntY_; idxY++)
  {
    prec ly = static_cast<prec>(idxY * binSizeY_);
    prec uy = ly + binSizeY_;
    for (int idxX = 0; idxX < binCntX_; idxX++)
    {
      prec lx = static_cast<prec>(idxX * binSizeX_);
      prec ux = lx + binSizeX_;
      binStor_.emplace_back(idxX, idxY, lx, ly, ux, uy, targetDensity_);
      bins_.push_back(&binStor_.back());
    }
  }

  LOG_DEBUG("NumBins: {}", bins_.size());

  // only initialized once
  updateBinsNonPlaceArea();

  // initialize fft structrue based on bins
  fft_.init(binCntX_, binCntY_, binSizeX_, binSizeY_);
}

void BinGrid::updateBinsNonPlaceArea()
{
  for(auto& bin : bins_)
  {
    bin->setNonPlaceArea(0);
  }

  for(Instance* inst : die_->insts())
  {
    if(inst->isFixed())
    {
      auto pairX = getMinMaxIdxX(static_cast<prec>(inst->lx()),
                                 static_cast<prec>(inst->ux()));
      auto pairY = getMinMaxIdxY(static_cast<prec>(inst->ly()),
                                 static_cast<prec>(inst->uy()));
      for(int i = pairX.first; i < pairX.second; i++)
      {
        for(int j = pairY.first; j < pairY.second; j++)
        {
          Bin* bin = bins_[j * binCntX_ + i];

          // Note that nonPlaceArea should have scale-down with
          // target density. 
          // See MS-replace paper
          //
          double overlap = getOverlapArea(bin, inst);
          overlap *= bin->targetDensity();
          bin->addNonPlaceArea(static_cast<prec>(overlap));
        }
      }
    }
  }
}


// Core Part
void BinGrid::updateBinsGCellDensityArea()
{
  // clear the Bin-area info
  for(auto& bin : bins_)
  {
    bin->setInstPlacedArea(0);
    bin->setFillerArea(0);
  }

  for(auto& cell : gCells_)
  {
    // The following function is critical runtime hotspot for global placer.
    if(cell->isInstance() && cell->isMacro())
      addMacroBinArea(cell);
    else if(cell->isInstance())
      addStdCellBinArea(cell);
    else
      addFillerBinArea(cell);
  }

  overflowArea_ = 0;

  // update density and overflowArea 
  // for nesterov use and FFT library
  for(auto& bin : bins_)
  {
    double binArea = bin->binArea();
    double cellArea = bin->instPlacedArea() + bin->fillerArea() + bin->nonPlaceArea();
    prec density = static_cast<prec>(cellArea / (binArea * bin->targetDensity()));
    bin->setDensity(density);

    overflowArea_ 
      += std::max(
          0.0,
          (bin->instPlacedArea() + bin->nonPlaceArea()) - (binArea * bin->targetDensity()));
  }
}

void BinGrid::addFillerBinArea(const GCell* gcell)
{
  std::pair<int, int> pairX = getMinMaxIdxX(gcell->lx(), gcell->ux());
  std::pair<int, int> pairY = getMinMaxIdxY(gcell->ly(), gcell->uy());
  for(int i = pairX.first; i < pairX.second; i++)
  {
    for(int j = pairY.first; j < pairY.second; j++)
    {
      Bin* bin = bins_[j * binCntX_ + i];
      double overlapArea = getOverlapDensityArea(bin, gcell);
      overlapArea *= gcell->densityScale();
      bin->addFillerArea(static_cast<prec>(overlapArea));
    }
  }
}

void BinGrid::addStdCellBinArea(const GCell* gcell)
{
  std::pair<int, int> pairX = getMinMaxIdxX(gcell->lx(), gcell->ux());
  std::pair<int, int> pairY = getMinMaxIdxY(gcell->ly(), gcell->uy());
  for(int i = pairX.first; i < pairX.second; i++)
  {
    for(int j = pairY.first; j < pairY.second; j++)
    {
      Bin* bin = bins_[j * binCntX_ + i];
      double overlapArea = getOverlapDensityArea(bin, gcell);
      overlapArea *= gcell->densityScale();
      bin->addInstPlacedArea(static_cast<prec>(overlapArea));
    }
  }
}

void BinGrid::addMacroBinArea(const GCell* gcell)
{
  std::pair<int, int> pairX = getMinMaxIdxX(gcell->lx(), gcell->ux());
  std::pair<int, int> pairY = getMinMaxIdxY(gcell->ly(), gcell->uy());
  for(int i = pairX.first; i < pairX.second; i++)
  {
    for(int j = pairY.first; j < pairY.second; j++)
    {
      Bin* bin = bins_[j * binCntX_ + i];
      double overlapArea = getOverlapDensityArea(bin, gcell);
      overlapArea *= (gcell->densityScale() * bin->targetDensity());
      bin->addInstPlacedArea(static_cast<prec>(overlapArea));
    }
  }
}

void BinGrid::addMacroBinAreaWithTheta(const GCell* gcell)
{
  std::pair<int, int> pairX = getMinMaxIdxX(gcell->lx(), gcell->ux());
  std::pair<int, int> pairY = getMinMaxIdxY(gcell->ly(), gcell->uy());
  prec scale = gcell->densityScale() * bellShape1(gcell->theta());
  for(int i = pairX.first; i < pairX.second; i++)
  {
    for(int j = pairY.first; j < pairY.second; j++)
    {
      Bin* bin = bins_[j * binCntX_ + i];
      double overlapArea = getOverlapDensityArea(bin, gcell);
      overlapArea *= (scale * bin->targetDensity());
      bin->addInstPlacedArea(static_cast<prec>(overlapArea));
    }
  }

  prec rlx = gcell->cx() - 0.5f * gcell->dy();
  prec rux = gcell->cx() + 0.5f * gcell->dx();
  prec rly = gcell->cy() - 0.5f * gcell->dy();
  prec ruy = gcell->cy() + 0.5f * gcell->dy();
  pairX = getMinMaxIdxX(rlx, rux);
  pairY = getMinMaxIdxY(rly, ruy);
  scale = gcell->densityScale() * bellShape2(gcell->theta());
  for(int i = pairX.first; i < pairX.second; i++)
  {
    for(int j = pairY.first; j < pairY.second; j++)
    {
      Bin* bin = bins_[j * binCntX_ + i];
      prec overlapLx = std::max(rlx, bin->lx());
      prec overlapLy = std::max(rly, bin->ly());
      prec overlapUx = std::min(rux, bin->ux());
      prec overlapUy = std::min(ruy, bin->uy());
      if(overlapLx < overlapUx && overlapLy < overlapUy)
      {
        prec overlapArea = (overlapUx - overlapLx) * (overlapUy - overlapLy);
        overlapArea *= (scale * bin->targetDensity());
        bin->addInstPlacedArea(static_cast<prec>(overlapArea));
      }
    }
  }
}

std::pair<int, int> BinGrid::getMinMaxIdxX(prec lx1, prec ux1)
{
  int lowerIdx = static_cast<int>(std::floor((lx1 - lx_) / binSizeX_));
  int upperIdx = static_cast<int>(std::ceil((ux1 - lx_) / binSizeX_));
  return std::make_pair(std::max(lowerIdx, 0), std::min(upperIdx, binCntX_));
}

std::pair<int, int> BinGrid::getMinMaxIdxY(prec ly1, prec uy1)
{
  int lowerIdx = static_cast<int>(std::floor((ly1 - ly_) / binSizeY_));
  int upperIdx = static_cast<int>(std::ceil((uy1 - ly_) / binSizeY_));
  return std::make_pair(std::max(lowerIdx, 0), std::min(upperIdx, binCntY_));
}

void BinGrid::addGCell(GCell* gc)
{
  gCells_.push_back(gc);
  totalCellArea_ += gc->dx() * gc->dy();
}

void BinGrid::updateDensityForceBin()
{
  // copy density to utilize FFT
  for(Bin* bin : bins_)
    fft_.updateDensity(bin->x(), bin->y(), bin->density());  

  // do FFT
  fft_.doFFT();

  // update electroPhi and electroForce
  // update sumPhi_ for nesterov loop
  sumPhi_ = 0;
  for(Bin* bin : bins_)
  {
    auto eForcePair = fft_.getElectroForce(bin->x(), bin->y());
    bin->setElectroForceX(eForcePair.first);
    bin->setElectroForceY(eForcePair.second);

    prec electroPhi = fft_.getElectroPhi(bin->x(), bin->y());
    bin->setElectroPhi(electroPhi);

    sumPhi_ += electroPhi 
      * static_cast<prec>(bin->nonPlaceArea() 
          + bin->instPlacedArea() + bin->fillerArea());
  }
}

void BinGrid::updateGCellDensityScaleAndSize()
{
  // update densitySize and densityScale in each gCell
  for(auto& gCell : gCells_)
  {
    prec scaleX = 0, scaleY = 0;
    prec densitySizeX = 0, densitySizeY = 0;
    if(gCell->dx() < REPLACE_SQRT2 * binSizeX_)
    {
      scaleX = static_cast<prec>(gCell->dx()) / static_cast<prec>(REPLACE_SQRT2 * binSizeX_);
      densitySizeX = static_cast<prec>(REPLACE_SQRT2 * binSizeX_);
    }
    else
    {
      scaleX = 1.0f;
      densitySizeX = gCell->dx();
    }

    if(gCell->dy() < REPLACE_SQRT2 * binSizeY_)
    {
      scaleY = static_cast<prec>(gCell->dy()) / static_cast<prec>(REPLACE_SQRT2 * binSizeY_);
      densitySizeY = static_cast<prec>(REPLACE_SQRT2 * binSizeY_);
    }
    else
    {
      scaleY = 1.0f;
      densitySizeY = gCell->dy();
    }

    gCell->setSize(densitySizeX, densitySizeY);
    gCell->setDensityScale(scaleX * scaleY);
  }
}


////////////////////////////////////////////////
// NesterovBaseVars
NesterovBaseVars::NesterovBaseVars() 
  : targetDensity(1.0), 
    minAvgCut(0.1), maxAvgCut(0.9),
    binCntX(0), binCntY(0),
    minWireLengthForceBar(-300),
    isSetBinCntX(0), isSetBinCntY(0) 
{
}

void  NesterovBaseVars::reset()
{
  targetDensity = 1.0;
  minAvgCut = 0.1;
  maxAvgCut = 0.9;
  isSetBinCntX = isSetBinCntY = 0;
  binCntX = binCntY = 0;
  minWireLengthForceBar = -300;
}


////////////////////////////////////////////////
// NesterovBase 

NesterovBase::NesterovBase()
  : pb_(nullptr), sumPhi_()
{
}

NesterovBase::NesterovBase(
    NesterovBaseVars nbVars, 
    std::shared_ptr<PlacerBase> pb)
    : NesterovBase()
{
  nbVars_ = nbVars;
  pb_ = pb;
  init();
}

void NesterovBase::init()
{
  LOG_TRACE("start NesterovBase::init");
  std::unordered_map<Instance*, GCell*> gCellMap;
  std::unordered_map<Pin*, GPin*> gPinMap;
  std::unordered_map<Net*, GNet*> gNetMap;
  std::unordered_map<Die*, BinGrid*> binGridMap;

  // create bin grid
  binGridStor_.reserve(pb_->dies().size());
  for(Die* die : pb_->dies())
  {
    // each die should own one bin grid
    binGridStor_.emplace_back(die);

    // send param into binGrid structure
    if(nbVars_.isSetBinCntX)
    {
      binGridStor_.back().setBinCntX(nbVars_.binCntX);
    }
    if(nbVars_.isSetBinCntY)
    {
      binGridStor_.back().setBinCntY(nbVars_.binCntY);
    }
    binGridStor_.back().setTargetDensity(nbVars_.targetDensity);

    // apply the bin setting
    binGridStor_.back().initBins();
  }

  // gCellStor init
  gCellStor_.reserve(pb_->insts().size());
  for(Instance* inst : pb_->insts())
  {
    if(!inst->isFixed())
    {
      GCell myGCell(inst);
      gCellStor_.push_back(myGCell);
    }
  }

  // TODO: 
  // at this moment, GNet and GPin is equal to
  // Net and Pin

  // gPinStor init
  gPinStor_.reserve(pb_->pins().size());
  for(auto& pin : pb_->pins())
  {
    GPin myGPin(pin);
    gPinStor_.push_back(myGPin);
  }

  LOG_TRACE("in NesterovBase::init gNetStor init");
  // gNetStor init
  gNetStor_.reserve(pb_->nets().size());
  for(auto& net : pb_->nets())
  {
    GNet myGNet(net);
    gNetStor_.push_back(myGNet);
  }

  // binGrid ptr init
  binGrids_.reserve(binGridStor_.size());
  for(auto& bg : binGridStor_)
  {
    binGrids_.push_back(&bg);
    binGridMap.emplace(bg.die(), &bg);
  }

  LOG_TRACE("in NesterovBase::init create filler for each die");
  // create filler for each die
  std::unordered_map<BinGrid*, std::pair<int, int>> fillerIdxMap;
  for(BinGrid* bg : binGrids_)
  {
    int startIdx = static_cast<int>(gCellStor_.size());
    initFillerGCells(bg);
    int endIdx = static_cast<int>(gCellStor_.size());
    fillerIdxMap.emplace(bg, std::make_pair(startIdx, endIdx));
  }

  LOG_TRACE("in NesterovBase::init gCell ptr init");
  // gCell ptr init
  gCells_.reserve(gCellStor_.size());
  for(auto& gCell : gCellStor_)
  {
    gCells_.push_back(&gCell);
    if(gCell.isInstance())
      gCellMap.emplace(gCell.instance(), &gCell);
  }

  LOG_TRACE("in NesterovBase::init gPin ptr init");
  // gPin ptr init
  gPins_.reserve(gPinStor_.size());
  for(auto& gPin : gPinStor_)
  {
    gPins_.push_back(&gPin);
    gPinMap.emplace(gPin.pin(), &gPin);
  }

  // gNet ptr init
  gNets_.reserve(gNetStor_.size());
  for(auto& gNet : gNetStor_)
  {
    gNets_.push_back(&gNet);
    gNetMap.emplace(gNet.net(), &gNet);
  }

  // gCellStor_'s pins_ fill
  for(auto& gCell : gCellStor_)
  {
    if(gCell.isInstance())
    {
      for(auto& pin : gCell.instance()->pins())
      {
        auto it = gPinMap.find(pin);
        gCell.addGPin(it == gPinMap.end() ? nullptr : it->second);
      }
    }
  }

  // gPinStor_' GNet and GCell fill
  for(auto& gPin : gPinStor_)
  {
    auto iit = gCellMap.find(gPin.pin()->instance());
    auto nit = gNetMap.find(gPin.pin()->net());
    gPin.setGCell(iit == gCellMap.end() ? nullptr : iit->second);
    gPin.setGNet(nit == gNetMap.end() ? nullptr : nit->second);
  } 

  // gNetStor_'s GPin fill
  for(auto& gNet : gNetStor_)
  {
    for(auto& pin : gNet.net()->pins())
    {
      // NOTE: What if a pin is on a fixed instances?
      auto it = gPinMap.find(pin);
      gNet.addGPin(it == gPinMap.end() ? nullptr : it->second);
    }
  }

  LOG_DEBUG("FillerInit: NumGCells: {}", gCells_.size());
  LOG_DEBUG("FillerInit: NumGNets: {}", gNets_.size());
  LOG_DEBUG("FillerInit: NumGPins: {}", gPins_.size());

  LOG_TRACE("in NesterovBase::init binGrids_ gCell fill");
  // binGrids_ gCell fill
  for(BinGrid* bg : binGrids_)
  {
    // assign inst gCell to binGrid
    for(Instance* inst : bg->die()->insts())
    {
      if(!inst->isFixed())
      {
        auto it = gCellMap.find(inst);
        GCell* gCell = (it == gCellMap.end() ? nullptr : it->second);
        bg->addGCell(gCell);
        gCell->setBinGrid(bg);
      }
    }
    // assign filler gCell to binGrid
    auto idx = fillerIdxMap[bg];
    for(int i = idx.first; i < idx.second; i++)
    {
      bg->addGCell(&gCellStor_[i]);
      gCellStor_[i].setBinGrid(bg);
    }

    // update densitySize and densityScale in each gCell
    bg->updateGCellDensityScaleAndSize();
  }
  LOG_TRACE("finish NesterovBase::init");
}


// virtual filler GCells
void NesterovBase::initFillerGCells(BinGrid* bg)
{
  // if die doesn't have any instance, we will not create any filler.
  if(bg->placeInstanceCount() == 0)
    return;

  double dxSum = 0.0, dySum = 0.0;
  int numStdcell = 0;
  for(Instance* inst : bg->die()->insts())
  {
    if(!inst->isFixed())
    {
      dxSum += inst->dx();
      dySum += inst->dy();
      numStdcell++;
    }
  }
  prec avgDx = static_cast<prec>(dxSum / numStdcell);
  prec avgDy = static_cast<prec>(dySum / numStdcell);

  double coreArea = (double)bg->die()->coreDx() * bg->die()->coreDy();

  // nonPlaceInstsArea should not have targetDensity downscaling!!! 
  double whiteSpaceArea = coreArea - bg->fixedInstanceArea();

  // TODO density screening
  double movableArea = whiteSpaceArea * nbVars_.targetDensity;
  
  double totalFillerArea = movableArea - bg->placeStdCellArea() 
                          - bg->placeMacroArea() * nbVars_.targetDensity;

  if(totalFillerArea < 0)
  {
    LOG_ERROR("Filler area is negative!!\n"
              "\tPlease put higher target density or\n"
              "\tRe-floorplan to have enough coreArea");
  }

  assert(avgDx > 0);
  assert(avgDy > 0);
  int fillerCnt = static_cast<int>(totalFillerArea / ((double)avgDx * avgDy));

  LOG_DEBUG("FillerInit: CoreArea: {}", coreArea);
  LOG_DEBUG("FillerInit: WhiteSpaceArea: {}", whiteSpaceArea);
  LOG_DEBUG("FillerInit: MovableArea: {}", movableArea);
  LOG_DEBUG("FillerInit: TotalFillerArea: {}", totalFillerArea);
  LOG_DEBUG("FillerInit: NumFillerCells: {}", fillerCnt);
  LOG_DEBUG("FillerInit: FillerCellArea: {}", (double)avgDx * avgDy);
  LOG_DEBUG("FillerInit: FillerCellSize: {}, {}", avgDx, avgDy); 

  // 
  // mt19937 supports huge range of random values.
  // rand()'s RAND_MAX is only 32767.
  //
  mt19937 randVal(0);
  for(int i = 0; i < fillerCnt; i++)
  {
    // instability problem between g++ and clang++!
    auto randX = randVal();
    auto randY = randVal();

    // place filler cells on random coordi and
    // set size as avgDx and avgDy
    prec cx = static_cast<prec>(randX % bg->die()->coreDx() + bg->die()->coreLx());
    prec cy = static_cast<prec>(randY % bg->die()->coreDy() + bg->die()->coreLy());
    GCell myGCell(cx, cy, avgDx, avgDy);
    gCellStor_.push_back(myGCell);
  }
}

// gcell update
void NesterovBase::updateGCellLocation(std::vector<Point>& coordis)
{
  for(auto& coordi : coordis)
  {
    int idx = static_cast<int>(&coordi - &coordis[0]);
    gCells_[idx]->setLocation( coordi.x, coordi.y );
  }
  for (BinGrid* bg : binGrids_)
  {
    bg->updateBinsGCellDensityArea();
  }
}

// gcell update
void NesterovBase::updateGCellCenterLocation(std::vector<Point>& coordis)
{
  for(auto& coordi : coordis)
  {
    int idx = static_cast<int>(&coordi - &coordis[0]);
    gCells_[idx]->setCenterLocation( coordi.x, coordi.y );
  }
  for (BinGrid* bg : binGrids_)
  {
    bg->updateBinsGCellDensityArea();
  }
}

double  NesterovBase::overflowArea() const
{
  double area = 0;
  for(BinGrid* bg : binGrids_)
  {
    area += bg->overflowArea();
  }
  return area;
}

prec NesterovBase::sumPhi() const
{
  return sumPhi_;
}

prec NesterovBase::targetDensity() const
{
  return nbVars_.targetDensity;
}

void NesterovBase::updateDensityCoordiLayoutInside(GCell* gCell)
{
  prec targetLx = gCell->lx();
  prec targetLy = gCell->ly();

  BinGrid* bg = gCell->binGrid();
  targetLx = std::max(std::min(targetLx, bg->ux() - gCell->dx()), bg->lx());
  targetLy = std::max(std::min(targetLy, bg->uy() - gCell->dy()), bg->ly());
  gCell->setLocation(targetLx, targetLy);
}

prec NesterovBase::getDensityCoordiLayoutInsideX(GCell* gCell, prec cx)
{
  prec adjVal = cx;
  BinGrid* bg = gCell->binGrid();
  adjVal = std::max(adjVal, bg->lx() + 0.5f * gCell->dx());
  adjVal = std::min(adjVal, bg->ux() - 0.5f * gCell->dx());
  return adjVal;
}

prec NesterovBase::getDensityCoordiLayoutInsideY(GCell* gCell, prec cy)
{
  prec adjVal = cy;
  BinGrid* bg = gCell->binGrid();
  adjVal = std::max(adjVal, bg->ly() + 0.5f * gCell->dy());
  adjVal = std::min(adjVal, bg->uy() - 0.5f * gCell->dy());
  return adjVal;
}

// 
// WA force cals - wlCoeffX / wlCoeffY
//
// * Note that wlCoeffX and wlCoeffY is 1/gamma 
// in ePlace paper.
void NesterovBase::updateWireLengthForceWA(
    prec wlCoeffX, prec wlCoeffY)
{
  // clear all WA variables.
  for(auto& gNet : gNets_)
  {
    gNet->clearWaVars();
  }
  for(auto& gPin : gPins_)
  {
    gPin->clearWaVars();
  }

  for(auto& gNet : gNets_)
  {
    gNet->updateBox();

    for(auto& gPin : gNet->gPins())
    {
      prec expMinX = (gNet->lx() - gPin->cx()) * wlCoeffX; 
      prec expMaxX = (gPin->cx() - gNet->ux()) * wlCoeffX;
      prec expMinY = (gNet->ly() - gPin->cy()) * wlCoeffY;
      prec expMaxY = (gPin->cy() - gNet->uy()) * wlCoeffY;

      // min x
      if(expMinX > nbVars_.minWireLengthForceBar)
      {
        gPin->setMinExpSumX( fastExp(expMinX) );
        gNet->addWaExpMinSumX( gPin->minExpSumX() );
        gNet->addWaXExpMinSumX( gPin->cx() 
            * gPin->minExpSumX() );
      }
      
      // max x
      if(expMaxX > nbVars_.minWireLengthForceBar)
      {
        gPin->setMaxExpSumX( fastExp(expMaxX) );
        gNet->addWaExpMaxSumX( gPin->maxExpSumX() );
        gNet->addWaXExpMaxSumX( gPin->cx() 
            * gPin->maxExpSumX() );
      }
     
      // min y 
      if(expMinY > nbVars_.minWireLengthForceBar)
      {
        gPin->setMinExpSumY( fastExp(expMinY) );
        gNet->addWaExpMinSumY( gPin->minExpSumY() );
        gNet->addWaYExpMinSumY( gPin->cy() 
            * gPin->minExpSumY() );
      }
      
      // max y
      if(expMaxY > nbVars_.minWireLengthForceBar)
      {
        gPin->setMaxExpSumY( fastExp(expMaxY) );
        gNet->addWaExpMaxSumY( gPin->maxExpSumY() );
        gNet->addWaYExpMaxSumY( gPin->cy() 
            * gPin->maxExpSumY() );
      }
    }
    //cout << gNet->lx() << " " << gNet->ly() << " "
    //  << gNet->ux() << " " << gNet->uy() << endl;
  }
}

// get x,y WA Gradient values with given GCell
Point NesterovBase::getWireLengthGradientWA(GCell* gCell, prec wlCoeffX, prec wlCoeffY)
{
  Point gradientPair;

  for(auto& gPin : gCell->gPins())
  {
    auto tmpPair = getWireLengthGradientPinWA(gPin, wlCoeffX, wlCoeffY);
    gradientPair.x += tmpPair.x;
    gradientPair.y += tmpPair.y;
  }

  // return sum
  return gradientPair;
}

Point NesterovBase::getWireLengthGradientWAWithTheta(GCell* gCell, prec wlCoeffX, prec wlCoeffY, prec& gradTheta)
{
  Point gradientPair;
  gradTheta = 0.0f;
  
  prec cosTheta = std::cos(gCell->theta());
  prec sinTheta = std::sin(gCell->theta());
  for(auto& gPin : gCell->gPins())
  {
    auto tmpPair = getWireLengthGradientPinWA(gPin, wlCoeffX, wlCoeffY);
    gradientPair.x += tmpPair.x;
    gradientPair.y += tmpPair.y;
    gradTheta += tmpPair.x * (-gPin->offsetCx() * sinTheta - gPin->offsetCy() * cosTheta);
    gradTheta += tmpPair.y * ( gPin->offsetCx() * cosTheta - gPin->offsetCy() * sinTheta);
  }

  // return sum
  return gradientPair;
}

// get x,y WA Gradient values from GPin
// Please check the JingWei's Ph.D. thesis full paper, 
// Equation (4.13)
//
// You can't understand the following function
// unless you read the (4.13) formula
Point NesterovBase::getWireLengthGradientPinWA(GPin* gPin, prec wlCoeffX, prec wlCoeffY)
{
  prec gradientMinX = 0, gradientMinY = 0;
  prec gradientMaxX = 0, gradientMaxY = 0;

  // min x
  if( gPin->hasMinExpSumX() ) {
    // from Net.
    prec waExpMinSumX = gPin->gNet()->waExpMinSumX();
    prec waXExpMinSumX = gPin->gNet()->waXExpMinSumX();

    gradientMinX = 
      ( waExpMinSumX * ( gPin->minExpSumX() * ( 1.0 - wlCoeffX * gPin->cx()) ) 
          + wlCoeffX * gPin->minExpSumX() * waXExpMinSumX )
        / ( waExpMinSumX * waExpMinSumX );
  }
  
  // max x
  if( gPin->hasMaxExpSumX() ) {
    
    prec waExpMaxSumX = gPin->gNet()->waExpMaxSumX();
    prec waXExpMaxSumX = gPin->gNet()->waXExpMaxSumX();
    
    gradientMaxX = 
      ( waExpMaxSumX * ( gPin->maxExpSumX() * ( 1.0 + wlCoeffX * gPin->cx()) ) 
          - wlCoeffX * gPin->maxExpSumX() * waXExpMaxSumX )
        / ( waExpMaxSumX * waExpMaxSumX );

  }

  // min y
  if( gPin->hasMinExpSumY() ) {
    
    prec waExpMinSumY = gPin->gNet()->waExpMinSumY();
    prec waYExpMinSumY = gPin->gNet()->waYExpMinSumY();

    gradientMinY = 
      ( waExpMinSumY * ( gPin->minExpSumY() * ( 1.0 - wlCoeffY * gPin->cy()) ) 
          + wlCoeffY * gPin->minExpSumY() * waYExpMinSumY )
        / ( waExpMinSumY * waExpMinSumY );
  }
  
  // max y
  if( gPin->hasMaxExpSumY() ) {
    
    prec waExpMaxSumY = gPin->gNet()->waExpMaxSumY();
    prec waYExpMaxSumY = gPin->gNet()->waYExpMaxSumY();
    
    gradientMaxY = 
      ( waExpMaxSumY * ( gPin->maxExpSumY() * ( 1.0 + wlCoeffY * gPin->cy()) ) 
          - wlCoeffY * gPin->maxExpSumY() * waYExpMaxSumY )
        / ( waExpMaxSumY * waExpMaxSumY );
  }

  return Point(gradientMinX - gradientMaxX, 
      gradientMinY - gradientMaxY);
}


Point NesterovBase::getWireLengthPreconditioner(GCell* gCell)
{
  return Point(gCell->gPins().size(), gCell->gPins().size());
}

Point NesterovBase::getDensityPreconditioner(GCell* gCell)
{
  prec areaVal = static_cast<prec>(gCell->dx()) 
    * static_cast<prec>(gCell->dy());

  return Point(areaVal, areaVal);
}

// get GCells' electroForcePair
// i.e. get DensityGradient with given GCell
Point NesterovBase::getDensityGradient(GCell* gCell)
{
  // find 
  BinGrid* bg = gCell->binGrid();

  auto pairX = bg->getMinMaxIdxX(gCell->lx(), gCell->ux());
  auto pairY = bg->getMinMaxIdxY(gCell->ly(), gCell->uy());
  
  Point electroForce;

  for(int i = pairX.first; i < pairX.second; i++)
  {
    for(int j = pairY.first; j < pairY.second; j++)
    {
      Bin* bin = bg->bins()[j * bg->binCntX() + i];
      prec overlapArea 
        = getOverlapDensityArea(bin, gCell) * gCell->densityScale();

      electroForce.x += overlapArea * bin->electroForceX();
      electroForce.y += overlapArea * bin->electroForceY();
    }
  }
  return electroForce;
}

Point NesterovBase::getDensityGradientWithTheta(GCell* gCell, prec& gradTheta)
{
  BinGrid* bg = gCell->binGrid();
  Point electroForce;
  gradTheta = 0.0f;
  
  auto pairX = bg->getMinMaxIdxX(gCell->lx(), gCell->ux());
  auto pairY = bg->getMinMaxIdxY(gCell->ly(), gCell->uy());
  prec bell = bellShape1(gCell->theta());
  prec derive = bellShapeDerive1(gCell->theta());
  for(int i = pairX.first; i < pairX.second; i++)
  {
    for(int j = pairY.first; j < pairY.second; j++)
    {
      Bin* bin = bg->bins()[j * bg->binCntX() + i];
      prec overlapArea = getOverlapDensityArea(bin, gCell) * gCell->densityScale();
      electroForce.x += bell * overlapArea * bin->electroForceX();
      electroForce.y += bell * overlapArea * bin->electroForceY();
      gradTheta += derive * overlapArea * bin->electroPhi();
    }
  }

  prec rlx = gCell->cx() - 0.5f * gCell->dy();
  prec rux = gCell->cx() + 0.5f * gCell->dx();
  prec rly = gCell->cy() - 0.5f * gCell->dy();
  prec ruy = gCell->cy() + 0.5f * gCell->dy();
  pairX = bg->getMinMaxIdxX(rlx, rux);
  pairY = bg->getMinMaxIdxY(rly, ruy);
  bell = bellShape2(gCell->theta());
  derive = bellShapeDerive2(gCell->theta());
  for(int i = pairX.first; i < pairX.second; i++)
  {
    for(int j = pairY.first; j < pairY.second; j++)
    {
      Bin* bin = bg->bins()[j * bg->binCntX() + i];
      prec overlapLx = std::max(rlx, bin->lx());
      prec overlapLy = std::max(rly, bin->ly());
      prec overlapUx = std::min(rux, bin->ux());
      prec overlapUy = std::min(ruy, bin->uy());
      if(overlapLx < overlapUx && overlapLy < overlapUy)
      {
        prec overlapArea = (overlapUx - overlapLx) * (overlapUy - overlapLy) * gCell->densityScale();
        electroForce.x += bell * overlapArea * bin->electroForceX();
        electroForce.y += bell * overlapArea * bin->electroForceY();
        gradTheta += derive * overlapArea * bin->electroPhi();
      }
    }
  }
  return electroForce;
}

Point NesterovBase::getDensityGradientLocal(GCell *gCell, prec alpha, prec beta, prec& cellDelta)
{
  BinGrid* bg = gCell->binGrid();
  auto pairX = bg->getMinMaxIdxX(gCell->lx(), gCell->ux());
  auto pairY = bg->getMinMaxIdxY(gCell->ly(), gCell->uy());
  Point grad;

  prec binArea = bg->binSizeX() * bg->binSizeY();
  prec invTotalCellArea = static_cast<prec>(1.0 / bg->totalCellArea());

  for(int i = pairX.first; i < pairX.second; i++)
  {
    for(int j = pairY.first; j < pairY.second; j++)
    {
      Bin* bin = bg->bins()[j * bg->binCntX() + i];
      prec binCellArea = bin->density() * binArea;
      if(binCellArea > binArea)
        cellDelta += beta * (binCellArea - binArea) * invTotalCellArea;

      prec overlapLx = std::max(bin->lx(), gCell->lx());
      prec overlapLy = std::max(bin->ly(), gCell->ly());
      prec overlapUx = std::min(bin->ux(), gCell->ux());
      prec overlapUy = std::min(bin->uy(), gCell->uy());
      prec shareArea = 0.0f, dAdx = 0.0f, dAdy = 0.0f;
      if(overlapLx < overlapUx && overlapLy < overlapUy)
      {
        shareArea = (overlapUx - overlapLx) * (overlapUy - overlapLy);
        if(bin->lx() < gCell->lx() && bin->ux() < gCell->ux())
          dAdx = -(overlapUy - overlapLy);
        else if(bin->lx() > gCell->lx() && bin->ux() > gCell->ux())
          dAdx = (overlapUy - overlapLy);
        if(bin->ly() < gCell->ly() && bin->uy() < gCell->uy())
          dAdy = -(overlapUx - overlapLx);
        else if(bin->ly() > gCell->ly() && bin->uy() > gCell->uy())
          dAdy = (overlapUx - overlapLx);
      }
      prec commonDiv = alpha / binArea;
      prec nu = fastExp(commonDiv * (binCellArea - binArea));
      prec commonVal = commonDiv * bin->electroPhi() * nu * binCellArea;
      prec commonMul = nu * shareArea;

      grad.x += dAdx * commonVal;
      grad.x += commonMul * bin->electroForceX();
      grad.y += dAdy * commonVal;
      grad.y += commonMul * bin->electroForceY();
    }
  }
  return grad;
}

Point NesterovBase::getDensityGradientLocalWithTheta(GCell *gCell, prec alpha, prec beta, prec& cellDelta, prec& gradTheta)
{
  BinGrid* bg = gCell->binGrid();
  auto pairX = bg->getMinMaxIdxX(gCell->lx(), gCell->ux());
  auto pairY = bg->getMinMaxIdxY(gCell->ly(), gCell->uy());
  Point grad;

  prec binArea = bg->binSizeX() * bg->binSizeY();
  prec invTotalCellArea = static_cast<prec>(1.0 / bg->totalCellArea());
  prec bell = bellShape1(gCell->theta());
  prec derive = bellShapeDerive1(gCell->theta());
  for(int i = pairX.first; i < pairX.second; i++)
  {
    for(int j = pairY.first; j < pairY.second; j++)
    {
      Bin* bin = bg->bins()[j * bg->binCntX() + i];
      prec binCellArea = bin->density() * binArea;
      if(binCellArea > binArea)
        cellDelta += beta * (binCellArea - binArea) * invTotalCellArea;

      prec overlapLx = std::max(bin->lx(), gCell->lx());
      prec overlapLy = std::max(bin->ly(), gCell->ly());
      prec overlapUx = std::min(bin->ux(), gCell->ux());
      prec overlapUy = std::min(bin->uy(), gCell->uy());
      prec shareArea = 0.0f, dAdx = 0.0f, dAdy = 0.0f;
      if(overlapLx < overlapUx && overlapLy < overlapUy)
      {
        shareArea = (overlapUx - overlapLx) * (overlapUy - overlapLy);
        if(bin->lx() < gCell->lx() && bin->ux() < gCell->ux())
          dAdx = -(overlapUy - overlapLy);
        else if(bin->lx() > gCell->lx() && bin->ux() > gCell->ux())
          dAdx = (overlapUy - overlapLy);
        if(bin->ly() < gCell->ly() && bin->uy() < gCell->uy())
          dAdy = -(overlapUx - overlapLx);
        else if(bin->ly() > gCell->ly() && bin->uy() > gCell->uy())
          dAdy = (overlapUx - overlapLx);
      }
      prec commonDiv = alpha / binArea;
      prec nu = fastExp(commonDiv * (binCellArea - binArea));
      prec commonVal = commonDiv * bin->electroPhi() * nu * binCellArea;
      prec commonMul = nu * shareArea;

      grad.x += bell * dAdx * commonVal;
      grad.x += bell * commonMul * bin->electroForceX();
      grad.y += bell * dAdy * commonVal;
      grad.y += bell * commonMul * bin->electroForceY();
      gradTheta += derive * nu * shareArea * bin->electroPhi();
    }
  }

  prec rlx = gCell->cx() - 0.5f * gCell->dy();
  prec rux = gCell->cx() + 0.5f * gCell->dx();
  prec rly = gCell->cy() - 0.5f * gCell->dy();
  prec ruy = gCell->cy() + 0.5f * gCell->dy();
  pairX = bg->getMinMaxIdxX(rlx, rux);
  pairY = bg->getMinMaxIdxY(rly, ruy);
  bell = bellShape2(gCell->theta());
  derive = bellShapeDerive2(gCell->theta());
  for(int i = pairX.first; i < pairX.second; i++)
  {
    for(int j = pairY.first; j < pairY.second; j++)
    {
      Bin* bin = bg->bins()[j * bg->binCntX() + i];
      prec binCellArea = bin->density() * binArea;
      if(binCellArea > binArea)
        cellDelta += beta * (binCellArea - binArea) * invTotalCellArea;

      prec overlapLx = std::max(bin->lx(), rlx);
      prec overlapLy = std::max(bin->ly(), rly);
      prec overlapUx = std::min(bin->ux(), rux);
      prec overlapUy = std::min(bin->uy(), ruy);
      prec shareArea = 0.0f, dAdx = 0.0f, dAdy = 0.0f;
      if(overlapLx < overlapUx && overlapLy < overlapUy)
      {
        shareArea = (overlapUx - overlapLx) * (overlapUy - overlapLy);
        if(bin->lx() < rlx && bin->ux() < rux)
          dAdx = -(overlapUy - overlapLy);
        else if(bin->lx() > rlx && bin->ux() > rux)
          dAdx = (overlapUy - overlapLy);
        if(bin->ly() < rly && bin->uy() < ruy)
          dAdy = -(overlapUx - overlapLx);
        else if(bin->ly() > rly && bin->uy() > ruy)
          dAdy = (overlapUx - overlapLx);
      }
      prec commonDiv = alpha / binArea;
      prec nu = fastExp(commonDiv * (binCellArea - binArea));
      prec commonVal = commonDiv * bin->electroPhi() * nu * binCellArea;
      prec commonMul = nu * shareArea;

      grad.x += bell * dAdx * commonVal;
      grad.x += bell * commonMul * bin->electroForceX();
      grad.y += bell * dAdy * commonVal;
      grad.y += bell * commonMul * bin->electroForceY();
      gradTheta += derive * nu * shareArea * bin->electroPhi();
    }
  }
  return grad;
}

// Density force cals
void NesterovBase::updateDensityForceBin()
{
  sumPhi_ = 0.0f;
  for(BinGrid* bg : binGrids_)
  {
    bg->updateDensityForceBin();
    sumPhi_ += bg->sumPhi();
  }
}

double NesterovBase::hpwl()
{
  double hpwl = 0;
  for(auto& gNet : gNets_)
  {
    gNet->updateBox();
    hpwl += gNet->hpwl();
  }
  return hpwl;
}

prec NesterovBase::overflow() const
{
  double overflowArea = 0;
  double placeStdCellArea = 0;
  double placeMacroArea = 0;
  for(BinGrid* bg : binGrids_)
  {
    overflowArea += bg->overflowArea();
    placeStdCellArea += bg->placeStdCellArea();
    placeMacroArea += bg->placeMacroArea();
  }
  return static_cast<prec>(overflowArea / (placeStdCellArea + placeMacroArea * nbVars_.targetDensity));
}

// int64_t is recommended, but prec is 2x fast
static prec getOverlapDensityArea(const  Bin* bin, const GCell* cell)
{
  prec rectLx = max(bin->lx(), cell->lx()), 
       rectLy = max(bin->ly(), cell->ly()),
       rectUx = min(bin->ux(), cell->ux()), 
       rectUy = min(bin->uy(), cell->uy());
  
  if( rectLx >= rectUx || rectLy >= rectUy )
    return 0;
  else
    return (rectUx - rectLx) * (rectUy - rectLy);
}


static prec getOverlapArea(const Bin* bin, const  Instance* inst)
{
  prec rectLx = max(bin->lx(), static_cast<prec>(inst->lx())),
       rectLy = max(bin->ly(), static_cast<prec>(inst->ly())),
       rectUx = min(bin->ux(), static_cast<prec>(inst->ux())), 
       rectUy = min(bin->uy(), static_cast<prec>(inst->uy()));

  if( rectLx >= rectUx || rectLy >= rectUy )
    return 0;
  else
    return (rectUx - rectLx) * (rectUy - rectLy);
}

// 
// https://codingforspeed.com/using-faster-exponential-approximation/
static prec fastExp(prec a)
{
  a = 1.0 + a / 1024.0;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  return a;
}

static prec bellShape(prec theta_i, prec theta)
{
  prec delta_abs = (1.f / (2.f * PI)) * std::abs(theta_i - theta);
  if (delta_abs <= 0.125f)
    return 1.f - 32.f * delta_abs * delta_abs;
  else if (delta_abs <= 0.25f)
    return 32.f * (delta_abs - 0.25f) * (delta_abs - 0.25f);
  else
    return 0.f;
}

static prec bellShape1(prec theta_i)
{
  return bellShape(theta_i, 0.0f) + bellShape(theta_i, PI) + bellShape(theta_i, 2.0f * PI);
}

static prec bellShape2(prec theta_i)
{
  return bellShape(theta_i, 0.5f * PI) + bellShape(theta_i, 1.5f * PI);
}

static prec bellShapeDerive(prec theta_i, prec theta)
{
  prec delta = (1.f / (2.f * PI)) * (theta_i - theta);
  prec delta_abs = std::abs(delta);
  if (delta_abs <= 0.125f)
    return (-16.f / PI) * delta;
  else if (delta_abs <= 0.25f)
    return (32.f / PI) * delta + ((delta > 0) ? (-8.f / PI) : (8.f / PI));
  else
    return 0.f;
}

static prec bellShapeDerive1(prec theta_i)
{
  return bellShapeDerive(theta_i, 0.0f) + bellShapeDerive(theta_i, PI) + bellShapeDerive(theta_i, 2.0f * PI);
}

static prec bellShapeDerive2(prec theta_i)
{
  return bellShapeDerive(theta_i, 0.5f * PI) + bellShapeDerive(theta_i, 1.5f * PI);
}

}
