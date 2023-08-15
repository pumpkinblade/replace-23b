#include "nesterovBase.h"
#include "placerBase.h"
#include "fft.h"
#include "log.h"

#include <iostream>
#include <random>
#include <algorithm>


#define REPLACE_SQRT2 1.414213562373095048801L

namespace replace {

using namespace std;

static int 
fastModulo(const int input, const int ceil);

static int64_t
getOverlapArea(Bin* bin, Instance* inst);

// Note that
// int64_t is ideal in the following function, but
// runtime is doubled compared with prec.
//
// Choose to use "prec" only in the following functions
static int64_t 
getOverlapDensityArea(Bin* bin, GCell* cell);

static prec
fastExp(prec exp);


////////////////////////////////////////////////
// GCell 

GCell::GCell() 
    : inst_(nullptr), isMacro_(false),
      lx_(0), ly_(0), ux_(0), uy_(0),
      dLx_(0), dLy_(0), dUx_(0), dUy_(0),
      densityScale_(0), gradientX_(0), gradientY_(0)
{
}

GCell::GCell(Instance* inst) 
    : GCell()
{
  setInstance(inst);
}

GCell::GCell(int cx, int cy, int dx, int dy) 
    : GCell() 
{
  dLx_ = lx_ = cx - dx/2;
  dLy_ = ly_ = cy - dy/2;
  dUx_ = ux_ = cx + dx/2;
  dUy_ = uy_ = cy + dy/2; 
}

void GCell::setInstance(Instance* inst)
{
  inst_ = inst;
  // density coordi has the same center points.
  dLx_ = lx_ = inst->lx();
  dLy_ = ly_ = inst->ly();
  dUx_ = ux_ = inst->ux();
  dUy_ = uy_ = inst->uy();
  isMacro_ = inst->isMacro();
}

void GCell::addGPin(GPin* gPin)
{
  gPins_.push_back(gPin);
}

void GCell::setLocation(int lx, int ly)
{
  ux_ = lx + (ux_ - lx_);
  uy_ = ly + (uy_ - ly_);
  lx = lx_;
  ly = ly_;

  for(auto& gPin: gPins_)
  {
    gPin->updateLocation(this);
  }
}

void GCell::setCenterLocation(int cx, int cy)
{
  const int halfDx = dx()/2;
  const int halfDy = dy()/2;

  lx_ = cx - halfDx;
  ly_ = cy - halfDy;
  ux_ = cx + halfDx;
  uy_ = cy + halfDy;

  for(auto& gPin: gPins_) {
    gPin->updateLocation(this);
  }
}

// changing size and preserve center coordinates
void GCell::setSize(int dx, int dy)
{
  const int centerX = cx();
  const int centerY = cy();

  lx_ = centerX - dx/2;
  ly_ = centerY - dy/2;
  ux_ = centerX + dx/2;
  uy_ = centerY + dy/2;
}

void GCell::setDensityLocation(int dLx, int dLy)
{
  dUx_ = dLx + (dUx_ - dLx_);
  dUy_ = dLy + (dUy_ - dLy_);
  dLx_ = dLx;
  dLy_ = dLy;

  // assume that density Center change the gPin coordi
  for(auto& gPin: gPins_) {
    gPin->updateDensityLocation(this);
  }
}

void GCell::setDensityCenterLocation(int dCx, int dCy)
{
  const int halfDDx = dDx()/2;
  const int halfDDy = dDy()/2;

  dLx_ = dCx - halfDDx;
  dLy_ = dCy - halfDDy;
  dUx_ = dCx + halfDDx;
  dUy_ = dCy + halfDDy;
  

  // assume that density Center change the gPin coordi
  for(auto& gPin: gPins_) {
    gPin->updateDensityLocation(this);
  }
}

// changing size and preserve center coordinates
void GCell::setDensitySize(int dDx, int dDy)
{
  const int dCenterX = dCx();
  const int dCenterY = dCy();

  dLx_ = dCenterX - dDx/2;
  dLy_ = dCenterY - dDy/2;
  dUx_ = dCenterX + dDx/2;
  dUy_ = dCenterY + dDy/2;
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
  lx_ = ly_ = INT_MAX;
  ux_ = uy_ = INT_MIN;

  for(auto& gPin : gPins_) {
    lx_ = std::min(gPin->cx(), lx_);
    ly_ = std::min(gPin->cy(), ly_);
    ux_ = std::max(gPin->cx(), ux_);
    uy_ = std::max(gPin->cy(), uy_);
  } 
}

int64_t GNet::hpwl() const
{
  return (int64_t)(ux_ - lx_) + (uy_ - ly_);
}

void GNet::clearWaVars() {
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
  cx_ = pin->cx();
  cy_ = pin->cy();
  offsetCx_ = pin->offsetCx();
  offsetCy_ = pin->offsetCy();
}

void GPin::setCenterLocation(int cx, int cy)
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

void GPin::updateDensityLocation(const GCell* gCell)
{
  cx_ = gCell->dCx() + offsetCx_;
  cy_ = gCell->dCy() + offsetCy_;
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

Bin::Bin(int x, int y, int lx, int ly, int ux, int uy, prec targetDensity) 
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
      isSetBinCntX_(0), isSetBinCntY_(0)
{
}

BinGrid::BinGrid(Die* die)
   : BinGrid()
{
  setDie(die);
}

BinGrid::BinGrid(BinGrid&& other)
{
  *this = std::move(other);
}

BinGrid& BinGrid::operator=(BinGrid&& other)
{
  std::swap(binStor_, other.binStor_);
  std::swap(bins_, other.bins_);
  std::swap(lx_, other.lx_);
  std::swap(ly_, other.ly_);
  std::swap(ux_, other.ux_);
  std::swap(uy_, other.uy_);
  std::swap(binCntX_, other.binCntX_);
  std::swap(binCntY_, other.binCntY_);
  std::swap(binSizeX_, other.binSizeX_);
  std::swap(binSizeY_, other.binSizeY_);
  std::swap(targetDensity_, other.targetDensity_);
  std::swap(overflowArea_, other.overflowArea_);
  std::swap(isSetBinCntX_, other.isSetBinCntX_);
  std::swap(isSetBinCntY_, other.isSetBinCntY_);
  std::swap(sumPhi_, other.sumPhi_);
  std::swap(ly_, other.ly_);
  std::swap(die_, other.die_);
  std::swap(gCells_, other.gCells_);
  std::swap(fft_, other.fft_);

  return *this;
}

void BinGrid::setDie(Die* die)
{
  die_ = die;
  lx_ = die->coreLx();
  ly_ = die->coreLy();
  ux_ = die->coreUx();
  uy_ = die->coreUy();

  placeInstCnt_ = 0;
  placeStdCellArea_ = 0;
  placeMacroArea_ = 0;
  fixedInstArea_ = 0;
  for(Instance* inst : die_->insts())
  {
    if(inst->isFixed())
    {
      fixedInstArea_ += (int64_t)(inst->dx()) * inst->dy();
    }
    else
    {
      placeInstCnt_++;
      if(inst->isMacro())
        placeMacroArea_ += (int64_t)(inst->dx()) * inst->dy();
      else
        placeStdCellArea_ += (int64_t)(inst->dx()) * inst->dy();
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
  int64_t totalBinArea 
    = static_cast<int64_t>(ux_ - lx_) 
    * static_cast<int64_t>(uy_ - ly_);

  int foundBinCnt = 2;
  if (placeInstCnt_ != 0)
  {
    // assume no fixed instance
    int64_t avgPlaceInstArea = (placeMacroArea_ + placeStdCellArea_) / placeInstCnt_;

    int64_t idealBinArea = std::round(static_cast<prec>(avgPlaceInstArea) / targetDensity_);
    int idealBinCnt = totalBinArea / idealBinArea; 
    
    LOG_DEBUG("TargetDensity: {}", targetDensity_);
    LOG_DEBUG("AveragePlaceInstArea: {}", avgPlaceInstArea);
    LOG_DEBUG("IdealBinArea: {}", idealBinArea);
    LOG_DEBUG("IdealBinCnt: {}", idealBinCnt);
    LOG_DEBUG("TotalBinArea: {}", totalBinArea);

    // find binCnt: 2, 4, 8, 16, 32, 64, ...
    // s.t. binCnt^2 <= idealBinCnt <= (binCnt*2)^2.
    for(foundBinCnt = 2; foundBinCnt <= 1024; foundBinCnt *= 2)
    {
      if( foundBinCnt * foundBinCnt <= idealBinCnt 
          && 4 * foundBinCnt * foundBinCnt > idealBinCnt )
      {
        break;
      }
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
  
  binSizeX_ = (ux_ - lx_) / binCntX_;
  binSizeY_ = (uy_ - ly_) / binCntY_;
  
  LOG_DEBUG("BinSize: {}, {}", binSizeX_, binSizeY_);

  // initialize binStor_, bins_ vector
  binStor_.reserve(binCntX_ * binCntY_);
  bins_.reserve(binCntX_ * binCntY_);
  int x = lx_, y = ly_;
  int idxX = 0, idxY = 0;
  for (idxY = 0; idxY < binCntY_; idxY++)
  {
    int sizeY = (idxY == binCntY_ - 1) ? (uy_ - y) : binSizeY_;
    x = lx_;
    for (idxX = 0; idxX < binCntX_; idxX++)
    {
      int sizeX = (idxX == binCntX_ - 1) ? (ux_ - x) : binSizeX_;
      binStor_.emplace_back(idxX, idxY, x, y, x + sizeX, y + sizeY, targetDensity_);
      bins_.push_back(&binStor_.back());
      x += sizeX;
    }
    y += sizeY;
  }
  //for(auto& bin : binStor_)
  //{

  //  int sizeX = (x + binSizeX_ > ux_) ? ux_ - x : binSizeX_;
  //  int sizeY = (y + binSizeY_ > uy_) ? uy_ - y : binSizeY_;

  //  //cout << "idxX: " << idxX << " idxY: " << idxY 
  //  //  << " x:" << x << " y:" << y 
  //  //  << " " << x+sizeX << " " << y+sizeY << endl;
  //  bin = Bin(idxX, idxY, x, y, x+sizeX, y+sizeY, targetDensity_);
  //  
  //  // move x, y coordinates.
  //  x += binSizeX_;
  //  idxX += 1;

  //  if( x >= ux_ )
  //  {
  //    y += binSizeY_;
  //    x = lx_; 
  //    
  //    idxY ++;
  //    idxX = 0;
  //  }

  //  bins_.push_back( &bin );
  //}

  LOG_DEBUG("NumBins: {}", bins_.size());

  // only initialized once
  updateBinsNonPlaceArea();

  // initialize fft structrue based on bins
  std::unique_ptr<FFT> fft(new FFT(binCntX_, binCntY_, binSizeX_, binSizeY_));
  fft_ = std::move(fft);
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
      std::pair<int, int> pairX = getMinMaxIdxX(inst);
      std::pair<int, int> pairY = getMinMaxIdxY(inst);
      for(int i = pairX.first; i < pairX.second; i++)
      {
        for(int j = pairY.first; j < pairY.second; j++)
        {
          Bin* bin = bins_[ j * binCntX_ + i ];

          // Note that nonPlaceArea should have scale-down with
          // target density. 
          // See MS-replace paper
          //
          bin->addNonPlaceArea(
            getOverlapArea(bin, inst) * bin->targetDensity());
        }
      }
    }
  }
}


// Core Part
void BinGrid::updateBinsGCellDensityArea()
{
  // clear the Bin-area info
  for(auto& bin : bins_) {
    bin->setInstPlacedArea(0);
    bin->setFillerArea(0);
  }

  for(auto& cell : gCells_) {
    std::pair<int, int> pairX 
      = getDensityMinMaxIdxX(cell);
    std::pair<int, int> pairY 
      = getDensityMinMaxIdxY(cell);

    // The following function is critical runtime hotspot 
    // for global placer.
    //
    if(cell->isInstance())
    {
      // macro should have 
      // scale-down with target-density
      if(cell->isMacro())
      {
        for(int i = pairX.first; i < pairX.second; i++)
        {
          for(int j = pairY.first; j < pairY.second; j++)
          {
            Bin* bin = bins_[j * binCntX_ + i];
            int64_t overlapArea = getOverlapDensityArea(bin, cell);
            overlapArea = static_cast<int64_t>(overlapArea * cell->densityScale() * bin->targetDensity());
            bin->addInstPlacedArea(overlapArea);
          }
        }
      }
      // normal cells
      else
      {
        for(int i = pairX.first; i < pairX.second; i++)
        {
          for(int j = pairY.first; j < pairY.second; j++)
          {
            Bin* bin = bins_[j * binCntX_ + i];
            prec overlapArea = getOverlapDensityArea(bin, cell);
            overlapArea = static_cast<int64_t>(overlapArea * cell->densityScale());
            bin->addInstPlacedArea(overlapArea); 
          }
        }
      }
    }
    else
    {
      for(int i = pairX.first; i < pairX.second; i++)
      {
        for(int j = pairY.first; j < pairY.second; j++)
        {
          Bin* bin = bins_[j * binCntX_ + i];
          prec overlapArea = getOverlapDensityArea(bin, cell);
          overlapArea = static_cast<int64_t>(overlapArea * cell->densityScale());
          bin->addFillerArea(overlapArea);
        }
      }
    }
  }

  overflowArea_ = 0;

  // update density and overflowArea 
  // for nesterov use and FFT library
  for(auto& bin : bins_)
  {
    prec binArea = static_cast<prec>(bin->binArea());
    prec cellArea = static_cast<prec>(bin->instPlacedArea() + bin->fillerArea() + bin->nonPlaceArea());
    prec density = cellArea / (binArea * bin->targetDensity());
    bin->setDensity(density);
    if (std::isnan(density))
    {
      LOG_ERROR("NaN density at bin {} {}", bin->x(), bin->y());
      LOG_INFO("binArea : {}, cellArea: {}, targetDensity: {}", binArea, cellArea, bin->targetDensity());
      LOG_INFO("bin size x: {}, bin size y: {}", bin->dx(), bin->dy());
      LOG_INFO("instPlaceArea: {}, fillerArea: {}, nonPlaceArea: {}", bin->instPlacedArea(), bin->fillerArea(), bin->nonPlaceArea());
    }
    else if (std::isinf(density))
    {
      LOG_ERROR("INF density at bin {} {}", bin->x(), bin->y());
      LOG_INFO("binAarea : {}, cellArea: {}, targetDensity: {}", binArea, cellArea, bin->targetDensity());
      LOG_INFO("binArea : {}, cellArea: {}, targetDensity: {}", binArea, cellArea, bin->targetDensity());
      LOG_INFO("bin size x: {}, bin size y: {}", bin->dx(), bin->dy());
      LOG_INFO("instPlaceArea: {}, fillerArea: {}, nonPlaceArea: {}", bin->instPlacedArea(), bin->fillerArea(), bin->nonPlaceArea());
    }

    overflowArea_ 
      += std::max((prec)0,
          static_cast<prec>(bin->instPlacedArea() + bin->nonPlaceArea())
          - (static_cast<prec>(binArea) * bin->targetDensity()));
  }
}


std::pair<int, int> BinGrid::getDensityMinMaxIdxX(GCell* gcell)
{
  int lowerIdx = (gcell->dLx() - lx())/binSizeX_;
  int upperIdx = 
    (fastModulo((gcell->dUx() - lx()), binSizeX_) == 0)
    ? (gcell->dUx() - lx()) / binSizeX_ 
    : (gcell->dUx() - lx()) / binSizeX_ + 1;
  //int lowerIdx = (gcell->dLx() - lx()) / binSizeX_;
  //int upperIdx = (gcell->dUx() - lx() + binSizeX_ - 1) / binSizeX_;
  return std::make_pair(std::max(lowerIdx, 0), std::min(upperIdx, binCntX_));
}

std::pair<int, int> BinGrid::getDensityMinMaxIdxY(GCell* gcell)
{
  int lowerIdx = (gcell->dLy() - ly())/binSizeY_;
  int upperIdx =
    (fastModulo((gcell->dUy() - ly()), binSizeY_) == 0)
    ? (gcell->dUy() - ly()) / binSizeY_ 
    : (gcell->dUy() - ly()) / binSizeY_ + 1;
  //int lowerIdx = (gcell->dLy() - ly()) / binSizeY_;
  //int upperIdx = (gcell->dUy() - ly() + binSizeY_ - 1) / binSizeY_;
  return std::make_pair(std::max(lowerIdx, 0), std::min(upperIdx, binCntY_));
}

std::pair<int, int> BinGrid::getMinMaxIdxX(Instance* inst)
{
  int lowerIdx = (inst->lx() - lx()) / binSizeX_;
  int upperIdx = 
    (fastModulo((inst->ux() - lx()), binSizeX_) == 0)
    ? (inst->ux() - lx()) / binSizeX_ 
    : (inst->ux() - lx()) / binSizeX_ + 1;
  //int lowerIdx = (inst->lx() - lx()) / binSizeX_;
  //int upperIdx = (inst->ux() - lx() + binSizeX_ - 1) / binSizeX_;
  return std::make_pair(std::max(lowerIdx, 0), std::min(upperIdx, binCntX_));
}

std::pair<int, int> BinGrid::getMinMaxIdxY(Instance* inst)
{
  int lowerIdx = (inst->ly() - ly()) / binSizeY_;
  int upperIdx = 
    (fastModulo((inst->uy() - ly()), binSizeY_) == 0)
    ? (inst->uy() - ly()) / binSizeY_ 
    : (inst->uy() - ly()) / binSizeY_ + 1;
  //int lowerIdx = (inst->ly() - ly()) / binSizeY_;
  //int upperIdx = (inst->uy() - ly() + binSizeY_ - 1) / binSizeY_;
  return std::make_pair(std::max(lowerIdx, 0), std::min(upperIdx, binCntY_));
}

void BinGrid::addGCell(GCell* gc)
{
  gCells_.push_back(gc);
}

void BinGrid::updateDensityForceBin()
{
  // copy density to utilize FFT
  for(Bin* bin : bins_)
  {
    fft_->updateDensity(
      bin->x(), bin->y(), bin->density());  
  }

  // do FFT
  fft_->doFFT();

  // update electroPhi and electroForce
  // update sumPhi_ for nesterov loop
  sumPhi_ = 0;
  for(Bin* bin : bins_)
  {
    auto eForcePair = fft_->getElectroForce(bin->x(), bin->y());
    bin->setElectroForceX(eForcePair.first);
    bin->setElectroForceY(eForcePair.second);

    prec electroPhi = fft_->getElectroPhi(bin->x(), bin->y());
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
      scaleX = static_cast<prec>(gCell->dx()) 
        / static_cast<prec>(REPLACE_SQRT2 * binSizeX_);
      densitySizeX = REPLACE_SQRT2 
        * static_cast<prec>(binSizeX_);
    }
    else
    {
      scaleX = 1.0f;
      densitySizeX = gCell->dx();
    }

    if(gCell->dy() < REPLACE_SQRT2 * binSizeY_)
    {
      scaleY = static_cast<prec>(gCell->dy()) 
               / static_cast<prec>(REPLACE_SQRT2 * binSizeY_);
      densitySizeY = REPLACE_SQRT2 * static_cast<prec>(binSizeY_);
    }
    else
    {
      scaleY = 1.0f;
      densitySizeY = gCell->dy();
    }

    gCell->setDensitySize(densitySizeX, densitySizeY);
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
    //binGridMap_.emplace(bg.die(), &bg);
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

  // // extract average dx/dy in range (10%, 90%)
  // vector<int> dxStor;
  // vector<int> dyStor;
  // dxStor.reserve(die->placeInsts().size());
  // dyStor.reserve(die->placeInsts().size());
  // for(Instance* inst : die->placeInsts())
  // {
  //   dxStor.push_back(inst->dx());
  //   dyStor.push_back(inst->dy());
  // }
  
  // // sort
  // std::sort(dxStor.begin(), dxStor.end());
  // std::sort(dyStor.begin(), dyStor.end());

  // // average from (10 - 90%) .
  // int64_t dxSum = 0, dySum = 0;
  // int minIdx = static_cast<int>(dxStor.size()*0.05);
  // int maxIdx = static_cast<int>(dxStor.size()*0.95);
  // assert(maxIdx - minIdx > 0);
  // for(int i=minIdx; i<maxIdx; i++)
  // {
  //   dxSum += dxStor[i];
  //   dySum += dyStor[i];
  // }

  // // the avgDx and avgDy will be used as filler cells' 
  // // width and height
  // int avgDx = static_cast<int>(dxSum / (maxIdx - minIdx));
  // int avgDy = static_cast<int>(dySum / (maxIdx - minIdx));

  double dxSum = 0.0, dySum = 0.0;
  int numStdcell = 0;
  for(Instance* inst : bg->die()->insts())
  {
    if(!inst->isFixed() && !inst->isMacro())
    {
      dxSum += inst->dx();
      dySum += inst->dy();
      numStdcell++;
    }
  }
  int avgDx = static_cast<int>(dxSum / numStdcell);
  int avgDy = static_cast<int>(dySum / numStdcell);

  int64_t coreArea = (int64_t)bg->die()->coreDx() * bg->die()->coreDy();

  // nonPlaceInstsArea should not have targetDensity downscaling!!! 
  int64_t whiteSpaceArea = coreArea - bg->fixedInstanceArea();

  // TODO density screening
  int64_t movableArea = static_cast<int64_t>(whiteSpaceArea * nbVars_.targetDensity);
  
  int64_t totalFillerArea = movableArea - bg->placeStdCellArea() 
                          - static_cast<int64_t>(bg->placeMacroArea() * nbVars_.targetDensity);

  if(totalFillerArea < 0)
  {
    LOG_ERROR("Filler area is negative!!\n"
              "\tPlease put higher target density or\n"
              "\tRe-floorplan to have enough coreArea");
  }

  assert(avgDx> 0);
  assert(avgDy> 0);
  int fillerCnt = static_cast<int>(totalFillerArea / ((int64_t)avgDx * avgDy));

  LOG_DEBUG("FillerInit: CoreArea: {}", coreArea);
  LOG_DEBUG("FillerInit: WhiteSpaceArea: {}", whiteSpaceArea);
  LOG_DEBUG("FillerInit: MovableArea: {}", movableArea);
  LOG_DEBUG("FillerInit: TotalFillerArea: {}", totalFillerArea);
  LOG_DEBUG("FillerInit: NumFillerCells: {}", fillerCnt);
  LOG_DEBUG("FillerInit: FillerCellArea: {}", static_cast<int64_t>(avgDx*avgDy));
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
    GCell myGCell(
        randX % bg->die()->coreDx() + bg->die()->coreLx(), 
        randY % bg->die()->coreDy() + bg->die()->coreLy(),
        avgDx, avgDy);
    gCellStor_.push_back(myGCell);
  }
}

//GCell* NesterovBase::placerToNesterov(Instance* inst)
//{
//  auto gcPtr = gCellMap_.find(inst);
//  auto result = (gcPtr == gCellMap_.end()) ? nullptr : gcPtr->second;
//  // assert(result != nullptr);
//  return result;
//}
//
//GNet* NesterovBase::placerToNesterov(Net* net)
//{
//  auto gnPtr = gNetMap_.find(net);
//  auto res = (gnPtr == gNetMap_.end()) ? nullptr : gnPtr->second;
//  // assert(res != nullptr);
//  return res;
//}
//
//GPin* NesterovBase::placerToNesterov(Pin* pin)
//{
//  auto gpPtr = gPinMap_.find(pin);
//  return (gpPtr == gPinMap_.end()) ? nullptr : gpPtr->second;
//}
//
//BinGrid* NesterovBase::placerToNesterov(Die* die)
//{
//  auto bgPtr = binGridMap_.find(die);
//  return (bgPtr == binGridMap_.end()) ? nullptr : bgPtr->second;
//}

// gcell update
void NesterovBase::updateGCellLocation(
  std::vector<Point>& coordis)
{
  for(auto& coordi : coordis)
  {
    int idx = &coordi - &coordis[0];
    gCells_[idx]->setLocation( coordi.x, coordi.y );
  }
}

// gcell update
void NesterovBase::updateGCellCenterLocation(
    std::vector<Point>& coordis)
{
  for(auto& coordi : coordis)
  {
    int idx = &coordi - &coordis[0];
    gCells_[idx]->setCenterLocation( coordi.x, coordi.y );
  }
}

void NesterovBase::updateGCellDensityCenterLocation(
    std::vector<Point>& coordis)
{
  for(auto& coordi : coordis)
  {
    int idx = &coordi - &coordis[0];
    gCells_[idx]->setDensityCenterLocation( 
        coordi.x, coordi.y );
  }
  for(BinGrid* bg : binGrids_)
  {
    bg->updateBinsGCellDensityArea();
  }
}

int64_t  NesterovBase::overflowArea() const {
  int64_t area = 0;
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

void  NesterovBase::updateDensityCoordiLayoutInside(GCell* gCell)
{
  prec targetLx = gCell->dLx();
  prec targetLy = gCell->dLy();

  BinGrid* bg = gCell->binGrid();
  if( targetLx < bg->lx() )
  {
    targetLx = bg->lx();
  }
  if(targetLy < bg->ly())
  {
    targetLy = bg->ly();
  }
  if(targetLx + gCell->dDx() > bg->ux())
  {
    targetLx = bg->ux() - gCell->dDx();
  }
  if(targetLy + gCell->dDy() > bg->uy())
  {
    targetLy = bg->uy() - gCell->dDy();
  }
  gCell->setDensityLocation(targetLx, targetLy);
}

prec NesterovBase::getDensityCoordiLayoutInsideX(
  GCell* gCell, prec cx)
{
  prec adjVal = cx;
  // TODO will change base on each assigned binGrids.
  // 
  BinGrid* bg = gCell->binGrid();
  if(cx - gCell->dDx()/2 < bg->lx())
  {
    adjVal = bg->lx() + gCell->dDx() / 2;
  }
  if(cx + gCell->dDx()/2 > bg->ux())
  {
    adjVal = bg->ux() - gCell->dDx() / 2;
  }
  return adjVal;
}

prec NesterovBase::getDensityCoordiLayoutInsideY(
  GCell* gCell, prec cy)
{
  prec adjVal = cy;
  // TODO will change base on each assigned binGrids.
  // 
  BinGrid* bg = gCell->binGrid();
  if(cy - gCell->dDy()/2 < bg->ly())
  {
    adjVal = bg->ly() + gCell->dDy()/2;
  }
  if(cy + gCell->dDy()/2 > bg->uy())
  {
    adjVal = bg->uy() - gCell->dDy()/2; 
  }

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
  return Point(gCell->gPins().size(),
     gCell->gPins().size());
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

  auto pairX = bg->getDensityMinMaxIdxX(gCell);
  auto pairY = bg->getDensityMinMaxIdxY(gCell);
  
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

int64_t NesterovBase::hpwl()
{
  int64_t hpwl = 0;
  for(auto& gNet : gNets_)
  {
    gNet->updateBox();
    hpwl += gNet->hpwl();
  }
  return hpwl;
}

prec NesterovBase::overflow() const
{
  int64_t overflowArea = 0;
  int64_t placeStdCellArea = 0;
  int64_t placeMacroArea = 0;
  for(BinGrid* bg : binGrids_)
  {
    overflowArea += bg->overflowArea();
    placeStdCellArea += bg->placeStdCellArea();
    placeMacroArea += bg->placeMacroArea();
  }
  return static_cast<prec>(overflowArea) 
         / static_cast<prec>(placeStdCellArea 
                              + placeMacroArea * nbVars_.targetDensity);
}


// https://stackoverflow.com/questions/33333363/built-in-mod-vs-custom-mod-function-improve-the-performance-of-modulus-op
static int  fastModulo(const int input, const int ceil)
{
  return input >= ceil? input % ceil : input;
}

// int64_t is recommended, but prec is 2x fast
static int64_t getOverlapDensityArea(Bin* bin, GCell* cell)
{
  int rectLx = max(bin->lx(), cell->dLx()), 
      rectLy = max(bin->ly(), cell->dLy()),
      rectUx = min(bin->ux(), cell->dUx()), 
      rectUy = min(bin->uy(), cell->dUy());
  
  if( rectLx >= rectUx || rectLy >= rectUy ) {
    return 0;
  }
  else {
    //return static_cast<prec>(rectUx - rectLx) 
    //  * static_cast<prec>(rectUy - rectLy);
    return (int64_t)(rectUx - rectLx) * (rectUy - rectLy);
  }
}


static int64_t getOverlapArea(Bin* bin, Instance* inst)
{
  int rectLx = max(bin->lx(), inst->lx()), 
      rectLy = max(bin->ly(), inst->ly()),
      rectUx = min(bin->ux(), inst->ux()), 
      rectUy = min(bin->uy(), inst->uy());

  if( rectLx >= rectUx || rectLy >= rectUy ) {
    return 0;
  }
  else {
    return static_cast<int64_t>(rectUx - rectLx) 
      * static_cast<int64_t>(rectUy - rectLy);
  }
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

}
