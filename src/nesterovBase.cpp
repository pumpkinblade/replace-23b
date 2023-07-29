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
// runtime is doubled compared with float.
//
// Choose to use "float" only in the following functions
static float 
getOverlapDensityArea(Bin* bin, GCell* cell);

static float
fastExp(float exp);


////////////////////////////////////////////////
// GCell 

GCell::GCell() 
  : lx_(0), ly_(0), ux_(0), uy_(0),
  dLx_(0), dLy_(0), dUx_(0), dUy_(0),
  densityScale_(0), 
  gradientX_(0), gradientY_(0) {}


GCell::GCell(Instance* inst) 
  : GCell() {
  setInstance(inst);
}

GCell::GCell(std::vector<Instance*>& insts) 
  : GCell() {
  setClusteredInstance(insts);
}

GCell::GCell(int cx, int cy, int dx, int dy) 
  : GCell() {
  dLx_ = lx_ = cx - dx/2;
  dLy_ = ly_ = cy - dy/2;
  dUx_ = ux_ = cx + dx/2;
  dUy_ = uy_ = cy + dy/2; 
  setFiller();
}

GCell::~GCell() {
  insts_.clear();
}

void
GCell::setInstance(Instance* inst) {
  insts_.push_back(inst);
  // density coordi has the same center points.
  dLx_ = lx_ = inst->lx();
  dLy_ = ly_ = inst->ly();
  dUx_ = ux_ = inst->ux();
  dUy_ = uy_ = inst->uy();
}

Instance*
GCell::instance() const {
  return *insts_.begin();
}

void
GCell::addGPin(GPin* gPin) {
  gPins_.push_back(gPin);
}


// do nothing
void
GCell::setFiller() {
}

void
GCell::setClusteredInstance(std::vector<Instance*>& insts) {
  insts_ = insts;
}

void
GCell::setMacroInstance() {
  isMacroInstance_ = true;
}

void
GCell::setStdInstance() {
  isMacroInstance_ = false; 
}

void
GCell::setLocation(int lx, int ly) {
  ux_ = lx + (ux_ - lx_);
  uy_ = ly + (uy_ - ly_);
  lx = lx_;
  ly = ly_;

  for(auto& gPin: gPins_) {
    gPin->updateLocation(this);
  }
}

void
GCell::setCenterLocation(int cx, int cy) {
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
void
GCell::setSize(int dx, int dy) {
  const int centerX = cx();
  const int centerY = cy();

  lx_ = centerX - dx/2;
  ly_ = centerY - dy/2;
  ux_ = centerX + dx/2;
  uy_ = centerY + dy/2;
}

void
GCell::setDensityLocation(int dLx, int dLy) {
  dUx_ = dLx + (dUx_ - dLx_);
  dUy_ = dLy + (dUy_ - dLy_);
  dLx_ = dLx;
  dLy_ = dLy;

  // assume that density Center change the gPin coordi
  for(auto& gPin: gPins_) {
    gPin->updateDensityLocation(this);
  }
}

void
GCell::setDensityCenterLocation(int dCx, int dCy) {
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
void
GCell::setDensitySize(int dDx, int dDy) {
  const int dCenterX = dCx();
  const int dCenterY = dCy();

  dLx_ = dCenterX - dDx/2;
  dLy_ = dCenterY - dDy/2;
  dUx_ = dCenterX + dDx/2;
  dUy_ = dCenterY + dDy/2;
}

void
GCell::setDensityScale(float densityScale) {
  densityScale_ = densityScale;
}

void
GCell::setGradientX(float gradientX) {
  gradientX_ = gradientX;
}

void
GCell::setGradientY(float gradientY) {
  gradientY_ = gradientY;
}

bool 
GCell::isInstance() const {
  return (insts_.size() == 1);
}

bool
GCell::isClusteredInstance() const {
  return (insts_.size() > 0);
}

bool
GCell::isFiller() const {
  return (insts_.size() == 0);
}

bool
GCell::isMacroInstance() const {
  return isMacroInstance_; 
}

bool
GCell::isStdInstance() const {
  return !isMacroInstance_;
}

////////////////////////////////////////////////
// GNet

GNet::GNet()
  : lx_(0), ly_(0), ux_(0), uy_(0),
  customWeight_(1), weight_(1),
  waExpMinSumX_(0), waXExpMinSumX_(0),
  waExpMaxSumX_(0), waXExpMaxSumX_(0),
  waExpMinSumY_(0), waYExpMinSumY_(0),
  waExpMaxSumY_(0), waYExpMaxSumY_(0),
  isDontCare_(0) {}

GNet::GNet(Net* net) : GNet() {
  nets_.push_back(net);
}

GNet::GNet(std::vector<Net*>& nets) : GNet() {
  nets_ = nets;
}

GNet::~GNet() {
  gPins_.clear();
  nets_.clear();
}

Net* 
GNet::net() const { 
  return *nets_.begin();
}

void
GNet::setCustomWeight(float customWeight) {
  customWeight_ = customWeight; 
}

void
GNet::addGPin(GPin* gPin) {
  gPins_.push_back(gPin);
}

void
GNet::updateBox() {
  lx_ = ly_ = INT_MAX;
  ux_ = uy_ = INT_MIN;

  for(auto& gPin : gPins_) {
    lx_ = std::min(gPin->cx(), lx_);
    ly_ = std::min(gPin->cy(), ly_);
    ux_ = std::max(gPin->cx(), ux_);
    uy_ = std::max(gPin->cy(), uy_);
  } 
}

int64_t
GNet::hpwl() {
  return static_cast<int64_t>((ux_ - lx_) + (uy_ - ly_));
}

void
GNet::clearWaVars() {
  waExpMinSumX_ = 0;
  waXExpMinSumX_ = 0;

  waExpMaxSumX_ = 0;
  waXExpMaxSumX_ = 0;
    
  waExpMinSumY_ = 0;
  waYExpMinSumY_ = 0;

  waExpMaxSumY_ = 0;
  waYExpMaxSumY_ = 0;
}


void
GNet::setDontCare() {
  isDontCare_ = 1;
}

bool
GNet::isDontCare() {
  return (gPins_.size() == 0) || (isDontCare_ == 1);
}

////////////////////////////////////////////////
// GPin 

GPin::GPin()
  : gCell_(nullptr), gNet_(nullptr),
  offsetCx_(0), offsetCy_(0),
  cx_(0), cy_(0),
  maxExpSumX_(0), maxExpSumY_(0),
  minExpSumX_(0), minExpSumY_(0),
  hasMaxExpSumX_(0), hasMaxExpSumY_(0), 
  hasMinExpSumX_(0), hasMinExpSumY_(0) {}

GPin::GPin(Pin* pin)
  : GPin() {
  pins_.push_back(pin);
  cx_ = pin->cx();
  cy_ = pin->cy();
  offsetCx_ = pin->offsetCx();
  offsetCy_ = pin->offsetCy();
}

GPin::GPin(std::vector<Pin*> & pins) {
  pins_ = pins;
}

GPin::~GPin() {
  gCell_ = nullptr;
  gNet_ = nullptr;
  pins_.clear();
}

Pin* 
GPin::pin() const {
  return *pins_.begin();
}

void
GPin::setGCell(GCell* gCell) {
  gCell_ = gCell;
}

void
GPin::setGNet(GNet* gNet) {
  gNet_ = gNet;
}

void
GPin::setCenterLocation(int cx, int cy) {
  cx_ = cx;
  cy_ = cy;
}

void
GPin::clearWaVars() {
  hasMaxExpSumX_ = 0;
  hasMaxExpSumY_ = 0;
  hasMinExpSumX_ = 0;
  hasMinExpSumY_ = 0;
    
  maxExpSumX_ = maxExpSumY_ = 0;
  minExpSumX_ = minExpSumY_ = 0;
}

void
GPin::setMaxExpSumX(float maxExpSumX) {
  hasMaxExpSumX_ = 1;
  maxExpSumX_ = maxExpSumX;
}

void
GPin::setMaxExpSumY(float maxExpSumY) {
  hasMaxExpSumY_ = 1;
  maxExpSumY_ = maxExpSumY;
}

void
GPin::setMinExpSumX(float minExpSumX) {
  hasMinExpSumX_ = 1;
  minExpSumX_ = minExpSumX;
}

void
GPin::setMinExpSumY(float minExpSumY) {
  hasMinExpSumY_ = 1;
  minExpSumY_ = minExpSumY;
}


void
GPin::updateLocation(const GCell* gCell) {
  cx_ = gCell->cx() + offsetCx_;
  cy_ = gCell->cy() + offsetCy_;
}

void
GPin::updateDensityLocation(const GCell* gCell) {
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
  electroForceX_(0), electroForceY_(0) {}

Bin::Bin(int x, int y, int lx, int ly, int ux, int uy, float targetDensity) 
  : Bin() {
  x_ = x;
  y_ = y;
  lx_ = lx; 
  ly_ = ly;
  ux_ = ux;
  uy_ = uy;
  targetDensity_ = targetDensity;
}

Bin::~Bin() {
  x_ = y_ = 0;
  lx_ = ly_ = ux_ = uy_ = 0;
  nonPlaceArea_ = instPlacedArea_ = fillerArea_ = 0;
  electroPhi_ = electroForceX_ = electroForceY_ = 0;
  density_ = targetDensity_ = 0;
}


const int64_t 
Bin::binArea() const { 
  return static_cast<int64_t>( dx() ) 
    * static_cast<int64_t>( dy() );
}

float
Bin::density() const {
  return density_;
}

float
Bin::targetDensity() const {
  return targetDensity_;
}

float
Bin::electroForceX() const {
  return electroForceX_;
}

float
Bin::electroForceY() const {
  return electroForceY_;
}

float
Bin::electroPhi() const {
  return electroPhi_;
}

void
Bin::setDensity(float density) {
  density_ = density;
}

void
Bin::setTargetDensity(float density) {
  targetDensity_ = density;
}

void
Bin::setElectroForce(float electroForceX, float electroForceY) {
  electroForceX_ = electroForceX;
  electroForceY_ = electroForceY;
}

void
Bin::setElectroPhi(float phi) {
  electroPhi_ = phi;
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
}

void BinGrid::setTargetDensity(float density)
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
  if (die_->placeInsts().size() != 0)
  {
    // assume no fixed instance
    int64_t averagePlaceInstArea = 0;
    for(Instance* inst : die_->placeInsts())
    {
      int64_t area = (int64_t)inst->dx() * inst->dy();
      averagePlaceInstArea += area;
    }
    averagePlaceInstArea /= die_->placeInsts().size();

    int64_t idealBinArea = 
      std::round(static_cast<float>(averagePlaceInstArea) / targetDensity_);
    int idealBinCnt = totalBinArea / idealBinArea; 
    
    LOG_INFO("TargetDensity: {}", targetDensity_);
    LOG_INFO("AveragePlaceInstArea: {}", averagePlaceInstArea);
    LOG_INFO("IdealBinArea: {}", idealBinArea);
    LOG_INFO("IdealBinCnt: {}", idealBinCnt);
    LOG_INFO("TotalBinArea: {}", totalBinArea);

    // find binCnt: 2, 4, 8, 16, 32, 64, ...
    // s.t. binCnt^2 <= idealBinCnt <= (binCnt*2)^2.
    for(foundBinCnt = 2; foundBinCnt <= 1024; foundBinCnt *= 2) {
      if( foundBinCnt * foundBinCnt <= idealBinCnt 
          && 4 * foundBinCnt * foundBinCnt > idealBinCnt ) {
        break;
      }
    }
  }
  // setBinCntX_;
  if( !isSetBinCntX_ ) {
    binCntX_ = foundBinCnt;
  }

  // setBinCntY_;
  if( !isSetBinCntY_ ) {
    binCntY_ = foundBinCnt;
  }


  LOG_INFO("BinCnt: {}, {}", binCntX_, binCntY_);
  
  binSizeX_ = ceil(
      static_cast<float>((ux_ - lx_))/binCntX_);
  binSizeY_ = ceil(
      static_cast<float>((uy_ - ly_))/binCntY_);
  
  LOG_INFO("BinSize: {}, {}", binSizeX_, binSizeY_);

  // initialize binStor_, bins_ vector
  binStor_.resize(binCntX_ * binCntY_);
  bins_.reserve(binCntX_ * binCntY_);
  int x = lx_, y = ly_;
  int idxX = 0, idxY = 0;
  for(auto& bin : binStor_) {

    int sizeX = (x + binSizeX_ > ux_)? 
      ux_ - x : binSizeX_;
    int sizeY = (y + binSizeY_ > uy_)? 
      uy_ - y : binSizeY_;

    //cout << "idxX: " << idxX << " idxY: " << idxY 
    //  << " x:" << x << " y:" << y 
    //  << " " << x+sizeX << " " << y+sizeY << endl;
    bin = Bin(idxX, idxY, x, y, x+sizeX, y+sizeY, targetDensity_);
    
    // move x, y coordinates.
    x += binSizeX_;
    idxX += 1;

    if( x >= ux_ ) {
      y += binSizeY_;
      x = lx_; 
      
      idxY ++;
      idxX = 0;
    }

    bins_.push_back( &bin );
  }

  LOG_INFO("NumBins: {}", bins_.size());

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

  for(Instance* inst : die_->fixedInsts())
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
    if( cell->isInstance() ) {
      // macro should have 
      // scale-down with target-density
      if( cell->isMacroInstance() ) {
        for(int i = pairX.first; i < pairX.second; i++) {
          for(int j = pairY.first; j < pairY.second; j++) {
            Bin* bin = bins_[ j * binCntX_ + i ];
            bin->addInstPlacedArea( 
                getOverlapDensityArea(bin, cell) 
                * cell->densityScale() * bin->targetDensity() ); 
          }
        }
      }
      // normal cells
      else if( cell->isStdInstance() ) {
        for(int i = pairX.first; i < pairX.second; i++) {
          for(int j = pairY.first; j < pairY.second; j++) {
            Bin* bin = bins_[ j * binCntX_ + i ];
            bin->addInstPlacedArea( 
                getOverlapDensityArea(bin, cell) 
                * cell->densityScale() ); 
          }
        }
      }
    }
    else if( cell->isFiller() ) {
      for(int i = pairX.first; i < pairX.second; i++) {
        for(int j = pairY.first; j < pairY.second; j++) {
          Bin* bin = bins_[ j * binCntX_ + i ];
          bin->addFillerArea( 
              getOverlapDensityArea(bin, cell) 
              * cell->densityScale() ); 
        }
      }
    }
  }  

  overflowArea_ = 0;

  // update density and overflowArea 
  // for nesterov use and FFT library
  for(auto& bin : bins_) {
    int64_t binArea = bin->binArea(); 
    bin->setDensity( 
        ( static_cast<float> (bin->instPlacedArea())
          + static_cast<float> (bin->fillerArea()) 
          + static_cast<float> (bin->nonPlaceArea()) )
        / static_cast<float>(binArea * bin->targetDensity()));

    overflowArea_ 
      += std::max(0.0f, 
          static_cast<float>(bin->instPlacedArea()) 
          + static_cast<float>(bin->nonPlaceArea())
          - (binArea * bin->targetDensity()));

  }
}


std::pair<int, int>
BinGrid::getDensityMinMaxIdxX(GCell* gcell) {
  int lowerIdx = (gcell->dLx() - lx())/binSizeX_;
  int upperIdx = 
   ( fastModulo((gcell->dUx() - lx()), binSizeX_) == 0)? 
   (gcell->dUx() - lx()) / binSizeX_ 
   : (gcell->dUx() - lx()) / binSizeX_ + 1;
  return std::make_pair(lowerIdx, upperIdx);
}

std::pair<int, int>
BinGrid::getDensityMinMaxIdxY(GCell* gcell) {
  int lowerIdx = (gcell->dLy() - ly())/binSizeY_;
  int upperIdx =
   ( fastModulo((gcell->dUy() - ly()), binSizeY_) == 0)? 
   (gcell->dUy() - ly()) / binSizeY_ 
   : (gcell->dUy() - ly()) / binSizeY_ + 1;

  return std::make_pair(lowerIdx, upperIdx);
}


std::pair<int, int>
BinGrid::getMinMaxIdxX(Instance* inst) {
  int lowerIdx = (inst->lx() - lx()) / binSizeX_;
  int upperIdx = 
   ( fastModulo((inst->ux() - lx()), binSizeX_) == 0)? 
   (inst->ux() - lx()) / binSizeX_ 
   : (inst->ux() - lx()) / binSizeX_ + 1;

  return std::make_pair(
      std::max(lowerIdx, 0),
      std::min(upperIdx, binCntX_));
}

std::pair<int, int>
BinGrid::getMinMaxIdxY(Instance* inst) {
  int lowerIdx = (inst->ly() - ly()) / binSizeY_;
  int upperIdx = 
   ( fastModulo((inst->uy() - ly()), binSizeY_) == 0)? 
   (inst->uy() - ly()) / binSizeY_ 
   : (inst->uy() - ly()) / binSizeY_ + 1;
  
  return std::make_pair(
      std::max(lowerIdx, 0),
      std::min(upperIdx, binCntY_));
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
    bin->setElectroForce(eForcePair.first, eForcePair.second);

    float electroPhi = fft_->getElectroPhi(bin->x(), bin->y());
    bin->setElectroPhi(electroPhi);

    sumPhi_ += electroPhi 
      * static_cast<float>(bin->nonPlaceArea() 
          + bin->instPlacedArea() + bin->fillerArea());
  }
}

void BinGrid::updateGCellDensityScaleAndSize()
{
  // update densitySize and densityScale in each gCell
  for(auto& gCell : gCells_)
  {
    float scaleX = 0, scaleY = 0;
    float densitySizeX = 0, densitySizeY = 0;
    if(gCell->dx() < REPLACE_SQRT2 * binSizeX_)
    {
      scaleX = static_cast<float>(gCell->dx()) 
        / static_cast<float>(REPLACE_SQRT2 * binSizeX_);
      densitySizeX = REPLACE_SQRT2 
        * static_cast<float>(binSizeX_);
    }
    else
    {
      scaleX = 1.0f;
      densitySizeX = gCell->dx();
    }

    if(gCell->dy() < REPLACE_SQRT2 * binSizeY_)
    {
      scaleY = static_cast<float>(gCell->dy()) 
               / static_cast<float>(REPLACE_SQRT2 * binSizeY_);
      densitySizeY = REPLACE_SQRT2 * static_cast<float>(binSizeY_);
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
  isSetBinCntX(0), isSetBinCntY(0) {}



void 
NesterovBaseVars::reset() {
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
  : pb_(nullptr), sumPhi_(0) {}

NesterovBase::NesterovBase(
    NesterovBaseVars nbVars, 
    std::shared_ptr<PlacerBase> pb)
  : NesterovBase() {
  nbVars_ = nbVars;
  pb_ = pb;
  init();
}

NesterovBase::~NesterovBase() {
  pb_ = nullptr;
}

void NesterovBase::init()
{
  LOG_TRACE("start NesterovBase::init");
  // gCellStor init
  gCellStor_.reserve(pb_->placeInsts().size());
  for(Instance* inst : pb_->placeInsts())
  {
    GCell myGCell(inst);
    // Check whether the given instance is macro or not
    if(inst->isMacro())
    {
      myGCell.setMacroInstance();
    }
    else
    {
      myGCell.setStdInstance();
    }
    gCellStor_.push_back(myGCell);
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

  LOG_TRACE("in NesterovBase::init create filler for each die");
  // create filler for each die
  std::unordered_map<Die*, std::pair<int, int>> fillerIdxMap;
  for(Die* die : pb_->dies())
  {
    int startIdx = static_cast<int>(gCellStor_.size());
    initFillerGCells(die);
    int endIdx = static_cast<int>(gCellStor_.size());
    fillerIdxMap.emplace(die, std::make_pair(startIdx, endIdx));
  }

  LOG_TRACE("in NesterovBase::init gCell ptr init");
  // gCell ptr init
  gCells_.reserve(gCellStor_.size());
  for(auto& gCell : gCellStor_)
  {
    gCells_.push_back(&gCell);
    if(gCell.isInstance())
    {
      gCellMap_[gCell.instance()] = &gCell;
    }
  }

  LOG_TRACE("in NesterovBase::init gPin ptr init");
// gPin ptr init
  gPins_.reserve(gPinStor_.size());
  for(auto& gPin : gPinStor_)
  {
    gPins_.push_back(&gPin);
    gPinMap_[gPin.pin()] = &gPin;
  }

  // gNet ptr init
  gNets_.reserve(gNetStor_.size());
  for(auto& gNet : gNetStor_)
  {
    gNets_.push_back(&gNet);
    gNetMap_[gNet.net()] = &gNet;
  }

  // gCellStor_'s pins_ fill
  for(auto& gCell : gCellStor_)
  {
    if( gCell.isFiller())
    {
      continue;
    }

    for(auto& pin : gCell.instance()->pins())
    {
      gCell.addGPin(placerToNesterov(pin));
    }
  }

  // gPinStor_' GNet and GCell fill
  for(auto& gPin : gPinStor_)
  {
    gPin.setGCell( 
        placerToNesterov(gPin.pin()->instance()));
    gPin.setGNet(
        placerToNesterov(gPin.pin()->net()));
  } 

  // gNetStor_'s GPin fill
  for(auto& gNet : gNetStor_)
  {
    for(auto& pin : gNet.net()->pins())
    {
      // NOTE: What if a pin is on a fixed instances?
      gNet.addGPin(placerToNesterov(pin));
    }
  }

  LOG_INFO("FillerInit: NumGCells: {}", gCells_.size());
  LOG_INFO("FillerInit: NumGNets: {}", gNets_.size());
  LOG_INFO("FillerInit: NumGPins: {}", gPins_.size());

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

  // binGrid ptr init
  binGrids_.reserve(binGridStor_.size());
  for(auto& bg : binGridStor_)
  {
    binGrids_.push_back(&bg);
    binGridMap_.emplace(bg.die(), &bg);
  }

  LOG_TRACE("in NesterovBase::init binGrids_ gCell fill");
  // binGrids_ gCell fill
  for(BinGrid* bg : binGrids_)
  {
    // assign inst gCell to binGrid
    for(Instance* inst : bg->die()->placeInsts())
    {
      GCell* gCell = placerToNesterov(inst);
      bg->addGCell(gCell);
      gCell->setBinGrid(bg);
    }
    // assign filler gCell to binGrid
    auto idx = fillerIdxMap[bg->die()];
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
void NesterovBase::initFillerGCells(Die* die)
{
  // if die does have any instance, we will not create any filler.
  if(die->placeInsts().size() == 0)
    return;

  // extract average dx/dy in range (10%, 90%)
  vector<int> dxStor;
  vector<int> dyStor;
  dxStor.reserve(die->placeInsts().size());
  dyStor.reserve(die->placeInsts().size());
  for(Instance* inst : die->placeInsts())
  {
    dxStor.push_back(inst->dx());
    dyStor.push_back(inst->dy());
  }
  
  // sort
  std::sort(dxStor.begin(), dxStor.end());
  std::sort(dyStor.begin(), dyStor.end());

  // average from (10 - 90%) .
  int64_t dxSum = 0, dySum = 0;
  int minIdx = static_cast<int>(dxStor.size()*0.05);
  int maxIdx = static_cast<int>(dxStor.size()*0.95);
  assert(maxIdx - minIdx > 0);
  for(int i=minIdx; i<maxIdx; i++)
  {
    dxSum += dxStor[i];
    dySum += dyStor[i];
  }

  // the avgDx and avgDy will be used as filler cells' 
  // width and height
  int avgDx = static_cast<int>(dxSum / (maxIdx - minIdx));
  int avgDy = static_cast<int>(dySum / (maxIdx - minIdx));

  // double dxSum = 0.0, dySum = 0.0;
  // int numStdcell = 0;
  // for(Instance* inst : die->insts())
  // {
  //   if(!inst->isMacro())
  //   {
  //     dxSum += inst->dx();
  //     dySum += inst->dy();
  //     numStdcell++;
  //   }
  // }
  // int avgDx = static_cast<int>(dxSum / numStdcell);
  // int avgDy = static_cast<int>(dySum / numStdcell);

  int64_t coreArea = (int64_t)die->coreDx() * die->coreDy();

  // nonPlaceInstsArea should not have targetDensity downscaling!!! 
  int64_t whiteSpaceArea = coreArea - die->fixedInstsArea();

  // TODO density screening
  int64_t movableArea = static_cast<int64_t>(whiteSpaceArea * nbVars_.targetDensity);
  
  int64_t totalFillerArea = movableArea 
    - static_cast<int64_t>(die->placeStdcellsArea())
    - static_cast<int64_t>(die->placeMacrosArea() * nbVars_.targetDensity);

  if(totalFillerArea < 0)
  {
    LOG_ERROR("Filler area is negative!!\n"
              "\tPlease put higher target density or\n"
              "\tRe-floorplan to have enough coreArea");
  }

  assert(avgDx> 0);
  assert(avgDy> 0);
  int fillerCnt = static_cast<int>(totalFillerArea / ((int64_t)avgDx * avgDy));

  LOG_INFO("FillerInit: CoreArea: {}", coreArea);
  LOG_INFO("FillerInit: WhiteSpaceArea: {}", whiteSpaceArea);
  LOG_INFO("FillerInit: MovableArea: {}", movableArea);
  LOG_INFO("FillerInit: TotalFillerArea: {}", totalFillerArea);
  LOG_INFO("FillerInit: NumFillerCells: {}", fillerCnt);
  LOG_INFO("FillerInit: FillerCellArea: {}", static_cast<int64_t>(avgDx*avgDy));
  LOG_INFO("FillerInit: FillerCellSize: {}, {}", avgDx, avgDy); 

  // 
  // mt19937 supports huge range of random values.
  // rand()'s RAND_MAX is only 32767.
  //
  mt19937 randVal(0);
  for(int i = 0; i < fillerCnt; i++) {

    // instability problem between g++ and clang++!
    auto randX = randVal();
    auto randY = randVal();

    // place filler cells on random coordi and
    // set size as avgDx and avgDy
    GCell myGCell(
        randX % die->coreDx() + die->coreLx(), 
        randY % die->coreDy() + die->coreLy(),
        avgDx, avgDy);
    gCellStor_.push_back(myGCell);
  }
}

GCell* NesterovBase::placerToNesterov(Instance* inst)
{
  auto gcPtr = gCellMap_.find(inst);
  auto result = (gcPtr == gCellMap_.end()) ? nullptr : gcPtr->second;
  assert(result != nullptr);
  return result;
}

GNet* NesterovBase::placerToNesterov(Net* net)
{
  auto gnPtr = gNetMap_.find(net);
  auto res = (gnPtr == gNetMap_.end()) ? nullptr : gnPtr->second;
  assert(res != nullptr);
  return res;
}

GPin* NesterovBase::placerToNesterov(Pin* pin)
{
  auto gpPtr = gPinMap_.find(pin);
  return (gpPtr == gPinMap_.end()) ? nullptr : gpPtr->second;
}

BinGrid* NesterovBase::placerToNesterov(Die* die)
{
  auto bgPtr = binGridMap_.find(die);
  return (bgPtr == binGridMap_.end()) ? nullptr : bgPtr->second;
}

// gcell update
void NesterovBase::updateGCellLocation(
  std::vector<FloatPoint>& coordis)
{
  for(auto& coordi : coordis)
  {
    int idx = &coordi - &coordis[0];
    gCells_[idx]->setLocation( coordi.x, coordi.y );
  }
}

// gcell update
void NesterovBase::updateGCellCenterLocation(
    std::vector<FloatPoint>& coordis)
{
  for(auto& coordi : coordis)
  {
    int idx = &coordi - &coordis[0];
    gCells_[idx]->setCenterLocation( coordi.x, coordi.y );
  }
}

void NesterovBase::updateGCellDensityCenterLocation(
    std::vector<FloatPoint>& coordis)
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

float NesterovBase::sumPhi() const
{
  return sumPhi_;
}

float NesterovBase::targetDensity() const
{
  return nbVars_.targetDensity;
}

void  NesterovBase::updateDensityCoordiLayoutInside(GCell* gCell)
{
  float targetLx = gCell->dLx();
  float targetLy = gCell->dLy();

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

float NesterovBase::getDensityCoordiLayoutInsideX(
  GCell* gCell, float cx)
{
  float adjVal = cx;
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

float NesterovBase::getDensityCoordiLayoutInsideY(
  GCell* gCell, float cy)
{
  float adjVal = cy;
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
    float wlCoeffX, float wlCoeffY)
{
  // clear all WA variables.
  for(auto& gNet : gNets_) {
    gNet->clearWaVars();
  }
  for(auto& gPin : gPins_) {
    gPin->clearWaVars();
  }

  for(auto& gNet : gNets_) {
    gNet->updateBox();

    for(auto& gPin : gNet->gPins()) {
      float expMinX = (gNet->lx() - gPin->cx()) * wlCoeffX; 
      float expMaxX = (gPin->cx() - gNet->ux()) * wlCoeffX;
      float expMinY = (gNet->ly() - gPin->cy()) * wlCoeffY;
      float expMaxY = (gPin->cy() - gNet->uy()) * wlCoeffY;

      // min x
      if(expMinX > nbVars_.minWireLengthForceBar) {
        gPin->setMinExpSumX( fastExp(expMinX) );
        gNet->addWaExpMinSumX( gPin->minExpSumX() );
        gNet->addWaXExpMinSumX( gPin->cx() 
            * gPin->minExpSumX() );
      }
      
      // max x
      if(expMaxX > nbVars_.minWireLengthForceBar) {
        gPin->setMaxExpSumX( fastExp(expMaxX) );
        gNet->addWaExpMaxSumX( gPin->maxExpSumX() );
        gNet->addWaXExpMaxSumX( gPin->cx() 
            * gPin->maxExpSumX() );
      }
     
      // min y 
      if(expMinY > nbVars_.minWireLengthForceBar) {
        gPin->setMinExpSumY( fastExp(expMinY) );
        gNet->addWaExpMinSumY( gPin->minExpSumY() );
        gNet->addWaYExpMinSumY( gPin->cy() 
            * gPin->minExpSumY() );
      }
      
      // max y
      if(expMaxY > nbVars_.minWireLengthForceBar) {
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
FloatPoint
NesterovBase::getWireLengthGradientWA(GCell* gCell, float wlCoeffX, float wlCoeffY) {
  FloatPoint gradientPair;

  for(auto& gPin : gCell->gPins()) {
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
FloatPoint
NesterovBase::getWireLengthGradientPinWA(GPin* gPin, float wlCoeffX, float wlCoeffY) {

  float gradientMinX = 0, gradientMinY = 0;
  float gradientMaxX = 0, gradientMaxY = 0;

  // min x
  if( gPin->hasMinExpSumX() ) {
    // from Net.
    float waExpMinSumX = gPin->gNet()->waExpMinSumX();
    float waXExpMinSumX = gPin->gNet()->waXExpMinSumX();

    gradientMinX = 
      ( waExpMinSumX * ( gPin->minExpSumX() * ( 1.0 - wlCoeffX * gPin->cx()) ) 
          + wlCoeffX * gPin->minExpSumX() * waXExpMinSumX )
        / ( waExpMinSumX * waExpMinSumX );
  }
  
  // max x
  if( gPin->hasMaxExpSumX() ) {
    
    float waExpMaxSumX = gPin->gNet()->waExpMaxSumX();
    float waXExpMaxSumX = gPin->gNet()->waXExpMaxSumX();
    
    gradientMaxX = 
      ( waExpMaxSumX * ( gPin->maxExpSumX() * ( 1.0 + wlCoeffX * gPin->cx()) ) 
          - wlCoeffX * gPin->maxExpSumX() * waXExpMaxSumX )
        / ( waExpMaxSumX * waExpMaxSumX );

  }

  // min y
  if( gPin->hasMinExpSumY() ) {
    
    float waExpMinSumY = gPin->gNet()->waExpMinSumY();
    float waYExpMinSumY = gPin->gNet()->waYExpMinSumY();

    gradientMinY = 
      ( waExpMinSumY * ( gPin->minExpSumY() * ( 1.0 - wlCoeffY * gPin->cy()) ) 
          + wlCoeffY * gPin->minExpSumY() * waYExpMinSumY )
        / ( waExpMinSumY * waExpMinSumY );
  }
  
  // max y
  if( gPin->hasMaxExpSumY() ) {
    
    float waExpMaxSumY = gPin->gNet()->waExpMaxSumY();
    float waYExpMaxSumY = gPin->gNet()->waYExpMaxSumY();
    
    gradientMaxY = 
      ( waExpMaxSumY * ( gPin->maxExpSumY() * ( 1.0 + wlCoeffY * gPin->cy()) ) 
          - wlCoeffY * gPin->maxExpSumY() * waYExpMaxSumY )
        / ( waExpMaxSumY * waExpMaxSumY );
  }

  return FloatPoint(gradientMinX - gradientMaxX, 
      gradientMinY - gradientMaxY);
}


FloatPoint NesterovBase::getWireLengthPreconditioner(GCell* gCell)
{
  return FloatPoint(gCell->gPins().size(),
     gCell->gPins().size());
}

FloatPoint NesterovBase::getDensityPreconditioner(GCell* gCell)
{
  float areaVal = static_cast<float>(gCell->dx()) 
    * static_cast<float>(gCell->dy());

  return FloatPoint(areaVal, areaVal);
}

// get GCells' electroForcePair
// i.e. get DensityGradient with given GCell
FloatPoint NesterovBase::getDensityGradient(GCell* gCell)
{
  // find 
  BinGrid* bg = gCell->binGrid();

  auto pairX = bg->getDensityMinMaxIdxX(gCell);
  auto pairY = bg->getDensityMinMaxIdxY(gCell);
  
  FloatPoint electroForce;

  for(int i = pairX.first; i < pairX.second; i++)
  {
    for(int j = pairY.first; j < pairY.second; j++)
    {
      Bin* bin = bg->bins()[j * bg->binCntX() + i];
      float overlapArea 
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

int64_t NesterovBase::getHpwl()
{
  int64_t hpwl = 0;
  for(auto& gNet : gNets_)
  {
    gNet->updateBox();
    hpwl += gNet->hpwl();
  }
  return hpwl;
}

void NesterovBase::reset()
{ 
  pb_ = nullptr;
  nbVars_.reset();
}



// https://stackoverflow.com/questions/33333363/built-in-mod-vs-custom-mod-function-improve-the-performance-of-modulus-op
static int 
fastModulo(const int input, const int ceil) {
  return input >= ceil? input % ceil : input;
}

// int64_t is recommended, but float is 2x fast
static float 
getOverlapDensityArea(Bin* bin, GCell* cell) {
  int rectLx = max(bin->lx(), cell->dLx()), 
      rectLy = max(bin->ly(), cell->dLy()),
      rectUx = min(bin->ux(), cell->dUx()), 
      rectUy = min(bin->uy(), cell->dUy());
  
  if( rectLx >= rectUx || rectLy >= rectUy ) {
    return 0;
  }
  else {
    return static_cast<float>(rectUx - rectLx) 
      * static_cast<float>(rectUy - rectLy);
  }
}


static int64_t
getOverlapArea(Bin* bin, Instance* inst) {
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
static float
fastExp(float a) {
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
