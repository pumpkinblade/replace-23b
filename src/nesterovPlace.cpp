#include "placerBase.h"
#include "nesterovBase.h"
#include "nesterovPlace.h"
#include "log.h"
#include "plot.h"
#include <glpk.h>

using namespace std;

namespace replace {

static prec getDistance(vector<Point>& a, vector<Point>& b);
static prec getSecondNorm(vector<Point>& a);

static void getRotBox(const GCell *cell, bool rot, double &lx, double &ux, double &ly, double &uy);
static double getRotOverlap(const GCell* cell1, const  GCell* cell2, bool rot1, bool rot2);
static void ilpSolve(const std::vector<GCell*>& macros, std::vector<bool>& isRot);
static Orientation findNearestOrientation(prec theta);

NesterovPlaceVars::NesterovPlaceVars()
  : maxNesterovIter(2000), 
  maxBackTrack(10),
  initDensityPenalty(0.00008),
  initWireLengthCoef(0.25),
  targetOverflow(0.1),
  minPhiCoef(0.95),
  maxPhiCoef(1.05),
  minPreconditioner(1.0),
  initialPrevCoordiUpdateCoef(100),
  referenceHpwl(446000000),
  useLocalDensity(false),
  initAlpha(1e-12f),
  initBeta(1e-11f),
  useTheta(false) {}

NesterovPlace::NesterovPlace() 
  : nb_(nullptr), npVars_(), 
    wireLengthGradSum_(0), 
    densityGradSum_(0),
    stepLength_(0),
    densityPenalty_(0),
    baseWireLengthCoef_(0), 
    wireLengthCoefX_(0), 
    wireLengthCoefY_(0),
    prevHpwl_(0),
    isDiverged_(false) {}

NesterovPlace::NesterovPlace(NesterovPlaceVars npVars, std::shared_ptr<NesterovBase> nb) 
    : NesterovPlace() 
{
  npVars_ = npVars;
  nb_ = nb;
}

void NesterovPlace::init() {
  LOG_TRACE("NesterovInit Begin");

  size_t gCellSize = nb_->gCells().size();
  curSLPCoordi_.resize(gCellSize, Point());
  curSLPSumGrads_.resize(gCellSize, Point());

  nextSLPCoordi_.resize(gCellSize, Point());
  nextSLPSumGrads_.resize(gCellSize, Point());
  
  prevSLPCoordi_.resize(gCellSize, Point());
  prevSLPSumGrads_.resize(gCellSize, Point());

  curCoordi_.resize(gCellSize, Point());
  nextCoordi_.resize(gCellSize, Point());

  if(npVars_.useLocalDensity)
  {
    localAlpha_ = npVars_.initAlpha;
    localBeta_ = npVars_.initBeta;
    curCellDelta_.resize(gCellSize, 0.0f);
    nextCellDelta_.resize(gCellSize, 0.0f);
  }

  if(npVars_.useTheta)
  {
    for(int i = 0; i < gCellSize; i++)
    {
      if(nb_->gCells()[i]->isInstance() && nb_->gCells()[i]->isMacro())
        macroIndices_.push_back(i);
    }
    curTheta_.resize(macroIndices_.size(), 0.0f);
    nextTheta_.resize(macroIndices_.size(), 0.0f);
    curSLPTheta_.resize(macroIndices_.size(), 0.0f);
    curSLPSumGradTheta_.resize(macroIndices_.size(), 0.0f);
    nextSLPTheta_.resize(macroIndices_.size(), 0.0f);
    nextSLPSumGradTheta_.resize(macroIndices_.size(), 0.0f);
    prevSLPTheta_.resize(macroIndices_.size(), 0.0f);
    prevSLPSumGradTheta_.resize(macroIndices_.size(), 0.0f);
  }

  for(int i = 0; i < gCellSize; i++)
  {
    nb_->updateDensityCoordiLayoutInside(nb_->gCells()[i]);
    curSLPCoordi_[i] 
      = prevSLPCoordi_[i] 
      = curCoordi_[i] 
      = Point(nb_->gCells()[i]->cx(), nb_->gCells()[i]->cy());
  }

  if(npVars_.useTheta)
  {
    for(int idx : macroIndices_)
    {
      curSLPTheta_[idx]
        = prevSLPTheta_[idx]
        = curTheta_[idx]
        = nb_->gCells()[idx]->theta();
    }
  }

  // bin update
  nb_->updateGCellCenterLocation(curSLPCoordi_);
  
  prevHpwl_ = nb_->hpwl();

  LOG_DEBUG("InitialHPWL: {}", prevHpwl_);

  // FFT update
  nb_->updateDensityForceBin();

  prec avgBinSizeX = 0.f; 
  prec avgBinSizeY = 0.f;
  for(BinGrid* bg : nb_->binGrids())
  {
    avgBinSizeX += static_cast<prec>(bg->binSizeX());
    avgBinSizeY += static_cast<prec>(bg->binSizeY());
  }
  avgBinSizeX /= static_cast<prec>(nb_->binGrids().size());
  avgBinSizeY /= static_cast<prec>(nb_->binGrids().size());
  baseWireLengthCoef_  = npVars_.initWireLengthCoef / (0.5f * (avgBinSizeX + avgBinSizeY));

  LOG_DEBUG("BaseWireLengthCoef: {}", baseWireLengthCoef_);
  
  sumOverflow_ = nb_->overflow();

  LOG_DEBUG("InitSumOverflow: {}", sumOverflow_);

  updateWireLengthCoef(sumOverflow_);

  LOG_DEBUG("WireLengthCoef: {}, {}", wireLengthCoefX_, wireLengthCoefY_);

  // WL update
  nb_->updateWireLengthForceWA(wireLengthCoefX_, wireLengthCoefY_);
 
  // fill in curSLPSumGrads_
  updateGradients(curSLPSumGrads_);

  if( isDiverged_ ) {
    return;
  }

  // approximately fill in prevSLPCoordi_ to calculate lc vars
  updateInitialPrevSLPCoordi();

  // bin, FFT, wlen update with prevSLPCoordi.
  nb_->updateGCellCenterLocation(prevSLPCoordi_);
  nb_->updateDensityForceBin();
  nb_->updateWireLengthForceWA(wireLengthCoefX_, wireLengthCoefY_);
  
  // update previSumGrads_
  updateGradients(prevSLPSumGrads_);

  if(isDiverged_)
    return;
  
  LOG_DEBUG("WireLengthGradSum: {}", wireLengthGradSum_);
  LOG_DEBUG("DensityGradSum", densityGradSum_);

  densityPenalty_ 
    = (wireLengthGradSum_ / densityGradSum_ )
    * npVars_.initDensityPenalty; 
  
  LOG_DEBUG("InitDensityPenalty: {}", densityPenalty_);
  
  sumOverflow_ = nb_->overflow();
  
  LOG_DEBUG("PrevSumOverflow: {}", sumOverflow_);
  
  stepLength_  
    = getStepLength (prevSLPCoordi_, prevSLPSumGrads_, curSLPCoordi_, curSLPSumGrads_);


  LOG_DEBUG("InitialStepLength: {}", stepLength_);
  LOG_TRACE("NesterovInit End");
}

// to execute following function,
// 
// nb_->updateGCellDensityCenterLocation(coordi); // bin update
// nb_->updateDensityForceBin(); // bin Force update
//  
// nb_->updateWireLengthForceWA(wireLengthCoefX_, wireLengthCoefY_); // WL update
//
void NesterovPlace::updateGradients(std::vector<Point>& sumGrads)
{
  wireLengthGradSum_ = 0;
  densityGradSum_ = 0;
  localDensityGradSum_ = 0;
  prec gradSum = 0;
  prec cellDeltaSum = 0;
  LOG_DEBUG("Density Penalty: {}", densityPenalty_);

  for(size_t i=0; i<nb_->gCells().size(); i++)
  {
    GCell* gCell = nb_->gCells().at(i);

    Point wireLengthGrad = nb_->getWireLengthGradientWA(gCell, wireLengthCoefX_, wireLengthCoefY_);
    wireLengthGradSum_ += fabs(wireLengthGrad.x);
    wireLengthGradSum_ += fabs(wireLengthGrad.y);

    Point densityGrad = nb_->getDensityGradient(gCell); 
    densityGradSum_ += fabs(densityGrad.x);
    densityGradSum_ += fabs(densityGrad.y);

    sumGrads[i].x = wireLengthGrad.x + densityPenalty_ * densityGrad.x;
    sumGrads[i].y = wireLengthGrad.y + densityPenalty_ * densityGrad.y;

    if(npVars_.useLocalDensity)
    {
      prec delta = curCellDelta_[i];
      Point lgrad = nb_->getDensityGradientLocal(gCell, localAlpha_, localBeta_, delta);
      nextCellDelta_[i] = delta;
      localDensityGradSum_ += fabs(lgrad.x);
      localDensityGradSum_ += fabs(lgrad.y);
      cellDeltaSum += delta;

      sumGrads[i].x += delta * lgrad.x;
      sumGrads[i].y += delta * lgrad.y;
    }

    Point wireLengthPreCondi = nb_->getWireLengthPreconditioner(gCell);
    Point densityPrecondi = nb_->getDensityPreconditioner(gCell);
    Point sumPrecondi(
        wireLengthPreCondi.x + densityPenalty_ * densityPrecondi.x,
        wireLengthPreCondi.y + densityPenalty_ * densityPrecondi.y);

    sumPrecondi.x = std::max(sumPrecondi.x, npVars_.minPreconditioner);
    sumPrecondi.y = std::max(sumPrecondi.y, npVars_.minPreconditioner);
    
    sumGrads[i].x /= sumPrecondi.x;
    sumGrads[i].y /= sumPrecondi.y; 

    gradSum += fabs(sumGrads[i].x) + fabs(sumGrads[i].y);
  }
  
  if(npVars_.useLocalDensity)
  {
    cellDeltaSum /= nb_->gCells().size();
    LOG_DEBUG("CellDeltaAverage: {}", cellDeltaSum);
  }

  LOG_DEBUG("WireLengthGradSum: {}", wireLengthGradSum_);
  LOG_DEBUG("DensityGradSum: {}", densityGradSum_);
  if(npVars_.useLocalDensity)
    LOG_DEBUG("LocalDensityGradSum: {}", localDensityGradSum_);
  LOG_DEBUG("GradSum: {}", gradSum);

  // divergence detection
  isDiverged_ = (std::isnan(gradSum) || std::isinf(gradSum));
}

void
NesterovPlace::doNesterovPlace(string placename) {
  LOG_TRACE("start NesterovPlace::doNesterovPlace");

  // if replace diverged in init() function, 
  // replace must be skipped.
  if( isDiverged_ ) {
    LOG_ERROR("RePlAce diverged. Please tune the parameters again");
    return;
  }

  if (placename != "" && placename[placename.length() - 1] != '_') {
    placename.push_back('_');
  }
  Plot::plot(nb_.get(), PlotNesterovType::GCell, "./plot/cell", placename + "cell_0");
  Plot::plot(nb_.get(), PlotNesterovType::Arrow, "./plot/arrow", placename + "arrow_0");

  // backTracking variable.
  prec curA = 1.0;

  // divergence detection
  prec minSumOverflow = 1e30;
  double hpwlWithMinSumOverflow = 1e30;

  // dynamic adjustment of max_phi_coef
  bool isMaxPhiCoefChanged = false;

  // diverge error handling
  string divergeMsg = "";
  int divergeCode = 0;

  // Core Nesterov Loop
  for(int i=0; i<npVars_.maxNesterovIter; i++) {
    LOG_DEBUG("Iter: {}", i+1);
    
    prec prevA = curA;

    // here, prevA is a_(k), curA is a_(k+1)
    // See, the papers' Algorithm 4 section
    //
    curA = (1.0 + sqrt(4.0 * prevA * prevA + 1.0)) * 0.5;

    // coeff is (a_k -1) / ( a_(k+1)) in paper.
    prec coeff = (prevA - 1.0)/curA;
    
    LOG_DEBUG("PreviousA: {}", prevA);
    LOG_DEBUG("CurrentA: {}", curA);
    LOG_DEBUG("Coefficient: {}", coeff);
    LOG_DEBUG("StepLength: {}", stepLength_);

    // Back-Tracking loop
    int numBackTrak = 0;
    for(numBackTrak = 0; numBackTrak < npVars_.maxBackTrack; numBackTrak++) {
      
      // fill in nextCoordinates with given stepLength_
      for(size_t k=0; k<nb_->gCells().size(); k++) {
        Point nextCoordi(
          curSLPCoordi_[k].x + stepLength_ * curSLPSumGrads_[k].x,
          curSLPCoordi_[k].y + stepLength_ * curSLPSumGrads_[k].y );

        Point nextSLPCoordi(
          nextCoordi.x + coeff * (nextCoordi.x - curCoordi_[k].x),
          nextCoordi.y + coeff * (nextCoordi.y - curCoordi_[k].y));

        GCell* curGCell = nb_->gCells()[k];

        nextCoordi_[k] 
          = Point( 
              nb_->getDensityCoordiLayoutInsideX( 
                curGCell, nextCoordi.x),
              nb_->getDensityCoordiLayoutInsideY(
                curGCell, nextCoordi.y));
        
        nextSLPCoordi_[k]
          = Point(
              nb_->getDensityCoordiLayoutInsideX(
                curGCell, nextSLPCoordi.x),
              nb_->getDensityCoordiLayoutInsideY(
                curGCell, nextSLPCoordi.y));
      }
 

      nb_->updateGCellCenterLocation(nextSLPCoordi_);
      nb_->updateDensityForceBin();
      nb_->updateWireLengthForceWA(wireLengthCoefX_, wireLengthCoefY_);

      updateGradients(nextSLPSumGrads_);

      // NaN or inf is detected in WireLength/Density Coef 
      if( isDiverged_ ) {
        break;
      }
  
      prec newStepLength  
        = getStepLength (curSLPCoordi_, curSLPSumGrads_, nextSLPCoordi_, nextSLPSumGrads_);
     
      LOG_DEBUG("NetStepLength: {}", newStepLength);

      if( newStepLength > stepLength_ * 0.95) {
        stepLength_ = newStepLength;
        break;
      }
      else {
        stepLength_ = newStepLength;
      } 
    }

    LOG_DEBUG("NumBackTrak: {}", numBackTrak+1);

    // dynamic adjustment for
    // better convergence with
    // large designs 
    if( !isMaxPhiCoefChanged && sumOverflow_ 
        < 0.35f ) {
      isMaxPhiCoefChanged = true;
      npVars_.maxPhiCoef *= 0.99;
    }

    // usually, maxBackTrack should be 1~3
    // 10 is the case when
    // all of cells are not moved at all.
    if( npVars_.maxBackTrack == numBackTrak ) {
      divergeMsg = "RePlAce divergence detected. \n";
      divergeMsg += "        Please decrease init_density_penalty value";
      divergeCode = 3;
      isDiverged_ = true;
    } 

    if( isDiverged_ ) {
      break;
    }

    updateNextIter();

    // For JPEG Saving
    // debug

    if( i == 0 || (i+1) % 20 == 0 ) {
      LOG_DEBUG("[NesterovSolve] Iter: {} overflow: {} HPWL: {}", i+1, sumOverflow_, prevHpwl_);
      Plot::plot(nb_.get(), PlotNesterovType::GCell, "./plot/cell", placename + "cell_" + std::to_string(i+1));
      Plot::plot(nb_.get(), PlotNesterovType::Arrow, "./plot/arrow", placename + "arrow_" + std::to_string(i+1));
    }

    if( minSumOverflow > sumOverflow_ ) {
      minSumOverflow = sumOverflow_;
      hpwlWithMinSumOverflow = prevHpwl_; 
    }

    // diverge detection on
    // large max_phi_cof value + large design 
    //
    // 1) happen overflow < 20%
    // 2) Hpwl is growing
    //
    if( sumOverflow_ < 0.3f 
        && sumOverflow_ - minSumOverflow >= 0.02f
        && hpwlWithMinSumOverflow * 1.2f < prevHpwl_ ) {
      divergeMsg = "RePlAce divergence detected. \n";
      divergeMsg += "        Please decrease max_phi_cof value";
      divergeCode = 4;
      isDiverged_ = true;
      break;
    }

    // minimum iteration is 50
    if( i > 50 && sumOverflow_ <= npVars_.targetOverflow) {
      LOG_DEBUG("[NesterovSolve] Finished with Overflow: {}", sumOverflow_);
      break;
    }
  }
 
  // in all case including diverge, 
  // PlacerBase should be updated. 
  updatePlacerBase();
  Plot::plot(nb_.get(), PlotNesterovType::GCell, "./plot/cell", placename + "cell_end");
  Plot::plot(nb_.get(), PlotNesterovType::Arrow, "./plot/arrow", placename + "arrow_end");

  if( isDiverged_ ) { 
    LOG_ERROR("{} : Code `{}`", divergeMsg, divergeCode);
  }
}

void
NesterovPlace::updateWireLengthCoef(prec overflow) {
  if( overflow > 1.0 ) {
    wireLengthCoefX_ = wireLengthCoefY_ = 0.1;
  }
  else if( overflow < 0.1 ) {
    wireLengthCoefX_ = wireLengthCoefY_ = 10.0;
  }
  else {
    wireLengthCoefX_ = wireLengthCoefY_ 
      = 1.0 / pow(10.0, (overflow-0.1)*20 / 9.0 - 1.0);
  }

  wireLengthCoefX_ *= baseWireLengthCoef_;
  wireLengthCoefY_ *= baseWireLengthCoef_;
  LOG_DEBUG("NewWireLengthCoef: {}", wireLengthCoefX_);
}

void
NesterovPlace::updateInitialPrevSLPCoordi() {
  for(size_t i=0; i<nb_->gCells().size(); i++) {
    GCell* curGCell = nb_->gCells()[i];


    prec prevCoordiX 
      = curSLPCoordi_[i].x + npVars_.initialPrevCoordiUpdateCoef 
      * curSLPSumGrads_[i].x;
  
    prec prevCoordiY
      = curSLPCoordi_[i].y + npVars_.initialPrevCoordiUpdateCoef
      * curSLPSumGrads_[i].y;
    
    Point newCoordi( 
      nb_->getDensityCoordiLayoutInsideX( curGCell, prevCoordiX),
      nb_->getDensityCoordiLayoutInsideY( curGCell, prevCoordiY) );

    prevSLPCoordi_[i] = newCoordi;
  } 
}

void
NesterovPlace::updateNextIter() {
  // swap vector pointers
  std::swap(prevSLPCoordi_, curSLPCoordi_);
  std::swap(prevSLPSumGrads_, curSLPSumGrads_);
  
  std::swap(curSLPCoordi_, nextSLPCoordi_);
  std::swap(curSLPSumGrads_, nextSLPSumGrads_);

  std::swap(curCoordi_, nextCoordi_);

  if(npVars_.useLocalDensity)
  {
    std::swap(curCellDelta_, nextCellDelta_);
  }

  sumOverflow_ = nb_->overflow();

  LOG_DEBUG("Gradient: {}", getSecondNorm(curSLPSumGrads_));
  LOG_DEBUG("Phi: {}", nb_->sumPhi());
  LOG_DEBUG("Overflow: {}", sumOverflow_);

  updateWireLengthCoef(sumOverflow_);
  double hpwl = nb_->hpwl();
  
  LOG_DEBUG("PreviousHPWL: {}", prevHpwl_);
  LOG_DEBUG("NewHPWL: {}", hpwl);

  prec phiCoef = getPhiCoef( 
      static_cast<prec>(hpwl - prevHpwl_) 
      / npVars_.referenceHpwl );
  
  prevHpwl_ = hpwl;
  densityPenalty_ *= phiCoef;

  if(npVars_.useLocalDensity)
  {
    localAlpha_ *= phiCoef;
    localBeta_ *= phiCoef;
  }

  LOG_DEBUG("PhiCoef: {}", phiCoef);
}

prec
NesterovPlace::getStepLength(
    std::vector<Point>& prevSLPCoordi_,
    std::vector<Point>& prevSLPSumGrads_,
    std::vector<Point>& curSLPCoordi_,
    std::vector<Point>& curSLPSumGrads_ ) {

  prec coordiDistance 
    = getDistance(prevSLPCoordi_, curSLPCoordi_);
  prec gradDistance 
    = getDistance(prevSLPSumGrads_, curSLPSumGrads_);

  LOG_DEBUG("CoordinateDistance: {}", coordiDistance);
  LOG_DEBUG("GradientDistance: {}", gradDistance);

  return coordiDistance / gradDistance;
}

prec
NesterovPlace::getPhiCoef(prec scaledDiffHpwl) {
  LOG_DEBUG("Input ScaleDiffHPWL", scaledDiffHpwl);

  prec retCoef 
    = (scaledDiffHpwl < 0)? 
    npVars_.maxPhiCoef: 
    npVars_.maxPhiCoef * pow( npVars_.maxPhiCoef, scaledDiffHpwl * -1.0 );
  retCoef = std::max(npVars_.minPhiCoef, retCoef);
  return retCoef;
}

void
NesterovPlace::updatePlacerBase() {
  for(auto& gCell : nb_->gCells()) {
    if(gCell->isInstance()) {
      int lx = static_cast<int>(std::round(gCell->lx()));
      int ly = static_cast<int>(std::round(gCell->ly()));
      gCell->instance()->setLocation(lx, ly);
    }
  }
}

void NesterovPlace::determinMacroOrient()
{
  for (BinGrid* bg : nb_->binGrids())
  {
    const Technology& tech = *bg->die()->tech();
    std::vector<GCell*> macros;
    for (int idx : macroIndices_)
    {
      GCell* gcell = nb_->gCells()[idx];
      if(gcell->binGrid() != bg)
        continue;
      Orientation ori = findNearestOrientation(gcell->theta());
      int i = static_cast<int>(ori);
      if(std::abs(LEGAL_THETA[idx] - gcell->theta()) < 0.1 * PI)
      {
        gcell->setThetaNoUpdatePin(LEGAL_THETA[i]);
        gcell->instance()->setOrientSize(tech, ori);
      }
      else
      {
        macros.push_back(gcell);
      }
    }
    std::vector<bool> isRot(macros.size(), false);
    ilpSolve(macros, isRot);
    for (int u = 0; u < macros.size(); u++)
    {
      if (isRot[u])
      {
        if (std::abs(macros[u]->theta() - 0.5 * PI) < std::abs(macros[u]->theta() - 1.5 * PI))
        {
          macros[u]->setThetaNoUpdatePin(0.5 * PI);
          macros[u]->instance()->setOrientSize(tech, Orientation::R90);
        }
        else
        {
          macros[u]->setThetaNoUpdatePin(1.5 * PI);
          macros[u]->instance()->setOrientSize(tech, Orientation::R270);
        }
      }
      else
      {
        if (std::abs(macros[u]->theta()) < std::abs(macros[u]->theta() - PI)
            || std::abs(macros[u]->theta() - 2.0 * PI) < std::abs(macros[u]->theta() - PI))
        {
          macros[u]->setThetaNoUpdatePin(0.0);
          macros[u]->instance()->setOrientSize(tech, Orientation::R0);
        }
        else
        {
          macros[u]->setThetaNoUpdatePin(PI);
          macros[u]->instance()->setOrientSize(tech, Orientation::R180);
        }
      }
    }
  }
}

static prec getDistance(vector<Point>& a, vector<Point>& b)
{
  prec sumDistance = 0.0f;
  for(size_t i=0; i<a.size(); i++)
  {
    sumDistance += (a[i].x - b[i].x) * (a[i].x - b[i].x);
    sumDistance += (a[i].y - b[i].y) * (a[i].y - b[i].y);
  }

  return sqrt(sumDistance / (2.0 * a.size()));
}

static prec getSecondNorm(vector<Point>& a)
{
  prec norm = 0;
  for(auto& coordi : a)
    norm += coordi.x * coordi.x + coordi.y * coordi.y;
  return sqrt( norm / (2.0*a.size()) ); 
}

static void getRotBox(const GCell *cell, bool rot, double &lx, double &ux, double &ly, double &uy)
{
  if (rot)
  {
    double cx = cell->cx();
    double dx = cell->dx();
    double cy = cell->cy();
    double dy = cell->dy();
    lx = cx - 0.5 * dy;
    ux = cx + 0.5 * dy;
    ly = cy - 0.5 * dx;
    uy = cy + 0.5 * dx;
  }
  else
  {
    lx = cell->lx();
    ux = cell->ux();
    ly = cell->ly();
    uy = cell->uy();
  }
}

static double getRotOverlap(const GCell* cell1, const  GCell* cell2, bool rot1, bool rot2)
{
  double lx1, ux1, ly1, uy1;
  getRotBox(cell1, rot1, lx1, ux1, ly1, uy1);
  double lx2, ux2, ly2, uy2;
  getRotBox(cell2, rot2, lx2, ux2, ly2, uy2);

  double lx = std::max(lx1, lx2);
  double ux = std::min(ux1, ux2);
  double ly = std::max(ly1, ly2);
  double uy = std::min(uy1, uy2);

  if (lx < ux && ly < uy)
    return (ux - lx) * (uy - ly);
  else
    return 0.0;
}

static void ilpSolve(const std::vector<GCell*>& macros, std::vector<bool>& isRot)
{
  int numMacros = static_cast<int>(macros.size());
  std::vector<double> coeff(numMacros * (numMacros + 1) / 2, 0.0);
  for (int u = 0, j = numMacros; u < numMacros; u++)
  {
    for (int v = u + 1; v < numMacros; v++)
    {
      double phi_uv = getRotOverlap(macros[u], macros[v], false, false);
      double phi_uRv = getRotOverlap(macros[u], macros[v], true, false);
      double phi_uvR = getRotOverlap(macros[u], macros[v], false, true);
      double phi_uRvR = getRotOverlap(macros[u], macros[v], true, true);
      coeff[u] += (phi_uRv - phi_uv);
      coeff[v] += (phi_uvR - phi_uv);
      coeff[j++] = (phi_uv - phi_uRv - phi_uvR + phi_uRvR);
    }
  }

  glp_prob *lp = glp_create_prob();
  glp_set_prob_name(lp, "orientation");
  glp_set_obj_dir(lp, GLP_MIN);
  // cols
  int n = static_cast<int>(coeff.size());
  char name[16];
  glp_add_cols(lp, n);
  for (int u = 0, j = numMacros; u < numMacros; u++)
  {
    std::sprintf(name, "r%d", u + 1);
    glp_set_col_name(lp, u + 1, name);
    glp_set_obj_coef(lp, u + 1, coeff[u]);
    glp_set_col_kind(lp, u + 1, GLP_BV);
    for (int v = u + 1; v < numMacros; v++)
    {
      std::sprintf(name, "r%d_%d", u + 1, v + 1);
      glp_set_col_name(lp, j + 1, name);
      glp_set_obj_coef(lp, j + 1, coeff[j]);
      glp_set_col_kind(lp, j + 1, GLP_BV);
      j++;
    }
  }
  // rows
  int m = 2 * n;
  glp_add_rows(lp, m);
  for (int u = 0, i = 0, j = numMacros; u < numMacros; u++)
  {
    for (int v = u + 1; v < numMacros; v++)
    {
      int ind[1 + 3];
      double val[1 + 3];

      glp_set_row_bnds(lp, i + 1, GLP_UP, 0.0, 0.0);
      ind[1] = u + 1;
      ind[2] = j + 1;
      val[1] = -1.0;
      val[2] = 1.0;
      glp_set_mat_row(lp, i + 1, 2, ind, val);

      glp_set_row_bnds(lp, i + 2, GLP_UP, 0.0, 1.0);
      ind[1] = u + 1;
      ind[2] = v + 1;
      ind[3] = j + 1;
      val[1] = 1.0;
      val[2] = 1.0;
      val[3] = -1.0;
      glp_set_mat_row(lp, i + 2, 3, ind, val);

      i += 2;
      j++;
    }
  }

  glp_simplex(lp, NULL);
  glp_intopt(lp, NULL);
  for (int u = 0; u < numMacros; u++)
    isRot[u] = static_cast<bool>(glp_mip_col_val(lp, u + 1));
  glp_delete_prob(lp);
}

static Orientation findNearestOrientation(prec theta)
{
  double bestDist = std::abs(theta);
  Orientation bestOri = Orientation::R0;
  for(int idx = 1; idx < 5; idx++)
  {
    double dist = std::abs(theta - LEGAL_THETA[idx]);
    if(dist < bestDist)
    {
      bestDist = dist;
      bestOri = static_cast<Orientation>(idx);
    }
  }
  return bestOri;
}

}
