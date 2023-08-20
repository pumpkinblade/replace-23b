#include "placerBase.h"
#include "nesterovBase.h"
#include "nesterovPlace.h"
#include "log.h"
#include <iostream>
using namespace std;

#include "plot.h"

namespace replace {

static prec
getDistance(vector<Point>& a, vector<Point>& b);

static prec
getSecondNorm(vector<Point>& a);

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
  referenceHpwl(446000000) {}

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
  curSLPWireLengthGrads_.resize(gCellSize, Point());
  curSLPDensityGrads_.resize(gCellSize, Point());
  curSLPSumGrads_.resize(gCellSize, Point());

  nextSLPCoordi_.resize(gCellSize, Point());
  nextSLPWireLengthGrads_.resize(gCellSize, Point());
  nextSLPDensityGrads_.resize(gCellSize, Point());
  nextSLPSumGrads_.resize(gCellSize, Point());
  
  prevSLPCoordi_.resize(gCellSize, Point());
  prevSLPWireLengthGrads_.resize(gCellSize, Point());
  prevSLPDensityGrads_.resize(gCellSize, Point());
  prevSLPSumGrads_.resize(gCellSize, Point());

  curCoordi_.resize(gCellSize, Point());
  nextCoordi_.resize(gCellSize, Point());

  for(int i = 0; i < gCellSize; i++)
  {
    nb_->updateDensityCoordiLayoutInside(nb_->gCells()[i]);
    curSLPCoordi_[i] 
      = prevSLPCoordi_[i] 
      = curCoordi_[i] 
      = Point(nb_->gCells()[i]->cx(), nb_->gCells()[i]->cy());
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
 
  // fill in curSLPSumGrads_, curSLPWireLengthGrads_, curSLPDensityGrads_ 
  updateGradients(
      curSLPSumGrads_, curSLPWireLengthGrads_,
      curSLPDensityGrads_);

  if( isDiverged_ ) {
    return;
  }

  // approximately fill in 
  // prevSLPCoordi_ to calculate lc vars
  updateInitialPrevSLPCoordi();

  // bin, FFT, wlen update with prevSLPCoordi.
  nb_->updateGCellCenterLocation(prevSLPCoordi_);
  nb_->updateDensityForceBin();
  nb_->updateWireLengthForceWA(wireLengthCoefX_, wireLengthCoefY_);
  
  // update previSumGrads_, prevSLPWireLengthGrads_, prevSLPDensityGrads_
  updateGradients(
      prevSLPSumGrads_, prevSLPWireLengthGrads_,
      prevSLPDensityGrads_);
  
  if( isDiverged_ ) {
    return;
  }
  
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
void
NesterovPlace::updateGradients(
    std::vector<Point>& sumGrads,
    std::vector<Point>& wireLengthGrads,
    std::vector<Point>& densityGrads) {

  wireLengthGradSum_ = 0;
  densityGradSum_ = 0;

  prec gradSum = 0;

  LOG_DEBUG("Density Penalty: {}", densityPenalty_);

  for(size_t i=0; i<nb_->gCells().size(); i++) {
    GCell* gCell = nb_->gCells().at(i);
    wireLengthGrads[i] = nb_->getWireLengthGradientWA(
        gCell, wireLengthCoefX_, wireLengthCoefY_);
    densityGrads[i] = nb_->getDensityGradient(gCell); 

    // Different compiler has different results on the following formula.
    // e.g. wireLengthGradSum_ += fabs(~~.x) + fabs(~~.y);
    //
    // To prevent instability problem,
    // I partitioned the fabs(~~.x) + fabs(~~.y) as two terms.
    //
    wireLengthGradSum_ += fabs(wireLengthGrads[i].x);
    wireLengthGradSum_ += fabs(wireLengthGrads[i].y);
      
    densityGradSum_ += fabs(densityGrads[i].x);
    densityGradSum_ += fabs(densityGrads[i].y);

    sumGrads[i].x = wireLengthGrads[i].x + densityPenalty_ * densityGrads[i].x;
    sumGrads[i].y = wireLengthGrads[i].y + densityPenalty_ * densityGrads[i].y;

    Point wireLengthPreCondi 
      = nb_->getWireLengthPreconditioner(gCell);
    Point densityPrecondi
      = nb_->getDensityPreconditioner(gCell);

    Point sumPrecondi(
        wireLengthPreCondi.x + densityPenalty_ * densityPrecondi.x,
        wireLengthPreCondi.y + densityPenalty_ * densityPrecondi.y);

    if( sumPrecondi.x <= npVars_.minPreconditioner ) {
      sumPrecondi.x = npVars_.minPreconditioner;
    }

    if( sumPrecondi.y <= npVars_.minPreconditioner ) {
      sumPrecondi.y = npVars_.minPreconditioner; 
    }
    
    sumGrads[i].x /= sumPrecondi.x;
    sumGrads[i].y /= sumPrecondi.y; 

    gradSum += fabs(sumGrads[i].x) + fabs(sumGrads[i].y);
  }
  
  LOG_DEBUG("WireLengthGradSum: {}", wireLengthGradSum_);
  LOG_DEBUG("DensityGradSum: {}", densityGradSum_);
  LOG_DEBUG("GradSum: {}", gradSum);

  // divergence detection on 
  // Wirelength / density gradient calculation
  if( isnan(wireLengthGradSum_) || isinf(wireLengthGradSum_) ||
      isnan(densityGradSum_) || isinf(densityGradSum_) ) {
    isDiverged_ = true;
  }
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
  Plot::plot(nb_.get(), PlotNesterovType::Bin, "./plot/bin", placename + "bin_0");
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

      updateGradients(nextSLPSumGrads_, nextSLPWireLengthGrads_, nextSLPDensityGrads_ );

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

    if( i == 0 || (i+1) % 100 == 0 ) {
      LOG_DEBUG("[NesterovSolve] Iter: {} overflow: {} HPWL: {}", i+1, sumOverflow_, prevHpwl_);
      Plot::plot(nb_.get(), PlotNesterovType::GCell, "./plot/cell", placename + "cell_" + std::to_string(i+1));
      Plot::plot(nb_.get(), PlotNesterovType::Bin, "./plot/bin", placename + "bin_" + std::to_string(i+1));
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
      cout << "[NesterovSolve] Finished with Overflow: " << sumOverflow_ << endl;
      break;
    }
  }
 
  // in all case including diverge, 
  // PlacerBase should be updated. 
  updatePlacerBase();

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
  std::swap(prevSLPWireLengthGrads_, curSLPWireLengthGrads_);
  std::swap(prevSLPDensityGrads_, curSLPDensityGrads_);
  std::swap(prevSLPSumGrads_, curSLPSumGrads_);
  
  std::swap(curSLPCoordi_, nextSLPCoordi_);
  std::swap(curSLPWireLengthGrads_, nextSLPWireLengthGrads_);
  std::swap(curSLPDensityGrads_, nextSLPDensityGrads_);
  std::swap(curSLPSumGrads_, nextSLPSumGrads_);

  std::swap(curCoordi_, nextCoordi_);

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

static prec
getDistance(vector<Point>& a, vector<Point>& b) {
  prec sumDistance = 0.0f;
  for(size_t i=0; i<a.size(); i++) {
    sumDistance += (a[i].x - b[i].x) * (a[i].x - b[i].x);
    sumDistance += (a[i].y - b[i].y) * (a[i].y - b[i].y);
  }

  return sqrt( sumDistance / (2.0 * a.size()) );
}

static prec
getSecondNorm(vector<Point>& a) {
  prec norm = 0;
  for(auto& coordi : a) {
    norm += coordi.x * coordi.x + coordi.y * coordi.y;
  }
  return sqrt( norm / (2.0*a.size()) ); 
}


}