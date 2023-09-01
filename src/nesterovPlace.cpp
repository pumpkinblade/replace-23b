#include "placerBase.h"
#include "nesterovBase.h"
#include "nesterovPlace.h"
#include "macroLegalizer.h"
#include "log.h"
#include "plot.h"
#include <glpk.h>
#include <random>

using namespace std;

namespace replace
{
  static double getDistance(vector<Point> &a, vector<Point> &b);
  static double getSecondNorm(vector<Point> &a);

  static void getRotBox(const GCell *cell, bool rot, double &lx, double &ux, double &ly, double &uy);
  static double getRotOverlap(const GCell *cell1, const GCell *cell2, bool rot1, bool rot2);
  static void ilpSolve(const std::vector<GCell *> &macros, std::vector<bool> &isRot);
  static Orientation findNearestOrientation(double theta);
  static double getThetaInsideTwoPi(double theta);

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
        useTheta(false)
  {
  }

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
        isDiverged_(false)
  {
  }

  NesterovPlace::NesterovPlace(NesterovPlaceVars npVars, std::shared_ptr<NesterovBase> nb)
      : NesterovPlace()
  {
    npVars_ = npVars;
    nb_ = nb;
  }

  void NesterovPlace::init()
  {
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

    if (npVars_.useLocalDensity)
    {
      localAlpha_ = npVars_.initAlpha;
      localBeta_ = npVars_.initBeta;
      prevCellDelta_.resize(gCellSize, 0.0);
      curCellDelta_.resize(gCellSize, 0.0);
      nextCellDelta_.resize(gCellSize, 0.0);
    }

    useThetaMacroCount = 0;
    if (npVars_.useTheta)
    {
      while (useThetaMacroCount < gCellSize)
      {
        if (nb_->gCells()[useThetaMacroCount]->isMacro())
          useThetaMacroCount++;
        else
          break;
      }
      curTheta_.resize(useThetaMacroCount, 0.0);
      nextTheta_.resize(useThetaMacroCount, 0.0);
      curSLPTheta_.resize(useThetaMacroCount, 0.0);
      curSLPSumGradTheta_.resize(useThetaMacroCount, 0.0);
      nextSLPTheta_.resize(useThetaMacroCount, 0.0);
      nextSLPSumGradTheta_.resize(useThetaMacroCount, 0.0);
      prevSLPTheta_.resize(useThetaMacroCount, 0.0);
      prevSLPSumGradTheta_.resize(useThetaMacroCount, 0.0);
      // precalculate the preconditional of theta
      wireLengthPrecondiTheta_.resize(useThetaMacroCount, 0.0);
      densityPrecondiTheta_.resize(useThetaMacroCount, 0.0);
      for (int i = 0; i < useThetaMacroCount; i++)
      {
        double wireLengthPrecondi = nb_->getWireLengthPreconditionerTheta(nb_->gCells()[i]);
        double densityPrecondi = nb_->getDensityPreconditionerTheta(nb_->gCells()[i]);
        wireLengthPrecondiTheta_[i] = wireLengthPrecondi;
        densityPrecondiTheta_[i] = densityPrecondi;
        curSLPTheta_[i] = prevSLPTheta_[i] = curTheta_[i] = getThetaInsideTwoPi(nb_->gCells()[i]->theta());
      }
    }

    for (int i = 0; i < gCellSize; i++)
    {
      bool useTheta = npVars_.useTheta && nb_->gCells()[i]->isMacro();
      nb_->updateDensityCoordiLayoutInside(nb_->gCells()[i], useTheta);
      curSLPCoordi_[i] = prevSLPCoordi_[i] = curCoordi_[i] = Point(nb_->gCells()[i]->cx(), nb_->gCells()[i]->cy());
    }

    // bin update
    if (npVars_.useTheta)
      nb_->updateGCellCenterLocationWithTheta(curSLPCoordi_, curSLPTheta_);
    else
      nb_->updateGCellCenterLocation(curSLPCoordi_);

    // get hpwl
    nb_->updateNetsBox();
    prevHpwl_ = nb_->getHpwl();
    LOG_DEBUG("InitialHPWL: {}", prevHpwl_);

    // FFT update
    nb_->updateDensityForceBin();

    double avgBinSizeX = 0.f;
    double avgBinSizeY = 0.f;
    for (BinGrid *bg : nb_->binGrids())
    {
      avgBinSizeX += bg->binSizeX();
      avgBinSizeY += bg->binSizeY();
    }
    avgBinSizeX /= static_cast<double>(nb_->binGrids().size());
    avgBinSizeY /= static_cast<double>(nb_->binGrids().size());
    baseWireLengthCoef_ = npVars_.initWireLengthCoef / (0.5f * (avgBinSizeX + avgBinSizeY));
    LOG_DEBUG("BaseWireLengthCoef: {}", baseWireLengthCoef_);

    sumOverflow_ = nb_->getOverflow();
    LOG_DEBUG("InitSumOverflow: {}", sumOverflow_);
    updateWireLengthCoef(sumOverflow_);
    nb_->updateWireLengthForceWA(wireLengthCoefX_, wireLengthCoefY_);
    LOG_DEBUG("WireLengthCoef: {}, {}", wireLengthCoefX_, wireLengthCoefY_);

    // fill in curSLPSumGrads_
    updateGradients(curSLPSumGrads_, curSLPSumGradTheta_, curCellDelta_, curCellDelta_);

    if (isDiverged_)
      return;

    // approximately fill in prevSLPCoordi_ to calculate lc vars
    updateInitialPrevSLPCoordi();

    // bin, FFT, wlen update with prevSLPCoordi.
    if (npVars_.useTheta)
      nb_->updateGCellCenterLocationWithTheta(prevSLPCoordi_, prevSLPTheta_);
    else
      nb_->updateGCellCenterLocation(prevSLPCoordi_);
    nb_->updateDensityForceBin();
    nb_->updateWireLengthForceWA(wireLengthCoefX_, wireLengthCoefY_);

    // update previSumGrads_
    updateGradients(prevSLPSumGrads_, prevSLPSumGradTheta_, prevCellDelta_, prevCellDelta_);

    if (isDiverged_)
      return;

    LOG_DEBUG("WireLengthGradSum: {}", wireLengthGradSum_);
    LOG_DEBUG("DensityGradSum", densityGradSum_);

    densityPenalty_ = (wireLengthGradSum_ / densityGradSum_) * npVars_.initDensityPenalty;
    LOG_DEBUG("InitDensityPenalty: {}", densityPenalty_);

    sumOverflow_ = nb_->getOverflow();
    LOG_DEBUG("PrevSumOverflow: {}", sumOverflow_);

    stepLength_ = getStepLength(prevSLPCoordi_, prevSLPSumGrads_, curSLPCoordi_, curSLPSumGrads_);
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
  void NesterovPlace::updateGradients(
      std::vector<Point> &sumGrads,
      std::vector<double> &sumGradTheta,
      std::vector<double> &refCellDeltas,
      std::vector<double> &outCellDeltas)
  {
    wireLengthGradSum_ = 0;
    densityGradSum_ = 0;
    localDensityGradSum_ = 0;
    wireLengthGradSumTheta_ = 0;
    densityGradSumTheta_ = 0;
    localDensityGradSumTheta_ = 0;
    double gradSum = 0;
    double gradSumTheta = 0;
    double cellDeltaSum = 0;

    LOG_DEBUG("Density Penalty: {}", densityPenalty_);

    for (int i = 0; i < useThetaMacroCount; i++)
    {
      GCell *gCell = nb_->gCells().at(i);

      // wirelength
      Point wireLengthGrad;
      double wireLengthGradTheta;
      wireLengthGrad = nb_->getWireLengthGradientWAWithTheta(gCell, wireLengthCoefX_, wireLengthCoefY_, wireLengthGradTheta);
      wireLengthGradSumTheta_ += fabs(wireLengthGradTheta);
      wireLengthGradSum_ += fabs(wireLengthGrad.x);
      wireLengthGradSum_ += fabs(wireLengthGrad.y);

      // density
      Point densityGrad;
      double densityGradTheta;
      densityGrad = nb_->getDensityGradientWithTheta(gCell, densityGradTheta);
      densityGradSumTheta_ += fabs(densityGradTheta);
      densityGradSum_ += fabs(densityGrad.x);
      densityGradSum_ += fabs(densityGrad.y);

      // local density
      Point localDensityGrad;
      double localDensityGradTheta = 0;
      double delta;
      if (npVars_.useLocalDensity)
      {
        delta = refCellDeltas[i];
        localDensityGrad = nb_->getLocalDensityGradientWithTheta(gCell, localAlpha_, localBeta_, delta, localDensityGradTheta);
        localDensityGradSumTheta_ += fabs(localDensityGradTheta);
        localDensityGradSum_ += fabs(localDensityGrad.x);
        localDensityGradSum_ += fabs(localDensityGrad.y);
        cellDeltaSum += delta;
      }

      // sumGrad
      if (npVars_.useLocalDensity)
      {
        sumGrads[i].x = wireLengthGrad.x + densityPenalty_ * densityGrad.x + delta * localDensityGrad.x;
        sumGrads[i].y = wireLengthGrad.y + densityPenalty_ * densityGrad.y + delta * localDensityGrad.y;
        sumGradTheta[i] = wireLengthGradTheta + densityPenalty_ * densityGradTheta + delta * localDensityGradTheta;
        outCellDeltas[i] = delta;
      }
      else
      {
        sumGrads[i].x = wireLengthGrad.x + densityPenalty_ * densityGrad.x;
        sumGrads[i].y = wireLengthGrad.y + densityPenalty_ * densityGrad.y;
        sumGradTheta[i] = wireLengthGradTheta + densityPenalty_ * densityGradTheta;
      }

      // xy precondi
      Point wireLengthPrecondi = nb_->getWireLengthPreconditioner(gCell);
      Point densityPrecondi = nb_->getDensityPreconditioner(gCell);
      Point sumPrecondi(
          wireLengthPrecondi.x + densityPenalty_ * densityPrecondi.x,
          wireLengthPrecondi.y + densityPenalty_ * densityPrecondi.y);
      if (npVars_.useLocalDensity)
      {
        Point localDensityPrecondi = nb_->getLocalDensityPreconditioner(gCell);
        sumPrecondi.x += delta * localDensityPrecondi.x;
        sumPrecondi.y += delta * localDensityPrecondi.y;
      }
      sumPrecondi.x = std::max(sumPrecondi.x, npVars_.minPreconditioner);
      sumPrecondi.y = std::max(sumPrecondi.y, npVars_.minPreconditioner);
      sumGrads[i].x /= sumPrecondi.x;
      sumGrads[i].y /= sumPrecondi.y;
      gradSum += fabs(sumGrads[i].x) + fabs(sumGrads[i].y);

      // theta precondi
      double wireLenghtPrecondiTheta = wireLengthPrecondiTheta_[i];
      double densityPrecondiTheta = densityPrecondiTheta_[i];
      double sumPrecondiTheta = wireLenghtPrecondiTheta + densityPenalty_ * densityPrecondiTheta;
      if (npVars_.useLocalDensity)
      {
        sumPrecondiTheta += delta * nb_->getLocalDensityPreconditionerTheta(gCell);
      }
      sumPrecondiTheta = std::max(sumPrecondiTheta, npVars_.minPreconditioner);
      sumGradTheta[i] /= sumPrecondiTheta;
      gradSumTheta += fabs(sumGradTheta[i]);
    }

    for (int i = useThetaMacroCount; i < nb_->gCells().size(); i++)
    {
      GCell *gCell = nb_->gCells().at(i);

      // wirelength
      Point wireLengthGrad;
      wireLengthGrad = nb_->getWireLengthGradientWA(gCell, wireLengthCoefX_, wireLengthCoefY_);
      wireLengthGradSum_ += fabs(wireLengthGrad.x);
      wireLengthGradSum_ += fabs(wireLengthGrad.y);

      // density
      Point densityGrad;
      densityGrad = nb_->getDensityGradient(gCell);
      densityGradSum_ += fabs(densityGrad.x);
      densityGradSum_ += fabs(densityGrad.y);

      // local density
      Point localDensityGrad;
      double delta;
      if (npVars_.useLocalDensity)
      {
        delta = refCellDeltas[i];
        localDensityGrad = nb_->getLocalDensityGradient(gCell, localAlpha_, localBeta_, delta);
        localDensityGradSum_ += fabs(localDensityGrad.x);
        localDensityGradSum_ += fabs(localDensityGrad.y);
        cellDeltaSum += delta;
      }

      if (npVars_.useLocalDensity)
      {
        sumGrads[i].x = wireLengthGrad.x + densityPenalty_ * densityGrad.x + delta * localDensityGrad.x;
        sumGrads[i].y = wireLengthGrad.y + densityPenalty_ * densityGrad.y + delta * localDensityGrad.y;
        outCellDeltas[i] = delta;
      }
      else
      {
        sumGrads[i].x = wireLengthGrad.x + densityPenalty_ * densityGrad.x;
        sumGrads[i].y = wireLengthGrad.y + densityPenalty_ * densityGrad.y;
      }

      Point wireLengthPrecondi = nb_->getWireLengthPreconditioner(gCell);
      Point densityPrecondi = nb_->getDensityPreconditioner(gCell);
      Point sumPrecondi(
          wireLengthPrecondi.x + densityPenalty_ * densityPrecondi.x,
          wireLengthPrecondi.y + densityPenalty_ * densityPrecondi.y);
      sumPrecondi.x = std::max(sumPrecondi.x, npVars_.minPreconditioner);
      sumPrecondi.y = std::max(sumPrecondi.y, npVars_.minPreconditioner);
      if (npVars_.useLocalDensity)
      {
        Point localDensityPrecondi = nb_->getLocalDensityPreconditioner(gCell);
        sumPrecondi.x += delta * localDensityPrecondi.x;
        sumPrecondi.y += delta * localDensityPrecondi.y;
      }
      sumGrads[i].x /= sumPrecondi.x;
      sumGrads[i].y /= sumPrecondi.y;
      gradSum += fabs(sumGrads[i].x) + fabs(sumGrads[i].y);
    }

    if (npVars_.useLocalDensity)
    {
      cellDeltaSum /= nb_->gCells().size();
      LOG_DEBUG("CellDeltaAverage: {}", cellDeltaSum);
    }

    LOG_DEBUG("WireLengthGradSum: {}", wireLengthGradSum_);
    LOG_DEBUG("DensityGradSum: {}", densityGradSum_);
    if (npVars_.useLocalDensity)
      LOG_DEBUG("LocalDensityGradSum: {}", localDensityGradSum_);
    LOG_DEBUG("GradSum: {}", gradSum);
    if (npVars_.useTheta)
    {
      LOG_DEBUG("WireLengthGradSumTheta: {}", wireLengthGradSumTheta_);
      LOG_DEBUG("DensityGradSumTheta: {}", densityGradSumTheta_);
      if (npVars_.useLocalDensity)
        LOG_DEBUG("LocalDensityGradSumTheta: {}", localDensityGradSumTheta_);
      LOG_DEBUG("GradSumTheta: {}", gradSumTheta);
    }

    // divergence detection
    isDiverged_ = (std::isnan(gradSum) || std::isinf(gradSum) || std::isnan(gradSumTheta) || std::isinf(gradSumTheta));
  }

  void NesterovPlace::doNesterovPlace(string placename)
  {
    LOG_TRACE("start NesterovPlace::doNesterovPlace");

    // if replace diverged in init() function,
    // replace must be skipped.
    if (isDiverged_)
    {
      LOG_ERROR("RePlAce diverged. Please tune the parameters again");
      return;
    }

    if (placename != "" && placename[placename.length() - 1] != '_')
    {
      placename.push_back('_');
    }
    Plot::plot(nb_.get(), PlotNesterovType::GCell, "./plot/cell", placename + "cell_0");
    Plot::plot(nb_.get(), PlotNesterovType::Arrow, "./plot/arrow", placename + "arrow_0");

    // backTracking variable.
    double curA = 1.0;

    // divergence detection
    double minSumOverflow = 1e30;
    double hpwlWithMinSumOverflow = 1e30;

    // dynamic adjustment of max_phi_coef
    bool isMaxPhiCoefChanged = false;

    // diverge error handling
    string divergeMsg = "";
    int divergeCode = 0;

    // Core Nesterov Loop
    for (int i = 0; i < npVars_.maxNesterovIter; i++)
    {
      LOG_DEBUG("Iter: {}", i + 1);

      double prevA = curA;

      // here, prevA is a_(k), curA is a_(k+1)
      // See, the papers' Algorithm 4 section
      //
      curA = (1.0 + sqrt(4.0 * prevA * prevA + 1.0)) * 0.5;

      // coeff is (a_k -1) / ( a_(k+1)) in paper.
      double coeff = (prevA - 1.0) / curA;

      LOG_DEBUG("PreviousA: {}", prevA);
      LOG_DEBUG("CurrentA: {}", curA);
      LOG_DEBUG("Coefficient: {}", coeff);
      LOG_DEBUG("StepLength: {}", stepLength_);

      // Back-Tracking loop
      int numBackTrak = 0;
      for (numBackTrak = 0; numBackTrak < npVars_.maxBackTrack; numBackTrak++)
      {

        // fill in nextCoordinates with given stepLength_
        for (int k = 0; k < nb_->gCells().size(); k++)
        {
          GCell *curGCell = nb_->gCells()[k];

          Point nextCoordi(
              curSLPCoordi_[k].x + stepLength_ * curSLPSumGrads_[k].x,
              curSLPCoordi_[k].y + stepLength_ * curSLPSumGrads_[k].y);

          Point nextSLPCoordi(
              nextCoordi.x + coeff * (nextCoordi.x - curCoordi_[k].x),
              nextCoordi.y + coeff * (nextCoordi.y - curCoordi_[k].y));

          nextCoordi_[k] = Point(
              nb_->getDensityCoordiLayoutInsideX(
                  curGCell, nextCoordi.x, npVars_.useTheta && curGCell->isMacro()),
              nb_->getDensityCoordiLayoutInsideY(
                  curGCell, nextCoordi.y, npVars_.useTheta && curGCell->isMacro()));

          nextSLPCoordi_[k] = Point(
              nb_->getDensityCoordiLayoutInsideX(
                  curGCell, nextSLPCoordi.x, npVars_.useTheta && curGCell->isMacro()),
              nb_->getDensityCoordiLayoutInsideY(
                  curGCell, nextSLPCoordi.y, npVars_.useTheta && curGCell->isMacro()));
        }

        for (int k = 0; k < useThetaMacroCount; k++)
        {
          double nextTheta = curSLPTheta_[k] + stepLength_ * curSLPSumGradTheta_[k];
          double nextSLPTheta = nextTheta + coeff * (nextTheta - curTheta_[k]);
          nextTheta_[k] = getThetaInsideTwoPi(nextTheta);
          nextSLPTheta_[k] = getThetaInsideTwoPi(nextSLPTheta);
        }

        if (npVars_.useTheta)
          nb_->updateGCellCenterLocationWithTheta(nextSLPCoordi_, nextTheta_);
        else
          nb_->updateGCellCenterLocation(nextSLPCoordi_);
        nb_->updateDensityForceBin();
        nb_->updateWireLengthForceWA(wireLengthCoefX_, wireLengthCoefY_);

        updateGradients(nextSLPSumGrads_, nextSLPSumGradTheta_, curCellDelta_, nextCellDelta_);

        // NaN or inf is detected in WireLength/Density Coef
        if (isDiverged_)
        {
          break;
        }

        double newStepLength = getStepLength(curSLPCoordi_, curSLPSumGrads_, nextSLPCoordi_, nextSLPSumGrads_);

        LOG_DEBUG("NetStepLength: {}", newStepLength);

        if (newStepLength > stepLength_ * 0.95)
        {
          stepLength_ = newStepLength;
          break;
        }
        else
        {
          stepLength_ = newStepLength;
        }
      }

      LOG_DEBUG("NumBackTrak: {}", numBackTrak + 1);

      // dynamic adjustment for
      // better convergence with
      // large designs
      if (!isMaxPhiCoefChanged && sumOverflow_ < 0.35f)
      {
        isMaxPhiCoefChanged = true;
        npVars_.maxPhiCoef *= 0.99;
      }

      // usually, maxBackTrack should be 1~3
      // 10 is the case when
      // all of cells are not moved at all.
      if (npVars_.maxBackTrack == numBackTrak)
      {
        divergeMsg = "RePlAce divergence detected. \n";
        divergeMsg += "        Please decrease init_density_penalty value";
        divergeCode = 3;
        isDiverged_ = true;
      }

      if (isDiverged_)
      {
        break;
      }

      updateNextIter();

      // For JPEG Saving
      // debug

      if (i == 0 || (i + 1) % 20 == 0)
      {
        LOG_DEBUG("[NesterovSolve] Iter: {} overflow: {} HPWL: {}", i + 1, sumOverflow_, prevHpwl_);
        Plot::plot(nb_.get(), PlotNesterovType::GCell, "./plot/cell", placename + "cell_" + std::to_string(i + 1));
        Plot::plot(nb_.get(), PlotNesterovType::Arrow, "./plot/arrow", placename + "arrow_" + std::to_string(i + 1));
      }

      if (minSumOverflow > sumOverflow_)
      {
        minSumOverflow = sumOverflow_;
        hpwlWithMinSumOverflow = prevHpwl_;
      }

      // diverge detection on
      // large max_phi_cof value + large design
      //
      // 1) happen overflow < 20%
      // 2) Hpwl is growing
      //
      if (sumOverflow_ < 0.3f && sumOverflow_ - minSumOverflow >= 0.02f && hpwlWithMinSumOverflow * 1.2f < prevHpwl_)
      {
        divergeMsg = "RePlAce divergence detected. \n";
        divergeMsg += "        Please decrease max_phi_cof value";
        divergeCode = 4;
        isDiverged_ = true;
        break;
      }

      // minimum iteration is 50
      if (i > 50 && sumOverflow_ <= npVars_.targetOverflow)
      {
        LOG_DEBUG("[NesterovSolve] Finished with Overflow: {}", sumOverflow_);
        break;
      }
    }

    // in all case including diverge,
    // PlacerBase should be updated.
    updatePlacerBase();
    Plot::plot(nb_.get(), PlotNesterovType::GCell, "./plot/cell", placename + "cell_end");
    Plot::plot(nb_.get(), PlotNesterovType::Arrow, "./plot/arrow", placename + "arrow_end");

    if (isDiverged_)
    {
      LOG_ERROR("{} : Code `{}`", divergeMsg, divergeCode);
    }
  }

  void NesterovPlace::updateWireLengthCoef(double overflow)
  {
    if (overflow > 1.0)
    {
      wireLengthCoefX_ = wireLengthCoefY_ = 0.1;
    }
    else if (overflow < 0.1)
    {
      wireLengthCoefX_ = wireLengthCoefY_ = 10.0;
    }
    else
    {
      wireLengthCoefX_ = wireLengthCoefY_ = 1.0 / pow(10.0, (overflow - 0.1) * 20 / 9.0 - 1.0);
    }

    wireLengthCoefX_ *= baseWireLengthCoef_;
    wireLengthCoefY_ *= baseWireLengthCoef_;
    LOG_DEBUG("NewWireLengthCoef: {}", wireLengthCoefX_);
  }

  void NesterovPlace::updateInitialPrevSLPCoordi()
  {
    for (int i = 0; i < useThetaMacroCount; i++)
    {
      GCell *curGCell = nb_->gCells()[i];
      double prevCoordiX = curSLPCoordi_[i].x + npVars_.initialPrevCoordiUpdateCoef * curSLPSumGrads_[i].x;
      double prevCoordiY = curSLPCoordi_[i].y + npVars_.initialPrevCoordiUpdateCoef * curSLPSumGrads_[i].y;
      Point newCoordi(
          nb_->getDensityCoordiLayoutInsideX(curGCell, prevCoordiX, true),
          nb_->getDensityCoordiLayoutInsideY(curGCell, prevCoordiY, true));
      prevSLPCoordi_[i] = newCoordi;

      double prevTheta = curSLPTheta_[i] + npVars_.initialPrevCoordiUpdateCoef * curSLPSumGradTheta_[i];
      prevSLPTheta_[i] = getThetaInsideTwoPi(prevTheta);
    }

    for (int i = useThetaMacroCount; i < nb_->gCells().size(); i++)
    {
      GCell *curGCell = nb_->gCells()[i];
      double prevCoordiX = curSLPCoordi_[i].x + npVars_.initialPrevCoordiUpdateCoef * curSLPSumGrads_[i].x;
      double prevCoordiY = curSLPCoordi_[i].y + npVars_.initialPrevCoordiUpdateCoef * curSLPSumGrads_[i].y;
      Point newCoordi(
          nb_->getDensityCoordiLayoutInsideX(curGCell, prevCoordiX, false),
          nb_->getDensityCoordiLayoutInsideY(curGCell, prevCoordiY, false));
      prevSLPCoordi_[i] = newCoordi;
    }
  }

  void NesterovPlace::updateNextIter()
  {
    // swap vector pointers
    std::swap(prevSLPCoordi_, curSLPCoordi_);
    std::swap(prevSLPSumGrads_, curSLPSumGrads_);
    std::swap(curSLPCoordi_, nextSLPCoordi_);
    std::swap(curSLPSumGrads_, nextSLPSumGrads_);
    std::swap(curCoordi_, nextCoordi_);

    if (npVars_.useLocalDensity)
    {
      std::swap(prevCellDelta_, curCellDelta_);
      std::swap(curCellDelta_, nextCellDelta_);
    }

    if (npVars_.useTheta)
    {
      std::swap(prevSLPTheta_, curSLPTheta_);
      std::swap(prevSLPSumGradTheta_, curSLPSumGradTheta_);
      std::swap(curSLPTheta_, nextSLPTheta_);
      std::swap(curSLPSumGradTheta_, nextSLPSumGradTheta_);
    }

    sumOverflow_ = nb_->getOverflow();

    LOG_DEBUG("Gradient: {}", getSecondNorm(curSLPSumGrads_));
    LOG_DEBUG("Phi: {}", nb_->sumPhi());
    LOG_DEBUG("Overflow: {}", sumOverflow_);

    updateWireLengthCoef(sumOverflow_);
    nb_->updateNetsBox();
    double hpwl = nb_->getHpwl();

    LOG_DEBUG("PreviousHPWL: {}", prevHpwl_);
    LOG_DEBUG("NewHPWL: {}", hpwl);

    double phiCoef = getPhiCoef(static_cast<double>(hpwl - prevHpwl_) / npVars_.referenceHpwl);

    prevHpwl_ = hpwl;
    densityPenalty_ *= phiCoef;

    if (npVars_.useLocalDensity)
    {
      localAlpha_ *= phiCoef;
      localBeta_ *= phiCoef;
    }

    LOG_DEBUG("PhiCoef: {}", phiCoef);
  }

  double NesterovPlace::getStepLength(
      std::vector<Point> &prevSLPCoordi_,
      std::vector<Point> &prevSLPSumGrads_,
      std::vector<Point> &curSLPCoordi_,
      std::vector<Point> &curSLPSumGrads_)
  {

    double coordiDistance = getDistance(prevSLPCoordi_, curSLPCoordi_);
    double gradDistance = getDistance(prevSLPSumGrads_, curSLPSumGrads_);

    LOG_DEBUG("CoordinateDistance: {}", coordiDistance);
    LOG_DEBUG("GradientDistance: {}", gradDistance);

    return coordiDistance / gradDistance;
  }

  double NesterovPlace::getPhiCoef(double scaledDiffHpwl)
  {
    LOG_DEBUG("Input ScaleDiffHPWL", scaledDiffHpwl);

    double retCoef = (scaledDiffHpwl < 0) ? npVars_.maxPhiCoef : npVars_.maxPhiCoef * pow(npVars_.maxPhiCoef, scaledDiffHpwl * -1.0);
    retCoef = std::max(npVars_.minPhiCoef, retCoef);
    return retCoef;
  }

  void NesterovPlace::updatePlacerBase()
  {
    LOG_TRACE("Begin updatePlacerBase");
    if (npVars_.useTheta)
      determinMacroOrient();
    LOG_TRACE("End updatePlacerBase");
    for (auto &gCell : nb_->gCells())
    {
      if (gCell->isInstance())
      {
        int cx = static_cast<int>(std::round(gCell->cx()));
        int cy = static_cast<int>(std::round(gCell->cy()));
        gCell->instance()->setCenterLocation(cx, cy);
      }
    }
  }

  void NesterovPlace::determinMacroOrient()
  {
    for (BinGrid *bg : nb_->binGrids())
    {
      const Technology &tech = *bg->die()->tech();
      std::vector<GCell*> macros, macrosILP;
      std::vector<bool> rot, rotILP;

      for (int i = 0; i < useThetaMacroCount; i++)
      {
        GCell *gcell = nb_->gCells()[i];
        if (gcell->binGrid() != bg)
          continue;
        Orientation ori = findNearestOrientation(gcell->theta());
        int idx = static_cast<int>(ori);
        if (std::abs(LEGAL_THETA[idx] - gcell->theta()) < 0.1 * PI)
        {
          macros.push_back(gcell);
          rot.push_back((ori == Orientation::R90) || (ori == Orientation::R270));
        }
        else
        {
          macros.push_back(gcell);
          rot.push_back(false);
          macrosILP.push_back(gcell);
          rotILP.push_back(false);
        }
      }

      LOG_TRACE("Begin ILP Solve");
      ilpSolve(macrosILP, rotILP);
      LOG_TRACE("Finish ILP Solve");

      LOG_DEBUG("macros.size() = {}, rot.size() = {}, macrosILP.size = {}, rotILP.size() = {}", 
                macros.size(), rot.size(), macrosILP.size(), rotILP.size());
      for(int i = 0, j = 0; i < rot.size() && j < rotILP.size(); i++)
      {
        LOG_DEBUG("i = {}, j = {}", i, j);
        if (macros[i] == macrosILP[j])
        {
          rot[i] = rotILP[j];
          ++j;
        }
      }

      MacroLegalizer lg;
      std::mt19937 rng(114514);
      std::bernoulli_distribution bn(0.2);
      const int iter = 100000;
      std::vector<std::pair<Orientation, Orientation>> preferOri(macros.size());
      std::vector<Instance*> macroInstances(macros.size());
      for(int i = 0; i < macros.size(); i++)
      {
        if(std::abs(macros[i]->theta()) < std::abs(macros[i]->theta() - PI) || 
          std::abs(macros[i]->theta() - 2 * PI) < std::abs(macros[i]->theta() - PI))
          preferOri[i].first = Orientation::R0;
        else
          preferOri[i].first = Orientation::R180;

        if(std::abs(macros[i]->theta() - 0.5 * PI) < std::abs(macros[i]->theta() - 1.5 * PI))
          preferOri[i].second = Orientation::R90;
        else
          preferOri[i].second = Orientation::R270;

        macroInstances[i] = macros[i]->instance();
      }

      for(int i = 0; i < iter; i++)
      {
        for(int i = 0; i < rot.size(); i++)
        {
          // set macro's location by gcell's location
          int cx = static_cast<int>(macros[i]->cx());
          int cy = static_cast<int>(macros[i]->cy());
          macroInstances[i]->setCenterLocation(cx, cy);
          
          // apply rotation
          if(rot[i])
            macroInstances[i]->setOrientSize(tech, preferOri[i].first);
          else
            macroInstances[i]->setOrientSize(tech, preferOri[i].second);
        }
        
        // try do post legalization
        bool ok = lg.postLegalize(macroInstances, bg->die());
        if(ok)
          break;
        else
        {
          // each macro has 20% probability to rotate
          for(int i = 0; i < rot.size(); i++)
            rot[i] = bn(rng) ^ rot[i];
        }
      }
    }
  }

  static double getDistance(vector<Point> &a, vector<Point> &b)
  {
    double sumDistance = 0.0f;
    for (size_t i = 0; i < a.size(); i++)
    {
      sumDistance += (a[i].x - b[i].x) * (a[i].x - b[i].x);
      sumDistance += (a[i].y - b[i].y) * (a[i].y - b[i].y);
    }

    return sqrt(sumDistance / (2.0 * a.size()));
  }

  static double getSecondNorm(vector<Point> &a)
  {
    double norm = 0;
    for (auto &coordi : a)
      norm += coordi.x * coordi.x + coordi.y * coordi.y;
    return sqrt(norm / (2.0 * a.size()));
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

  static double getRotOverlap(const GCell *cell1, const GCell *cell2, bool rot1, bool rot2)
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

  static void ilpSolve(const std::vector<GCell *> &macros, std::vector<bool> &isRot)
  {
    int numMacros = static_cast<int>(macros.size());
    if (numMacros == 0)
      return;
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

  static Orientation findNearestOrientation(double theta)
  {
    double bestDist = std::abs(theta);
    Orientation bestOri = Orientation::R0;
    for (int idx = 1; idx < 5; idx++)
    {
      double dist = std::abs(theta - LEGAL_THETA[idx]);
      if (dist < bestDist)
      {
        bestDist = dist;
        bestOri = static_cast<Orientation>(idx);
      }
    }
    return bestOri;
  }

  static double getThetaInsideTwoPi(double theta)
  {
    return std::fmod(std::fmod(theta, 2.0f * PI) + 2.0f * PI, 2.0f * PI);
  }

}
