#include "replace.h"
#include "initialPlace.h"
#include "nesterovPlace.h"
#include "placerBase.h"
#include "nesterovBase.h"
#include "abacusLegalizer.h"
#include "abaxLegalizer.h"
#include "macroLegalizer.h"
#include <iostream>
#include <memory>
#include "log.h"
#include "fft.h"
#include "plot.h"

namespace replace
{

  using namespace std;

  Replace::Replace(float targetDensity)
      : pb_(nullptr), 
        initialPlaceMaxIter_(20),
        initialPlaceMinDiffLength_(1500),
        initialPlaceMaxSolverIter_(100),
        initialPlaceMaxFanout_(200),
        initialPlaceNetWeightScale_(800),
        nesterovPlaceMaxIter_(2000),
        binGridCntX_(0), binGridCntY_(0),
        overflow_(0.1f), density_(targetDensity),
        initDensityPenalityFactor_(0.00008f),
        initWireLengthCoef_(0.25f),
        minPhiCoef_(0.95f), maxPhiCoef_(1.05f),
        referenceHpwl_(446000000),
        incrementalPlaceMode_(false),
        useLocalDensity_(false),
        initLocalAlpha_(1e-12f), initLocalBeta_(1e-11f),
        useTheta_(false)
  {
  }

  Replace::~Replace()
  {
    reset();
  }

  void Replace::init()
  {
  }

  void Replace::reset()
  {
    initialPlaceMaxIter_ = 20;
    initialPlaceMinDiffLength_ = 1500;
    initialPlaceMaxSolverIter_ = 100;
    initialPlaceMaxFanout_ = 200;
    initialPlaceNetWeightScale_ = 800;

    nesterovPlaceMaxIter_ = 2000;
    binGridCntX_ = binGridCntY_ = 0;
    overflow_ = 0;
    density_ = 0;
    initDensityPenalityFactor_ = 0.0001f;
    initWireLengthCoef_ = 1.0f;
    minPhiCoef_ = 0.95f;
    maxPhiCoef_ = 1.05f;
    referenceHpwl_ = 446000000;

    incrementalPlaceMode_ = false;
    // verbose_ = 0;
  }

  void Replace::setPlacerBase(const std::shared_ptr<PlacerBase>& pb)
  {
    pb_ = pb;
  }

  void Replace::doInitialPlace()
  {
    InitialPlaceVars ipVars;
    ipVars.minIter = initialPlaceMinIter_;
    ipVars.maxIter = initialPlaceMaxIter_;
    ipVars.minDiffLength = initialPlaceMinDiffLength_;
    ipVars.maxSolverIter = initialPlaceMaxSolverIter_;
    ipVars.maxFanout = initialPlaceMaxFanout_;
    ipVars.netWeightScale = initialPlaceNetWeightScale_;
    ipVars.incrementalPlaceMode = incrementalPlaceMode_;

    std::unique_ptr<InitialPlace> ip(new InitialPlace(ipVars, pb_));
    ip->doBicgstabPlace();
  }

  void Replace::doNesterovPlace(string placename)
  {
    LOG_TRACE("start Replace::doNesterovPlace");
    NesterovBaseVars nbVars;
    nbVars.targetDensity = density_;

    if (binGridCntX_ != 0)
    {
      nbVars.isSetBinCntX = 1;
      nbVars.binCntX = binGridCntX_;
    }

    if (binGridCntY_ != 0)
    {
      nbVars.isSetBinCntY = 1;
      nbVars.binCntY = binGridCntY_;
    }

    std::shared_ptr<NesterovBase> nb = std::make_shared<NesterovBase>(nbVars, pb_);

    NesterovPlaceVars npVars;

    npVars.minPhiCoef = minPhiCoef_;
    npVars.maxPhiCoef = maxPhiCoef_;
    npVars.referenceHpwl = referenceHpwl_;
    npVars.initDensityPenalty = initDensityPenalityFactor_;
    npVars.initWireLengthCoef = initWireLengthCoef_;
    npVars.targetOverflow = overflow_;
    npVars.maxNesterovIter = nesterovPlaceMaxIter_;
    npVars.useLocalDensity = useLocalDensity_;
    npVars.initAlpha = initLocalAlpha_;
    npVars.initBeta = initLocalBeta_;
    npVars.useTheta = useTheta_;

    NesterovPlace np(npVars, nb);
    np.init();
    np.doNesterovPlace(placename);
  }

  void Replace::doAbacusLegalization()
  {
    LOG_TRACE("start Replace::doAbacusLegalization");
    AbacusLegalizerVars algVars;
    algVars.weightOpt = AbacusLegalizerVars::One;

    // make_unique is C++14 std
    // alg_ = std::make_unique<AbacusLegalizer>(algVars, pb_);
    AbacusLegalizer alg(algVars, pb_);
    alg.doLegalization();
    LOG_TRACE("end Replace::doAbacusLegalization");
  }

  void Replace::doAbaxLegalization()
  {
    LOG_TRACE("start Replace::doAbaxLegalization");
    AbaxLegalizer abax(pb_);
    abax.doLegalization();
    LOG_TRACE("end Replace::doAbaxLegalization");
  }

  void Replace::doMacroLegalization()
  {
    LOG_TRACE("start Replace::doMacroLegalization");

    MacroLegalizerVars mlgVars;
    mlgVars.maxPostLegalizeIter = 10000;

    MacroLegalizer mlg(mlgVars, pb_);
    mlg.doLegalization();

    LOG_TRACE("end Replace::doMacorLegalization");
  }

  void Replace::modifyTerminal()
  {
    // terminal resizing with spacing
    int spaceLeft = pb_->terminalSpacing();
    int spaceRight = pb_->terminalSpacing() - spaceLeft;
    for(Instance* inst : pb_->die("terminal")->insts())
    {
      inst->setBox(
        inst->lx() - spaceLeft,
        inst->ly() - spaceLeft,
        inst->ux() + spaceRight,
        inst->uy() + spaceRight
      );
    }

    // terminal die core area
    int coreLx = pb_->die("terminal")->coreLx() + spaceRight;
    int coreLy = pb_->die("terminal")->coreLy() + spaceRight;
    int coreUx = pb_->die("terminal")->coreUx() - spaceLeft;
    int coreUy = pb_->die("terminal")->coreUy() - spaceLeft;
    pb_->die("terminal")->setCoreBox(coreLx, coreLy, coreUx, coreUy);

    // set row Params
    int rowWidth = pb_->die("terminal")->coreDx();
    int rowHeight = pb_->terminalSizeX() + pb_->terminalSpacing();
    int repeatCount = pb_->die("terminal")->coreDy() / rowHeight;
    int rowStartX = pb_->die("terminal")->coreLx();
    int rowStartY = pb_->die("terminal")->coreLy();
    pb_->die("terminal")->setRowParams(rowStartX, rowStartY, rowWidth, rowHeight, repeatCount);
  }

  void Replace::recoverTerminal()
  {
    int spaceLeft = pb_->terminalSpacing();
    int spaceRight = pb_->terminalSpacing() - spaceLeft;
    for(Instance* inst : pb_->die("terminal")->insts())
    {
      inst->setBox(
        inst->lx() + spaceLeft,
        inst->ly() + spaceLeft,
        inst->ux() - spaceRight,
        inst->uy() - spaceRight
      );
    }

    // terminal die core area
    int coreLx = pb_->die("terminal")->coreLx() - spaceRight;
    int coreLy = pb_->die("terminal")->coreLy() - spaceRight;
    int coreUx = pb_->die("terminal")->coreUx() + spaceLeft;
    int coreUy = pb_->die("terminal")->coreUy() + spaceLeft;
    pb_->die("terminal")->setCoreBox(coreLx, coreLy, coreUx, coreUy);
  }

  void Replace::setInitialPlaceMinIter(int iter)
  {
    initialPlaceMinIter_ = iter;
  }

  void Replace::setInitialPlaceMaxIter(int iter)
  {
    initialPlaceMaxIter_ = iter;
  }

  void Replace::setInitialPlaceMinDiffLength(int length)
  {
    initialPlaceMinDiffLength_ = length;
  }

  void Replace::setInitialPlaceMaxSolverIter(int iter)
  {
    initialPlaceMaxSolverIter_ = iter;
  }

  void Replace::setInitialPlaceMaxFanout(int fanout)
  {
    initialPlaceMaxFanout_ = fanout;
  }

  void Replace::setInitialPlaceNetWeightScale(float scale)
  {
    initialPlaceNetWeightScale_ = scale;
  }

  void Replace::setNesterovPlaceMaxIter(int iter)
  {
    nesterovPlaceMaxIter_ = iter;
  }

  void Replace::setNesterovPlaceUseLocalDensity(bool on)
  {
    useLocalDensity_ = on;
  }

  void Replace::setNesterovInitLocalAlpha(float alpha)
  {
    initLocalAlpha_ = alpha;
  }

  void Replace::setNesterovInitLocalBeta(float beta)
  {
    initLocalBeta_ = beta;
  }

  void Replace::setNesterovUseTheta(bool on)
  {
    useTheta_ = on;
  }

  void Replace::setBinGridCntX(int binGridCntX)
  {
    binGridCntX_ = binGridCntX;
  }

  void Replace::setBinGridCntY(int binGridCntY)
  {
    binGridCntY_ = binGridCntY;
  }

  void Replace::setTargetOverflow(float overflow)
  {
    overflow_ = overflow;
  }

  void Replace::setTargetDensity(float density)
  {
    density_ = density;
  }

  void Replace::setInitDensityPenalityFactor(float penaltyFactor)
  {
    initDensityPenalityFactor_ = penaltyFactor;
  }

  void Replace::setInitWireLengthCoef(float coef)
  {
    initWireLengthCoef_ = coef;
  }

  void Replace::setMinPhiCoef(float minPhiCoef)
  {
    minPhiCoef_ = minPhiCoef;
  }

  void Replace::setMaxPhiCoef(float maxPhiCoef)
  {
    maxPhiCoef_ = maxPhiCoef;
  }

  void Replace::setReferenceHpwl(float refHpwl)
  {
    referenceHpwl_ = refHpwl;
  }

  void Replace::setIncrementalPlaceMode(bool mode)
  {
    incrementalPlaceMode_ = mode;
  }
}
