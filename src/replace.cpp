#include "replace.h"
#include "initialPlace.h"
#include "nesterovPlace.h"
#include "placerBase.h"
#include "nesterovBase.h"
#include "abacusLegalizer.h"
#include "macroLegalizer.h"
#include <iostream>
#include <memory>
#include "log.h"

namespace replace
{

  using namespace std;

  Replace::Replace(float targetDensity)
      : pb_(nullptr), nb_(nullptr),
        ip_(nullptr), np_(nullptr), alg_(nullptr),
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
        incrementalPlaceMode_(false)
        // verbose_(0)
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
    ipVars.maxIter = initialPlaceMaxIter_;
    ipVars.minDiffLength = initialPlaceMinDiffLength_;
    ipVars.maxSolverIter = initialPlaceMaxSolverIter_;
    ipVars.maxFanout = initialPlaceMaxFanout_;
    ipVars.netWeightScale = initialPlaceNetWeightScale_;
    ipVars.incrementalPlaceMode = incrementalPlaceMode_;

    std::unique_ptr<InitialPlace> ip(new InitialPlace(ipVars, pb_));
    ip_ = std::move(ip);
    ip_->doBicgstabPlace();
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

    nb_ = std::make_shared<NesterovBase>(nbVars, pb_);

    NesterovPlaceVars npVars;

    npVars.minPhiCoef = minPhiCoef_;
    npVars.maxPhiCoef = maxPhiCoef_;
    npVars.referenceHpwl = referenceHpwl_;
    npVars.initDensityPenalty = initDensityPenalityFactor_;
    npVars.initWireLengthCoef = initWireLengthCoef_;
    npVars.targetOverflow = overflow_;
    npVars.maxNesterovIter = nesterovPlaceMaxIter_;

    std::unique_ptr<NesterovPlace> np(new NesterovPlace(npVars, pb_, nb_));
    np_ = std::move(np);

    np_->doNesterovPlace(placename);
  }

  void Replace::doAbacusLegalization()
  {
    LOG_TRACE("start Replace::doAbacusLegalization");
    AbacusLegalizerVars algVars;
    algVars.weightOpt = AbacusLegalizerVars::One;

    // make_unique is C++14 std
    // alg_ = std::make_unique<AbacusLegalizer>(algVars, pb_);
    alg_.reset(new AbacusLegalizer(algVars, pb_));
    alg_->doLegalization();
    LOG_TRACE("end Replace::doAbacusLegalization");
  }

  void Replace::doMacroLegalization()
  {
    LOG_TRACE("start Replace::doMacroLegalization");

    MacroLegalizerVars mlgVars;
    mlgVars.maxPostLegalizeIter = 10000;

    mlg_.reset(new MacroLegalizer(mlgVars, pb_));
    mlg_->doLegalization();

    LOG_TRACE("end Replace::doMacorLegalization");
  }

  void Replace::doSAMacroLegalization()
  {
    LOG_TRACE("start Replace::doSAMacroLegalization");

    MacroLegalizerVars mlgVars;
    mlgVars.maxPostLegalizeIter = 1000;

    mlg_.reset(new MacroLegalizer(mlgVars, pb_));
    mlg_->doLegalization();

    LOG_TRACE("end Replace::doMacorLegalization");
  }

  void
  Replace::setInitialPlaceMaxIter(int iter)
  {
    initialPlaceMaxIter_ = iter;
  }

  void
  Replace::setInitialPlaceMinDiffLength(int length)
  {
    initialPlaceMinDiffLength_ = length;
  }

  void
  Replace::setInitialPlaceMaxSolverIter(int iter)
  {
    initialPlaceMaxSolverIter_ = iter;
  }

  void
  Replace::setInitialPlaceMaxFanout(int fanout)
  {
    initialPlaceMaxFanout_ = fanout;
  }

  void
  Replace::setInitialPlaceNetWeightScale(float scale)
  {
    initialPlaceNetWeightScale_ = scale;
  }

  void Replace::setNesterovPlaceMaxIter(int iter)
  {
    nesterovPlaceMaxIter_ = iter;
  }

  void
  Replace::setBinGridCntX(int binGridCntX)
  {
    binGridCntX_ = binGridCntX;
  }

  void
  Replace::setBinGridCntY(int binGridCntY)
  {
    binGridCntY_ = binGridCntY;
  }

  void
  Replace::setTargetOverflow(float overflow)
  {
    overflow_ = overflow;
  }

  void
  Replace::setTargetDensity(float density)
  {
    density_ = density;
  }

  void
  Replace::setInitDensityPenalityFactor(float penaltyFactor)
  {
    initDensityPenalityFactor_ = penaltyFactor;
  }

  void
  Replace::setInitWireLengthCoef(float coef)
  {
    initWireLengthCoef_ = coef;
  }

  void
  Replace::setMinPhiCoef(float minPhiCoef)
  {
    minPhiCoef_ = minPhiCoef;
  }

  void
  Replace::setMaxPhiCoef(float maxPhiCoef)
  {
    maxPhiCoef_ = maxPhiCoef;
  }

  void
  Replace::setReferenceHpwl(float refHpwl)
  {
    referenceHpwl_ = refHpwl;
  }

  void
  Replace::setIncrementalPlaceMode(bool mode)
  {
    incrementalPlaceMode_ = mode;
  }
}
