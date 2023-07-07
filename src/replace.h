#ifndef __REPLACE_HEADER__
#define __REPLACE_HEADER__

#include <memory>

namespace replace
{
  class PlacerBase;
  class NesterovBase;

  class InitialPlace;
  class NesterovPlace;
  class AbacusLegalizer;

  class Replace
  {
  public:
    Replace();
    ~Replace();

    void init();
    void reset();

    void setPlacerBase(const std::shared_ptr<PlacerBase>& pb);

    void doInitialPlace();
    void doNesterovPlace();
    void applyGlobalPlacement();
    void doAbacusLegalization();

    // Initial Place param settings
    void setInitialPlaceMaxIter(int iter);
    void setInitialPlaceMinDiffLength(int length);
    void setInitialPlaceMaxSolverIter(int iter);
    void setInitialPlaceMaxFanout(int fanout);
    void setInitialPlaceNetWeightScale(float scale);

    void setNesterovPlaceMaxIter(int iter);

    void setBinGridCntX(int binGridCntX);
    void setBinGridCntY(int binGridCntY);

    void setTargetDensity(float density);
    void setTargetOverflow(float overflow);
    void setInitDensityPenalityFactor(float penaltyFactor);
    void setInitWireLengthCoef(float coef);
    void setMinPhiCoef(float minPhiCoef);
    void setMaxPhiCoef(float maxPhiCoef);

    // HPWL: half-parameter wire length.
    void setReferenceHpwl(float deltaHpwl);

    void setIncrementalPlaceMode(bool mode);

  private:
    std::shared_ptr<PlacerBase> pb_;
    std::shared_ptr<NesterovBase> nb_;

    std::unique_ptr<InitialPlace> ip_;
    std::unique_ptr<NesterovPlace> np_;
    std::unique_ptr<AbacusLegalizer> alg_;

    int initialPlaceMaxIter_;
    int initialPlaceMinDiffLength_;
    int initialPlaceMaxSolverIter_;
    int initialPlaceMaxFanout_;
    float initialPlaceNetWeightScale_;

    int nesterovPlaceMaxIter_;
    int binGridCntX_;
    int binGridCntY_;
    float overflow_;
    float density_;
    float initDensityPenalityFactor_;
    float initWireLengthCoef_;
    float minPhiCoef_;
    float maxPhiCoef_;
    float referenceHpwl_;

    bool incrementalPlaceMode_;
  };
}

#endif
