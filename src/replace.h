#ifndef __REPLACE_HEADER__
#define __REPLACE_HEADER__

#include <memory>
#include <string>
using std::string;

namespace replace
{
  class PlacerBase;
  class NesterovBase;

  class InitialPlace;
  class NesterovPlace;
  class AbacusLegalizer;
  class LayerAssignmenter;
  class MacroLegalizer;

  class Replace
  {
  public:
    // targetDesntiy can be set between 0.7~1, or 1~2 if we place cells of 
    // two layers on one layer
    Replace(float targetDensity);
    ~Replace();

    void init();
    void reset();

    void setPlacerBase(const std::shared_ptr<PlacerBase>& pb);

    void doInitialPlace();
    // placename is only used for LOG and plotted picture names to distinguish
    // the placements before and after partition
    void doNesterovPlace(string placename = "");
    void doAbacusLegalization();
    void doAbaxLegalization();
    void doMacroLegalization();

    void modifyTerminal();
    void recoverTerminal();

    // Initial Place param settings
    void setInitialPlaceMinIter(int iter);
    void setInitialPlaceMaxIter(int iter);
    void setInitialPlaceMinDiffLength(int length);
    void setInitialPlaceMaxSolverIter(int iter);
    void setInitialPlaceMaxFanout(int fanout);
    void setInitialPlaceNetWeightScale(float scale);

    void setNesterovPlaceMaxIter(int iter);
    void setNesterovPlaceUseLocalDensity(bool on);
    void setNesterovInitLocalAlpha(float alpha);
    void setNesterovInitLocalBeta(float beta);
    void setNesterovUseTheta(bool on);

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

    void setMacroPostMaxIter(int iter);

  private:
    std::shared_ptr<PlacerBase> pb_;

    int initialPlaceMinIter_;
    int initialPlaceMaxIter_;
    int initialPlaceMinDiffLength_;
    int initialPlaceMaxSolverIter_;
    int initialPlaceMaxFanout_;
    float initialPlaceNetWeightScale_;
    bool incrementalPlaceMode_;

    int nesterovPlaceMaxIter_;
    int binGridCntX_;
    int binGridCntY_;
    float overflow_;      // target density overflow
    float density_;       // target density
    float initDensityPenalityFactor_;
    float initWireLengthCoef_;
    float minPhiCoef_;
    float maxPhiCoef_;
    float referenceHpwl_;
    bool useLocalDensity_;
    float initLocalAlpha_;
    float initLocalBeta_;
    bool useTheta_;

    int macroPostMaxIter_;
  };
}

#endif
