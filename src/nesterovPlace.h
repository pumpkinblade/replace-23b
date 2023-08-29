#ifndef __REPLACE_NESTEROV_PLACE__
#define __REPLACE_NESTEROV_PLACE__

#include "point.h"
#include <memory>
#include <vector>
#include <string>
using std::string;

namespace replace
{
  class PlacerBase;
  class Instance;
  class NesterovBase;

  class NesterovPlaceVars
  {
  public:
    int maxNesterovIter;
    int maxBackTrack;
    double initDensityPenalty;          // INIT_LAMBDA
    double initWireLengthCoef;          // base_wcof
    double targetOverflow;              // overflow
    double minPhiCoef;                  // pcof_min
    double maxPhiCoef;                  // pcof_max
    double minPreconditioner;           // MIN_PRE
    double initialPrevCoordiUpdateCoef; // z_ref_alpha
    double referenceHpwl;               // refDeltaHpwl
    // local density
    bool useLocalDensity;
    double initAlpha;
    double initBeta;
    // theta
    bool useTheta;
    NesterovPlaceVars();
  };

  class NesterovPlace
  {
  public:
    NesterovPlace();
    NesterovPlace(NesterovPlaceVars npVars, std::shared_ptr<NesterovBase> nb);
    ~NesterovPlace() = default;

    void init();
    void doNesterovPlace(string placename = "");

    void updateGradients(std::vector<Point> &sumGrads,
                         std::vector<double> &sumGradTheta,
                         std::vector<double> &refCellDeltas,
                         std::vector<double> &outCellDeltas);

    void updateWireLengthCoef(double overflow);

    void updateInitialPrevSLPCoordi();

    double getStepLength(
        std::vector<Point> &prevCoordi_,
        std::vector<Point> &prevSumGrads_,
        std::vector<Point> &curCoordi_,
        std::vector<Point> &curSumGrads_);

    void updateNextIter();
    double getPhiCoef(double scaledDiffHpwl);

    void updatePlacerBase();

    void determinMacroOrient();

  private:
    std::shared_ptr<NesterovBase> nb_;
    NesterovPlaceVars npVars_;

    // SLP is Step Length Prediction.
    //
    // y_st, y_dst, y_wdst, w_pdst
    std::vector<Point> curSLPCoordi_;
    std::vector<Point> curSLPSumGrads_;

    // y0_st, y0_dst, y0_wdst, y0_pdst
    std::vector<Point> nextSLPCoordi_;
    std::vector<Point> nextSLPSumGrads_;

    // z_st, z_dst, z_wdst, z_pdst
    std::vector<Point> prevSLPCoordi_;
    std::vector<Point> prevSLPSumGrads_;

    // x_st and x0_st
    std::vector<Point> curCoordi_;
    std::vector<Point> nextCoordi_;

    double wireLengthGradSum_;
    double densityGradSum_;
    double localDensityGradSum_;

    // alpha
    double stepLength_;

    // opt_phi_cof
    double densityPenalty_;

    // base_wcof
    double baseWireLengthCoef_;

    // wlen_cof
    double wireLengthCoefX_;
    double wireLengthCoefY_;

    // phi is described in ePlace paper.
    double sumPhi_;
    double sumOverflow_;

    // half-parameter-wire-length
    double prevHpwl_;

    bool isDiverged_;

    // local density
    double localAlpha_;
    double localBeta_;
    std::vector<double> prevCellDelta_;
    std::vector<double> curCellDelta_;
    std::vector<double> nextCellDelta_;

    // theta
    int useThetaMacroCount;
    std::vector<double> curTheta_;
    std::vector<double> nextTheta_;
    std::vector<double> curSLPTheta_;
    std::vector<double> curSLPSumGradTheta_;
    std::vector<double> nextSLPTheta_;
    std::vector<double> nextSLPSumGradTheta_;
    std::vector<double> prevSLPTheta_;
    std::vector<double> prevSLPSumGradTheta_;
    std::vector<double> wireLengthPrecondiTheta_;
    std::vector<double> densityPrecondiTheta_;
    double wireLengthGradSumTheta_;
    double densityGradSumTheta_;
    double localDensityGradSumTheta_;
  };
}

#endif
