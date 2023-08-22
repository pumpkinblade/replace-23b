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

class NesterovPlaceVars {
  public:
  int maxNesterovIter;
  int maxBackTrack;
  prec initDensityPenalty; // INIT_LAMBDA
  prec initWireLengthCoef; // base_wcof
  prec targetOverflow; // overflow
  prec minPhiCoef; // pcof_min
  prec maxPhiCoef; // pcof_max
  prec minPreconditioner; // MIN_PRE
  prec initialPrevCoordiUpdateCoef; // z_ref_alpha
  prec referenceHpwl; // refDeltaHpwl
  // local density
  bool useLocalDensity;
  prec initAlpha;
  prec initBeta;
  // theta
  bool useTheta;
  NesterovPlaceVars();
};

class NesterovPlace {
public:
  NesterovPlace();
  NesterovPlace(NesterovPlaceVars npVars, std::shared_ptr<NesterovBase> nb);
  ~NesterovPlace() = default;

  void init();
  void doNesterovPlace(string placename = "");

  void updateGradients(std::vector<Point>& sumGrads, std::vector<prec>& sumGradTheta);

  void updateWireLengthCoef(prec overflow);

  void updateInitialPrevSLPCoordi();

  prec getStepLength(
      std::vector<Point>& prevCoordi_,
      std::vector<Point>& prevSumGrads_,
      std::vector<Point>& curCoordi_,
      std::vector<Point>& curSumGrads_ );

  void updateNextIter();
  prec getPhiCoef(prec scaledDiffHpwl);

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

  prec wireLengthGradSum_;
  prec densityGradSum_;
  prec localDensityGradSum_;

  // alpha
  prec stepLength_;

  // opt_phi_cof
  prec densityPenalty_;

  // base_wcof
  prec baseWireLengthCoef_;

  // wlen_cof
  prec wireLengthCoefX_;
  prec wireLengthCoefY_;

  // phi is described in ePlace paper.
  double sumPhi_;
  prec sumOverflow_;

  // half-parameter-wire-length
  double prevHpwl_;

  bool isDiverged_;

  // local density
  prec localAlpha_;
  prec localBeta_;
  std::vector<prec> curCellDelta_;
  std::vector<prec> nextCellDelta_;

  // theta
  std::vector<int> macroIndices_;
  std::vector<prec> curTheta_;
  std::vector<prec> nextTheta_;
  std::vector<prec> curSLPTheta_;
  std::vector<prec> curSLPSumGradTheta_;
  std::vector<prec> nextSLPTheta_;
  std::vector<prec> nextSLPSumGradTheta_;
  std::vector<prec> prevSLPTheta_;
  std::vector<prec> prevSLPSumGradTheta_;
  std::vector<prec> wireLengthPrecondiTheta_;
  std::vector<prec> densityPrecondiTheta_;
  prec wireLengthGradSumTheta_;
  prec densityGradSumTheta_;
  prec localDensityGradSumTheta_;
};
}

#endif
