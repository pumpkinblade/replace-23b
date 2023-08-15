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
  NesterovPlaceVars();
};

class NesterovPlace {
public:
  NesterovPlace();
  NesterovPlace(NesterovPlaceVars npVars, std::shared_ptr<NesterovBase> nb);
  ~NesterovPlace() = default;

  void init();
  void doNesterovPlace(string placename = "");

  void updateGradients(
      std::vector<Point>& sumGrads,
      std::vector<Point>& wireLengthGrads,
      std::vector<Point>& densityGrads );

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

private:
  std::shared_ptr<NesterovBase> nb_;
  NesterovPlaceVars npVars_;

  // SLP is Step Length Prediction.
  //
  // y_st, y_dst, y_wdst, w_pdst
  std::vector<Point> curSLPCoordi_;
  std::vector<Point> curSLPWireLengthGrads_;
  std::vector<Point> curSLPDensityGrads_;
  std::vector<Point> curSLPSumGrads_;

  // y0_st, y0_dst, y0_wdst, y0_pdst
  std::vector<Point> nextSLPCoordi_;
  std::vector<Point> nextSLPWireLengthGrads_;
  std::vector<Point> nextSLPDensityGrads_;
  std::vector<Point> nextSLPSumGrads_;

  // z_st, z_dst, z_wdst, z_pdst
  std::vector<Point> prevSLPCoordi_;
  std::vector<Point> prevSLPWireLengthGrads_;
  std::vector<Point> prevSLPDensityGrads_;
  std::vector<Point> prevSLPSumGrads_;

  // x_st and x0_st
  std::vector<Point> curCoordi_;
  std::vector<Point> nextCoordi_;

  prec wireLengthGradSum_;
  prec densityGradSum_;

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
  prec sumPhi_;
  prec sumOverflow_;

  // half-parameter-wire-length
  int64_t prevHpwl_;

  prec isDiverged_;
};
}

#endif
