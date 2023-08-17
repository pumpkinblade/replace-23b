#ifndef __REPLACE_FFT__
#define __REPLACE_FFT__

#include <vector>
#include "point.h"

namespace replace {

class FFT {
  public:
    FFT();
    ~FFT() = default;

    void init(int binCntX, int binCntY, prec binSizeX, prec binSizeY);

    // input func
    void updateDensity(int x, int y, prec density);

    // do FFT
    void doFFT();

    // returning func
    std::pair<prec, prec> getElectroForce(int x, int y);
    prec getElectroPhi(int x, int y);

  private:
    // 2D array; width: binCntX_, height: binCntY_;
    // No hope to use Vector at this moment...
    //prec** binDensity_;
    //prec** electroPhi_;
    //prec** electroForceX_;
    //prec** electroForceY_;

    std::vector<prec> binDensityStor_;
    std::vector<prec*> binDensity_;

    std::vector<prec> electroPhiStor_;
    std::vector<prec*> electroPhi_;

    std::vector<prec> electroForceXStor_;
    std::vector<prec*> electroForceX_;

    std::vector<prec> electroForceYStor_;
    std::vector<prec*> electroForceY_;

    // cos/sin table (prev: w_2d)
    // length:  max(binCntX, binCntY) * 3 / 2
    std::vector<prec> csTable_;

    // wx. length:  binCntX_
    std::vector<prec> wx_;
    std::vector<prec> wxSquare_;

    // wy. length:  binCntY_
    std::vector<prec> wy_;
    std::vector<prec> wySquare_;

    // work area for bit reversal (prev: ip)
    // length: round(sqrt( max(binCntX_, binCntY_) )) + 2
    std::vector<int> workArea_;

    int binCntX_;
    int binCntY_;
    prec binSizeX_;
    prec binSizeY_;
};

}

#endif
