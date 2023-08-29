#ifndef __REPLACE_FFT__
#define __REPLACE_FFT__

#include <vector>
#include "point.h"

namespace replace {

class FFT {
  public:
    FFT();
    ~FFT() = default;

    void init(int binCntX, int binCntY, double binSizeX, double binSizeY);

    // input func
    void updateDensity(int x, int y, double density);

    // do FFT
    void doFFT();

    // returning func
    std::pair<double, double> getElectroForce(int x, int y);
    double getElectroPhi(int x, int y);

  private:
    // 2D array; width: binCntX_, height: binCntY_;
    // No hope to use Vector at this moment...
    //double** binDensity_;
    //double** electroPhi_;
    //double** electroForceX_;
    //double** electroForceY_;

    std::vector<double> binDensityStor_;
    std::vector<double*> binDensity_;

    std::vector<double> electroPhiStor_;
    std::vector<double*> electroPhi_;

    std::vector<double> electroForceXStor_;
    std::vector<double*> electroForceX_;

    std::vector<double> electroForceYStor_;
    std::vector<double*> electroForceY_;

    // cos/sin table (prev: w_2d)
    // length:  max(binCntX, binCntY) * 3 / 2
    std::vector<double> csTable_;

    // wx. length:  binCntX_
    std::vector<double> wx_;
    std::vector<double> wxSquare_;

    // wy. length:  binCntY_
    std::vector<double> wy_;
    std::vector<double> wySquare_;

    // work area for bit reversal (prev: ip)
    // length: round(sqrt( max(binCntX_, binCntY_) )) + 2
    std::vector<int> workArea_;

    int binCntX_;
    int binCntY_;
    double binSizeX_;
    double binSizeY_;
};

}

#endif
