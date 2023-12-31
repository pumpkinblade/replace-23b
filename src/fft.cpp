#include <cstdlib>
#include <cmath>
#include <cfloat>

#include <iostream>

#include "fft.h"
#include "log.h"

#define REPLACE_FFT_PI 3.141592653589793238462L 

//
// The following FFT library came from
// http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html
//
//
/// 1D FFT ////////////////////////////////////////////////////////////////
void cdft(int n, int isgn, double *a, int *ip, double *w);
void ddct(int n, int isgn, double *a, int *ip, double *w);
void ddst(int n, int isgn, double *a, int *ip, double *w);

/// 2D FFT ////////////////////////////////////////////////////////////////
void cdft2d(int, int, int, double **, double *, int *, double *);
void rdft2d(int, int, int, double **, double *, int *, double *);
void ddct2d(int, int, int, double **, double *, int *, double *);
void ddst2d(int, int, int, double **, double *, int *, double *);
void ddsct2d(int n1, int n2, int isgn, double **a, double *t, int *ip, double *w);
void ddcst2d(int n1, int n2, int isgn, double **a, double *t, int *ip, double *w);

/// 3D FFT ////////////////////////////////////////////////////////////////
void cdft3d(int, int, int, int, double ***, double *, int *, double *);
void rdft3d(int, int, int, int, double ***, double *, int *, double *);
void ddct3d(int, int, int, int, double ***, double *, int *, double *);
void ddst3d(int, int, int, int, double ***, double *, int *, double *);
void ddscct3d(int, int, int, int isgn, double ***, double *, int *, double *);
void ddcsct3d(int, int, int, int isgn, double ***, double *, int *, double *);
void ddccst3d(int, int, int, int isgn, double ***, double *, int *, double *);

namespace replace {

FFT::FFT()
    : binCntX_(0), binCntY_(0), binSizeX_(0), binSizeY_(0)
{
}

void 
FFT::init(int binCntX, int binCntY, double binSizeX, double binSizeY) {
  binCntX_ = binCntX;
  binCntY_ = binCntY;
  binSizeX_ = binSizeX;
  binSizeY_ = binSizeY;

  binDensityStor_.resize(binCntX_ * binCntY_, 0);
  electroPhiStor_.resize(binCntX_ * binCntY_, 0);
  electroForceXStor_.resize(binCntX_ * binCntY_, 0);
  electroForceYStor_.resize(binCntX_ * binCntY_, 0);

  binDensity_.resize(binCntX_, nullptr);
  electroPhi_.resize(binCntX_, nullptr);
  electroForceX_.resize(binCntX_, nullptr);
  electroForceY_.resize(binCntX_, nullptr);
  for (int i = 0; i < binCntX_; i++)
  {
    binDensity_[i] = &binDensityStor_[i * binCntY_];
    electroPhi_[i] = &electroPhiStor_[i * binCntY_];
    electroForceX_[i] = &electroForceXStor_[i * binCntY_];
    electroForceY_[i] = &electroForceYStor_[i * binCntY_];
  }

  csTable_.resize( std::max(binCntX_, binCntY_) * 3 / 2, 0 );

  wx_.resize( binCntX_, 0 );
  wxSquare_.resize( binCntX_, 0);
  wy_.resize( binCntY_, 0 );
  wySquare_.resize( binCntY_, 0 );

  workArea_.resize( round(sqrt(std::max(binCntX_, binCntY_))) + 2, 0 );
 
  for(int i=0; i<binCntX_; i++) {
    wx_[i] = REPLACE_FFT_PI * static_cast<double>(i) 
      / static_cast<double>(binCntX_);
    wxSquare_[i] = wx_[i] * wx_[i]; 
  }

  for(int i=0; i<binCntY_; i++) {
    wy_[i] = REPLACE_FFT_PI * static_cast<double>(i)
      / static_cast<double>(binCntY_) 
      * static_cast<double>(binSizeY_) 
      / static_cast<double>(binSizeX_);
    wySquare_[i] = wy_[i] * wy_[i];
  }
}

void
FFT::updateDensity(int x, int y, double density) {
  binDensity_[x][y] = density;
}

std::pair<double, double> 
FFT::getElectroForce(int x, int y) {
  return std::make_pair(
      electroForceX_[x][y],
      electroForceY_[x][y]);
}

double
FFT::getElectroPhi(int x, int y) {
  return electroPhi_[x][y]; 
}

using namespace std;

void
FFT::doFFT() {
  ddct2d(binCntX_, binCntY_, -1, binDensity_.data(), 
      NULL, (int*) &workArea_[0], (double*)&csTable_[0]);
  
  for(int i = 0; i < binCntX_; i++) {
    binDensity_[i][0] *= 0.5;
  }

  for(int i = 0; i < binCntY_; i++) {
    binDensity_[0][i] *= 0.5;
  }

  for(int i = 0; i < binCntX_; i++) {
    for(int j = 0; j < binCntY_; j++) {
      binDensity_[i][j] *= 4.0 / binCntX_ / binCntY_;
    }
  }

  for(int i = 0; i < binCntX_; i++) {
    double wx = wx_[i];
    double wx2 = wxSquare_[i];

    for(int j = 0; j < binCntY_; j++) {
      double wy = wy_[j];
      double wy2 = wySquare_[j];

      double density = binDensity_[i][j];
      double phi = 0;
      double electroX = 0, electroY = 0;

      if(i == 0 && j == 0) {
        phi = electroX = electroY = 0.0f;
      }
      else {
        //////////// lutong
        //  denom =
        //  wx2 / 4.0 +
        //  wy2 / 4.0 ;
        // a_phi = a_den / denom ;
        ////b_phi = 0 ; // -1.0 * b / denom ;
        ////a_ex = 0 ; // b_phi * wx ;
        // a_ex = a_phi * wx / 2.0 ;
        ////a_ey = 0 ; // b_phi * wy ;
        // a_ey = a_phi * wy / 2.0 ;
        ///////////
        phi = density / (wx2 + wy2);
        electroX = phi * wx;
        electroY = phi * wy;
      }
      electroPhi_[i][j] = phi;
      electroForceX_[i][j] = electroX;
      electroForceY_[i][j] = electroY;
    }
  }
  // Inverse DCT
  ddct2d(binCntX_, binCntY_, 1, 
      electroPhi_.data(), NULL, 
      (int*) &workArea_[0], (double*) &csTable_[0]);
  ddsct2d(binCntX_, binCntY_, 1, 
      electroForceX_.data(), NULL, 
      (int*) &workArea_[0], (double*) &csTable_[0]);
  ddcst2d(binCntX_, binCntY_, 1, 
      electroForceY_.data(), NULL, 
      (int*) &workArea_[0], (double*) &csTable_[0]);
}

}
