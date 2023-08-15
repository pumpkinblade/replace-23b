#include <cstdlib>
#include <cmath>
#include <cfloat>

#include <iostream>

#include "fft.h"
#include "log.h"

#define REPLACE_FFT_PI 3.141592653589793238462L 

namespace replace {


FFT::FFT()
    : binCntX_(0), binCntY_(0), binSizeX_(0), binSizeY_(0) {}

FFT::FFT(int binCntX, int binCntY, int binSizeX, int binSizeY)
    : binCntX_(binCntX), binCntY_(binCntY), binSizeX_(binSizeX), binSizeY_(binSizeY) 
{
  init();   
}

void
FFT::init() {
  //binDensity_ = new prec*[binCntX_];
  //electroPhi_ = new prec*[binCntX_];
  //electroForceX_ = new prec*[binCntX_];
  //electroForceY_ = new prec*[binCntX_];

  //for(int i=0; i<binCntX_; i++) {
  //  binDensity_[i] = new prec[binCntY_];
  //  electroPhi_[i] = new prec[binCntY_];
  //  electroForceX_[i] = new prec[binCntY_];
  //  electroForceY_[i] = new prec[binCntY_];

  //  for(int j=0; j<binCntY_; j++) {
  //    binDensity_[i][j] 
  //      = electroPhi_[i][j] 
  //      = electroForceX_[i][j] 
  //      = electroForceY_[i][j] 
  //      = 0.0f;  
  //  }
  //}

  binDensityStor_.resize(binCntX_ * binCntY_, (prec)0);
  electroPhiStor_.resize(binCntX_ * binCntY_, (prec)0);
  electroForceXStor_.resize(binCntX_ * binCntY_, (prec)0);
  electroForceYStor_.resize(binCntX_ * binCntY_, (prec)0);

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
    wx_[i] = REPLACE_FFT_PI * static_cast<prec>(i) 
      / static_cast<prec>(binCntX_);
    wxSquare_[i] = wx_[i] * wx_[i]; 
  }

  for(int i=0; i<binCntY_; i++) {
    wy_[i] = REPLACE_FFT_PI * static_cast<prec>(i)
      / static_cast<prec>(binCntY_) 
      * static_cast<prec>(binSizeY_) 
      / static_cast<prec>(binSizeX_);
    wySquare_[i] = wy_[i] * wy_[i];
  }
}

void
FFT::updateDensity(int x, int y, prec density) {
  binDensity_[x][y] = density;
}

std::pair<prec, prec> 
FFT::getElectroForce(int x, int y) {
  return std::make_pair(
      electroForceX_[x][y],
      electroForceY_[x][y]);
}

prec
FFT::getElectroPhi(int x, int y) {
  return electroPhi_[x][y]; 
}

using namespace std;

void
FFT::doFFT() {
  ddct2d(binCntX_, binCntY_, -1, binDensity_.data(), 
      NULL, (int*) &workArea_[0], (prec*)&csTable_[0]);
  
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
    prec wx = wx_[i];
    prec wx2 = wxSquare_[i];

    for(int j = 0; j < binCntY_; j++) {
      prec wy = wy_[j];
      prec wy2 = wySquare_[j];

      prec density = binDensity_[i][j];
      prec phi = 0;
      prec electroX = 0, electroY = 0;

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
      (int*) &workArea_[0], (prec*) &csTable_[0]);
  ddsct2d(binCntX_, binCntY_, 1, 
      electroForceX_.data(), NULL, 
      (int*) &workArea_[0], (prec*) &csTable_[0]);
  ddcst2d(binCntX_, binCntY_, 1, 
      electroForceY_.data(), NULL, 
      (int*) &workArea_[0], (prec*) &csTable_[0]);
}

}
