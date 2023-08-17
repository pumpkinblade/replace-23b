#ifndef __REPLACE_COORDI__
#define __REPLACE_COORDI__

namespace replace {

  class FloatPoint
  {
  public:
    float x;
    float y;
    FloatPoint() : x(0.f), y(0.f) {}
    FloatPoint(float x, float y) : x(x), y(y) {}
  };

  class IntPoint
  {
  public:
    int x;
    int y;
    IntPoint() : x(0), y(0) {}
    IntPoint(int x, int y) : x(x), y(y) {}
  };

  class DoublePoint
  {
  public:
    double x;
    double y;
    DoublePoint() : x(0.), y(0.) {}
    DoublePoint(double x, double y) : x(x), y(y) {}
  };
}

using prec = double;
using Point = replace::DoublePoint;

#endif
