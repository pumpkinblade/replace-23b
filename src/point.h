#ifndef __REPLACE_COORDI__
#define __REPLACE_COORDI__

class FloatPoint
{
public:
  float x;
  float y;
  FloatPoint() : x(0.f), y(0.f) {}
  FloatPoint(float x, float y) : x(x), y(y) {}
};

class FloatPoint3
{
public:
  float x;
  float y;
  float z;
  FloatPoint3() : x(0.f), y(0.f), z(0.f) {}
  FloatPoint3(float x, float y, float z) : x(x), y(y), z(z) {}
};

class DoublePoint
{
public:
  double x;
  double y;
  DoublePoint() : x(0.), y(0.) {}
  DoublePoint(double x, double y) : x(x), y(y) {}
};

class DoublePoint3
{
public:
  double x;
  double y;
  double z;
  DoublePoint3() : x(0.), y(0.), z(0.) {}
  DoublePoint3(double x, double y, double z) : x(x), y(y), z(z) {}
};

using prec = double;
using Point = ::DoublePoint;
using Point3 = ::DoublePoint3;

constexpr double PI = 3.1415926535897932;
constexpr double LEGAL_THETA[5] = {0.0, 0.5 * PI, 1.0 * PI, 1.5 * PI, 2.0 * PI};
enum class Orientation : unsigned char { R0, R90, R180, R270, R360 };

#endif
