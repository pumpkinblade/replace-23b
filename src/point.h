#ifndef __REPLACE_COORDI__
#define __REPLACE_COORDI__

class Point
{
public:
  double x;
  double y;
  Point() : x(0.), y(0.) {}
  Point(double x, double y) : x(x), y(y) {}
};

constexpr double PI = 3.1415926535897932;
constexpr double SQRT2 = 1.4142135623730951;
constexpr double LEGAL_THETA[5] = {0.0, 0.5 * PI, 1.0 * PI, 1.5 * PI, 2.0 * PI};
enum class Orientation : unsigned char { R0, R90, R180, R270, R360 };

#endif
