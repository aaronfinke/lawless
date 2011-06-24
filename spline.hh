// spline.hh

#ifndef SCALA_SPLINE
#define SCALA_SPLINE

#include <vector>
typedef std::pair<float,float> RPair;

namespace scala
{
  class Spline
  // Translated from Fortran code in truncate
  {
  public:
    Spline() : n(0) {}
    // construct from set of x,y pairs
    Spline(const std::vector<RPair>& xyin);

    float Interpolate(const float& xx) const;

  private:
    int n;
    std::vector<double> x; 
    std::vector<double> y; 
    std::vector<double> y2;  // second derivatives of interpolating
			      // function
  };
}

#endif
