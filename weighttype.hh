// weighttype.hh

#ifndef WEIGHTTYPE_HEADER
#define WEIGHTTYPE_HEADER

#include <string>

namespace scala{

  class WeightType {
  public:
    // Weight type for averaging
    //   UNIT        unit weights
    //   VARIANCE    weight = 1/variance
    //   SQRTSCALE   weight = sqrt(g)  g = 1/scale
    //   SCALE       weight = g        g = 1/scale
    enum AverageWeightType {UNIT, VARIANCE, SQRTSCALE, SCALE};

    static std::string formatWeightType (const AverageWeightType& weighttype);

  };
}

#endif
