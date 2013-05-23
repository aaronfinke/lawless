// weighttype.hh

#ifndef WEIGHTTYPE_HEADER
#define WEIGHTTYPE_HEADER

namespace scala{

  class WeightType {
  public:
    // Weight type for averaging
    //   UNIT        unit weights
    //   VARIANCE    weight = 1/variance
    //   SQRTSCALE   weight = 1/sqrt(g)  g = 1/scale
    enum AverageWeightType {UNIT, VARIANCE, SQRTSCALE};
  };
}

#endif
