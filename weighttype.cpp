// weighttype.cpp

#include "weighttype.hh"
namespace scala {
  // ------------------------------------------------------------
  //! return formatted version of weight
  std::string WeightType::formatWeightType
  (const WeightType::AverageWeightType& weighttype)
  // Weight type for averaging
  //   UNIT        unit weights
  //   VARIANCE    weight = 1/variance
  //   SQRTSCALE   weight = 1/sqrt(g)  g = 1/scale
  {
    if (weighttype == WeightType::UNIT) {return "unit weights";}
    if (weighttype == WeightType::VARIANCE) {return "variance weights";}
    if (weighttype == WeightType::SQRTSCALE) {return "SquareRoot(scale) weights";}
    if (weighttype == WeightType::SCALE) {return "scale weights";}
    return "";
  }
}
