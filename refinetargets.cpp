// refinetargets.cpp
//
// Labels for possible refine target functions R,
//   based on some weighted difference
//  D = [ w (yobs - f(ycalc) ]
//  
//    QUADRATIC     least-squares R = Sum (D^2)
//    LNCOSH        R = ln cosh (D)   (recommended by Garib Murshudov)
//

#include "refinetargets.hh"

namespace scala {
  std::string RefineTargets::format(const REFINETARGETTYPES& type)
  {
    if (type == QUADRATIC) {
      return "quadratic";
    } else if (type == LNCOSH) {
      return "ln(cosh())";
    }
    return "";
  }
}
