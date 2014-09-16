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
  // ---------------------------------------------------------
  // maximum value of argument to cosh(x) before using approximation
  const double RefineTargets::MAXCOSHARG = 3.0;
  // if wd < MINCOSHARG, 2nd derivative = 1
  const double RefineTargets::MINCOSHARG = 0.01;
  // ---------------------------------------------------------
  std::string RefineTargets::format(const REFINETARGETTYPES& type)
  {
    if (type == QUADRATIC) {
      return "quadratic";
    } else if (type == LNCOSH) {
      return "ln(cosh())";
    }
    return "";
  }
  // ---------------------------------------------------------
}
