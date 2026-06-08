// refinetargets.hh
//
// Labels for possible refine target functions R,
//   based on some weighted difference
//  D = [ w (yobs - f(ycalc) ]
//  
//    QUADRATIC     least-squares R = Sum (D^2)
//    LNCOSH        R = ln cosh (D)   (recommended by Garib Murshudov)
//
//  
//  Function to minimise is either:
//    1) least squares, ie R = Sum [(wd)^2]  if targettype == RefineTargets::QUADRATIC, or
//    
//    2) R = Sum(ln(cosh(d)) as recommended by Garib Murshudov as more robust
//       targettype = RefineTargets::LNCOSH
//
//  First (gradient) and second (Hessian) derivatives:
//  1) Quadratic
//     dR/dpi      = Sum[ 2 w^2 d dd/dpi ]
//     d2R/dpidpj  = Sum[ 2 w^2 dd/dpi dd/dpj ]
//  2) ln cosh
//     dR/dpi      = Sum[ w tanh(wd) dd/dpi ]
//     d2R/dpidpj  = Sum[ (w/cosh^2(wd)) dd/dpi dd/dpj ] (strictly) ...
//            ... but Garib says this is better replaced by
//     d2R/dpidpj  = Sum[ w (1/wd) tanh(wd) dd/dpi dd/dpj ]
//     Note that if |wd| > ~3, ln(cosh(wd)) ~= |wd|, tanh(wd) ~= +-1
//               if |wd| very small (or 0), (1/wd) tanh(wd) = +1


#ifndef REFINETARGETS_HEADER
#define REFINETARGETS_HEADER

#include <string>

namespace scala {
  class RefineTargets {
  public:
    enum REFINETARGETTYPES {QUADRATIC, LNCOSH};

    static std::string format(const REFINETARGETTYPES& type);

    // maximum value of argument to cosh(x) before using approximation
    static const double MAXCOSHARG; //  = 3.0;
    // if wd < MINCOSHARG, 2nd derivative = 1
    static const double MINCOSHARG; //  = 0.01;


  };
}
#endif
