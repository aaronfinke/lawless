// simpleminimise.hh

// A simple Gauss-Newton minimiser, using Clipper classes

#ifndef SIMPLEMINIMISE_HEADER
#define SIMPLEMINIMISE_HEADER

#include <vector>

// Clipper
#include <clipper/clipper.h>

namespace SimpleMinimise {
  // ---------------------------------------------------------
  class TGH {  // Target Gradient Hessian
  public:
    TGH(){}
    TGH(const double& Target, const std::vector<double>& Gradient,
	const clipper::Matrix<double>& Hessian)
      : target(Target), gradient(Gradient), H(Hessian) {}

    double target;
    std::vector<double> gradient;
    clipper::Matrix<double> H;    // Hessian
  };
  // ---------------------------------------------------------
  class FitBase {
  public:

    FitBase(){}

    virtual int Nparameters() = 0;
    virtual int Ndata() = 0;
    
    virtual void ApplyShifts(const std::vector<double> shifts) = 0;

    virtual TGH TargetGradientHessian() = 0;
  };
  // ---------------------------------------------------------
  class DampedGaussNewton {
    // Uses Clipper routines to iterate a Gauss-Newton algorithm, with simple damping
  public:
    DampedGaussNewton(){}
    DampedGaussNewton(FitBase& fitstuff,
		      const int& Ncycles,
		      const double& tolerance,
		      const double& damp);

    void run(FitBase& fitstuff,
	     const int& Ncycles,
	     const double& tolerance,
	     const double& damp);

  private:
  };
}
#endif
