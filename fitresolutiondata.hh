//  fitresolutiondata.hh
//

//  Fit a function of resolution to some binned data, using Clipper

#ifndef FITRESOLUTIONDATA_HEADER
#define FITRESOLUTIONDATA_HEADER

#include <vector>
#include "simpleminimise.hh"
#include "refinetargets.hh"
#include "radialfunction.hh"

// Clipper
#include <clipper/clipper.h>

using namespace SimpleMinimise;

namespace scala {
  // ---------------------------------------------------------
  class ResolutionData {
  public:
    ResolutionData(){}
    ResolutionData(const double& S, const double& V, const double& W)
      : s(S), v(V), w(W) {}

    double s;  // 1/d^2
    double v;  // a value at position s
    double w;  // sqrt(weight), <= 0 to ignore
  };
  // ---------------------------------------------------------
  // ---------------------------------------------------------
  class FitResolutionData : public FitBase
  {
  public:
    FitResolutionData(){}
    FitResolutionData(const std::vector<ResolutionData>& Rdata,
   		      RadialBase& Radialfunction);

    // set true to use quadratic residual, else ln(cosh(d))
    void setQuadratic(const bool& quadratic);

    int Nparameters() const {return npar;}
    int Ndata() const {return rdata->size();}
    
    int NvalidData() const;  // number of data with weight > 0

    void ApplyShifts(const std::vector<double> shifts);
    
    void settargetType(const RefineTargets::REFINETARGETTYPES& Targettype)
    {targettype = Targettype;}
    
    // returns target etc
    TGH TargetGradientHessian(const bool& dogradient=true);

    // just target
    double Target();

    std::vector<double> parameters() const;
    void SetParameters(const std::vector<double>& params);

  private:
    const std::vector<ResolutionData>* rdata;
    RadialBase* radialfunction;
    std::vector<double> params;
    int npar;   // number of parameters
    std::vector<int> nrefpar; // number of contributions to each parameter
    
    double target;
    int nobs;
    
    // default to use ln(cosh(d)) residual (LNCOSH), else QUADRATIC
    RefineTargets::REFINETARGETTYPES targettype;
  };
}
#endif
