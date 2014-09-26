// radialfunction.hh
//
// radial functions,
// eg of resolution s = 4(sin theta/lambda)**2


#ifndef RADIALFUNCTION_HEADER
#define RADIALFUNCTION_HEADER

#include <vector>

#include "range.hh"
#include "scala_util.hh"

namespace scala {
  //---------------------------------------------------------------
  class RadialBase {
  public:
    RadialBase(){}
    RadialBase(const std::vector<double>& params){}

    virtual void init(const std::vector<double>& params) = 0;

    virtual void SetParameters(const std::vector<double>& params) = 0;

    //!  value of function
    virtual double value(const double& s) const = 0;

    //! derivatives
    virtual std::vector<double> deriv(const double& s) const = 0;

    //! return s corresponding to function value v
    virtual double inverse(const double& v) const = 0;

    virtual int Nparameters() const = 0;

    virtual std::vector<double> parameters() const = 0;

    //    virtual std::vector<double> getLargeShifts() const = 0;

    virtual std::string format() const {return "";}

    virtual std::string formatfunction() const {return "";}

  };
  //---------------------------------------------------------------
  class RadialTanhFunction : public RadialBase {
    // A suitable function for modelling the fall-off of
    //   CC(1/2) vs. s^2, suggested by Ed Pozharski
    // CC(s) = (1/2)(1-tanh((s-d0)/r))*dcc - dcc + 1
    //   where d0 is the value of s when CC = 0.5, and r is the steepness
    //   of fall-off, dcc is an (optional) offset to allow a fit to
    //   negative values
    // No theoretical justification!
  public:
    RadialTanhFunction(){}

    //! construct from 2 or 3 parameters
    RadialTanhFunction(const std::vector<double>& params);
    //! initialise from 2 or 3 parameters
    void init(const std::vector<double>& params);

    void SetParameters(const std::vector<double>& params);

    //!  value of function
    double value(const double& s) const; 

    //! derivatives
    std::vector<double> deriv(const double& s) const;

    //! return s corresponding to function value v
    double inverse(const double& v) const;

    int Nparameters() const;

    std::vector<double> parameters() const;

    //    std::vector<double> getLargeShifts() const;

    std::string format() const;

    std::string formatfunction() const;


  private:
    double d0;
    double r;
    double dcc;  // range of CC, usually = 1.0 for CC in range 1.0 -> 0.0
    int npar;  // 2 or 3
  };

}
#endif
