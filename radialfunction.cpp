// radialfunction.cpp

#include "radialfunction.hh"

namespace scala {
  // ------------------------------------------------------------
  RadialTanhFunction::RadialTanhFunction(const std::vector<double>& params)
  // A suitable function for modelling the fall-off of
  //   CC(1/2) vs. s^2, suggested by Ed Pozharski
  // CC(s) = (1/2)(1-tanh((s-d0)/r))*dcc - dcc + 1
  //   where d0 is the value of s when CC = 0.5, and r is the steepness
  //   of fall-off, dcc is an (optional) offset to allow a fit to
  //   negative values
  // No theoretical justification!
  {
    init(params);
  }
  // ------------------------------------------------------------
  void RadialTanhFunction::init(const std::vector<double>& params)
  {
    d0 = params[0];
    r = params[1];
    dcc = 1.0;
    npar = params.size();
    if (npar == 3) {
      dcc = params[2];
    }
  }
  // ------------------------------------------------------------
  void RadialTanhFunction::SetParameters(const std::vector<double>& params)
  {init(params);}
  // ------------------------------------------------------------
  //!  value of function
  double RadialTanhFunction::value(const double& s) const
  {
    double z = (s-d0)/r;
    return 0.5*(1.0 - tanh(z))*dcc - dcc + 1.0;;
  }
  // ------------------------------------------------------------
  //! derivatives
  std::vector<double> RadialTanhFunction::deriv(const double& s) const
  {
    double z = (s-d0)/r;
    // d(tanh(z))/dz = 1/cosh^2(z)
    double dvdz = -0.5*dcc/(cosh(z)*cosh(z));
    std::vector<double> dvdp(npar);
    // dz/dd0 = -1/r
    dvdp[0] = -dvdz/r;
    // dz/dr = -z/r
    dvdp[1] = dvdp[0] * z;
    if (npar == 3) {
      dvdp[2] = 0.5*(1.0 - tanh(z)) - 1.0;
    }
    return dvdp;
  }
  // ------------------------------------------------------------
  double RadialTanhFunction::inverse(const double& v) const
  //! return s corresponding to function value v, or -1.0 if invalid
  {
    // tanh(z) = 1 - 2v
    double a = 1.0 - 2.0*v;
    if (npar == 3) {
      a = 1.0 - 2.0*(v + dcc - 1.0)/dcc;
    }
    if (std::abs(a) > 0.99999) {
      //      std::cout << "inverse " <<a<<" "<<v<<" "<<dcc<<"\n";
      return -1.0;
    }
    //    ASSERT (std::abs(a) < 1.0);
    double z = atanh(a);
    //^
    //    std::cout << "Inverse: "<<" v = "<<v <<
    //      " z= "<<z<<" a " <<a<<" d0 "<<d0<<" r "<<r;
    //      if (npar == 3) {
    //  std::cout << " dcc " <<dcc;
    //      }
    //    std::cout <<std::endl;
    //^-
    return z*r + d0;
  }
  // ------------------------------------------------------------
  int RadialTanhFunction::Nparameters() const
  {
    return npar;
  }
  // ------------------------------------------------------------
  std::vector<double> RadialTanhFunction::parameters() const
  {
    std::vector<double> p(npar);
    p[0] = d0;
    p[1] = r;
    if (npar == 3) {
      p[2] = dcc;
    }
    return p;
  }
  // ------------------------------------------------------------
  std::string RadialTanhFunction::format() const
  {
    std::string s = "d0 = " + clipper::String(d0) +
      ", r = " + clipper::String(r);
    if (npar == 3) {
      s += ", dcc = "+clipper::String(dcc);
    }
    return s;
  }
  // ------------------------------------------------------------
  std::string RadialTanhFunction::formatfunction() const
  {
    std::string s;
    if (npar == 2) {
      s = "(1/2)(1 - tanh((s - d0)/r))";
    } else if (npar == 3) {
      s = "(1/2)(1 - tanh((s - d0)/r))*dcc - dcc + 1";
    }
    return s;
  }
  // ------------------------------------------------------------
}
