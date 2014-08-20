//
//   linearlsq.cpp
//
// Linear least squares class
//
// Weighted least squares
//
// Uses clipper matrix inversion: documentation indicates that this
// is suitable for < 20 parameters. No check is made on this
//


#include "linearlsq.hh"

#include <assert.h>
#define ASSERT assert

namespace scala {
  //--------------------------------------------------------------
  LinearLSQ::LinearLSQ(const int Nparam) {
    clear(Nparam);
  }
  //--------------------------------------------------------------
  void LinearLSQ::clear(const int Nparam) {
    npar = Nparam;
    AA = clipper::Matrix<double>(npar, npar, 0.0);
    ATy = std::vector<double>(npar, 0.0);
    nobs = 0;
  }
  //--------------------------------------------------------------
  // Add in observation:
  //   y   observed value
  //   x   measurement vector (the first constant
  //       parameter is explicit and should = 1.0,
  //       so x.size() == Nparam)
  //   w   sqrt(weight)
  void LinearLSQ::add(const double& y,
           const std::vector<double> x, const double& w) {
    ASSERT (int(x.size()) == npar);
    nobs++;
    for (int l=0;l<npar;l++) {
      ATy[l] += w * x[l] * y;          // contribution ot [A]T y
      for (int m=0;m<npar;m++) {
        AA(l,m) += w * x[l] * x[m];
      }
    }
  }
  //--------------------------------------------------------------
  // Return solution = parameter vector
  std::vector<double> LinearLSQ::solve() {
    if (nobs <= 1) return std::vector<double>(npar, 0.0);
    return AA.solve(ATy);
  }
  //--------------------------------------------------------------
}
