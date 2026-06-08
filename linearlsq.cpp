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
      ATy[l] += w * x[l] * y;          // contribution to [A]T y
      for (int m=0;m<npar;m++) {
        AA(l,m) += w * x[l] * x[m];
      }
    }
  }
  //--------------------------------------------------------------
  // Return solution = parameter vector
  std::vector<double> LinearLSQ::solve() {
    if (nobs <= 1) return std::vector<double>(npar, 0.0);
    // Scale normal matrix by diagonal U
    std::vector<double> U(npar);    // diagonal of scaling matrix
    std::vector<double> Uinv(npar); // diagonal of inverse scaling matrix
    clipper::Matrix<double> UB(npar,npar);
    for (int i=0;i<npar;++i) {
      ASSERT (AA(i,i) > 0.0);  // positive definite
      U[i] = sqrt(AA(i,i));
      Uinv[i] = 1.0/U[i];
    }
    // Scaled gradient vector
    std::vector<double> UAT(npar, 0.0);
    for (int j=0;j<npar;++j) {
      UAT[j] = ATy[j] * Uinv[j]; // scaled gradient
      for (int i=0;i<npar;++i) {
        UB(i,j) = AA(i,j) * Uinv[i] * Uinv[j];
      }}

    //^^
    /*
    for (int l=0;l<npar;l++) {
      printf("Par %2d, b %8.1f  :", l, UAT[l]);
      for (int m=0;m<npar;m++) {
        printf("  %8.1f", UB(l,m));
      }
      printf("\n");
    }
    clipper::Matrix<double> BB = UB;
    std::vector<double> ev = BB.eigen();
    for (size_t k=0; k<ev.size(); k++) {
      std::cout <<"  "<<ev[k];
    }
    std::cout <<"\n";
    //-^^
    */
    std::vector<double> pp = UB.solve(UAT);  // scaled parameters
    for (int i=0;i<npar;++i) {
      pp[i] *= Uinv[i]; // unscale
    }
    return pp;
  }
  //--------------------------------------------------------------
}
