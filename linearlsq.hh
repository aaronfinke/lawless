//
//   linearlsq.hh
//
// Linear least squares class
//
// Weighted least squares
//
// Uses clipper matrix inversion: documentation indicates that this
// is suitable for < 20 parameters. No check is made on this
//

#ifndef LINEARLSQ_HEADER
#define LINEARLSQ_HEADER

#include <vector>

// Clipper
#include <clipper/clipper.h>

namespace scala {

  class LinearLSQ {
  public:
    LinearLSQ(){}
    LinearLSQ(const int Nparam);

    void clear(const int Nparam); // Clear totals    

    // Add in observation:
    //   y   observed value
    //   x   measurement vector (the firat constant
    //       parameter is explicit and should = 1.0,
    //       so x.size() == Nparam)
    //   w   sqrt(weight)
    void add(const double& y,
	     const std::vector<double> x, const double& w);

    // Return solution = parameter vector
    std::vector<double> solve();

    int Nobs() const {return nobs;}

  private:
    int npar;
    clipper::Matrix<double> AA;  // normal matrix [A]T [A]
    std::vector<double>     ATy; // [A]T y
    int nobs;
  };

}

#endif
