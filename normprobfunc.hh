// normprobfunc.hh
//
// Inverse normal probability class
//
//  For a Gaussian distribution,
//  the fraction of data lying within x sds of mean is
//  P(x) (1/2pi) Integral(-x->x) exp(-a^2/2) da
//
// This class tabulates P(x), then provides a reverse lookup
// to get x given P(x)
//

#ifndef NPRBFUNC_HEADER
#define NPRBFUNC_HEADER

#include <vector>


class NormalProbability
{
public:
  NormalProbability();

  // fractional rank FracRank = (n - 2i + 1)/n
  // < 0 if i > n/2
  // rank is from 1 -> Nobs
  double ExpectedDelta(const int& rank, const int& Nobs);
  double ExpectedDelta(const double& FracRank);

  void Print();

private:
  static const int ntab;   // do not change this
  std::vector<double> prob;

  static const double TabInterval;  // 4/511
};

#endif
