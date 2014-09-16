// normprobfunc.cpp
// Inverse normal probability class
//
//  For a Gaussian distribution,
//  the fraction of data lying within x sds of mean is
//  P(x) (1/2pi) Integral(-x->x) exp(-a^2/2) da
//
// This class tabulates P(x), then provides a reverse lookup
// to get x given P(x)
//
// Translated from Fortran in Scala, probably originally from
// David Smith & Lynne Howell program NORMAN
// Note that the internal table is hard-wired to 512 elements

#include <iostream>
#include <cmath>
#include <math.h>


#include "normprobfunc.hh"

//---------------------------------------------------------
const double NormalProbability::TabInterval = double(4)/double(511);
const int NormalProbability::ntab = 512; // do not change this

NormalProbability::NormalProbability()
{
  prob.resize(ntab);

  //  Build probability integral for X = 0 to 4.0 in
  //     increments of 0.007827789.  (= 4/511)
  double p = 1.0;
  double den[28]; // indexed from 1
  double pi = 4.0 * atan(1.0);
  double TwoOvRoottwoPi = 2.0/sqrt(2.0*pi); // 0.7978845608

  for(int i=1;i<=27;i++)
    {
      p = p * i;
      den[i] = log10(p * (2 * i + 1) * pow(2.0,i));
    }

  int id = 0;
  prob[0] = 0.0;

  for (int i=1;i<ntab;i++)   // from 1->511
    {
      id = id + 1;
      double del = id;
      double p = del * TabInterval;
      double sum = 1.0;
      for (int j=1;j<=27;j++)
        {
          int k = (j+1)%2 - j%2;
          double psum = (2 * j) * log10(p) - den[j];
          if (psum + 7.0 < 0.0) break;
          psum = pow(10.0,psum);
          sum = sum + k * psum;
        }
      prob[i] = TwoOvRoottwoPi * p * sum;
    }
}
//---------------------------------------------------------
double NormalProbability::ExpectedDelta(const int& rank, const int& Nobs)
// fractional rank FracRank = (n - 2i + 1)/n
// < 0 if i > n/2
{
  return ExpectedDelta((double(Nobs)-2.0*rank+1.0)/double(Nobs));
}
//---------------------------------------------------------
double NormalProbability::ExpectedDelta(const double& FracRank)
// fractional rank FracRank = (n - 2i + 1)/n
// < 0 if i > n/2
{
  // interpolate from probability table

  double px = std::abs(FracRank);

  int i = 256;   // = ntab/2 - 1
  int id = i/2;
  //  start of binary search loop
  for (int k=0;k<8;k++)
    {
      double del = px - prob[i-1];
         if (del > 0.0)
           {
             i = i + id;
           }
         else
           {
             i = i - id;
           }
         id = id / 2;
    }
  if (px < prob[i-1]) i = i - 1;
  double del = (px - prob[i-1]) / (prob[i] - prob[i-1]);
  double x   = TabInterval * (double(i-1) + del); // 4/511
  if (FracRank > 0.0) x = -x;
  return x;
}
//---------------------------------------------------------
void NormalProbability::Print()
{
  for (int i=0;i<ntab;i++)
    {
      std::cout << i << " " << i*TabInterval << " " << prob[i] << "\n";
    }
}
//---------------------------------------------------------
