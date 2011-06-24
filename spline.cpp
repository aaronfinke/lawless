#include <iostream>
#include "spline.hh"


namespace scala
{
  Spline::Spline(const std::vector<RPair>& xyin)
  {
    n = xyin.size();
    x.resize(n);
    y.resize(n);
    for (int i=0;i<n;i++)
      {
	x[i] = xyin[i].first;
	y[i] = xyin[i].second;
      }
    // Calculate 2nd derivatives
    y2.resize(n);
    std::vector<double> u(n);

    y2[0] = 0.0;
    u[0] = 0.0;

    for (int i=1;i<n-1;i++)
      {
	double sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
	double p = sig*y2[i-1] + 2.0;
	y2[i] = (sig - 1.0)/p;
	u[i] = (6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/
		(x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1])/p;
      }
    double qn=0.0;
    double un=0.0;
    y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.);

    for (int i=n-2;i>=0;i--)
      {
	y2[i] = y2[i]*y2[i+1]+u[i];
      }
  }
  float Spline::Interpolate(const float& xx) const
  {
    int lo = 0;
    int hi = n-1;
    while (hi-lo > 1)
      {
	int k=(hi+lo)/2;
	if (x[k] > xx)
	  hi = k;
	else
	  lo = k;
      }
    double h = x[hi] - x[lo];
    if (h == 0.0) return y[hi];
    double a = (x[hi]-xx)/h;
    double b = (xx-x[lo])/h;
    return a*y[lo] + b*y[hi] +
      ((a*a*a - a)*y2[lo] + (b*b*b-b)*y2[hi])*h*h/6.0;
  }
}
