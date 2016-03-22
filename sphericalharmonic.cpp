//
//  sphericalharmonic.cpp
//

#include <math.h>
#include "sphericalharmonic.hh"
#include "util.hh"


//----------------------------------------------------------------------------
SphericalHarmonic::SphericalHarmonic()     // construct with default order 6
{
  lmaxeven = 6;
  lmaxodd  = 5;
  SetUp(Max(lmaxeven, lmaxodd));
}
//----------------------------------------------------------------------------
SphericalHarmonic::SphericalHarmonic(const int& lmaxEven,
                                     const int& lmaxOdd)
//  lmaxEven, lmaxOdd   maximum orders for even & odd terms
{
  lmaxeven = (lmaxEven/2)*2;   // force even
  lmaxodd  = ((lmaxOdd+1)/2)*2-1;  // force odd
  //^  std::cout << "even, odd " << lmaxeven << " " << lmaxodd << "\n";
  SetUp(Max(lmaxeven, lmaxodd));
}
//----------------------------------------------------------------------------
void SphericalHarmonic::SetUp(const int& maxOrder)
{
  // Set up log factorial table
  lmax = maxOrder;
  fact.resize(2*lmax+1);
  fact[0] = 0.;
  for (int l=1;l<=2*lmax;++l) {
    fact[l] = fact[l-1] + log(double(l));
  }
  //^  std::cout << "SphericalHarmonic::SetUp, lmax " << lmax << " " << fact.size() << "\n"; //^
}
//----------------------------------------------------------------------------
int SphericalHarmonic::Nterms() const
//  Number of Yml terms excluding 00
{
 return lmaxeven*(lmaxeven+3)/2+1 + lmaxodd*(lmaxodd+3)/2;
}
//----------------------------------------------------------------------------
std::vector<double> SphericalHarmonic::Ylm(const double& theta,
                                           const double& phi) const
//
// Calculate Ylm(theta, phi)
//
// On entry:
//  theta               colatitude, in range 0 -> pi
//  phi (radians)       longitude, 0 -> 2pi
//  lmaxEven, lmaxOdd   maximum orders for even & odd terms
//
// Returns ylm  vector of spherical harmonic parameters
//  number of even or odd terms = lmax/2(lmax+3) + 1
//
// Normalisation factor = 1/(4*pi).
// Spherical harmonics are Hermitian antisymmetric i.e.
// S(L,-M) = (-1)^M * CONJG( S(L,M) )
{
  std::vector<double> ylm;

  double plm1, plm2, plm, phm;
  double sintheta = sin(theta);
  double costheta = cos(theta);
  int m1;
  double temp;

  std::vector<double> sinthetaPowm(lmax+1);  // sin theta ^ m, m 0, lmax
  sinthetaPowm[1] = sintheta;
  for (int m=2;m<=lmax;++m) {
    sinthetaPowm[m] = sinthetaPowm[m-1] * sintheta;
  }

  // Loop l  (exclude 00 term)
  for (int l=1;l<=lmax;++l) {
    if ((l%2 == 0) ? (l <= lmaxeven) : (l <= lmaxodd)) {
      // loop l
      for (int m=l;m>=0;--m) {
        if (m == l) {
          plm1=0.;
          plm2=exp(.5*fact[2*l]-fact[l]-l*log(2.0))*sqrt(double(2*l+1));
        } else {
          m1=m+1;
          temp=plm1 * sintheta * sintheta;
          plm1=plm2;
          plm2=(2.*m1*costheta*plm2-sqrt(double((l-m1)*(l+m1+1)))*temp)/
            sqrt(double((l+m1)*(l-m)));
        }
        //    std::cout << "l,m,plm2 " << l << " " << m << " " << plm2 << "\n";
        if (m > 0) {
          //      double plm=plm2*pow(sintheta, double(m));
          plm=plm2*sinthetaPowm[m];
           phm=m*phi;
          ylm.push_back(plm*cos(phm));
          ylm.push_back(plm*sin(phm));
        } else {
          ylm.push_back(plm2);
        }
      }
    }
  }
  return ylm;
}
