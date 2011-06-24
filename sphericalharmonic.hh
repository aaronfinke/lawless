//
// sphericalharmonic.hh
//
// Class to calculate spherical harmonics
//
// Originally taken from some Fortran code by Ian Tickle

#ifndef SPHERICALHARMONIC_HEADER
#define SPHERICALHARMONIC_HEADER

#include <iostream>
#include <vector>
#include "util.hh"

class SphericalHarmonic
{
public:
  SphericalHarmonic();     // construct with default order 6
  //  lmaxEven, lmaxOdd   maximum orders for even & odd terms
  SphericalHarmonic(const int& lmaxEven,
		    const int& lmaxOdd);

  //  Number of Yml terms
  int Nterms() const;
  int MaxOrder() const {return Max(lmaxeven, lmaxodd);}

  //!
  void Check() const {std::cout << "%%%% SphericalHarmonic: lmax " << lmax<<"\n";}

  // Calculate Ylm(theta,phi) coefficients
  //
  // On entry:
  //  theta (radians)     colatitude, in range 0 -> pi
  //  phi (radians)       longitude, 0 -> 2pi
  //  lmaxEven, lmaxOdd   maximum orders for even & odd terms
  //  
  // Returns ylm  vector of spherical harmonic parameters
  //  number of even or odd terms = lmax/2(lmax+3) + 1
  // Note that the first entry is Y00, a constant term
  //
  // Normalisation factor = 1/(4*pi).
  //  Spherical harmonics are Hermitian antisymmetric i.e.
  //    Ylm(L,-M) = (-1)**M * CONJG( Ylm(L,M) )
  // Number of terms returned = 
  std::vector<double> Ylm(const double& theta,
			  const double& phi) const;

private:
  void SetUp(const int& maxOrder);
  int lmax, lmaxeven, lmaxodd; 

  std::vector<double> fact;
};

#endif
