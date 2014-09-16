
#include "samplegaussian.hh"
#include "scala_util.hh"
using namespace scala;

//--------------------------------------------------------------
float SampleGaussian::Get(const float& Mean, const float& SD)
{
  return SD*Get() + Mean;
}

//--------------------------------------------------------------
float SampleGaussian::Get()
// generate a random number with mean 0 & unit variance,
// using mysterious Box-Muller transformation
//   see Numerical Recipes
{
  if (gotone)
    {
      gotone = false;
      return previous;  // return buffered value
    }
  float v1, v2;
  float rsq = +1000000.;
  while (rsq > 1.)
    {
      // two random numbers in range -1 to +1
      v1 = 2.0*FRandom(1.0) - 1.0;
      v2 = 2.0*FRandom(1.0) - 1.0;
      // Radius^2
      rsq = v1*v1 + v2*v2;
    }
  float fac = sqrt(-2.*log(rsq)/rsq);
  previous = v2*fac;  // store one for later
  gotone = true;
  return v1*fac;
}
