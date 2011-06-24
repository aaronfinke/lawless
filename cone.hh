// cone.hh
//
// A class to define cones around three orthogonal axes in reciprocal space
// 

#ifndef CONE_HEADER
#define CONE_HEADER

#include "hkl_datatypes.hh"

namespace scala {
  class Cone
  {
  public:
    Cone();
    Cone(const double& angledegrees);

    //! set cone angle in degrees
    void SetConeAngle(const double& angledegrees);
    //! return cone angle in degrees
    double ConeAngle() const;

    // return nearest axis index, if within cone around that axis, -1 if outside
    int Axis(const Hkl& hkl, const Rtype& invresolsq,
	     const Scell& cell);

  private:
    double coneangle; // radians
  };
} // namespace scala


#endif
