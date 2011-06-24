// cone.cpp
//
// A class to define cones around three orthogonal axes in reciprocal space
// 

// Clipper
#include <clipper/clipper.h>

#include "cone.hh"

namespace scala {
// ------------------------------------------------------------
  Cone::Cone()  {
    //    const double CONEDEGREES = 30.0;  // default angle
    const double CONEDEGREES = 20.0;  // default angle
    coneangle = clipper::Util::d2rad(CONEDEGREES);
  }
// ------------------------------------------------------------
  Cone::Cone(const double& angledegrees) {
    SetConeAngle(angledegrees);
  }
  // ------------------------------------------------------------
  //! set cone angle in degrees
  void Cone::SetConeAngle(const double& angledegrees)
  {
    coneangle = clipper::Util::d2rad(angledegrees);
  }  
  // ------------------------------------------------------------
  //! return cone angle in degrees
  double Cone::ConeAngle() const
  {
    return clipper::Util::rad2d(coneangle);
  }
  // ------------------------------------------------------------
  int Cone::Axis(const Hkl& hkl, const Rtype& invresolsq,
		 const Scell& cell)
  // return nearest reciprocal axis index,
  // if within cone around that axis, -1 if outside
  {
    DVect3 x = hkl.orth(cell);  // orthogonalised coordinate
    // dot products of x with a*,b*,c* vectors
    DVect3 xabc = x * cell.Bmat();
    double dstar = sqrt(invresolsq);  // 1/d
    std::vector<Dtype> reccell = cell.ReciprocalCell();
    int jaxis = -1;  // nearest axis, if with cone
    double minang = 10000.0; // angle
    for (int i=0;i<3;++i) { //angles from each axis
      // dot products have length = d* a*, d* b*, d* c*
      double ang = acos(Min(1.0,std::abs(xabc[i])/(dstar*reccell[i])))  ;
      if (ang < coneangle && ang < minang) {
	jaxis = i;
	minang = ang;
      }
    }
    return jaxis;
  }
}
