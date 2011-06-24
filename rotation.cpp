/* rotation.cpp: fundamental data types for the clipper libraries */
// Moved into scala namesapce Phil Evans
//C Copyright (C) 2000-2004 Kevin Cowtan and University of York
//C Copyright (C) 2000-2005 Kevin Cowtan and University of York
//L
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA


#include "rotation.hh"


namespace scala {


template<int T> Euler<T>::Euler( const Rotation& rot )
{
  int r1, r2, r3, s; params( r1, r2, r3, s );
  ftype ca, cb, cg, sa, sb, sg;
  if ( r1 == r3 ) {
    ftype w, x, y, z;
    w = rot.w();
    if      ( r1==2 && r2==1 && r3==2 ) { x= rot.x(); y= rot.y(); z= rot.z(); }
    else if ( r1==1 && r2==0 && r3==1 ) { x= rot.z(); y= rot.x(); z= rot.y(); }
    else if ( r1==0 && r2==2 && r3==0 ) { x= rot.y(); y= rot.z(); z= rot.x(); }
    else if ( r1==2 && r2==0 && r3==2 ) { x=-rot.y(); y= rot.x(); z= rot.z(); }
    else if ( r1==1 && r2==2 && r3==1 ) { x=-rot.x(); y= rot.z(); z= rot.y(); }
    else if ( r1==0 && r2==1 && r3==0 ) { x=-rot.z(); y= rot.y(); z= rot.x(); }
    else Message::message( Message_fatal( "Rotation::euler() type invalid" ) );
    cb = 1.0 - 2.0 * (x*x + y*y);
    sb = 2.0 * sqrt( (x*x + y*y) * (w*w + z*z) );
    if ( sb > 0.0001 ) {
      ca = 2.0 * (x*z + w*y);
      sa = 2.0 * (y*z - w*x);
      cg = 2.0 * (w*y - x*z);
      sg = 2.0 * (y*z + w*x);
    } else {
      ca = 1.0;
      sa = 0.0;
      cg = cb;
      sg = 2.0*(y*z + w*x);
    }
    if ( s==1 ) (*this) = Euler<T>( atan2(sg,cg), atan2(sb,cb), atan2(sa,ca) );
    else        (*this) = Euler<T>( atan2(sa,ca), atan2(sb,cb), atan2(sg,cg) );
  } else {
    ftype w, x, y, z;
    w = rot.w();
    if      ( r1==0 && r2==1 && r3==2 ) { x= rot.x(); y= rot.y(); z= rot.z(); }
    else if ( r1==2 && r2==0 && r3==1 ) { x= rot.z(); y= rot.x(); z= rot.y(); }
    else if ( r1==1 && r2==2 && r3==0 ) { x= rot.y(); y= rot.z(); z= rot.x(); }
    else if ( r1==1 && r2==0 && r3==2 ) { x=-rot.y(); y= rot.x(); z= rot.z(); }
    else if ( r1==0 && r2==2 && r3==1 ) { x=-rot.x(); y= rot.z(); z= rot.y(); }
    else if ( r1==2 && r2==1 && r3==0 ) { x=-rot.z(); y= rot.y(); z= rot.x(); }
    else Message::message( Message_fatal( "Rotation::euler() type invalid" ) );
    if ( s == 0 ) x = -x;
    sa = 2.0 * ( x*y + z*w );
    ca = x*x - y*y - z*z + w*w;
    sb = 2.0 * ( y*w - x*z );
    cb = sqrt( 1.0 - sb*sb );
    sg = 2.0 * ( z*y + x*w );
    cg = z*z - x*x - y*y + w*w;
    if ( ( 1 + r1 - r2 + s ) % 3 == 0 ) sg = -sg;
    (*this) = Euler<T>( atan2(sg,cg), atan2(sb,cb), atan2(sa,ca) );
  }
}

template<int T> Rotation Euler<T>::rotation() const
{
  int r1, r2, r3, s; params( r1, r2, r3, s );
  ftype x[3];

  x[0] = x[1] = x[2] = 0.0;
  x[r1] = sin(0.5*alpha());
  Rotation rot1( cos(0.5*alpha()), x[0], x[1], x[2] );

  x[0] = x[1] = x[2] = 0.0;
  x[r2] = sin(0.5*beta() );
  Rotation rot2( cos(0.5*beta() ), x[0], x[1], x[2] );

  x[0] = x[1] = x[2] = 0.0;
  x[r3] = sin(0.5*gamma());
  Rotation rot3( cos(0.5*gamma()), x[0], x[1], x[2] );

  if ( s == 1 ) return rot1*(rot2*rot3);
  else          return rot3*(rot2*rot1);
}

template<int T> String Euler<T>::format() const
{
  char xyz[] = {'X','Y','Z'}; char rs[] = {'r','s'};
  int r1, r2, r3, s; params( r1, r2, r3, s );
  return String("Euler") + xyz[r1] + xyz[r2] + xyz[r3] + rs[s] + " = ("+String(Util::rad2d(alpha()))+","+String(Util::rad2d(beta()))+","+String(Util::rad2d(gamma()))+")";
}

template<int T> void Euler<T>::params( int& r1, int& r2, int& r3, int& s )
{
  r1 = (      ((T>>3)&3)     ) % 3;
  r2 = ( r1 + ((T>>2)&1) + 1 ) % 3;
  r3 = ( r2 + ((T>>1)&1) + 1 ) % 3;
  s  = ( T&1 );
}


String Euler_ccp4::format() const
{ return "Euler = ("+String(Util::rad2d(alpha()))+","+String(Util::rad2d(beta()))+","+String(Util::rad2d(gamma()))+")"; }


String FormatV3(const Vec3<ftype>& v)
// format axis vector for printing
{
  String s ="["+String(v[0],6,3)+", "+String(v[1],6,3)+", "+String(v[2],6,3)+"]";
  String ss;
  // Strip out spaces
  for (size_t i=0;i<s.size();i++)
    if (s[i] != ' ') ss.push_back(s[i]);
  return ss;
}

Euler_explicit::Euler_explicit(const Vec3<ftype>& e1, const ftype& phi1,
			       const Vec3<ftype>& e2, const ftype& phi2,
			       const Vec3<ftype>& e3, const ftype& phi3)
      : e1_(e1), e2_(e2), e3_(e3), phi1_(phi1), phi2_(phi2), phi3_(phi3)
{
  // Normalise direction cosines
  e1_ = e1_.unit();
  e2_ = e2_.unit();
  e3_ = e3_.unit();
}

String Euler_explicit::format() const
{ return "Explicit Euler ="+
    String(Util::rad2d(phi1()))+" "+FormatV3(e1())+";"+
    String(Util::rad2d(phi2()))+" "+FormatV3(e2())+";"+
    String(Util::rad2d(phi3()))+" "+FormatV3(e3())
    ; }

String Polar_ccp4::format() const
{ return "Polar = ("+String(Util::rad2d(omega()))+","+String(Util::rad2d(phi()))+","+String(Util::rad2d(kappa()))+")"; }

Rotation::Rotation( const Euler_ccp4& euler )
{
  // not optimised
  ftype a1 = 0.5 * euler.alpha();
  ftype a2 = 0.5 * euler.beta();
  ftype a3 = 0.5 * euler.gamma();
  Rotation r1( cos(a1), 0.0, 0.0, sin(a1) );
  Rotation r2( cos(a2), 0.0, sin(a2), 0.0 );
  Rotation r3( cos(a3), 0.0, 0.0, sin(a3) );
  *this = r3*(r2*r1);   // gamma, then beta, then alpha
}

//! constructor: from Euler_explicit
Rotation::Rotation( const Euler_explicit& euler_explicit )
{
  // rotate first r3, then r2, then r1
  Rotation r1(euler_explicit.e1(), euler_explicit.phi1());
  Rotation r2(euler_explicit.e2(), euler_explicit.phi2());
  Rotation r3(euler_explicit.e3(), euler_explicit.phi3());
  *this = r3*(r2*r1);
}


Rotation::Rotation( const Polar_ccp4& polar )
{
  w_ = cos( 0.5 * polar.kappa() );
  x_ = sin( 0.5 * polar.kappa() ) * cos( polar.phi() ) * sin( polar.omega() );
  y_ = sin( 0.5 * polar.kappa() ) * sin( polar.phi() ) * sin( polar.omega() );
  z_ = sin( 0.5 * polar.kappa() ) *                      cos( polar.omega() );
}

//! constructor from direction cosines and angle
Rotation::Rotation( const Vec3<ftype>& axis, const ftype& kappa )
{
  Vec3<ftype> ax = axis.unit();
  w_ = cos( 0.5 * kappa );
  x_ = sin( 0.5 * kappa ) * ax[0];
  y_ = sin( 0.5 * kappa ) * ax[1];
  z_ = sin( 0.5 * kappa ) * ax[2];
}

Rotation::Rotation( const Mat33<>& m )
{
  ftype tr = m(0,0) + m(1,1) + m(2,2) + 1.0;

  // check the diagonal
  if ( tr > 1.0e-8 ) {
    ftype s( sqrt(tr) );
    w_ = s * 0.5;
    s = 0.5 / s;
    x_ = s * ( m(2,1) - m(1,2) );
    y_ = s * ( m(0,2) - m(2,0) );
    z_ = s * ( m(1,0) - m(0,1) );
  } else {
    if ( m(0,0) > m(1,1) && m(0,0) > m(2,2) ) {
      ftype s( sqrt(1.0 + m(0,0) - m(1,1) - m(2,2) ) );
      x_ = 0.5 * s;
      if ( s != 0.0 ) s = 0.5 / s;
      w_ = s * ( m(2,1) - m(1,2) );
      y_ = s * ( m(0,1) + m(1,0) );
      z_ = s * ( m(0,2) + m(2,0) );
    } else if ( m(1,1) > m(2,2) ) {
      ftype s( sqrt(1.0 + m(1,1) - m(2,2) - m(0,0) ) );
      y_ = 0.5 * s;
      if ( s != 0.0 ) s = 0.5 / s;
      w_ = s * ( m(0,2) - m(2,0) );
      z_ = s * ( m(1,2) + m(2,1) );
      x_ = s * ( m(1,0) + m(0,1) );
    } else {
      ftype s( sqrt(1.0 + m(2,2) - m(0,0) - m(1,1) ) );
      z_ = 0.5 * s;
      if ( s != 0.0 ) s = 0.5 / s;
      w_ = s * ( m(1,0) - m(0,1) );
      x_ = s * ( m(2,0) + m(0,2) );
      y_ = s * ( m(2,1) + m(1,2) );
    }
  }
}

/*! The normalisation is performed in-place. If a rotation becomes
  significantly denormalised, the conversion methods will
  fail. Therefore it may be safer to call this before a conversion. */
const Rotation& Rotation::norm()
{
  ftype s = w_*w_+x_*x_+y_*y_+z_*z_;
  if ( s < 1.0e-12 ) {
    w_ = 1.0;
    x_ = y_ = z_ = 0.0;
  } else {
    s = 1.0/sqrt(s);
    w_ *= s;
    x_ *= s;
    y_ *= s;
    z_ *= s;
  }
  return *this;
}

/*! If beta ~= 0, then alpha is set to zero.
  \return The Euler_ccp4 angles. */
Euler_ccp4 Rotation::euler_ccp4() const
{
  ftype ca, cb, cg, sa, sb, sg;
  cb = 1.0 - 2.0 * (x_*x_ + y_*y_);
  sb = 2.0 * sqrt( (x_*x_ + y_*y_) * (w_*w_ + z_*z_) );
  if ( sb > 0.0001 ) {
    ca = 2.0 * (x_*z_ + w_*y_);
    sa = 2.0 * (y_*z_ - w_*x_);
    cg = 2.0 * (w_*y_ - x_*z_);
    sg = 2.0 * (y_*z_ + w_*x_);
  } else {
    ca = 1.0;
    sa = 0.0;
    cg = cb;
    sg = 2.0*(y_*z_ + w_*x_);
  }
  return Euler_ccp4( atan2(sa,ca), atan2(sb,cb), atan2(sg,cg) );
}

/*! If omega ~= 0, then phi is set to zero.
  \return The Polar_ccp4 angles. */
Polar_ccp4 Rotation::polar_ccp4() const
{
  ftype om, ph, ka;
  om = ph = ka = 0.0;
  if ( w_ < 0.999999 ) {
    ftype r = sqrt( x_*x_ + y_*y_ );
    om = atan2( r, z_ );
    if ( r > 0.000001 ) ph = atan2( y_, x_ );
    ka = 2.0*acos( w_ );
  }
  return Polar_ccp4( om, ph, ka );
}

/*! The resulting rotation matrix would commonly be used to construct
  a clipper::RTop_orth.
  \return The rotation matrix. */
Mat33<> Rotation::matrix() const
{
  ftype xx( 2.0*x_*x_ );
  ftype yy( 2.0*y_*y_ );
  ftype zz( 2.0*z_*z_ );
  ftype xy( 2.0*x_*y_ );
  ftype yz( 2.0*y_*z_ );
  ftype xz( 2.0*z_*x_ );
  ftype wx( 2.0*w_*x_ );
  ftype wy( 2.0*w_*y_ );
  ftype wz( 2.0*w_*z_ );
  return Mat33<>( 1.0-yy-zz, xy-wz,     xz+wy,
		       xy+wz,     1.0-xx-zz, yz-wx,
		       xz-wy,     yz+wx,     1.0-xx-yy );
}

  ftype AngleInRange(const ftype& angle)
  {
    // Put angle in range -pi to +pi
    if (angle >  Util::pi()) return (angle - Util::twopi());
    if (angle < -Util::pi()) return (angle + Util::twopi());
    return angle;
  }

  bool Rotation::euler_explicit(Euler_explicit& euler_explicit,
				const int& SolutionNumber) const
//!< return Euler_explicit angles given axes; return false if no solution
//!< SolutionNumber = 1 or 2 for the two solutions
//
//  If the angle between e1 & e2 or e2 & e3 are < 90 degrees, then
//  there are matrices which have no solution
{
  //  Note that the algorithm for decomposing a rotation into three arbitrary
  //  rotations was described by David Thomas
  //  (Thomas, D.J.(1990) Acta Cryst. A46, 321-343)
  //  and detailed by Gerard Bricogne
  //  (Proceedings of the CCP4 Study Weekend 23-24 January 1987)
  //  The code here is a essentially translation of a Fortran routine
  //  by Zbyszek Otwinowski  (solveu)
  //
  //  The notation used follows Bricogne
  //
  //     The Euler_explicit class represents three angles of rotation around
  //     specified axes e1, e2, e3 defined by their direction cosines
  //     ie rotate about e3, then e2, then e1
  //
  //     R =  R(e1, phi1) R(e2, phi2) R(e3, phi3) 
  //
  //     where R(e, phi) is the rotation matrix for a rotation by angle
  //     phi around axis e
  //----
  // Preliminaries:-
  // The rotation matrix itself
  Mat33<ftype> R = matrix();
  Vec3<ftype> e1 = euler_explicit.e1();
  Vec3<ftype> e2 = euler_explicit.e2();
  Vec3<ftype> e3 = euler_explicit.e3();

  // Dot products of the three vectors
  ftype e1e2 = e1 * e2;
  ftype e1e3 = e1 * e3;
  ftype e2e3 = e2 * e3;

  ftype e1e2e3 = e1 *
    Vec3<ftype>::cross(e2,e3);  // triple product

  // R e3
  Vec3<ftype> Re3 = R * e3;
  // e1 . Re3
  ftype e1Re3 = e1 * Re3;

  // We need a test vector perpendicular to e3: use e2 x e3
  Vec3<ftype> u = Vec3<ftype>::cross(e2, e3);
  if (u*u < 1.0e-6)
    {
      // Fail if e2 & e3 are parallel
      return false;
    }
  u = u.unit();

  // ** Step 1 ** Calculation of phi2  (Bricogne equation (4))
  // e1.(R e3) = (e1.e2)(e2.e3) + {(e1.e3) - (e1.e2)(e2.e3)} cos(phi2)
  //           + (e1.e2 x e3) sin(phi2)
  // The coefficients of cos & sin phi2
  ftype cc = e1e3 - e1e2 * e2e3;
  ftype ss = e1e2e3;
  // Error if both are zero (indeterminate)
  if (std::abs(cc) < 1.0e-6 && std::abs(ss) < 1.0e-6) {return false;}
  // Normalise equation (4)
  ftype norm = sqrt(cc*cc + ss*ss);
  // rhs = e1.Re3 - (e1.e2)(e2.e3)
  ftype rhs = (e1Re3 - e1e2 * e2e3)/norm;
  // abs(rhs) should not be greater than 1.0
  if (std::abs(rhs) > 1.000002) {return false;}
  // Allow a small tolerance
  if (rhs > 1.0)
    {rhs = 1.0;}
  else if (rhs < -1.0)
    {rhs = -1.0;}
  cc /= norm;
  ss /= norm;
  // Solve  rhs  = cos(phi2) * cc + sin(phi2) * ss
  // using cos(a-b) = cos(a) cos(b) + sin(a) sin(b)
  //   where b = phi2
  ftype a = atan2(ss, cc);
  ftype amb = acos(rhs);  // +-(a-b)
  // Two solutions (from arc cos), phi2 = b = a -+(a-b) 
  //   1)  phi2 = a - amb
  //   2)  phi2 = a + amb
  // in range -pi to +pi
  // Note that if e1 == e3, ss = 0, a = 0 & phi2b = -phi2a
  ftype phi2a  = AngleInRange(a - amb);
  ftype phi2b  = AngleInRange(a + amb);
  // Choose solution 1 (larger phi2 ie positive) or
  // solution 2 (smaller (negative) phi2)
  //  make phi2a > phi2b  
  if (phi2a < phi2b) Util::swap<ftype>(phi2a, phi2b);
  ftype phi2 = phi2b;
  if (SolutionNumber == 1) phi2 = phi2a;   // pick one solution
  
  // ** Step 2 ** Calculation of phi1
  ftype phi1 = 0.0;
  Rotation R2(e2, phi2);  // R2 = R(e2, phi2)
  Vec3<ftype> v = R2.matrix() * e3;  // v = R2 e3 
  // w = R e3
  Vec3<ftype> w = Re3;
  // v1 = v - (v.e1) e1
  // w1 = w - (w.e1) e1
  Vec3<ftype> v1 = v - (v*e1)*e1;
  Vec3<ftype> w1 = w - (w*e1)*e1;
  norm = (v1*v1)*(w1*w1);
  // If norm = 0, rotations 1 & 3 are around same axis (for this phi2),
  // so any value for phi1 is OK (leave = 0.0)
  if (norm > 1.0e-8)
    {
      norm = sqrt(norm);
      // norm = sqrt((v1.v1)(w1.w1))
      // cos(phi1) = (v1.w1)/norm
      // sin(phi1) = (v1.w1 x e1)/norm
      phi1 = AngleInRange(atan2((v1*Vec3<ftype>::cross(w1, e1))/norm, (v1*w1)/norm));
    }
  
  // ** Step 3 ** Calculation of phi3
  // R3 = R(e3, phi3) = R * R(e1, -phi1) * R(e2, -phi2) 
  //   note rotations mulitply in reverse order to matrices
  Rotation R3 =  *(this) * (Rotation(e1, -phi1) * R2.inverse());
  Vec3<ftype> R3u = R3.matrix() * u;
  //  sin(phi3) = u. R3u x e3
  //  cos(phi3) = u . R3u
  ftype phi3 = atan2(u*Vec3<ftype>::cross(R3u, e3), u*R3u);
  euler_explicit.SetAngles(phi1, phi2, phi3);
  return true;
}


/*! Note: This multiplication operator combines rotations in the
  opposite order to that of matrices. Thus, the rotation which arises
  from applying rotation r1 followed by rotation r2 is given by
  r1*r2. Similarly, the rotation which arises from applying rotation
  r1 followed by rotation r2 and rotation r3 is given by
  r1*(r2*r3). */
Rotation operator* ( const Rotation& r1, const Rotation& r2 )
{
  return Rotation( r1.w_*r2.w_ - r1.x_*r2.x_ - r1.y_*r2.y_ - r1.z_*r2.z_,
		   r1.w_*r2.x_ + r1.x_*r2.w_ + r1.z_*r2.y_ - r1.y_*r2.z_,
		   r1.w_*r2.y_ + r1.y_*r2.w_ + r1.x_*r2.z_ - r1.z_*r2.x_,
		   r1.w_*r2.z_ + r1.z_*r2.w_ + r1.y_*r2.x_ - r1.x_*r2.y_ );
}

String Rotation::format() const
{ return "Quaternion wxyz = ("+String(w_)+","+String(x_)+","+String(y_)+","+String(z_)+")"; }


template class Euler< 0>;
template class Euler< 1>;
template class Euler< 2>;
template class Euler< 3>;
template class Euler< 4>;
template class Euler< 5>;
template class Euler< 6>;
template class Euler< 7>;
template class Euler< 8>;
template class Euler< 9>;
template class Euler<10>;
template class Euler<11>;
template class Euler<12>;
template class Euler<13>;
template class Euler<14>;
template class Euler<15>;
template class Euler<16>;
template class Euler<17>;
template class Euler<18>;
template class Euler<19>;
template class Euler<20>;
template class Euler<21>;
template class Euler<22>;
template class Euler<23>;


} // namespace clipper
