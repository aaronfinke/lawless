// matvec.hh

// Various utility routines for handling 3x3 matrices & vectors
// Some routines using CCtbx routines are in cctbx_utils
//
// By default all are <double> unless there are special reasons
// for eg float
// (mainly for interface to external libraries such as CCP4)

// For various reasons, I am using up to 4 different representations
// of 3x3 matrices (+ sgtbx::rt_mx
//
//  3x3 matrices are
//     ( 0,0  0,1  0,2 )
//     ( 1,0  1,1  1,2 )     element i,j
//     ( 2,0  2,1  2,2 )
//
// a) simple arrays double[3][3]   [i][j]             AMat33
// b) clipper::Mat33<double> (i,j)                    CMat33
// c) std::vector<double> M33(9) elements in order 00, 01, 02, 10 ...
//    M33(k)  k = i*3 + j                             VMat33
//
// Routines here are:-
//  1) converters                Set[type] (other_type)
//                               (not all combinations written yet!)
//  2) tests of identity         is_mat33_ident(Matrix)
//  3) format reindex matrix     FormatReindex_as_hkl(Matrix)
//  4) vector utilities          Modulus(Vec3<T>)
//                               IntVec(Vec3<double>)
//  5) LSQ general matrix to superimpose vector set y = [U] x
//     not necessarily orthonormal (no constraints)
//

#ifndef MATVEC_UTILS
#define MATVEC_UTILS

#include <cstring>

// Clipper
#include <clipper/clipper.h>
using clipper::Vec3;
using clipper::Mat33;
using clipper::RTop;
typedef clipper::Vec3<double> DVect3;
typedef clipper::Mat33<double> DMat33;

#include "csymlib.h"    // CCP4 symmetry stuff

#include <assert.h>
#define ASSERT assert

namespace MVutil {
  //--------------------------------------------------------------
  // Make clipper::Mat33 from various things
  //   2D array
  //     float -> double
  clipper::Mat33<double>  SetCMat33(const float m[3][3]);
  clipper::Mat33<double> SetCMat33(const float m[9]);
  //     double
  clipper::Mat33<double> SetCMat33(const double m[3][3]);
  // vector<double> (9)
  clipper::Mat33<double> SetCMat33(const std::vector<double>& m);
  // CCP4 symop
  clipper::Mat33<double> SetCMat33(CSym::ccp4_symop& op);
  // Diagonal matrix
  clipper::Mat33<double> SetDiagCMat33(const clipper::Vec3<double>& v);
  //--------------------------------------------------------------
  // Convert Mat33 to vector type
  std::vector<double> SetVMat33(const clipper::Mat33<double>& R);
  std::vector<double> SetVMat33(const float m[3][3]);
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  bool is_mat33_ident(const clipper::Mat33<float>& R);
  // Returns true if 3x3 matrix R is identity
  //--------------------------------------------------------------
  bool is_mat33_ident(const clipper::Mat33<double>& R);
  // Returns true if 3x3 matrix R is identity
  //--------------------------------------------------------------
  bool is_rtop_ident(const clipper::RTop<double>& R);
  // Returns true if 3x3 matrix R is identity
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  // Format reciprocal space matrix (eg reindex or symmetry)
  // as h,k,l
  std::string FormatSymop_as_hkl(const clipper::Symop& op,
				 const std::string& brackets = "[]");
  std::string FormatReindex_as_hkl(const clipper::RTop<double>& op,
  				   const std::string& brackets = "[]");
  //--------------------------------------------------------------
  template<class T> inline T Modulus( const clipper::Vec3<T>& v)
  { return sqrt(v*v); }
  //--------------------------------------------------------------
  clipper::Vec3<int> IntVec(const clipper::Vec3<double>& vector);
  //--------------------------------------------------------------
  template<class T> int LargestVectorElement(const clipper::Vec3<T>& v)
  // Return index 0,1,2 of largest element in vector
  {
    if (v[0] > v[1]) {
      if (v[0] > v[2]) {
	return 0; // v[0] biggest
      } else {
	return 2; // v(2] biggest
      }
    } else if (v[1] > v[2]) {
      return 1; // v[1] biggest
    }
    return 2;  // v(2] biggest     
  }
  //--------------------------------------------------------------
  // Generate LSQ general matrix transforming one set of 3-vectors into another
  // if the 3x3 matrix [U] transforms x -> y, with a possible error e
  // then e = [U] x - y
  //
  //  Return [U] & mean residual E
  //
  // Minimising E = 1/2N Sum(e.e)  ie least-squares gives
  //   [U] = [R] [S]^-1
  // where
  //   [S]ij  =  Sum(xi xj)    i, j = 1,3 for 3-vector
  //   [R]ij  =  Sum(yi xj)
  // On entry:
  //  x, y    corresponding vector lists
  //  AvoidPlanar  if true, generate a few orthogonal vectors to avoid the ambiguities
  //          introduced in all x vectors are coplanar (not included in residual E)
  //          but only if needed (only valid for orthogonal frames)
  DMat33 LsqTransform(const std::vector<DVect3>& x, const std::vector<DVect3>& y,
		      double& E, const bool& AvoidPlanar=false);
  //--------------------------------------------------------------
  // Get axis direction for rotation matrix R
  // Returns normalised axis direction
  clipper::Vec3<double> AxisDirection(const clipper::Mat33<double>& R);
  //--------------------------------------------------------------
}

#endif

