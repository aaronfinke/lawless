// \file cctbx_utils.hh

// Various utility routines for handling 3x3 matrices & vectors
//   using CCtbx routines
// For non-cctbx things see mtavec_utils
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
// cctbx::scitbx::mat3<double> (i,j)               SMat33
//
// Routines here are:-
//  1) converters                Set[type] (other_type)
//                               (not all combinations written yet!)
//  3) format reindex matrix     FormatReindex_as_hkl(Matrix)
//

#ifndef CCTBX_UTILS
#define CCTBX_UTILS

#include <cstring>

// Clipper
#include <clipper/clipper.h>
using clipper::Vec3;
using clipper::Mat33;
using clipper::RTop;
typedef clipper::Vec3<double> DVect3;
typedef clipper::Mat33<double> DMat33;

///#include "csymlib.h"    // CCP4 symmetry stuff

#include <cctbx/crystal/symmetry.h>
#include <cctbx/sgtbx/lattice_symmetry.h>
#include <cctbx/uctbx/fast_minimum_reduction.h>

#include <assert.h>
#define ASSERT assert

using namespace cctbx;

namespace MVutil {
  //--------------------------------------------------------------
  // vector(9) ->  matrix 3x3
  scitbx::mat3<double> SetSMat33(const std::vector<double>& Rmatrix);
  // Make scitbx::vec3 from clipper::Vec3 <double>
  scitbx::vec3<double> SetSvec3(const clipper::Vec3<double>& v3);
  //--------------------------------------------------------------
  // Construct rt_mx from vector(9) matrix 3x3
  // Default denominators are for crystal symmetry matrices
  sgtbx::rt_mx SetRtMx(scitbx::mat3<double>& cdmat,
		       const int r_den=1,
		       const int t_den=sgtbx::sg_t_den);
  //--------------------------------------------------------------
  sgtbx::rt_mx SetRtMx(const std::vector<double>& Rmatrix,
		       const int r_den=1,
		       const int t_den=sgtbx::sg_t_den);
  //--------------------------------------------------------------
  sgtbx::rt_mx SetRtMx(const clipper::RTop<double>& RT,
		       const int r_den=1,
		       const int t_den=sgtbx::sg_t_den);
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  // Format reciprocal space matrix (eg reindex or symmetry)
  // as h,k,l
  std::string FormatReindex_as_hkl(const scitbx::mat3<double>& op); //!< for cctbx
}

#endif

