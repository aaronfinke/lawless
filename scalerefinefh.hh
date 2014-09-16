//
// scalerefinefh.hh
//
// Main scaling routine , traditional filtered Newton method (Fox-Holmes)
//
//

#ifndef SCALEREFINEFH_HEADER
#define SCALEREFINEFH_HEADER

#include "hkl_unmerge.hh"
#include "scalemodel.hh"
#include "controls.hh"
#include "Output.hh"
#include "fmat.h"


namespace scala {
// ---------------------------------------------------------
class Hessian
{
public:
  Hessian():Np(0){}
  Hessian(const TNT::Fortran_Matrix<floatType>& H);

  //   Note: Hessian H is addressed by fortran-style indexing from 1
  inline floatType operator()(int i, int j) {return hessian(i,j);}
  inline const floatType operator()(int i, int j) const {return hessian(i,j);}

  TNT::Vector<floatType> Ubs();

  // Filtered inverse rescaled Hessian
  //   MinFiltered  minimum number of eigenvalues filtered
  TNT::Fortran_Matrix<floatType> FilteredInverse(const int& MinFiltered);

  std::string format() const;

private:
  TNT::Fortran_Matrix<floatType> hessian;
  int Np;  // dimension
  TNT::Vector<floatType> ubs;  // scaling vector
};

// ---------------------------------------------------------
  void ScaleRefineFH(hkl_unmerge_list& hkl_list, ScaleModel& AllScales,
		     const all_controls& controls, const int& Ncycles,
		     phaser_io::Output& output);
}


#endif

