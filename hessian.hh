//
// hessian.hh
//

#ifndef HESSIAN_HEADER
#define HESSIAN_HEADER

#include "Minimizer.h"
#include "fmat.h"


namespace scala {
// ---------------------------------------------------------
class Hessian
{
public:
  Hessian():Np(0){}
  Hessian(const TNT::Fortran_Matrix<floatType>& H);

  //   Note: Hessian H is addressed by fortran-style indexing from 1
  //  inline floatType operator()(int i, int j) {return hessian(i,j);}
  //  inline const floatType operator()(int i, int j) const {return hessian(i,j);}

  //TNT::Vector<floatType> Ubs();

  // Filtered inverse rescaled Hessian
  //   MinFiltered  minimum number of eigenvalues filtered
  TNT::Fortran_Matrix<floatType> FilteredInverse(const int& MinFiltered);

  std::string format() const;

private:
  TNT::Fortran_Matrix<floatType> hessian;
  int Np;  // dimension
  std::vector<floatType> ubs;  // scaling vector

  void Ubs();

};

}
#endif

