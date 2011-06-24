//(c) 2000-2007 Cambridge University Technical Services Ltd
//All rights reserved
#ifndef __RefineBase2Class__
#define __RefineBase2Class__
#include "phaser_types.h"
#include "Output.h"
#include "Errors.h"
#include "fmat.h"
#include "RefineBase.h"
#include "vec.h"

namespace phaser {

class RefineBase2 : public RefineBase 
{
  public: 
    RefineBase2() {};
    virtual ~RefineBase2() {};

    double finiteDiffGradient(TNT::Vector<double>&,double);
    double finiteGDiffHessian(TNT::Fortran_Matrix<double>&,double);
    double finiteGDiffDiagHessian(TNT::Fortran_Matrix<double>&,double);
    double finiteFDiffHessian(TNT::Fortran_Matrix<double>&,double);
    double finiteFDiffDiagHessian(TNT::Fortran_Matrix<double>&,double);
};

} //phaser
#endif
