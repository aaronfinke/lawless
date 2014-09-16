//  refinereferencescale.hh
//
// 
//  Refine scale of reference intensities to test set, using BFGS minimiser

#ifndef REFINEREFERENCESCALE_HEADER
#define REFINEREFERENCESCALE_HEADER

#include <assert.h>
#define ASSERT assert

#include "mergedlist.hh"
#include "referencescalemodel.hh"
#include "timer.hh"
#include "refinetargets.hh"

#include "phaser_types.hh"
#include "ProtocolBase.h"
#include "RefineBase2.h"
#include "Minimizer.h"

#include "ProtocolScale.hh"


namespace scala {
  // ---------------------------------------------------------
  class RefineReferenceScale  : public phaser::RefineBase2 {
  public:
    RefineReferenceScale(){}

    RefineReferenceScale
    (clipper::HKL_data<clipper::data32::I_sigI>& Isigi,
     const hkl_merge& Hklmergelist,
     ReferenceScaleModel& Referencescalemodel,
     const double& Meanintensity,
     const int& Nprocs);

    // set true to use quadratic residual, else ln(cosh(d))
    void setQuadratic(const bool& Quadratic);
    // Set target type
    void setTargetType(const RefineTargets::REFINETARGETTYPES& targetType);
    // Return target type
    RefineTargets::REFINETARGETTYPES targetType() const {return targettype;}

    floatType    targetFn();  //this is where main body goes
    floatType    gradientFn(TNT::Vector<floatType>&);
    floatType    hessianFn(TNT::Fortran_Matrix<floatType>&,bool&);

    void         applyShift(TNT::Vector<floatType>&);
    void         logCurrent(outStream,Output&);
    bool1D       getRefineMask(protocolPtr);

    std::string  whatAmI(int&);
    std::vector<bounds>    getLowerBounds();
    std::vector<bounds>    getUpperBounds();
    TNT::Vector<floatType> getRefinePars();
    TNT::Vector<floatType> getLargeShifts();
    int Nobservations() const {return nref;}

  private:
    const clipper::HKL_data<clipper::data32::I_sigI>* isigi;  // reflection data
    const hkl_merge*     hklmergelist;
    ReferenceScaleModel* referencescalemodel;

    int npar;
    int nprocs;
    int nref;

    double meanintensity;

    std::vector<double> params;    // the parameters
    TNT::Vector<floatType> gradient;
    floatType target;
    bool gradientOK;

    // default to use ln(cosh(d)) residual (LNCOSH), else QUADRATIC
    RefineTargets::REFINETARGETTYPES targettype;

    void TargetGradientHessian(bool DoGradient,
			       bool DoHessian,
			       TNT::Fortran_Matrix<floatType>& H);
    void TargetGradientHessianU(bool DoGradient,
			       bool DoHessian,
			       TNT::Fortran_Matrix<floatType>& H);
    void TargetGradientHessianM(bool DoGradient,
			       bool DoHessian,
			       TNT::Fortran_Matrix<floatType>& H);


  };
}
#endif
