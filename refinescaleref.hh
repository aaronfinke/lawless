//  refinescaleref.hh
//

//  Scaling against reference set, using BFGS minimiser

#ifndef REFINESCALEREF_HEADER
#define REFINESCALEREF_HEADER

#include <vector>

#include "phaser_types.hh"
#include "ProtocolBase.h"
#include "RefineBase2.h"
#include "Minimizer.h"

#include "ProtocolScale.hh"
#include "hkl_unmerge.hh"
#include "referencelist.hh"
#include "scalemodel.hh"
#include "sdmodel.hh"

using namespace phaser;

namespace scala {
class RefineScaleRef : public phaser::RefineBase2
{
public:
  RefineScaleRef(){}
  RefineScaleRef(const hkl_unmerge_list& Hkl_list,
		 const ReferenceList& Hklreflist,
		 ScaleModel& Scalemodel,
		 const SDmodel& SDM,
		 const int& Nprocs);

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
  std::vector<int> Nobservations() const {return nrefpar;}
  int NobservationsAll() const {return nobs;}

private:
  ScaleModel* scalemodel;   // all the scales
  const SDmodel* sdmodel;
  int npar;   // number of parameters
  std::vector<double> params;    // the parameters
  std::vector<int> nrefpar; // number of contributions to each parameter
  
  const hkl_unmerge_list* hkl_list;  // reflection data
  const ReferenceList* hklreflist;   // reference list


  bool gradientOK;  // true if there is an up-to-date gradient

  TNT::Vector<floatType> gradient;
  floatType target;
  int nobs;

  int nprocs;
  // weight type:
  //  > 0  1/(var(Ihl) + var(Iref))
  //  < 0  1/var(Ihl)
  //  else = 1.0
  int weighttype;

  void TargetGradientHessian(bool DoGradient,
			     bool DoHessian,
			     TNT::Fortran_Matrix<floatType>& H);

  double vweight(const double& sd, const IsigI& isigiref) const;
};
}
#endif
