//  refinescale.hh
//

//  Main Fox-Holmes scaling with all parameters, using BFGS minimiser

#ifndef REFINESCALE_HEADER
#define REFINESCALE_HEADER

#include <vector>

#include "phaser_types.hh"
#include "ProtocolBase.h"
#include "RefineBase2.h"
#include "Minimizer.h"

#include "ProtocolScale.hh"
#include "hkl_unmerge.hh"
#include "scalemodel.hh"

using namespace phaser;

namespace scala {
class RefineScale : public phaser::RefineBase2
{
public:
  RefineScale(){}
  RefineScale(const hkl_unmerge_list& Hkl_list,
			   ScaleModel& Scalemodel);

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

private:
  ScaleModel* scalemodel;   // all the scales
  int npar;   // number of parameters
  std::vector<double> params;    // the parameters
  std::vector<int> nrefpar; // number of contributions to each parameter
  
  const hkl_unmerge_list* hkl_list;  // reflection data

  bool gradientOK;  // true if there is an up-to-date gradient

  TNT::Vector<floatType> gradient;
  floatType target;

  void TargetGradientHessian(bool DoGradient,
			     bool DoHessian,
			     TNT::Fortran_Matrix<floatType>& H);

};
}
#endif
