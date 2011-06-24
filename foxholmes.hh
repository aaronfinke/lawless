//  FoxHolmes.hh
//

//  Simple-minded Fox-Holmes scaling, using BFGS minimiser

#ifndef FOXHOLMES_HEADER
#define FOXHOLMES_HEADER

#include <vector>

#include "phaser_types.hh"
#include "ProtocolBase.h"
#include "RefineBase2.h"
#include "Minimizer.h"

#include "initialscales.hh"
#include "ProtocolScale.hh"

using namespace phaser;

namespace scala { 
class FoxHolmes : public phaser::RefineBase2
{
public:
  FoxHolmes(){}
  FoxHolmes(const scala::InitialData& Data);

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
  

private:
  int npar;
  std::vector<double> scales;
  const scala::InitialData* data;
  bool gradientOK;  // true if there is an up-to-date gradient
  TNT::Vector<floatType> gradient;
  floatType target;

  int MeanI(const std::vector<DPair>& y,
	    double& mnI, double& sumwg2) const;
  void TargetGradientHessian(bool DoGradient,
			     bool DoHessian,
			     TNT::Fortran_Matrix<floatType>& H);

};
}
#endif
