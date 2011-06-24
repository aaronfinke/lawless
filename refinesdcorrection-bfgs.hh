//  refinesdcorrection.hh
//

// sd' = SDfac * Sqrt(sd^2 + (SDadd * <I>)^2)
// Refine SDfac and SDadd

#ifndef REFINE_SDCORR_HEADER
#define REFINE_SDCORR_HEADER

#include <vector>

#include "hkl_unmerge.hh"
#include "sdmodel.hh"
#include "intensitybin.hh"
#include "selectedobservations.hh"
#include "sdanalysis.hh"

#include "phaser_types.hh"
#include "ProtocolBase.h"
#include "RefineBase2.h"

using namespace phaser;
namespace scala {
  //---------------------------------------------------------------
  // Refine SD correction model using BGFS minimiser
  void RefineSDcorrectionFactors(SDmodel& SDM,
		      const hkl_unmerge_list& hkl_list,
		      const all_controls& controls,
		      IntensityBin& Irange,
		      const double& tolerance,
		      const int&  max_cycles,
		      phaser_io::Output& output);

  //---------------------------------------------------------------
class RefineSDCorrection : public phaser::RefineBase2
{
public:
  RefineSDCorrection() : DEBUG(false){}
  RefineSDCorrection(SDmodel& sdm,
		     const hkl_unmerge_list& Hkl_list,
		     IntensityBin& Irange,
		     const bool& Anomalous);
  void init(SDmodel& sdm,
	    const hkl_unmerge_list& Hkl_list,
	    IntensityBin& Irange,
	    const bool& Anomalous);

  floatType    targetFn();
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
  SDmodel* SDM;                      // SD model (may change)
  const hkl_unmerge_list* hkl_list;  // reflection data (const, doesn't change)
  IntensityBin* irange;
  bool anomalous;
  std::vector<bool> combine;  // combine fulls & partials, for each run
  mutable SDanalysis sdanal;

  int npar; // number of parameters

  bool gradientOK;  // true if there is an up-to-date gradient
  TNT::Vector<floatType> gradient;
  floatType target;
  bool targetOK;

  bool DEBUG;

  void TargetGradientHessian(bool DoGradient,
			     bool DoHessian,
			     TNT::Fortran_Matrix<floatType>& H);

};
}
#endif
