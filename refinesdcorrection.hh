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
#include "refinetargets.hh"

namespace scala {
  //---------------------------------------------------------------
  // Refine SD correction model using BGFS minimiser
  SDanalysis RefineSDcorrectionFactors(SDmodel& SDM,
				       const hkl_unmerge_list& hkl_list,
				       const all_controls& controls,
				       IntensityBin& Irange,
				       const double& tolerance, const double& rtolerance,
				       const int&  max_cycles,
				       phaser_io::Output& output);
  //---------------------------------------------------------------
  SDanalysis SumsforSDcorrection(const SDmodel& SDM,
				 const hkl_unmerge_list& hkl_list,
				 const bool& anomalous,
				 IntensityBin& irange);

  //---------------------------------------------------------------
  //! returns residuals
  TargetResiduals UpdateParameters(SDmodel& SDM, const SDanalysis& sdanal,
				   const double& tolerance, const double& damp,
				   const bool& Noupdate);

}
#endif
