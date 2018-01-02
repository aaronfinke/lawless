//  linearsdcorrection.hh
//

// sd' = SDfac * Sqrt(sd^2 + SdB* <I> + (SDadd * <I>)^2)
// Refine SDfac, SdB and SDadd, linear model
//
// DO NOT USE, doesn't work (yet)


#ifndef LINEAR_SDCORR_HEADER
#define LINEAR_SDCORR_HEADER

#include <vector>

#include "hkl_unmerge.hh"
#include "sdmodel.hh"
#include "intensitybin.hh"
#include "selectedobservations.hh"
#include "sdanalysis.hh"
#include "refinetargets.hh"
#include "linearlsq.hh"

namespace scala {
  // ---------------------------------------------------------
  // Refine SD correction model using linear model
  class LinearSDcorrection {
  public:
    LinearSDcorrection(SDmodel& SDM,
		       const hkl_unmerge_list& hkl_list,
		       const all_controls& controls,
		       IntensityBin& Irange,
		       const double& tolerance, const double& rtolerance,
		       const int&  max_cycles,
		       phaser_io::Output& output);
    void init(SDmodel& SDM,
	      const hkl_unmerge_list& hkl_list,
	      const all_controls& controls,
	      IntensityBin& Irange,
	      const double& tolerance, const double& rtolerance,
	      const int&  max_cycles,
	      phaser_io::Output& output);

  private:
    int Ndatasets;
    int npargroups;   // number of parameter sets for run, full/partial

    std::vector<LinearLSQ> linearlsq;
    LinearLSQ linearlsq1;

    void SumsforSDcorrection(const SDmodel& SDM,
			     const hkl_unmerge_list& hkl_list,
			     const bool& anomalous,
			     IntensityBin& irange);

    // Add in contribution to LSQ objects for this observation set
    void addContribution(SelectedObservations& selobs,
			 const double& Iav,
			 const SDmodel& SDM);
    //---------------------------------------------------------------
    //! returns residuals
    void UpdateParameters(SDmodel& SDM);
    
  };
}
#endif
