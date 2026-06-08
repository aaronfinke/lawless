// analysesd.hh

#ifndef ANALYSESD_HEADER
#define ANALYSESD_HEADER

#include "hkl_unmerge.hh"
#include "sdmodel.hh"
#include "Output.hh"
#include "score_datatypes.hh"
#include "intensitybin.hh"
#include "normalise.hh"
#include "selectedobservations.hh"
#include "normalprobanal.hh"

namespace scala
{
  // ------------------------------------------------------------
  // On entry:
  //  SDM        SD correction model
  //  hkl_list   reflection data
  //  controls   
  //  NormRes    normalisation
  //  firstAnalysis
  //     = 0  two analyses expected, after rough scaling and main scaling
  //     = +1    same, after first analysis
  //     = -1 only one analysis expected ie first & last
  void AnalyseSD(SDmodel& SDM, hkl_unmerge_list& hkl_list,
		 all_controls& controls,
		 const Normalise& NormRes,
		 const int& firstAnalysis,
		 phaser_io::Output& output);

  void AnalyseNormalProbability(SDmodel& SDM, hkl_unmerge_list& hkl_list,
				all_controls& controls, const bool& plot,
				phaser_io::Output& output);
  // ------------------------------------------------------------
  class SDMdataNumbers {  // mostly just a struct to return data
  public:
    SDMdataNumbers(){}
    SDMdataNumbers(const int& Nsets) : nsets(Nsets) {
      enoughdata = true;
      nfnp.resize(nsets);
    }

    int Number(const int& iset) const {  // total number for this set
      return nfnp.at(iset).first + nfnp.at(iset).second;
    }

    bool enoughdata;
    int nsets;
    std::vector<std::pair<int,int> > nfnp;  // number of fulls/partials for each set

  };
  // ------------------------------------------------------------
  SDMdataNumbers UpdateSDMfromNPlot(SDmodel& SDM, const hkl_unmerge_list& hkl_list,
				    const all_controls& controls, const bool& fixup,
				    const bool& update,
				    phaser_io::Output& output);
  // ------------------------------------------------------------
  int AccumulateNormProb(const SDmodel& SDM, const hkl_unmerge_list& hkl_list,
			  const bool& Anomalous,
			 std::vector<NormalProbAnal>& normalprobanal);
  //---------------------------------------------------------------
  // for debuggery
  void MakeSDplot(SDmodel& SDM,
		  hkl_unmerge_list& hkl_list,
		  const all_controls& controls,
		  const Normalise& NormRes,
		  phaser_io::Output& output);
}

#endif
