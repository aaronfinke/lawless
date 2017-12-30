// analyseanom.hh

#ifndef ANALYSEANOM_HEADER
#define ANALYSEANOM_HEADER

#include "hkl_unmerge.hh"
#include "sdmodel.hh"
#include "Output.hh"
#include "controls.hh"
#include "normalprobanal.hh"

namespace scala
{
  class AnalyseAnom {

  public:
    AnalyseAnom(){}

  // Analyse anomalous differences & adjust anomalous rejection criterion (in controls)
  // Returns mid-slopes of DelAnom normal probability plot for each dataset
    AnalyseAnom(const hkl_unmerge_list& hkl_list,
		const SDmodel& SDM,
		all_controls& controls,
		const ResoRange& ResRange,
		const bool& plot,
		phaser_io::Output& output);


    // Returns mid-slopes of DelAnom normal probability plot for each dataset
    std::vector<double> Slopes() const {return slopes;}

    // for each dataset for each resolution bin
    std::vector<std::vector<MeanSD> > RmsDelAnom() const {return rmsdelanom;}


  private:
    int ndatasets;
    std::vector<double> slopes;      // for each dataset
    int nresbin;                    // number of resolution bins
    std::vector<std::vector<MeanSD> > rmsdelanom; // for each dataset for each resolution bin

    int AccumulateDelanomNormProb(const SDmodel& SDM,
				  const hkl_unmerge_list& hkl_list,
				  const OutlierControl& outliercontrol,
				  const ResoRange& ResRange,
				  std::vector<NormalProbAnal>& normalprobanal);
  };
}

#endif
