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
  // Analyse anomalous differences & adjust anomalous rejection criterion (in controls)
  // Returns mid-slopes of DelAnom normal probability plot for each dataset
  std::vector<float> AnalyseAnom(const hkl_unmerge_list& hkl_list,
				 const SDmodel& SDM,
				 all_controls& controls,
				 const bool& plot,
				 phaser_io::Output& output);

  int AccumulateDelanomNormProb(const SDmodel& SDM,
				const hkl_unmerge_list& hkl_list,
				const OutlierControl& outliercontrol,
				std::vector<NormalProbAnal>& normalprobanal);
}

#endif
