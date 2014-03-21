// statistics.h
//
// Print all merging statistics etc
//

#ifndef STATISTICS_HEADER
#define STATISTICS_HEADER

#include "scalemodel.hh"
#include "hkl_unmerge.hh"
#include "Output.hh"
#include "sdmodel.hh"
#include "controls.hh"
#include "Output.hh"
#include "normalise.hh"
#include "anomdistribution.hh"
#include "summarystatistics.hh"
#include "referencelist.hh"
 
namespace scala {
  //
  // Statistics for within a dataset datasetIndex
  //
  //  On entry:
  //   AllScales    scale model
  //   hkl_list     list with scales applied & outliers rejected
  //   SDM          Sd correction model
  //   controls     all controls
  //   datasetIndex dataset index
  //   NormRes      normalisation object, over all data (no run/batch dependence)
  //   anomProbSlope slope of anomalous normal probability plot
  //   hklreflist   reference data for analysis, if present
  //   output
  //  Returns summary statistics for this dataset 
  SummaryStatistics Statistics(const ScaleModel& AllScales,
		  const hkl_unmerge_list& hkl_list,
		  SDmodel& SDM,
		  const all_controls& controls, const int& datasetIndex, 
		  const ResoRange& ResRange,
		  const Normalise& NormRes,
		  const AnomDistribution& anomDistribution,
		  const float& anomProbSlope,
		  const ReferenceList& hklreflist,
		  phaser_io::Output& output);
}

#endif
