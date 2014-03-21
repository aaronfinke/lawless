// selectscalingreflections.hh
//
// Mark reflections to be used for scaling with accept/reject flags
//
// Possible criteria:
//
//  1) accept if all (accepted) observations have I/sd'(I) > IovSDmin
//     sd'(I) after correction by SDmodel
//
//  2) accept if |E^2| > E2min after normalisation
//  

#ifndef SELECTSCALINGREFLECTIONS_HEADER
#define SELECTSCALINGREFLECTIONS_HEADER

#include "hkl_unmerge.hh"
#include "sdmodel.hh"
#include "scalemodel.hh"

namespace scala {
  // On entry:
  //   hkl_list    reflection list, scales applied if needed
  //   SDM         sd correction model
  //   AllScales   scale model just used to set Phi bins for I/sd cutoff
  //   IovSDmin    minimum value for  <I>/sd'(<I>), == 0 no test
  //                < 0 negative value from default, to be reset here 
  //   E2min       minimum |E^2|, <= 0 no test
  //   E2max       maximum |E^2|, <= 0 no test
  //
  // On exit:
  //   hkl_list    reflection list, reflection accept flags updated
  //
  // returns number of reflections rejected, and nskip

  std::pair<int,int> SelectScalingReflections(hkl_unmerge_list& hkl_list,
					      const SDmodel& SDM,
					      const ScaleModel& AllScales,
					      float& IovSDmin,
					      const float& E2min, const float& E2max);
}

#endif
