// observationstatuscontrol.hh



#ifndef OBSERVATIONSTATUSCONTROL_HEADER
#define OBSERVATIONSTATUSCONTROL_HEADER

#include "hkl_unmerge.hh"

namespace scala {

  //  Set or unset multi-lattice overlap status flags for all observations
  //  Other flags are left unaltered
  // if Exclude == true, then set flag to exclude all multiple observations
  // if Exclude == false, then set flag to include multiple observations
  // return number of overlapped rejections
  int SetOverlapFlags(const bool& Exclude, hkl_unmerge_list& hkl_list);

  // Set to use runs flagged in userun
  void SetRunsToUse(const std::vector<bool>& userun, hkl_unmerge_list& hkl_list);

  // Clear all observation status flags back to the ObservationFlag setting
  // ie clear outlier & Emax status flags
  // the ObservationFlag setting is left unaltered
  // Don't change RejectOverlap flag
  void ClearObsStatus(hkl_unmerge_list& hkl_list);
}

#endif
