// applyscales.h
//
// Testing for reading refecltions & getting scales
//

#ifndef APPLYSCALES_HEADER
#define APPLYSCALES_HEADER

#include "scalemodel.hh"
#include "hkl_unmerge.hh"

namespace scala {
  //! Apply scales to all data, return mean(I) in resolution limits
 double ApplyScales(const ScaleModel& AllScales, hkl_unmerge_list& hkl_list);
}

#endif
