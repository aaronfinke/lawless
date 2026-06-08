//
// scalerefinefh.hh
//
// Main scaling routine , traditional filtered Newton method (Fox-Holmes)
//
//

#ifndef SCALEREFINEFH_HEADER
#define SCALEREFINEFH_HEADER

#include "hkl_unmerge.hh"
#include "scalemodel.hh"
#include "controls.hh"
#include "Output.hh"

namespace scala {

// ---------------------------------------------------------
  void ScaleRefineFH(hkl_unmerge_list& hkl_list, ScaleModel& AllScales,
		     const all_controls& controls, const int& Ncycles,
		     phaser_io::Output& output);
}


#endif

