//
// scalerefine.hh
//
// Main scaling routine
//

#ifndef SCALEREFINE_HEADER
#define SCALEREFINE_HEADER

#include "hkl_unmerge.hh"
#include "scalemodel.hh"
#include "controls.hh"
#include "Output.hh"
#include "sdmodel.hh"

namespace scala {
// ---------------------------------------------------------
  void ScaleRefine(hkl_unmerge_list& hkl_list, ScaleModel& AllScales,
		   const SDmodel& SDM,
		   const all_controls& controls, const int& Ncycles,
		   const bool& print, phaser_io::Output& output);
}

#endif

