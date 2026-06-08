// selectsdcorrreflections.hh
//
// Mark reflections with accept/reject flags to be used
// for SD correction optimisation
//
// 1) no singletons
// 2) Roughly evenly distributed of |E^2| in lower intensity bins
// 3) all higher intensities E^2 > E2min
//  

#ifndef SELECTSDCORRREFLECTIONS_HEADER
#define SELECTSDCORRREFLECTIONS_HEADER

#include "hkl_unmerge.hh"
#include "controls.hh"
#include "normalise.hh"

namespace scala {
  // ------------------------------------------------------------
  // On entry:
  //   hkl_list    reflection list, scales applied if needed
  //   Nbintarget  target minimum number of reflections / intensity bin
  //
  // On exit:
  //   hkl_list    reflection list, reflection accept flags updated
  //
  // returns number of reflections accepted
  int SelectSDcorrReflections(hkl_unmerge_list& hkl_list,
			      const all_controls& controls,
                              const int& Nbintarget,
			      const Normalise& NormRes);

}
#endif
