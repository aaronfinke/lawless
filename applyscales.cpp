// applyscales.cpp
//
// calculate & store all scales
//

#include "scalemodel.hh"
#include "hkl_unmerge.hh"

namespace scala {
  void ApplyScales(const ScaleModel& AllScales, hkl_unmerge_list& hkl_list)
  {
    reflection this_refl;
    observation this_obs;
    hkl_list.rewind();

    // * * * * Loop reflections
    // loop all reflections unconditionally
    for (int jref=0;jref<hkl_list.num_reflections();++jref) {
      this_refl = hkl_list.get_reflection(jref);
      Rtype invresolsq = this_refl.invresolsq();
      //  Loop all observations
      for (int i=0;i<this_refl.num_observations();++i) {
	this_obs = this_refl.get_observation(i);
	AllScales.ScaleObs(this_obs, invresolsq);  // apply scale to observation
	this_refl.replace_observation(this_obs);
      }
      hkl_list.replace_reflection(this_refl); // store updated reflection
    }
  }
}
