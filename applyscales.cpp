// applyscales.cpp
//
// calculate & store all scales
//

#include "scalemodel.hh"
#include "hkl_unmerge.hh"
#include "scala_util.hh"

namespace scala {
  double ApplyScales(const ScaleModel& AllScales, hkl_unmerge_list& hkl_list)
  //! Apply scales to all data, return mean(I) within resolution limits
  {
    reflection this_refl;
    observation this_obs;
    hkl_list.rewind();
    ResoRange resrange = hkl_list.ResLimRange();  // resolution limits
    MeanSD meanI;  // mean scaled I within resolution limits

    // * * * * Loop reflections
    // loop all reflections unconditionally
    for (int jref=0;jref<hkl_list.num_reflections();++jref) {
      this_refl = hkl_list.get_reflection(jref);
      Rtype invresolsq = this_refl.invresolsq();
      bool inrange = (resrange.tbin(invresolsq) >= 0); // true if in reso limits

      //  Loop all observations
      for (int i=0;i<this_refl.num_observations();++i) {
	this_obs = this_refl.get_observation(i);
	AllScales.ScaleObs(this_obs, invresolsq);  // apply scale to observation
	this_refl.replace_observation(this_obs);
	if (inrange) {
	  meanI.Add(this_obs.kI());
	}
      }
      hkl_list.replace_reflection(this_refl); // store updated reflection
    }
    return meanI.Mean();  // mean scaled I
  }
}
