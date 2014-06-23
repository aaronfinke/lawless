// observationstatuscontrol.cpp
//
//

#include  "observationstatuscontrol.hh"


namespace scala {

// ------------------------------------------------------------
int SetOverlapFlags(const bool& Exclude, hkl_unmerge_list& hkl_list)
//  Set or unset multi-lattice overlap status flags for all observations
//  Other flags are left unaltered
// if Exclude == true, then set flag to exclude all multiple observations
// if Exclude == false, then set flag to include multiple observations
  {
    if (!hkl_list.MultiLattice()) {return 0;}  // ignore unless multilattice

    reflection this_refl;
    observation this_obs;
    hkl_list.rewind();
    int nrej = 0;

    // loop all reflections unconditionally
    for (int jref=0;jref<hkl_list.num_reflections();++jref) {
      this_refl = hkl_list.get_reflection(jref);
      // loop all observations, ignoring accept flag
      for (int lobs=0;lobs<this_refl.num_observations();++lobs) {
	this_obs = this_refl.get_observation(lobs);
	ObservationStatus status = this_obs.ObsStatus();
	if (Exclude && !this_obs.IsSingleton()) {
	  nrej++;
	  // reject multiples if Exclude 
	  status.SetRejectOverlap();
	} else { // else keep
	  status.UnsetRejectOverlap();  // keep
	}
	this_obs.UpdateStatus(status);
	this_refl.replace_observation(this_obs);
      }
      hkl_list.replace_reflection(this_refl);
    }
    hkl_list.SetExcludeOverlaps(Exclude);  // store record
    return nrej;
  }
// ------------------------------------------------------------
  void SetRunsToUse(const std::vector<bool>& userun,
		    hkl_unmerge_list& hkl_list)
  // Set to use runs flagged in userun
  {
    ASSERT (int(userun.size()) == hkl_list.num_runs());

    reflection this_refl;
    observation this_obs;
    hkl_list.rewind();

    // loop all reflections unconditionally
    for (int jref=0;jref<hkl_list.num_reflections();++jref) {
      this_refl = hkl_list.get_reflection(jref);
      // loop all observations, ignoring accept flag
      for (int lobs=0;lobs<this_refl.num_observations();++lobs) {
	this_obs = this_refl.get_observation(lobs);
	ObservationStatus status = this_obs.ObsStatus();
	if (userun[this_obs.run()]) {
	  status.UnsetRejectRun();   // Use
	} else {
	  status.SetRejectRun();     // Don't use
	}
	this_obs.UpdateStatus(status);
	this_refl.replace_observation(this_obs);
      }
      hkl_list.replace_reflection(this_refl);
    }
  }
// ------------------------------------------------------------
  void SetBatchesToUse(const std::vector<bool>& usebatch,
		    hkl_unmerge_list& hkl_list)
  // Set to reject batches flagged in usebatch
  // Updates:
  //   hkl_list.batches  flag rejected batches
  //   hkl_list.runlist  remove rejected batches from run
  //   all observations, eliminate all observations from rejected batches
  {
    ASSERT (int(usebatch.size()) == hkl_list.num_batches());

    bool somerejected = false; // are there any rejections?
    for (size_t ib=0; ib<usebatch.size(); ib++) { 
      if (!usebatch[ib]) {
	somerejected = true;
	break;
      }
    }
    if (!somerejected) {return;} // nothing to do

    std::vector<int> rejectedbatches;
    for (size_t ib=0; ib<usebatch.size(); ib++) { 
      if (!usebatch[ib]) {
	// reject this batch from batch list and run
	hkl_list.RejectBatchSerial(ib);
	int batchnum = hkl_list.batch(ib).num();
	rejectedbatches.push_back(batchnum);
      }
    }

    if (rejectedbatches.size() > 0) {
      int nrejb = hkl_list.EliminateBatches(rejectedbatches);
      hkl_list.prepare();
      hkl_list.sum_partials();
    }

    /*
    reflection this_refl;
    observation this_obs;
    hkl_list.rewind();

    // loop all reflections unconditionally
    for (int jref=0;jref<hkl_list.num_reflections();++jref) {
      this_refl = hkl_list.get_reflection(jref);
      // loop all observations, ignoring accept flag
      for (int lobs=0;lobs<this_refl.num_observations();++lobs) {
	this_obs = this_refl.get_observation(lobs);
	ObservationStatus status = this_obs.ObsStatus();
	int batch = this_obs.Batch(); // central batch
	int jbatch = hkl_list.batch_serial(batch);
	if (!usebatch[jbatch]) {
	  if (!this_obs.IsFull()) {
	    // sanity check
	    Message::message(Message_fatal
	       ("SetBatchesToUse: cannot reject individual batches for Partials"));
	  }
	  status.SetRejectBatch();     // Don't use
	  this_obs.UpdateStatus(status);
	  this_refl.replace_observation(this_obs);
	}
      }
      hkl_list.replace_reflection(this_refl);
    }
    */
  }
  // ------------------------------------------------------------
  void ClearObsStatus(hkl_unmerge_list& hkl_list, const bool& allflags)
  // Clear all observation status flags back to the ObservationFlag setting
  // ie clear outlier & Emax status flags
  // the ObservationFlag setting is left unaltered
  // Set reflection status to Accept
  // Don't change RejectOverlap, run & batch rejections flag unless allflags true
  {
    reflection this_refl;
    observation this_obs;
    hkl_list.rewind();

    // loop all reflections unconditionally
    for (int jref=0;jref<hkl_list.num_reflections();++jref) {
      this_refl = hkl_list.get_reflection(jref);
      // loop all observations, ignoring accept flag
      for (int lobs=0;lobs<this_refl.num_observations();++lobs) {
	this_obs = this_refl.get_observation(lobs);
	ObservationStatus status = this_obs.ObsStatus();
	bool overlapstatus = status.TestRejectOverlap();
	if (allflags) {
	  status.ResetStatusAll();
	} else {
	  // leaves ObservationFlag, resolution, overlap, run & batch rejections
	  status.ResetStatus();
	}
	this_obs.UpdateStatus(status);
	this_refl.replace_observation(this_obs);
      }
      this_refl.SetStatus(0);  // accept reflection
      hkl_list.replace_reflection(this_refl);
    }
  }
} // namespace scala
