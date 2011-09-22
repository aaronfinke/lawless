// hkl_controls.cpp
//
#include <iostream>

#include "hkl_controls.hh"

namespace scala {
  //--------------------------------------------------------------
  file_select::file_select()
    //        *********
    // Default values if no input
  {
    init();
  }
  //--------------------------------------------------------------
  file_select::file_select(phaser_io::InputAll input, const int& NumFileSeries)
  {
    init();
    if (input.Valid()) {
      // Overall resolution limits
      range_sel = input.getRESO();
      // batch exclusions
      batchexclude = input.BatchExclusions();
      // batch inclusions from explicit RUN specification
      batchinclude = input.RunBatches();
      // Check that fileSeries specified on selection commands match
      // specified files. Fail here if not
      batchexclude.CheckSeries(NumFileSeries); 
      nullResolutionfraction = input.NullResolutionfraction();
      nullNegativeReject = input.NullNegativeReject();
    }
  }
  //--------------------------------------------------------------
  void file_select::init()
  {
    // Reject totals
    nrej_reso = 0;
    nrej_mflag = 0;
    inputscale = 1.0;
    // fraction of maximum resolution to use in test for blank batches
    nullResolutionfraction = -1.0;
    // Threshold on proportion of negative reflections
    nullNegativeReject = -1.0;
  }
  //--------------------------------------------------------------
  void file_select::set_reslimits(const ResoRange& resrange)
  {
    range_sel = resrange;
  }
  //--------------------------------------------------------------
  bool file_select::in_reslimits(const Rtype& s) const
  {
    if (range_sel.tbin(s) >= 0) return true;
    return false;
  }
  //--------------------------------------------------------------
  bool file_select::accept_batch
  (const int& batch_number, const int& fileSeries) const
  // Return false if batch_number is in rejection lists
  //  NB for reject options specifying batch numbers _after_
  //  any renumbering (fileSeriesList == 0), no test will be done
  //  (ie always returns true) unless fileSeries == 0
  {
    bool accept = true;
    if (!batchinclude.Null()) {
      // We have inclusions from specified RUNs
      //  accept if within selection
      accept = batchinclude.InSelection(batch_number, fileSeries);
    }
    if (batchexclude.InSelection(batch_number, fileSeries)) {
      accept = false;
    }
    return accept;
  }
}
