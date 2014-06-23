// rejectbatches.cpp
//
// Reject batches if their scales are way out of line with the others
//

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;

#include "rejectbatches.hh"
#include "statistics.hh"
#include "scala_util.hh"
#include "string_util.hh"
#include "observationstatuscontrol.hh"

using phaser_io::LOGFILE;
using phaser_io::LXML;

RejectBatches::RejectBatches(hkl_unmerge_list& hkl_list,
			     ScaleModel& AllScales,
			     const all_controls& controls,
			     phaser_io::Output& output)
{
  // Sanity checks: only use this for batch scaling
  if (!AllScales.isAllBatch()) {
      Message::message(Message_fatal
       ("Cannot use REJECT BATCH unless BATCH scaling is used"));
  }

  float batchrejectfactor = controls.outlierScale.Reject(ALL).batchrejectfactor;
  if (batchrejectfactor < 0.00001) {return;}

  output.logTab(0, LOGFILE,
		"\nTesting Batch scales for unlikely values\n");
  output.logTabPrintf(0, LOGFILE,
    "Reject batches if scale is negative or greater than %7.3f times the median value\n",
		      batchrejectfactor);

  // Get all batch scale factors
  std::vector<float> scale0batch;
  std::vector<float> bfacbatch;
  std::vector<Batch> batches = hkl_list.Batches();
  BatchScales0(batches, -1, AllScales,
	       scale0batch, bfacbatch);
  std::vector<float> scale0batch_sort = scale0batch;
  std::sort(scale0batch_sort.begin(), scale0batch_sort.end());
  float medianscale = Median(scale0batch_sort, -1);

  float threshold = medianscale * batchrejectfactor;
  output.logTabPrintf(0, LOGFILE, "\nMedian scale %8.3f, reject threshold %8.3f\n",
		      medianscale, threshold);

  int nbatches = batches.size();
  std::vector<bool> usebatch(nbatches, true);
  std::vector<int> negativebatches;
  int nrej = 0;
  for (size_t ib=0; ib<scale0batch.size(); ib++) { 
    if (scale0batch[ib] < 0.000001) {
      //      std::cout << "Negative scale, batch serial "<<ib<<"\n"; //^
      negativebatches.push_back(ib);
      usebatch[ib] = false;
    }
    if (scale0batch[ib] > threshold) { // scale too large
      usebatch[ib] = false;
      nrej++;
    }
  }

  std::string s;
  if (negativebatches.size() > 0) {
    output.logTab(0, LOGFILE,
		  "Batches with negative scales:");
    s = "";
    for (size_t ib=0; ib<negativebatches.size(); ib++) { 
      s += " "+StringUtil::itos(batches[negativebatches[ib]].num(), 7);
    }
    output.logTab(0, LOGFILE, StringUtil::WrapLine(s, 80, 0));
  }
  if (nrej > 0) {
    output.logTab(0, LOGFILE, 
		  "\nBatches reject for scales too large (Batch: scale)\n");
    s = "";
    for (size_t ib=0; ib<usebatch.size(); ib++) { 
      if (!usebatch[ib]) {
	int batchnumber = batches[ib].num();
	float scale = scale0batch[ib];
	s += " "+StringUtil::itos(batchnumber, 7) +
	  ": "+StringUtil::ftos(scale,7,3);
      }
    }
    output.logTab(0, LOGFILE, StringUtil::WrapLine(s, 80, 0));
  }

  // set flag to mark observations in rejected batches
  SetBatchesToUse(usebatch, hkl_list);
  
  // Update scale model
  std::vector<int> batchnumbers(nbatches);
  for (int ib=0;ib<nbatches;++ib) {
    batchnumbers[ib] = batches[ib].num();
  }
  AllScales.setBatchReject(usebatch, batchnumbers);

}
