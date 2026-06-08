// restore.cpp

#include "restore.hh"

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;

namespace scala {
  //--------------------------------------------------------------
  void ScaleRestoreError::error(const std::string& message)
  {
    clipper::Message::message(Message_fatal
                              ("RESTORE error: incompatible save file: "+message));
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  void RunsFromSavefile::init(Fileread& FR, const std::vector<Run>& runlist)
  // Read run block from save file, get run definitions for runs which are wanted
  {
    FR.ReadTag("Nruns");
    svnruns = FR.Int();  // number of runs in save file
    nruns = runlist.size();  // number of runs in runlist = number to be used
    // for each run in save file, set to run index in runlist, or -1 if not used
    runsfromsavefile.assign(svnruns, -1);
    nrfound = 0;

    // Loop runs in save file
    for (int ir=0;ir<svnruns;++ir) {
      FR.ReadTag("RunNumber");
      int runnum = FR.Int();
      runnum = runnum; // dummy
      FR.ReadTag("Run");
      if (FR.GetTag() != "V1") {  // version check
        clipper::Message::message(Message_fatal
                                  ("RESTORE incompatible run version"));
      }
      FR.Skip(); // skip "{"
      FR.ReadTag("Batch_number_list");
      int nbat = FR.Int();
      std::vector<int> batchnumbers = FR.IntVec(nbat);
      if (!FR.CheckEnd()) {
        clipper::Message::message(Message_warn
                                  ("RESTORE unexpected tag "+FR.Tag()));
      }
      // Does this run match any in runlist?
      // return index in runlist, -1 if not found
      int irun = RunNotFound(runlist, batchnumbers);
      if (irun >= 0) {
        // Build list of runs from save file which should be kept
        runsfromsavefile[ir] = irun;
        nrfound++;
      }
    } // end loop runs in save file
// runsfromsavefile now contains the index in runlist for each run in save file
  }
  //--------------------------------------------------------------
  // Does this run match any in runlist?
  // return index in runlist, -1 if not found
  int RunsFromSavefile::RunNotFound(const std::vector<Run>& runlist,
                                    const std::vector<int>& batchnumbers) const
  {
    int nbn = batchnumbers.size();
    for (size_t ir=0;ir<runlist.size();++ir) { // loop runs
      std::vector<int> bl = runlist[ir].BatchList();
      if (int(bl.size()) == nbn) { // same size
        if (bl == batchnumbers) { // all same
          return ir;
        }}
    }
    return -1;
  }
  //--------------------------------------------------------------
}
