// restore.hh
//
// classes for RESTORE

#ifndef RESTORE_HEADER
#define RESTORE_HEADER

#include "fileread.hh"
#include "runthings.hh"

namespace scala {
  class RunsFromSavefile {
    //! extract run definitions from save (dump) file
  public:
    RunsFromSavefile(){}
    RunsFromSavefile(Fileread& FR, const std::vector<Run>& runlist)
    {init(FR, runlist);}

    void init(Fileread& FR, const std::vector<Run>& runlist);

    int NumberRunsInSaveFile() const {return svnruns;}

    int NumberRunsFound() const {return nrfound;}

    // for each run in save file, set to run index in runlist, or -1 if not used
    std::vector<int> Runsfromsavefile() const {return runsfromsavefile;}

  private:
    int svnruns;  // number of runs in save file
    int nrfound;  // number accepted from save file

    int nruns;    // number of runs in runlist = number to be used
    // for each run in save file, set to run index in runlist, or -1 if not used
    std::vector<int> runsfromsavefile;


    // Does this run match any in runlist?
    // return index in runlist, -1 if not found
    int RunNotFound(const std::vector<Run>& runlist,
		    const std::vector<int>& batchnumbers) const;

  };
} // namespace scala
#endif
