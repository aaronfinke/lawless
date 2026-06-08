//
//  batchgroup.hh
//
// Optional mapping of Batch number to batch analysis serial, grouping
// batches to sensible ranges (eg ~1 degree) for analysis

#ifndef BATCHGROUP_HEADER
#define BATCHGROUP_HEADER

#include "hkl_datatypes.hh"
#include "controls.hh"
#include "hkl_unmerge.hh"
#include "scalemodel.hh"

#include <string>

namespace scala {
  //--------------------------------------------------------------
  class Batchgroup {
    // A class to map actual batch number 
  public:
    Batchgroup() : nbatches(-1), ngroup(0) {}

    // allscales             scale model, for maximum primary scale smoothing range
    // analysiscontrols  for batchgroupwidth
    // hkl_list          address store for batch number -> batch serial lookup
    // datasetIndex      < 0 for all datasets, see statistics
    Batchgroup(const ScaleModel& allscales,
	       const AnalysisControls& analysiscontrols,
	       const hkl_unmerge_list& hkl_list,
	       const int& datasetIndex);

    void init(const ScaleModel& allscales,
	      const AnalysisControls& analysiscontrols,
	      const hkl_unmerge_list& hkl_list,
	      const int& datasetIndex);

    // Dummy for testing, just set ngroups
    //  Assume batch numbers are 1 to Ngroups
    Batchgroup(const int& Ngroups);

    // return number of groups, = 0 if unset
    int numberofgroups() const {
      if (nbatches <= 0) {return ngroup;}
      return batchindex.size();
    }

    // return group number for batch number
    int batchgroup(const int& batchnumber) const;

    // return (central) batch serial for group
    int batchserial(const int& jgroup) const;

    // return group width
    double groupWidth() const {return batchgroupwidth;}

    // return number of batches in each group
    int numberInGroup() const {return numberinGroup;}

    int numGroupSmooth() const {return numgroupsmooth;}

    // Return batch list
    std::vector<Batch> batches() const;

    std::vector<Run> runlist() const {
      return hkl_list_pointer->RunList();
    }

    // return dataset index
    int getdatasetindex() const {return datasetindex;}

    // return run for group
    int groupRun(const int& igroup) const
    {return grouprun.at(igroup);}

    // Number of batch groups in run
    int nGroupsinRun(const int& irun) const
    {return ngroupsinrun.at(irun);}

    // List of groups in run
    std::vector<int> groupsinrun(const int& irun) const;

    // format for printing
    std::string format() const;

    // format for XML
    std::string formatXML() const;

  private:
    int nbatches;
    int ngroup;
    double batchgroupwidth;
    int numberinGroup;
    const hkl_unmerge_list* hkl_list_pointer;

    int numgroupsmooth; // smoothing across batches

    // group number for each batch serial
    std::vector<int> groupindex;
    // First batch serial in each group (most groups of numberinGroup)
    std::vector<int> batchindex;
    std::vector<int> nbatchingroup; // number of batches in each group
    std::vector<int> ngroupsinrun; // number of groups in each run
    std::vector<int> grouprun;     // run for each group
    bool validPhi;
    bool differentPhiRanges;  // true if not all batches have same range
    // true if each batch is its own group, either because
    // specified (batchgroupwidth <= 0) or invalid Phi
    bool individualbatches;
    int datasetindex;

    // Given a list of contiguous batch numbers, assign them to groups
    //   returns number of groups
    int makeGroups(const std::vector<int>& batchnumberlist,
		   const int& irun);

  };
  //--------------------------------------------------------------
  class Xbreaks {
  public:
    Xbreaks(){}
    Xbreaks(const std::vector<Batch>& batches,
	    const int& datasetIndex);
    // datasetIndex = -1 for all datasets

    Xbreaks(const Batchgroup& batchgroup,
	    const int& datasetIndex);
    // datasetIndex = -1 for all datasets

    void init(const std::vector<Batch>& batches,
	      const int& datasetIndex);
    
    std::vector<Range> get_breaks() const {return breaks;}
    IntRange get_batchnumberrange() const {return validbatchnumbers;}
    
  private:
    std::vector<Range> breaks;
    IntRange validbatchnumbers;  // 1st and last actual accepted batch numbers, for x-axis range
  };
  //--------------------------------------------------------------
  //  symbol size  = 0 if "too many" points
  int GetSymbolSizeforNpoints(const int& npoints);
}
#endif
