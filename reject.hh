//
// reject.hh
// 
// Control of volatile acceptance flags on observations
// Each observation has:
// 1) an ObservationFlag
//   This comes from the input (MTZ) file from Mosflm & flags various error
//   conditions. It is static as read from the file
// 2) an ObservationStatus
//   This is volatile, & includes
//   (a) a flag based on the ObservationFlag, but with conditional acceptance
//   (b) outlier flag
//   (c) Emax flag
//   (d) flag for too strong or too weak for the scaling pass
//   

#ifndef REJECT_HEADER
#define REJECT_HEADER

#include "hkl_unmerge.hh"
#include "controls.hh"
#include "sdmodel.hh"
#include "plotfiles.hh"
#include "normalise.hh"

namespace scala {
  // ------------------------------------------------------------
  class WriteRogues
  // Write rogues entry
  {
  public:
    WriteRogues() :rogues(NULL) {}  // Null entry, no output

    // Open ROGUES file & write header if Start true
    // Open ROGUESPLOT file & write header if Plot true
    // title & maximum resolution d* = lambda/d
    //  outliercontrol   parameters for rejection
    WriteRogues(const std::string& filename,
		const bool& Start, const bool& Plot,
		const bool& multilattice,
		const std::string& title,
		const float& dstarMax, const float& wavelength,
		const Rings& icerings,
		const OutlierControl& outliercontrol,
		const bool& xmgraceoutput);

    bool Open() const {return (rogues != NULL);}

    // Write one rogue reflection to file
    void RogueReflection(const reflection& this_refl,
			 const std::vector<float>& deviations,
			 const Normalise& NormRes);

    void End() {
      if (rogueplot.IsPlot()) {rogueplot.End();}
    }

    std::string formatXML() const {return rogueplot.formatXML();}
    
  private:
    FILE* rogues;
    RoguePlot rogueplot;
  };
  // ------------------------------------------------------------
  //  Clear outlier status flags for all observations
  //  Other flags are left unaltered
  void ClearOutlierFlags(hkl_unmerge_list& hkl_list);
  // ------------------------------------------------------------
  // Check for outliers in all reflections, using parameters in outliercontrol,
  // and set status flags as required on each observation
  //
  //  SDM              current sd correction model
  //  anomOn           true to do main outlier check only within the
  //                   I+ and I- sets
  //  outliercontrol   parameters for rejection
  //  RoguesList       optional rogues list output
  //  
  void RejectOutlier(hkl_unmerge_list& hkl_list,
		     const SDmodel& SDM, const Normalise& NormRes,
		     const bool& anomOn,
		     const OutlierControl& outliercontrol,
		     WriteRogues& RoguesList);
  // ------------------------------------------------------------
  // Returns counts of flagged outliers within I+/-, between +/- and on Emax
  std::vector<int> CountOutliers(const hkl_unmerge_list& hkl_list);
  // ------------------------------------------------------------
  // Returns counts of flagged outliers as XML
  std::string CountOutliersXML(const std::vector<int>& nrejs);
  // ------------------------------------------------------------
  // Returns counts of flagged outliers within I+/-, between +/- and on Emax
  // return[0] number of rejects [1] number on I+- [2] number on Emax
  // On exit:
  //  rejectedbatch  count of rejected reflections for each batch
  //  rejecteddataset count of rejected reflections for each dataset
  std::vector<int> CountOutliers(const hkl_unmerge_list& hkl_list,
				 std::vector<int>& rejectedbatch,
				 std::vector<int>& rejecteddataset);

}
#endif
