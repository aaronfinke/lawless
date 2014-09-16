// optimisecombine.hh
//
// 

#ifndef OPTIMISECOMBINE_HEADER
#define OPTIMISECOMBINE_HEADER

#include <vector>

#include "hkl_unmerge.hh"
#include "score_datatypes.hh"
#include "range.hh"

namespace scala {

  class OptimiseCombine {
    //! Choose or optimise combination of intensities from summation integration and profile-fitting, if present

    // Mosflm produces 2 estimates of each intensity, from just summing pixels (Isum) or from
    // profile fitting (Ipr). In some cases Isum is better for strong intensities while Ipr is always
    // better for weak ones. In this, case the best option is probably to use a combination depending
    // on the "raw" intensity Iraw, ie the intensity (of summed partials) back-corrected for the LP
    // factor (ie divided by the LP factor in the file from Mosflm) (INTENSITIES COMBINE)
    //
    // then  I = w Ipr + (1-w) Isum
    // where w = 1/(1 + (Isumraw/Imid^Ipower))
    //       Isumraw = Isum/LP       Ipower = 3 by default (not changed here)
    //
    // Optimisation is based on the best overall Rmeas
    // 
    // The static class SelectI is set to choose one or the combination 
  public:
    OptimiseCombine():bothIpresent(false){}

    //! Construct and set SelectI to the "optimum"
    OptimiseCombine(const hkl_unmerge_list& hkl_list,
		    phaser_io::Output& output);
    void init(const hkl_unmerge_list& hkl_list,
	      phaser_io::Output& output);

    bool IsOptimised() const {return bothIpresent;}

  private:
    // Set imid value (= 0 Isum, < 0 Ipr, > 0 combine), set Rmeas by resolution &
    // return overall Rmeas
    double TestValue(const double& imid, const hkl_unmerge_list& hkl_list,
		     std::vector<Rfactor>& Rmeas);

    void PrintR(const double& flag,
		const double& score, const std::vector<Rfactor>& Rmeas,
		phaser_io::Output& output);

    // Return score, Rfactors vs. resolution    
    std::vector<Rfactor> GetScores(const hkl_unmerge_list& hkl_list) const;
    // Return Mean(Iraw)
    double RawIntensityDistribution(const hkl_unmerge_list& hkl_list) const;

    // private data
    bool bothIpresent;     // true if both intensity estimates are present in the file
    double meanI;          // mean raw intensity
    double Itop;           // top of search range

    ResoRange ResRange;             // resolution range & bins
    int nresbin;                    // number of resolution bins
    std::vector<Rfactor> RmeasBest; // best Rmeas score by resolution
    double Rmbest;                  // best Rmeas overall
    
  }; // class OptimiseCombine
}

#endif
