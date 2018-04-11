//
// secondaryscalestats.hh
//
// Statistics about the applied and calculated secondary scales
//

#ifndef SECONDARYSCALESTATS_HEADER
#define SECONDARYSCALESTATS_HEADER

#include <vector>
#include "Output.hh"
#include "hkl_unmerge.hh"
#include "scalemodel.hh"
#include "histogram.hh"

#include <clipper/clipper.h>


namespace scala {

  class SecondaryScaleStats {
  public:
    SecondaryScaleStats(){}
    void init(const ScaleModel& allScales);

    // go through list to get secondary scale statistics
    void getScaleStats(hkl_unmerge_list& hkl_list,
		       const bool& onlyUseSingletons);

  // Print secondary corrections as 2D array
    void PrintSecondaryCorrections(phaser_io::Output& output) const;

  private:
    double angleinterval;  // angle interval for analysis on theta and phil

    // Table of calculated corrections
    int nsecscales;
    // Phi values
    int nphi;
    std::vector<double> phivalues;
    // Theta values
    int ntheta;
    std::vector<double> thetavalues;
    std::vector<clipper::Array2d<double> > secscales;
    // count of number in each theta, phi bin
    std::vector<clipper::Array2d<int> > sscount;
    Range secscalerange;
    Histogram secscalehisto; // histogram of values

    const ScaleModel*  allscales;

    void setupHistogram();
    void PrintHistogram(phaser_io::Output& output) const;

    std::string secscaletype(const std::string& sectype) const;

  };
}
#endif
