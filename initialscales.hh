//
// initialscales.hh
//
// Initial rough scaling by making average intensities equal
//

#ifndef INITIALSCALES_HEADER
#define INITIALSCALES_HEADER

#include "hkl_unmerge.hh"
#include "scalemodel.hh"
#include "controls.hh"
#include "Output.hh"


typedef std::pair<float,float> FPair;
typedef std::pair<double,double> DPair;

namespace scala {

class InitialData 
// Class to store initial scaling data, for FoxHolmes scaling
{
public:
  InitialData();
  InitialData(const clipper::Array2d<double>& AvI);

  int Npar() const {return np;}    // N valid rotation ranges 
  int NparAll() const {return npall;} // total N rotation ranges

  // return for each rotation range true if there are data, else false
  std::vector<bool> validRanges() const {return validranges;}

  // number of occupied resolution ranges for each rotation range
  std::vector<int> rangeCount() const {return rangecount;}

  // Data are in 2D array AvI(rotation, resolution)
  // For each resolution bin, we want to make all the <Irot> equal over
  //   all rotation bins  ie <I(rot,reso)>/scales(rot)  constant for all "rot"

  // return one observation, I, sigma pair, length np
  // Return false if end of data
  bool ObsArray(std::vector<DPair>& obs) const;

private:
  int npall;  // Number of parameters = nrotranges
  int np;     // Number of active parameters <= nrotranges
  // Data is a 2D array(nrotranges, nresbins)
  const clipper::Array2d<double>*  avi;
  std::vector<bool> validranges;
  std::vector<int> rangecount; // number of occupied resolution ranges for each rotation range
  mutable int next;
};


  class InitialScales {
  public:
    InitialScales() : status(-1) {}
    // Get initial estimates of primary scales, from making intensity
    // averages equal
    InitialScales(hkl_unmerge_list& hkl_list, ScaleModel& AllScales,
		  const all_controls& controls,
		  phaser_io::Output& output);

    // return true if there seems to be enough data to refine scales
    // If minimum_overlap <= 0.0, no check is made
    bool enoughData(const double& minimum_overlap,
		    const int& maximum_gap) const;


    // Report overlap status information to XML
    void reportOverlapXML(phaser_io::Output& output) const;

  private:
    std::vector<int> numobsrotrange; // number of observations for rotrange
    std::vector<double> fractionaloverlapbyrotrange;
    double averageoverlap;
    mutable double overlapthreshold;  // minimum allowed fractional overlap
    double minimumoverlap;            // minimum  fractional overlap
    mutable int allowedgap;                   // maximum allowed gap
    std::vector<double> gscales;  // inverse scales (g)

    // = -1, never called, = 0 no scales determined, > 0 n scales determined
    int status;

  };
}


#endif
