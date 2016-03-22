//
// comparesds.hh
//
// Compare sample sd(I) estimates with corrected input sd(I) estimates
//

#ifndef COMPARESDS_HEADER
#define COMPARESDS_HEADER

#include <vector>

#include "range.hh"
#include "selectedobservations.hh"
#include "intensitybin.hh"

namespace scala {

  class CompareSDs {
  public:
    CompareSDs(){}

    CompareSDs(const ResoRange& ResRange, const IntensityBin& Irange);
    void init(const ResoRange& ResRange, const IntensityBin& Irange);

    void add(SelectedObservations& selobs,
	     const int& mres, const int& mint); // add in contributions

    void printByResolution(phaser_io::Output& output) const;

    void printByIntensity(const IntensityBin& irange,
			  phaser_io::Output& output) const;


  private:
    ResoRange resrange;
    std::vector<MeanValue> meanSDratioreso; // mean(SD/sampleSD) by resolution
    std::vector<MeanValue> meanSDratioInt;  // mean(SD/sampleSD) by intensity


  };
}

#endif
