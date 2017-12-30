
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
#include "batchgroup.hh"
#include "score_datatypes.hh"

namespace scala {

  class CompareSDs {
  public:
    CompareSDs(){}

    CompareSDs(const ResoRange& ResRange, const IntensityBin& Irange,
	       const int& Minimumsample, const Batchgroup& batchgroup);
    void init(const ResoRange& ResRange, const IntensityBin& Irange,
	      const int& Minimumsample, const Batchgroup& batchgroup);

    void add(SelectedObservations& selobs,
	     const int& mres, const int& mint); // add in contributions

    void printByResolution(phaser_io::Output& output) const;

    void printByIntensity(const IntensityBin& irange,
			  phaser_io::Output& output) const;

    void printByBatch(const int& datasetIndex,
		      phaser_io::Output& output) const;

    //! return minimum sample
    int MinimumSample() const {return minimumsample;}

  private:
    ResoRange resrange;
    const Batchgroup * pbatchgroup;

    std::vector<MeanValue> meanSDratioreso; // mean(SD/sampleSD) by resolution
    std::vector<MeanValue> meanSDratioInt;  // mean(SD/sampleSD) by intensity
    std::vector<correl_coeff> CCSDreso;     // CC(SD,sampleSD) by resolution
    std::vector<correl_coeff> CCSDint;      // CC(SD,sampleSD) by intensity

    std::vector<MeanValue> meanChiSqreso; // mean(ChiSq-Variance) by resolution
    std::vector<MeanValue> meanChiSqInt;  // mean(ChiSq-Variance) by intensity
    std::vector<MeanValue> meanChiSqSamplereso; // mean(ChiSq-Sample) by resolution
    std::vector<MeanValue> meanChiSqSampleInt;  // mean(ChiSq-Sample) by intensity

    // Batch group statistics: mean(ChiSq) = Mean(Delta2^2)
    std::vector<MeanValue> meanChiSqBtGp;       // mean(ChiSq-Variance) by batch group
    std::vector<MeanValue> meanChiSqSampleBtGp; // mean(ChiSq-Sample) by batch group

    int minimumsample;  // sample SD for more this number of observations

    MeanValue csq1;
    MeanValue csq2;
    MeanValue csq3;



  };
}

#endif
