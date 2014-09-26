// applyscales.h
//
// Testing for reading reflections & getting scales
//

#ifndef APPLYSCALES_HEADER
#define APPLYSCALES_HEADER

#include "scalemodel.hh"
#include "hkl_unmerge.hh"
#include "Output.hh"

namespace scala {
  //----------------------------------------------------------------
  class UnusualObservation {
    // record of an observation which has eg the maximum correction from sd(k)
  public:

    UnusualObservation() : batch(-1) {}

    UnusualObservation(const observation& this_obs, const Rtype& Invresolsq);
    void init(const observation& this_obs, const Rtype& Invresolsq);

    std::string format() const;

  private:
    Hkl hkl;
    Rtype invresolsq;
    int batch; // batch number
    Rtype gscale;
    Rtype varg;  // Var(g)
    Rtype I;     // scaled intensity
    Rtype sigI0; // unscaled sigI
    Rtype sigI;  // sigI after scaling and correction by sd(k)
  };
  //----------------------------------------------------------------
  class ApplyScales {
  public:
    ApplyScales() : maxrelsddiff(-1000.0) {}

    //! Apply scales to all data, return mean(I) in resolution limits
    // if onlyUseSingletons true, do not attempt to apply scales to overlaps
    ApplyScales(const ScaleModel& AllScales, hkl_unmerge_list& hkl_list,
		const bool& onlyUseSingletons);

    void scale(const ScaleModel& AllScales, hkl_unmerge_list& hkl_list,
	       const bool& onlyUseSingletons);

    // return mean scaled I
    double meanI() const {return meani.Mean();}

    // Mean sd(1/g)
    std::vector<MeanValue> meanSDk() const {return meansdk;}
    // Mean sigIcorrected/sigI
    std::vector<MeanValue> meanRelSDdiff() const {return meanrelsddiff;}

    void print(phaser_io::Output& output) const;

  private:
    MeanValue meani;  // mean scaled I within resolution limits
    //  in resolution bins ...
    ResoRange resrange;
    std::vector<MeanValue> meansdk;  // Mean sd(1/g) * g
    std::vector<MeanValue> meanrelsddiff;  // Mean sigIcorrected/sigI
    // sigIcorrected = sd' allows for sd(1/g)

    Rtype maxrelsddiff; // Maximum (sd' - sd)/sd'
    Rtype maxsdk;
    UnusualObservation maxsdkobs;
  };

}

#endif
