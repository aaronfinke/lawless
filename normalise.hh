// normalise.hh
// Phil Evans August 2003-2004 etc
// Simplified August 2018

#ifndef SCALA_NORMALISE
#define SCALA_NORMALISE

#include "hkl_unmerge.hh"
#include "spline.hh"
#include "icering.hh"
#include "score_datatypes.hh"
#define ASSERT assert
#include <assert.h>

namespace scala {
  //--------------------------------------------------------------
  class Normalise
  // Normalise intensities with respect to resolution
  //
  // Use binned <I> to get normalisation object
  // Binned of resolution (from ResRange) over all batches and runs
  // Optionally output plot to file lnI.plot (if Printlevel > 0)
  //
  // MinIsigRatio   minimum I/sigI ratio on averaged data
  //                ranges beyond this threshold get reset
  //                if < 0, no reset
  //
  {
  public:
    Normalise() : valid(false) {}

    Normalise(const hkl_unmerge_list& hkl_list,
	      const double& MinIsigRatio,
	      Rings& Icerings,
	      const int PrintLevel);


    //  true if there is a valid normalisation factor at this resolution
    bool validResolution(const float& sSqr) const; //? reject

    // Total average correction, multiplying scale
    float Corr(const float& sSqr) const;

    // Correction factors
    float apply(const float& I, const float& sSqr) const;
    IsigI apply(const IsigI& Is, const float& sSqr) const;

    // <I> overall input data
    double Imean() const {return imean;}  //? Ibinning
    // maximum intensity I
    double Imax() const {return imax;}

    void dump() const;

  private:
    bool valid; // true for valid normalisation

    ResoRange resorange;  // internal resolution range, with more bins
    int nrbin;  // number of resolution bins

    std::vector<double> mnsSqr;    // <sSqr> by resolution bin
    std::vector<double> mnI;       //  <I> by resolution bin
    std::vector<double> sdI;       // corresponding sd(<I>)
    std::vector<int> mcount;

    //^^
    std::vector<MeanVariance> wmnI;
    std::vector<MeanValue> unwmnI;


    double imean; // overall <I>
    double imax;  // maximum intensity

    Spline bincorr;
    float E2max;  // Maximum allowed E**2

    // Internal methods
    void setstores();

    void store(const int& ibin, const double& mnsSqr,
	       const MeanVariance& mnv);

    // Weak high resolution bins are unreliable, so (pending a better method)
    // replace <I> by a value extrapolated from the last accepted bin
    // Very crude!!
    void resetweak(const double& minIsigRatio);

    // return I/sigI for resolution bin
    double iovsig(const int& ibin) const;

  };  // class Normalise


}

#endif
