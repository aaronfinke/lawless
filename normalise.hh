// normalise.hh
// Phil Evans August 2003-2004 etc
// Simplified August 2018

#ifndef SCALA_NORMALISE
#define SCALA_NORMALISE

#include "hkl_unmerge.hh"
#include "spline.hh"
#include "icering.hh"
#include "score_datatypes.hh"
#include "median.hh"

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
  // PrintLevel  if > 0, dump to norm.plot
  //
  {
  public:
    Normalise() : valid(false), useAniso(false) {}

    Normalise(const hkl_unmerge_list& hkl_list,
	      const double& MinIsigRatio,
	      Rings& Icerings,
	      const clipper::U_aniso_frac& u_aniso_frac,
	      const int PrintLevel,
	      bool final=false);

    // Store u_aniso_frac tensor (scaled), may be null
    //  tensor as returned from AnisotropicAnalysis.U_aniso_frac()
    // Stored here scaled by -twoPi^2 for use in quadratic form
    void setAniso(const clipper::U_aniso_frac& u_aniso_frac);
    // Disable anisotropy correction
    void noAniso() {useAniso = false;}

    // this should be called after setAniso or noAniso
    void init(const hkl_unmerge_list& hkl_list,
	      const double& MinIsigRatio,
	      Rings& Icerings,
	      const int PrintLevel,
	      bool final=false);

    //  true if there is a valid normalisation factor at this resolution
    bool validResolution(const float& sSqr) const; //? reject
    
    // Total average correction, multiplying scale
    // This is the scale needed to bring an individual I to match the average
    float Corr(const float& sSqr,
	       const DVect3& rhkl=DVect3()) const;
    // anisotropic part of correction, multiplying scale
    // This is the scale needed to bring an individual I to match the average
    float anisoCorr(const DVect3& rhkl) const;

    // Apply correction factors to get E^2 from I
    //  The relevant <I> for hkl is anisoCorr/<I>, so divide by that
    float apply(const float& I, const float& sSqr,
		const DVect3& rhkl=DVect3()) const;
    IsigI apply(const IsigI& Is, const float& sSqr,
		const DVect3& rhkl=DVect3()) const;



    // <I> overall input data
    double Imean() const {return imean;}  //? Ibinning
    // maximum intensity I (weighted mean)
    double Imax() const {return imax;}

    // mean/median ratio, averaged over some low resolution bins
    double mmratio() const;

  private:
    bool valid; // true for valid normalisation

    ResoRange resorange;  // internal resolution range, with more bins
    int nrbin;  // number of resolution bins

    std::vector<double> mnsSqr;  // <sSqr> by resolution bin
    //  <I> by resolution bin from trimmed range, aniso corrected if useAniso true
    std::vector<double> mnI;
    std::vector<double> sdI;     // corresponding sd(<I>) from weights
    std::vector<double> medianI; //  median(I) by resolution bin
    std::vector<double> anisomedI; //  median(I) by resolution bin aniso corrected
    std::vector<int> mcount;     // number used for mnI

    double imean; // overall <I>
    double imax;  // maximum intensity

    Spline bincorr;
    float E2max;  // Maximum allowed E**2

    // Scaled by twoPi^2 for use in quadratic form,
    // and by 2.0 to allow for scaling intensities rather than amplitudes
    clipper::U_aniso_frac u_aniso_frac_scaled;
    bool useAniso;

    // Internal methods
    void setstores();

    void store(const int& ibin, const double& sSqr,
	       const MeanVariance& mnv,
	       Median<float>& medI,
	       Median<float>& anisoMedI);

    // Weak high resolution bins are unreliable, so (pending a better method)
    // replace <I> by a value extrapolated from the last accepted bin
    // Very crude!!   NOT USED NOW, replaced by fitLogCurve

    void resetweak(const double& minIsigRatio);

    // fit exponential curve to data beyond resolution of ~3A
    std::vector<double> fitLogCurve(const double& minIsigRatio);


    // return I/sigI for resolution bin
    double iovsig(const int& ibin) const;

    void dump(const std::string& filename="");

  };  // class Normalise


}

#endif
