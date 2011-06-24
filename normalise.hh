// normalise.hh
// Phil Evans August 2003-2004 etc

#ifndef SCALA_NORMALISE
#define SCALA_NORMALISE

#include "hkl_merged_list.hh"
#include "hkl_unmerge.hh"
#include "spline.hh"
#include "icering.hh"
#include "score_datatypes.hh"
#include "linearlsq.hh"
#define ASSERT assert
#include <assert.h>

namespace scala
{
//--------------------------------------------------------------
  class BinSums {
  public:
    BinSums(const Rtype& Time0=0.0)
      : sum_I(0.0), sum_sSqr(0.0), N(0), sum_time(0.0),
	max_time(0.0)
    {time0=Time0;}  // construct & zero

    void clear(const Rtype& Time0=0.0)
    {sum_I=0.0;sum_sSqr=0.0;N=0;sum_time=0.0;
      time0=Time0;max_time=0.0;}

    // Add in I, s^2, time
    void add(const Rtype& I, const Rtype& sSqr, const Rtype& Time=0.0);

    // <I>
    double MeanI() const;
    // <s^2>
    double MeanSSqr() const;
    // N
    int num() const {return N;}
    // <time>
    double MeanTime() const;
    // time0
    double Time0() const {return time0;}
    // Maximum time
    double MaxTime() const {return max_time;}

  private:
    // sums (and later means) for <I> and <s^2>
    double sum_I;
    double sum_sSqr;
    int N;
    // "Time" information
    double sum_time;
    double time0;
    double max_time;
  };
  //--------------------------------------------------------------
  class MeanIsdIsSqr {
    // <I>, <sdI> and <sSqr> eg for a resolution bin
  public:
    MeanIsdIsSqr()  : MnI(0.0), SdMnI(0.0), MnSdI(0.0), SdMnSdI(0.0), MnsSqr(0.0), N(0) {}

    double MnI;
    double SdMnI;
    double MnSdI;
    double SdMnSdI;
    double MnsSqr;
    int N;
  };
  //--------------------------------------------------------------
  class BfactorModel {
    // Model for Bfactor type conversion factor
    //
    // This model is for a Wilson-type scale factor, but with the
    // B-factor a linear function of "time" (optionally fixed)
    //
    // c(s^2, t) = c0 exp[(B0 + B*t)*s2]
    //
    // Three parameters: c0 (determinated as ln(c0)); B0; B
  public:
    BfactorModel()
      : valid(false), time0(0.0),  dtime(0.0), nparam(3) {}
    // Constant = true for B constant with time
    BfactorModel(const bool& Constant) : valid(false), dtime(0.0) {
      Setup(Constant);}
    
    // Constant = true for B constant with time
    void Setup(const bool& Constant, const double& Time0=0.0);

    int Nparam() const {return nparam;}
    
    // Store parameters
    void Set(const std::vector<double>& params);

    // Mark as invalid (no data)
    void Invalid() {valid = false;}
    bool Valid() const {return valid;}

    // Return measurement vector = (1.0  sSqr,  sSqr*t)
    std::vector<double> MeasurementVector(const double& sSqr, const double& t);

    // constant term
    double C0() const {return exp(lnC0);}
    // B0
    double Bfac() const {return B0;}

    // Return parameters
    std::vector<double> Params() const;

    // time length
    double DTime() const {return dtime;}  // get
    double& DTime() {return dtime;}  // set

    // Return dividing correction factor for resolution sSqr &
    // "time" (actual time, not relative to start
    double Factor(const double& sSqr, const double& t=0.0) const;
    // time relative to start
    double FactorT(const double& sSqr, const double& t=0.0) const;
    
  private:
    bool valid;
    double lnC0;   // ln(c0)
    double B0;
    double B;
    double time0;
    double dtime;  //  time length
    int nparam;
  };
  //--------------------------------------------------------------
  class Normalise
  // Normalise intensities with respect to resolution
  //
  // Conversion factor to convert I to E is done by run
  // Within each run, the dividing factor is a function of s2 = sSqr = 1/d^2 and
  // "time" (typically relative batch number is a proxy for time)
  //
  //    c(s2, t) = c0 exp[(B0 + B*t)*s2]
  // ie the "B-factor" is a linear function of time t
  //
  // For determination, this is linearised by taking logs
  //    ln(c) = ln(c0) + B0*s2 + B*s2*t
  //
  // Normalisation is then E = (1/c) * I
  //  
  //
  // Type +1  just B-factor
  // Type +2  B-factor + smoothed bins (over all runs
  //
  //   a) construct from vector of bin average <I> & bin average <s^2>
  //   (binned on resolution & "time")
  //   b) calculate & return overall B-factor
  //         Boverall = .calcB()  (this done in the constructor)
  //   c) smooth bins (spline)
  //   d) apply correction
  //         Icorr = .apply(I, sSqr)
  //
  // If number of "runs" == 1 == number of batches, then a single normalisation
  // scheme is used for all data, ie the data sre assumed to be already scaled.
  // In this case, the function "applAvg" must be used to apply the normalisation
  {
  public:
    Normalise() : validRun(false), validAll(false) {type=-1;}

    Normalise(const std::vector<clipper::Array2d<BinSums> >& sums_st);

    // set max E^2
    void SetE2max(const float& EEmax) {E2max = EEmax;}

    // return SclCorr
    float SclCor() const {return AvgFactor.C0();}

    // return average B correction
    float Bcorr() const {return AvgFactor.Bfac();}

    // Total correction factor (multiplying factor)
    float Corr(const float& sSqr,
		const int& irun, const float& time) const;

    std::vector<BfactorModel> BfactorCorr() const {return Bfactors;}


    // return normalised intensity or I/sigI
    // "time" here is actual time or batch, not relative to start of run:
    // offset to make relative to start of run is done internally
    float apply(const float& I, const float& sSqr,
		const int& irun, const float& time) const;
    IsigI apply(const IsigI& Is, const float& sSqr,
		const int& irun, const float& time) const;
    // Average factors
    float applyAvg(const float& I, const float& sSqr) const;
    IsigI applyAvg(const IsigI& Is, const float& sSqr) const;

    // <I> overall input data
    float Imean() const {return imean;}
    // maximum intensity I
    float Imax() const {return imax;}

    //! Store <I>, <sdI>, <sSqr> in resolution bins
    void StoreMeanIsdIsSqr(const std::vector<MeanIsdIsSqr>& Meanisdissqr) {
       meanisdissqr = Meanisdissqr;
    }
    //! return <I>, <sdI>, <sSqr> in resolution bins
    std::vector<MeanIsdIsSqr>  GetMeanIsdIsSqr() const {return meanisdissqr;}

    bool NotTooLarge(const float& EE) const {
      return (EE < E2max);
    }

    bool NotTooLarge(const IsigI& Is) const {
      return (Is.I() < E2max);
    }

  private:
    bool validRun; // true for valid normalisation by run & batch (function "apply")
    bool validAll; // true for valid overall normalisation (no run/batch dependence)(function "applyAvg")
    int type;
    float bcmin;    // minimum spline correction
    int Nbins;
    int Nruns;
    std::vector<LinearLSQ> RunFactor;    // linear fit for each run
    std::vector<BfactorModel> Bfactors;  // B-factors etc for each run
    float SclCorr;  // Average constant term
    BfactorModel AvgFactor;   // Average factor (valid for t=0)

    float imean; // overall <I>
    float imax;  // maximum intensity

    std::vector<BinSums> sSqrmnI;
    std::vector<RPair> sSqrmnIcorr;
    Spline bincorr;
    float E2max;  // Maximum allowed E**2

    std::vector<MeanIsdIsSqr> meanisdissqr; // <I>, <sdI>, <sSqr> in resolution bins

  };  // class Normalise

// Set up intensity normalisation object NormRes
//
// Use binned <I> to get normalisation object
// Optionally output plot to file lnI.plot (if Printlevel > 0)
// Optionally, if MinIsigRatio > 0.0
//   return updated ResoRange, setting resolution cutoffs and bins
//   (if MinIsigRatio < 0, don't change resolution range)
//
// If Overall true, generate a single normalisation for all runs, ie
// assume the data is already scaled
//   
// MinIsigRatio   minimum I/sigI ratio on unaveraged data
//                if < 0, do not reset resolution range
Normalise SetNormalise(const hkl_unmerge_list& ref_list,
		       const double& MinIsigRatio,
		       const bool& Overall,
		       ResoRange& ResRange,
		       Rings& Icerings,
		       const int PrintLevel);

// Set up intensity normalisation object NormRes from merged data
//
// Use binned <I> to get normalisation object
// Optionally, if MinIsigRatio > 0.0
//   return updated ResoRange, setting resolution cutoffs and bins
//
// MinIsigRatio   minimum I/sigI ratio on unaveraged data
Normalise SetNormaliseMerged(const hkl_merged_list& ref_list,
			     const double& MinIsigRatio,
			     ResoRange& ResRange);

}


#endif
