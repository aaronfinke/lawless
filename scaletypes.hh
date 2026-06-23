// scaletypes.hh
//
//  Scale types for scale model
//
//  

#ifndef SCALETYPES_HEADER
#define SCALETYPES_HEADER

#include <vector>
///#include "InputAll.hh"
///#include "keywords_aimless.hh"
///#include "hkl_unmerge.hh"
#include "sphericalharmonic.hh"
#include "tie.hh"
#include "fileread.hh"
#include "tile.hh"
#include "hash.hh"


namespace scala {
  //--------------------------------------------------------------
  class SmoothedValue
  {
    // Provides a smoothed moving average value over a range of "normalised"
    // coordinate
    //
    // For normalised coordinate z which runs over the range 0 to zmax, we
    // have values at positions -0.5, +0.5, 1.5, ... Nint(zmax)+0.5 (ie extended
    // beyond actual range at both ends). Number of points Np = Ns + 2
    // Normalised coordinate z is derived from raw coordinate x as
    //    z  =  (x - x0)/Dx
    // Dx (= spacing) is adjusted to make an integral multiple into x range
    // ie Dx = (x1 - x0) / Ns
    // Raw coordinates range from x0 to x1 corresponding to z = 0 -> Ns
    // 
    // The moving average is taken over Nav points
    //   (at present limited to 3, 4 or 5)
    // The weight at each point w(i) = exp[- d^2] where d = (z - z(i))/sigma
    //
    // Special cases:
    //  1) If number of intervals = 0, then value is a single constant
    //  2) If number of sample points Np < Nav, then Nav = Np
    //     Note that since Np = Ns + 2, even if Ns = 1, Np = 3, so Nav can be = 3
    //  3) Ends: always uses Nav points, abutting ends if necessary
    //  
  public:
    SmoothedValue(): nvalues(0) {}
    // Construct from range of raw unnormalised coordinate & number of sample intervals
    // Set smoothing values to defaults, Nav = 3
    SmoothedValue(const Range& Xrange, const int& Ns);
    // Set smoothing values: number of points, sigma1
    // If sigma < 0, set to "optimum" (!)
    // (or at least a suitable) value from Naverage
    // Naverage defaults to 3
    void SetSmoothing(const int& Naverage, const double& Sigma);
    // Store values (nvalues points)
    void StoreValues(const std::vector<float>& Values);
    void StoreValues(const std::vector<double>& Values);
    void StoreValue(const double& value); // set all to same value

    // Return number of values
    int Nvalues() const {return nvalues;}
    // Return number of sample intervals
    int Nsample() const {return nsample;}
    // Return all values
    std::vector<double> Values() const {return values;}
    // Return interpolated value at point (original unnormalised coordinate)
    double Value(const double& x) const;
    // Return interpolated value at point, plus weights at each point,
    // for original unnormalised coordinate
    double ValueWeight(const double& x, std::vector<double>& weight,
		      double& sumweight) const;
    // Return number of points averaged
    int Naverage() const {return naverage;}
    // Return sigma factor
    double Sigma() const {return sigma;}
    // Return spacing
    double Spacing() const {return spacing;}
    // Return positions
    std::vector<double> Positions() const {return positions;}

    // format for save/restore
    std::string FormatSave() const;
    // restore
    void Restore(Fileread& FR);

  private:
    int nvalues;      // Number of value points
    double xmin, xmax; // actual range of unnormalised coordinate
    double x0;         // position of 1st point in unnormalised coordinate (= xmin)
    double spacing;    // interval on raw coordinate
    int nsample;      // number of samples in x range
    std::vector<double> values;     // the values
    std::vector<double> positions;  // the positions of the values
				   //   z(i), i = 0, nvalues-1
    int naverage;     // Number of points to smooth
    double half_nav;   //   0.5 * naverage
    double sigma;      // smoothing weight
  };  // class SmoothedValue
  //--------------------------------------------------------------
  class PrimaryScale
  {
    // 1) Smooth scaling on phi (or equivalent)
    // To calculate the scale for an observation at "rotation angle" phi,
    // uses the class SmoothedValue
    // 
    // We have have nscales scale factors at intervals of scalespacing,
    // including one extra on each end of actual range (unless nscales == 1)
    //  at positions z(i) = -0.5, +0.5, 1.5, ... (nscales-1.5) in the "normalised"
    //  coordinate 
    // Normalised coordinate z = (phi - phi0)/scalespacing
    // scalespacing will be adjusted to make an integral number of intervals in the range
    // Then
    //  inverse scale g  = Sum(i) [ w(i) Scale(i) ] / Sum(i)[w(i)]
    //      for i=0,nscales-1
    // and w(i) depends on the distance between the coordinate of the
    // observation z & the coordinates for each scale factor z(i)
    //        d(i) = (z - z(i))/ SDwz    where SDwz is a "SD" factor
    // Use the 3, 4 or 5 points closest to z
    //        w(i) = exp(-d(i)^2)
    //
    // Partial derivatives:
    //   dg/Scale(i) = w(i) / Sum(i)[w(i)] 
    //
    // Smoothing is controlled by two parameters, which are set globally
    // (static)
    //  i)  NavgScale  number of points in moving average,
    //      == 3, 4 or 5 [default 3]
    //      The actual value used will be forced to be <= nscales
    //  ii) SDwz  "sigma" value for Gaussian weight: default set depending
    //      on NavgScale
    //
    // 2) Batch scaling (or if nscales == 1)
    //   g(j) = Scale(j)
    //   dg(j)/Scale(i) = 1 if j == i, else  = 0
    //
  public:
    PrimaryScale() : nscales(0) {}

    // Construct smooth scaling from number of scale intervals
    // Note that actual number nscales will be NscalesIntervals+2
    //   unless NscalesIntervals = 0, in which nscales = 1
    PrimaryScale(const int& NscaleIntervals, const Range& phirange);

    // Construct smooth scaling from scale spacing 
    //  scalespacing will be adjusted to give an integral number of intervals
    // Always at least 2 scales
    PrimaryScale(const double& scaleSpacing,
		 const Range& phirange);

    // Construct batch scaling from list of batch numbers (from run)
    PrimaryScale(const std::vector<int>& BatchNumbers);

    // Return list of ties: sdtie is sd for weight, idx0 is index to first global
    // parameter for setting ties, since they refer to the global parameter index
    std::vector<Tie> Ties(const double& sdtie, const int& idx0);

    // Store scales vector (length nscales)
    void StoreScales(const std::vector<double>& Scales);
    // Store Nobs vector (length nscales)
    void StoreNobservations(const std::vector<int>& Nobs);
    // Retrieve scales (length nscales)
    std::vector<double> Scales() const;
    // Retrieve number of contributions (length nscales)
    std::vector<int> Nobservations() const;

    // Number of scales, == 0 if null
    int Number() const {return nscales;}
    // Number of intervals, minimum 1
    int Nintervals() const {return Max(1,nscaleintervals);}
    // Interval
    double Spacing() const {return scalespacing;}

    // Smooth scaling:
    // Scale at position phi
    double Scale(const double& phi) const;
    // Return scale and derivative vector at position phi
    void ScaleDeriv(const double& phi, double& scale, std::vector<double>& dgdp) const;
    // Return scale and derivative vector at position phi
    double ScaleDeriv(const double& phi, std::vector<double>& dgdp) const;

    // Batch scaling:
    // Scale at batch number batch
    double Scale(const int& batch) const;
    // Return scale and derivative vector at batch number batch
    void ScaleDeriv(const int& batch, double& scale, std::vector<double>& dgdp) const;
    // Return scale and derivative vector at batch number batch
    double ScaleDeriv(const int& batch, std::vector<double>& dgdp) const;

    // true if Batch mode
    bool IsBatchScale() const {return batchscale;}

    // valid scale at batch number batchnum
    bool ValidScale(const int& batchnum) const;

    // index number in this batch list for batch number batchnum, omitting rejects
    int batchSerialIndex(const int& batchnum) const;

    // set ib'th batch if Usebatch[ib] false, batchnumbers corresponding
    void setBatchReject(const std::vector<bool>& Usebatch,
			const std::vector<int>& batchnumbers);

    // Format layout for printing
    std::string format() const;

    // Format all information into a labelled save format for later restoration
    std::string FormatSave() const;

    // restore
    void Restore(Fileread& FR);

    // Set global smoothing parameters
    // NB these must be set before constructing an object, otherwise they
    // have no effect, & the default values will be used
    static void SetNavgScale(const int& NavgScale);
    static int NavgScale() {return navgscale;}
    //  Set SDwz < 0 to use default value
    static void SetSDwz(const double& SDwz) {sdwz = SDwz;}
    static double SDwz() {return sdwz;}

  private:
    void init(const Range& phirange);
    void ScaleDeriv(const bool& Deriv, const double& phi,
		    double& scale, std::vector<double>& dgdp) const;
    void ScaleDeriv(const bool& Deriv, const int& batch,
		    double& scale, std::vector<double>& dgdp) const;

    bool batchscale;     // batch | smooth
    int nscaleintervals; //BS number of scale intervals, = 0 for single scale
    int nscales;         //BS actual number of scale parameters
    double scalespacing; //S
    double phi0;         //S
    hash_table batch_lookup;  //B batch lookup for batch scaling

    SmoothedValue smoothscale;  //S for smoothed scales
    std::vector<double> batchscales;  //B batch scales
    std::vector<int> nobsPar;   // BS number of observations for each scale
    bool allbatches;            // true if all batches are used
    std::vector<bool> usebatch; // flag to use this batch
    std::vector<int> batchscaleindex; // if !allbatches, index into scale list for this batch
    std::vector<int> scalebatchindex; // if !allbatches, index into batch list for this scale

    static int navgscale;   //S  Number of points in moving average, 3, 4 or 5
    static double sdwz;     //S "SD" for weighting
  };
  //--------------------------------------------------------------
  class RelativeBfactor
  {
    // 1) Smooth B-factor on time (or equivalent, usually phi as a proxy)
    // To calculate the Bfactor for an observation at "time" time, uses the
    // class SmoothedValue
    //
    // We have have nbfac factors at intervals of bfacspacing, including one extra
    // on each end of actual range (unless nbfac == 1)
    //   at positions z(i) = -0.5, +0.5, 1.5, ... (nbfac-1.5) in the "normalised"
    //  coordinate 
    // Normalised coordinate z = (time - time0)/scalespacing
    // Then B  = Sum(i) [ w(i) B(i) ] / Sum(i)[w(i)]  for i=0,nbfac-1
    //      gB = exp [ 2 s2 B ] where s2 = (sin theta/lambda)^2  
    // and w(i) depends on the distance between the coordinate of the observation z &
    // the coordinates for each B factor z(i)
    //        d(i) = (z - z(i)) / SDwt    where SDwt is a "SD" factor
    // Use the 3, 4 or 5 points closest to z
    //        w(i) = exp(-d(i)^2)
    //
    // Partial derivatives:
    //   dg/B(i) = 2 s2 gB w(i) / Sum(i)[w(i)]
    //   
    // Smoothing is controlled by two parameters, which are set globally (static)
    //  i)  NavgBfac  number of points in moving average, == 3, 4 or 5 [default 3]
    //      The actual value used will be forced to be <= nbfac
    //  ii) SDwt  "sigma" value for Gaussian weight: default set depending on
    //      NavgBfac
    //
    // 2) Batch scaling
    //   g(j) = B(j)
    //   dg(j)/B(i) = 2 s2 gB  if j == i, else  = 0
    //
  public:
    // Type of B-factor:
    //   NONE     none
    //   BATCH    for each batch
    //   SMOOTH   at intervals, smoothed with Gaussian interpolation
    enum RelativeBtype {NONE, BATCH, SMOOTH, DECAY};

    RelativeBfactor();  // null constructor

    // Construct smooth B-factors from number of Bfactor intervals
    // Note that actual number nbfac will be NbfacIntervals+2
    //   unless NbfacIntervals = 0, in which nbfac = 1
    RelativeBfactor(const int& NbfacIntervals,
		    const Range& timerange);

    // Construct smooth B-factors from spacing 
    //  spacing will be adjusted to give an integral number of intervals
    RelativeBfactor(const double& bfacSpacing, const Range& timerange);

    // Construct batch B-factors from  list of batch numbers (from run)
    RelativeBfactor(const std::vector<int>& BatchNumbers);

    // Return list of ties: sdtie is sd for weight, idx0 is index to first global
    // parameter for setting ties, since they refer to the global parameter index
    std::vector<Tie> Ties(const double& sdtie, const int& idx0) const;

    // Return list of ties: sdtie is sd for weight, idx0 is index to first global
    // parameter for setting ties, since they refer to the global parameter index
    // Tie B-factor to zero
    std::vector<Tie> ZeroTies(const double& sdtie, const int& idx0) const;

    // Store B-factor vector
    void StoreBfactors(const std::vector<double>& Bfacs);
    // Store Nobs vector (length nbfac)
    void StoreNobservations(const std::vector<int>& Nobs);
    // Retrieve B-factors
    std::vector<double> Bfactors() const;
    // Retrieve number of contributions (length nbfac)
    std::vector<int> Nobservations() const;

    // Number of B-factors, == 0 if null
    int Number() const {return nbfac;}
    // Interval
    double Spacing() const {return bfacspacing;}

    // Smooth B-factors:
    // B-factor at position time for 4(sin theta/lambda)^2 = invresolsq
    double BfactorScale(const double& time, const double& invresolsq) const;
    // Return B-factor and derivative vector at position time
    void BfactorScaleDeriv(const double& time, const double& invresolsq,
		      double& bfac, std::vector<double>& dgdp) const;
    // Return B-factor and derivative vector at position time
    double BfactorScaleDeriv(const double& time, const double& invresolsq,
		       std::vector<double>& dgdp) const;
    // Bfactor value at time "time"
    double BfactorValue(const double& time);

    // Batch B-factors:
    // B-factor at batch number batch
    double BfactorScale(const int& batch, const double& invresolsq) const;
    // Return B-factor scale and derivative vector at batch number batch
    void BfactorScaleDeriv(const int& batch, const double& invresolsq,
		      double& gbfac, std::vector<double>& dgdp) const;
    // Return B-factor scale and derivative vector at batch number batch
    double BfactorScaleDeriv(const int& batch, const double& invresolsq,
			    std::vector<double>& dgdp) const;
    // Bfactor value for batch
    double BfactorValueB(const int& batch) const;

    // true if Batch mode
    bool IsBatchBfactor() const {return batchbfac;}

    // valid B-factor at batch number batchnum
    bool ValidBfactor(const int& batchnum) const;

    // index number in this batch list for batch number batchnum, omitting rejects
    int batchSerialIndex(const int& batchnum) const;

    // set ib'th batch if Usebatch[ib] false, batchnumbers corresponding
    void setBatchReject(const std::vector<bool>& Usebatch,
			const std::vector<int>& batchnumbers);

    std::string format() const;

    // Format all information into a labelled save format for later restoration
    std::string FormatSave() const;

    // restore
    void Restore(Fileread& FR);

    static void SetNavgBfac(const int& NavgBfac);
    static int NavgBfac() {return navgbfac;}
    static void SetSDwt(const double& SDwt) {sdwt = SDwt;}
    static double SDwt() {return sdwt;}

  private:
    void init(const Range& timerange);
    void BfactorScaleDeriv(const bool& Deriv, const double& time, const double& invresolsq,
		    double& gbfac, std::vector<double>& dgdp) const;
    void BfactorScaleDeriv(const bool& Deriv,
		      const int& batch, const double& invresolsq,
		      double& gbfac, std::vector<double>& dgdp) const;

    bool batchbfac;     // batch | smooth
    int nbfacintervals;
    int nbfac;
    double bfacspacing;
    double time0;
    hash_table batch_lookup;  // batch lookup for batch Bfactors

    SmoothedValue smoothB;        // smooth Bfactors
    std::vector<double> bfactors;  // batch B-factors
    std::vector<int> nobsPar;  // number of observations for each Bfactors
    bool allbatches;            // true if all batches are used
    std::vector<bool> usebatch; // flag to use this batch
    std::vector<int> batchbfacindex; // if !allbatches, index into scale list for this batch
    std::vector<int> bfacbatchindex; // if !allbatches, index into batch list for this scale

    static int navgbfac;  //  Number of points in moving average, 3, 4 or 5
    static double sdwt;    // "SD" for weighting
  };
  //--------------------------------------------------------------
  class SecondaryScale
  {
    // Secondary beam (absorption) correction
    // Parameterised as sum of spherical harmonics

  public:
    enum SecondaryScaleType {NONE, SECONDARY, ABSORPTION};

    SecondaryScale() : secscltype(NONE), lmax(0) {}
    SecondaryScale(const SecondaryScaleType& secSclType,
		   const int& lmax, const int& lmaxodd,
		   const int& pole);

    // Return list of ties: sdtie is sd for weight, idx0 is index to first global
    // parameter for setting ties, since they refer to the global parameter index
    std::vector<Tie> Ties(const double& sdtie, const int& idx0);

    // Coordinate type
    SecondaryScaleType Type() const {return secscltype;}

    // Store coefficient vector (length ncoeffs)
    void StoreCoefficients(const std::vector<double>& Sphcoefficients);
    // Store Nobs vector (length ncoeffs)
    void StoreNobservations(const std::vector<int>& Nobs);
    // Retrieve coefficient vector (length ncoeffs)
    std::vector<double> Coefficients() const {return sphcoefficients;}
    // Retrieve number of contributions (length ncoeffs)
    std::vector<int> Nobservations() const {return nobsPar;}

    int Number() const {return ncoeffs;}

    // Return scale for secondary beam directions polar angles thetap, phip
    double Scale(const double& thetap, const double& phip) const;
    // Return scale & derivatives for secondary beam directions polar angles thetap, phip
    void ScaleDeriv(const double& thetap, const double& phip,
		    double& scale, std::vector<double>& dgdp) const;
    // Return scale & derivatives for secondary beam directions polar angles thetap, phip
    double ScaleDeriv(const double& thetap, const double& phip,
		    std::vector<double>& dgdp) const;

    std::string format() const;

    // Return h, k, l for pole = 1,2,3, else "none",
    std::string formatPole() const;

    //!
    void Check() const {sphHarmonic.Check();} //!

    // Format all information into a labelled save format for later restoration
    std::string FormatSave() const;

    // restore
    void Restore(Fileread& FR);

  private:
    SecondaryScaleType secscltype;
    SphericalHarmonic  sphHarmonic;  // spherical harmonic function
    int lmax;  // default initialised to 0
    int lmaxodd;  // maximum odd order < lmax
    int pole;    // for ABSORPTION, crystal mode, = 1,2, or 3 for h,k,l
    int ncoeffs; // Number of coefficents
    std::vector<double> sphcoefficients;  // coefficients for each spherical harmonic
    std::vector<int> nobsPar;  //  number of observations for each coefficient

  };
  //--------------------------------------------------------------
  class ScaleSpecification
  // The specification from one SCALES command (in case of multiple runs)
  {
  public:
    ScaleSpecification() : run(-1), isdefault(true),
			   batch(false), nscales(-1), spacing(5.0),
			   nbfac(-1), bspacing(20.0),
			   sec_abs(scala::SecondaryScale::SECONDARY),
			   lmax(4), lmaxodd(3), pole(-1),
			   ntilex(-1), ntiley(-1),
			   detectorscaletype(DetectorScale::NONE) {}

    enum ParameterSDusage {NONE, DIAGONAL, COVARIANCE};

    void dump() const;

    void SetConstant(const int& irun=-1); // SCALES CONSTANT

    int run;     // Run number for this specification, = -1 for all runs
    bool isdefault; // true if this is the default, ie not explicit  

    bool batch;  // true for batch mode
    int nscales; // Number of scales, = -1 for spacing specified
    float spacing; // ROTATION SPACING

    int nbfac;   // number of Bfactors, = 0 OFF, = -1 spacing specified
    float bspacing; // BROTATION SPACING

    SecondaryScale::SecondaryScaleType sec_abs; // NONE, SECONDARY, ABSORPTION
    int lmax;    // Order for secondary|absorption correction (must be even)
    int lmaxodd; //   maximum order for odd terms, < lmax
    int pole;  // for ABSORPTION, = 1,2,3 for h,k,l, = -1 unspecified, = 0 SECONDARY
    
    int ntilex, ntiley;  // number of tiles, ntilex < 0 for no correction
    DetectorScale::DetectorScaleType detectorscaletype;
  };
  //--------------------------------------------------------------
  class WavelengthChebyshevScale
  // Chebyshev polynomial wavelength normalization for Laue data
  // Scale w(lambda) = f(lambda) / f(lambda_ref) where
  //   f(lambda) = sum_{k=0}^{degree} c_k T_k(z(lambda))
  //   z(lambda) = (2*lambda - lam_min - lam_max) / (lam_max - lam_min)
  // Supports up to MAXRANGES wavelength ranges, each with its own polynomial
  {
  public:
    static const int MAXRANGES = 5;

    struct WavelengthRange {
      double lam_min, lam_max;
      int degree;
      int offset;  // offset into global coefficients vector
    };

    WavelengthChebyshevScale() : ncoeffs(0), lambda_ref(0.0), ref_range(-1) {}

    // Construct from list of ranges and reference wavelength
    WavelengthChebyshevScale(const std::vector<WavelengthRange>& Ranges,
                              const double& LambdaRef);

    // Number of Chebyshev coefficients (sum of degree+1 over all ranges)
    int Number() const {return ncoeffs;}
    // Number of ranges
    int Nranges() const {return int(ranges.size());}

    // Return reference wavelength
    double LambdaRef() const {return lambda_ref;}

    // Store all coefficients (length ncoeffs)
    void StoreCoefficients(const std::vector<double>& Coeffs);
    // Retrieve all coefficients
    std::vector<double> Coefficients() const {return coeffs;}

    // Return normalization scale for wavelength lambda
    // Returns 1.0 if lambda is out of all ranges
    double Scale(const double& lambda) const;

    // Return scale and derivatives d(scale)/d(c_j) for all j
    // Returns 1.0 and zero derivatives if lambda is out of range
    double ScaleDeriv(const double& lambda, std::vector<double>& dgdp) const;

    // Return true if any ranges are defined
    bool IsActive() const {return (ncoeffs > 0);}

    // Return range struct for range index ir
    WavelengthRange Range(const int& ir) const {return ranges.at(ir);}

    // Return true if parameter index idx_in_wav (0-based within wavelength block)
    // is the c[0] (constant term) of any range — these need a positive lower bound
    bool IsConstantTerm(const int& idx_in_wav) const {
      for (int ir = 0; ir < int(ranges.size()); ++ir) {
        if (idx_in_wav == ranges[ir].offset) return true;
      }
      return false;
    }

    // Format normalization table for log output: w(lambda) at npoints per range
    std::string PrintNormalization(const int& npoints = 10) const;

    // XML representation of the fit: ranges, coefficients, sampled w(lambda)
    std::string asXML() const;

    // format for save/restore
    std::string FormatSave() const;
    void Restore(Fileread& FR);

  private:
    // Clenshaw evaluation: sum_{k=0}^{n} c[k] T_k(z)
    double chebeval(const std::vector<double>& c, const double& z) const;

    // ASCII line plot of w(lambda) over range irange (for log output)
    std::string AsciiPlot(const int& irange, const int& width = 60,
                          const int& height = 15) const;

    // Evaluate Chebyshev basis T_k(z) for k=0..degree into T
    void chebbasis(const int& degree, const double& z,
                   std::vector<double>& T) const;

    // Return range index for wavelength, -1 if outside all ranges
    int rangeIndex(const double& lambda) const;

    // Map wavelength to z in [-1,1] for given range index
    double mapToZ(const double& lambda, const int& irange) const;

    // Evaluate f(lambda) using range irange
    double evalRange(const double& lambda, const int& irange) const;

    std::vector<WavelengthRange> ranges;
    std::vector<double> coeffs;   // all coefficients, concatenated across ranges
    int ncoeffs;
    double lambda_ref;
    int ref_range;  // range index for lambda_ref (-1 if not set)
    double f_ref;   // f(lambda_ref), cached after StoreCoefficients
  };
  //--------------------------------------------------------------
  class WavelengthGPRScale
  // Gaussian-process wavelength normalization for Laue data.
  //
  // Non-parametric alternative to WavelengthChebyshevScale. Unlike the
  // Chebyshev model this is NOT refined inside the BFGS scaling loop: it is
  // fitted once as a pre-pass and applied thereafter as a fixed multiplicative
  // correction (a precomputed lookup table), so it contributes no parameters to
  // the scale-model parameter vector and has zero derivatives.
  //
  // Fit procedure (see Fit()):
  //   1) bin per-observation log-ratios  y = log(I_obs / <I>_symmetry)  by wavelength
  //   2) form a heteroscedastic training set (bin centre, bin-mean, SEM^2)
  //   3) fit a zero-mean Gaussian process in LOG space with a squared-exponential
  //      (or Matern-3/2) kernel; the length scale is chosen by maximising the
  //      log marginal likelihood (unless fixed by the user)
  //   4) evaluate the posterior mean g(lambda) on a dense uniform grid (lookup table)
  // Scale:  ws(lambda) = exp(g(lambda) - g(lambda_ref))   (>0 by construction)
  {
  public:
    enum KernelType {SQEXP, MATERN32};

    // User-supplied controls (from LAUE NORMGPR keyword)
    struct GPRControl {
      bool enabled;
      double lam_min, lam_max; // wavelength range; <=0 => take limits from data
      double lengthscale;      // GP length scale; <=0 => optimise by marginal likelihood
      int nbins;               // number of wavelength bins for the training set
      KernelType kernel;
      GPRControl() : enabled(false), lam_min(0.0), lam_max(0.0),
                     lengthscale(0.0), nbins(50), kernel(SQEXP) {}
    };

    WavelengthGPRScale() : active(false), lam_min(0.0), lam_max(0.0),
                           lambda_ref(0.0), g_ref(0.0), grid_step(0.0),
                           ngrid(0), kernel(SQEXP), lengthscale(0.0),
                           sigf(0.0), signoise(0.0), ntrain(0) {}

    // Fit the GP from per-observation samples; returns number of training bins
    // used (0 => fit failed / inactive). Diagnostic text appended to fitlog.
    int Fit(const std::vector<double>& lambdas,
            const std::vector<double>& logratios,
            const std::vector<double>& weights,
            const GPRControl& ctrl, const double& LambdaRef,
            std::string& fitlog);

    bool IsActive() const {return active;}

    double LambdaRef() const {return lambda_ref;}

    // Normalization scale ws(lambda); returns 1.0 if inactive or out of range
    double Scale(const double& lambda) const;

    // Format normalization table for log output
    std::string PrintNormalization(const int& npoints = 12) const;

    // XML representation of the fit: range, hyperparameters, sampled w(lambda)
    std::string asXML() const;

    // format for save/restore
    std::string FormatSave() const;
    void Restore(Fileread& FR);

  private:
    double kernelValue(const double& la, const double& lb) const;

    // ASCII line plot of w(lambda) over the fitted range (for log output)
    std::string AsciiPlot(const int& width = 60, const int& height = 15) const;

    bool active;
    double lam_min, lam_max;   // fitted wavelength range
    double lambda_ref, g_ref;  // reference wavelength and g(lambda_ref)
    double grid_step;          // uniform spacing of the lookup grid
    int ngrid;                 // number of lookup points
    std::vector<double> grid_g; // posterior-mean log-scale on the grid
    KernelType kernel;         // kernel used (for reporting)
    double lengthscale;        // fitted/used length scale (for reporting)
    double sigf;               // signal sdev (for reporting)
    double signoise;           // mean training noise sdev (for reporting)
    int ntrain;                // number of training bins used (for reporting)
  };  // class WavelengthGPRScale
  //--------------------------------------------------------------
  class LinkSpecs {
    // Links for SURFACE parameters
  public:

    LinkSpecs() : nlinks(0), linkall(false), nunlinks(0) {}

    // true is any [un]links are set
    bool isSet() const {return ((nlinks != 0) || (nunlinks != 0));}

    void setLinkAll() {linkall = true; nlinks = -1; links.clear();}
    void setUnlinkAll() {unlinkall = true; nunlinks = -1; unlinks.clear();}

    bool linkAll() const {return linkall;}
    bool unlinkAll() const {return unlinkall;}

    void addLink(const std::pair<int,int> runs2) {
      links.push_back(runs2); nlinks = links.size()+1;
    }
    void addUnlink(const std::pair<int,int> runs2) {
      unlinks.push_back(runs2); nunlinks = unlinks.size()+1;
    }

    int Nlinks() const {return nlinks;}
    int Nunlinks() const {return nunlinks;}

    std::vector<std::pair<int,int> > Links() const {return links;}
    std::vector<std::pair<int,int> > Unlinks() const {return unlinks;}

  private:
    // nlinks = 0, no links, = -1 link all
    int nlinks;
    std::vector<std::pair<int,int> > links;    // for each pair, link run(first) -> run(second)
    bool linkall;
    // nunlinks = 0, no unlinks, = -1 unlink all
    int nunlinks;
    std::vector<std::pair<int,int> > unlinks;  // for each pair, unlink run(first) -> run(second)
    bool unlinkall;
  };
}

#endif
