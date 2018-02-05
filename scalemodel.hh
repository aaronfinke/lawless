// scalemodel.hh
//
//  Scale model for reflection list
//
//  The inverse scale factor ghl for a given observation (or part) is the product of
//  number of parts
//
//    ghl = g(primary_beam_direction) g(time) g(secondary_beam_direction) ...
//
//  For each run, there is a primary beam component, and an index into lists for
//  each of the other components, so that multiple runs may share these (eg secondary)
//  

#ifndef SCALEMODEL_HEADER
#define SCALEMODEL_HEADER

#include <vector>
#include "InputAll.hh"
#include "hkl_unmerge.hh"
#include "scaletypes.hh"
#include "Output.hh"
#include "tie.hh"

#include "Minimizer.h"
#include "fmat.h"


namespace scala {
  //--------------------------------------------------------------
  class ValidScaleModel {
    // For each run, does the data have suffcient information to support each
    // type of scale model?
  public:
    ValidScaleModel(): nruns(0) {}
    ValidScaleModel(hkl_unmerge_list& hkl_list);
    void init(hkl_unmerge_list& hkl_list);

    //! validity flag for primary beam corrections, for irun'th run
    bool ValidPrimary(const int& irun) const
    {return validprimary.at(irun);}

    //! validity flag for batch corrections, for irun'th run
    bool ValidBatch(const int& irun) const
    {return validbatch.at(irun);}

    //! validity flag for secondary beam corrections, for irun'th run
    bool ValidSecondary(const int& irun) const
    {return validsecondary.at(irun);}

    //! validity flag for tile corrections, for irun'th run
    bool ValidTile(const int& irun) const
    {return validtile.at(irun);}

  private:
    // Validity flags for each run
    int nruns;
    std::vector<bool> validprimary;   // smooth primary scale (& B-factor)
    std::vector<bool> validbatch;     // batch number
    std::vector<bool> validsecondary; // secondary beam, ie full geometry
    std::vector<bool> validtile;      // tile information
  }; // class ValidScaleModel
  //--------------------------------------------------------------
  class ScaleModel 
  {
  public:
    // types of parameter
    enum ScaleParameterType {NONE, SCALE, BFACTOR, SECONDARY, TILE};
    // String representation
    static std::string ScaleParameterTypeString(const ScaleParameterType& type);

    ScaleModel():status(0){}
    // Construct from input commands and reflection list
    //  All secondary beam direction in hkl_list will be calculated if needed
    ScaleModel(const phaser_io::InputAll& input,
	       hkl_unmerge_list& hkl_list,
	       phaser_io::Output& output);

    // Initialise from input commands and reflection list
    //  All secondary beam direction in hkl_list will be calculated if needed
    void init(const phaser_io::InputAll& input,
	      hkl_unmerge_list& hkl_list,
	      phaser_io::Output& output);

    // call at beginning of each refinement cycle to clear negative sec scale flag
    void clearNegativeSecScaleFlag();

    // return > 0 if a negative secondary scale occurred, and clear the flag
    int NegativeSecScaleOccurred();

    // return > 0 if a negative secondary scale occurred on last (previous) cycle
    int NegativeSecScaleLast() {return negativeSecScaleLast;}

    void SetConstant(hkl_unmerge_list& hkl_list,
		     phaser_io::Output& output);  // set SCALE CONSTANT for all runs

    // If true, allow tile corrections to vary azimuthally
    // if false, force to be radially symmetric
    // returns true if anything has changed
    bool symmetricTiles(const bool& symmetric);

    // switch secondary scales On (true) or Off (false)
    // returns true if anything has changed
    bool switchSecondaryScales(const bool& on);

    // Set reject list for batches, relevant for BATCH scale mode only (fail if not)
    void setBatchReject(const std::vector<bool>& usebatch,
			const std::vector<int>& batchnumbers);

    // Get vector of parameters
    std::vector<double> GetParameters() const;
    // get type for all parameters
    std::vector<ScaleParameterType> GetParameterType() const;
    // get type for a parameter
    ScaleParameterType GetParameterType(const int& Ipar) const;

    //! return detector scale number & scale index into detector parameter list
    std::pair<int,int> DetectorParameterNumber(const int& Ipar) const;

    // get lower bound for parameter, depending on type: return false if unbounded
    bool GetLowerBound(const int& Ipar, double& Lower) const;
    // get upper bound for parameter, depending on type: return false if unbounded
    bool GetUpperBound(const int& Ipar, double& Upper) const;
    // get "large shift" value for parameter, depending on type
    double GetLargeShift(const int& Ipar) const;

    // Clear observation counts at beginning of cycle, mainly relevant for tiles
    void clearCounts();

    // Set all parameters from vector and count of number of contributions
    void SetParameters(const std::vector<float>& params, const std::vector<int>& Nobs,
		       const bool& donormalise=true);
    void SetParameters(const std::vector<double>& params, const std::vector<int>& Nobs,
		       const bool& donormalise=true);

    // Normalise scales & B-factors
    void NormaliseParameters();

    // Fix up secondary scales if there is a negative scale for any observation
    void fixupSecondaryScales();

    // Set initial primary scales, eg from InitialScales
    // also store count of number of observations
    void SetInitialScales(const std::vector<double>& gscales,
			  const std::vector<int>& numobsrotrange);

    // Return true if model is refinable, ie not just one scale and one B-factor
    bool IsRefinable() const;
    // Total number of parameters
    int Nparameters() const {return nparameters;}
    //  Number of primary scale parameters
    int Nscales() const {return nprimaryscale;}
    //  Number of Bfactor parameters
    int NBfactors() const {return nbfactors;}
    // Number of secondary scale parameters
    int Nsecondary() const {return nsecondaryscale;}
    // Number of normalisation parameters: usually 2 (k,B);
    //   1 if fixed Bfactors for any run; = 0 if not refineable
    int Nfilter() const;

    // Return true if BATCH scales for all runs
    bool isAllBatch() const;

    // return vector of batch scale factors, empty if not batch
    std::vector<double> getBatchScales() const;

    // Return primary scale for specified run
    PrimaryScale primary_scale(const int& irun) const {return primary_scales.at(irun);}
    // Return Bfactor for specified run
    RelativeBfactor Bfactor(const int& irun) const {return  relative_bfactors.at(irun);}

    // return number of ranges in run, for batch scale = Nbatches, else number of scales-1
    int Nrange(const int& irun) const {return primary_scales.at(irun).Nintervals();}

    // Scale observation, returns scale applied
    double ScaleObs(observation& obs, const Rtype& invresolsq,
		    const bool& onlyUseSingletons=true) const;
      
    // Scale observation, returns scale applied and
    // partial derivative vector d(ghl)/dp
    double ScaleObs(observation& obs, const Rtype& invresolsq,
		   std::vector<double>& dghldp) const;

    // Return restraint target, and optionally gradient & Hessian contributions
    double TieValues(const bool& DoGradient, const bool& DoHessian,
		     const std::vector<double>& params,
		     std::vector<double>& dRdpi,
		     std::vector<TieHessian>& Htie);

    // Print scale layout
    void PrintLayout(phaser_io::Output& output) const;
    // Print all scale parameters
    void PrintScales(phaser_io::Output& output) const;
    // Print secondary corrections
    void PrintSecondaryCorrections(phaser_io::Output& output) const;
    

    //! Write image[s] for each detector scale
    void WriteImage(const std::string imagefilename,
		    phaser_io::Output& output) const;

    //!
    void Check() const {if (nsecscales > 0) secondary_scales[0].Check();}

    //! Dump scale model to file
    //    void Save(const std::string& dumpfilename,
    //	      const std::vector<Run>& runlist) const;

    //! Format scalemodel for save/restore
    // NB ties are not saved
    std::string FormatSave(const std::vector<Run>& runlist) const;

    //! Restore from file
    void Restore(const std::string& restorefilename,
		 const std::vector<Run>& runlist);

    //! store parameter variance information
    void storeParameterVariances(const TNT::Fortran_Matrix<floatType>& H,
				 const double& wd2in,
				 const int& nminusm);
    
    //! turn off parameter variances
    void ignoreParameterVariances() {nfreedom = 0;}

    //! true if we have parameter variances swtiched on
    bool haveParameterVariances() const {return (nfreedom > 0);}

    //! Usage of parameter variances, = {NONE, DIAGONAL, COVARIANCE};
    scala::ScaleSpecification::ParameterSDusage parameterSDusage() const
    {return parametersdusage;}

    // return false if number of variance parameters is not same as nparameters
    // OK (true) if no variance used
    bool checkVarianceNumbers() const;

  private:
    // Setup from scale specifications and reflection list
    // Sets pole, for ABSORPTION, = 1,2,3 for h,k,l, = -1 unspecified, = 0 SECONDARY
    void setup(const std::vector<ScaleSpecification>& scaleSpecs,
	       const LinkSpecs& linkspecs,
	       hkl_unmerge_list& hkl_list,
	       phaser_io::Output& output);

    std::string SetupScale(const int& irun,
			   const ScaleSpecification& scaleSpec,
			   const Run& run,
			   const ValidScaleModel&  validscalemodel,
			   phaser_io::Output& output);

    void CountParameters();

    // Returns index into scale specification list for run index irun
    //  returns -1 if not found
    int scaleSpecIndex(const int& irun,
		       const std::vector<ScaleSpecification>& scaleSpecs,
		       const std::vector<Run>& runList) const;

    // Status:
    //  = 0  unset
    //  > 0  set and OK
    //  < 0  insufficient information (ie not IsRefinable)
    int status;

    // Run information
    std::vector<int> runnumbers;
    // Index to Run number from lattice number, for multilattice data
    std::vector<int> idxrunlattice;

    // Primary beam things
    int nruns;  // number of runs == number of primary models
    std::vector<PrimaryScale> primary_scales;        // for each run
    // index to 1st primary scale parameter for each run
    std::vector<int> idxrun_primary_scales;
    std::vector<RelativeBfactor> relative_bfactors;  // for each run
    // index to 1st B-factor parameter for each run
    std::vector<int> idxrun_bfactors;

    // Secondary beam things
    std::vector<SecondaryScale> secondary_scales;  // secondary models
    int nsecscales;                                // number of secondary models <= nruns
    std::vector<int> sec_scale_index_run;  // index into secondary scale list for each run
					   // set up by [UN]LINK or by default
    // index to 1st secondary parameter for each run
    std::vector<int> idxrun_secondary;
    // for ABSORPTION, pole = 1,2,3 for h,k,l, = -1 unspecified, = 0 SECONDARY
    int pole;
    // initially = 0 on each cycle (call to clearNegativeSecScaleFlag),
    // count > 0 if secondary scale for observation is negative
    mutable int negativeSecScale;  // within each refine cycle
    mutable int negativeSecScaleLast;  // on previous (last) cycle
    // ... or at any time 
    int negativeSecScaleOccurred;

    // Scaling by tile (or other detector scale)
    //   typically only one scale set, unless different runs are from different detectors
    std::vector<DetectorType> detectortypes;   // for each run
    std::vector<DetectorScale> detector_scales;    // the tile scales
    std::vector<int> detector_scale_index_run; // which scale for each run?
    int ndetscales;  // number of detector scales
    // index to 1st detector parameter for each run
    std::vector<int> idxrun_detector;

    // Ties
    int nties;
    std::vector<Tie> ties;
    // counts for each type
    int nties_rot, nties_bfac, nties_zerob, nties_surf, nties_tiles;

    double sd_rotation;
    double sd_bfactor;
    double sd_zerob;
    double sd_surface;

    // For CCD tiles:
    //  tie_tile[0] for r
    //  tie_tile[1] for w
    //  tie_tile[2] for A
    //  tie_tile[3] for x0, y0
    //  tie_tile[4] for Fourier coefficients
    //  tie_tile[5] for r target
    //  tie_tile[6] for w target
    std::vector<double> tie_tile; // other types, no ties

    int nparameters;      // Number of parameters
    int nprimaryscale;    //  Number of primary scale parameters
			  //  maybe = 0 if one CONSTANT run
    int nbfactors;        //  Number of B-factor parameters
    int nsecondaryscale;  //  Number of secondary scale parameters
    int ntilescale;       //  Number of tile scale parameters

    // Normalisation:
    int scalenormrun;    // run for scale  normalisation
    int scalenormbatch;  // batch number for scale  normalisation, -1 for 1st
			 //  after construction, batch serial number in run

    int bfacnormrun;     // run number for B-factor normalisation
    int bfacnormbatch;   // batch number for B-factor  normalisation, -1 for best
			 //  after construction, batch serial number in run
    bool normalisebfac;  // usually true to normalise Bfactor,
			 // false if no Bfactor refinement for at least one run

    // Variance/covariance information
    clipper::Array2d<double> VC;  // variance/covariance matrix for parameters
    std::vector<double> varpar;   // parameter variances from VC diagonal
    double wd2;                   // Sum(w Del^2)
    // degrees of freedom (n-m), initially = -1, = 0 if not set yet
    int nfreedom;

    //  {NONE, DIAGONAL, COVARIANCE};
    //  != NONE to use parameter SDs in sig(I) scaling
    scala::ScaleSpecification::ParameterSDusage parametersdusage; 

    // Print wrapping lines:
    //   line 1, values v (double), label t1
    //   line 2, values v2 (int),   label t2  (optional, if size > 0)
    //   line 3, values n (int),    label t3
    //   fw  field width
    //   fd  number of decimal points for v
    static std::string PrintWrappingLines(const std::vector<double>& v,
				   const std::string& t1,
				   const std::vector<double>& v2,
				   const std::string& t2,
				   const std::vector<int>& n,
				   const std::string& t3,
				   const int& fw,
				   const int& fd);

    static std::string PrintWrappingLinesWithSD(const std::vector<double>& v,
					 const std::string& t1,
					 const std::vector<double>& sds,
					 const int& fw,
					 const int& fd);

    // extract nsd SDs from variance, beginning at index idxsd
    std::vector<double> extractSDs(const int& idxsd, const size_t& nsd) const;



    // return index in runlist, -1 if not found
    int RunNotFound(const std::vector<Run>& runlist,
		    const std::vector<int>& batchnumbers) const;

    // Returns scale for observation, using scale set jscale (== irun for main observation
    double ScaleFactor(const int& jscale,
		       const observation& obs, const Rtype& invresolsq) const;

    // Returns scale for observation, using scale set jscale (== irun for main observation
    // and partial derivative vector d(ghl)/dp
    double ScaleFactorDeriv(const int& jscale,
			    observation& obs,
			    const Rtype& invresolsq,
			    std::vector<double>& dghldp) const;

    void SetupTies(const phaser_io::InputAll& input);
    void SetupTies();

    // set automatic TILE settings if appropriate
    // modifies scaleSpecs
    void autoTiles(std::vector<scala::ScaleSpecification>& scaleSpecs,
		   hkl_unmerge_list& hkl_list,
		   phaser_io::Output& output);

    // Variance(gscale)
    double VarScale(const std::vector<double>& dghldp) const;

    void processLinks(const LinkSpecs& linkspecs,
		      const std::vector<Run>& runlist,
		      phaser_io::Output& output);

    // returns run indices for both ends of the link, first = -1 if not found
    std::pair<int,int> checkLink(const std::pair<int,int>& link,
				 const std::vector<Run>& runlist) const;

    // list of run indices which share scale index sclidx
    std::vector<int> runsWithScaleIndex(const int& sclidx) const;

    // get number of secondary scales from sec_scale_index_run, check that all
    // are present
    int numberSecondaryScale() const;

    // format as "Set <scaleset>, run[s]: <runnumbers>"
    std::string formatSecondaryrunset(const int& scaleset) const;

  }; // class ScaleModel 
}

#endif
