// sdmodel.hh

// List of SDcorrection models, for each run

#ifndef SDMODEL_HEADER
#define SDMODEL_HEADER

#include "hkl_datatypes.hh"
#include "sdctypes.hh"
#include "hkl_unmerge.hh"
#include "selectedobservations.hh"
#include "InputAll.hh"

namespace scala
{
//--------------------------------------------------------------
  class SDties {
  public:
    // Ties, same for all runs
    // tietype = 0 no tie, = -1 defaults, = +1 set from input,
    //         = +2 similarity tie ie target is average value over runs (nruns >1)
    int tietype;
    std::vector<double> targets;   // 3 targets
    std::vector<double> sdtargets;  // ... and their SDs (= 0 no target)
  };
//--------------------------------------------------------------
  class SDmodel
  // A set of SDcorrections for each run, and for partials & fulls separately
  // Optionally, the same values may be used for all runs, but separate
  // values are still stored for each run.
  // Also, for runs with few or no fulls or partials, values for the other
  // full/partial class may be used
  {
  public:
    SDmodel() {init();}

    void init();

    //! Set flag to refine or not
    // +1 refine, 0 no refine, -1 linear fit
    void SetRefine(const int& refineFlag) {refine = refineFlag;}
    //! return flag to refine or not
    int Refine() const {return refine;}

    void SetDamp(const double& Damp) {damp = Damp;} //!< set damp factor
    double Damp() const {return damp;} //!< return damp factor

    //! Set flag = true to fix SDb = 0.0
    void SetNoSDb(const bool& NoSDb);
    //! return No SDb flag
    bool NoSDb() const {return nosdb;}
    //! Set flag = true to force same parameters for all runs
    void SetAllRunsSame(const bool& flag);
    //! return flag = true if the same parameters are used for all runs
    bool AllRunsSame() const {return allrunssame;}

    //! Reset all minimum & maximum values
    void ResetRange();  

    //! Add SD correction for one runs: runs must be added in sequentially!
    void AddRun(const int& runNum,
		const SDcorrection& SDCfull, const SDcorrection& SDCpartial);

    //! Store ties for all SD corrections
    // tietype = 0 no tie, = -1 defaults, = +1 set from parameters
    void SetTies(const int& tietype,
		 const std::vector<double>& targets,
		 const std::vector<double>& sdtargets);
    //! Store one tie for all SD corrections, for parameter ipar (0-2)
    void ResetTie(const int& ipar,
		  const double& Target,
		  const double& SDtarget);
    //! format tie information
    std::string formatTie() const;

    //! true if there are active ties
    bool restrained() const;

    //! Set similarity targets from current parameters (initialising)
    void SetTargetsFromAverageParameters();

    //! return tie settings
    SDties Ties() const {return ties;}
    //! reset tie settings
    void ResetTies(const SDties& Ties) {ties = Ties;SetTies();}
    //! return true if Similarity ties (tietype = +2)
    bool IsSimilarityTie() const {return (ties.tietype == +2);}

    //! Number of runs (actual number)
    int Nruns() const {return sdc_full_run.size();}
    //! Number of sets, = Nruns or 1 if AllRunsSame
    int Nsets() const {return nsets;}

    int Nparams() const; //!< Number of parameters

    //! Number of parameter groups, depends on allrunsame flag and full/partial usage
    int Ngroups() const;

    //! Set index list for each run
    void SetIdxParam();

    // RunIndex is index from 0
    // 0 <= RunIndex < Nruns
    //! Return SDcorrection for fulls for given run
    SDcorrection SDCfull(const int& RunIndex) const
    {return sdc_full_run.at(RunIndex);}
    //! Return SDcorrection for partials for given run
    SDcorrection SDCpartial(const int& RunIndex) const
    {return sdc_partial_run.at(RunIndex);}

    //! Set SDcorrection for fulls for given run
    void setSDCfull(const int& RunIndex, const SDcorrection& sdc)
    {sdc_full_run.at(RunIndex) = sdc;}
    //! Set SDcorrection for partials for given run
    void setSDCpartial(const int& RunIndex, const SDcorrection& sdc)
    {sdc_partial_run.at(RunIndex) = sdc;}

    int NparamsFull(const int& RunIndex) const //!< Number of parameters for fulls
    {return sdc_full_run.at(RunIndex).Nparams();}
    int NparamsPartial(const int& RunIndex) const //!< Number of parameters for fulls
    {return sdc_partial_run.at(RunIndex).Nparams();}

    int NparamsPerGroup() const {return nparamspergroup;}

    //! use partial correction for fulls 
    // if few == true, there are few fulls, if false none
    void FullAsPartial(const int& RunIndex, const bool& few);

    //! use full correction for partials 
    // if few == true, there are few partials, if false none
    void PartialAsFull(const int& RunIndex, const bool& few);

    //! return true if at least one run has few fulls||partials
    bool SomeFewRuns() const;

    //! Update SDC for run RunIndex by multiplying SDfac values
    //    UpdateFull     update factor for fulls
    //    UpdatePartial  update factor for partials
    void UpdateFactor(const int& RunIndex,const bool& copy,
		      const float& UpdateFull, const float& UpdatePartial);

    //! Reset SDadd for all runs
    void SetSDadd(const double& SDadd);

    //! Apply appropriate SD correction to all valid observations in reflection
    //! Return vector of uncorrected sigmas, for all observations,
    //!  even unselected ones
    std::vector<float> CorrectReflection(reflection& Ref) const;

    //! Apply appropriate SD correction to all observations in reflection, for ROGUES
    void CorrectAllReflection(reflection& Ref) const;

    //! Return overall minimum & maximum values
    void GetSDcorrectionRanges(float& minSDcorrFulls, float& maxSDcorrFulls,
			       float& minSDcorrPartials, float& maxSDcorrPartials) const;

    // Weighting scheme for averaging observations in SD correction
    //! Set weight
    void SetWeight(const WeightType::AverageWeightType& weightType)
    {weighttype = weightType;}

    //! Set variance weights
    void SetVarianceWeights() {SetWeight(WeightType::VARIANCE);}

    //! Set SqrtScale weights
    void SetSqrtScaleWeights() {SetWeight(WeightType::SQRTSCALE);}

    //! Set Scale weights
    void SetScaleWeights() {SetWeight(WeightType::SCALE);}

    //! Set unitweights
    void SetUnitWeights() {SetWeight(WeightType::UNIT);}

    WeightType::AverageWeightType Weight() const {return weighttype;}

    // = = =  API for refinement
    //! Get vector of parameters
    std::vector<double> GetParameters() const;

    //! Get vector of parameter "shifts" ie starting values for varying
    //!  the parameters in simplex optimisation
    std::vector<double> GetShifts(const double& scale) const;

    //! Set all parameters from vector
    void SetParameters(const std::vector<double>& params);

    //! Set all parameters from vector, parameterupdated true if this parameter has been updated
    void SetParameters(const std::vector<double>& params,
		       const std::vector<bool>& parameterupdated);

    // return vector elements for each observation in Selobs
    // each element is vector of elements for each parameter
    //  elements for each parameter are d(delta(iobs))/dp(k)
    std::vector <std::vector<double> > GetDerivatives
    (SelectedObservations& Selobs,
     const std::vector<float>& sigmaI) const;
    // scale-weighted <I>
    std::vector <std::vector<double> > GetDerivativesscalewt
    (SelectedObservations& Selobs,
     const std::vector<float>& sigmaI) const;

    //! lower bounds for each parameter (0.0 means no bound)
    std::vector<double> LowerBounds() const;
    //! upper bounds for each parameter (0.0 means no bound)
    std::vector<double> UpperBounds() const;
    //! large shifts  for each parameter
    std::vector<double> LargeShifts();

    //! Return restraint R2 for each parameter group
    std::vector<double> GetRestraintR() const;
    //! return derivatives dR2/dp and Hessian for each parameter group
    void GetRestraintDerivatives
    (std::vector <std::vector<double> >& dr2dp,
     std::vector <clipper::Array2d<double> >& H) const;

    //! return the "coordinate" of the iobs'th entry in selobs
    // Also return parameter group index idxgroup
    //   Iav = average I
    std::vector<double> coordinate(const SelectedObservations& selobs,
				   const int& iobs, const double& Iav,
				   int& idxgroup) const;

    //! return the corrected SD (from parameters) for the iobs'th entry in selobs
    //   Iav = average I
    double sdCorrected(const SelectedObservations& selobs,
		       const int& iobs, const double& Iav) const;

    // = = =  

    //! return use flag (only||few + fulls||partials)
    // = 0 normal, both fulls & partials; < 0 treat full as partial; > 0 treat partial as full
    //   = +1 only fulls    = +2 few partials
    //   = -1 only partials = -2 few fulls
    int UseFlag(const int& RunIndex) const {return usetype.at(RunIndex);}

    //! return use flag (only||few + fulls||partials)
    // = 0 normal, both fulls & partials; < 0 treat full as partial; > 0 treat partial as full
    //   = +1 only fulls    = +2 few partials
    //   = -1 only partials = -2 few fulls
    std::vector<int> UseFlags() const {return usetype;}

    //! set true if sample SD is to be used in the final averaging
    void SetSampleSD(const bool& samplesd,
		     const int& Minimumsample)
    {sampleSD = samplesd; minimumsample = Minimumsample;}

    //! return true if sample SD is to be used in the final averaging
    bool SampleSD() const {return sampleSD;}

    //! return minimum sample
    int MinimumSample() const {return minimumsample;}

    //! return formatted only||few + fulls||partials information for all runs
    std::string formatFullPartialInfo() const;

    //! return formatted use flag (only||few + fulls||partials)
    std::string formatUseFlag(const int& RunIndex) const;

    //! return formatted values
    std::string format() const;

    //! as XML
    std::string asXML() const;

    //! return formatted weight information
    std::string formatWeightType() const;

    //! format for Dump
    std::string FormatSave() const;

    void Restore(const std::string& restorefilename,
		 const std::vector<Run>& runlist);

    void dump() const;

  private:
    std::vector<SDcorrection> sdc_full_run;     // for fulls
    std::vector<SDcorrection> sdc_partial_run;  // for partials
    std::vector<int> runnumbers;  // actual run numbers
    std::vector<int> idxfullparam;       // index to first parameter for each run, fulls
    std::vector<int> idxpartialparam;    // index to first parameter for each run, partials
    // index to parameter group for each run, for fulls, partials
    std::vector<std::pair<int,int> > idxparamgroups;
    int nparamspergroup; // number of paramters in each group

    // = 0 normal, both fulls & partials; < 0 treat full as partial; > 0 treat partial as full
    //   = +1 only fulls    = +2 few partials
    //   = -1 only partials = -2 few fulls
    std::vector<int> usetype;

    bool nosdb;  // if true fix SDb = 0
    bool allrunssame;  // if true use same parameters for all runs
    // refine = +1 refine (non-linear), 0 no refine, -1 linear fit
    int refine;
    int nsets;   // number of unique runs, = number of runs or 1 if allrunssame
    bool sampleSD;  // true if sample SD is to be used in the final averaging
    int minimumsample;  // ... with more than this number of observations

    double damp;  // damp factor for refinement

    // Ties, same for all runs
    SDties ties; // just a data structure for save and restore
    //    int tietype;     // = 0 no tie, = -1 defaults, = +1 set from input
    //         = +2 similarity tie ie target is average value over runs (nruns >1)
    //    std::vector<double> targets;   // 3 targets
    //    std::vector<double> sdtargets;  // ... and their SDs (= 0 no target)

    // Weighting scheme for averaging observations in SD correction
    WeightType::AverageWeightType weighttype;   // type of weighting for average

    void SetTies();

    // add 2 or 3 parameters into averages    
    void AddToAverages(std::vector<MeanValue>& averagerealparameters,
		       const std::vector<double>& realparameters) const;

    double ISa(const SDcorrection& sdc) const;
    // ISa = 1/(Sdfac*SDadd)   =~ (I/sig(I))asymtotic for large I
    // see K.Diederichs, Acta Cryst. D66,733

    //! return a reference to the relevant SDC model (run, full/partial
    // Also return parameter group index idxgroup
    const SDcorrection& SDCmodel(const SelectedObservations& selobs,
				 const int& iobs,
				 int& idxgroup) const;

  };
  //--------------------------------------------------------------
  //! Create SDmodel for each run from input or defaults
  SDmodel CreateSDmodel(const phaser_io::InputAll& input,
			const std::vector<Run>& runlist,
			const bool& setnull=false);
  //--------------------------------------------------------------
  SDmodel CreateSDmodel(const std::vector<Run>& runlist,
			const bool& setnull=false);
  // Create SDmodel for each run from  defaults (no input)
  // mainly for testing
//-------------------------------------------------------------
  void SetSdmFullPartialFlags(const Run::FullsAndPartials& FandP,
			      SDmodel& SDM, const int& irun);
}
// ---------------------------------------------------------
class TargetResiduals {
public:
  TargetResiduals() :  quadratic(true), converged(false),
		       R1(0.0), R1lsq(0.0), R2(0.0) {}
  void Add(const double& r1, const double& r2)
  {R1 += r1;R2 += r2;}
  void AddLsq(const double& r1lsq)
  {R1lsq += r1lsq;}
  
  double R()      // == R1 + R2
  {return R1 + R2;}
  
  bool quadratic; // true if LSQ
  bool converged;
  std::vector<bool> updated;   // true if updated, for each class
  //! Just some residuals
  double R1;     // optimised residual, LSQ or ln cosh
  double R1lsq;  // LSQ residual, may be equal to R1
  double R2;     // restraints
};
#endif
