// sdmodel.hh

// List of SDcorrection models, for each run

#ifndef SDMODEL_HEADER
#define SDMODEL_HEADER

#include "hkl_datatypes.hh"
#include "sdctypes.hh"
#include "hkl_unmerge.hh"
#include "selectedobservations.hh"

namespace scala
{
//--------------------------------------------------------------
  class SDmodel
  // A set of SDcorrections for each run, and for partials & fulls separately
  // Optionally, the same values may be used for all runs, but separate
  // vlues are still stored for each run.
  // Also, for runs with few or no fulls or partials, values for the other
  // full/partial class may be used
  {
  public:
    SDmodel() {init();}

    void init();

    //! Set flag to refine or not
    void SetRefine(const bool& refineFlag) {refine = refineFlag;}
    //! return flag to refine or not
    bool Refine() const {return refine;}

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
    //! format tie information
    std::string formatTie() const;

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

    int NparamsFull(const int& RunIndex) const //!< Number of parameters for fulls
    {return sdc_full_run.at(RunIndex).Nparams();}
    int NparamsPartial(const int& RunIndex) const //!< Number of parameters for fulls
    {return sdc_partial_run.at(RunIndex).Nparams();}

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

    //! Apply appropriate SD correction to all valid observations in reflection
    //! Return vector of uncorrected sigmas, for all observations,
    //!  even unselected ones
    std::vector<float> CorrectReflection(reflection& Ref) const;

    //! Return overall minimum & maximum values
    void GetSDcorrectionRanges(float& minSDcorrFulls, float& maxSDcorrFulls,
			       float& minSDcorrPartials, float& maxSDcorrPartials) const;

    // = = =  API for refinement
    //! Get vector of parameters
    std::vector<double> GetParameters() const;

    //! Get vector of parameter "shifts" ie starting values for varying
    //!  the parameters in simplex optimisation
    std::vector<double> GetShifts(const double& scale) const;

    //! Set all parameters from vector
    void SetParameters(const std::vector<double>& params);

    // return vector elements for each observation in Selobs
    // each element is vector of elements for each parameter
    //  elements for each parameter are d(delta(iobs))/dp(k)
    std::vector <std::vector<double> > GetDerivatives
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

    //! return formatted only||few + fulls||partials information for all runs
    std::string formatFullPartialInfo() const;

    //! return formatted use flag (only||few + fulls||partials)
    std::string formatUseFlag(const int& RunIndex) const;

    //! return formatted values
    std::string format() const;

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

    // = 0 normal, both fulls & partials; < 0 treat full as partial; > 0 treat partial as full
    //   = +1 only fulls    = +2 few partials
    //   = -1 only partials = -2 few fulls
    std::vector<int> usetype;

    bool nosdb;  // if true fix SDb = 0
    bool allrunssame;  // if true use same parameters for all runs
    bool refine; // if true refine parameters
    int nsets;   // number of unique runs, = number of runs or 1 if allrunssame

    double damp;  // damp factor for refinement
    // Ties, same for all runs
    int tietype;     // = 0 no tie, = -1 defaults, = +1 set from input
    std::vector<double> targets;   // 3 targets
    std::vector<double> sdtargets;  // ... and their SDs (= 0 no target)

    void SetTies();


  };
  //--------------------------------------------------------------
  //! Create SDmodel for each run from input or defaults
  SDmodel CreateSDmodel(const phaser_io::InputAll& input,
			const std::vector<Run>& runlist);
//-------------------------------------------------------------
  void SetSdmFullPartialFlags(const Run::FullsAndPartials& FandP,
			      SDmodel& SDM, const int& irun);
}

#endif
