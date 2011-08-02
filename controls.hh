// controls.hh

#ifndef SCALA_CONTROLS
#define SCALA_CONTROLS

#include "hkl_datatypes.hh"
#include "observationflags.hh"
#include "eprob.hh"
#include "range.hh"

namespace scala
{
  //  
  //------------------------------------------------------------
  class run_controls
  // Run definition 
  {
  public:
    run_controls();
    void SetStatus(const int& Status) {RunStatus = Status;}
    int Status() const {return RunStatus;}
    // true if set up (auto or specified)
    bool Set() const {return (RunStatus >= 0);}
    // return true if run specifications were given on input, 
    // irrespective of whether they have been imposed yet or not
    bool Explicit() const {return ((RunStatus == -2) || (RunStatus == +1));}

    void StoreRunBatchSelection(const BatchSelection& runselection)
    {batchrangesruns = runselection; RunStatus = -2;}

    //! return run number for batch, -1 if not in list, 0 if no selections
    int RunNumber(const int& batchnumber);

    //! return list of run numbers specified
    std::vector<int> RunNumberList() const;

    //! return list of batch ranges for specified run
    std::vector<IntRange> BatchRanges(const int& RunNumber)
    {return batchrangesruns.BatchRanges(RunNumber);}

    //! set list of RunNumber, resolution range from input
    void StoreResoByRun(const std::vector<std::pair<int,ResoRange> > Run_resolution_ranges);
    //! return list of RunNumber (NOT Run index!(), resolution range, if set 
    std::vector<std::pair<int,ResoRange> > GetResoByRun() const
    {return run_resolution_ranges;}
    //! clear all resolution by run selections
    void ClearResoByRun() {run_resolution_ranges.clear();} // don't reset resobyrun
    //! return true is any resolution-by-run ranges have ever been set
    bool IsResoByRun() const {return resobyrun;}

  private:
    // Run definition status
    //  -2  explicitly defined from input, but not yet done
    //  -1  initial value, nothing set, default to auto setting
    //   0  auto-setting done
    //  +1  explicitly defined from input and set
    int RunStatus;

    // Batch selections for explicitly defined runs
    BatchSelection batchrangesruns;

    std::vector<std::pair<int,ResoRange> > run_resolution_ranges; // run ID number and reso range
    bool resobyrun;  // true if there have ever been run_resolution_ranges set, even if cleared
  }; // run_controls
  //------------------------------------------------------------
  class partial_controls
  // Selection & treatment of partials 
  //  accept_fract_min, accept_fract_max
  //                  acceptable range of total fraction_calc to be
  //                  accepted as a complete partial
  //  correct_fract_min   minimum total fraction_calc for
  //                  observation to be scaled (if not complete)
  //                  = 0.0 don't 
  //  check           true if Mpart flags should be checked for
  //                  consistency
  //  maxgap          maximum accepted gap in batch number,
  //                  usually = 0 ie no gap
  {
  public:
    partial_controls();
    partial_controls(const double& FrMin, const double& FrMax,
		     const double& FrCorr, const bool& Check,
		     const int& MaxGap);

    std::string format() const;
    
    Rtype& accept_fract_min() {return accept_fract_min_;}  // Set
    Rtype accept_fract_min() const {return accept_fract_min_;} // Get
    Rtype& accept_fract_max() {return accept_fract_max_;}  // Set
    Rtype accept_fract_max() const {return accept_fract_max_;} // Get
    Rtype& correct_fract_min() {return correct_fract_min_;}  // Set
    Rtype correct_fract_min() const {return correct_fract_min_;} // Get
    bool& check() {return check_;}
    bool check() const {return check_;}
    int& maxgap() {return maxgap_;}
    int maxgap() const {return maxgap_;}

    void Clear(); // clear totals
    
    //  rejections because of gaps
    void IncrementNrejGap() {nrejgap++;}
    int NrejGap() const {return nrejgap;}
    //  rejections because total fraction too small
    void IncrementNrejFractionTooSmall() {nrejfractionsmall++;}
    int NrejFractionTooSmall() const {return nrejfractionsmall;}
    //  rejections because total fraction too large
    void IncrementNrejFractionTooLarge() {nrejfractionlarge++;}
    int NrejFractionTooLarge() const {return nrejfractionlarge;}

  private:
    Rtype accept_fract_min_, accept_fract_max_;
    Rtype correct_fract_min_;
    bool check_;
    int maxgap_;
    int nrejgap;               //  rejections because of gaps
    int nrejfractionsmall;     //  rejections because total fraction too small
    int nrejfractionlarge;     //  rejections because total fraction too large
  }; // partial_controls
  //------------------------------------------------------------
class col_controls
// Column selection options
//  (a) profile-fitted or integrated intensity column
{
public:
  col_controls() : prfpresent(false), SelectIcolFlag(-1), imid(-1.0), IpowerComb(3) {}

  // Set column selection flags
  // see class observation_part (hkl_unmerge.hh) for implementation
  //   INTEGRATED   integrated intensity I                  (SelectIcolFlag=0)
  //   PROFILE      profile-fitted intensity IPR            (SelectIcolFlag=-1)
  //   COMBINE      weighted mean after combining partials  (SelectIcolFlag=+1)
  void SetIcolFlag(const int& IcolFlag,  const double& Imid, const int Ipower=3);

  // Set up observation_part class (static method) with current flags
  // conditional on Ipr present or absent (col_Ipr == 0 absent)
  void SetupColSelection(const int col_Ipr=1);

  // Return selection flag
  int IcolFlag() const {return SelectIcolFlag;}

  //! return true if Imid was set from command input
  bool IsImidSet() const {return (imid > 0.1);}

  //! return true if file contains both profile-fitted and summation Is
  bool BothIpresent() const {return prfpresent;}

private:
  bool prfpresent;     // true if file contains both profile-fitted and summation Is
  int SelectIcolFlag;
  int IpowerComb;
  double imid;
}; // col_controls
 //------------------------------------------------------------
class AnalysisControls
// Things to control analysis of data
//  - resolution binning
//  - intensity binning
//  - anisotropy analysis
{
public:
  AnalysisControls() : nresobins(10), nibins(10) {}
  AnalysisControls(const int& Nresobins, const int& Nibins, const double& Coneangle,
		   const double& MinimumHalfdatasetCC,
		   const double& MinimumIoverSigma, const double& MinimumBatchIoverSigma,
		   const int& Nbatchsmooth=1)
    : nresobins(Nresobins), nibins(Nibins), coneangledegrees(Coneangle),
      minimumhalfdatasetcc(MinimumHalfdatasetCC), minimumioversigma(MinimumIoverSigma),
      minimumbatchioversigma(MinimumBatchIoverSigma), nbatchsmooth(Nbatchsmooth)
{}

  int NresoBins() const {return nresobins;}
  int NiBins() const {return nibins;}
  double ConeAngle() const {return coneangledegrees;}
  double MinimumHalfdatasetCC() const {return minimumhalfdatasetcc;}
  double MinimumIoverSigma() const {return minimumioversigma;}
  double MinimumBatchIoverSigma() const {return minimumbatchioversigma;}
  int NbatchSmooth() const {return nbatchsmooth;}
  void SetNbatchSmooth(const int& Nbatchsmooth) {nbatchsmooth = Nbatchsmooth;}

private:
  int nresobins;    // number of resolution bins
  int nibins;       // number of intensity bins
  double coneangledegrees;  // cone angle 
  double minimumhalfdatasetcc;
  double minimumioversigma;
  double minimumbatchioversigma;
  int nbatchsmooth; // number of batches over which to smooth statistics
};
//------------------------------------------------------------
class RejectFlags {
public:

  enum Reject2Policy {REJECT, KEEP, REJECTLARGER, REJECTSMALLER};

  RejectFlags(){}
  RejectFlags(const float& Sdrej, const float& Sdrej2,
	      const Reject2Policy& Rej2policy)
    : sdrej(Sdrej), sdrej2(Sdrej2), rej2policy(Rej2policy) {}

  std::string format() const;

  std::string formatReject2Policy() const;

  float sdrej;        // SD multiplier for outlier rejection
  float sdrej2;       //  special for two observations
  //  enum Reject2Policy {REJECT, KEEP, REJECTLARGER, REJECTSMALLER};
  Reject2Policy rej2policy;  // what to do with 2 deviant observations

};
//------------------------------------------------------------
class OutlierControl
{
public:
  OutlierControl(); // Set defaults: COMBINE
  OutlierControl(const int& Ndatasets); // Set defaults: COMBINE

  bool& Combine() {return combine;}       // Set
  bool Combine() const {return combine;}  // Get

  bool Anom() const; // return true if rejection is set between I+ & I- for all datasets

  void SetNdatasets(const int& Ndatasets);

  // rejection criteria, within I+, I-  or between I+ & I-
  // Set
  RejectFlags& Reject(const AnomalousClass& selclass, const int& dts_index=0);
  // Get
  RejectFlags Reject(const AnomalousClass& selclass, const int& dts_index=0) const;

  void SetEmax(const float& Emax); //!< set Emax (acentric)
  EProb EMaxTest() const {return emaxtest;}

private:
  bool combine;                // true for outlier checks between datasets
  RejectFlags reject;          // main rejection criteria, within I+, I- set
  // for each dataset
  std::vector<RejectFlags> rejectanom;      // between I+ & I-
  EProb emaxtest;  
}; // OutlierControl
//=================================================================
  //! Controls on scale refinement  
  class RefineControl
  {
  public:
    RefineControl(); // set defaults

    bool& BFGS() {return bfgs;} // set
    bool BFGS() const {return bfgs;} // get

    int& Ncyc1() {return ncyc1;}
    int Ncyc1() const {return ncyc1;}

    int& Ncycles() {return ncycles;}
    int Ncycles() const {return ncycles;}

    float& Converge() {return converge;}
    float Converge() const {return converge;}

    float& IovSDmin() {return iovsdmin;}
    float IovSDmin() const {return iovsdmin;}

    float& E2min() {return e2min;}
    float E2min() const {return e2min;}

    float& E2max() {return e2max;}
    float E2max() const {return e2max;}

  private:
    bool bfgs;           // or Fox-Holmes if false
    int  ncyc1;          // number of 1st stage cycles
    int  ncycles;        // number of main stage cycles
    float converge;      // convergence limit (multiplier of sd)
    float iovsdmin;      // <I>/sd'(<I>) limit for 1st pass scaling
    float e2min;         // |E^2| limit for 2nd pass scaling
    float e2max;         // |E^2| maximum limit for 2nd pass scaling
  };
//=================================================================
class DatasetControl
// selection of datasets etc
{
public:
  DatasetControl() :basedataset(-1){}

  int BaseDataset() const {return basedataset;}
private:
  int basedataset;
};
//=================================================================
class PolarisationControl {
  //! data for update of polarisation correction, for XDS/INTEGRATE files
public:
  PolarisationControl() : set(false), polarisationfactor(0.0) {}

  void SetFactor(const double& Polarisationfactor) {
    polarisationfactor = Polarisationfactor;
    set = true;

  }
  double Factor() const {return polarisationfactor;}
  bool IsSet() const {return set;}

  // Default value for synchrotron
  static double Default() {return +0.98;}

private:
  bool set;                  // true if the value has been explictly set from input
  double polarisationfactor; // = 0 for unpolarised, eg in-house source
};
//=================================================================
class all_controls
// all controls to store in hkl list
//  - run controls
//  - partial controls
//  - outlier controls
//  - anomalous on flag
//
// Just a public data structure
{
public:
  all_controls(){};

  run_controls runs;
  partial_controls partials;
  ObservationFlagControl observationflagcontrol; // acceptable flags
  AnalysisControls analysis;
  OutlierControl outlierScale;
  OutlierControl outlierMerge;
  RefineControl refinecontrol;	
  bool Anomalous;   // true if "anomalous on"
  bool AnomalousSDcorr;   // true to separate I+ & I- for SD correction
  DatasetControl datasetcontrol;
  PolarisationControl polarisationcontrol;
}; // all_controls
}  // namespace scala
#endif
