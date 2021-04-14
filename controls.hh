// controls.hh

#ifndef SCALA_CONTROLS
#define SCALA_CONTROLS

#include "hkl_datatypes.hh"
#include "observationflags.hh"
#include "eprob.hh"
#include "range.hh"
#include "runthings.hh"
#include "weighttype.hh"

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

    //! store explicit run selection options
    void StoreRunBatchSelection(const RunSelection& runselection);

    //! return run number for batch, -1 if not in list, 0 if no selections
    int RunNumber(const int& batchnumber);

    //! return list of run numbers specified
    std::vector<int> RunNumberList() const;

    //! return list of batch ranges for specified run, or empty range for BYFILE or AUTO
    std::vector<IntRange> BatchRanges(const int& RunNumber);

    //! Return true if in selection, for given file filenum
    bool InSelection(const int& runnum,
		     const int& batchnum,
		     const int& originalbatchnum,
		     const int& filenum) const;

    // return true if run specifications were given on input, 
    // irrespective of whether they have been imposed yet or not
    bool Explicit() const;

    //! Return true if BYFILE
    bool Byfile() const
    {return (runsettype == RunSelection::BYFILE);}

    //! Return true if AUTO
    bool Auto() const
    {return (runsettype == RunSelection::AUTO);}

    //! Store use run flags
    void StoreUseRun(const std::vector<bool>& userun) {
      userun_ = userun;}
    //! Return use run flags
    std::vector<bool> UseRun() const {return userun_;}

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
    RunSelection::RunSetType runsettype;

    std::vector<bool> userun_;  // runs to use, see ChooseRuns & SetRunsToUse

    std::vector<std::pair<int,ResoRange> > run_resolution_ranges; // run ID number and reso range
    bool resobyrun;  // true if there have ever been run_resolution_ranges set, even if cleared
  }; // run_controls
  //-----------------------------------------------------------
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

    // nopartials true if we know there are no partials
    void setNoPartials(const bool& noPartials);
    bool noPartials() const {return nopartials;}

  private:
    Rtype accept_fract_min_, accept_fract_max_;
    Rtype correct_fract_min_;
    bool check_;
    int maxgap_;
    int nrejgap;               //  rejections because of gaps
    int nrejfractionsmall;     //  rejections because total fraction too small
    int nrejfractionlarge;     //  rejections because total fraction too large
    bool nopartials;           // true if we know there are no partials
  }; // partial_controls
  //------------------------------------------------------------
class col_controls
// Column selection options
//  (a) profile-fitted or integrated intensity column
{
public:
  col_controls() : prfpresent(false), SelectIcolFlag(0), IpowerComb(3), imid(-1.0) {}

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
//  - batch binning
{
public:
  AnalysisControls() : nresobins(10), nibins(10), batchgroupwidth(-1.0) {}
  AnalysisControls(const int& Nresobins, const int& Nibins, const double& Coneangle,
		   const double& MinimumHalfdatasetCC,
		   const double& MinimumHalfdatasetAnomCC,
		   const double& MinimumIoverSigma, const double& MinimumBatchIoverSigma,
		   const double& Smoothstatisticsrange, const double& Batchgroupwidth,
		   const bool& Detectoranalysis)
    : nresobins(Nresobins), nibins(Nibins), coneangledegrees(Coneangle),
      minimumhalfdatasetcc(MinimumHalfdatasetCC),
      minimumhalfdatasetanomcc(MinimumHalfdatasetAnomCC),
      minimumioversigma(MinimumIoverSigma),
      minimumbatchioversigma(MinimumBatchIoverSigma), 
      smoothstatisticsrange(Smoothstatisticsrange),
      batchgroupwidth(Batchgroupwidth),
      detectoranalysis(Detectoranalysis)
{}

  int NresoBins() const {return nresobins;}
  int NiBins() const {return nibins;}
  double ConeAngle() const {return coneangledegrees;}
  double MinimumHalfdatasetCC() const {return minimumhalfdatasetcc;}
  double MinimumHalfdatasetAnomCC() const {return minimumhalfdatasetanomcc;}
  double MinimumIoverSigma() const {return minimumioversigma;}
  double MinimumBatchIoverSigma() const {return minimumbatchioversigma;}
  double SmoothStatisticsRange() const {return smoothstatisticsrange;}
  double Batchgroupwidth() const {return batchgroupwidth;}
  void  SetBatchgroupwidth(const double& Batchgroupwidth) {batchgroupwidth = Batchgroupwidth;}
  bool DetectorAnalysis() const {return detectoranalysis;}

private:
  int nresobins;    // number of resolution bins
  int nibins;       // number of intensity bins
  double coneangledegrees;  // cone angle 
  double minimumhalfdatasetcc;
  double minimumhalfdatasetanomcc;
  double minimumioversigma;
  double minimumbatchioversigma;

  // phi range over which to smooth statistics
  double smoothstatisticsrange;

  // for analysis, batches may be grouped by batchgroupwidth degrees,
  // <= 0.0 to use actual batch number
  double batchgroupwidth;
  bool detectoranalysis;
};
//------------------------------------------------------------
class RejectFlags {
public:

  enum Reject2Policy {REJECT, KEEP, REJECTLARGER, REJECTSMALLER};

  RejectFlags(){}
  RejectFlags(const float& Sdrej, const float& Sdrej2,
	      const Reject2Policy& Rej2policy,
	      const float& batchRejectFactor=-1.0)
    : sdrej(Sdrej), sdrej2(Sdrej2), rej2policy(Rej2policy),
      batchrejectfactor(batchRejectFactor) {}

  std::string format() const;

  std::string formatReject2Policy() const;

  std::string formatXML(const std::string& tag,
			const std::string& attribute="") const;

  float sdrej;        // SD multiplier for outlier rejection
  float sdrej2;       //  special for two observations
  //  enum Reject2Policy {REJECT, KEEP, REJECTLARGER, REJECTSMALLER};
  Reject2Policy rej2policy;  // what to do with 2 deviant observations

  // For Batch mode rejection, valid only for batch scaling eg for serial data:
  //  if > 0, reject batches with scales > batchrejectfactor * medianscale
  //  ( and negatives)
  float batchrejectfactor;

};
//------------------------------------------------------------
class OutlierControl
{
public:

  // Outlier rejection: 
  //  No outlier check, [just on discrepancies, just on Emax,] Both [default]
  enum OutlierPolicy {NOREJECT, REJECTDIFFERENCES, REJECTEMAX, REJECTBOTH};

  OutlierControl(); // Set defaults: COMBINE
  OutlierControl(const int& Ndatasets); // Set defaults: COMBINE

  void init(const int& Ndatasets);

  bool& Combine() {return combine;}       // Set
  bool Combine() const {return combine;}  // Get

 // return true if rejection is set between I+ & I- for all datasets
  bool Anom() const {return anomreject;}

  // Weighting scheme for averaging observations in outlier testing
  WeightType::AverageWeightType weightType() const {return weighttype;}

  void SetNdatasets(const int& Ndatasets); // copy outliercontrols from 1st dataset

  // rejection criteria, within I+, I-  or between I+ & I-
  // Set rejection flags:
  // if selclass == BOTH, then either set a specific dataset,
  //   or all if dts_index<0
  // AnomalousClass = ALL   ignore differences between I+ & I-
  //                = BOTH  for test between I+ & I-
  //                = IPLUS, IMINUS  I+ or I-
  void SetReject(const RejectFlags& flags,
		 const AnomalousClass& selclass, const int& dts_index=-1);
  // Get
  RejectFlags Reject(const AnomalousClass& selclass, const int& dts_index=0) const;

  void SetEmax(const float& Emax); //!< set Emax (acentric)
  EProb EMaxTest() const {return emaxtest;}
  // true if there is an Emax test
  bool isEmaxTest() const {return !emaxtest.Null();}

  void SetOutlierPolicy(const OutlierPolicy& Outlierpolicy) {outlierpolicy = Outlierpolicy;}
  OutlierPolicy GetOutlierPolicy() const {return outlierpolicy;}

  std::string formatXML() const;  // format for XML

private:
  bool combine;                // true for outlier checks between datasets
  RejectFlags reject;          // main rejection criteria, within I+, I- set
  // for each dataset
  std::vector<RejectFlags> rejectanom;      // between I+ & I-
  bool anomreject;             // true if any rejection between I+ & I-
  EProb emaxtest;  
  // Weighting scheme for averaging observations in outlier testing
  WeightType::AverageWeightType weighttype;   // type of weighting for average
  OutlierPolicy outlierpolicy;

}; // OutlierControl
//=================================================================
  //! Controls on scale refinement  
  class RefineControl
  {
  public:
    RefineControl(); // set defaults

    // method = +1 BFGS, = 0 scale to reference, = -1 Fox-Holmes
    void setMethodBFGS() {method = +1;}
    void setMethodReference() {method = 0;}
    void setMethodFH() {method = -1;}
    bool BFGS() const {return (method>0);} // get
    bool Reference() const {return (method==0);} // get
    bool FH() const {return (method<0);} // get

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

    int& Nprocs() {return nprocs;}
    int Nprocs() const {return nprocs;}

    // Set number of processors from number or fraction (<1)
    void SetNprocs(const float& fproc);

    std::string format() const;

  private:
    // method = +1 BFGS, = 0 scale to reference, = -1 Fox-Holmes
    int method;
    int  ncyc1;          // number of 1st stage cycles
    int  ncycles;        // number of main stage cycles
    float converge;      // convergence limit (multiplier of sd)
    float iovsdmin;      // <I>/sd'(<I>) limit for 1st pass scaling
    float e2min;         // |E^2| limit for 2nd pass scaling
    float e2max;         // |E^2| maximum limit for 2nd pass scaling
    int nprocs;          // number of processors to use
    int maxprocs;        // maximum number
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
class PolarizationControl {
  //! data for update of polarization correction, for XDS/INTEGRATE files
  //
  // see Kahn et al (1982) J Appl Cryst 15, 330-337
  // Since this is used for XDS data, the internal conventions here are
  // as in XDS, except for E'pi (eprimepi)

public:
  PolarizationControl() : setfraction(false), setnormal(false),
			  polarizationfraction(-1.0) {}
  void setFraction(const double& Polarizationfraction)
  // Polarizationfactor = 0.5 for unpolarised, eg in-house source,
  //    ~+0.99 for a synchrotron
  //  Note that this is the definition as in XDS
  {
    polarizationfraction = Polarizationfraction;
    setfraction = true;
  }

  void setDirection(const clipper::Vec3<double>& PN)
  // PN = polarization normal ie perpendicular to synchrotron plane (as in XDS)
  {
    Pn = PN;
    setnormal = true;
  }

  double Factor() const {
    if (setfraction) {
      if (polarizationfraction == 0.0) { // explicitly turned off
	return polarizationfraction;
      }
      return 2.0*(polarizationfraction - 0.5);
    } else {
      // not set, return default
      return Default();
    }
  }

  double Fraction() const {
    return polarizationfraction;
  }

  clipper::Vec3<double> Normal() const {return Pn;}

  // component of electric vector in synchrotron plane, = Pn x s0, unit vector
  //   Cambridge frame
  void setEprimePi(const clipper::Vec3<double>& Eprimepi) {
    eprimepi = Eprimepi;}

  //   Cambridge frame
  clipper::Vec3<double> EprimePi() const {
    return eprimepi;
  }

  bool IsFractionSet() const {return setfraction;}
  bool IsNormalSet() const {return setnormal;}

  // Default value for synchrotron
  static double Default() {return +0.99;}

private:
  bool setfraction;  // true if fraction set
  bool setnormal;    // true if normal set
  // = 0.5 for unpolarised, eg in-house source, ~+0.99 for a synchrotron
  //  Note that this is the definition as in XDS
  // The electrical field vector of the incident beam is found in
  // the x,z-plane of the laboratory coordinate system with a
  // probability of polarizationfraction
  double polarizationfraction ;
  clipper::Vec3<double> Pn;         // polarization normal, XDS frame
  clipper::Vec3<double> eprimepi;
};
//=================================================================
// Controls for analysis and selection of anomalous
class AnomalousControl {
public:
  AnomalousControl();

  // All public
  bool FlagInput;   // true if explicit "anomalous" command given
  bool Anomalous;   // true if "anomalous on"
  bool AnomalousSDcorr;   // true to separate I+ & I- for SD correction
  
  // At present, anomalous scattering is considered to be present if any one of
  // the following is true (defaults in brackets):
  //  1) Anomplot slope > anomslopethreshold (1.3)
  //  2) CCanom > anomCCthreshold (0.3) in more than anomNbinthreshold bins (2)
  //  3) RCRanom > anomRCRthreshold (1.3) in more than anomNbinthreshold bins (2)
  double anomslopethreshold;
  double anomCCthreshold;
  double anomRCRthreshold;
  int anomNbinthreshold;
};
//=================================================================
class PlotControl
// Controls correlation plots and scatterplots 
{
 public:
  PlotControl() : xmgraceoutput(true) {}

  bool xmgraceoutput;  // true to output .xmgr files

 };
//=================================================================
class all_controls
// all controls to store in hkl list
//  - run controls
//  - partial controls
//  - outlier controls
//  - anomalous controls (on flag)
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
  AnomalousControl anomalouscontrol;
  DatasetControl datasetcontrol;
  PolarizationControl polarizationcontrol;
  PlotControl plotcontrol;
}; // all_controls
}  // namespace scala
#endif
