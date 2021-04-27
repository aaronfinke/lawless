// keyword_aimless.hh

#ifndef KEYWORDS_HEADER
#define KEYWORDS_HEADER

#include "CCP4base.hh"
#include "hkl_datatypes.hh"
#include "range.hh"
#include "controls.hh"
#include "sdctypes.hh"
#include "observationflags.hh"
#include "globalcontrols.hh"
#include "scaletypes.hh"
#include "weighttype.hh"

using scala::ResoRange;

namespace scala {
  class ScaleSpecification;
}

namespace phaser_io {
  //--------------------------------------------------------------
  class ANOMALOUS: public InputBase, virtual public CCP4base
    //  ANOMALOUS [ON] [OFF]
  {
  public:
    ANOMALOUS();
    virtual ~ANOMALOUS() {}
    Token_value parse(std::istringstream&);

    void setANOMALOUS(const bool& OnOff=true) {anomalous = OnOff;}
    bool getANOMALOUS() const {return anomalous;}
    // return true if the anomalous flag was explicitly given
    bool AnomalousFlagInput() const {return given;}

    void analyse(){}

  private:
    bool anomalous;
    bool given;      // true if ON or OFF were given
  };
  //--------------------------------------------------------------
  class SCALES: public InputBase, virtual public CCP4base
    // SCALES [RUN <irun>]
    // [BATCH || ROTATION [<nscales> || SPACING <spacing>]
    // BFACTOR ON || OFF BROTATION [<nscales> || SPACING <spacing>]
    // [SECONDARY  [<Lmax>]]
    // [ABSORPTION [<Lmax>] [POLE [h|k|l]]]
    // [CONSTANT]
    // [[NO]TILE [<Ntilex> [<Ntiley>]] [CCD[x] | FLAT | PIXEL]]
  {
  public:
    SCALES();
    virtual ~SCALES() {}
    Token_value parse(std::istringstream&);

    void analyse(){}

    std::vector<scala::ScaleSpecification> getScaleSpecifications() const
    {return specs;}

    bool noTile() const {return notile;}  // true if explicit NOTILE is given
    
  private:
    int nspecs;
    // 1st specification is always the default one (run -1)
    std::vector<scala::ScaleSpecification> specs;
    bool notile;  // true if explicit NOTILE is given
  };
  //--------------------------------------------------------------
  class RUNSET : public InputBase, virtual public CCP4base
  {
    // Syntax:
    //   RUN <runnumber> BATCH [FILE|SERIES  <Jfile>]  <b1> TO <b2>
    //   
    //   *** Not done *** RUN <runnumber> DATASET <datasetname> | <crystalname>/<datasetname>
    //   
  public:
    RUNSET();
    virtual ~RUNSET() {}
    Token_value parse(std::istringstream&);
    // Return list of batch ranges associated with a run and selection type
    scala::RunSelection Runsetselection() const {return runsetselection;}
    // Return list of batch ranges associated with a run
    scala::BatchSelection RunBatches() const
    {return runsetselection.batchranges;}

    void analyse(){}
  private:
    scala::RunSelection runsetselection;
  };
  //--------------------------------------------------------------
  class RESO : public InputBase, virtual public CCP4base
  {
    // Syntax:
    //  RESOlution [RUN <Irun>]  <high> [<low>]  either order
    //  RESOlution [RUN <Irun>] HIGH <high>
    //  RESOlution [RUN <Irun>] LOW  <low>
    //
    // Resolution by run may override earlier run commands 
  public:
    RESO();
    virtual ~RESO() {}
    Token_value parse(std::istringstream&);

    void setRESO(const ResoRange& reso_range);
    ResoRange getRESO() const;

    std::vector<std::pair<int,ResoRange> > GetResoByRun() const {return run_resolution_ranges;}

    void analyse(){}

    bool ValidRESO() const {return Valid;}

  private:
    ResoRange input_resolution_range;
    std::vector<std::pair<int,ResoRange> > run_resolution_ranges; // run ID number and reso range
    bool Valid;
  };
  //--------------------------------------------------------------
  class PARTIALS: public InputBase, virtual public CCP4base
  {
    // Syntax: PARTials [TEST <lower_limit> <upper_limit>] [CORRECT <minimum_fraction>] [NOCHECK] [[NO]GAP <gap limit>]
    //   
  public:
    PARTIALS();
    virtual ~PARTIALS() {}
    Token_value parse(std::istringstream&);


    void setFracLimMin(const double& MinLim) 
    {minfraclim = MinLim;}
    void setFracLimMax(const double& MaxLim) 
    {maxfraclim = MaxLim;}

    void setFracScaleLim(const double& MinLim)
    {minscalefrac = MinLim;}

    double getFracLimMin() const  {return minfraclim;}
    double getFracLimMax() const  {return maxfraclim;}

    double getSclMinLim() const {return minscalefrac;}
    bool getCheck() const {return check;}

    void setMaxGap(const int& MaxGap) {maxgap = MaxGap;}  // set maxgap
    int getMaxGap() const {return maxgap;} // return maxgap

    void analyse(){}

  private:
    double minfraclim, maxfraclim;
    double minscalefrac;
    bool check;
    int maxgap;
  };
  //--------------------------------------------------------------
  class EXCLUDE : public InputBase, virtual public CCP4base
  {
    // Syntax:
    //   EXCLUDE DATASET <datasetname> | <crystalname>/<datasetname>
    //     exclude dataset
    //   EXCLUDE BATCH [FILE|SERIES  <Jfile>]  <b1> <b2> <b3> ... | <b1> TO <b2>
    //     exclude batch list or range: if FILE key present, b1 etc
    //     refer to original (file) batch numbers (SERIES is the equivalent for a
    //     series of files defined with a wild-card)
    //   
  public:
    EXCLUDE();
    virtual ~EXCLUDE() {}
    Token_value parse(std::istringstream&);
    scala::BatchSelection BatchExclusions() const {return batchexclude;}

    void analyse(){}
  private:
    scala::BatchSelection batchexclude;
  };
  //--------------------------------------------------------------
  class BINS : public InputBase, virtual public CCP4base
    // BINS RESOLUTION <NresoBins> INTENSITY <NintensityBins>
  {
  public:
    BINS();
    virtual ~BINS(){}
    Token_value parse(std::istringstream&);

    int getResoBins() const {return nrbins;}
    int getIntBins() const {return nibins;}

    void analyse(){}

  private:
    int nrbins, nibins;
  };
  //--------------------------------------------------------------
  class REJECT : public InputBase, virtual public CCP4base
    // Set outlier flags
    //
    // Sub-keywords:
    //
    //  [SCALE|MERGE]  use these values for scaling|merging steps
    //            if not specified, use for both
    //  [COMBINE] compare observations across all datasets [default]
    //  [SEPARATE]  outlier checks only within datasets
    //   sdrej    sd multiplier for maximum deviation from scale-weighted mean I
    //  [sdrej2]  special value for reflections measured twice
    //  [REJECT|KEEP|SMALLER|LARGER]] these flags control what to do in the 
    //            merging stage if there are two observations which disagree
    //            In scaling, both are always rejected.
    //            REJECT      reject both
    //            KEEP        keep both [default]
    //            LARGER      reject larger
    //            SMALLER     reject smaller
    //  [ALL sdreja [sdrej2a]]
    //            if ALL is present and the ANOMALOUS option is on, then
    //            test all data (ie all I+ & I-) values against these limits
    //            after testing within the I+ & I- set
    //            NOALL (== ALL 0) switches off this test
    //            Only applies to merging step (scaling step checks all anyway)
    //  EMAX <Emax> maximum normalised F accepted
    //  BATCH  <batchrejectfactor>
    //          if > 0, reject batches with scales >
    //             batchrejectfactor * medianscale
    //           (and negatives)
    //           valid only for batch scaling eg for serial data
  {
  public:
    REJECT();
    virtual ~REJECT(){}
    Token_value parse(std::istringstream&);
    
    void analyse();

    scala::OutlierControl GetOutlierControlsScale() const
    {return outliercontrolsscale;}

    scala::OutlierControl GetOutlierControlsMerge() const
    {return outliercontrolsmerge;}

    // true if this command has been given
    bool isrejectset() const {return set;}

  private:
    bool set;   // true if set explicitly from coammnds
    scala::OutlierControl outliercontrolsscale;
    scala::OutlierControl outliercontrolsmerge;
  };
  //--------------------------------------------------------------
  class TIE : public InputBase, virtual public CCP4base
    // Define parameter restraints ("ties")
    //
    // TIE <parameter> <sd> [sd2>]
    //  <parameter> = SURFACE  for SECONDARY or ABSORPTION
    //              = ROTATION for primary scale parameters (eg BATCH)
    //              = BFACTOR  for B-factors
    //              = ZEROB    for B-factors tied to B = 0
    //              = TILE     for tile correction parameters (5 sds)
    //              = TARGETTILE targets for tile correction r & w
    //
    // SD parameters defaulted to -1 if no restraint
  {
  public:
    TIE();
    virtual ~TIE(){}
    Token_value parse(std::istringstream&);
    
    void analyse(){}
    
    // Return values
    float TIE_sd_surface()  const {return tiesd_surface;}
    float TIE_sd_rotation() const {return tiesd_rotation;}
    float TIE_sd_bfactor()  const {return tiesd_bfactor;}
    float TIE_sd_zerob()    const {return tiesd_zerob;}
    std::vector<double> TIE_tile  ()   const {return ties_tile;}

  private:
    float tiesd_surface;  // sd for SURFACE, < 0.0 for no restraint
    float tiesd_rotation; // sd for ROTATION, < 0.0 for no restraint
    float tiesd_bfactor;  // sd for BFACTOR, < 0.0 for no restraint
    float tiesd_zerob;    // sd for ZEROB, < 0.0 for no restraint
    // sds for TILE, < 0.0 for no restraint, 5 values for CCD tiles
    // + target values for r, w, ie 7 values
    std::vector<double> ties_tile;

  };
  //--------------------------------------------------------------
  class NAME : public InputBase, virtual public CCP4base
  {
    // Syntax: NAME PROJECT <Pname> CRYSTAL <Xname> DATASET <Dname>
  public:
    NAME();
    virtual ~NAME() {}
    Token_value parse(std::istringstream&);


    static void setNAME(const scala::PxdName& PXDname);
    static scala::PxdName getNAME() {return pxdname;}
    void analyse(){}

  private:
    static scala::PxdName pxdname;
  };
  //--------------------------------------------------------------
  class REFINE : public InputBase, virtual public CCP4base
  {
    // Syntax: REFINE BFGS|FH CYCLE [Ncyc1>] <Ncycles>
    //            SELECT <IovSDmin> <E2min> [<E2Max>]
    //            PARALLEL [AUTO] | <nproc> | <fproc>
    //    CYCLE
    //         Ncyc1   number of cycles in 1st stage [default 2]
    //         Ncycles number of cycles in main scaling [default 10]
    //    BFGS use BFGS minimiser
    //    FH   use Fox_Holmes least-squares
    // Selection criteria for scaling:
    //    IovSDmin   <I>/sd'(<I>) limit for 1st pass scaling
    //    E2min      |E^2| limit for 2nd pass scaling
    //    E2max      |E^2| limit for 2nd pass scaling
    // If OpenMP is enabled:
    //    PARALLEL  number of processors to use in scaling, or
    //              fraction of available processors to use, or
    //              AUTO determine a "best" number of processors to use
    //    If PARALLEL is specified without an argument, then AUTO is assumed
  public:
    REFINE();
    virtual ~REFINE() {}
    Token_value parse(std::istringstream&);

    scala::RefineControl RefineControl() const {return refinecontrol;}

    void analyse(){}

  private:
    scala::RefineControl refinecontrol;
  };
  //--------------------------------------------------------------
  class TITLE : public InputBase, virtual public CCP4base
  {
    // Syntax: TITLE <title>
    // Title for run
  public:
    TITLE();
    virtual ~TITLE() {}
    Token_value parse(std::istringstream&);

    std::string Title() const {return title;}
    void analyse(){}

  private:
    std::string title;
  };
  //--------------------------------------------------------------
  class ONLYMERGE : public InputBase, virtual public CCP4base
  {
    // Syntax: ONLYMERGE
    //   No scaling
  public:
    ONLYMERGE();
    virtual ~ONLYMERGE() {}
    Token_value parse(std::istringstream&);

    bool Onlymerge() const {return onlymerge;}
    void analyse(){}

  private:
    bool onlymerge;
  };
  //--------------------------------------------------------------
  class BLANK : public InputBase, virtual public CCP4base
    //
    // Read parameters for detecting blank images
    // Syntax:
    //  BLANK [RESO[FRACTION] <fraction>] [NEGATIVE <negativefraction>]
    //    <fraction>   blank image test will be done out to this fraction
    //                 of the maximum resolution, default 0.5
    //    <negativefraction> fraction of reflections with negative intensity
    //                 above which the image is considered to be blank
    //                 Default 0.3
    //   
  {
  public:
    BLANK();
    virtual ~BLANK() {}
    Token_value parse(std::istringstream&);
    float NullResolutionfraction() const {return nullResolutionfraction;}
    float NullNegativeReject() const {return nullNegativeReject;}
    
    void analyse(){}
    
  private:
    float nullResolutionfraction;
    float nullNegativeReject;
    bool isset;
  };
  //--------------------------------------------------------------
  class SDCORRECTION : public InputBase, virtual public CCP4base
    //
    // Read parameters for SD correction
    // Syntax:
    //   SDCORRECTION [[NO]REFINE] [INDIVIDUAL|SAME] [FIXSDB]
    //     [RUN <RunNumber>] [FULL | PARTIAL] <SdFac> [<SdB>] <SdAdd>
    //     DAMP <dampfactor>
    //     TIE [<parameter> <value> <sd>] | NOTIE
    //  <parameter> is "SdFac" "SdB" or "SdAdd" (case insensitive)
    //     SIMILAR <sd1> [<sd2>] <sd3>   for SDfac, [SDb,] SDadd 
    //     WEIGHT VARIANCE | UNIT | SQRTSCALE  set weighting scheme
    //        for averaging Ih in calculating deviations
    //        VARIANCE  w = 1/var(I)  [default]
    //        UNIT      w = 1
    //        SQRTSCALE w = 1/sqrt(g) = sqrt(scale)
    //        SCALE     w = g = 1/scale
    //        SAMPLESD use sample SD in final averaging
    //        
    // NB LINEAR and GAUSSIAN options don't work - do not use
  {
  public:
    SDCORRECTION();
    virtual ~SDCORRECTION() {}
    Token_value parse(std::istringstream&);
    void analyse();

    //! return != 0 if parameters should be refined
    // refine = +1 refine (non-linear), 0 no refine, -1 linear fit
    int SDC_Refine() const {return refine;}
    //! return true if SDC REFINE is explicitly set
    bool SDC_RefineSet() const {return refine_set;}
    //! return true if the same values should be used for all runs
    bool SDC_AllRunsSame() const {return allsame;}
    //! return true if fixed SdB = 0.0
    bool SDC_NoSDb() const {return fixsdb;}

    //! return number of input corrections, = 0 none, = -1 overall
    int SDC_NumberInput() const;
    //! return list of given SD corrections, full, partial pairs
    std::vector<std::pair<scala::SDcorrection,scala::SDcorrection> >
       SDC_SDcorrections() const {return sdinput;}
    //! return corresponding run numbers, = -1 for all runs
    std::vector<int> SDC_RunNumbers() const {return runnumbers;}

    //! Return damp factor
    double SDCdamp() const {return damp;}

    //! return tietype, target & SDs, = 0 no tie, = -1 defaults, = +1 set here,
    //   = +2 similarity
    int SDCties(std::vector<double>& Targets,
		std::vector<double>& SDtarget) const;

    //! return weight type
    scala::WeightType::AverageWeightType SDCweightType() const
      {return weighttype;}

    //! return sampleSD flag
    bool SampleSD() const {return sampleSD;}

    //! return set flag, true if this command has been given
    bool isSDcorrectionSet() const {return set;}

  private:
    // refine = +1 refine (non-linear), 0 no refine, -1 linear fit
    int refine;
    bool refine_set; // true if explicit refine flag set
    bool allsame; // true for all runs the same
    bool fixsdb;  // true for fixed SdB = 0.0
    // list of given SD corrections (full, partial)
    std::vector<std::pair<scala::SDcorrection,scala::SDcorrection> > sdinput;
    std::vector<int> runnumbers;       // corresponding run numbers, = -1 all  
    double damp;
    int tietype;     // = 0 no tie, = -1 defaults, = +1 set here, = +2 similar
    std::vector<double> targets;   // 3 targets
    std::vector<double> sdtargets;  // ... and their SDs (= 0 no target)
    scala::WeightType::AverageWeightType weighttype;
    bool sampleSD;
    bool set;  // true if this command has been given
  };
  //--------------------------------------------------------------
  class INTENSITIES : public InputBase, virtual public CCP4base
    //
    // Read parameters for selecting intensities
    // Syntax:
    // INTENSITIES [SUMMATION | PROFILE | COMBINE [<Imid>] [POWER <Ipower>] ]

  {
  public:
    INTENSITIES();
    virtual ~INTENSITIES() {}
    Token_value parse(std::istringstream&);

    int GetIcolFlag() const {return selecticolflag;}
    double GetCombineImid() const {return imid;}
    int GetCombinePower() const {return ipowercomb;}

    void analyse(){}

  private:
    int selecticolflag;
    int ipowercomb;
    double imid;
  };
  //--------------------------------------------------------------
  class KEEP : public InputBase, virtual public CCP4base
    //
    // Read options to keep flagged observations, for merging step only
    // Syntax:
    // KEEP [OVERLOADS|BGRATIO <bgratio_max>|PKRATIO <pkratio_max>|
    //       GRADIENT <bg_gradient_max>|EDGE]
  {
  public:
    KEEP();
    virtual ~KEEP() {}
    Token_value parse(std::istringstream&);

    //! return observation flag controls
    scala::ObservationFlagControl Observationflagcontrol() const
    {return observationflagcontrol;}

    void analyse(){}

  private:
    scala::ObservationFlagControl observationflagcontrol;
  };
  //--------------------------------------------------------------
  class OUTPUT : public InputBase, virtual public CCP4base
    //
    // Read options to define output reflection files
    // Syntax:
    //  OUTPUT [MTZ] [NO]MERGED | UNMERGED [SPLIT | TOGETHER]
    //        [POLISH MERGED | UNMERGED]
    //        [ORIGINAL | REDUCED]
  {
  public:
    OUTPUT();
    virtual ~OUTPUT() {}
    Token_value parse(std::istringstream&);

    scala::OutputControls Outputcontrols() const {return outputcontrols;}

    void analyse(){}

  private:
    scala::OutputControls outputcontrols;
  };
  //--------------------------------------------------------------
  class DUMP : public InputBase, virtual public CCP4base
    //
    //  Dump scale model from file [default SCALES]
    // Syntax:
    //  DUMP <filename>
  {
  public:
    DUMP();
    virtual ~DUMP() {}
    Token_value parse(std::istringstream&);

    std::string  DumpFileName() const {return dumpfilename;}

    void analyse(){}

  private:
    std::string dumpfilename;
  };
  //--------------------------------------------------------------
  class RESTORE : public InputBase, virtual public CCP4base
    //
    //  Restore scale model from file [default SCALES]
    // Syntax:
    //  RESTORE <filename>
  {
  public:
    RESTORE();
    virtual ~RESTORE() {}
    Token_value parse(std::istringstream&);
    //! true if restore set
    bool Restore() const {return restore;}
    //! Filename, default "SCALES"
    std::string  RestoreFileName() const {return restorefilename;}

    void analyse(){}

  private:
    bool restore;
    std::string restorefilename;
  };
  //--------------------------------------------------------------
  class ANALYSIS : public InputBase, virtual public CCP4base
    //
    //  Controls for ANALYSIS analysis
    // Syntax:
    //  ANALYSIS  CONE <angle> 
    //    CCMINIMUM <MinimumHalfdatasetCC>
    //    CCANOMMINIMUM <MinimumHalfdatasetAnomCC>
    //    ISIGMINIMUM <MinimumIoverSigma>
    //    BATCHISIGMINIMUM <MinimumBatchIoverSigma>
    //    SMOOTHSTATISTICS <SmoothStatisticsRange>
    //    GROUPBATCH  <BatchGroupRange>
    //
    // Cone angle is the half-angle (degrees) for cones around each reciprocal axis
    // MinimumHalfdatasetCC  minimum CC for resolution warning
    // MinimumHalfdatasetAnomCC  minimum CCanom for resolution warning
    // MinimumIoverSigma          minimum <<I>/sd(<I>)> for resolution warning
    // MinimumBatchIoverSigma     minimum <I/sd(I)> for resolution warning by batch, from unmerged I
    // SmoothStatisticsRange angle in degrees over which (roughly) to smooth
    //            batch statistics, <0 to default to automatic setting,
    // BatchGroupRange        phi range for grouping batches in analysis
    //
    //   ANALYSIS [NO]DETECTOR
    //     Analyse (or not) scales on the detector, images written to DETECTORIMAGE
    //
  {
  public:
    ANALYSIS();
    virtual ~ANALYSIS() {}
    Token_value parse(std::istringstream&);

    double ConeAngle() const {return coneangledegrees;}
    double MinimumHalfdatasetCC() const {return minimumhalfdatasetcc;}
    double MinimumHalfdatasetAnomCC() const {return minimumhalfdatasetanomcc;}
    double MinimumIoverSigma() const {return minimumioversigma;}
    double MinimumBatchIoverSigma() const {return minimumbatchioversigma;}
    double SmoothStatisticsRange() const {return smoothstatisticsrange;}
    double BatchGroupRange() const {return batchgrouprange;}
    bool DetectorAnalysis() const {return detector;}  // detector analysis flag

    void analyse(){}

  private:
    double coneangledegrees;
    double minimumhalfdatasetcc;
    double minimumhalfdatasetanomcc;
    double minimumioversigma;
    double minimumbatchioversigma;
    double smoothstatisticsrange;
    double batchgrouprange;
    // Detector analysis things
    bool detector;
  };
  //--------------------------------------------------------------
  class INITIAL : public InputBase, virtual public CCP4base
    //
    //  Controls for INITIAL scale
    // Syntax:
    //  INITIAL MEAN   set all initial scales from mean intensities [default]
    //  INITIAL UNITY  set all initial scales to unity
    //  INITIAL MINIMUM_OVERLAP  read minimum overlap fraction
    //  INITIAL MAXIMUM_GAP  read maximum contiguous "gaps" in rotation ranges
  {
  public:
    INITIAL();
    virtual ~INITIAL() {}
    Token_value parse(std::istringstream&);

    bool InitialUnity() const {return unity;}

    double Minimum_overlap() const {return minimum_overlap;}

    int Maximum_gap() const {return maximum_gap;}

    void analyse(){}

  private:
    bool unity;
    double minimum_overlap;
    int maximum_gap;
  };
  //--------------------------------------------------------------
  class XMLOUT : public InputBase, virtual public CCP4base
  {
    // Syntax: XMLOUT <filename>
  public:
    XMLOUT();
    virtual ~XMLOUT() {}
    Token_value parse(std::istringstream&);


    void setXMLOUT(const std::string& Name) {name = Name;}
    std::string getXMLOUT() const {return name;}
    void analyse(){}

  private:
    std::string name;
  };
  //--------------------------------------------------------------
  class ROGUES : public InputBase, virtual public CCP4base
  {
    // Syntax: ROGUES <filename>
  public:
    ROGUES();
    virtual ~ROGUES() {}
    Token_value parse(std::istringstream&);


    void setROGUES(const std::string& Name) {name = Name;}
    std::string getROGUES() const {return name;}
    void analyse(){}

  private:
    std::string name;
  };
  //--------------------------------------------------------------
  class HKLIN : public InputBase, virtual public CCP4base
  {
    // Syntax: HKLIN <filename>
  public:
    HKLIN();
    virtual ~HKLIN() {}
    Token_value parse(std::istringstream&);


    void setHKLIN(const std::string& Name) {name = Name;}
    std::string getHKLIN() const {return name;}
    void analyse(){}

  private:
    std::string name;
  };
  //--------------------------------------------------------------
  class HKLOUT : public InputBase, virtual public CCP4base
  {
    // Syntax: HKLOUT <filename>
  public:
    HKLOUT();
    virtual ~HKLOUT() {}
    Token_value parse(std::istringstream&);


    void setHKLOUT(const std::string& Name) {name = Name;}
    std::string getHKLOUT() const {return name;}
    void analyse(){}

  private:
    std::string name;
  };
  //--------------------------------------------------------------
  class UNMERGEDOUT : public InputBase, virtual public CCP4base
  {
    // Syntax: UNMERGEDOUT <filename>
  public:
    UNMERGEDOUT();
    virtual ~UNMERGEDOUT() {}
    Token_value parse(std::istringstream&);


    void setUNMERGEDOUT(const std::string& Name) {name = Name;}
    std::string getUNMERGEDOUT() const {return name;}
    void analyse(){}

  private:
    std::string name;
  };
  //--------------------------------------------------------------
  class HKLREF : public InputBase, virtual public CCP4base
  {
    // Syntax: HKLREF <filename>
  public:
    HKLREF();
    virtual ~HKLREF() {}
    Token_value parse(std::istringstream&);


    void setHKLREF(const std::string& Name) {name = Name;}
    std::string getHKLREF() const {return name;}
    void analyse(){}

  private:
    std::string name;
  };
  //--------------------------------------------------------------
  class LABREF : public InputBase, virtual public CCP4base
  {
    // Syntax: LABREF [F|I = ] <F|Ilabel> [[SIGF|I = ] <sigF|Ilabel>]
  public:
    LABREF();
    virtual ~LABREF() {}
    Token_value parse(std::istringstream&);

    void setLABREF_I(const std::string& lab_I) {FIlabel = lab_I;}
    void setLABREF_sigI(const std::string& lab_sigI) {sigFIlabel = lab_sigI;}
    std::string getLABREF_I() const {return FIlabel;}
    std::string getLABREF_sigI() const {return sigFIlabel;}
    void analyse(){}

  private:
    std::string FIlabel;
    std::string sigFIlabel;

  };
  //--------------------------------------------------------------
  class XYZIN : public InputBase, virtual public CCP4base
  {
    // Syntax: XYZIN <filename>
  public:
    XYZIN();
    virtual ~XYZIN() {}
    Token_value parse(std::istringstream&);


    void setXYZIN(const std::string& Name) {name = Name;}
    std::string getXYZIN() const {return name;}
    void analyse(){}

  private:
    std::string name;
  };
  //--------------------------------------------------------------
  class USESDPARAMETER : public InputBase, virtual public CCP4base
  {
    // Syntax: USESDPARAMETER [NO | DIAGONAL | COVARIANCE]
    //    (default DIAGONAL if not explicit)
  public:
    USESDPARAMETER();
    virtual ~USESDPARAMETER() {}
    Token_value parse(std::istringstream&);

    void setUSESDPARAMETER
    (const scala::ScaleSpecification::ParameterSDusage& Use_sd_parameter)
    {parametersdusage = Use_sd_parameter;}
    scala::ScaleSpecification::ParameterSDusage getUSESDPARAMETER() const
    {return parametersdusage;}
    void analyse(){}

  private:
    //  {NONE, DIAGONAL, COVARIANCE};
    scala::ScaleSpecification::ParameterSDusage parametersdusage; 

  };
  //--------------------------------------------------------------
  class LINK : public InputBase, virtual public CCP4base
  {
    // Syntax: [UN]LINK [SURFACE] [ALL] | <run2> TO <run1>
  public:
    LINK();
    virtual ~LINK() {}
    Token_value parse(std::istringstream&);

    scala::LinkSpecs getLINKs() const
    {return linkspecs;}
    void analyse(){}

  private:
    scala::LinkSpecs linkspecs; 

};
  //--------------------------------------------------------------
  class PLOT : public InputBase, virtual public CCP4base
  {
    // Option to suppress XMgrace file output 
    //    NORMPLOT, ANOMPLOT, ROGUEPLOT, CORRELPLOT
    //    
    // Syntax: PLOT NOXMGR | XMGR
  public:
    PLOT();
    virtual ~PLOT() {}
    Token_value parse(std::istringstream&);

    bool XMGRoutput() const {return xmgroutput;}
    void setXMGRoutput(const bool& xmgrout) {xmgroutput = xmgrout;}
    void analyse(){}

  private:
    bool xmgroutput;

};
} // phaser_io

#endif
