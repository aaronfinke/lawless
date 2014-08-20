// controls.cpp

#include <sstream>
#if _OPENMP
#include <omp.h>
#endif

#include "controls.hh"
#include "hkl_unmerge.hh"
#include "string_util.hh"

namespace scala
{
  //------------------------------------------------------------
  run_controls::run_controls() : RunStatus(-1), resobyrun(false) {}
 //------------------------------------------------------------
  //! return run number for batch, -1 if not in list, 0 if no selections
  int run_controls::RunNumber(const int& batchnumber)
  {
    if (batchrangesruns.Null()) return 0; // no selections
    return batchrangesruns.FlagNumber(batchnumber);
  }
 //------------------------------------------------------------
  std::vector<int> run_controls::RunNumberList() const
  // List of run numbers specified
  {
    return  batchrangesruns.UniqueFlags();
  }
 //------------------------------------------------------------
  //! set list of RunNumber, resolution range from input
  void  run_controls::StoreResoByRun
  (const std::vector<std::pair<int,ResoRange> > Run_resolution_ranges)
  {
    run_resolution_ranges = Run_resolution_ranges;
    resobyrun = (run_resolution_ranges.size() > 0);
  }
 //------------------------------------------------------------
  partial_controls::partial_controls()
    //    : accept_fract_min(0.95), accept_fract_max(1.05),
    //      correct_fract_min(0.0), check(true)
  {
    accept_fract_min_ = 0.95;
    accept_fract_max_ = 1.05;
    correct_fract_min_ = 0.0;
    check_ = true;
    maxgap_ = 0;
    nopartials = false;  // assume partials
  }
  //------------------------------------------------------------
  partial_controls::partial_controls(const double& FrMin, const double& FrMax,
                                     const double& FrCorr, const bool& Check,
                                     const int& MaxGap)
  {
    accept_fract_min_ = FrMin;
    accept_fract_max_ = FrMax;
    correct_fract_min_ = FrCorr;
    check_ = Check;
    maxgap_ = MaxGap;
    nopartials = false;  // assume partials
  }
  //------------------------------------------------------------
  void partial_controls::Clear()
  {
    nrejgap = 0;
    nrejfractionsmall = 0;
    nrejfractionlarge = 0;
  }
  //------------------------------------------------------------
  // nopartials true if we know there are no partials
  void partial_controls::setNoPartials(const bool& noPartials)
  {
    nopartials = noPartials;
    if (nopartials) {Clear();}
  }
  //------------------------------------------------------------
  std::string partial_controls::format() const
  {
    std::string s;
    if (nopartials) {
      {s += FormatOutput::logTab
          (0,"\nNo partially recorded reflections\n");}
    } else {
      if (check_)
        {s += FormatOutput::logTab
            (0,"\nHandling of partials:\n  MPART flags are checked\n");}
      else
        {s += FormatOutput::logTab
            (0,"\nHandling of partials:\n  MPART flags are not checked\n");}
      s += FormatOutput::logTabPrintf
        (0,"  Summed partials accepted if total fraction is between %5.2f & %5.2f\n",
         accept_fract_min_, accept_fract_max_);
      if (correct_fract_min_ > 0.001)
        {s += FormatOutput::logTabPrintf(0,
                                         "  Incomplete partials scaled by 1/fraction_calc if fraction is > %5.2f\n",
                                         correct_fract_min_);
        }
      if (maxgap_ > 0)
        {s += FormatOutput::logTabPrintf(0,"  Partials with up to %2d missing parts in the middle will be accepted\n",
                                         maxgap_);
        }
    }
    return s;
  }
  //------------------------------------------------------------
  void col_controls::SetIcolFlag(const int& IcolFlag, const double& Imid,
                                 const int Ipower)
  // Set column selection flags
  // see class observation_part (hkl_unmerge.hh) for implementation
  //   INTEGRATED   integrated intensity I                  (SelectIcolFlag=0)
  //   PROFILE      profile-fitted intensity IPR            (SelectIcolFlag=-1)
  //   COMBINE      weighted mean after combining partials  (SelectIcolFlag=+1)
  {
    SelectIcolFlag = IcolFlag;
    IpowerComb = Ipower;
    if (IcolFlag > 0) {
      SelectIcolFlag = 1;
      imid = Imid;
    }
  }
  //------------------------------------------------------------
  void col_controls::SetupColSelection(const int col_Ipr)
  // Set up observation_part class (static method) with current flags
  // conditional on Ipr present or absent (col_Ipr == 0 absent)
  {
    if (col_Ipr <= 0) {
      SelectIcolFlag = 0;
      prfpresent = false;
    } else {
      prfpresent = true; // both profile-fitted and summation intensities
    }
    int Iflag = SelectIcolFlag;
    SelectI::SetIcolFlag(Iflag, imid, IpowerComb);
    SelectI::SetIprPresent(prfpresent);
  }
  //------------------------------------------------------------
  std::string RejectFlags::formatReject2Policy() const
  {
    //  enum Reject2Policy {REJECT, KEEP, REJECTLARGER, REJECTSMALLER};
    if (rej2policy == REJECT) {
      return "REJECT";
    } else if (rej2policy == KEEP) {
      return "KEEP";
    } else if (rej2policy == REJECTLARGER) {
      return "REJECTLARGER";
    } else if (rej2policy == REJECTSMALLER) {
      return "REJECTSMALLER";
    }
    return "";
  }
  //------------------------------------------------------------
  std::string RejectFlags::format() const
  {
    return "Reflections measured 3 or more times: "+
      clipper::String(sdrej,3,3)+
      " maximum deviation from weighted mean of all other observations\n"+
      "Reflections measured twice: "+
      clipper::String(sdrej2,3,3)+" maximum deviation from weighted mean\n"+
      "   Policy for deviant reflections measured twice: "+
      formatReject2Policy()+"\n";
  }
  //------------------------------------------------------------
  OutlierControl::OutlierControl()
  // Set defaults
  {
    combine = true;  // default combine
    reject = RejectFlags(6.0, 6.0, RejectFlags::KEEP);
    int ndatasets = 1;
    rejectanom.assign(ndatasets, RejectFlags(9.0, 9.0, RejectFlags::KEEP));
    anomreject = true;
    emaxtest.init(10.0);
  }
  //------------------------------------------------------------
  OutlierControl::OutlierControl(const int& Ndatasets)
  // Set defaults, if Ndatasets <= 0 clear anomalous flags
  {
    combine = true;  // default combine
    reject = RejectFlags(6.0, 6.0, RejectFlags::KEEP);
    emaxtest.init(10.0);
    SetNdatasets(Ndatasets);
  }
  //------------------------------------------------------------
  void OutlierControl::SetNdatasets(const int& Ndatasets)
  // Copy rejectanom from 1st dataset or set if no first dataset
  // if Ndatasets <= 0 clear anomalous flags
  {
    if (Ndatasets <= 0) {
      anomreject = false;
      rejectanom.clear();
      return;
    }
    RejectFlags rejflags(9.0, 9.0, RejectFlags::KEEP);
    if (rejectanom.size() > 0) {
      rejflags = rejectanom.at(0);
    }
    rejectanom.assign(Ndatasets, rejflags);
    anomreject = true;
  }
  //------------------------------------------------------------
  // rejection criteria, within I+, I-  or between I+ & I-
  void OutlierControl::SetReject(const RejectFlags& flags,
                                 const AnomalousClass& selclass, const int& dts_index)
  // Set
  // if selclass == BOTH, dts_index may be = -1 for all datasets, in which case set all
  {
    if (selclass == BOTH) {
      if (dts_index >= 0) {
        rejectanom.at(dts_index) = flags;
      } else {
        int ndts = rejectanom.size();
        rejectanom.assign(ndts, flags); // set all datasets
      }
    } else {
      reject = flags;
    }
    anomreject = false;
    for (size_t id=0;id<rejectanom.size();++id) {
      if (rejectanom[id].sdrej > 0.0) anomreject = true;
    }
  }
  //------------------------------------------------------------
  RejectFlags OutlierControl::Reject(const AnomalousClass& selclass, const int& dts_index) const
  // Get
  {
    if (selclass == BOTH) {
      if (dts_index > int(rejectanom.size())-1) {
        std::cout << "OutlierControl::Reject " << dts_index <<" " << rejectanom.size() <<"\n"; //^
      }
      return rejectanom.at(Max(0,dts_index));
    } else {return reject;}
  }
  //------------------------------------------------------------
  void OutlierControl::SetEmax(const float& Emax) //!< set Emax (acentric)
  {
    emaxtest.init(Emax);
  }
  //------------------------------------------------------------
  RefineControl::RefineControl() {
    bfgs = true;
    ncyc1 = 2;
    ncycles = 10;
    converge = 0.3;
    iovsdmin = -3.0;
    e2min = 0.8;
    e2max = 5.0;
    nprocs = 1;
    maxprocs = 1;
  }
  //------------------------------------------------------------
  // Set number of processors from number or fraction (<1)
  void RefineControl::SetNprocs(const float& fproc)
  {
    // Ignore if OpenMP not enabled
    nprocs = 1;
#if _OPENMP
    //  Find out how many we are allowed
    maxprocs = 0;
    //if (getenv("OMP_NUM_THREADS") != NULL) {
    //  std::stringstream(std::string(getenv("OMP_NUM_THREADS"))) >> maxprocs;
    //}
    maxprocs = omp_get_num_procs();
    if (fproc > 0.99) {
      // Explicitly set
      nprocs = Min(Nint(fproc), maxprocs); // reset to maximum if greater
    } else if (fproc >= 0.0) {
      // fraction of maximum, at least 1
      nprocs = Max(1, Nint(fproc*maxprocs));
    } else if (fproc < 0.0) {
      nprocs = -1;  // set later (AUTO)
    }
#endif
  }
  //------------------------------------------------------------
  std::string RefineControl::format() const
  {
    std::string s;
    if (nprocs == 1) {
      s = "Refinement stages will use a single processor\n";
    } else if (nprocs > 1) {
      s = FormatOutput::logTabPrintf(0,
   "Number of processors used for refinement stages = %2d of maximum %2d\n",
                                     nprocs, maxprocs);
    } else {
      s =
     "Number of processors for refinement stages will be determined from number of observations\n";
    }
    return s;
  }
  //------------------------------------------------------------
  AnomalousControl::AnomalousControl()
  // Set defaults FIXME allow user input
  {
    Anomalous = false;   // true if "anomalous on"
    AnomalousSDcorr = false;   // true to separate I+ & I- for SD correction

    // At present, anomalous scattering is considered to be present if any one of
    // the following is true (defaults in brackets):
    //  1) Anomplot slope > anomslopethreshold (1.3)
    //  2) CCanom > anomCCthreshold (0.3) in more than anomNbinthreshold bins (2)
    //  3) RCRanom > anomRCRthreshold (1.3) in more than anomNbinthreshold bins (2)
    anomslopethreshold = 1.3;
    anomCCthreshold = 0.3;
    anomRCRthreshold = 1.3;
    anomNbinthreshold = 2;
  }
  //------------------------------------------------------------
}
