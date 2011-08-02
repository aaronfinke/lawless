// controls.cpp

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
  }
  //------------------------------------------------------------
  void partial_controls::Clear()
  {
    nrejgap = 0;
    nrejfractionsmall = 0;
    nrejfractionlarge = 0;
  }
  //------------------------------------------------------------
  std::string partial_controls::format() const
  {
    std::string s;
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
    emaxtest.init(10.0);
  }
  //------------------------------------------------------------
  OutlierControl::OutlierControl(const int& Ndatasets)
  // Set defaults
  {
    combine = true;  // default combine
    reject = RejectFlags(6.0, 6.0, RejectFlags::KEEP);
    int ndatasets = Ndatasets;
    rejectanom.assign(ndatasets, RejectFlags(9.0, 9.0, RejectFlags::KEEP));
    emaxtest.init(10.0);
  }
  //------------------------------------------------------------
  void OutlierControl::SetNdatasets(const int& Ndatasets)
  // Set defaults
  {
    rejectanom.assign(Ndatasets, RejectFlags(9.0, 9.0, RejectFlags::KEEP));
  }
  //------------------------------------------------------------
  // rejection criteria, within I+, I-  or between I+ & I-
  RejectFlags& OutlierControl::Reject(const AnomalousClass& selclass, const int& dts_index)
  // Set
  // dts_index may be = -1 for all datasets, in which case use 0
  {
    if (selclass == BOTH) {return rejectanom.at(Max(0,dts_index));}
    else {return reject;}
  }
  //------------------------------------------------------------
  RejectFlags OutlierControl::Reject(const AnomalousClass& selclass, const int& dts_index) const
  // Get
  {
    if (selclass == BOTH) {return rejectanom.at(Max(0,dts_index));}
    else {return reject;}
  }
  //------------------------------------------------------------
  bool OutlierControl::Anom() const
  // return true if rejection is set between I+ & I- for all datasets
  {
    for (size_t id=0;id<rejectanom.size();++id) {
      if (rejectanom[id].sdrej <= 0.0) return false;
    }
    return true;
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
  }
  //------------------------------------------------------------
}
