// sdanalysis.hh

#ifndef SDANALYSIS_HEADER
#define SDANALYSIS_HEADER

#include <vector>

#include "sdmodel.hh"
#include "score_datatypes.hh"
#include "intensitybin.hh"
#include "selectedobservations.hh"
#include "file_util.hh" //^^ dumping option

namespace scala
{
  // ------------------------------------------------------------
  class SDparameterGroup
  //! a subclass of observations for SD analysis, and indices of parameters
  //  For SD correction, observations are grouped by run (optionally) and by
  //  full/partial
  //  This class handles one of these groups, and the sums over all observations
  //  in that bin class needed for analysis by intensity
  {
  public:
    SDparameterGroup(){}
    void init(const int& Numintbins, const int& Nparam);

    // Add in delta for intensity bin
    void AddDelta(const int& intBin, const float& delta);

    // Add in derivative for intensity bin and full/partial
    // kpl is local parameter number
    void AddDerivative(const int& intBin, const int& kpl,
		       const float& dddp, const float& delta);
    
    //! return sd(delta), = 0 if none
    std::vector<double> SDdelta() const;

    //! return MeanSD(delta) for each intensity bin
    std::vector<MeanSD> MeanDelta() const;

    // return d(sigma(delta[inb]))/dp[k] for intensity bins inb
    // outer vector length nparams (k), inner length numintbins (inb)
    std::vector<std::vector<double> > dSigmaDeltaDp() const;

    int Nparam() const {return nparams;}

    // to sum totals
    SDparameterGroup& operator += (const SDparameterGroup& other);
    friend SDparameterGroup& operator + (const SDparameterGroup& a, const SDparameterGroup& b);

  private:
    int numintbins; // number of intensity bins
    std::vector<MeanSD> mndelta; // by intensitybin
    int nparams;    // number of (local) parameters for partial derivative components, 0 none

    // in the following, outer vector length nparams, inner length numintbins
    std::vector<std::vector<double> > sumddeltadp;
    std::vector<std::vector<double> > sumddelta2dp;
    std::vector<std::vector<int> >  nsumddeltadp;
  };
  // ------------------------------------------------------------
  // ------------------------------------------------------------
  class SDanalysis {
  public:
    SDanalysis() : numintensitybins(0), dumpfile(NULL) {}
    // Construct and clear
    // Derivatives = false, no derivative components accumulated
    SDanalysis(const IntensityBin& Irange, const SDmodel& SDM,
	       const bool& Allsamerun, const bool& Derivatives=false);
    // also full/partial
    // reset and clear
    void init (const IntensityBin& Irange, const SDmodel& SDM,
	       const bool& Allsamerun, const bool& Derivatives=false);
    // also full/partial
    // just clear
    void clear();

    //^ dumping for analysis
    void SetDump(const std::string& dumpfilename);

    //! returns true if object is empty
    bool Empty() const; 

    //! return true if there are fulls and they were used (cf useflag)
    bool FullsUsed(const int& irun) const;

    //! return true if there are Partials and they were used (cf useflag)
    bool PartialsUsed(const int& irun) const;

    //! return true if allsame run
    bool AllSameRun() const {return allsamerun;}

    //! set maximum allowed delta, just to eliminate idiocies
    void SetDeltaLimit(const float& deltaLimit) {deltalimit = deltaLimit;}

    //! add in delta, for intensity bin intBin, run irun, full == false for partial
    bool AddDelta(const float& delta, const int& intBin,
		  const int& irun, const bool& full);  //!< returns true if OK

    // Add in delta1 contributions for intensity bin intBin
    //  delta1 = (Ihl - <Ih>others)/sqrt(SD(Ihl)^2 + SD(Ihothers)^2)
    //   where <Ih>others is the average over all observations of reflection h,
    //   excluding Ihl itself and SD(Ihothers) is its ESD
    //   cf AddSelobsDelta2 which uses delta2 definition
    void AddSelobsDelta(SelectedObservations& selobs,
			const int& intBin);

    // Add in delta2 contributions for intensity bin intBin
    //  delta2 = sqrt(n/n-1) (Ihl - <Ih>)/SD(Ihl)
    //   where <Ih> is the average over all observations of reflection h,
    //   including Ihl itself
    //   cf AddSelobsDelta which uses delta1 definition
    void AddSelobsDelta2(SelectedObservations& selobs,
			 const int& intBin);

    // Add in to sums for derivatives
    //  ddeltadp[iobs][k] is d(delta(iobs))/dp(k) for the iobs'th observation in selobs
    //  p(k) is the k'th parameter of nparams
    void AddDerivatives(SelectedObservations& selobs,
			const int& intBin,
			const std::vector <std::vector<double> >& ddeltadp);

    //    void PrintTable(const std::vector<Run>& runlist,
    //		    const SDmodel& SDM,
    //		    const int& datasetIndex, const PxdName& dataset_pxd,
    //		    phaser_io::Output& output);

    // Return one element
    MeanSD GetMeanSD(const int& intBin,
		  const int& irun, const bool& full) const;
    // Return one element
    MeanSD GetMeanSD(const int& intBin,
		  const int& irun, const int& fullpart) const;
    // Return vector for all intensity bins
    std::vector<MeanSD> GetMeanSD(const int& irun, const bool& full) const;

    //! return number of observations in each intensity bin, over all runs and full/partials
    std::vector<int> NumberinIntbins() const;

    //! return number of observations in each  class
    std::vector<int> NumberinClass() const;

    //! return total number of observations, fulls, partials
    std::pair<int,int> Number() const;

    //! return intensity range
    IntensityBin Irange() const {return irange;}

    //! return number of intensity bins
    int NumberIntensityBins() const {return numintensitybins;}

    //! return number of runs
    int NumberRuns() const {return numruns;}

    //! return number of runs used
    int NumberRunsUsed() const {return numrunsused;}

    //! return number of parameter classes
    int NumberOfParameterGroups() const {return numparametergroups;}

    //! return bin class for parameter class and intensity bin
    int BinClass(const int& jpc, const int& mint) const;

    //! return number of parameters in parameter class
    int Nparam(const int& paramclass) const;

    //! return index in complete parameter list to first parameter in class
    int IdxParam(const int& paramclass) const;

    //! SDs for each "bin class" jc
    std::vector<double> SDdelta() const;

    //! MeanSD for deltas for all bin classes
    std::vector<MeanSD> MeanDelta() const;

    // SDs and their derivatives w.r.t. parameters p(k) for
    // each "bin class" jc
    // A bin class jc is intensity/run/full|partial
    //   returns dsigDeldp  [jc][k] partial derivatives d(sigma(delta(jc)))/dp(k) for
    //              bin class jc, parameter k
    std::vector<std::vector<std::vector<double> > > Derivatives() const;
    //! corresponding numbers
    std::vector<std::vector<int> > NumberInDerivatives() const;

  private:
    std::string pxdname;
    int numintensitybins;
    IntensityBin irange;
    int numruns;   // actual number
    int numrunsused; // = 1 if allsamerun
    bool allsamerun;// true if all runs stored as one
    bool fulls;     // true if any fulls
    bool partials;  // true if any partials
    int nparams;  // over all classes

    // = 0 normal, both fulls & partials; < 0 treat full as partial; > 0 treat partial as full
    //   = +1 only fulls    = +2 few partials
    //   = -1 only partials = -2 few fulls
    std::vector<int> useflag;  // for each run

    std::vector<SDparameterGroup> sdparametergroup; // for each parameter group
    int numparametergroups; // number of parameter classes
    std::vector<int> grouprunnumber;  // run number for each parameter group
    std::vector<bool> groupfull;  // true if group is for fulls, false if for partials
    std::vector<int> rungroupnumberfulls;    // group number for fulls in each run
    std::vector<int> rungroupnumberpartials; // group number for partials in each run

    int nbinclasses;  // number of analysis bin classes, <= numrunsused*2 * numintensitybins

    // local parameter index for each global parameter
    std::vector<int> idxlocal;
    // index into global parameter list for 1st parameter of each parameter group
    std::vector<int> idxparametergroup;

    //! parameter group number for run and full/partial
    int ParameterGroup(const int& irun, const bool& full) const;

    float deltalimit;  // maximum allowed delta, just to eliminate idiocies

    FILE* dumpfile; 

  };
  // ------------------------------------------------------------
  void PrintSDanalysis(const SDanalysis& sdanal1, const SDanalysis& sdanal2,
		       const RejectFlags& rejflags,
		       const IntensityBin& Irange,
		       const std::vector<Run>& runlist,
		       const SDmodel& SDM,
		       const int& datasetIndex, const PxdName& dataset_pxd,
		       const bool& fullprint,
		       phaser_io::Output& output);
  // ------------------------------------------------------------
  void PrintSDanalysisTable(const RejectFlags& rejflags,
			    const IntensityBin& Irange,
			    const std::string& ttitle,
			    const int& fullpartial,
			    const bool& outer,
			    const std::vector<std::vector<MeanSD> >& msdanal,
			    phaser_io::Output& output);
  // ------------------------------------------------------------
}

#endif
