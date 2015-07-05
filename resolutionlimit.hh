// resolutionlimit.hh
//
// Class to determine resolution "limit"
//

#ifndef RESOLUTIONLIMIT_HEADER
#define RESOLUTIONLIMIT_HEADER

#include <vector>

#include "range.hh"
#include "scala_util.hh"
#include "radialfunction.hh"
#include "fitresolutiondata.hh"


namespace scala {
  //---------------------------------------------------------------
  class ResolutionLimit {
    //! A high resolution limit is determined from a score list for
    //! each resolution bin, by finding the point at which the score
    //! drops below the limit. If the score in the highest resolution
    //! bin is above the limit, then the maximum resolution is used.
    //! Linear interpolation is used between bin values
    //
  public:
    ResolutionLimit() : limit(0.0), highres(0.0), status(-3), 
			sufficientdata(false), nrej(0), fittype(NONE) {}

    // Fittype: NONE   no curve-fitting (eg for mnI/sigI)
    //          TANH   curve fit to tanh function ~= (0.5*(1-tanh((s-d0)r)
    //          LINEAR straight-line fit (eg for all ~0 or negative)
    enum FitType {NONE, TANH, LINEAR};

    //! construct from score list, resolution range and minimum score
    ResolutionLimit(const std::vector<double> score,
		    const ResoRange& ResRange,
		    const double& Limit,
		    const FitType& Fittype);
    //! construct from score list, resolution range and minimum score
    void init(const std::vector<double> score,
	      const ResoRange& ResRange,
	      const double& Limit,
	      const FitType& Fittype);
    //! construct from MeanSD list (score = Mean), resolution range and minimum score
    ResolutionLimit(const std::vector<MeanSD> mnsd,
		    const ResoRange& ResRange,
		    const double& Limit,
		    const FitType& Fittype);
    //! initialise from MeanSD list (score = Mean), resolution range and minimum score
    void init(const std::vector<MeanSD> mnsd,
	      const ResoRange& ResRange,
	      const double& Limit,
	      const FitType& Fittype);

    //! return high resolution limit in A
    double HighResolution() const {return highres;}

    //! return score limit used
    double Limit() const {return limit;}

    //! return fitted value at s = 1/d^2
    double fitvalue(const double& s) const;

    // =  NONE, no function fit
    // =  TANH, fit radial tanh function (for CC(1/2))
    FitType Fittype() const  {return fittype;}

    //! status = -1 all scores below limit
    //!          +1 all scores above limit
    //!           0 intermediate limit set
    //!          -2 unset because of no data
    //!          -3 unset at all
    int Status() const { return status;}

    //! true if there is a any data
    bool valid() const
    {return status != -3;}

    //! true if there were enough data for a fit
    bool sufficientData() const {return sufficientdata;}

    // anomalous == true for assessment of anomalous signal (changes wording)
    std::string format(const bool& anomalous) const;

    std::string formatparameters() const;

  private:
    double limit;   // score limit used
    double highres; // high resolution limit, = 0.0 if undefined
    // status = -1 all scores below limit
    //          +1 all scores above limit
    //           0 intermediate limit set
    //          -2 unset because of no data
    //          -3 unset at all 
    int status;
    bool sufficientdata;  // true if there were enough data for a fit

    mutable int nrej;

    // =  NONE, no function fit
    // =  TANH, fit radial tanh function (for CC(1/2))
    // =  LINEAR, fit straight-line (for CCanom)
    FitType fittype;

    RadialTanhFunction radialfunction;

    double slope, intercept;  // for straight-line fit

    double fit(const std::vector<double>& score,
	       const ResoRange& ResRange);

    int rejectoutliers(std::vector<ResolutionData>& data,
		       const double& reject) const;

    // return true if there are not enough data points to determine
    bool testinsufficientdata(const std::vector<double> score) const;

  };
  //---------------------------------------------------------------
} // namespace scala 
#endif
