// resolutionlimit.hh
//
// Class to determine resolution "limit"
//

#ifndef RESOLUTIONLIMIT_HEADER
#define RESOLUTIONLIMIT_HEADER

#include <vector>

#include "range.hh"
#include "scala_util.hh"

namespace scala {
  class ResolutionLimit {
    //! A high resolution limit is determined from a score list for
    //! each resolution bin, by finding the point at which the score
    //! drops below the limit. If the score in the highest resolution
    //! bin is above the limit, then the maximum resolution is used.
    //! Linear interpolation is used between bin values
    //
  public:
    ResolutionLimit() : status(-2){}
    //! construct from score list, resolution range and minimum score
    ResolutionLimit(const std::vector<double> score,
		    const ResoRange& ResRange,
		    const double& Limit);
    //! construct from score list, resolution range and minimum score
    void init(const std::vector<double> score,
	      const ResoRange& ResRange,
	      const double& Limit);
    //! construct from MeanSD list (score = Mean), resolution range and minimum score
    ResolutionLimit(const std::vector<MeanSD> mnsd,
		    const ResoRange& ResRange,
		    const double& Limit);
    //! initialise from MeanSD list (score = Mean), resolution range and minimum score
    void init(const std::vector<MeanSD> mnsd,
	      const ResoRange& ResRange,
	      const double& Limit);

    //! return high resolution limit in A
    double HighResolution() const {return highres;}

    //! return score limit used
    double Limit() const {return limit;}

    //! status = -1 all scores below limit
    //!          +1 all scores above limit
    //!           0 intermediate limit set
    //!          -2 unset
    int Status() const { return status;}

  private:
    double limit;   // score limit used
    double highres; // high resolution limit
    // status = -1 all scores below limit
    //          +1 all scores above limit
    //           0 intermediate limit set
    //          -2 unset
    int status;
  };


} // namespace scala 
#endif
