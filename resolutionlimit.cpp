// resolutionlimit.cpp
//
// Class to determine resolution "limit"
//

#include "resolutionlimit.hh"

#define ASSERT assert
#include <assert.h>

namespace scala {
  // ------------------------------------------------------------
  ResolutionLimit::ResolutionLimit(const std::vector<double> score,
				   const ResoRange& ResRange,
				   const double& Limit)
    : status(-2)
  //! construct from score list, resolution range and minimum score
  {
    init(score, ResRange, Limit);
  }
  // ------------------------------------------------------------
  void ResolutionLimit::init(const std::vector<double> score,
			     const ResoRange& ResRange,
			     const double& Limit) 
  //! initialise from score list, resolution range and minimum score
  {
    ASSERT (int(score.size()) == ResRange.Nbins());
    //! A high resolution limit is determined from a score list for
    //! each resolution bin, by finding the point at which the score
    //! drops below the limit. If the score in the highest resolution
    //! bin is above the limit, then the maximum resolution is used.
    //! Linear interpolation is used between bin values
    limit = Limit;
    if (score.back() >= limit) {
      highres = ResRange.ResHigh();
      status = +1;
      return;
    }
    int nbins = ResRange.Nbins();
    int i1 = -1;
    for (int i=nbins-1;i>=0;--i) { // loop backwards
      if (score[i] >= limit) {
	i1 = i;
	break;
      }
    }
    if (i1 < 0) {
      // All bins below limit, use lowest bin & set status
      highres = ResRange.boundsA(0).second;
      status = -1;
      return;
    }
    // linear interpolate on 1/d^2 between bins i1 and i1+1
    ASSERT (score[i1+1] <= score[i1]);  // just checking
    double den = score[i1]-score[i1+1];
    double f = 1.0;
    if (den > 0.0) { // trap for divide by 0
      f = (score[i1]-limit)/den;
    }
    highres = ResRange.middle(i1) +
      f * (ResRange.middle(i1+1) - ResRange.middle(i1));
    if (highres > 0.0) highres = 1.0/sqrt(highres);
    status = 0;
  }
  // ------------------------------------------------------------
  ResolutionLimit::ResolutionLimit(const std::vector<MeanSD> mnsd,
				   const ResoRange& ResRange,
				   const double& Limit)
  //! construct from MeanSD list (score = Mean), resolution range and minimum score
  {
    init (mnsd, ResRange, Limit);
  }
  // ------------------------------------------------------------
  void ResolutionLimit::init(const std::vector<MeanSD> mnsd,
			     const ResoRange& ResRange,
			     const double& Limit)
  //! initialise from MeanSD list (score = Mean), resolution range and minimum score
  {
    std::vector<double> score(mnsd.size());
    for (size_t i=0;i<mnsd.size();++i) {
      score[i] = mnsd[i].Mean();
    }
    init (score, ResRange, Limit);
  }
}
