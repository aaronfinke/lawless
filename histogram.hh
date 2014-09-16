// histogram.hh
//
// A histogram class
//

#ifndef HISTOGRAM_HEADER
#define HISTOGRAM_HEADER

#include <vector>
#include "range.hh"

namespace scala {
class Histogram {
public:
  Histogram(){} //!< empty constructor
 
  Histogram(const Range& range);  //!< construct from range & Nbins
  void init(const Range& range);  //!< initialise from range & Nbins

  void SetNbins(const int& Nbins);  //!< store number of bins

  void add(const double& val);   // add value into histogram

  //! return counts in each bin
  std::vector<int> Counts() const {return counts;}
  //! return range
  Range HRange() const {return limits;}
  //! return numbers below and above limits
  std::pair<int, int> Noutside() const {return noutside;}
  //! return actual range
  Range ActualRange() const {return actualrange;}
  //! return total counts
  int total() const;

private:
  Range limits;
  Range actualrange;
  int nbins;
  std::vector<int> counts;
  std::pair<int, int> noutside;  // numbers below and above limits
};
}
#endif
