// histogram.cpp


#include "histogram.hh" 
namespace scala {
// ------------------------------------------------------------
Histogram::Histogram(const Range& range)  //!< construct from range & Nbins
{
  init(range);
}
// ------------------------------------------------------------
void Histogram::init(const Range& range)  //!< initialise from range & Nbins
{
  limits = range;
  actualrange.clear();
  nbins = limits.Nbins();
  counts.assign(nbins,0);
}
// ------------------------------------------------------------
void Histogram::SetNbins(const int& Nbins)  //!< store number of bins
{
  nbins = Nbins;
  counts.assign(nbins,0);
}
 // ------------------------------------------------------------
void Histogram::add(const double& val)   // add value into histogram
{
  int m = limits.tbin(val);
  if (m < 0) {
    if (val < limits.min()) {
      noutside.first++;
    } else if (val > limits.max()) {
      noutside.second++;
    }
  } else {
    counts[m]++;
  }
  actualrange.update(val);
}
// ------------------------------------------------------------
//! return total counts
int Histogram::total() const
{
  int total = 0;
  for (size_t i=0; i<counts.size(); i++) { 
    total += counts[i];
  }
  return total;
}
// ------------------------------------------------------------
}
