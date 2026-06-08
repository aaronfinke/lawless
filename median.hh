// median.hh

#ifndef MEDIAN_HEADER
#define MEDIAN_HEADER

#include <algorithm>

//======================================================================
template<class T> class Median {
  //! Median, interquartile range, etc
  public:
  Median():sorted_(false), count_(0){}

  // Construct from list, sort unless sorted = true
  Median(const std::vector<T>& f, 
	 const bool& sorted=false);

  // - - - -
  void init(const std::vector<T>& f, 
	    const bool& sorted=false);

  void clear() {v.clear(); count_=0;}

  // - - - -
  void add(const T& val)
  {v.push_back(val);}

  // - - - -
  void sort();

  // - - - -
  T median() {return valueatfraction(0.5);}

  // - - - -
  T interquartilerange() {
    return std::abs(valueatfraction(0.75) - valueatfraction(0.25));
  }

  T valueatfraction(const double& frac);

  // Mean of values between lowerfrac and upperfrac, set count
  T meanofrange(const double& upperfrac, const double& lowerfrac=0.0);

  // MAD Median Absolute Deviation
  T MAD();

  // Count
  int count() const {return count_;}

private:
  std::vector<T> v;
  bool sorted_;

  // Values valid after meanofrange
  int count_;   // = 0 before meanofrange called
};
//--------------------------------------------------------------
template <class T> Median<T>::Median(const std::vector<T>& f, 
	       const bool& sorted)
{init(f, sorted);}
//--------------------------------------------------------------
template <class T> void Median<T>::init(const std::vector<T>& f, 
			    const bool& sorted)
{
  v = f;
  count_ = 0;
  if (!sorted) {
    std::sort(v.begin(), v.end());
  }
  sorted_ = true;
}
//--------------------------------------------------------------
template <class T> void Median<T>::sort()
{
  std::sort(v.begin(), v.end());
  sorted_ = true;
}
//-------------------------------------------- ------------------
template <class T> T Median<T>::valueatfraction(const double& frac)
// get score at fraction, eg 0.5 for median
{
  int n = v.size();
  if (n == 0) {
    return 0.0;
  }
  if (!sorted_) {sort();}
  T score = 0.0;
  T rindex = frac * (n-1);
  int i1 = int(rindex); // lower index
  T remainder = rindex - i1;
  T tolerance = 0.5*frac;
  
  if (std::abs(remainder) < tolerance) {
    // on i1
    return v[i1];
  }
  if (i1+1 >= n) {
    return v[i1];
  }
  // Linear interpolation
  return v[i1]*(1.0-remainder) + v[i1+1]*remainder;
}
//-------------------------------------------- ------------------
// Mean of values between lowerfrac and upperfrac
template <class T> T Median<T>::meanofrange(const double& upperfrac, const double& lowerfrac)
{
  int n = v.size();
  if (n == 0) {
    return 0.0;
  }
  if (!sorted_) {sort();}

  T rindex = lowerfrac * (n-1);
  int i1 = int(rindex); // lower index
  rindex = upperfrac * (n-1);
  int i2 = int(rindex); // upper index
  
  T mean = 0.0;
  count_ = 0;
  for (int i=i1;i<=i2;++i) {
    mean += v[i];
    count_++;
  }
  mean = mean/T(i2-i1+1);
  return mean;
}
//-------------------------------------------- ------------------
// MAD Median Absolute Deviation
template <class T> T Median<T>::MAD()
{
  int n = v.size();
  if (n == 0) {
    return 0.0;
  }
  if (!sorted_) {sort();}

  T med = median();
  Median<T> deviations;
  for (int i=0;i<n;++i) {
    // | v[i] - median |
    deviations.add(std::abs(v[i] - med));
  }
  // Median of absolute deviaitons from median
  T mad = deviations.median();
  return mad;
}

#endif
