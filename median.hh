// median.hh

#ifndef MEDIAN_HEADER
#define MEDIAN_HEADER

#include <algorithm>

//======================================================================
template<class T> class Median {
  //! Median, interquartile range, etc
  public:
  Median(){}

  // Construct from list, sort unless sorted = true
  Median(const std::vector<T>& f, 
	 const bool& sorted=false);

  // - - - -
  void init(const std::vector<T>& f, 
	    const bool& sorted=false);

  // - - - -
  T median() const {return valueatfraction(0.5);}

  // - - - -
  T interquartilerange() const {
    return std::abs(valueatfraction(0.75) - valueatfraction(0.25));
  }

  T valueatfraction(const double& frac) const;

private:
  std::vector<T> v;

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
  if (!sorted) {
    std::sort(v.begin(), v.end());
  }
}
//--------------------------------------------------------------
template <class T> T Median<T>::valueatfraction(const double& frac) const
// get score at fraction, eg 0.5 for median
{
  T score = 0.0;
  int n = v.size();
  if (n == 0) {
    return 0.0;
  }
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
  return 0.5*(v[i1] + v[i1+1]);
}

#endif
