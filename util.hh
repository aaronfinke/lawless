#ifndef PHIL_UTIL
#define PHIL_UTIL
// Simple utility functions, in unnamed namespace
// Self-contained
//   (ie no .cpp file, everything is in this header file)

#include <cmath>
#include <climits>
#include <vector>
#include <cstdlib>

#include <clipper/clipper.h>

namespace 
{
  template<class T> inline static int Nint( const T& a ) { return int( rint( a ) ); }
//--------------------------------------------------------------
    //! max
  //    template<class T> inline static T Max(const T& a, const T& b)
    template<class T1, class T2> inline static T1 Max(const T1& a, const T2& b)
	{ return (a > b) ? a : b; }
//--------------------------------------------------------------
    //! min
  //    template<class T> inline static T Min(const T& a, const T& b)
    template<class T1, class T2> inline static T1 Min(const T1& a, const T2& b)
      { return (a < b) ? a : b; }
//--------------------------------------------------------------
  // Close(a,b[,tol])  true if a == b within tolerance
  template<class T1, class T2> inline static bool
  Close(const T1& a, const T2& b,
	const T1& tol=1.0e-6)
  { return std::abs(a-b)<=tol; }
//--------------------------------------------------------------
  template<class T> inline bool IsInList(const std::vector<T>& list,
					 const T& item)
    //returns true if item is in list
  {
    return (std::find(list.begin(),list.end(),item) != list.end());	    
  }
//--------------------------------------------------------------
  template<class T> bool EqualList(const std::vector<T>& list1,
				   const std::vector<T>& list2)
    // Returns true if lists contain the same items, in any order
  {
    if (list1.size() != list2.size()) return false;
    for (int i=0;i<list2.size();i++)
      {
	if ( ! IsInList<T>(list1,list2[i])) return false;
      }
    return true;
  }
  //--------------------------------------------------------------
  int LargestCommonFactor(const clipper::Vec3<int>& vec)
  // find largest common factor of 3-vector
  {
    int small = INT_MAX;
    for (int i=0;i<3;++i) {
      if (std::abs(vec[i]) < small) small = std::abs(vec[i]); // smallest
    }
    // Test all factors from small to 1
    for (int jf=small;jf>=1;--jf) {
      bool found = true;
      for (int i=0;i<3;++i) {
	if (vec[i]%jf != 0) {found = false; break;}
      }
      if (found) {
	return jf;
      }
    }
    // Not found
    return 1;
  }
//--------------------------------------------------------------
}


#endif
