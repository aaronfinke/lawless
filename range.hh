// range.hh
//
// Binning & range classes

#ifndef SCALA_RANGE
#define SCALA_RANGE

#include <cstdlib>

#include "clipper/clipper.h"

#include "util.hh"

typedef std::pair<float,float> RPair;

namespace scala
{
  class IntRange;
  //========================================================================
  class Range
  {
  public:
    // Default constructor, set min & max to extremes ready
    // for update, Nbin:=0
    Range();
    // Rfirst & Rlast will be swapped if necessary so that Rfirst < Rlast
    //  unless Ascending == false
    Range(const double& Rfirst, const double& Rlast,
	  const bool& Ascending=true, const int& Nbin=1);

    Range(const IntRange& intrange); // construct from an IntRange

    virtual ~Range(){}
    
    // Allow negative range (first > last), set ascending = false
    void AllowDescending();

    // Rfirst & Rlast will be swapped if necessary so that Rfirst < Rlast
    //  unless ascending == false
    void SetRange(const double& Rfirst, const double& Rlast,
		  const bool& Ascending=true, const int& Nbin=0);
    void SetNbin(const int& Nbin);

    // Clear min & max prior to update
    void clear();
    // Update min & max  (min is first_!)
    virtual void update(const double& value);

    // Set or return first or last
    double& first() {return first_;} // set
    double first() const {return first_;} // get
    double& last() {return last_;} // set
    double last() const {return last_;} // get
    // Return min or max
    double min() const {return Min(first_, last_);} // get
    double max() const {return Max(first_, last_);} // get

    // Return absolute value of range (last - first)
    double AbsRange() const {return std::abs(last_ - first_);}

    // Return bin number in range 0->(Nbin-1) (reset to edge bin if outside)
    int bin(const double& value) const;
    // Return bin number in range 0->(Nbin-1) or -1 (reject) if outside range
    int tbin(const double& value) const;

    // Return number of bins
    int Nbins() const {return Nbin_;}
    // Lower & upper bounds of bin
    RPair bounds(const int& bin) const;
    // Middle of bin
    float middle(const int& bin) const;

    // Return true if valid (ie set)
    bool Valid() const {return valid;}

    // format
    std::string format() const;

    // Returns maximum range
    Range MaxRange(const Range& other) const;

  private:
    double first_, last_;
    mutable int Nbin_;
    mutable double width;
    mutable double tolerance; // fraction of bin width as tolerance
    bool ascending;  // false if negative width is allowed 
    bool valid; // true if set

    void init();
    void CheckWidth() const;
  };

  

  class ResoRange : public Range
  {
  public:
    ResoRange();
    // Construct from low, high in A, number of observations
    //   (not reflections)
    ResoRange(const double& lowreso, const double& hireso,
	      const int& Nobs=0);
    // construct from resolution range in s = 1/d^2
    ResoRange(const Range& range);
    // Initialise from Range in 1/d^2
    void init_range(const Range& range);

    // true if explicitly set
    bool isSet() const {return set;}
    // Mark as set (from update) and initialise
    void Set();

    // Reset resolution range
    void SetRange(const double& lowreso, const double& hireso);
    void SetRange(const double& lowreso, const double& hireso,
		  const int& Nobs);
    void SetRange(const int& Nobs);

    // Extend range by small tolerance
    void ExtendRange();

    // Unconditionally set number of bins
    void  SetNbins(const int& NumBin);
    int Nbins() const {return Nbin;}

    // Force bin width irrespective of Nobservations
    void SetWidth(const double& width);
    // Return bin width
    double Width() const {return delta_sSqr;}

    // Return limits
    double ResLow() const;   // in A
    double ResHigh() const;  // in A
    double SResLow() const;  // in 1/d^2 = 4(sin theta/lambda)^2
    double SResHigh() const; //
    Range RRange() const;    // in A
    Range SRange() const;    // in 1/d^2 = 4(sin theta/lambda)^2

    // Middle of bin (in A)
    double middleA(const int& bin) const;
    // Middle of bin (in 1/d^2), use middle()
    // limits of bin (in A)
    RPair boundsA(const int& bin) const;
    // limits of bin (in 1/d^2)
    RPair boundsS(const int& bin) const;
    // Range of bin
    ResoRange BinRange(const int& bin) const;
    // Returns maximum range
    ResoRange MaxRange(const ResoRange& other) const;
    // Returns minimum range
    ResoRange MinRange(const ResoRange& other) const;

    std::string format() const;

  private:
    static const double LowDef;  // Default low resolution
    static const double HiDef;   //         high

    bool set;  // true if range explicitly set

    int Nbin;

    double LowReso, HiReso;
    double sSqrmin, sSqrmax;

    int MinNbin, MaxNbin; 
    int MinNrefBin, MaxNrefBin;
    int Nobservations;
    double delta_sSqr;

    void init();
    void init_reso();
    void init_bins();

  };
  //========================================================================
  class IntRange
  {
    // Range class for integers  (eg batch numbers)
    // no binning
  public:
    // Default constructor, set min & max to extremes ready
    // for update
    IntRange();
    // Imin & Imax will be swapped if necessary so that Imin < Imax
    IntRange(const int& Imin, const int& Imax);

    // Imin & Imax will be swapped if necessary so that Imin < Imax
    void SetRange(const int& Imin, const int& Imax);

    // Clear min & max prior to update
    void clear();
    // Update min & max (in loop)
    void update(const int& value);
    void update(const IntRange& range);

    // Return min or max
    int min() const;
    int max() const;

    // Return true if in range 
    bool InRange(const int& value) const;

    // Returns maximum range
    IntRange MaxRange(const IntRange& other) const;

    int midrange() const {return (min_+max_)/2;}

    // Return absolute value of range (max - min)
    int AbsRange() const {return std::abs(max_ - min_);}

    // test equality
    bool operator == (const IntRange& other)
    {return ((min_==other.min_) && (max_==other.max_));}
    bool operator != (const IntRange& other)
    {return ((min_!=other.min_) || (max_!=other.max_));}

    // true if ranges overlap
    bool Overlap(const IntRange& other) const;

  private:
    int min_, max_;

    void init();
  };

}
#endif
