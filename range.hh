// range.hh
//
// Binning & range classes

#ifndef SCALA_RANGE
#define SCALA_RANGE

#include "clipper/clipper.h"

#include "util.hh"

typedef std::pair<float,float> RPair;

namespace scala
{
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
	  const bool& Ascending=true, const int& Nbin=0);

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

    // format
    std::string format() const;

  private:
    double first_, last_;
    mutable int Nbin_;
    mutable double width;
    mutable double tolerance; // fraction of bin width as tolerance
    bool ascending;  // false if negative width is allowed 

    void init();
    void CheckWidth() const;
  };

  

  class ResoRange : public Range
  {
  public:
    ResoRange();
    // Construct from low, high in A, number of observations
    //   (not reflections)
    ResoRange(const float& lowreso, const float& hireso,
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
    void SetRange(const float& lowreso, const float& hireso);
    void SetRange(const float& lowreso, const float& hireso,
		  const int& Nobs);
    void SetRange(const int& Nobs);

    // Unconditionally set number of bins
    void  SetNbins(const int& NumBin);
    int Nbins() const {return Nbin;}

    // Force width irrespective of Nobservations
    void SetWidth(const float& width);

    // Return limits
    float ResLow() const;   // in A
    float ResHigh() const;  // in A
    float SResLow() const;  // in 1/d^2 = 4(sin theta/lambda)^2
    float SResHigh() const; //
    Range RRange() const;    // in A
    Range SRange() const;    // in 1/d^2 = 4(sin theta/lambda)^2

    // Middle of bin (in A)
    float middleA(const int& bin) const;
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

  private:
    static const float LowDef;  // Default low resolution
    static const float HiDef;   //         high

    bool set;  // true if range explicitly set

    int Nbin;

    float LowReso, HiReso;
    float sSqrmin, sSqrmax;

    int MinNbin, MaxNbin; 
    int MinNrefBin, MaxNrefBin;
    int Nobservations;
    float delta_sSqr;

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

    // Return min or max
    int min() const;
    int max() const;

    // Return true if in range 
    bool InRange(const int& value) const;

  private:
    int min_, max_;

    void init();
  };

}
#endif
