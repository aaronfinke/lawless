// range.cpp

#include "range.hh"

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;

namespace scala 
{
  //--------------------------------------------------------------
  Range::Range()
    : first_(+1000000), last_(-1000000), Nbin_(0),
      width(0.0), tolerance(0.001), ascending(true)  {}
  //--------------------------------------------------------------
    // Rfirst & Rlast will be swapped if necessary so that Rfirst < Rlast
    //  unless ascending == false
  Range::Range(const double& Rfirst, const double& Rlast,
	       const bool& Ascending, const int& Nbin)
    : first_(Rfirst), last_(Rlast), Nbin_(Nbin), ascending(Ascending)
  {
    init();
  }
  //--------------------------------------------------------------
  void Range::SetRange(const double& Rfirst, const double& Rlast,
		       const bool& Ascending, const int& Nbin)
  {
    first_ = Rfirst;
    last_ = Rlast;
    Nbin_ = Nbin;
    ascending = Ascending;
    init();
  }
  //--------------------------------------------------------------
  void Range::init()
  {
    // Allow first > last if ascending false
    if (ascending) {
      if (last_ < first_) {
	double a = last_;
	last_ = first_;
	first_ = a;
      }}
    width = 0.0;
    tolerance = 0.0001;
  }
  //--------------------------------------------------------------
    // Allow negative range (first > last), set ascending = false
  void Range::AllowDescending()
  {
    ascending = false;
  }
  //--------------------------------------------------------------
  void Range::SetNbin(const int& Nbin)
  {
    Nbin_ = Nbin;
  }
  //--------------------------------------------------------------
  void Range::clear()
    // Reset first & last to extremes ready for update
  {
    first_ = +1.0e+10;
    last_ = -first_;
  }
  //--------------------------------------------------------------
  void Range::update(const double& value)
  {
    first_ = Min(first_, value);
    last_ = Max(last_, value);
  }
  //--------------------------------------------------------------
  void Range::CheckWidth() const
  {
    if (width == 0.0) {
      if (Nbin_ <= 0)
	clipper::Message::message(Message_fatal("Range: Nbin <= 0 "
						  +clipper::String(Nbin_)));
      width = (double(last_ - first_))/double(Nbin_);
      if (ascending && width <= 0.0) {
	clipper::Message::message(Message_fatal("Range: width <= 0  "

						+clipper::String(width)));
      } else if (width == 0.0) {
	// If width == 0, reset Nbin to 1
	Nbin_ = 1;
      }
    }
    if (width > 0.0) {tolerance = std::abs(tolerance);}
    else if (width < 0.0) {tolerance = -std::abs(tolerance);}
  }
  //--------------------------------------------------------------
  // Return bin number in range 0->(Nbin-1)
  int Range::bin(const double& value) const
  {
    CheckWidth();
    if (Nbin_ <= 1) return 0;
    return Max(Min(int((double(value-first_)/width) - tolerance), Nbin_-1),0);
  }
  //--------------------------------------------------------------
  // Return bin number in range 0->(Nbin-1) or -1 (reject) if outside range
  int Range::tbin(const double& value) const
  {
    CheckWidth();
    if (Nbin_ <= 1) return 0;
    int n = int((double(value-first_)/width) - tolerance);
    if (width > 0.0) {
      if (value < first_ || value > last_) n = -1;
    } else if (width < 0.0) {
      if (value > first_ || value < last_) n = -1;
    }
    return n;
  }
  //--------------------------------------------------------------
  RPair Range::bounds(const int& bin) const
  // Lower & upper bounds of bin
  {
    CheckWidth();
    return RPair(float(bin)*width+first_,float(bin+1)*width+first_);
  }
  //--------------------------------------------------------------
  // Middle of bin
  float Range::middle(const int& bin) const
  {
    CheckWidth();
    return (float(bin)+0.5)*width+first_;
  }
  //--------------------------------------------------------------
  // format
  std::string Range::format() const
  {
    CheckWidth();
    return "Range: first "+clipper::String(first_)+
      " last "+clipper::String(last_)+" Nbin "+clipper::String(Nbin_)+
      " width "+clipper::String(width)+"\n";
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  const float ResoRange::LowDef = 10000.;  // Default low resolution
  const float ResoRange::HiDef = 0.001;   //         high

  //--------------------------------------------------------------
  ResoRange::ResoRange()
    : Range(), Nobservations(0)
  {
    LowReso = LowDef;
    HiReso = HiDef;
    init_reso();
    init_bins();
    set = false;
  }
  //--------------------------------------------------------------
  ResoRange::ResoRange(const float& lowreso, const float& hireso,
		       const int& Nobs)
    : Range(), LowReso(lowreso), HiReso(hireso), Nobservations(Nobs)
  {
    init_reso();
    init_bins();
  }
  //--------------------------------------------------------------
  ResoRange::ResoRange(const Range& range)
    : Range(), Nobservations(0)
    // Construct from Range in 1/d^2
  {
    init_range(range);
  }
  //--------------------------------------------------------------
  void ResoRange::init_range(const Range& range)
    // Initialise from Range in 1/d^2
  {
    sSqrmin = range.min();
    sSqrmax = range.max();
    LowReso = LowDef;
    HiReso = HiDef;
    if (sSqrmin > 0.0)
      LowReso = 1./sqrt(range.min());
    if (sSqrmax > 0.0)
      HiReso = 1./sqrt(range.max());
    init_bins();
    init();
  }
  //--------------------------------------------------------------
  // Mark as set (from update)
  void ResoRange::Set() {
    sSqrmin = Range::min();
    sSqrmax = Range::max();
    set = true;
    init_range(*this);
  }
  //--------------------------------------------------------------
  void ResoRange::SetRange(const float& lowreso, const float& hireso)
  {
    LowReso = lowreso;
    HiReso = hireso;
    init_reso();
    init();
  }
  //--------------------------------------------------------------
  void ResoRange::SetRange(const float& lowreso, const float& hireso,
			   const int& Nobs)
  {
    Nobservations = Nobs;
    SetRange(lowreso, hireso);
  }
  //--------------------------------------------------------------
  void ResoRange::SetRange(const int& Nobs)
  {
    Nobservations = Nobs;
    init();
  }
  //--------------------------------------------------------------
  void ResoRange::SetWidth(const float& width)
  {
    // Force width irrespective of Nobservations
    delta_sSqr = width;
    Nobservations = 0;
    init();
  }
  //--------------------------------------------------------------
  void ResoRange::init_bins()
  {
    // Set up default guides
    MinNbin = 5;
    MaxNbin = 30;
    MinNrefBin = 200;
    MaxNrefBin = 3000;
    delta_sSqr = 0.012;
    init();
  }
  //--------------------------------------------------------------
  void  ResoRange::init_reso()
  {
    if (LowReso < HiReso)
      std::swap(LowReso, HiReso);
    
    if (LowReso <= 0.0) 
      LowReso = LowDef;
    sSqrmin = 1./(LowReso*LowReso);
    if (HiReso <= 0.0) 
      HiReso = HiDef;
    sSqrmax = 1./(HiReso*HiReso);
    set = true;
  }
  //--------------------------------------------------------------
  void  ResoRange::SetNbins(const int& NumBin)
  {
    Nbin = NumBin;
    Range::SetRange(sSqrmin, sSqrmax, true, Nbin);
    set = true;
  }
  //--------------------------------------------------------------
  void  ResoRange::init()
  {
    // try bin width
    int n = int((sSqrmax-sSqrmin)/delta_sSqr);

    // Tune by number of observations
    if (Nobservations > 0 && n > 0) {
      if (Nobservations/n < MinNrefBin) n = Max(1,Nobservations/MinNrefBin);
      if (Nobservations/n > MaxNrefBin) n = Nobservations/MaxNrefBin;
    }
  
    // Overall minimum & maximum
    Nbin = Min(Max(n, MinNbin), MaxNbin);
    Range::SetRange(sSqrmin, sSqrmax, true, Nbin);
    set = true;
  }
  //--------------------------------------------------------------
  float  ResoRange::ResLow() const
  {
    return 1./sqrt(Range::min());
  }
  //--------------------------------------------------------------
  float  ResoRange::ResHigh() const
  {
    return 1./sqrt(Range::max());
  }

  //--------------------------------------------------------------
  float  ResoRange::SResLow() const
  {
    return Range::min();
  }
  //--------------------------------------------------------------
  float  ResoRange::SResHigh() const
  {
    return Range::max();
  }
  //--------------------------------------------------------------
  Range ResoRange::RRange() const
  {
    return Range(1./sqrt(Range::max()),1./sqrt(Range::min()));
  }
  //--------------------------------------------------------------
  Range ResoRange::SRange() const
  {
    return Range(min(),max());
  }
  //--------------------------------------------------------------
  // Middle of bin (in A)
  float ResoRange::middleA(const int& bin) const
  {
    return 1./sqrt(Range::middle(bin));
  }
  //--------------------------------------------------------------
  // limits of bin (in A)
  RPair ResoRange::boundsA(const int& bin) const
  {
    return RPair(1./sqrt(Range::bounds(bin).first),
		 1./sqrt(Range::bounds(bin).second));
  }
  //--------------------------------------------------------------
  // limits of bin (in 1/d^2)
  RPair ResoRange::boundsS(const int& bin) const
  {
    return RPair(Range::bounds(bin).first,
		 Range::bounds(bin).second);
  }
  //--------------------------------------------------------------
  // Range of bin
  ResoRange ResoRange::BinRange(const int& bin) const
  {
    return ResoRange(1./sqrt(Range::bounds(bin).first),
		     1./sqrt(Range::bounds(bin).second));
  }
  //--------------------------------------------------------------
  // Returns maximum range
  ResoRange ResoRange::MaxRange(const ResoRange& other) const
  {
    int nobs = Nobservations + other.Nobservations;
    if (set) {
      if (other.set) {
	// Both set
	return ResoRange(Max(ResLow(), other.ResLow()),
			 Min(ResHigh(), other.ResHigh()),
			 nobs);
      } else {
	// other not set
	return *this;
      }
    } else {
      // this not set
      return other;  // if not set
    }
  }
  //--------------------------------------------------------------
  // Returns minumum range
  ResoRange ResoRange::MinRange(const ResoRange& other) const
  {
    int nobs = Nobservations + other.Nobservations;
    if (set) {
      if (other.set) {
	// Both set
	return ResoRange(Min(ResLow(), other.ResLow()),
			 Max(ResHigh(), other.ResHigh()),
			 nobs);
      } else {
	// other not set
	return *this;
      }
    } else {
      // this not set
      return other;  // if not set
    }
  }
  //--------------------------------------------------------------
  IntRange::IntRange() {clear();}
  //--------------------------------------------------------------
  IntRange::IntRange(const int& Imin, const int& Imax)
    : min_(Imin), max_(Imax)
  {
    init();
  }
  //--------------------------------------------------------------
  void IntRange::SetRange(const int& Imin, const int& Imax)
  {
    min_ = Imin;
    max_ = Imax;
    init();
  }
  //--------------------------------------------------------------
  void IntRange::init()
  {
    if (max_ < min_)
      {
	int a = max_;
	max_ = min_;
	min_ = a;
      }
  }
  //--------------------------------------------------------------
  void IntRange::clear()
    // Reset min & max to extremes ready for update
  {
    min_ = 400000000;
    max_ = -min_;
  }
  //--------------------------------------------------------------
  void IntRange::update(const int& value)
  {
    min_ = Min(min_, value);
    max_ = Max(max_, value);
  }
  //--------------------------------------------------------------
  int IntRange::min() const
  {
    return min_;
  }
  //--------------------------------------------------------------
  int IntRange::max() const
  {
    return max_;
  }
  //--------------------------------------------------------------
  // Return true if in range 
  bool IntRange::InRange(const int& value) const
  {
    if (value >= min_) {
      if (value <= max_) {
	return true;
      }}
    return false;      
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
}
