// observationflags.cpp

#include <iostream>
#include <stdio.h>

#include "observationflags.hh"
#include "util.hh"
#include "string_util.hh"

namespace scala
{
//--------------------------------------------------------------
    const int ObservationFlag::FLAG_BGRATIO =      1;
    const int ObservationFlag::FLAG_PKRATIO =      2;
    const int ObservationFlag::FLAG_TOONEGATIVE =  4;
    const int ObservationFlag::FLAG_GRADIENT =     8;
    const int ObservationFlag::FLAG_OVERLOAD =    16;
    const int ObservationFlag::FLAG_EDGE =        32;
//--------------------------------------------------------------
  ObservationFlag::ObservationFlag(const ObservationFlag& flag)
  {
    bitflags = flag.bitflags;
    bgpk = flag.bgpk;
  }
//--------------------------------------------------------------
  // Return unpacked values of BGratio, PKratio, gradient
  void ObservationFlag::BgPkValues
  (float& BGratio, float& PKratio, float& Gradient) const
  {
    int igood = int(bgpk);
    Gradient = bgpk - float(igood);
    int ipkr = igood/100;
    BGratio = 0.1*(igood - ipkr*100);
    PKratio = 0.01*ipkr;
  }
//--------------------------------------------------------------
  float PackFlagValues(const float& BGratio,
		       const float& PKratio, const float& Gradient)
  {
    float pv = Min(Gradient, 0.99);
    int iv = Min(Nint(10.*BGratio), 99) + Nint(100.*PKratio)*100;
    return pv + float(iv);
  }		 
//--------------------------------------------------------------
  // Combine with another flags
  void ObservationFlag::AddFlag(const ObservationFlag& other)
  {
    // set bits from either source
    bitflags = bitflags | other.bitflags;
    // Get maximum values
    float bgr, pkr, grd;
    BgPkValues(bgr, pkr, grd);
    float bgr2, pkr2, grd2;
    other.BgPkValues(bgr2, pkr2, grd2);
    bgr = Max(bgr, bgr2);
    pkr = Max(pkr, pkr2);
    grd = Max(grd, grd2);
    bgpk = PackFlagValues(bgr, pkr, grd);
  }
//--------------------------------------------------------------
  //! return brief formatted version of which flags are set
  std::string ObservationFlag::format() const
  {
    std::string s = "      ";
    if (TestBGratio()) {s[0] = 'B';}
    if (TestPKratio()) {s[1] = 'P';}
    if (TestTooNeg()) {s[2] = 'N';}
    if (TestGradient()) {s[3] = 'G';}
    if (TestOverload()) {s[4] = 'O';}
    if (TestEdge()) {s[5] = 'E';}
    return s;
  }
//--------------------------------------------------------------
//--------------------------------------------------------------
  ObservationFlagControl::ObservationFlagControl()
  {Init();}
//--------------------------------------------------------------
  ObservationFlagControl::ObservationFlagControl
  (const float& BGrlim,  const float& PKrlim,  const float& Gradlim,
   const bool& AcceptOverload, const bool& AcceptEdge)
    :  bgrlim(BGrlim), pkrlim(PKrlim), grdlim(Gradlim),
       acceptoverload(AcceptOverload), acceptedge(AcceptEdge)
  {
    // Clear counts
    Clear();
  }
//--------------------------------------------------------------
  void ObservationFlagControl::Init()
  {
    bgrlim = -1.0;
    pkrlim = -1.0;
    grdlim = -1.0;
    acceptoverload = false;
    acceptedge = false;
    Clear();
    
  }
//--------------------------------------------------------------
  void ObservationFlagControl::Clear()
  {
    NBGratio = 0;
    NPKratio = 0;
    NTooNeg = 0;
    NGradient = 0;
    Noverload = 0;
    Nedge = 0;
    NaccBGratio = 0;
    NaccPKratio = 0;
    NaccTooNeg = 0;
    NaccGradient = 0;
    Naccoverload = 0;
    Naccedge = 0;

    MaxBGratio = 0.0;
    MaxPKratio = 0.0;
    MaxGradient = 0.0;
    
    MaxAccBGratio = 0.0;
    MaxAccPKratio = 0.0;
    MaxAccGradient = 0.0;
  }
//--------------------------------------------------------------
  void ObservationFlagControl::SetBGRlimit(const float& BGrlim)
  {
    bgrlim = BGrlim;
  }
//--------------------------------------------------------------
  void ObservationFlagControl::SetPKRlimit(const float& PKrlim)
  {
    pkrlim = PKrlim;
  }
//--------------------------------------------------------------
  void ObservationFlagControl::SetGradlimit(const float& Gradlim)
  {
    grdlim = Gradlim;
  }
//--------------------------------------------------------------
  void ObservationFlagControl::SetAcceptOverload()
  {
    acceptoverload = true;
  }
//--------------------------------------------------------------
  void ObservationFlagControl::SetAcceptEdge()
  {
    acceptedge = true;
  }
//--------------------------------------------------------------
  // Returns true is observation accepted, & count them
  bool ObservationFlagControl::IsAccepted(const ObservationFlag& flag)
  {
    if (flag.OK()) {
      // Acceptance due to no bit flags set, no counting
      return true;
    } else {
      // Some bit flags are set
      bool OK = true;
      float bgr, pkr, grd;
      flag.BgPkValues(bgr, pkr, grd);
      MaxBGratio = Max(bgr, MaxBGratio);
      MaxPKratio = Max(pkr, MaxPKratio);
      MaxGradient = Max(grd, MaxGradient);
      
      // BGratio
      if (flag.TestBGratio()) {
	NBGratio++;
	if (bgrlim > 0.0 && bgr < bgrlim) {
	  NaccBGratio++;
	  MaxAccBGratio = Max(bgr, MaxAccBGratio);
	} else
	  {OK = false;}
      }
      if (flag.TestPKratio()) {
	NPKratio++;
	if (pkrlim > 0.0 && pkr < pkrlim) {
	  NaccPKratio++;
	  MaxAccPKratio = Max(pkr, MaxAccPKratio);
	} else
	  {OK = false;}
      }
      // Gradient
      if (flag.TestGradient()) {
	NGradient++;
	if (grdlim > 0.0 && grd < grdlim) {
	  NaccGradient++;
	  MaxAccGradient = Max(grd, MaxAccGradient);
	} else
	  {OK = false;}
      }
      // TooNeg
      if (flag.TestTooNeg()) {
	NTooNeg++;
	OK = false;
      }
      // Overload
      if (flag.TestOverload()) {
	Noverload++;
	if (acceptoverload) 
	  {Naccoverload++;}
	else
	  {OK = false;}
      }
      // Edge
      if (flag.TestEdge()) {
	Nedge++;
	if (acceptedge) 
	  {Naccedge++;}
	else
	  {OK = false;}
      }
      return OK;
    }
  }
  //--------------------------------------------------------------
  std::string ObservationFlagControl::PrintCounts() const
  {
    std::string s;
    if (NBGratio+NPKratio+NTooNeg+NGradient+Noverload+Nedge == 0) return s;

    s += FormatOutput::logTab(0,
		   "\n\nNumbers of observations marked in the FLAG column");
    s += FormatOutput::logTab(0,
		   "By default all flagged observations are rejected");
    s += FormatOutput::logTab(0,
		   "Observations may be counted in more than one category\n\n");
    s += FormatOutput::logTab(0,
		   "                             Flagged  Accepted   Maximum   MaxAccepted");
    s += FormatOutput::logTabPrintf(0, "   BGratio too large       %8d%8d%12.3f%12.3f\n",
			 NBGratio, NaccBGratio, MaxBGratio, MaxAccBGratio);
    s += FormatOutput::logTabPrintf(0, "   PKratio too large       %8d%8d%12.3f%12.3f\n",
			 NPKratio, NaccPKratio, MaxPKratio, MaxAccPKratio);
    s += FormatOutput::logTabPrintf(0, "   Negative < 5sigma       %8d%8d\n",
			 NTooNeg, NaccTooNeg);
    s += FormatOutput::logTabPrintf(0, "   Gradient too large      %8d%8d%12.3f%12.3f\n",
			 NGradient, NaccGradient, MaxGradient, MaxAccGradient);
    s += FormatOutput::logTabPrintf(0, "   Profile-fitted overloads%8d%8d\n",
			 Noverload, Naccoverload);
    s += FormatOutput::logTabPrintf(0, "   Spots on edge           %8d%8d\n\n",
			 Nedge, Naccedge);
    return s;
  }
//--------------------------------------------------------------
  const unsigned int ObservationStatus::wordmask;  //  = 0xFFFF
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
}
