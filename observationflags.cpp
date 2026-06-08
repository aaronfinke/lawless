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
    const int ObservationFlag::FLAG_MISFIT =      64;
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
    if (TestMisfit()) {s[5] = 'X';}  // XDS Misfit
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
    acceptmisfit = false;
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
    Nmisfit = 0;
    NaccBGratio = 0;
    NaccPKratio = 0;
    NaccTooNeg = 0;
    NaccGradient = 0;
    Naccoverload = 0;
    Naccedge = 0;
    Naccmisfit = 0;

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
  void ObservationFlagControl::SetAcceptMisfit()
  {
    acceptmisfit = true;
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
      // Misfit
      if (flag.TestMisfit()) {
        Nmisfit++;
        if (acceptmisfit)
          {Naccmisfit++;}
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
    if (NBGratio+NPKratio+NTooNeg+NGradient+Noverload+Nedge+Nmisfit == 0) {
      return s;
    }
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
    s += FormatOutput::logTabPrintf(0, "   Spots on edge           %8d%8d\n",
                         Nedge, Naccedge);
    s += FormatOutput::logTabPrintf(0, "   XDS misfits (outliers)  %8d%8d\n\n",
                         Nmisfit, Naccmisfit);
    return s;
  }
  //--------------------------------------------------------------
  std::string ObservationFlagControl::XMLset(const int& nflagged, const int& naccepted,
                                             const float& maximumvalue,
                                             const float& maxaccepted) const
  {
    std::string s;
    s += StringUtil::MakeXMLtag("NumberFlagged", nflagged);
    s += StringUtil::MakeXMLtag("NumberAccepted", naccepted);
    if (maximumvalue > -9999.0) {
      s += StringUtil::MakeXMLtag("Maximum", maximumvalue, 8,3);
      s += StringUtil::MakeXMLtag("MaximumAccepted", maxaccepted,8,3);
    }
    return s+"\n";
  }
  //--------------------------------------------------------------
  std::string ObservationFlagControl::asXML() const
  {
    std::string s;
    if (NBGratio+NPKratio+NTooNeg+NGradient+Noverload+Nedge == 0) return s;

    s += "<ObservationFlags>\n";
    s += StringUtil::MakeXMLtag("BGratioTooLarge",
        XMLset(NBGratio, NaccBGratio, MaxBGratio, MaxAccBGratio), false);
    s += StringUtil::MakeXMLtag("PKratioTooLarge",
        XMLset(NPKratio, NaccPKratio, MaxPKratio, MaxAccPKratio), false);
    s += StringUtil::MakeXMLtag("TooNegative",
                XMLset(NTooNeg, NaccTooNeg), false);
    s += StringUtil::MakeXMLtag("GradientTooLarge",
        XMLset(NGradient, NaccGradient, MaxGradient, MaxAccGradient), false);
    s += StringUtil::MakeXMLtag("ProfileFittedOverloads",
                                XMLset(Noverload, Naccoverload), false);
    s += StringUtil::MakeXMLtag("Edge",
                                XMLset(Nedge, Naccedge), false);
    s += StringUtil::MakeXMLtag("Misfit",
                                XMLset(Nmisfit, Naccmisfit), false);
    s += "</ObservationFlags>\n";
    return s;
  }
//--------------------------------------------------------------
  const unsigned int ObservationStatus::wordmask;  //  = 0xFFFF
  //--------------------------------------------------------------
  void ObservationStatus::ResetStatus()
  // Clear all flags except ObsFlag & resolution flag,
  //   overlap, run & batch rejections
  {
    unsigned int mask =  OBSSTAT_FLAG | OBSSTAT_RESOLUTION |
      OBSSTAT_OVERLAP | OBSSTAT_RUN | OBSSTAT_BATCH;
    bitflags &= mask;
  }
  //--------------------------------------------------------------
  void ObservationStatus::ResetStatusAll()
  // Clear all flags except ObsFlag
  {
    unsigned int mask =  OBSSTAT_FLAG;
    bitflags &= mask;
  }
//--------------------------------------------------------------
  // true is OK or outlier or > Emax (ie suitable for Rogues file)
  // rejected or kept
  bool ObservationStatus::IsOKforRogues() const {
    if (bitflags == 0) return true;
    const unsigned int ROGUES_FLAG =
      OBSSTAT_OUTLIER & OBSSTAT_OUTLIERANOM &
      OBSSTAT_EMAX & OBSSTAT_EMAX_OK &
      OBSSTAT_STRONG & OBSSTAT_WEAK & OBSSTAT_DEVIANT;
    return (bitflags & ROGUES_FLAG) == 0;
  }
  //--------------------------------------------------------------
    //  accepted if no status bits are set (except DEVIANT & EMAX_OK)
  bool ObservationStatus::IsAccepted() const
  {
    if (bitflags == 0) {
      return true;
    }
    const unsigned int ACCEPT_MASK =
      wordmask - (OBSSTAT_DEVIANT | OBSSTAT_EMAX_OK);
    return ((bitflags & ACCEPT_MASK) == 0) ;
  }
  //--------------------------------------------------------------
  // add in status
  ObservationStatus ObservationStatus::mergestatus(const ObservationStatus& status) const
  {
    unsigned int bitstatus = (bitflags | status.Bitflags());
    return bitstatus;
  }
  //--------------------------------------------------------------
  // Format for debugging
  std::string ObservationStatus::format() const
  {
    std::string s;
    if (IsAccepted()) {
      s = "Accepted";
      return s;
    }
    if (TestObsFlag()) {s += "|OBSSTAT_FLAG";}
    if (TestResolution()) {s += "|OBSSTAT_RESOLUTION";}
    if (TestOutlier()) {s += "|OBSSTAT_OUTLIER";}
    if (TestOutlierAnom()) {s += "|OBSSTAT_OUTLIERANOM";}
    if (TestEmax()) {s += "|OBSSTAT_EMAX";}
    if (TestEmaxOK()) {s += "|OBSSTAT_EMAX_OK";}
    if (TestTooStrong()) {s += "|OBSSTAT_STRONG";}
    if (TestTooWeak()) {s += "|OBSSTAT_WEAK";}
    if (TestRejectOverlap()) {s += "|OBSSTAT_OVERLAP";}
    if (TestRejectRun()) {s += "|OBSSTAT_RUN";}
    if (TestRejectBatch()) {s += "|OBSSTAT_BATCH";}
    if (TestDeviant()) {s += "|OBSSTAT_DEVIANT";}
    return s;
  }
//--------------------------------------------------------------
} // namespace scala
