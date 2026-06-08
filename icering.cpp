// icering.cpp

#include <assert.h>
#define ASSERT assert

#include "icering.hh"
#include "string_util.hh"

using phaser_io::LOGFILE;
using phaser_io::LXML;

namespace scala
{
  //--------------------------------------------------------------
  void Rings::setIceRings(const int& ringtype)
  {
    nrings = 0;
    rings.clear();
    listtype = ringtype;
    if (ringtype == 1) {
      // Set up default ice rings: 3.90, 3.67, 3.44A
      // clear any existing ones first
      //
      // Harry's list - note the widths are variable
      //   3.94 - 3.86
      //   3.705 - 3.635
      //   3.465 - 3.415
      //   2.690 - 2.650
      //   2.265 - 2.235
      //   2.090 - 2.070
      //   1.960 - 1.940
      //   1.925 - 1.915
      //   1.895 - 1.885
      //   1.730 - 1.720
      
      // use a constant width in reciprocal space
      // rings should get wider at higher resolution, but they
      // probably get weaker as well
      const double RWIDTH = 0.005;
      // resolution in A, full width in d* 1/A
      AddRing(3.8996, RWIDTH);
      AddRing(3.6697, RWIDTH);
      AddRing(3.4398, RWIDTH);
      AddRing(2.6699, RWIDTH);
      AddRing(2.2499, RWIDTH);
      AddRing(2.0800, RWIDTH);
      AddRing(1.9499, RWIDTH);
      AddRing(1.9200, RWIDTH);
      AddRing(1.8900, RWIDTH);
      AddRing(1.7250, RWIDTH);
    } else if (ringtype == 2) {
      // This is a list from Clemens Vonrhein, as used in autoPROC

      AddRing(3.895194,0.00081690);
      AddRing(3.663871,0.00097591);
      AddRing(3.439251,0.00103830);
      AddRing(2.667950,0.00152917);
      AddRing(2.248670,0.00190315);
      AddRing(2.068321,0.00232813);
      AddRing(1.947494,0.00245315);
      AddRing(1.916274,0.00256933);
      AddRing(1.882113,0.00275378);
      AddRing(1.719354,0.00336036);
      AddRing(1.522158,0.00399225);
      AddRing(1.471920,0.00395105);
      AddRing(1.443011,0.00392707);
      AddRing(1.370488,0.00512793);
      AddRing(1.298091,0.00523005);
      AddRing(1.260481,0.00547263);
      AddRing(1.221687,0.01447748);
      AddRing(1.168801,0.01540575);
      AddRing(1.123103,0.01259269);
      AddRing(1.072467,0.03340847);
      AddRing(1.037006,0.01969054);
      AddRing(1.009667,0.01538879);
      AddRing(0.987219,0.01554805);
      AddRing(0.973464,0.01795071);
      AddRing(0.964942,0.01729539);
      AddRing(0.940017,0.02181351);
      AddRing(0.921414,0.02055424);
      AddRing(0.915060,0.01816664);
      AddRing(0.904326,0.01841565);
      AddRing(0.890075,0.04032559);
      AddRing(0.869309,0.02100578);
      AddRing(0.849285,0.03303786);
      AddRing(0.839363,0.02110332);
      AddRing(0.828018,0.02473026);
      AddRing(0.810936,0.02677606);
      AddRing(0.797286,0.02178282);
      AddRing(0.778083,0.03175621);
      AddRing(0.763379,0.02602994);
      AddRing(0.749456,0.03634337);
      AddRing(0.742847,0.02517499);
      AddRing(0.735067,0.04878942);
      AddRing(0.721629,0.03395375);
      AddRing(0.710740,0.02991279);
      AddRing(0.705646,0.02937006);
      AddRing(0.698276,0.03911991);
      AddRing(0.688685,0.03232857);
      AddRing(0.680085,0.03318906);
      AddRing(0.673595,0.03415765);
      AddRing(0.667806,0.03397932);
      AddRing(0.658283,0.04865441);
      AddRing(0.650295,0.03708961);
      AddRing(0.639763,0.05178115);
      AddRing(0.632089,0.04751414);
      AddRing(0.624413,0.04386679);
      AddRing(0.620587,0.03849147);
      AddRing(0.615737,0.05208592);
      AddRing(0.607883,0.04665315);
      AddRing(0.588122,0.05190859);
    }
  }
  //--------------------------------------------------------------
  void Rings::CheckRing(const int& Iring) const
  {
    ASSERT (Iring < nrings && Iring >= 0);
  }
  //--------------------------------------------------------------
  // Resolution in A, width in 1/d^2 units
  void Rings::AddRing(const double& Resolution, const double& width)
  {
    rings.push_back(IceRing(Resolution, width));
    nrings++;
  }
  //--------------------------------------------------------------
  // Copy rejected rings only
  void Rings::CopyRejRings(const Rings& other)
  {
    rings.clear();
    for (int i=0;i<other.nrings;i++) {
      if (other.rings[i].Reject())
        {rings.push_back(other.rings[i]);}
    }
    nrings = rings.size();
    listtype = other.listtype;
  }
  //--------------------------------------------------------------
  // Clear list
  void Rings::Clear()
  {
    nrings = 0;
    rings.clear();
  }
  //--------------------------------------------------------------
  // If in ring, returns ring number (0,n-1), else = -1
  int Rings::InRing(const double& invresolsq) const
  {
    for (size_t i=0;i<rings.size();i++) {
      if (rings[i].InRing(invresolsq))
        {return i;}
    }
    return -1;
  }
  //--------------------------------------------------------------
  void Rings::ClearSums()
  {
    for (size_t i=0;i<rings.size();i++) {
      rings[i].ClearSums();
    }
  }
  //--------------------------------------------------------------
  void Rings::ClearCounts() const
  {
    for (size_t i=0;i<rings.size();i++) {
      rings[i].ClearCount();
    }
  }
  //--------------------------------------------------------------
  void Rings::AddObs(const int& Iring, const IsigI& I_sigI,
                     const double& invresolsq)
  {
    CheckRing(Iring);
    rings[Iring].AddObs(I_sigI, invresolsq);
  }
 //--------------------------------------------------------------
  void Rings::SetReject(const int& Iring)
  {
    CheckRing(Iring);
    rings[Iring].SetReject();
  }
 //--------------------------------------------------------------
  void Rings::SetReject(const int& Iring, const bool& Rej)
  {
    CheckRing(Iring);
    rings[Iring].SetReject(Rej);
  }
  //--------------------------------------------------------------
  void Rings::SetRejectAll(const bool& Rej)
  {
    for (size_t i=0;i<rings.size();i++) {
      rings[i].SetReject(Rej);
    }
  }
 //--------------------------------------------------------------
  bool Rings::Reject(const int& Iring) const
  {
    CheckRing(Iring);
    return rings[Iring].Reject();
  }
 //--------------------------------------------------------------
  double Rings::MeanI(const int& Iring) const
  {
    CheckRing(Iring);
    return rings[Iring].MeanI();
  }
  //--------------------------------------------------------------
  double Rings::MeanSigI(const int& Iring) const
  {
    CheckRing(Iring);
    return rings[Iring].MeanSigI();
  }
  //--------------------------------------------------------------
  double Rings::MeanSSqr(const int& Iring) const
  {
    CheckRing(Iring);
    return rings[Iring].MeanSSqr();
  }
  //--------------------------------------------------------------
  int Rings::N(const int& Iring) const
  {
    CheckRing(Iring);
    return rings[Iring].N();
  }
  //--------------------------------------------------------------
  int Rings::NinRing(const int& Iring) const
  {
    CheckRing(Iring);
    return rings[Iring].NinRing();
  }
  //--------------------------------------------------------------
  void Rings::report(phaser_io::Output& output) const
  {
    output.logTab(0,LXML,"<IceRings>");
    std::string message;
    if (nrings == 0) {
      if (listtype > 0) {
	message = "Ice rings were rejected from scaling, but kept in statistics and output";
      } else {
	message = "No ice ring rejections from scaling, statistics or output file";
      }
    } else {
	message = "Ice rings were rejected from scaling, statistics and output";
    }
    output.logTab(0,LOGFILE, message);
    output.logTab(0,LXML,
		    StringUtil::MakeXMLtag("Message", message));
    output.logTab(0,LXML,
		  StringUtil::MakeXMLtag("RingType", listtype));
    if (nrings > 0) {
      int total = 0;
      for (size_t i=0; i<rings.size(); i++) { 
	total += rings[i].NinRing();
      }
      output.logTabPrintf(0,LOGFILE,
	   "Number of reflections rejected in ice rings = %6d\n",
			  total);
      output.logTab(0,LXML,
		    StringUtil::MakeXMLtag("TotalRejected", total));
      output.logTab(0,LOGFILE, "For each ring:");
      for (size_t i=0; i<rings.size(); i++) { 
	if (rings[i].NinRing() > 0) {
	  output.logTab(0,LXML,"<Ring>");
	  double resolution = 1.0/rings[i].Dstar();
	  output.logTabPrintf(1,LOGFILE,
			      "Resolution %7.3f Number rejected %6d\n",
			      resolution, rings[i].NinRing());
	  output.logTab(0,LXML,
			StringUtil::MakeXMLtag("RingIndex", int(i)));
	  output.logTab(0,LXML,
			StringUtil::MakeXMLtag("RingResolution", resolution));
	  output.logTab(0,LXML,
			StringUtil::MakeXMLtag("RingWidth", rings[i].width(),
					       12,8));
	  output.logTab(0,LXML,
			StringUtil::MakeXMLtag("NumberRejected", rings[i].NinRing()));
	  output.logTab(0,LXML,"</Ring>");
	}
      }
    }
    output.logTab(0,LXML,"</IceRings>");
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  IceRing::IceRing(const double& Resolution, const double& width)
  {
    ring_invressqr = 1./(Resolution*Resolution);
    halfwidth_invressqr = 0.5*width;
    ClearSums();
    reject = false;
  }
  //--------------------------------------------------------------
  bool IceRing::InRing(const double& invresolsq) const
  {
    if (Close<double,double>(invresolsq,
                             ring_invressqr,halfwidth_invressqr))  {
      ninring++;
      return true;
    }
    return false;
  }
  //--------------------------------------------------------------
  void IceRing::ClearSums()
  {
    sum_I = 0.0;
    sum_sigI = 0.0;
    sum_sSqr = 0.0;
    nI = 0;
    ninring = 0;
  }
  //--------------------------------------------------------------
  void IceRing::ClearCount() const
  {
    ninring = 0;  // mutable
  }
  //--------------------------------------------------------------
  void IceRing::AddObs(const IsigI& I_sigI, const double& invresolsq)
  {
    sum_I += I_sigI.I();
    sum_sigI += I_sigI.sigI();
    sum_sSqr += invresolsq;
    nI++;
  }
  //--------------------------------------------------------------
  double IceRing::MeanI() const
  {
    if (nI > 0)
      {return sum_I/nI;}
    else
      {return 0.0;}
  }
  //--------------------------------------------------------------
  double IceRing::MeanSigI() const
  {
    if (nI > 0)
      {return sum_sigI/nI;}
    else
      {return 0.0;}
  }
  //--------------------------------------------------------------
  double IceRing::MeanSSqr() const
  {
    if (nI > 0)
      {return sum_sSqr/nI;}
    else
      {return 0.0;}
  }
  //--------------------------------------------------------------
}
