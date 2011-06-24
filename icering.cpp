// icering.cpp

#include <assert.h>
#define ASSERT assert

#include "icering.hh"

namespace scala
{
  //--------------------------------------------------------------
  void Rings::DefaultIceRings()
  {
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

    nrings = 0;
    rings.clear();
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
  //--------------------------------------------------------------
  IceRing::IceRing(const double& Resolution, const double& width)
  {
    ring_invressqr = 1./(Resolution*Resolution);
    halfwidth_invressqr = 0.5*width;
    ClearSums();
  }
  //--------------------------------------------------------------
  bool IceRing::InRing(const double& invresolsq) const
  {
    if (Close<double,double>(invresolsq,
			     ring_invressqr,halfwidth_invressqr))
      {return true;}
    return false;
  }
  //--------------------------------------------------------------
  void IceRing::ClearSums()
  {
    sum_I = 0.0;
    sum_sigI = 0.0;
    sum_sSqr = 0.0;
    nI = 0;
    reject = false;
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
