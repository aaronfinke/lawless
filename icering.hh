// icering.hh


#ifndef ICERING_HEADER
#define ICERING_HEADER

#include <vector>

#include "hkl_datatypes.hh"

namespace scala
{
  //--------------------------------------------------------------
  class IceRing
  {
  public:
    // Resolution in A, width in 1/d^2 units
    IceRing(const double& Resolution, const double& width);

    // Resolution of ring: return centre of ring as d* = 1/d
    double Dstar() const {return sqrt(ring_invressqr);}

    // If in ring, returns true
    bool InRing(const double& invresolsq) const;
    // Clear intensity sums
    void ClearSums();
    // Add in IsigI
    void AddObs(const IsigI& I_sigI, const double& invresolsq);

    void SetReject() {reject=true;}
    void SetReject(const bool& Rej) {reject=Rej;}

    // Results
    double MeanI() const;
    double MeanSigI() const;
    double MeanSSqr() const;
    int N() const {return nI;}

    // Reject flag, true means reject reflection in this range
    bool Reject() const {return reject;}

  private:
    double ring_invressqr;  // centre of ring in 1/d^2
    double halfwidth_invressqr;  // halfwidth of ring in 1/d^2
    double sum_I;
    double sum_sigI;
    double sum_sSqr;
    int nI;
    bool reject;
  };
  //--------------------------------------------------------------
  class Rings
  {
  public:
    Rings() : nrings(0) {}

    // Resolution in A, width in 1/d^2 units
    void AddRing(const double& Resolution, const double& width);

    // Set up default ice rings: 3.90, 3.67, 3.44A etc
    void DefaultIceRings();

    // Clear list
    void Clear();

    int Nrings() const {return nrings;}

    // Resolution of iring'th ring as d*
    double Dstar(const int& iring) const
    {return rings.at(iring).Dstar();}

    // Copy rejected rings only
    void CopyRejRings(const Rings& other);

    // If in ring, returns ring number (0,n-1), else = -1 
    int InRing(const double& invresolsq) const;
    // Clear intensity sums
    void ClearSums();
    // Add in IsigI
    void AddObs(const int& Iring, const IsigI& I_sigI, const double& invresolsq);

    void SetReject(const int& Iring);
    void SetReject(const int& Iring, const bool& Rej);

    // Results
    double MeanI(const int& Iring) const;
    double MeanSigI(const int& Iring) const;
    double MeanSSqr(const int& Iring) const;
    // Reject flag, true means reject reflection in this range
    bool Reject(const int& Iring) const; 
    int N(const int& Iring) const;

  private:
    int nrings;
    std::vector<IceRing> rings;

    void CheckRing(const int& Iring) const;
  };

}
#endif
