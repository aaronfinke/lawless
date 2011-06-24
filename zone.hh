// zone.hh

#ifndef ZONE_HEADER
#define ZONE_HEADER

#include <clipper/clipper.h>
#include "hkl_datatypes.hh"
#include "range.hh"
#include "globalcontrols.hh"
#include "hkl_symmetry.hh"
#include "score_datatypes.hh"

namespace scala
{
  //--------------------------------------------------------------
  double IscoreVal(const IsigI& Isig);
  // Score function value from I,sig
  //  = I/sigI 
  //=================================================================
  class IndexIsigI
  {
  public:
    IndexIsigI(){}
    IndexIsigI(const int& j, const IsigI& Is) : index(j), Isig(Is) {}

    int index;
    IsigI Isig;
  };
  //=================================================================
  class OneDFourier
  // One-dimensional Fourier, just a few grid points (eg 0, 1/2)
  {
  public:
    OneDFourier(){}
    OneDFourier(const std::vector<int>& Ngrid);
    OneDFourier(const std::vector<double>& Xgrid);

    // Add into totals
    void AddRef(const int& j, const double& Val);
    // Return number of contributions
    int Nobs() const;
    std::vector<double> FourierVal() const;


  private:
    std::vector<int> ngrid;   // grid intervals are 1/ngrid
    std::vector<double> xgrid;
    std::vector<double> Sumx;
    int nobs;
  };
  //=================================================================
  class Zone
  {
  public:
    Zone() {}  // dummy constructor
    // Construct for screw axis
    //  Axis == a,b,c
    //  Nfold is order of axis
    // Nfold negated to test only the 1/Nfold point
    // (sets validpoint array)
    //   this allows testing of 4(1) but not 4(2)
    // LGsym is Laue group symmetry
    Zone(const std::string& Axis, const int& Nfold,
	 const hkl_symmetry& LGsym); 
    // Construct for glide plane
    //  GlidePlane == a,b,c,110,101,101  (normal to mirror)
    //  Glide = a,b,c,n,d  glide direction
    // LGsym is Laue group symmetry
    Zone(const std::string& GlidePlane, const std::string& Glide,
	 const hkl_symmetry& LGsym);

    // return true if hkl is in defined zone
    bool InZone(const Hkl& hkl) const;
    // Return zone type, true if axis, false if glide
    bool Axis() const {return axis;}
    // Return angular order
    int Order() const {return order;}
    // Format zone, in "new" frame
    // RefToNew is reindex operator from "reference" (lattice) frame
    // (stored with StoreReindex function) to some new frame in
    // which to format the zone
    // if Sub > 1, the screw component is order/Sub
    std::string formatNewFrame(const ReindexOp& RefToNew, const int Sub=1) const;
    std::string formatRefFrame(const int Sub=1) const;
    // Format in Laue group frame
    std::string formatLGFrame(const int Sub=1) const;
    // Format condition, eg "0kl: h+k = 2n"
    // For axis, Sub is grid point
    std::string FormatConditionNewFrame(const ReindexOp& RefToNew,
				const int Sub=1) const;
    std::string FormatConditionRefFrame(const int Sub=1) const;
    std::string FormatConditionLGFrame(const int Sub=1) const;

    // Axis direction in reference frame
    std::string Direction() const;

    // Return true if two zones apply to the same hkl zone
    // Glides only, eg b(a), c(a), n(a) all apply to 0kl
    bool SameZone(const Zone& other) const;

    // Add into totals
    void AddRef(const Hkl& hkl, const IsigI& Isig, const float& sSqr);

    // Store SD of control reflections for Z-score
    // Adjust to avoid zeroes, but set flag to indicate this
    void StoreMeanSD(const std::vector<double>& Mean, const std::vector<double>& SD);

    // Store reindex from some standard reference frame, for comparison of zones
    void StoreReindex(const ReindexOp& reindexin);
    ReindexOp Reindex() const {return reindexmat;}
    // Store overall probability if valid, ie except for 4- & 6-fold screws
    //  fail if invalid
    void StoreProb(const double& ProbYes);

    // Return Laue group symmetry
    hkl_symmetry LGsymm() const {return lgsymm;}


    // Return Fourier values at each grid point
    std::vector<double> FourierVal() const;

    // Return probability at each grid point
    std::vector<double> p() const;
    // Return overall probability if valid, ie except for 4- & 6-fold screws
    //  fail if invalid
    double Prob() const;
    // Return grid points
    std::vector<int> Ngrid() const {return ngrid;}
    // Return number of grid points
    int NgridPoints() const {return npoint;}
    // Return number of contributions
    int Nobs() const;
    // Inverse resolution range
    Range InvResoRange() const {return InvResRange;}
    // Return minimum & maximum indices
    std::pair<int,int> IndexRange() const 
    {return std::pair<int,int>(minIndx, maxIndx);}
    // Return list of indices
    std::vector<int> Indices() const;
    // Return valid flag
    bool Valid() const {return (valid);}
    //    bool Valid() const {return (valid && active);}
    bool ValidPoint(const int& ig) const {return validpoint.at(ig);}
    // Return SD of control reflections for Z-score
    double GetSD(const int& i) const {return controlsd[i];}
    // Return true if data are systematically missing eg all odd indices
    bool PrunedData() const {
      // Return true if any are true
      for (int i=0;i<npoint;++i) {
	if (prunedData[i]) {return true;}
      }
      return false;
    }
    // Return true if data are systematically missing eg all odd indices
    bool PrunedData(const int& i) const {return prunedData[i];}

    // Return two test hkl which will be systematically absent
    // if in this zone, for glide planes
    // RefToNew is reindex operator from "reference" (lattice) frame
    // (stored with StoreReindex function) to some new frame in
    // which to return the indices
    std::vector<Hkl> GlideTestHkl(const ReindexOp& RefToNew) const;
    // Return test hkl's which will be systematically absent
    // if in this zone, for axes pricipal only, not 110)
    // One test hkl for each grid point, test in reverse order
    std::vector<Hkl> AxisTestHkl(const ReindexOp& RefToNew,
		    std::vector<std::vector<int> >& screwabsence) const;

    std::vector<IndexIsigI> IndexedData() const {return IndxIsigI;}

   // weighted average I, sigma
    std::vector<IsigI> AverageIsigI() const { return averageIsigI;}
    // average I,sig after neighbour "correction"
    std::vector<IsigI> AdjustedIsigI() const {return adjustedIsigI;}


    // Fraction of intensity on screw axis to subtract from neighbouring reflection
    static void SetNeighbourFraction(const float& nf) {neighbourFraction = nf;}
    static float GetNeighbourFraction() {return neighbourFraction;}

    void dump() const;

    // For sorting, all glides after axes, then on order
    friend bool operator < (const Zone& a,const Zone& b)
      // true if a before b in sort order
    {
      if (a.axis) {
	if (!b.axis) return true;  // b is glide, a is axis
      } else {
	if (b.axis) return false;    // a is glide, b is axis
      }
      // Same, compare "order"
      return a.order > b.order;
    }

    // Compare two zones, return true if equivalent after allowing for reindexing
    bool Compare(const Zone& other) const;

  private:
    void init(const bool& disable);
    void TestValidIndices() const;
    void CalcResults() const;
    Hkl CondReference() const;
    std::vector<double> GetFourierValues() const;
    std::vector<IsigI> AverageI(const std::vector<IndexIsigI>& IdxIs) const;

    clipper::Vec3<double> DirectionRefFrame() const;
    // Get zone direction in reference frame
    //  vector returned is axis direction or zone normal (for glide plane)



    // Permutation operator to convert input hkl to internal
    // standard  h'T = hT [P]
    // Axis: 00l
    // Principle glide: hk0
    // Diagonal glide: hhl
    ReindexOp permute;  
    // Reindex from lattice frame to internal standard [P'] = [H] [P]
    ReindexOp permuteIndex;

    // Index along axis = [hc].qc
    //  where hc is hkl in constructor frame &
    //  cond  = qc
    //  cond is stored here permuted by [P] permutation
    //       to convert external hkl
    clipper::Vec3<int> cond;
    //  qcond = ql = [H] qc    where [H] is the reindexing operator from
    //                         external reference (lattice) frame
    //                         Note that cond = qc is integral, while qcond = ql
    //                         may be half-integral if the lattice is centred
    clipper::Vec3<double> qcond;

    // Direction is axis direction or glide plane
    //  "a", "b" "c", "110", "101", or "011"
    std::string direction;
    std::string glide;
    int order;
    bool axis;   // true if axis, false if glideplane
    bool diagonal;  // true if glide zone is diagonal ie 110
    Range InvResRange;
    std::vector<IndexIsigI> IndxIsigI;
    int minIndx;
    int maxIndx;
    // next 2 vectors indexed by axial index 0 -> maxIndx
    mutable std::vector<IsigI> averageIsigI;   // weighted average I, sigma
    mutable std::vector<IsigI> adjustedIsigI;  // average I,sig after neighbour "correction"

    mutable std::vector<bool> validpoint;  // Permanent flag for each point
                                   // this will set to false to avoid testing
                                   // this grid point, eg so as not to test
                                   // 4(2) axes in an I or F lattice

    // Totals for Fourier analysis
    std::vector<int> ngrid;   // grid intervals are 1/ngrid
                              //  eg for 6-fold axis, ngrid = (1,2,3,6)
    mutable OneDFourier fsum;                 // 1D Fourier at grid intervals
    std::vector<double> controlsd;   // Standard deviation of control score for each grid point
    std::vector<double> controlmean; //  ... and its mean
    int npoint;   // ngrid.size()
    bool ZeroControlSD;             // at least one control SD was zero before resetting
    std::vector<bool> UnitControlMean;           // at least one control Mean == 1.0 and control SD == 0.0
    static double MinControlSD;

    // fraction of neighbouring reflection to subtract from intensity
    static float neighbourFraction;

    // Reindex operator from reference cell to constructor frame
    // for comparison of zones
    ReindexOp reindexmat;
    bool NonIdentityReindex;   // true if reindexmat is not identity
    hkl_symmetry lgsymm;       // Laue group symmetry


    bool singleprob;           // true if only the last probability (Fourier)
			       // point is useful. This is the case for all glides
			       // and 2- & 3-fold screws
    // Things set in method CalcResults
    mutable bool valid;    // true if we have enough data
    mutable bool results;
    mutable std::vector<double> Pfor;
    mutable double prob_yes;   // unique probability values if singleprob = true
    mutable std::vector<bool> prunedData;  // true if indices of reflection data are
                                         // systematically missing (eg only even orders)
                                         // so that the zone is indeterminate

  };   // Zone
  //================================================================
  class SysAbsScore
  //  Class to accumulate & calculate total "score" for systematic absences
  //  in a putative spacegroup
  //  Each condition (eg A) provide a probability p(A) which is
  //  multiplied in either as p(A) or p(!A) = 1-p(A) depending on the
  //  putative spacegroup
  //  
  {
  public:
    SysAbsScore() : possible(true), TotalProb(1.0), n(0) {}
    SysAbsScore(const PossibleSpaceGroup& PSG)
      :  possible(true), TotalProb(1.0), n(0), PossibleSG(PSG) {}

    // Add in contribution (only if non-zero)
    void ProbYes(const double& p) {if (p >= 0.0) {TotalProb *= p; n++;}}
    // Add in 1-contribution 
    void ProbNo(const double& p) {if (p >= 0.0) {TotalProb *= 1.0-p; n++;}}
    // Set as impossible
    void Impossible() {possible=false;} 
    // Set as possible
    void Possible() {possible=true;}
    // Set condition string
    void SetCondition(const std::string& cnd) {condition = cnd;}
    // Set condition string
    void SetConditionLG(const std::string& cnd) {conditionLG = cnd;}
    // Possible spacegroup
    void SetSGposs(const PossibleSpaceGroup& PSG) {PossibleSG = PSG;}
    // Normalise by multiplying factor
    void Normalise(const double& fac) {TotalProb *= fac;}
    // Add zone number to list
    void AddZone(const int& izone) {zonesingroup.push_back(izone);}
    // Return zone list
    std::vector<int> ZoneList() const {return zonesingroup;}
    // Return true if zone iz is in group
    bool IsZoneInGroup(const int& iz) const;


    // return total probability
    double TotalProbability() const {return TotalProb;}
    // Possible flag
    bool IsPossible() const {return possible;}
    // Condition string
    std::string Condition() const {return condition;}
    // Condition string Laue group frame
    std::string ConditionLG() const {return conditionLG;}

    PossibleSpaceGroup SGposs() const {return PossibleSG;}

  private:
    bool possible;
    double TotalProb;
    int n;
    std::string condition;
    std::string conditionLG;   // Laue group frame
    PossibleSpaceGroup PossibleSG;
    std::vector<int> zonesingroup;   // list of zone numbers in this group
  };  // SysAbsScore

}

#endif
