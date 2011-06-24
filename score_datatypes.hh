// score_datatypes.hh
//
// Data type definitions



#ifndef SCORE_DATATYPES_HEADER
#define SCORE_DATATYPES_HEADER

// Clipper
#include <clipper/clipper.h>
using clipper::Vec3;
using clipper::Mat33;
using clipper::String;
using clipper::Metric_tensor;
typedef clipper::Vec3<double> DVect3;
typedef clipper::Mat33<double> DMat33;
typedef clipper::Vec3<float> FVect3;
typedef clipper::Mat33<float> FMat33;

#include "cmtzlib.h"    // CCP4 MTZlib headers (namespace CMtz)
#include "csymlib.h"    // CCP4 symmetry stuff
#include "matvec_utils.hh"  // Matrix & vector utilities
#include "Output.hh"
#include "hkl_datatypes.hh"
#include "globalcontrols.hh"

#include "util.hh"

typedef float  Rtype;
typedef double Dtype;

typedef std::pair<int,int> IndexPair;
typedef std::pair<float,float> RPair;
typedef std::pair<double,double> DPair;



namespace scala
{
  //--------------------------------------------------------------
  class LinearFit
  // Very simple-minded linear regression
  {
  public:
    LinearFit()
      : sumw(0.0),sumwx(0.0),sumwy(0.0),sumwxx(0.0),sumwxy(0.0),np(0)
    {maxy=0.0;miny=0.0;}

    // Add in one contribution
    void add(const float& x, const float& y, const float& w);
    // Return slope, intercept
    RPair result() const;
    // Return slope only, fixed intercept b (eg = 0)
    double slope(const float& b=0.0) const;
    // Return number of points
    int Number() const {return np;}

  private:
    double sumw, sumwx, sumwy, sumwxx, sumwxy;
    float  miny,maxy;
    int np;
  };

  //--------------------------------------------------------------
  class IsigNorm {
    // I/sigI & normalisation factor (multiplying)
  public:
    IsigNorm(){}
    IsigNorm(const IsigI& Is, const float& f)
      : Isig(Is), fac(f) {}

    IsigI Isig;
    float fac;  // dividing correction
  };
  //--------------------------------------------------------------
  //  typedef std::pair<IsigI, IsigI> IsPair;
  class IsPair {
  public:
    // a pair of I/sigI with run & "time" (batch) information

    IsPair(){}
    IsPair(const IsigI& Is1_in, const float& f1_in,
	   const IsigI& Is2_in, const float& f2_in)
      : IsN1(Is1_in, f1_in), IsN2(Is2_in, f2_in) {}

    void Set(const IsigI& Is1_in, const float& f1_in,
	     const IsigI& Is2_in, const float& f2_in)
    {
      IsN1.Isig = Is1_in;
      IsN1.fac = f1_in;
      IsN2.Isig = Is2_in;
      IsN2.fac = f2_in;
    }

    IsigI Is1() const {return IsN1.Isig;}
    //    float f1() const {return IsN1.fac;}

    IsigI Is2() const {return IsN2.Isig;}
    //    float f2() const {return IsN2.fac;}
    // Normalised
    IsigI Is1norm() const {return IsN1.Isig.scaleIs(IsN1.fac);}
    IsigI Is2norm() const {return IsN2.Isig.scaleIs(IsN2.fac);}

  private:
    IsigNorm IsN1, IsN2;
  };

  //--------------------------------------------------------------
  class IKode
  {
  public:
    IKode() {}
    IKode(const IsigI& Is_in, const int& kode_in, const float& invres,
	  const int& irun_in, const float& time_in)
      : Is(Is_in), sSqr(invres), irun(irun_in), time(time_in), kode(kode_in) {}
    IKode(const IsigI& Is_in, const int& kode_in, const float& invres)
      : Is(Is_in), sSqr(invres), irun(-1), time(-1.0), kode(kode_in) {}


    IsigI Is;  // intensity/sigI
    float sSqr;  // 1/d^2
    int irun;
    float time;
    int kode; // packed hkl code
  };
  //--------------------------------------------------------------
  bool operator < (const IKode& a,const IKode& b);
  // for sort on resolution

  //======================================================================
  class ValCount
  {
  public:
    ValCount(const Rtype v, const int n) : val(v), count(n) {}

    Rtype val;
    int count;
  };

  //======================================================================
  class correl_coeff
  {
  public:
    // Constructor, initialise all sums to zero
    correl_coeff()
      : sum_wx(0.0), sum_wx2(0.0), sum_wy(0.0),
	sum_wy2(0.0), sum_wxy(0.0), sum_w(0.0),
	sum_wwx(0.0), sum_wwy(0.0), sum_w2x(0.0), sum_w2y(0.0),
	n(0), nw(0) {};

    void zero(); // initialise all sums to zero

    // Add in contribution x, y with weight w
    void add(const Rtype& x, const Rtype& y, const Rtype& w=1.0);
    // Add in contribution x.I(), y.I() with weight w (for IsigI correlation)
    void add(const IsigI& x, const IsigI& y, const Rtype& w=1.0);
    // Add in contribution x, y with weight w &
    // sdx,y =  wx, wy
    void add(const Rtype& x, const Rtype& y,
	     const Rtype& w, const Rtype& wx, const Rtype& wy);
    // Number of contributions
    int Number() const;
    // set CC as dummy
    void SetCC();
 
    // Calculate correlation coefficient from sums
    ValCount result() const;

    void dump() const;

    correl_coeff& operator+=(const correl_coeff& other);
    friend correl_coeff& operator+ (const correl_coeff& a, const correl_coeff& b);

  private:
    double sum_wx, sum_wx2, sum_wy, sum_wy2, sum_wxy, sum_w;
    double sum_wwx, sum_wwy, sum_w2x, sum_w2y;
    int n, nw;
  };
  //======================================================================
  class MSdiff
  // Root Mean square weighted pair deviation
  // Sqrt (Sum( (I1-I2)^2 / (sigI1**2 + sigI2**2 [ + VarK])))
  // Two intensities with their sds: I1,sigI1; I2,sigI2
  // VarK   optional estimate of varance of relative scale
  {
  public:
    MSdiff() : sum_ndf2(0.0), sum_w(0.0), n_f(0) {}

    void zero();

    // add in contribution 
    void add(const IsigI& Is1, const IsigI& Is2, const double& w=1.0,
	     const double& VarK=0.0);
  
    ValCount result() const
    {
      if (n_f == 0)
	{return ValCount(0.0,0);}
      else
	{return ValCount(-sqrt(sum_ndf2/sum_w), n_f);}
    }
    
    MSdiff& operator+=(const MSdiff& other);
    friend MSdiff& operator+ (const MSdiff& a, const MSdiff& b);

  private:
    double sum_ndf2;
    double sum_w;
    int n_f;
  };
  //======================================================================
  class Rfactor
  {
  public:
    Rfactor() : sum_df(0.0), sum_f(0.0), n_f(0) {}

    // add in contribution with weight w
    void add(const double& df, const double& f, const double& w)
    {sum_df += fabs(df) * w; sum_f += f; n_f += 1;}
  
    ValCount result() const    {
      if (n_f == 0) {
	return ValCount(0.0,0);
      } else {
	return ValCount(sum_df/sum_f, n_f);
      }
    }
     Rtype R() const    {
      if (n_f == 0) {
	return 0.0;
      } else {
	return sum_df/sum_f;
      }
    }

    // Set dummy R = 0.0
    void SetR() {n_f=0;}

    Rfactor& operator+=(const Rfactor& other);
    friend Rfactor& operator+ (const Rfactor& a, const Rfactor& b);

  private:
    double sum_df, sum_f;
    int n_f;
  };

  //======================================================================
  class PairSet
  {
  public:
    PairSet() {};
    PairSet(const int& I1, const int& I2, const int& SymElement);    // weight = 1.0
    PairSet(const int& I1, const int& I2, const int& SymElement, 
	    const double& wt, const float& fac1, const float& fac2);


    // Add in pair if set empty or it matches existing pair(s)
    // matches if has same symmetry element & at least
    // one common observation
    // Returns true if added
    bool AddPair(const int& I1, const int& I2, const int& SymElement,
		 const float& fac1, const float& fac2);

    // Access
    int SymElement() const {return SymElmt;}
    std::vector<IndexPair> Pairs() const {return pairs;}
    std::vector<RPair> Fac() const {return fac12;}
    std::vector<int> ObsIndices() const {return obsindices;}
    double Weight() const {return weight;}
    int Nobs() const {return obsindices.size();}

    void dump()
    {
      std::cout << "PairSet: \n"
		<< "  SymEl = " << SymElmt << "\n"
		<< "  Npair = " << pairs.size() << "\n";
      for (size_t i=0;i<pairs.size();i++) std::cout << pairs[i].first << " - "
						 << pairs[i].second << "\n";
      std::cout << "  Nobs  = " << obsindices.size() << "\n";
      for (size_t i=0;i<obsindices.size();i++) std::cout << obsindices[i] << "\n";

    }


  private:
    std::vector<IndexPair> pairs;
    std::vector<RPair> fac12;       // multiplying scales
    std::vector<int> obsindices;
    double weight;
    int SymElmt;
  };
  //--------------------------------------------------------------
  template<class T> RPair Mean(const std::vector<T>& list)
  // Return mean & SD of list of numbers
  {
    int n = list.size();
    double sum_sc = 0.0;
    double sum_sc2 = 0.0;
    double sc;
    double sd;
    int count = 0;
    
    if (n<1) return RPair(0.0,0.0);
    if (n==1) 
      {
	sc = list[0];
	sd = 0.0;
      }
    else
      {
	for (int i = 0;i<n;i++)
	  {
	    sc = list[i];
	    sum_sc  += sc;
	    sum_sc2 += sc*sc;
	  }
	// Mean
	sc = sum_sc/double(n);
	// SD
	sd = sqrt((sum_sc2 - double(n)*sc*sc)/double(n-1));
      }
    
    return RPair(sc, sd);
  }
  //--------------------------------------------------------------
  template<class T> class Array3D
  // General 3D array class, indexed from 0 or from offset
  {
  public:
    // Dummy constructor 
    inline Array3D() : nu(0),nv(0),nw(0),u0(0),v0(0),w0(0), offset(false) {data.clear();}
    // Construct with size, offset 0,0,0
    inline Array3D(const int& iu, const int& iv, const int& iw)
      : nu(iu),nv(iv),nw(iw),  u0(0),v0(0),w0(0), offset(false)
    {data.resize(iu*iv*iw);}
    // Construct with size (iu,iv,iw) and offset (ku,kv,kw)
    inline Array3D(const int& iu, const int& iv, const int& iw,
		   const int& ku, const int& kv, const int& kw)
      : nu(iu),nv(iv),nw(iw),  u0(ku),v0(kv),w0(kw)
    {data.resize(iu*iv*iw);
      offset=false;if(ku!=0||kv!=0||kw!=0)offset = true;}

    void Offset(const int& ku, const int& kv, const int& kw)
    // Give index of first element (default 0,0,0)
    {u0=ku; v0=kv; w0=kw; offset=false;if(ku!=0||kv!=0||kw!=0)offset = true;}

    // Construct with size, offset 0,0,0
    inline void init(const int& iu, const int& iv, const int& iw, const T& val)
    {
      nu=iu; nv=iv; nw=iw;  u0=0; v0=0; w0=0;
      offset = false;
      data.assign(iu*iv*iw, val);
    }
    
    inline const T& operator ()(const int& i, const int& j, const int& k) const
    { 
      if (offset) return data[i-u0+nu*(j-v0+nv*(k-w0))];
      return data[i+nu*(j+nv*k)]; }      //!< get element
    
    inline T& operator ()(const int& i, const int& j, const int& k)
    {
      if (offset) return data[i-u0+nu*(j-v0+nv*(k-w0))];
      return data[i+nu*(j+nv*k)]; }      //!< set element
    
    inline const T& operator ()(clipper::Vec3<int> v) const
    { 
      if (offset) return data[v[0]-u0+nu*(v[1]-v0+nv*(v[2]-w0))];
      return data[v[0]+nu*(v[1]+nv*v[2])]; }      //!< get element
    
    inline T& operator ()(clipper::Vec3<int> v)
    {
      if (offset) return data[v[0]-u0+nu*(v[1]-v0+nv*(v[2]-w0))];
      return data[v[0]+nu*(v[1]+nv*v[2])]; }      //!< get element
    
    // Access as 1D array (raw index from 0 always, no offset)
    inline const T& operator ()(const int& i) const
    {return data[i];}      //!< get element
    
    inline T& operator ()(const int& i)
    {return data[i];}      //!< set element
    
    void Zero() {data.assign(nu*nv*nw, T(0));}
    
    int size() const {return nu*nv*nw;}

    bool Test(const int& i, const int& j, const int& k)
    // test indices: return true if in range
    { 
      if (i<u0 || i-u0>=nu || j<v0 || j-v0>=nv || k<w0 || k-w0>=nw) return false; 
      return true;
    }
    bool Test(const clipper::Vec3<int> v)
    // test indices: return true if in range
    { 
      if (v[0]<u0 || v[0]-u0>=nu || v[1]<v0 || v[1]-v0>=nv ||
	  v[2]<w0 || v[2]-w0>=nw) return false; 
      return true;
    }

  private:
    int nu,nv,nw;   // size of each dimension
    int u0,v0,w0;   // offsets, ie indices of first element
    bool offset;  // true if offset non-zero
    std::vector<T> data;
  };
//--------------------------------------------------------------
class ValidElementObserved
// Class to store data to determine if sufficient of the possible
// symmetry elements in the lattice or in the original (HKLIN)
// Laue group have been observed
{
public:
  ValidElementObserved()
    : numelementLattice(0),
      numelementpresent(0),
      numoriglgelement(0),
      numoriglgelementpresent(0){}

  // NumElementLattice,	  total number of symmetry elements
  //                          in lattice group (including identity)
  // NumElementPresent,	  number of lattice elements Observedd in data
  // NumOrigLGElement,	  number of elements in original laue group
  // NumOrigLGElementPresent  number Observedd in data
  ValidElementObserved(const int& NumElementLattice,
		       const int& NumElementPresent,
		       const int& NumOrigLGElement,
		       const int& NumOrigLGElementPresent)
    : numelementLattice(NumElementLattice),
      numelementpresent(NumElementPresent),
      numoriglgelement(NumOrigLGElement),
      numoriglgelementpresent(NumOrigLGElementPresent){}

  double ObservedInLatticeGroup() const
  {
    if (numelementLattice-1 <= 0) return 1.0;
    return (double(numelementpresent)/double(numelementLattice-1));
  }

  double ObservedInOrigGroup() const
  {
    if (numoriglgelement <= 0) return 1.0;
    return  (double(numoriglgelementpresent)/double(numoriglgelement));
  }

  bool Valid(const double ThresholdLattice = 0.49,
	     const double ThresholdOrig = 0.49) const
  // Returns true if a sufficient proportion of possible
  // symmetry elements have been observed
  {
    if (numelementLattice == 0) return false;   // no data
    //!+
    //    std::cout << "\n*-*-*-* NumSymop in original LG " << numoriglgelementpresent
    //	      << " of " << numoriglgelement << " total " <<  numelementpresent
    //	      << " of " << numelementLattice-1 << "\n\n";
    //!-
    return ((ObservedInLatticeGroup() >= ThresholdLattice) &&
	    (ObservedInOrigGroup() >= ThresholdOrig));
  }

private:
  // Numbers of elements
  int numelementLattice;       // total number of symmetry elements
  // in lattice group (including identity)
  int numelementpresent;     // number of lattice elements Observedd in data
  int numoriglgelement;        // number of elements in original laue group
  int numoriglgelementpresent; // number Observedd in data

}; // ValidElementObserved
  //================================================================
  class PossibleSpaceGroup
  {
  public:
    PossibleSpaceGroup(){}
    // Constructor arguments:
    //   Sname    spacegroup name (as found in Laue group frame)
    //   number   spacegroup number
    //   Rname    spacegroup name in standard setting
    //   Scond    absence conditions 
    //   ScondLG    absence conditions Laue group frame
    //   Sreindex reindex operator from found SG to standard setting
    //   Schiral  chiral flag
    //   probin   probability of systematic absence
    //   ZonesinGroup list of zone numbers in this group
    //
    PossibleSpaceGroup(const std::string& Sname, const int& number,
		       const std::string& Rname,
		       const std::string& Scond, const std::string& ScondLG,
		       const ReindexOp& Sreindex, const Chirality& Schiral,
		       const double& probin, const std::vector<int>& ZonesinGroup);

    void StoreAcenProb(const double& probacen);
    void StoreLaueGroupName(const std::string& LGname) {lauegroupname = LGname;;}
    void StorePointGroupName(const std::string& PGname) {pointgroupname = PGname;;}
    void StoreLaueGroupProb(const double& problg, const double& confidence);
    void StoreSGname(const std::string& SGname) {name = SGname;}
    void Accept() {accepted=true;}
    void StoreOrigReindex(const ReindexOp& LGreindex) {reindex_orig = LGreindex;}
    void StoreReindex(const ReindexOp& Reindex) {reindex = Reindex;}
    // Return zone list
    std::vector<int> ZoneList() const {return zonesingroup;}

    // Access
    // basic spacegroup name if SGname true, else reference space group name 
    std::string Name(const bool& nameSG=true) const;
    std::string Refname() const {return refname;};    // Reference spacegroup name
    int Sgnumber() const {return sgnumber;}           // Spacegroup number
    // List of reflection conditions
    //  eg "0kl: k+l = 2n, h00: h = 2n"
    // basic spacegroup name if SGname true, else reference space group name 
    std::string Condition(const bool& nameSG=true) const;
    // List of reflection conditions Laue group frame
    std::string ConditionLG() const  {return conditionLG;} 
    // If true, original reindex
    // else if false [default], total reindex operator from original to reference setting
    ReindexOp Reindex(const bool& Orig=false) const;
    // Original reindex
    ReindexOp ReindexOrig() const {return reindex_orig;}
    Chirality Chiral() const {return chiral;}         // Chirality
    // Total probability
    double Prob() const {return prob;} 
    // systematic absence probability
    double SysAbsProb() const {return sysabsprob;} 
    // centrosymmetric or non-centro probability
    double CentroProb() const {return centroprob;} 
    // Laue group probability
    double LaueGroupProb() const {return lauegroupprob;}
    // Laue group confidence
    double LaueGroupConfidence() const {return lauegroupconfidence;}
    std::string LaueGroupName() const {return lauegroupname;}
    std::string PointGroupName() const {return pointgroupname;}

    bool Accepted() const {return accepted;}
    // Return true if P,I or F orthorhombic (to allow for non-"reference" setting option
    bool IsPIForthorhombic() const;
    // Return true if group is I2 etc (mI class)
    bool IsI2() const;

    // Reindex rhombohedral H lattice to rhombohedral R setting
    //  OrigIsR  true if original lattice was R
    // Do nothing if not H lattice
    // Return true if done
    bool ReindexRhombohedral(const bool& OrigIsR);

  private:
    std::string name;         // Spacegroup name
    std::string refname;      // Reference spacegroup name
    int sgnumber;             // Spacegroup number
    std::string lauegroupname;  // Lauegroup name
    std::string pointgroupname;  // Pointgroup name
    std::string condition;    // List of reflection conditions
    //  eg "0kl: k+l = 2n, h00: h = 2n"
    std::string conditionLG;    // List of reflection conditions in Laue group frame
    // reindex operator from Laue group to reference setting
    ReindexOp reindex;
    // reindex operator original setting to Laue group
    ReindexOp reindex_orig;
    Chirality chiral;         // Chirality
    double prob;              // Total probability 
    double lauegroupprob;     // Laue group probability
    double lauegroupconfidence;     // Laue group confidence
    double sysabsprob;        // systematic absence probability
    double centroprob;        // probability of having the right centrosymmetricity
    bool accepted;
    std::vector<int> zonesingroup;   // list of zone numbers in this group

    // Compare names, for find
    friend bool operator == (const PossibleSpaceGroup& a,
			     const PossibleSpaceGroup& b)
    {return a.name == b.name;}
    friend bool operator != (const PossibleSpaceGroup& a,
			     const PossibleSpaceGroup& b)
    {return a.name != b.name;}
    // Sort order on Prob, largest first
    friend bool operator < (const PossibleSpaceGroup& a,
			    const PossibleSpaceGroup& b)
    {return a.prob > b.prob;}
  };  // class PossibleSpaceGroup
  //------------------------------------------------------
  // check all subgroups, reindex H lattices to R
  // Only use if original lattice type was R
  void ResetRhombohedralHexLatticestoR
  (std::vector<scala::PossibleSpaceGroup>& groups);
 }  // namespace scala
#endif

// end  scala_datatypes.h
