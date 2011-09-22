// \file hkl_symmetry.hh

#ifndef HKL_SYMMETRY_HEADER
#define HKL_SYMMETRY_HEADER

#include <vector>
// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_info;

#include "hkl_datatypes.hh"
#include "lattice.hh"

namespace scala {
  //--------------------------------------------------------------
  //! Change H to R or vv in space group name
  //! If 1st character of name is H or R, change to :
  //!   1. character after ":" if H or R
  //!   2. HorR if 'H' or 'R' [default 'H']
  //!   3. else leave
  //! eg  R 3 2 -> H 3 2
  //!     R 3 :H -> H 3
  std::string SGnameHtoR(const std::string& sgname, const char& HorR = ' ');
  //--------------------------------------------------------------
  //! Mostly same as (and derived from) clipper::Spacegroup

  //! but preserves the order of symmetry operators, as needed for
  //! unmerged MTZ files
  class SpaceGroup : public clipper::Spacegroup
		     // Clipper Spacegroup with symops in defined order
		     // ie same as in constructor
  {
  public:
    // Constructors
    SpaceGroup() : Nsymp(0) {}
    SpaceGroup(std::vector<clipper::Symop>& symops); //!< constructor from symops
    SpaceGroup(const std::string& spgname);          //!< constructor from name
    SpaceGroup(const int& SpgNumber);                //!< constructor from number
    SpaceGroup(const clipper::Spacegroup& ClpSG);    //!< constructor from clipper

    void init(std::vector<clipper::Symop>& symops);//!< initialise from symops 
    void init(const std::string& spgname);	   //!< initialise from name	  
    void init(const int& SpgNumber);		   //!< initialise from number 
    void init(const clipper::Spacegroup& ClpSG);   //!< initialise from clipper

    //! True if no valid space group has been constructed
    bool Null() const {return (Nsymp == 0);}

    clipper::Symop Symop(const int& symN) const
    {return csymops.at(symN);}  //!< return symN'th symop
    clipper::Symop RotSymop(const int& symN) const
    {return rotsymops.at(symN);} //!< return symN'th rotational symop
    clipper::Symop InvRotSymop(const int& symN) const
    {return invrotsymops.at(symN);}  //!< return symN'th inverse symop

    //! return the symop corresponding to Isym from put_in_asu
    clipper::Symop SymopFromIsym(const int& isym) const;

    //! return symop codes for rotation symmetry operators
    clipper::Symop_code RotSymopCode(const int& symN) const
    {return clipper::Symop_code(rotsymops.at(symN));}
    //! return inverse symop codes for rotation symmetry operators
    clipper::Symop_code InvRotSymopCode(const int& symN) const
    {return clipper::Symop_code(invrotsymops.at(symN));}

    //! return derived Patterson (Laue) group
    SpaceGroup PattersonGroup() const;
    //! return derived point group
    SpaceGroup PointGroup() const;

    //! return lattice type character
    char LatType() const {return lattype;}
    //! Return group with new lattice, keeping point group operators
    SpaceGroup NewLatticePointGroup(const char& Lattype) const;

    //! return all codes for ordered operators from clipper group
    std::vector<clipper::Symop_code> SymopCodes() const;

    bool IsSymopIdentity(const int& symN) const; //!< true if symN'th operator is identity

    std::string formatAllSymops_as_xyz() const; //!< real-space symops
    std::string formatAllSymops_as_hkl() const; //!< reciprocal symops

    //! Space group name
    std::string Symbol_hm() const {return spacegroupname;}
    //! Space group name
    std::string Symbol_hall() const;
    //! Space group number
    int Spacegroup_number() const {return spacegroupnumber;}
    //! CCP4 Space group number
    int CCP4_Spacegroup_number() const {return CCP4spacegroupnumber;}

    //! put hkl into asymmetric unit, return ISYM symmetry number (even for -h-k-l)
    clipper::HKL put_in_asu(const clipper::HKL& hkl, int& isym) const;
    //! get hkl from asymmetric unit, given ISYM symmetry number (even for -h-k-l)
    clipper::HKL get_from_asu(const clipper::HKL& hkl, const int& isym) const;

    //! Compare whole symops (including translations)
    friend bool operator == (const SpaceGroup& a,const SpaceGroup& b);
    friend bool operator != (const SpaceGroup& a,const SpaceGroup& b);

    //! Change basis, reindex
    void ChangeBasis(const scala::ReindexOp& reindex);
    
  private:
    int Nsymp;
    char lattype;  // 'H' for hexagonal-setting rhombohedral
    std::string spacegroupname;
    int spacegroupnumber;
    int CCP4spacegroupnumber;

    // The following operators are in constructor-order
    // (unlike clipper::Spacegroup which may reorder them
    // Note that if [Ri] is the rotation part of the i'th symmetry operator,
    // then in reciprocal space we have indices actual h(i) and reduced h(red)
    //   h(red) = h(i)T   [Ri]
    //   h(i)   = h(red)T [Ri]^-1 
    std::vector<clipper::Symop> csymops;      // full space group operators
    std::vector<clipper::Symop> rotsymops;    // rotation operators [Ri]
    std::vector<clipper::Symop> invrotsymops; // inverse rotation operators [Ri]^-1

    // Private functions
    void init();
    // Return lattice type extracted from space group name & symops
    void SetLatType();
    // CCP4 symop from clipper symop
    CSym::ccp4_symop MakeCCP4symop(const clipper::Symop& S) const;
    // Use CCP4 library routine to get name from operators
    std::string CCP4spaceGroupName() const;
    // Use CCP4 library routine to get numbers from operators
    void CCP4spaceGroupNumber();
    // Clipper symop from CCP4 symop
    clipper::Symop MakeClippersymop(const CSym::ccp4_symop& S) const;
    // Make symop string from CCP4 space group
    std::string AllSymopsfromCCP4(const CSym::CCP4SPG* ccp4sg) const;
  };  // end class SpaceGroup
  //--------------------------------------------------------------
  class SymElement {
    //! Symmetry element: an N-fold rotation 
  public:
    int Nfold;  // order of rotation
    clipper::Vec3<int> iaxis;  // rotation axis direction vector
    std::vector<int> symops;  // list of symop numbers representing this element
    SymElement() : Nfold(0) {};  // null constructor
    void print() const;
  };
  //--------------------------------------------------------------
  class hkl_symmetry
  //
  //! All required reciprocal-space symmetry stuff
  //
  //  Constructor
  //    hkl_symmetry(const scala::SpaceGroup& SG)
  //    hkl_symmetry(const std::string SpgName);
  //    hkl_symmetry(const int& SpgNumber);
  //    hkl_symmetry(const clipper::Spacegroup& ClpSG);
  //
  //  Reset symmetry
  //    void set_symmetry()
  //
  //  Constructor & set_symmetry setup tables to identify symmetry
  //  elements (N-fold rotations) and list corresponding symmetry
  //  operators
  //  Also sets up the operation corresponding to pairs of ISYM codes
  //  ( ISYM code = SymopNumber*2 + 1(hkl) or +2(-h-k-l))
  //
  //  An n-fold rotation axis corresponds to n operators including I
  //  I is a 1-fold rotation
  //
  //  int get_symelmt(isym1, isym2) returns symmetry element number
  //                               (from 0) of primitive operator
  //                                relating a pair ofobservations
  //                                with isym1 & isym2
  //
  {
  public:
    // Constructors
    hkl_symmetry();
    hkl_symmetry(const clipper::Spacegroup& cspgp); //!< constructor from clipper
    hkl_symmetry(const scala::SpaceGroup& cspgp); //!< constructor form SpaceGroup
    hkl_symmetry(const std::string SpgName); //!< constructor from name
    hkl_symmetry(const int& SpgNumber); //!< constructor from number

    //! true if object is null
    bool IsNull() const {return (Nsymp == 0);}
    //! store reindex operator
    void set_reindex(const ReindexOp& op) {ReindexToZeroFrame = op;}

    //! return SpaceGroup
    SpaceGroup GetSpaceGroup() const {return spaceGroup;}
    //! return clipper spacegroup
    clipper::Spacegroup ClipperGroup() const {return spaceGroup;}

    //! put hkl into asymmetric unit, return ISYM symmetry number (odd for h+)
    Hkl put_in_asu(const Hkl& hkl, int& isym) const; 
    //! put hkl into asymmetric unit, return ISYM symmetry number (odd for h+)
    clipper::HKL put_in_asu(const clipper::HKL& hkl, int& isym) const;
    //! get hkl from asymmetric unit, given ISYM symmetry number (odd for h+)
    Hkl get_from_asu(const Hkl& hkl, const int& isym) const;
    bool is_centric(const Hkl& hkl) const; //!< true if reflection centric
    float epsilon(const Hkl& hkl) const;  //!< epsilon, = 0 for systematic absences

    bool IsCentro() const {return centro;}  //!< true if space group centrosymmetric
    bool IsChiral() const {return chiral;}  //!< true if space group chiral

    //! true if hkl is present in lattice, false if systematically absent from lattice
    bool LatticePresent(const Hkl& hkl) const; 
    //! change basis, reindex
    void ChangeBasis(const scala::ReindexOp& reindex);

    // *****
    // ***  Symmetry element stuff
    //! return number of symmetry elements 
    int Nelement() const {return Nelement_;}
    //! return true if kelement'th element is identity
    bool IsElementIdent(const int& kelement) const;
    //! return axis direction for kelement'th symmetry element
    clipper::Vec3<int>  AxisDirection(const SymElement& elmt) const;
    //! make formatted version of symmetry element
    std::string format_element(const int& kelement) const;
    //! make XML version of symmetry element
    std::string XML_element(const int& kelement) const;
    void print_element(const int& idx) const; //!< print idx'th symmetry element (debug)
    void print_elements() const; //!< print all symmetry elements (debug)


    //! return number of operators in kelement'th symmetry element
    int NopInElement(const int& kelement) const;
    //! Return rotation part of j'th inverse symmetry operator of kelement'th symmetry element
    std::vector<double> SymopInElement(const int& j,
				       const int& kelement) const;
    //! Return lsym'th inverse symmetry operator of  kelement'th symmetry element
    clipper::Symop ClipperSymopInElement(const int& lsym,
					 const int& kelement) const;


    //! Return symmetry-element number (from 0->Nelement-1) relating observations with ISYM =  isym1,isym2
    /*! Isym = 2 * SymNumber + 1  for h+
     Isym = 2 * SymNumber + 2  for h-
     Ignore difference between h+ and h- for this purpose */
    int get_symelmt(const int& isym1, const int& isym2)  const;

    //! Return list of symmetry elements (from 0) relating original indices h1 & h2
    /*! This call is used for centric reflections where more than one
     symmetry operator relates two observations
     Ignore difference between h+ and h- for this purpose */
    std::vector<int> get_symelmt(const Hkl& h1, const Hkl& h2) const;

    // Return spacegroup names etc
    std::string symbol_Hall() const;        //!< return Hall symbol
    std::string symbol_xHM(const char& HorR=' ') const;  //!< return extended HM symbol
    std::string symbol() const {return spgname;} //!< return name from constructor or from xHM

    //! return lattice type (P, C, I, F, H, ?R)
    char lattice_type() const {return LatType;}
    //! return number of primitive symops, == 0 for null object
    int NsymP() const {return Nsymp;}
    //! return number of symops
    int Nsym() const {return spaceGroup.num_symops();}

    //! return all rotational symops in element
    std::vector<clipper::Symop> RotSymopsInElement(const int& kelement) const;
    //! return all primitive rotational symops
    std::vector<clipper::Symop> PrimRotSymops() const;


    //! return crystal system
    CrystalSystem CrysSys() const {return cryssys;}
    //! return cell constraint flags, vector(6)
    /*!  = -1  anything, = 0 constrained angle (=90 or 120)
      > 0 = j, same as j'th parameter */
    std::vector<int> CellConstraint() const;

    //! test equality with another symmetry object
    /*!  not equals_r & equals_rt test number of operators &  lattice type (centering) */
    bool equals_rt(const hkl_symmetry& other) const; //!< test all operators including translation
    bool equals_r(const hkl_symmetry& other) const; //!< test just rotation parts
    // true if spacegroup other is same ignoring translations
    // with operators in the same order
    bool equals_r_order(const hkl_symmetry& other) const;

    //!  simple test, ignore translations
    bool operator==(const hkl_symmetry& other) const {return equals_r(other);}
    //!  simple test, ignore translations
    bool operator!=(const hkl_symmetry& other) const {return !equals_r(other);}

    //! true if jel'th element of this symmetry object == jel2'th element of other
    bool equalElement(const int& jel, const hkl_symmetry& other, const int& jel2) const;

  private:
    int Nsymp;
    SpaceGroup spaceGroup; //  the actual space group
    std::string spgname;  // name from constructor or xHM name
    char LatType;
    int Nelement_;
    CrystalSystem cryssys;
    bool centro;
    bool chiral;

    std::vector<clipper::Symop_code> inv_rot_symcodes;
    std::vector<clipper::Symop_code> rot_symcodes;
    std::vector<SymElement> elements; // symops for each sub-set
    std::vector<int> element_index;  // element number for this symop (Nsymp)
    std::vector<int> symelmt;   // a 2D array Nelement x Nelement

    // This is the reciprocal space operator to convert back to
    // initial frame
    // Reindex operator H 
    // Applies to index h such that h'T = hT H
    ReindexOp ReindexToZeroFrame; // initially identity

    // private functions
    bool ElementEqual(const SymElement& a, const SymElement& b) const;
    void set_symmetry(); // set up symmetry stuff on initialisation
    // Get crystal system
    CrystalSystem CrysSysfromSGnumber(const int& SpaceGroupNumber);
    // return true if space group is chiral
    bool ChiralTest();
  }; 
  //--------------------------------------------------------------
}
#endif
