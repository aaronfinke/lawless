// pointgroup.hh
//
// This class uses CCtbx routines, but keeps them internal
// All imports & exports are done as standard types (std::vector etc)
// or scala types
// No CCtbx types are used in the external interface

#ifndef SYM_POINTGROUP
#define SYM_POINTGROUP

#include <cctbx/crystal/symmetry.h>
#include <cctbx/sgtbx/lattice_symmetry.h>
#include <cctbx/uctbx/fast_minimum_reduction.h>


using namespace cctbx;

#include <iostream>
#include "hkl_datatypes.hh"
#include "lattice.hh"
#include "Output.hh"
#include "controls.hh"
#include "zone.hh"
#include "score_datatypes.hh"

using scala::Chirality;

namespace CCtbxSym
{
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  //! Fix up rhombohedral space group names for CCTBX
  //  If lattice type  (1st character of name) is R or H, change to "R xxx :H||R"
  std::string CCTBX_SGsymbol_HorR(const std::string& SName);

  std::string PointGroupName(const std::string& SGname);
  // Get equivalent Laue (Patterson) name from spacegroup name
  std::string LaueGroupName(const std::string& SGname);
  void PrintCctbxSymops(const sgtbx::space_group& SpSgrp);
  void PrintChBOp(const sgtbx::change_of_basis_op& ChBasis);
  std::string UcellFormat(const uctbx::unit_cell uc);
  std::string ChangeBasisFormat_as_Reindex(const sgtbx::change_of_basis_op& ChB);
  // returns true if groups are enantiomers
  bool EnantiomorphicGroups(const std::string& SGname1, const std::string& SGname2);
  sgtbx::change_of_basis_op
      MakeChangeOfBasisOp(const scala::ReindexOp& reindex_op);
  // Convert change-of-basis operator to ReindexOp
  scala::ReindexOp SetReindexOp(const sgtbx::change_of_basis_op& ChB);
  std::string Reindex_as_xyz(const scala::ReindexOp& reindex);
  //--------------------------------------------------------------
  // Returns change of basis to "best" cell
  // only affects triclinic, monoclinic & orthorhombic cells
  sgtbx::change_of_basis_op GetBestCell(const sgtbx::space_group& Group,
					const uctbx::unit_cell& uccell,
					const int& AllowI2);
  //==================================================================
  //==================================================================
  class LatticeGroup
  // Lattice group corresponding to space or point group
  // includes crystal system
  {
  public:
    LatticeGroup(){};
    LatticeGroup(const sgtbx::space_group& Group);
    LatticeGroup(const sgtbx::space_group& Group, const char LatType);
    // accessors
    CrystalSystem crystal_system() const {return CrysSys;}
    sgtbx::space_group lattice_group() const {return LatGroup;}
    bool Valid() const {return valid;}
    
  private:
    void init(const sgtbx::space_group& Group,
	      const char LatType); // for constructors
    sgtbx::space_group SpaceGroup;
    sgtbx::space_group LatGroup;
    CrystalSystem CrysSys;
    bool valid;
  };
  //==================================================================
  //==================================================================
  class AlternativeBases
  // Vector of symmetry-related change of basis operators +
  //  "difference" in cell dimensions from some reference cell
  {
  public:
    AlternativeBases(){};
    AlternativeBases(const sgtbx::change_of_basis_op& cb_op,
		     const double& celldiff);
    AlternativeBases(const std::vector<sgtbx::change_of_basis_op>& cb_ops,
		     const double& celldiff);
    
    void AddOp(const sgtbx::change_of_basis_op& cb_op);
    
    double CellDiff() const {return CellDiff_;}
    
    int Nop() const {return CbOps.size();}
    sgtbx::change_of_basis_op FirstOp() const {return CbOps[0];}
    sgtbx::change_of_basis_op Op(const int& i) const {return CbOps[i];}
    // Return "simplest" operation:
    //  identity or one with smallest non-zero elements
    sgtbx::change_of_basis_op SimplestOp() const;
    
    friend bool operator < (const AlternativeBases& a, const AlternativeBases& b)
      // for sorting by rank on celldiff
    {return (a.CellDiff_ < b.CellDiff_);}
    
  private:
    std::vector<sgtbx::change_of_basis_op> CbOps;
    double CellDiff_;
  };
  
  //==================================================================
  //==================================================================
  class PointGroup
  // Laue PointGroup 
  {
  public:
    PointGroup(); // construct as P1
    PointGroup(const std::string& Name); // construct from point-group name
    // Construct as "LatType"1 (C1, I1, F1 ...)
    PointGroup(const char& LatticeType);
    // Construct from two 3x3 matrices representing two operators
    // also recording the element numbers from which it was
    // constructed
    PointGroup(const int& Kelement1,
	       const std::vector<double>& Rmatrix1,
	       const int& Kelement2,
	       const std::vector<double>& Rmatrix2,
	       const char& LatticeType);
    // Construct from PointGroup, with flag to select whether to convert to
    // equivalent LatticeGroup
    PointGroup(const bool SetLatticeGroup, const PointGroup& PG);
    // Add this symmetry element to list if it is present in pointgroup,
    // if it is return true
    bool AddElement(const int& Kelement,
		    const std::vector<double>& Rmatrix);

    // Store cell, and change-of-basis from cell frame
    // to this frame (from reindex operator)
    double SetCell(const std::vector<double>& cellin,
		   const scala::ReindexOp& reindex_op,
		   const int& AllowI2);

    // <<<< Access methods >>>>
    // Returns true if element is present in element list
    bool HasElement(const int& Kelement) const;
    bool HasElement(const std::vector<double>& Rmatrix) const;
    // Return list of symmetry elements 
    std::vector<int> Elements() const {return ElementNums;}
    // Return number of symmetry elements 
    int NElements() const {return ElementNums.size();}

    // Return Laue group name for reference setting (standard)
    std::string RefLGname() const;
    // Return Laue group name for constructor setting
    std::string LGname() const;
    // Return point group name for reference setting (standard)
    std::string RefPGname() const;

    // Return reindex operator from cell to standard setting
    // as 3x3 matrix (in vector in order 11, 12, 13, etc)
    // This operator C can be used to reindex old hkl as
    //   h' = h C   (row vectors h, h')
    scala::ReindexOp RefSGreindex() const;
    // as string
    std::string RefSGreindexFormat() const;
    // Reindex operator constructor to reference 
    scala::ReindexOp SGreindex() const;
    // Reindex operator original to constructor
    scala::ReindexOp SGreindexOrig() const;


    CrystalSystem crystal_system() const;

    // returns maximum angular deviation from imposing     
    // symmetry constraints
    double Delta() const {return delta;}

    // Return cell in standard (reference) frame
    std::vector<double> TransformedCell() const;

    void PrintAlternativeCells(phaser_io::Output& output,
			       const bool& OutputXML,
			       const float& max_delta,
			       const int& AllowI2) const;

    // Return a number indicating the "order" of the group,
    // number of operators, SG number etc
    // ie higher value is higher symmetry
    int OrderValue() const;

    // Equality of operators in original (cell) frame
    // (same rotation operators in reference frame &&
    //   same change-of-basis)
    bool Equals(const PointGroup& other) const;
    // Equality of operators in reference frame
    bool EqualsRef(const PointGroup& other) const;

    // Friends
    friend bool operator == (const PointGroup& a,const PointGroup& b)
    // Test only Laue group, not unitcell
    {return (a.RotGrp == b.RotGrp);}

    //    {return ((a.LaueGrp_ref == b.LaueGrp_ref) && (a.ChBasis_ref.c() == b.ChBasis_ref.c()));}

    friend bool operator != (const PointGroup& a,const PointGroup& b)
      // Test only Laue group, not unitcell
    {return (a.RotGrp != b.RotGrp);}
    //    {return (!(a.LaueGrp_ref == b.LaueGrp_ref) && (a.ChBasis_ref.c() == b.ChBasis_ref.c()));}

    friend std::vector<scala::ReindexOp>
        AlternativeIndexing(const PointGroup& PG,
			    const bool& strict,
			    const scala::Scell target_cell,
			    const float& max_delta,
			    const int& AllowI2);

    void dump() const;

    std::vector<scala::PossibleSpaceGroup> TestSpaceGroupList
    (std::vector<scala::Zone>& Zones,
     const Chirality& chiral) const;

    // List of compatible spacegroups
    std::vector<std::string> SpaceGroupList(const Chirality chiral=scala::CHIRAL,
					    const bool AllLattice = false) const;



    void PrintAllSpacegroups(phaser_io::Output& output,
			     const Chirality chiral=scala::CHIRAL,
			     const bool AllLattice = false) const;

    void SetCellDiff(const double& diff) {Cell_Diff = diff;}
    double CellDiff() const;

    char GetLatType() const;

    // Find list of cells close to target cell cell_target,
    //   within diff_tolerance
    // Return reindex operators
    // Also deviations
    std::vector<scala::ReindexOp>
    GetCloseCell(const scala::Scell& cell_target,
		 const double& diff_tolerance,
		 const bool& ExcludeIdentity, 
		 std::vector<double>& celldiff,
		 const int& AllowI2);

  private:
    // For constructors
    void init(const sgtbx::space_group& Pgroup, const char& LatticeType);
    sgtbx::space_group GetSpaceGroup() const {return RotGrp_ref;}

    scala::SysAbsScore TestPossible(const cctbx::sgtbx::space_group& sg,
			     std::vector<scala::Zone>& Zones) const;

    // There are three relevant basis frames:
    //  1) frame used in constructor, ie frame of symmetry operators
    //  2) frame in which the unit cell is defined, original frame _cell
    //  3) standard (reference) frame _ref

    // acentric pointgroup without lattice translations
    sgtbx::space_group RotGrp;       // Constructor basis
    sgtbx::space_group RotGrp_ref;   // reference frame
    sgtbx::space_group PntGrp_ref;   // RotGrp + lattice centering, reference frame
    // Laue group (including inversion & lattice centering) in
    // constructor (lattice) frame
    sgtbx::space_group_type LaueGrp_type;  // Constructor frame
    sgtbx::space_group LaueGrp_ref;     // reference frame
    sgtbx::space_group_type LaueGrp_ref_type;     // reference frame


    // Cell in original frame
    std::vector<double> input_cell;
    uctbx::unit_cell uccell;

    // Change of basis (Orig(cell)->Constructor)
    scala::ReindexOp CellReindexOp;
    sgtbx::change_of_basis_op ChBasis_cell;

    // Change of basis Original (cell) -> reference
    sgtbx::change_of_basis_op ChBasis;
    uctbx::unit_cell uccell_ref;

    // Constructor -> reference
    sgtbx::change_of_basis_op ChBasis_ref;

    // Reference -> "best"
    sgtbx::change_of_basis_op ChBasis_best;

    // Symmetry element list
    std::vector<int> ElementNums;

    char LatType; // 'R' not 'H' for rhombohedral, CCtbx convention

    double delta;

    double Cell_Diff;
  };
  //--------------------------------------------------------------
  // Get list of alternative basis sets for Laue group
  // subject to criteria set by flags
  //
  // On entry:
  //  Group       symmetry group
  //  ExcludeIdentity if true, exclude settings identical or
  //                   symmetry-related to initial setting
  //  AnyCell     if true, accept any cell (no test against refcell)
  //              if false, accept only cells within tolerance of
  //              given refcell, ranked by difference
  //  BestCell    if true, only accept "best" conventional cell
  //              starting from uccell
  //              Not relevent for symmetry > Orthorhombic
  //              definition of "best" depends on AnyCell
  //                AnyCell true  - best according to standard
  //                AnyCell false - closest to refcell
  //
  //  uccell      initial unit cell
  //  refcell     reference unit cell for comparison
  //  tolerance   tolerance on mean square base vector difference
  //  angular_tolerance  similarity allowed for
  //                     monoclinic & triclinic cells
  //  AllowI2      > 0 to add I2 settings to C2 (mC or mI)
  //
  // On exit:
  //  returns vector of vectors of change-of-basis operators
  //  grouped if they are related by the Group symmetry
  
  // Best Cell definitions:
  //   1) orthorhombic
  //      test all permutations of axes, standard has a<b<c
  //   3) monoclinic
  //      test affine transformations around unique axis
  //      best has minimum beta - 90
  //   3) triclinic
  //      get reduced cell
  //      
  //
  std::vector<AlternativeBases>
  GetAlternativeBases(const sgtbx::space_group& Group,
		      const bool ExcludeIdentity, bool AnyCell,
		      const bool BestCell,
		      const uctbx::unit_cell& uccell,
		      const uctbx::unit_cell& refcell,
		      const double& tolerance,
		      const double& angular_tolerance,
		      const int& AllowI2);

  //--------------------------------------------------------------

  // Return list of possible alternative indexing schemes
  // compatible with cell or with lattice symmetry
  //
  // (1) If strict = true
  //   return symmetry operators from symmetry elements 
  //   present in lattice group but not in pointgroup
  //   (max_delta is ignored)
  //   This can only happen for symmetries > orthorhombic
  // else (2) strict = false
  //   find accidental alternatives arising from special
  //   cell dimension relationships
  //   max_delta is angular tolerance for cell similarity
  //
  //  Returns:-
  //   List of reindex operators, including:-
  //     strict flag: true if "strict"
  //     deviations - cell deviations from initial cell
  //                  (maximum angle in degrees)
  //                  = 0.0 for strict settings
  //
  // Friend of class PointGroup
  std::vector<scala::ReindexOp> 
  AlternativeIndexing(const PointGroup& PGz1,
		      const bool& strict,
		      const scala::Scell target_cell,
		      const float& max_delta,
		      const int& AllowI2);

  // Get list of possible change of basis operators, depending
  // on symmetry group
  //
  // On entry:
  //  LatGroup     lattice group
  //  Group        symmetry group
  //  AllowI2      >0  to add I2 settings to C2 (mC or mI)
  //
  // 1) Symmetry > orthorhombic
  //    find operators in LatGroup which are not present in Group
  // 2) Symmetry == orthorhombic
  //    axis permutation operators (independent of Group)
  // 3) Symmetry == monoclinic or triclinic
  //    find all affine transformations which preserve point group
  //
  std::vector<sgtbx::change_of_basis_op> PossibleChBOp
  (const  LatticeGroup& LatGroup, const sgtbx::space_group& Group,
   const int& AllowI2);

}
#endif
