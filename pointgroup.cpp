// pointgroup.cpp

#include "pointgroup.hh"

//#include <cctbx/crystal/symmetry.h>
//#include <cctbx/sgtbx/lattice_symmetry.h>
#include <cctbx/sgtbx/find_affine.h>
#include <cctbx/sgtbx/rot_mx_info.h>
#include <cctbx/sgtbx/group_codes.h>

//#include <clipper/clipper.h>

#include "latsym.hh"
#include "hkl_symmetry.hh"
#include "getsubgroups.hh"
#include "zone.hh"
#include "scala_util.hh"
#include "string_util.hh"
#include "cctbx_utils.hh"

/*
#include "lattice.hh"
#include "pgscore.hh"
#include "hkl_merged_list.hh"
*/

using namespace cctbx;
using phaser_io::LOGFILE;
using phaser_io::LXML;


#include <assert.h>
#define ASSERT assert

namespace CCtbxSym
{
  //--------------------------------------------------------------
  std::string CCTBX_SGsymbol_HorR(const std::string& SName)
  // Fix up rhombohedral space group names for CCTBX
  //  If lattice type  (1st character of name) is R or H, change to "R xxx :H||R"
  {
    std::string name = StringUtil::Trim(SName);
    if (name.size() > 0) {
      // First character should be a lattice type
      char lattype = toupper(name[0]);
      if (lattype == 'H' || lattype == 'R') { // only do anything if H or R
	std::vector<std::string> parts = StringUtil::split(name, ":");
	if (parts.size() > 1) {
	  // extension :H or :R found
	  std::string ext = StringUtil::Trim(parts[1]);
	  if (!(ext == "H" || ext == "R")) {
	    // extension must be H or R
	    clipper::Message::message(clipper::Message_fatal
		("CCTBX_SGsymbol_HorR: Illegal space group name "+SName));
	  }
	} else {
	  // No extension, so add one, :R or :H
	  name += ":"+std::string(1,lattype);
	}
	// Change first character to R
	name = "R"+name.substr(1);
      }
    }
    return name;
  }
  //--------------------------------------------------------------
  // Get equivalent pointgroup name from spacegroup name
  std::string PointGroupName(const std::string& SGname)
  {
    sgtbx::space_group SG = sgtbx::space_group
      (sgtbx::space_group_symbols(CCTBX_SGsymbol_HorR(SGname)).hall());
    return SpaceGroupName(SG.build_derived_reflection_intensity_group(true).type(), 'H');
  }
  //--------------------------------------------------------------
  // Get equivalent Laue (Patterson) name from spacegroup name
  std::string LaueGroupName(const std::string& SGname)
  {
    sgtbx::space_group SG = sgtbx::space_group
      (sgtbx::space_group_symbols(CCTBX_SGsymbol_HorR(SGname)).hall());
    return SpaceGroupName(SG.build_derived_patterson_group().type(), 'H');
  }
  //--------------------------------------------------------------
  bool EnantiomorphicGroups(const std::string& SGname1, const std::string& SGname2)
  // returns true if groups are enantiomers
  {
    sgtbx::space_group SG1(sgtbx::space_group_symbols(SGname1).hall());
    sgtbx::space_group_type SGT1(SG1);
    if (!SGT1.is_enantiomorphic()) return false;
    sgtbx::space_group SG2(sgtbx::space_group_symbols(SGname2).hall());
    sgtbx::space_group_type SGT2(SG2);
    if (!SGT2.is_enantiomorphic()) return false;
    // change hand of SG2 & check for equality with SG1
    sgtbx::change_of_basis_op cb = SGT2.change_of_hand_op();
    sgtbx::space_group SG2m = SG2.change_basis(cb);
    if (SG2m == SG1)
      {return true;}
    return false;
  }
  //--------------------------------------------------------------
  std::string Reindex_as_xyz(const scala::ReindexOp& reindex_op)
  {
    return StringUtil::Strip("("+MakeChangeOfBasisOp(reindex_op).as_xyz()+")");
  }
  //--------------------------------------------------------------
  std::string DeColonName(const std::string& name)
    // remove any part of spacegroup name after ":" inclusive
    // Change 'R' to 'H'
  {
    std::string::size_type i = name.find(":");
    std::string s = name;
    if (i < name.size()) {
      // ":" found
	if (name.find("Hall") < name.size()) {
	  std::string::size_type j = name.find("(");
	  // Strip off starting "Hall:" & trailing "(" string 
	  if (name[j-1] == ' ') j--;
	  s = name.substr(i+2,j-i-2);
	} else {
	  if (name[i-1] == ' ') i--;
	  s = name.substr(0,i);
	}
    }
    return scala::SGnameHtoR(s,'H');
  }
  //--------------------------------------------------------------
  char CentringSymbol(const sgtbx::space_group& Group)
    // Extract conventional lattice symbol from spacegroup
    // For rhombohedral lattices, set it to 'R' even in rhombohedral
    // axis system (cctbx returns 'P' in this case)
  {
    char LT = 'P';
    if (Group.crystal_system() == sgtbx::crystal_system::trigonal)  {
      // Extract code from Hermann-Mauguin symbol, first non-blank
      std::string HMcode = sgtbx::space_group_symbols(Group.type().number())
	.hermann_mauguin();
      for (size_t i=0;i<HMcode.size();i++) {
	if (HMcode[i] != ' ') {
	  LT = HMcode[i];
	  break;
	}
      }
    } else {
      LT = Group.conventional_centring_type_symbol();
    }
    return LT;
  }
  //--------------------------------------------------------------
  bool GlidePresent(const sgtbx::space_group& sg, const scala::Zone& glidezone,
		    const scala::ReindexOp& LatToSG)
  // Return true if glidezone is present in spacegroup
  //   LatToSG is reindex operator from lattice (zone reference) frame to
  //   space group reference frame
  {
    std::vector<scala::Hkl> Chkl = glidezone.GlideTestHkl(LatToSG); // pair of test hkl's
    // Both reflections must be absent for glide to be present
    // in spacegroup
    bool sysabs = true;
    for (size_t i=0;i<Chkl.size();i++) {
      if (!sg.is_sys_absent
	  (miller::index<int>(Chkl[i].h(),Chkl[i].k(),Chkl[i].l())))
	{sysabs = false;}
    }
    return sysabs;
  }
  //--------------------------------------------------------------
  int ScrewPresent(const sgtbx::space_group& sg, const scala::Zone& axiszone,
		    const scala::ReindexOp& LatToSG)
    // Return > 0 if axiszone is present in spacegroup
    //  with value depending on which screw component found
    //  eg for 6-fold
    //   = 1 for 6(3), = 2 for 6(2|4), = 3 for 6(1|5)
  //   LatToSG is reindex operator from lattice (zone reference) frame to
  //   space group reference frame
  {
    // list of test hkl's and test conditions
    std::vector<std::vector<int> > screwabsence;
    std::vector<scala::Hkl> Chkl = axiszone.AxisTestHkl(LatToSG, screwabsence);
    int nr = Chkl.size();
    int ip = 0;
    int np = axiszone.NgridPoints();
    if (nr > 0) {
      ASSERT (nr == 2);
      ASSERT (int(screwabsence.size()) == np);
      for (ip=0;ip<np;ip++) {
	if (screwabsence[ip][0] >= 0) {
	  // skip if no test
	  bool accept = true;
	  for (int i=0;i<2;i++) {
	    //^
	    //^		    std::cout << "ScrewPresent:" << Chkl[i].h() << " "
	    //^      << Chkl[i].k() << " " << Chkl[i].l() <<"\n";
	    //^-
	    bool absent = sg.is_sys_absent
	      (miller::index<int>(Chkl[i].h(),Chkl[i].k(),Chkl[i].l()));
	    // if absent in spacegroup, must be flagged as absent ie = +1
	    // if present in spacegroup, must be flagged as present ie = 0
	    // else not found
	    if (!((absent && screwabsence[ip][i] > 0) ||
		  (!absent && screwabsence[ip][i] == 0))) accept = false;
	  }
	  if (accept) break;
	}
      }
    }
    if (ip < 0 || ip >= np) ip = 0;
    return ip;
  }
  //--------------------------------------------------------------
  sgtbx::change_of_basis_op GetBestCell(const sgtbx::space_group& Group,
					const uctbx::unit_cell& uccell,
					const int& AllowI2)
  // Returns change of basis to "best" cell
  // only affects triclinic, monoclinic & orthorhombic cells
  {
    bool ExcludeIdentity = false, AnyCell = true, BestCell = true;
    std::vector<AlternativeBases> CbOp_list = 
      GetAlternativeBases(Group, ExcludeIdentity, AnyCell, BestCell,
			  uccell, uccell, 3.0, 0.0001, AllowI2);
    // should return only one value
    return  CbOp_list[0].FirstOp();
  }
  //--------------------------------------------------------------
  LatticeGroup::LatticeGroup(const sgtbx::space_group& Group)
  {
    init (Group, CentringSymbol(Group));
  }
  //--------------------------------------------------------------
  LatticeGroup::LatticeGroup
  (const sgtbx::space_group& Group, const char LatType)
  {
    init (Group, LatType);    
  }
  //--------------------------------------------------------------
  void LatticeGroup::init(const sgtbx::space_group& Group,
			  const char LatType)
      // Determines lattice group corresponding to spacegroup
    // (from lookup) 
    // and crystal system CrysSys (eg MONOCLINIC)
    // International Tables spacegroup number
    //   195-230 cubic
    //   168-194 hexagonal
    //   143-167 trigonal
    //   75-142  tetragonal
    //   16-74   orthorhombic
    //   3-15    monoclinic
    //   1-2     triclinic
  {
    valid = false;
    SpaceGroup = Group;
    int SpaceGroupNumber = SpaceGroup.type().number();
    if ( ! scala::AllowedLatticeType(LatType))
      {
	clipper::Message::message(clipper::Message_fatal
		  ("LatticeGroup: Illegal lattice-type "+std::string(1,LatType)));
      }

    std::string LatGrpSymbol = "";

    if (SpaceGroupNumber > 194) {
      // Cubic
      CrysSys = CUBIC;
      if (LatType == 'P')
	LatGrpSymbol = "P m -3 m";
      else if (LatType == 'I')
	LatGrpSymbol = "I m -3 m";
      else if (LatType == 'F')
	LatGrpSymbol = "F m -3 m";
    } else if (SpaceGroupNumber > 142) {
      // Hexagonal or trigonal
      if (SpaceGroupNumber > 167)
	CrysSys = HEXAGONAL;
      else
	CrysSys = TRIGONAL;
      if (LatType == 'P')
	LatGrpSymbol = "P 6/m m m";
      else if (LatType == 'R' || LatType == 'H')
	LatGrpSymbol = "R -3 m :H";
    } else if (SpaceGroupNumber > 74) {
      // Tetragonal
      CrysSys = TETRAGONAL;
      if (LatType == 'P')
	LatGrpSymbol = "P 4/m m m";
      else if (LatType == 'I')
	LatGrpSymbol = "I 4/m m m";
    } else if (SpaceGroupNumber > 15) {
      // Orthorhombic
      CrysSys = ORTHORHOMBIC;
      if (LatType == 'P')
	LatGrpSymbol = "P m m m";
      else if (LatType == 'I')
	LatGrpSymbol = "I m m m";
      else if (LatType == 'C')
	LatGrpSymbol = "C m m m";
      else if (LatType == 'F')
	LatGrpSymbol = "F m m m";
    } else if (SpaceGroupNumber > 2) {
      CrysSys = MONOCLINIC;
      if (LatType == 'P')
	LatGrpSymbol = "P 1 2/m 1";
      else if (LatType == 'C')
	LatGrpSymbol = "C 1 2/m 1";
      else if (LatType == 'I')
	LatGrpSymbol = "I 1 2/m 1";
    } else {
      // Triclinic
      CrysSys = TRICLINIC;
      LatGrpSymbol = "P -1";
    }
    if (LatGrpSymbol != "") {
      LatGroup = sgtbx::space_group
	(sgtbx::space_group_symbols(LatGrpSymbol).hall());
      valid = true;
    }
  }
  //--------------------------------------------------------------
  AlternativeBases::AlternativeBases(const sgtbx::change_of_basis_op& cb_op,
				     const double& celldiff)
    : CbOps(1,cb_op), CellDiff_(celldiff) {}
  //--------------------------------------------------------------
  AlternativeBases::AlternativeBases
  (const std::vector<sgtbx::change_of_basis_op>& cb_ops,
   const double& celldiff)
    : CellDiff_(celldiff) {CbOps = cb_ops;}
  //--------------------------------------------------------------
  void AlternativeBases::AddOp(const sgtbx::change_of_basis_op& cb_op)
  {
    // Special for identity operator, make sure it is the first one
    if (cb_op.is_identity_op())
      CbOps.insert(CbOps.begin(), cb_op);
    else
      CbOps.push_back(cb_op);
  }
  //--------------------------------------------------------------
  sgtbx::change_of_basis_op AlternativeBases::SimplestOp() const
  {
    // Return "simplest" operation:
    //  identity or one with smallest non-zero elements
    if (CbOps[0].is_identity_op())
      return CbOps[0];
    float besttot = 1000000.;
    int minneg = 100000;
    int kbest = 0;
    float tol = 0.01;  // tolerance
    for (size_t k=0;k<CbOps.size();k++) {
      float total = 0.0;
      int nneg = 0;
      scitbx::mat3<double> op = CbOps[k].c_inv().r().as_double();
      // Criteria for "simplest"
      // 1) smallest sum of absolute value of all elements
      // 2) minimum number of negatives
      for (int i=0;i<3;i++) {
	for (int j=0;j<3;j++) {
	  total += std::abs(op(i,j));
	  if (op(i,j) < 0.0) nneg += 1;
	}
      }
      if (total < besttot-tol) {
	besttot = total;
	kbest = k;
	minneg = nneg;
      } else if (Close<float>(total, besttot, tol)) {
	if (nneg < minneg) {
	  besttot = total;
	  kbest = k;
	  minneg = nneg;
	}
      }
    }
    return CbOps[kbest];
  }
  //--------------------------------------------------------------
  void AppendUniqueChBasis
  (const sgtbx::change_of_basis_op& cb_op,
   const sgtbx::rt_mx& cb_op_mx,
   const double& celldiff,
   const sgtbx::space_group& Group,
   std::vector<AlternativeBases>& ChBasisVec)
    // Add cb_op into ChBasisVec (vector of AlternativeBases)
    //   note that the AlternativeBases class is itself a vector of ops
    // If it is equivalent by symmetry in Group (rotational part), add
    //  to equivalent vector, else start new vector
    // Don't store if identical to existing op
    // Special for P1: no check on existing operators, just append
    //
    // On entry:
    //  cb_op      change of basis operator
    //  cb_op_mx   equivalent rt_mx operator
    //  Group      symmetry group
    //  ChBasisVec vector to update
  {
    int kop = 0;
    bool equiv_op = false;
    bool ident = false;
    if (ChBasisVec.size() > 0 && Group.n_smx() > 1) {
      // do we have a symmetry-related version of this one already?
      for (size_t j=0;j<Group.n_smx();j++) {
	sgtbx::rt_mx S = Group.smx(j).multiply(cb_op_mx).
	  new_denominators(cb_op_mx);
	for (size_t k=0;k<ChBasisVec.size();k++) {
	  if (S ==  ChBasisVec[k].FirstOp().c().new_denominators(S)) {
	    if (j == 0) {
	      // Identical
	      ident = true;
	    }
	    equiv_op = true;
	    kop = k;
	    break;
	  }
	}
      }		
    }
    if (!ident) {
      if (equiv_op) {
	// Equivalent operator, accumulate in vector
	ChBasisVec[kop].AddOp(cb_op);
      } else {
	// New operator
	ChBasisVec.push_back(AlternativeBases(cb_op,celldiff));
      }
    }
  }
  //--------------------------------------------------------------
  void PrintChBOp(const sgtbx::change_of_basis_op& ChBasis)
  {
    printf("Change of basis: %s\n",ChBasis.as_xyz().c_str());

    printf("        Inverse: %s\n",ChBasis.as_hkl().c_str());
  }
  //--------------------------------------------------------------
  std::string UcellFormat(const uctbx::unit_cell uc)
  {
    clipper::String line;
    for (int i=0;i<6;i++) line += " "+clipper::String(uc.parameters()[i],6,5);
    return line;
  }
  //--------------------------------------------------------------
  void PrintCctbxSymops(const sgtbx::space_group& SpSgrp)
  {
    for (size_t i=0;i<SpSgrp.n_ltr();i++) {
      scitbx::vec3<double> vl = SpSgrp.ltr(i).as_double();
      printf("Lattice translation: ");
      for (int j=0;j<3;j++) printf(" %8.3f",vl[j]);
      printf("\n");
    }
    af::shared<sgtbx::rt_mx> symops = SpSgrp.all_ops();
    for (size_t i=0;i<symops.size();i++) {
      printf("Symmetry operator %3d: %s\n",int(i+1),symops[i].as_xyz().c_str());
    } 
  }   
  //--------------------------------------------------------------
  sgtbx::change_of_basis_op
  MakeChangeOfBasisOp(const scala::ReindexOp& reindex_op)
    // Return change_of_basis_op based on reindex operator (its inverse)  
  {
    // Reindex Op
    clipper::Mat33<double> R = reindex_op.rot();
    scitbx::mat3<double> CM(R(0,0), R(0,1), R(0,2),
			    R(1,0), R(1,1), R(1,2),
			    R(2,0), R(2,1), R(2,2));

    sgtbx::change_of_basis_op C(sgtbx::rt_mx(CM.inverse(),
					     scitbx::vec3<double>(0.0,0.0,0.0),
					     sgtbx::cb_r_den,sgtbx::cb_t_den));
    return C;
  }
  //--------------------------------------------------------------
  scala::ReindexOp SetReindexOp(const sgtbx::change_of_basis_op& ChB)
    // Convert change-of-basis operator to ReindexOp
  {
    scitbx::mat3<double> op = ChB.c_inv().r().as_double();
    clipper::Mat33<double> vop;
    for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
	vop(i,j) = op(i,j);
    
    scala::ReindexOp vp(vop);
    return vp;
  }
  //--------------------------------------------------------------
  std::string ChangeBasisFormat_as_Reindex(const sgtbx::change_of_basis_op& ChB)
  {
    scitbx::mat3<double> op = ChB.c_inv().r().as_double();
    return MVutil::FormatReindex_as_hkl(op);
  }
  //--------------------------------------------------------------
  sgtbx::space_group LaueGroup(const sgtbx::space_group& RotGrp,
			       const char& LatType)
    // Make Laue group from rotation group by adding inversion centre
    // and lattice centring
  {
    sgtbx::space_group Pgroup(RotGrp);
    Pgroup.expand_inv(sgtbx::tr_vec(0,0,0));  // add inversion centre
    if (LatType != ' ') {
      Pgroup.expand_conventional_centring_type(LatType);
    }
    Pgroup.make_tidy();
    return Pgroup;
  }
  //--------------------------------------------------------------
  bool OperatorInGroup(sgtbx::change_of_basis_op& ChBasis,
		       const sgtbx::space_group& SG)
    // Return true if change-of-basis operator is a member of 
    // the spacegroup
  {
    sgtbx::space_group::smx_array_type smx_list = SG.smx();
    // Denominator of rotation part of change-of-basis operator
    int den = ChBasis.c().r().den();

    for (size_t i=0;i<smx_list.size();i++)
      {
	if (smx_list[i].r().new_denominator(den) == ChBasis.c().r())
	  return true;
      }
    return false;
  }
  //--------------------------------------------------------------
  int MonoclinicUniqueAxis(const sgtbx::space_group& Group)
    // Returns unique axis for monoclinic spacegroup
    // = 0,1,2 for a,b,c
    // = -1 if not monoclinic
    // from Ralf Grosse-Kunstleve's lattice_symmetry.cpp
  {
    if (Group.n_smx() != 2) {return -1;}
    // Second symmetry operator (1st is identity)
    sgtbx::rot_mx_info two_fold_info(Group(1).r());
    // must be dyad
    ASSERT(scitbx::fn::absolute(two_fold_info.type()) == 2);
    // Axis direction (eigenvector), the element we want = +1
    // Assumed to be some permutation of (0,0,1)
    scitbx::vec3<int> const& ev = two_fold_info.ev();
    ASSERT(std::count(ev.begin(), ev.end(), 0) == 2);
    // Axis = 0,1,2
    return static_cast<int>(std::find(ev.begin(), ev.end(), 1) - ev.begin());
  }
  //--------------------------------------------------------------
  bool SymInGroup(const sgtbx::space_group& Group,
		  const sgtbx::rt_mx& cb_op_mx)
  {
    // Returns true if ChBoperator is in group (including identity)
    for (size_t j=0;j<Group.n_smx();j++)
      {
	if (cb_op_mx == Group.smx(j))
	  {
	    return true;
	  }
      }
    return false;
  }
  //--------------------------------------------------------------
  std::vector<sgtbx::change_of_basis_op> PossibleChBOp
  (const  LatticeGroup& LatGroup, const sgtbx::space_group& Group,
   const int& AllowI2)
  // Get list of possible change of basis operators, depending
  // on symmetry group
  //
  // On entry:
  //  LatGroup     lattice group
  //  Group        symmetry group
  //  AllowI2      +1 to add I2 settings to C2 (mC or mI)
  //
  // 1) Symmetry > orthorhombic
  //    find operators in LatGroup which are not present in Group
  // 2) Symmetry == orthorhombic
  //    axis permutation operators (independent of Group)
  // 3) Symmetry == monoclinic or triclinic
  //    find all affine transformations which preserve point group
  //    (plus I2 alternative to C2 if AllowI2 > 0)
  //
  {

    std::vector<sgtbx::change_of_basis_op> cb_ops;
    sgtbx::change_of_basis_op cb_op_def;  // to define correct denominators

    CrystalSystem CrysSys = LatGroup.crystal_system();

    if (CrysSys > ORTHORHOMBIC) {    
      // Higher symmetry than orthorhombic
      cb_ops.push_back(sgtbx::change_of_basis_op());  // add in identity
      // Test if there any symmetry operators which are in
      // the lattice group but not the Laue group 
      // These are potential alternative bases
      sgtbx::space_group LatGrp =  LatGroup.lattice_group();
      for (size_t i=1;i<LatGrp.n_smx();i++) {
	// Is this operator in Laue group?
	bool found = false;
	for (size_t j=1;j<Group.n_smx();j++) {
	  if (LatGrp.smx(i) == Group.smx(j)) {
	    found=true; break;
	  }
	}
	if (!found) {
	  // Add operator into vector as change of basis
	  cb_ops.push_back(sgtbx::change_of_basis_op(LatGrp.smx(i)));
	}
      }
    } else if (CrysSys == ORTHORHOMBIC) {
      // Orthorhombic
      cb_ops.push_back(sgtbx::change_of_basis_op());  // add in identity
      //  Use cubic symmetry to generate permutation operators
      sgtbx::space_group affine_group("P 4 3*");
      for (std::size_t i_smx=1;i_smx<affine_group.n_smx();i_smx++) {
	cb_ops.push_back(sgtbx::change_of_basis_op
			 (affine_group(i_smx)).new_denominators(cb_op_def));
      }
    } else if (CrysSys == MONOCLINIC || CrysSys == TRICLINIC) {
      int ChooseI2 = 0;
      // Is this centred monoclinic?
      if (CrysSys == MONOCLINIC && (CentringSymbol(Group) == 'C')) {
	if (AllowI2 > 0) {
	  ChooseI2 = +1;  // allow I2
	} else if (AllowI2 < 0) {
	  ChooseI2 = -1;  // force I2
	}
      }

      // List of change-of-basis operators which leave
      // spacegroup invariant and leave cell the same within tolerance

      const int affine_range=4;  // default = 2
      if (ChooseI2 >= 0) {
	sgtbx::find_affine affine(Group, affine_range);
	af::const_ref<sgtbx::rt_mx> affine_cb_mx = affine.cb_mx().const_ref();
	
	// Loop and store potential change-of-basis operators
	for(std::size_t i_cb_mx=0;i_cb_mx<affine_cb_mx.size();i_cb_mx++) {
	  cb_ops.push_back(sgtbx::change_of_basis_op
			   (affine_cb_mx[i_cb_mx]).new_denominators(cb_op_def));
	}
      }
      if (ChooseI2 !=0) {
	// Change of basis from C2 to I2
	sgtbx::space_group SG_C2 = sgtbx::space_group(sgtbx::space_group_symbols("C2/m"));
	sgtbx::space_group SG_I2 = sgtbx::space_group(sgtbx::space_group_symbols("I2/m"));
	sgtbx::change_of_basis_op I2toC2 = SG_I2.type().cb_op().new_denominators(cb_op_def);
	sgtbx::change_of_basis_op C2toI2 = I2toC2.inverse().new_denominators(cb_op_def);
	
	sgtbx::find_affine affineI2(Group.change_basis(C2toI2), affine_range);
	af::const_ref<sgtbx::rt_mx> affine_cb_mxI2 = affineI2.cb_mx().const_ref();
	
	// Loop and store potential change-of-basis operators
	for(std::size_t i_cb_mx=0;i_cb_mx<affine_cb_mxI2.size();i_cb_mx++)  {
	  sgtbx::change_of_basis_op new_cbop =
	    sgtbx::change_of_basis_op(affine_cb_mxI2[i_cb_mx]).new_denominators(cb_op_def);
	  // Post-multiply each operator in I2 frame by I2 to C2 transformation
	  cb_ops.push_back(new_cbop * I2toC2);
	}
      }
    }
    return cb_ops;
  }
  //--------------------------------------------------------------
  bool UpdateBestCell(const bool& sameGroup,
		      const CrystalSystem& CrysSys,
		      const uctbx::unit_cell test_cell,
		      const int& UniqueAxis,
		      const sgtbx::change_of_basis_op& cb_op,
		      const double& angular_tolerance,
		      sgtbx::change_of_basis_op& best_cb_op,
		      uctbx::unit_cell& best_cell)
  // If test_cell is closer to the "standard" cell than
  // best_cell, update best_cell := test_cell and best_cbop := cb_op
  // Monoclinic & orthorhombic (not C-centred, detected by sameGroup == false)
  //  only, no change otherwise
  // Return true if accepted
  {
    int compare = -1;
    if (CrysSys == ORTHORHOMBIC && sameGroup) {
      compare = best_cell.compare_orthorhombic(test_cell);
    } else if (CrysSys == MONOCLINIC) {
      compare = best_cell.compare_monoclinic
	(test_cell, UniqueAxis, angular_tolerance);
    }

    if (compare > 0) {
      //^
      //	std::cout << "Updating cell: " << CrysSys << " "
      //		  << sameGroup << " " << angular_tolerance << "\n"
      //		  << "Old best cell: " << UcellFormat(best_cell) <<"\n"
      //		  << "New best cell: " << UcellFormat(test_cell) <<"\n";
      //-!
      best_cb_op = cb_op;
      best_cell = test_cell;
    }
    return (compare > 0);
  }
  //--------------------------------------------------------------
  std::vector<AlternativeBases>
  GetAlternativeBases(const sgtbx::space_group& Group,
		      const bool ExcludeIdentity, bool AnyCell,
		      const bool BestCell,
		      const uctbx::unit_cell& uccell,
		      const uctbx::unit_cell& targetcell,
		      const double& tolerance,
		      const double& angular_tolerance,
		      const int& AllowI2)
    // Get list of alternative basis sets for Laue group
    // subject to criteria set by flags
    //
    // On entry:
    //  Group       symmetry group
    //  ExcludeIdentity if true, exclude settings identical or
    //                   symmetry-related to initial setting
    //  AnyCell     if true, accept any cell (no test against targetcell)
    //              if false, accept only cells within tolerance of
    //              given targetcell, ranked by difference
    //  BestCell    if true, only accept "best" conventional cell
    //              starting from uccell
    //              Not relevent for symmetry > Orthorhombic
    //              definition of "best" depends on AnyCell
    //                AnyCell true  - best according to standard
    //                AnyCell false - closest to targetcell
    //
    //  uccell      initial unit cell
    //  targetcell  target unit cell for comparison
    //  tolerance   tolerance on mean square base vector difference
    //  angular_tolerance  similarity allowed for
    //                     monoclinic & triclinic cells
    //  AllowI2      > 0 to add I2 settings to C2 (mC or mI)
    //               < 0 force I2 setting
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
    //      If AllowI2 > 0, test I2 setting as well as C2 for centred monoclinic
    //   3) triclinic
    //      get reduced cell
    //      
    //
  {
    //^
    //    std::cout << "GetAlternativeBases " << Group.type().hall_symbol() << "\n";
    //    std::cout << "Cell: " << UcellFormat(uccell) << "\n";
    //-!


    std::vector<sgtbx::change_of_basis_op> ChBops;
    std::vector<AlternativeBases> ChBasisVec;


    sgtbx::change_of_basis_op best_cb_op;

    sgtbx::space_group AltGroup = Group;

    LatticeGroup LatGroup(Group);
    CrystalSystem CrysSys = LatGroup.crystal_system();

    uctbx::unit_cell best_cell = uccell;
    uctbx::unit_cell test_cell = uccell;
    sgtbx::rt_mx cb_op_mx;

    double celldiff;  // difference between two cells

    int UniqueAxis = -1;

    if (CrysSys == MONOCLINIC) {
      // Unique axis = 0,1,2 for a,b,c, = -1 if triclinic
      UniqueAxis = MonoclinicUniqueAxis(Group);
    }

    bool GotPossibleChBOpList = false;


    // First get standard "best" cell in all cases where it is needed
    if (AnyCell) {
      if (CrysSys == TRICLINIC) {
	// triclinic, just get reduced cell
	best_cb_op = ReduceCell(test_cell);
      } else {
	//^ 
	//	    std::cout << "First operator:\n";
	//	    PrintChBOp(best_cb_op);
	//-!
	if (CrysSys == MONOCLINIC || CrysSys == ORTHORHOMBIC) {
	  // Get list of possible change of basis operators to test,
	  // depending on crystal system
	  // This is independent of cell, but will run over a range
	  // (+-affine_range) of cell offsets
	  ChBops = PossibleChBOp(LatGroup, Group, AllowI2);
	  //^
	  //      std::cout << "Number of ChBOp " << ChBops.size() << "\n";
	  //      for (int i=0;i< ChBops.size();++i) {
	  //	PrintChBOp(ChBops[i]);
	  //      }
	  //-!
	  GotPossibleChBOpList = true;
	  int op1 = 0;
	  if (AllowI2 < 0) {
	    // Force I2 option - apply first operator
	    best_cb_op = ChBops[0];
	    best_cell = best_cb_op.apply(uccell);
	    op1 = 1; // skip 1st op
	  }
	  //^
	  //		std::cout << "\n**** Base cell :" << UcellFormat(uccell) << "\n"
	  //			  << "Symmetry " << Group.type().hall_symbol() << "\n";
	  //		std::cout << "First " << " cell "
	  //			  << UcellFormat(best_cell) << "\n";
	  //		PrintChBOp(best_cb_op);
		//^-
		// Loop potential change-of-basis operators
	  for (size_t i=op1;i<ChBops.size();i++) {
	    sgtbx::change_of_basis_op cb_op = ChBops[i];
	    cb_op_mx = cb_op.c().new_denominators(Group.smx(0));
		  
	    AltGroup = Group.change_basis(cb_op);
	    bool sameGroup = (AltGroup == Group);
	    //		    if (AltGroup == Group) {
	    // Exclude operators which change group
	    //  (possibly never happens)
	    test_cell = cb_op.apply(uccell);
	    // Best standard cell
	    bool updated = UpdateBestCell
	      (sameGroup, CrysSys, test_cell, UniqueAxis, cb_op,
	       angular_tolerance,
	       best_cb_op, best_cell);
	    updated = updated;
	    //^
	    //		    if (updated) {
	    //		      std::cout <<"\nUpdated\n";
	    //		      std::cout << "Alternative " << i
	    //				<< " group "
	    //				<< AltGroup.type().universal_hermann_mauguin_symbol()
	    //				<< " lattice "<< AltGroup.conventional_centring_type_symbol()
	    //				<< " cell " << UcellFormat(cb_op.apply(uccell)) << "\n";
	    //		      PrintChBOp(cb_op);
	    //		    } else {
	    //		      std::cout << "Not accepted, groups "
	    //				<< Group.type().hall_symbol() << " : "
	    //				<< AltGroup.type().hall_symbol() << "\n";
	    //		    		    }
	    //-!
	  }
	}
      }
      // We now have a "best" cell
      if (BestCell) {
	// Store best operator as sole solution
	// that is all we need
	ChBasisVec.clear();
	// Exclude identity option
	cb_op_mx = best_cb_op.c().new_denominators(Group.smx(0));
	if (! (ExcludeIdentity && SymInGroup(Group, cb_op_mx)) )
	  {ChBasisVec.push_back(AlternativeBases(best_cb_op,0.0));}
	//^
	//	    std::cout << "**** Best cell " 
	//		      << UcellFormat(best_cb_op.apply(uccell)) << "\n";
	//	    PrintChBOp(best_cb_op);
	//^-
	return ChBasisVec;
      }
      best_cell = best_cb_op.apply(uccell);
    } else {
      best_cell = targetcell;
    }
    
    //^
    //    std::cout << "* * Input cell        " << UcellFormat(uccell) << "\n";
    //    std::cout << "* * Target(best) cell " << UcellFormat(best_cell) << "\n";
    // Get differences against either "best" cell or target cell
    // only get here if !BestCell

    // Get list of possible change of basis operators to test,
    // depending on crystal system, if we haven't got them already
    // This is independent of cell, but will run over a range
    // (+-affine_range) of cell offsets
    if (!GotPossibleChBOpList) {
      ChBops = PossibleChBOp(LatGroup, Group, AllowI2);
    }

    // Loop potential change-of-basis operators
    for(std::size_t i=0;i<ChBops.size();i++) {
      sgtbx::change_of_basis_op cb_op = ChBops[i];
      cb_op_mx = cb_op.c().new_denominators(Group.smx(0));

      if (! (ExcludeIdentity && SymInGroup(Group, cb_op_mx)) ) {
	// Not identity or symmetry-related, if excluded
	AltGroup = Group.change_basis(cb_op);
	if (AltGroup == Group) {
	  // Exclude operators which change group
	  test_cell = cb_op.apply(uccell);
	  // Measure of difference between cells
	  //  root mean square difference of bases
	  celldiff = sqrt(best_cell.bases_mean_square_difference(test_cell));
	  //^
	  //	std::cout << "* * Test cell " << UcellFormat(test_cell) << "\n";
	  if (AnyCell || celldiff < tolerance) {
	    //		    targetcell.is_similar_to
	    //		    (test_cell,length_tolerance,angular_tolerance))
	    //^
	    //		    std::cout << "**** Got it! " << UcellFormat(test_cell) << "\n";
	    //		    PrintChBOp(cb_op);
	    //^-
	    // Cell acceptable on AnyCell or similar
	    // or all acceptable cells (not BestCell)
	    AppendUniqueChBasis(cb_op, cb_op_mx, celldiff, Group, ChBasisVec);
	  }
	}
      }  // end if ExcludeIdentity
    }  // end loop operators
    // Sort solutions of deviation from best or target cell

    std::sort(ChBasisVec.begin(),ChBasisVec.end());

    // for "BestCell" option, just return first example
    if (ChBasisVec.size() > 0 && BestCell) ChBasisVec.resize(1);

    return ChBasisVec;
  }
  //--------------------------------------------------------------
  std::vector<AlternativeBases>
  ChangeBasesList(const std::vector<AlternativeBases>& CBlist,
		  const sgtbx::change_of_basis_op& cb_op)
  // Apply change-of-basis operator to all operators in list
  // to allow for change of frame
  // Post-multiply operators by cb_op
  // Also reorder operators if necessary to put identity operator
  // first in sublist
  {
    std::vector<AlternativeBases> CB_newlist(CBlist.size());
    if (CBlist.size() > 0) {
      for (size_t i=0;i<CBlist.size();i++) {
	int Nops = CBlist[i].Nop();
	std::vector<sgtbx::change_of_basis_op> CbOps;
	for (int j=0;j<Nops;j++) {
	  // Operator relative to constructor frame
	  sgtbx::change_of_basis_op cb_ij =
	    CBlist[i].Op(j).new_denominators(cb_op);
	  
	  sgtbx::change_of_basis_op CB = cb_ij *  cb_op;
	  if (CB.is_identity_op()) {
	    CbOps.insert(CbOps.begin(), CB);
	  } else {
	    CbOps.push_back(CB);
	  }
	}
	CB_newlist[i] = AlternativeBases(CbOps, CBlist[i].CellDiff());
      }
    }
    return CB_newlist;
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  PointGroup::PointGroup()
    // Default constructor as P1
  {
    init(sgtbx::space_group(), 'P');
  }
  //--------------------------------------------------------------
  PointGroup::PointGroup(const char& LatticeType)
    // Construct as "LatType"1 (C1, I1, F1 ...)
  {
    sgtbx::space_group Pgroup;
    init(Pgroup, LatticeType);
  }
  //--------------------------------------------------------------
  PointGroup::PointGroup(const std::string& Name)
    // Construct from pointgroup name
    // Remove any translations, don't add inversion
  {
    sgtbx::space_group_symbols sgsymbol;
    try {  // try name
      sgsymbol = sgtbx::space_group_symbols(CCTBX_SGsymbol_HorR(Name));
    }
    catch (cctbx::error) {
      // Fall back via number (and CCP4 libraries) eg for I 1 21 1
      sgsymbol = sgtbx::space_group_symbols(scala::SpaceGroup(Name).Spacegroup_number());
    }

    sgtbx::space_group Pgroup =
      sgtbx::space_group(sgsymbol.hall()).build_derived_reflection_intensity_group(false);
    //^
    //    std::cout << "PGinitName " << Name << " " << CCTBX_SGsymbol_HorR(Name) << " "
    //	      << sgtbx::space_group_symbols(CCTBX_SGsymbol_HorR(Name)).hall() << " "
    //	      << Pgroup.type().lookup_symbol() << " " << CentringSymbol(Pgroup) << "\n";
    init(Pgroup, CentringSymbol(Pgroup));
  }
  //--------------------------------------------------------------
  PointGroup::PointGroup(const int& Kelement1,
			 const std::vector<double>& Rmatrix1,
			 const int& Kelement2,
			 const std::vector<double>& Rmatrix2,
			 const char& LatticeType)
  {
    ASSERT(Rmatrix1.size() == 9);
    ASSERT(Rmatrix2.size() == 9);
  
    // record unique element numbers
    ElementNums.push_back(Kelement1);
    // Make rt_mx matrices
    // input operators are reciprocal space operators
    // invert to get real-space
    sgtbx::rt_mx R1 = MVutil::SetRtMx(Rmatrix1).inverse();

    sgtbx::rt_mx R2;
    if (Kelement2 != Kelement1) 
      {
	ElementNums.push_back(Kelement2);
	R2 = MVutil::SetRtMx(Rmatrix2).inverse();
      }

    //  initialise space-group, add  1st & 2nd operators
    sgtbx::space_group Pgroup;
    Pgroup.expand_smx(R1);
    if (Kelement2 != Kelement1) Pgroup.expand_smx(R2);
    Pgroup.make_tidy();
    init(Pgroup, LatticeType);
  }
  //--------------------------------------------------------------
  // Construct from PointGroup, with optional flag to convert to
  // equivalent LatticeGroup
  PointGroup::PointGroup(const bool SetLatticeGroup, const PointGroup& PG)
  {
    char LT = PG.GetLatType();
    sgtbx::space_group Pgroup = PG.GetSpaceGroup();
    if (SetLatticeGroup) {
      LatticeGroup LG(Pgroup, LT);
      // Lattice group is in the reference frame, so we need
      // to transform it into the frame of PG
      Pgroup = LG.lattice_group().change_basis(PG.ChBasis.inverse());
    }
    init(Pgroup, LT);
  }
  //--------------------------------------------------------------
  void PointGroup::init(const sgtbx::space_group& Pgroup, const char& LatticeType)
  {
    Cell_Diff = 0.0;
    // Force rotation group acentric (remove any centre)
    RotGrp = Pgroup.build_derived_acentric_group();

    LatType = LatticeType;
    if (LatType == 'H') LatType = 'R';

    // Laue group in constructor frame
    //  Add inversion & centering 
    sgtbx::space_group LaueGrp = LaueGroup(RotGrp, LatType);
    //^
    //    std::cout <<"LaueGrp\n";
    //    PrintCctbxSymops(LaueGrp);
    //^-
    // Reference setting
    LaueGrp_type = LaueGrp.type();
    // change of basis from current setting to reference setting
    ChBasis_ref = LaueGrp_type.cb_op();

    // change of basis from current setting to reference setting
    ChBasis = ChBasis_ref;
    input_cell.resize(6);
    // Don't initialise cell (takes too long to call SetCell!)

    // Rotation group in reference frame
    RotGrp_ref =
      RotGrp.change_basis(ChBasis_ref).build_derived_point_group();
    LaueGrp_ref = LaueGrp.change_basis(ChBasis_ref);
    LaueGrp_ref_type = sgtbx::space_group_type(LaueGrp_ref);
    LatType = CentringSymbol(LaueGrp_ref);
    // Point group, including lattice centering
    PntGrp_ref = RotGrp_ref;
    PntGrp_ref.expand_conventional_centring_type(LatType);

    //^
    //    std::cout << "\nPGinit " <<  LaueGrp_type.lookup_symbol() << " "
    //    	      <<  LaueGrp_ref_type.lookup_symbol() << "\n";
    //    std::cout << "   ChBasis_ref\n";
    //    PrintChBOp(ChBasis_ref);
    //^-
  }
  //--------------------------------------------------------------
  char PointGroup::GetLatType() const
  {
    return CentringSymbol(LaueGrp_ref);
  }
  //--------------------------------------------------------------
  bool PointGroup::AddElement(const int& Kelement,
			      const std::vector<double>& Rmatrix)
    // Test if this symmetry operator is present in pointgroup,
    // if it is return true & add element number to list
  {
    if (HasElement(Rmatrix)) {
      ElementNums.push_back(Kelement);
      return true;
    }
    return false;
  }
  //--------------------------------------------------------------
  bool PointGroup::HasElement(const std::vector<double>& Rmatrix) const
  // Test if this symmetry operator is present in pointgroup,
  // if it is return true
  {
    sgtbx::rt_mx R = MVutil::SetRtMx(Rmatrix).inverse();
    // Compare this with all the primitive operators of the
    // point group
    bool found = false;
    for (size_t k=0;k<RotGrp.n_smx();k++) {
      if (R == RotGrp(0, 0, k)) {
	found = true;
	break;
      }
    }
    return found;
  }
  //--------------------------------------------------------------
  bool PointGroup::HasElement(const int& Kelement) const
    // Returns true if element is present in element list
  {
    if (std::find(ElementNums.begin(), ElementNums.end(), Kelement)
	!=  ElementNums.end()) return true;
    return false;
  }
  //--------------------------------------------------------------
  double PointGroup::SetCell(const std::vector<double>& cellin,
			     const scala::ReindexOp& reindex_op,
			     const int& AllowI2)
    // Store cell, returns maximum angular deviation from imposing
    // symmetry constraints, and store change-of-basis from
    // cell frame to symmetry frame used to construct this object
    // Note the change of basis operator is inverse of reindex operator
    // 
    // reindex_op is reindexing from cell frame to
    // "constructor" (lattice) frame
  {
    ASSERT(cellin.size() == 6);
    input_cell = cellin;
    scitbx::af::double6 dcell;
    for (int i=0;i<6;i++) dcell[i] = cellin[i];
    uccell = uctbx::unit_cell(dcell);

    // Reindex Op from original cell frame to this constructor frame
    CellReindexOp = reindex_op;
    ChBasis_cell = MakeChangeOfBasisOp(CellReindexOp);

    ChBasis = ChBasis_ref * ChBasis_cell;
    uctbx::unit_cell uccell_chb = ChBasis.apply(uccell);
    //^
    //    std::cout << "\n>>> SetCell " << RefLGname() << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";

    // Try to get "best" setting of cell and new reference
    //  cb_op_best transforms reference to "best"
    sgtbx::change_of_basis_op cb_op_best = GetBestCell(LaueGrp_ref, uccell_chb, AllowI2);
    uccell_ref = cb_op_best.apply(uccell_chb);

    //^
    //    std::cout << "Reindex original->constructor ChBasis_cell): "
    //	      <<  ChangeBasisFormat_as_Reindex(ChBasis_cell) << "\n"
    //	      << "Reindex original->reference ChBasis): "
    //	      <<  ChangeBasisFormat_as_Reindex(ChBasis) << "\n"
    //	      << "Reindex constructor->reference ChBasis_ref): "
    //	      <<  ChangeBasisFormat_as_Reindex(ChBasis_ref) << "\n"
    //	      << "CellIn:  " << UcellFormat(uccell) << "\n"
    //	      << "CellChB: " << UcellFormat(uccell_chb) << "\n"
    //	      << "CellRef: " << UcellFormat(uccell_ref) << "\n";
    //-!

    // Test for change of symmetry, C2 to I2
    // Does change of basis operator change the Laue group?
    sgtbx::space_group LGnew = LaueGrp_ref.change_basis(cb_op_best);
    if (LGnew != LaueGrp_ref) {
      // Update stored groups, those in new reference frame
      RotGrp_ref = RotGrp_ref.change_basis(cb_op_best);
      PntGrp_ref = PntGrp_ref.change_basis(cb_op_best);
      LaueGrp_ref = LGnew;
      LaueGrp_ref_type = sgtbx::space_group_type(LaueGrp_ref);
      LatType = CentringSymbol(LaueGrp_ref);
      //^
      //                        std::cout << "\n==== SetCell: updated groups " << LatType << " "
      //                  		<< LaueGrp_ref_type.hall_symbol() << "\n";
      //                        std::cout << "Cell_ref: " << UcellFormat(uccell_ref) << "\n";
      //                        std::cout << "ChBasis:\n";
      //                        PrintChBOp(ChBasis);
      //                        std::cout << "ChBasis_ref:\n";
      //                        PrintChBOp(ChBasis_ref);
      //                        std::cout << "cb_op_best:\n";
      //                        PrintChBOp(cb_op_best);
      //^-
    }
    // ChBasis is operator for Original(cell) -> reference
    ChBasis = cb_op_best * ChBasis;
    ChBasis_ref = cb_op_best * ChBasis_ref;
    //    } else {
    //      // ChBasis is operator for Original(cell) -> reference
    //      ChBasis = cb_op_best * ChBasis;
    //    }
    // If unit cell and symmetry haven't changed, reset
    // change-of-basis to identity
    if (uccell_ref.is_similar_to(uccell)) {
      //	if (OperatorInGroup(ChBasis, RotGrp_ref))
      LatticeGroup LatGroup(RotGrp_ref);
      if (OperatorInGroup(ChBasis, LatGroup.lattice_group())) {
	ChBasis = sgtbx::change_of_basis_op();  // identity
      }
    }
    //^
    //    std::cout << "End of SetCell:\n"
    //	      << "  Reindex original->reference ChBasis): "
    //	      <<  ChangeBasisFormat_as_Reindex(ChBasis) << "\n";
    //^-
    // Maximum angular deviation of cell from that
    // required by rotation group
    delta = sgtbx::lattice_symmetry::
      find_max_delta(uccell_ref,
		     LatticeGroup(RotGrp_ref).lattice_group().build_derived_acentric_group());

    return delta;
  }
  //--------------------------------------------------------------
  void PointGroup::dump() const
  {
    printf("\n>>>>> PG dump <<<<<\n");
    show_space_group_type(LaueGrp_ref_type);
    printf("\n");
    // All possible symbols
    sgtbx::space_group_symbols SGsymbols(LaueGrp_ref_type.number());
    printf("HM symbol %s || extension: %c\n",SGsymbols.hermann_mauguin().c_str(),
	   SGsymbols.extension());

    printf("Change of basis Cell -> Constructor (ChBasis_cell, SGreindexOrig)\n");
    PrintChBOp(ChBasis_cell);
    std::cout << ">>> Cell to standard/reference (ChBasis, RefSGreindex):\n";
    PrintChBOp(ChBasis);
    std::cout << ">>> Constructor to standard/reference (ChBasis_ref, SGreindex):\n";
    PrintChBOp(ChBasis_ref);

    if (input_cell[0] > 0.0) {
      printf("Input cell: ");
      for (int i=0;i<6;i++) printf(" %6.2f",input_cell[i]);
      printf("\n");
      printf("Transformed cell: ");
      printf("%s",UcellFormat(uccell_ref).c_str());
      printf("\n\n");
    }

    // Symmetry elements
    printf("\nRotationGroup:\n");
    PrintCctbxSymops(RotGrp);

    printf("\nRotationGroup in reference frame:\n");
    PrintCctbxSymops(RotGrp_ref);

    printf("Laue Group in reference frame:\n");
    PrintCctbxSymops(LaueGrp_ref);
  }
  //--------------------------------------------------------------
  std::string PointGroup::RefLGname() const
    // Return Laue group name for reference setting (standard)
  {
    return SpaceGroupName(LaueGrp_ref_type, 'H');
  }
  //--------------------------------------------------------------
  std::string PointGroup::LGname() const
    // Return Laue group name for constructor setting
  {
    return SpaceGroupName(LaueGrp_type, 'H');
  }
  //--------------------------------------------------------------
  std::string PointGroup::RefPGname() const
    // Return point group name for reference setting (standard)
  {
    return SpaceGroupName(PntGrp_ref.type(), 'H');
  }
  //--------------------------------------------------------------
  scala::ReindexOp PointGroup::RefSGreindex() const
    // Return change-of-basis operator from cell to standard setting
    // as 3x3 matrix (ReindexOp type)
    // This operator C can be used to reindex old hkl as
    //   h' = h C   (row vectors h, h')
  {
    return SetReindexOp(ChBasis);
  }
  //--------------------------------------------------------------
  scala::ReindexOp PointGroup::SGreindex() const
    // Reindex operator constructor to reference 
  {
    return SetReindexOp(ChBasis_ref);
  }
  //--------------------------------------------------------------
  scala::ReindexOp PointGroup::SGreindexOrig() const
    // Reindex operator original to constructor
  {
    return SetReindexOp(ChBasis_cell);
  }
  //--------------------------------------------------------------
  std::string PointGroup::RefSGreindexFormat() const
    // Return change-of-basis operator to standard setting
    // as string
  {    
    return ChangeBasisFormat_as_Reindex(ChBasis);
  }
  //--------------------------------------------------------------
  std::vector<double> PointGroup::TransformedCell() const
  {
    std::vector<double> tcell(6);
    for (int i=0;i<6;i++) tcell[i] = uccell_ref.parameters()[i];
    return tcell;
  }
  //--------------------------------------------------------------
  int PointGroup::OrderValue() const
    // Return a number indicating the "order" of the group,
    // Just use spacegroup number
  {
    return LaueGrp_ref_type.number();
  }
  //--------------------------------------------------------------
  double PointGroup::CellDiff() const {return Cell_Diff;}
  //--------------------------------------------------------------
  std::vector<scala::ReindexOp>
  PointGroup::GetCloseCell(const scala::Scell& cell_target,
			   const double& diff_tolerance,
			   const bool& ExcludeIdentity, 
			   std::vector<double>& celldiff,
			   const int& AllowI2)
  // Return list of reindex operators & cell differences for
  // alternative indexing schemes which preserve Laue group
  // and have cell close to target.
  //
  // Reindex operators are in target cell frame
  //
  // On entry:
  //   cell_target      target cell in "original(cell)" frame
  //   diff_tolerance   tolerance
  //                     (root mean square difference of base vectors)
  //   ExcludeIdentity  true to exclude identity operator
  //   celldiff         returns cell differences
  //                    (root mean square difference of base vectors)
  //                     (matching reindex operators)
  {
    double diff_max = 10.0;
    scitbx::af::double6 dcell;
    for (int i=0;i<6;i++) dcell[i] = cell_target[i];
    // Target in original frame
    uctbx::unit_cell uccell_target = uctbx::unit_cell(dcell);

    // Get list of alternative bases in "reference" frame
    bool AnyCell = false, BestCell = false;
    double max_delta = 180.0;  // angle tolerance
    // List of bases 
    std::vector<AlternativeBases> CbOp_list = 
      GetAlternativeBases(LaueGrp_ref, ExcludeIdentity, AnyCell, BestCell,
			  uccell_ref, uccell_target,
			  diff_tolerance, max_delta, AllowI2);
    /*    //^
    std::cout << "\nLauegroup ";
    show_space_group_type(LaueGrp_ref_type);
    std::cout << "\nPGCellORIG" << UcellFormat(uccell) << "\n";
    std::cout << "\nCelltargetORIG" << UcellFormat(uccell_target) << "\n";
    std::cout << "\nPGCellREF" << UcellFormat(uccell_ref) << "\n";
    std::cout << "Nop " << CbOp_list.size() << "\n"
	      << "ChBasis original -> reference\n";
    PrintChBOp(ChBasis);   
    std::cout << "ChBasis original -> constructor\n";
    PrintChBOp(ChBasis_cell);   
    //^-    
    */

    double diff0 = -1.0;
    double diff  = -1.0;

    std::vector<AlternativeBases> CloseBases;

    // Things to return
    std::vector<scala::ReindexOp> CloseReindex;
    celldiff.clear();

    if (CbOp_list.size() >= 1) {
      // Change all operators to be relative original cell frame
      std::vector<AlternativeBases> CBlist =
	ChangeBasesList(CbOp_list, ChBasis);
      
      // Best cell difference
      diff0 = CBlist[0].CellDiff();

      //^
      //^std::cout << "CBlist size " << CBlist.size() << "\n";
      
      for (size_t k=0;k<CBlist.size();k++) {
	// Deviation from target
	diff = CBlist[k].CellDiff();
	
	/*//^!
	  std::cout << "***ChBop " << k << "\n";
	  PrintChBOp(CBlist[k].Op(0));
	  std::cout << " Cell" << UcellFormat(CBlist[k].Op(0).apply(uccell))
	  << " Diff: " << diff << "\n";
	  //^*/
	  
	  if (diff < diff_tolerance &&
	      (k == 0 || (diff-diff0) < diff_max)) {
	    // Reindex operator + diff, save for sorting
	    CloseBases.push_back
	      (AlternativeBases(CBlist[k].SimplestOp(), diff));
	  }
      }
      // Sort list & transfer for output
      std::sort(CloseBases.begin(), CloseBases.end());
      
      for (size_t k=0;k<CloseBases.size();k++) {
	CloseReindex.push_back(SetReindexOp(CloseBases[k].FirstOp()));
	celldiff.push_back(CloseBases[k].CellDiff());
      }
    }
    return CloseReindex;
  }
  //--------------------------------------------------------------
  void PointGroup::PrintAlternativeCells(phaser_io::Output& output,
					 const bool& OutputXML,
					 const float& max_delta,
					 const int& AllowI2) const
  {
    bool ExcludeIdentity = true, AnyCell = false, BestCell = false;
    double tolerance = 5.0; 
    std::vector<AlternativeBases> CbOp_list = 
      GetAlternativeBases(LaueGrp_ref, ExcludeIdentity, AnyCell, BestCell,
			  uccell_ref, uccell_ref, tolerance, max_delta, AllowI2);
    // change-of-basis operators are relative to reference cell

    if (CbOp_list.size() > 0) {
      // Change basis of all operators to cell frame
      std::vector<AlternativeBases> CBlist =
	ChangeBasesList(CbOp_list, ChBasis);
      
      for (size_t k=0;k<CBlist.size();k++) {
	uctbx::unit_cell uccell_alt =
	  CBlist[k].Op(0).apply(uccell);
	//  Similarity: tolerance parameters are relative length,
	//    absolute angle
	if (OutputXML)
	  {output.logTabPrintf(0,LXML,"<Alternative Number=\"%3d\">\n", k+1);}
	if (uccell_ref.is_similar_to(uccell_alt, 0.005, 0.1)) {
	  output.logTab(0,LOGFILE,"  Same cell");
	  for (int i=0;i<6;i++) output.logTab(0,LOGFILE,"      ");
	  output.logTab(0,LOGFILE,"         ");
	} else {
	  output.logTab(0,LOGFILE," Other cell");
	  std::vector<double> acell(6);
	  for (int i=0;i<6;i++) {
	    output.logTabPrintf(0,LOGFILE,"%6.1f",uccell_alt.parameters()[i]);
	    acell[i] = uccell_alt.parameters()[i];
	  }

	  output.logTabPrintf(0,LOGFILE," %7.2f ",CBlist[k].CellDiff());
	  if (OutputXML) {
	    output.logTab(1,LXML,scala::Scell(acell).xml());}
	}
	for (int j=0;j<CBlist[k].Nop();j++) {
	  if (j>0) {
	    output.logTabPrintf(0,LOGFILE,"\n                                                          ...");
	  }
	  output.logTabPrintf(0,LOGFILE," %s",
			  ChangeBasisFormat_as_Reindex(CBlist[k].Op(j)).c_str());
	  if (OutputXML) {
	    output.logTabPrintf(1,LXML,"<ReindexNumber> %4d\n", j+1);
	    output.logTab(1,LXML,SetReindexOp(CBlist[k].Op(j)).as_hkl_XML()+"\n");
	    output.logTab(1,LXML,SetReindexOp(CBlist[k].Op(j)).as_XML());
	    output.logTabPrintf(1,LXML,"</ReindexNumber>\n");
	  }
	}
	output.logTab(0,LOGFILE,"\n");
	if (OutputXML)
	  {output.logTab(0,LXML,"</Alternative>");}
      }
    }
  }
  //--------------------------------------------------------------
    CrystalSystem PointGroup::crystal_system() const
    {
      return LatticeGroup(LaueGrp_ref).crystal_system();
    }
  //--------------------------------------------------------------
  std::vector<std::string>
  PointGroup::SpaceGroupList(const Chirality chiral,
			     const bool AllLattice) const
    // Return list of all spacegroups compatible with Laue group (reference frame) 
    // If chiral = CHIRAL, list only chiral spacegroups
    // If AllLattice, include all spacegroups compatible with lattice group
  {
    char CLtype = CentringSymbol(LaueGrp_ref);
    std::vector<std::string> SGnames;

    cctbx::sgtbx::space_group LG_ref = LaueGrp_ref;
    cctbx::sgtbx::space_group LSG;

    if (AllLattice) {
      LatticeGroup LG(LaueGrp_ref, CLtype);
      if (! LG.Valid()) return SGnames;
      LG_ref = LG.lattice_group();
    }

    cctbx::sgtbx::space_group_symbol_iterator sgiter;
    std::string Name;
    for(;;) {
      cctbx::sgtbx::space_group_symbols symbol = sgiter.next();
      if (symbol.number() == 0) break;
      cctbx::sgtbx::space_group sg(symbol.hall());
      bool valid =
	(chiral == scala::CHIRAL && sg.is_chiral()) ||
	(chiral == scala::CENTROSYMMETRIC && sg.is_centric()) ||
	(chiral == scala::NONCHIRAL);
      char LT = CentringSymbol(sg);
      if (AllLattice) {
	LatticeGroup LG(sg, LT);
	if (! LG.Valid())
	  {valid = false;}
	else 
	  {LSG = LG.lattice_group();}
      } else {
	LSG = sg.build_derived_patterson_group();
      }
      if (valid) {
	if (LT == CLtype &&
	    LSG == LG_ref) {
	  // Name = sg.type().lookup_symbol();
	  Name = SpaceGroupName(sg.type(), 'H');
	  if (std::find(SGnames.begin(), SGnames.end(), Name) == SGnames.end())
	    SGnames.push_back(Name);
	}
      }
    }
    return SGnames;
  }
  //--------------------------------------------------------------
  scala::SysAbsScore PointGroup::TestPossible(const cctbx::sgtbx::space_group& sg,
					      std::vector<scala::Zone>& Zones) const
  {
    bool DEBUG = false;

    std::string Name = sg.type().lookup_symbol();
    int sgnumber = sg.type().number();
    //^
    if (DEBUG) std::cout << "\n\n***** SG: " << Name << "  " << sgnumber << "\n";
    // This spacegroup belongs to the right Laue group 
    // Loop systematic absence zones
    //?//    bool PossSG = true; // Start condition true, may then be falsified

    // Reindex operator from constructor (lattice) frame to
    // current reference frame
    scala::ReindexOp LtoS = SetReindexOp(ChBasis_ref);
    //  but we want the absence conditions in the space group
    //  standard frame
    sgtbx::change_of_basis_op SGref = sg.type().cb_op();
    scala::ReindexOp LtoSPGref = SetReindexOp(SGref * ChBasis_ref);

    if (DEBUG)
      {
	std::cout << "Reindex Lat to Laue group ref:  " << LtoS.as_hkl() << "\n";
	std::cout << "Reindex Lat to space group ref: " << LtoSPGref.as_hkl() << "\n";
      }

    std::string condition;
    std::string conditionLG;
    bool cond1 = true;
    bool cond1LG = true;
    scala::SysAbsScore Pyes;  // combined score

    for (size_t iz=0;iz<Zones.size();iz++)  {
      //^
      if (DEBUG) std::cout << "\nZone: " << Zones[iz].formatRefFrame();
      // Does this zone belong to this Laue group?
      if (RefLGname() == Zones[iz].LGsymm().symbol()) {
	// Yes
	// Does this zone have absences corresponding to
	// this spacegroup?
	if (!Zones[iz].Axis()) {
	  // This is a glide plane (implies !CHIRAL)
	  // (note all glides precede all axes)
	  if (GlidePresent(sg, Zones[iz], LtoS)) {
	    //^
	    if (DEBUG) 
	      std::cout << "   InSG, zone " << iz
			<< " Cond: "
			<< Zones[iz].FormatConditionNewFrame(LtoSPGref, 0)
			<< " Nobs " << Zones[iz].Nobs();
	    //^-
	    if (!cond1) condition += ", ";
	    if (!cond1LG) conditionLG += ", ";
	    cond1 = false;
	    cond1LG = false;
	    condition   += Zones[iz].FormatConditionNewFrame(LtoSPGref, 0);
	    conditionLG += Zones[iz].FormatConditionNewFrame(LtoS, 0);  // Laue group frame
	    Pyes.AddZone(iz+1);
	    // This glide is present in this spacegroup
	    // Only last status point is interesting
	    // ie 3rd (4n) point for d glide, else 2nd point (2n)
	    // Glide may be observed,  score positive
	    Pyes.ProbYes(Zones[iz].p().back());
	    //^
	    if (DEBUG) std::cout << "\nProbYes " << Zones[iz].p().back() 
				 << " Ptot " << Pyes.TotalProbability() << "\n";
	  } else {
	    //^
	    if (DEBUG) std::cout << "   NotInSG, zone " << iz << " Nobs " << Zones[iz].Nobs() << "\n";
	    // This glide is NOT present in this spacegroup
	    // Glide may be observed, score negative
	    Pyes.ProbNo(Zones[iz].p().back());
	    //^
	    if (DEBUG) std::cout << "ProbNo " << 1.0-Zones[iz].p().back()
				 << " Ptot " << Pyes.TotalProbability() << "\n";
	  }
	} else {
	  // This is an axis
	  // Does this spacegroup contain a screw along here?
	  int jscrew = ScrewPresent(sg, Zones[iz], LtoS);
	  
	  if (jscrew > 0) {
	    // This spacegroup contains the absence due to screw
	    //  (including screw obscured by glide)
	    //^
	    if (DEBUG) std::cout << "   InSG, zone " << iz << " " << jscrew << "  "
				 << Zones[iz].FormatConditionNewFrame(LtoSPGref, jscrew) << "  ";
	    //^-
	    if (!cond1)   condition += ", ";
	    if (!cond1LG) conditionLG += ", ";
	    cond1 = false;
	    condition += Zones[iz].FormatConditionNewFrame(LtoSPGref, jscrew);
	    cond1LG = false;
	    conditionLG += Zones[iz].FormatConditionNewFrame(LtoS, jscrew);
	    Pyes.AddZone(iz+1);
	    // Axis may be observed, score positive
	    Pyes.ProbYes(Zones[iz].p()[jscrew]);
	    //^
	    if (DEBUG) std::cout << "\nProbYes " << Zones[iz].p()[jscrew]
				 << " Ptot " << Pyes.TotalProbability()
				 << " Nobs " << Zones[iz].Nobs() << "\n";
	  } else {
	    // This spacegroup contains no screw
	    //^
	    if (DEBUG) std::cout << "   NotInSG, zone " << iz <<  "\n";
	    // Axis may be observed, multiply in score for no screw
	    // this is the 0'th element of the vector p
	    Pyes.ProbYes(Zones[iz].p()[0]);
	    //^
	    if (DEBUG) std::cout << "\nProbNo " << Zones[iz].p()[0]
				 << " Ptot " << Pyes.TotalProbability()
				 << " Nobs " << Zones[iz].Nobs() << "\n";
	  }
	}
      }
    }  // end zone loop
    Pyes.SetCondition(condition);
    Pyes.SetConditionLG(conditionLG);
        
    Chirality ch;
    if (sg.is_chiral()) ch = scala::CHIRAL;
    else if (sg.is_centric()) ch = scala::CENTROSYMMETRIC;
    else ch = scala::NONCHIRAL;

    Pyes.SetSGposs(scala::PossibleSpaceGroup
		    (DeColonName(Name),
		     sgnumber,
		     SpaceGroupName(sg.type(), 'H'),           // this should be reference name
		     Pyes.Condition(), Pyes.ConditionLG(),
		     //  Change of basis to reference frame
		     SetReindexOp(sg.type().cb_op()),
		     ch,
		     Pyes.TotalProbability(),
		     Pyes.ZoneList()));

    return Pyes;
  }
  //--------------------------------------------------------------
  std::vector<scala::PossibleSpaceGroup> PointGroup::TestSpaceGroupList
    (std::vector<scala::Zone>& Zones,
     const Chirality& chiral) const
  // Return list of all spacegroups compatible with Laue group (actual frame)
  // and scored zones
  // If chiral = CHIRAL, list only chiral spacegroups
  //
  // Zones are in some arbitrary subgroup frame, but with stored reindex
  // operator from reference "lattice" frame
  //
  // Each space group needs to be compared to each zone in the current class
  // (Laue group) reference frame
{
    //^
    bool DEBUG = false;

    char CLtype = CentringSymbol(LaueGrp_ref);
    std::vector<scala::PossibleSpaceGroup> SGposs;
    //^
    if (DEBUG) {
      std::cout << "SGposs.size " << SGposs.size() << "\n";
    }

    cctbx::sgtbx::space_group LG = LaueGrp_ref;
    cctbx::sgtbx::space_group LSG;
    scala::PossibleSpaceGroup SGp;

    cctbx::sgtbx::space_group_symbol_iterator sgiter;
    for(;;) {
      cctbx::sgtbx::space_group_symbols symbol = sgiter.next();
      if (symbol.number() == 0) break;
      cctbx::sgtbx::space_group sg(symbol.hall());
      bool valid =
	(chiral == scala::CHIRAL && sg.is_chiral()) ||
	(chiral == scala::CENTROSYMMETRIC && sg.is_centric()) ||
	(chiral == scala::NONCHIRAL);
      
      if (valid) {
	std::string Name = sg.type().lookup_symbol();
	int sgnumber = sg.type().number();
	scala::SysAbsScore Pcomb;
	char LT = CentringSymbol(sg);
	LSG = sg.build_derived_patterson_group();
	if (LT == CLtype &&
	    LSG == LG) {
	  Pcomb = TestPossible(sg, Zones);
	  if (Pcomb.IsPossible()) {
	    SGp = Pcomb.SGposs();
	    // This is a possible spacegroup, add to list
	    if (std::find(SGposs.begin(), SGposs.end(), SGp) == SGposs.end()) {
	      SGposs.push_back(SGp);
	      if (DEBUG) std::cout << "\n>> Spacegroup added: "
				   << SGposs.size() << " " << SGp.Name() << "\n";		      
	    }
	    //^+
	  } else {
	    if (DEBUG) std::cout << "\n>> Spacegroup not possible: "
				 << Name << "  " << sgnumber << "\n";
	  }
	  //^-
	  // SPECIAL for space group 205 p 21/a -3
	  // Try permuted version
	  if (sgnumber == 205) {
	    sgtbx::change_of_basis_op cb("l,k,-h");
	    sg = sg.change_basis(cb);
	    Pcomb = TestPossible(sg, Zones);
	    if (Pcomb.IsPossible()) {
	      scala::ReindexOp rindx = scala::ReindexOp("-l,k,h");
	      scala::PossibleSpaceGroup SGp2
		(Pcomb.SGposs().Name(true), SGp.Sgnumber(), SGp.Name(false),
		 SGp.Condition(), SGp.ConditionLG(),
		 rindx, SGp.Chiral(),
		 Pcomb.SGposs().Prob(), SGp.ZoneList());
	      SGposs.push_back(SGp2);
	    }
	  }  // end SPECIAL
	}  // end possible Laue group
      }
    }  // end sg loop
    if (DEBUG) {
      std::cout << "SGpass.size (2) " << SGposs.size() << "\n";
    }
    return SGposs;
  }
  //--------------------------------------------------------------
  bool PointGroup::Equals(const PointGroup& other) const
  // Equality of operators in original (cell) frame
  {
    return (RotGrp_ref == other.RotGrp_ref) &&
      (ChBasis.c() == other.ChBasis.c());
  }
  //--------------------------------------------------------------
  bool PointGroup::EqualsRef(const PointGroup& other) const
  // Equality of operators in reference frame
  {
    return (RotGrp_ref == other.RotGrp_ref);
  }
  //--------------------------------------------------------------
  void PointGroup::PrintAllSpacegroups(phaser_io::Output& output,
				       const Chirality chiral,
				       const bool AllLattice) const
  {
    std::vector<std::string> SGnames = SpaceGroupList(chiral, AllLattice);
    if (SGnames.size() > 0) {
      output.logTab(0,LOGFILE,"Possible spacegroups:\n");
      int nc = 0;
      for (size_t i=0;i<SGnames.size();i++) {
	if (nc > 80) {
	  output.logTab(0,LOGFILE,"\n");
	  nc = 0;
	}
	output.logTab(0,LOGFILE," <"+SGnames[i]+">");
	nc += SGnames[i].size() + 3;
      }
      output.logTab(0,LOGFILE,"\n");
    } else {
      output.logTab(0,LOGFILE,"No compatible spacegroups\n");
    }
  }
  //--------------------------------------------------------------
  std::vector<scala::ReindexOp>
  AlternativeIndexing(const PointGroup& PG,
		      const bool& strict,
		      const scala::Scell target_cell,
		      const float& max_delta,
		      const int& AllowI2)
    // Return list of possible alternative indexing schemes
    // compatible with cell or with lattice symmetry
    //
    // (1) If strict = true
    //   return symmetry operators from symmetry elements 
    //   present in lattice group but not in pointgroup
    //   (max_delta & target_cell are ignored)
    //   This can only happen for symmetries > orthorhombic
    // else (2) strict = false
    //   find accidental alternatives arising from special
    //   cell dimension relationships
    //   max_delta is tolerance for cell similarity
    //    to target_cell
    //
    //  Returns:-
    //   List of reindex operators, including:-
    //     strict flag: true if "strict"
    //     deviations - cell deviations from target cell
    //                  (maximum angle in degrees)
    //                  = 0.0 for strict settings
    //
    // Friend of class PointGroup
  {
    std::vector<scala::ReindexOp> ReindexList;

    if (strict) {
      if (LatticeGroup(PG.LaueGrp_ref).crystal_system() <= ORTHORHOMBIC)
	// strict option only valid for high symmetry
	{return ReindexList;}

      bool ExcludeIdentity = false, AnyCell = false, BestCell = false;
      std::vector<AlternativeBases> CbOp_list = 
	GetAlternativeBases(PG.LaueGrp_ref, ExcludeIdentity, AnyCell, BestCell,
			    PG.uccell_ref, PG.uccell_ref, 3.0, max_delta, AllowI2);
      // change-of-basis operators are relative to reference cell
      if (CbOp_list.size() > 0) {
	// Change basis of all operators to cell frame
	std::vector<AlternativeBases> CBlist =
	  ChangeBasesList(CbOp_list, PG.ChBasis);
	
	for (size_t k=0;k<CBlist.size();k++) {
	  ReindexList.push_back(SetReindexOp(CBlist[k].SimplestOp()));
	  ReindexList.back().SetStrict(true);
	}
      }
    } else {
      // not strict
      // Lattice group of pointgroup
      LatticeGroup LG(PG.LaueGrp_ref);
      // List of possible "strict" bases (high symmetry only)
      std::vector<scala::ReindexOp> StrictOps;
      if (LG.crystal_system() > ORTHORHOMBIC) {
	std::vector<sgtbx::change_of_basis_op> StrictCB
	  = PossibleChBOp(LG, PG.LaueGrp_ref, AllowI2);
	for (size_t i=0;i< StrictCB.size();i++)
	  {StrictOps.push_back(SetReindexOp(StrictCB[i]));}
      }
      // Maximum lattice symmetry 
      CCtbxSym::LatticeSymmetry lat(PG.input_cell, PG.LatType, AllowI2, max_delta);
      // Reindexing operator for "best" spacegroup
      //  from original -> lattice (best), inverse operator for hkl
      scala::ReindexOp reindex_op = MVutil::SetCMat33(lat.best_sg_reindex_op());
      // All subgroups
      std::vector<PointGroup> subgroups =
	scala::GetSubGroups(scala::hkl_symmetry(lat.best_spacegroup_symbol()));
      if (subgroups.size() > 0) {
	// Loop subgroups to find any which are the same as this
	for (size_t k=0;k<subgroups.size();k++) {
	  if (subgroups[k].LaueGrp_ref == PG.LaueGrp_ref) {
	    // Yes it is
	    double celldiff = subgroups[k].SetCell(PG.input_cell,reindex_op, AllowI2);
	    celldiff = celldiff;
	    //^ debug
	    /*		    std::cout << "\nAlternative found: reindex "
	      << subgroups[k].RefSGreindexFormat() << "\n";
	      std::cout << "Cell: ";
	      std::vector<double> tcell =
	      subgroups[k].TransformedCell();
	      for (int i=0;i<6;i++) std:: cout << " " << tcell[i];
	      std::cout << "\nInput Cell: ";
	      for (int i=0;i<6;i++) std:: cout << " " << PG.input_cell[i];
	      std::cout << "\n";*/
	    //^
	    
	    std::vector<double> cell_diffs;

	    std::vector<scala::ReindexOp> possible_reindex =
	      subgroups[k].GetCloseCell(target_cell,
					max_delta, false,
					cell_diffs, AllowI2);
	    
	    for (size_t i=0;i<possible_reindex.size();i++) {
	      // Is this a "strict" operator?
	      possible_reindex[i].
		SetStrict(IsInList<scala::ReindexOp>
			  (StrictOps,possible_reindex[i]));
	      possible_reindex[i].SetDeviation(cell_diffs[i]);
	      ReindexList.push_back(possible_reindex[i]);
	    }
	  }
	} 
      }
    }  
    std::sort(ReindexList.begin(),ReindexList.end());

    return ReindexList;
  }
}
