// latsym.cpp

#include "latsym.hh"
#include "pointgroup.hh"
#include <cstdio>

using namespace cctbx;

namespace CCtbxSym{
  //--------------------------------------------------------------
  bool centred_monoclinic(const int& sgn)
  // true if space group number sgn is C/I-centred monoclinic  
  {
    if (sgn == 5 || sgn == 8 || sgn == 9 || sgn == 12 || sgn == 15) return true;
    return false;
  }
  //--------------------------------------------------------------
  // Format spacegroup name in "standard" convention
  // but allow I2 setting of centred monoclinic
  // if ReferenceSetting true (default), then use the reference setting
  // rather than the actual setting
  std::string SpaceGroupName(const sgtbx::space_group_type& SGtype,
			     const char Rlattice, const bool& ReferenceSetting)
  {
    std::string symbol = SGtype.lookup_symbol();
    sgtbx::space_group_symbols sgs(symbol);
    std::string sgname = sgs.hermann_mauguin();
    //    std::cout << "SpaceGroupName " << sgname << "\n"; //^
    int sgn = SGtype.number();
    // if ReferenceSetting true, then use the reference setting rather than to actual
    // setting
    if (ReferenceSetting) {
      sgname = sgtbx::space_group_symbols(sgn).hermann_mauguin();
    }
    // Check for monoclinic centred lattice
    if (centred_monoclinic(sgn)) {
      sgname = sgtbx::space_group_symbols(SGtype.lookup_symbol()).hermann_mauguin();
    }
    // Look for "R" or "H"
    std::string::size_type i = sgname.find("R");
    if (i >= sgname.size())   // "R" not found
      i = sgname.find("H");      // look for H
    if (i < sgname.size())    // "R" or "H" found
      sgname[i] = Rlattice;
    return sgname;
  }
  //--------------------------------------------------------------
  // Format spacegroup name in "standard" convention
  std::string SpaceGroupName(const sgtbx::space_group& SG,
			     const char Rlattice, const bool& ReferenceSetting)
  {
    return SpaceGroupName(SG.type(), Rlattice, ReferenceSetting);
  }
  //--------------------------------------------------------------
  std::string SpaceGroupName(const std::string& SGname,
			     const char Rlattice, const bool& ReferenceSetting)
  {
    sgtbx::space_group SG = sgtbx::space_group
      (sgtbx::space_group_symbols(SGname).hall());
    return SpaceGroupName(SG.type(), Rlattice, ReferenceSetting);
  }
  //--------------------------------------------------------------
  std::string SpaceGroupName(const int& SpGpNumber,
			     const char Rlattice, const bool& ReferenceSetting)
  {
    std::string extension = "";
    if (Rlattice == 'R' || Rlattice == 'H') {extension = Rlattice;}
    sgtbx::space_group SG(sgtbx::space_group_symbols(SpGpNumber, extension));
    return SpaceGroupName(SG.type(), Rlattice, ReferenceSetting);
  }
  //--------------------------------------------------------------
  std::string SpaceGroup_xHM(const int& SpGpNumber,
			      const char Rlattice)
  // return extended Hermann-Mauguin symbol
  // for rhombohedral space groups,set lattice type to
  // Rlattice = 'R' or 'H'
  {
    sgtbx::space_group_symbols SGS(SpGpNumber);
    char ext = SGS.extension();
    std::string extension = "";
    if (ext == 'H')
      // extension = 'H' for rhombohedral space groups
      // reset to  'R' or 'H' from input
      {if (Rlattice == 'R' || Rlattice == 'H') {extension = Rlattice;}}
    SGS = sgtbx::space_group_symbols(SpGpNumber, extension);
    //    std::cout << "SGS: HM" << SGS.hermann_mauguin()
    //	      << " extension :" << SGS.extension()
    //	      << ": Hall " << SGS.hall()
    //	      << " EHM " << SGS.extended_hermann_mauguin()
    //	      << "\n";
    //**  call changed in cctbx version 2
   #if defined (CCTBX_VERSION) && (CCTBX_VERSION <= 2006)
    return  SGS.extended_hermann_mauguin();
   #else
    return  SGS.universal_hermann_mauguin();
   #endif
  }
  //--------------------------------------------------------------
  void
  show_space_group_type(const sgtbx::space_group_type& SG_type)
  {
    std::printf("%s (No. %d)", SpaceGroupName(SG_type).c_str(),
		SG_type.number());
  }

  //--------------------------------------------------------------
  void
  show_unit_cell(uctbx::unit_cell const& unit_cell)
  {
    for(std::size_t i=0;i<6;i++) {
      std::printf("%s%.6g%s",
		  (i < 1 ? "(" : " "),
		  unit_cell.parameters()[i],
		  (i < 5 ? "," : ")"));
    }
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  void SymPrint(const crystal::symmetry& sym)
  {
    std::printf("\n");
    std::printf("Unit cell: ");
    show_unit_cell(sym.unit_cell());
    std::printf("\n");
    std::printf("Space group: ");
    show_space_group_type(sgtbx::space_group_type(sym.space_group()));
    std::printf("\n");
    std::printf("\n");
  }
  //--------------------------------------------------------------
  void SGPrint(const sgtbx::space_group& spgrp)
  {
    std::printf("\n");
    std::printf("Space group: ");
    show_space_group_type(sgtbx::space_group_type(spgrp));
    std::printf("\n");
  }
  //--------------------------------------------------------------
  void PrintChBasisOp(const sgtbx::change_of_basis_op ch_op,
		      const std::string label)
  {
    std::printf("\nChange of basis operator for %s",label.c_str());
    std::printf(":  %s\n",ch_op.c().as_xyz().c_str());
    std::printf("Inverse change of basis operator for %s",label.c_str());
    std::printf(":  %s\n",ch_op.c_inv().as_xyz().c_str());
  }
  //--------------------------------------------------------------
  sgtbx::space_group ChangeBasis(const sgtbx::space_group& SG,
					const sgtbx::change_of_basis_op& CbOp)
  // Create new spacegroup with changed basis from
  // primitive operators only
  {
//^!    std::cout << "\nChangeBasis:";
//^!    PrintChBasisOp(CbOp,"");
//^!    std::cout << "  starting symmetry\n";
//^!    PrintCctbxSymops(SG);
    cctbx::sgtbx::space_group new_sg(false, cctbx::sgtbx::sg_t_den);;
    for (std::size_t i=1;i<SG.n_smx();i++) {
      cctbx::sgtbx::rt_mx sgsmx = SG.smx(i).new_denominators(CbOp.c());
      cctbx::sgtbx::rt_mx cinv = CbOp.c_inv();
      cctbx::sgtbx::rt_mx c = CbOp.c();
      cctbx::sgtbx::rt_mx c1 = sgsmx.multiply(cinv);
      cctbx::sgtbx::rt_mx c2 = c.multiply(c1);
      cctbx::sgtbx::rt_mx c3 = c2.new_denominators(new_sg.smx(0));
      //^      cctbx::sgtbx::rt_mx cbsmx = CbOp(sgsmx);
      //      cctbx::sgtbx::rt_mx cbsmx = CbOp.apply(sgsmx);
      new_sg.expand_smx(c3);
      //      new_sg.expand_smx(CbOp.apply(SG.smx(i)));
      //^      new_sg.expand_smx(CbOp(SG.smx(i)));
    }
    new_sg.make_tidy();
//^!    std::cout << "  changed symmetry\n";
//^!    PrintCctbxSymops(new_sg);
//^!    std::cout << " >> END ChangeBasis\n";
    return new_sg;
  }
  //--------------------------------------------------------------
  bool RhombohedralAxes(const std::vector<double>& unit_cell_dimensions)
    // True if not hexagonal setting (angles 90,90,120)
  {
    double tol =  0.2;  // tolerance on angles
    if (Close<double,double>(unit_cell_dimensions.at(3), 90.0, tol))
      if (Close<double,double>(unit_cell_dimensions.at(4), 90.0, tol))
	if (Close<double,double>(unit_cell_dimensions.at(5), 120.0, tol))
	  return false;
    return true;
  }
  //--------------------------------------------------------------
  sgtbx::change_of_basis_op ReduceCell(const uctbx::unit_cell& cell,
				       const sgtbx::change_of_basis_op RefOp)
  {
    uctbx::fast_minimum_reduction<> red(cell);
    return sgtbx::change_of_basis_op
      (sgtbx::rt_mx(sgtbx::rot_mx(red.r_inv(), 1)).inverse())
      .new_denominators(RefOp);
  }
  //--------------------------------------------------------------
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--------------------------------------------------------------
  // LatticeSymmetry class

  LatticeSymmetry::LatticeSymmetry(const std::vector<double>& unit_cell_dimensions,
				   const char lattice_type,
				   const int& AllowI2,
				   const double max_delta)
    : lattice_type_(lattice_type), max_delta_(max_delta)
  {
    // Symmetry object combines unit_cell object and
    // space-group (in this P1, C1, I1, F1 or R1)
    scitbx::af::double6 dcell;
    for (int i=0;i<6;i++) dcell[i] = unit_cell_dimensions[i];
    // CCtbx routines need lattice type H as R
    if (lattice_type_ == 'H') lattice_type_ = 'R';
    // but if the rhombohedral setting is really rhombohedral,
    // then it needs to be treated as P here, otherwise the z2p_op
    // is wrong
    if (lattice_type_ == 'R')
      {
	// Hexagonal setting has angles 90,90,120
	// Rhombohedral has alpha=beta=gamma
	if(RhombohedralAxes(unit_cell_dimensions))
	  lattice_type_ = 'P';
      }

    input_symmetry_ = crystal::symmetry(uctbx::unit_cell(dcell),
					sgtbx::space_group(std::string(1,lattice_type_)+" 1"));

    // get an operator to convert centred cell to primitive
    //  == I if P already
    sgtbx::change_of_basis_op z2p_op = input_symmetry_.space_group().z2p_op();

    // Reduce primitive cell
    uctbx::unit_cell p_cell = input_symmetry_.unit_cell().change_basis
      (z2p_op.c_inv().r().as_double());
    sgtbx::change_of_basis_op red_op = ReduceCell(p_cell, z2p_op);

    z2p_op = red_op * z2p_op;
    crystal::symmetry primitive_symmetry = input_symmetry_.change_basis(z2p_op);

    // Get highest symmetry compatible with lattice
    // (note the group_search overloads operator()
    //   to create spacegroup)
    // computes potential axes
    //**  call changed in cctbx version 2
    #if defined (CCTBX_VERSION) && (CCTBX_VERSION <= 2006)
        //** old version 1
	sgtbx::lattice_symmetry::group_search lattice_symmetry_group;
        cctbx::sgtbx::space_group lattice_group =
	  lattice_symmetry_group(primitive_symmetry.unit_cell(), max_delta);
    #else
 	//** for new cctbx **
	cctbx::sgtbx::space_group lattice_group =
	  cctbx::sgtbx::lattice_symmetry::group(primitive_symmetry.unit_cell(), max_delta);
    #endif

    lattice_group.make_tidy();
    
    lattice_symmetry_ = crystal::symmetry(primitive_symmetry.unit_cell(),
					  lattice_group);
     
    // Adjust unit cell to fit lattice symmetry
    //   (eg force angles = 90 etc)
    crystal::symmetry adjust_sym
      (lattice_group.average_unit_cell(primitive_symmetry.unit_cell()),
       lattice_group);

    // Convert adjusted cell to reference setting
    sgtbx::space_group_type adjust_sym_type(adjust_sym.space_group());
    sgtbx::change_of_basis_op cb_op_ref = adjust_sym_type.cb_op();
    best_symmetry_ = adjust_sym.change_basis(cb_op_ref);

    //^
    //    std::cout << "BestSymmetry " << best_symmetry_.space_group().type().hall_symbol() << "\n";
    //    std::cout << "Cell: " <<  UcellFormat(best_symmetry_.unit_cell()) << "\n";
    //^-
    //  Select "best" orthorhombic or monoclinic cell
    sgtbx::change_of_basis_op cb_op_opt;
    LatticeGroup LatG(best_symmetry_.space_group());
    if (LatG.crystal_system() == ORTHORHOMBIC || LatG.crystal_system() == MONOCLINIC)
      {
	cb_op_opt =
	  GetBestCell(best_symmetry_.space_group(), best_symmetry_.unit_cell(), AllowI2);
	best_symmetry_ = best_symmetry_.change_basis(cb_op_opt);
      }

    // Total basis transformation
    cb_op_inp_best =  cb_op_opt * cb_op_ref * z2p_op;

    // Use identity change-of-basis operator if possible
    if (best_symmetry_.unit_cell().is_similar_to(input_symmetry_.unit_cell())) {
      sgtbx::change_of_basis_op cb_op_corr = cb_op_inp_best.inverse();
      if (   best_symmetry_.change_basis(cb_op_corr).space_group()
             == best_symmetry_.space_group()) {
	cb_op_inp_best = cb_op_corr * cb_op_inp_best;
      }
    }
    // Maximum deviation from original cell
    delta = sgtbx::lattice_symmetry::
      find_max_delta(cb_op_inp_best.apply(input_symmetry_.unit_cell()),
		     best_symmetry_.space_group().build_derived_point_group());
  }

  //--------------------------------------------------------------
  void LatticeSymmetry::print() const
  { 
    std::printf("\n\n***********\n\n");
    std::printf("Maximum (lattice) Symmetry : ");
    SymPrint(lattice_symmetry_);
    std::printf("Best (reference) Symmetry : ");
    SymPrint(best_symmetry_);
    std::printf("          Reindex operator (input->best): %s\n",
		cb_op_inp_best.c_inv().as_xyz().c_str());
    std::printf("\n***********\n\n");
  }

  //--------------------------------------------------------------
  std::string LatticeSymmetry::best_spacegroup_symbol() const
  {
    return SpaceGroupName(best_symmetry_.space_group(), 'H');
  }
  //--------------------------------------------------------------
  std::vector<double> LatticeSymmetry::best_sg_reindex_op() const
  {
    scitbx::mat3<double> op = cb_op_inp_best.c_inv().r().as_double();
    std::vector<double> vop(9);
    int k=0;
    for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
	vop[k++] = op(i,j);
    
    return vop;
  }
  //--------------------------------------------------------------
  std::vector<double> LatticeSymmetry::best_unit_cell() const
  {
    std::vector<double> cell(6);
    for (int i=0;i<6;i++) cell[i] = best_symmetry_.unit_cell().parameters()[i];
    return cell;
  }
  //--------------------------------------------------------------
  std::vector<double> LatticeSymmetry::orig_unit_cell() const
    // Return original unit cell with changed basis
  {
    cctbx::crystal::symmetry orig_sym =
      input_symmetry_.change_basis(cb_op_inp_best);
    std::vector<double> cell(6);
    for (int i=0;i<6;i++) cell[i] = orig_sym.unit_cell().parameters()[i];
    return cell;
  }
} // namespace CCtbxSym
