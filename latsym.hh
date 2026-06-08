// latsym.hh

#ifndef SYM_LATSYM
#define SYM_LATSYM

#include <cctbx/crystal/symmetry.h>
#include <cctbx/sgtbx/lattice_symmetry.h>

namespace CCtbxSym
{
  //--------------------------------------------------------------
  // Format spacegroup name in "standard" convention
  // Second argument is symbol to use for rhombohedral lattice (R or H)
  // if ReferenceSetting true (default), then use the reference setting
  // rather than the actual setting
  std::string SpaceGroupName(const cctbx::sgtbx::space_group_type& SGtype,
			     const char Rlattice = 'R', const bool& ReferenceSetting=true);
  //--------------------------------------------------------------
  // Format spacegroup name in "standard" convention
  // Second argument is symbol to use for rhombohedral lattice (R or H)
  std::string SpaceGroupName(const cctbx::sgtbx::space_group& SG,
			     const char Rlattice = 'R', const bool& ReferenceSetting=true);
  //--------------------------------------------------------------
  // Format spacegroup name in "standard" convention
  // Second argument is symbol to use for rhombohedral lattice (R or H)
  std::string SpaceGroupName(const std::string& SGname,
			     const char Rlattice = 'R', const bool& ReferenceSetting=true);
  //--------------------------------------------------------------
  std::string SpaceGroupName(const int& SpGpNumber,
			     const char Rlattice = 'R', const bool& ReferenceSetting=true);
  //--------------------------------------------------------------
  // return extended Hermann-Mauguin symbol
  // for rhombohedral space groups,set lattice type to
  // Rlattice = 'R' or 'H'
  std::string SpaceGroup_xHM(const int& SpGpNumber,
			     const char Rlattice = 'H');
  //--------------------------------------------------------------
  bool RhombohedralAxes(const std::vector<double>& unit_cell_dimensions);
  //--------------------------------------------------------------
  void
  show_space_group_type(cctbx::sgtbx::space_group_type const& space_group_type);
  //--------------------------------------------------------------
  void
  show_unit_cell(cctbx::uctbx::unit_cell const& unit_cell);
  //--------------------------------------------------------------
  cctbx::sgtbx::change_of_basis_op ReduceCell(const cctbx::uctbx::unit_cell& cell,
  		 const cctbx::sgtbx::change_of_basis_op RefOp = cctbx::sgtbx::change_of_basis_op());
  //--------------------------------------------------------------
  class LatticeSymmetry{
  public:
    // Constructor from vector of cell dimensions (6 numbers)
    LatticeSymmetry(const std::vector<double>& unit_cell_dimensions,
		    const char lattice_type,
		    const int& AllowI2,
		    const double max_delta = 3.0);
    void print() const;

    // Return spacegroup symbol (Hall symbol) for "best" spacegroup
    std::string best_spacegroup_symbol() const;
    // return change-of-basis operator from input to "best"
    // cell as 3x3 matrix (in vector in order 11, 12, 13, etc)
    // This operator C can be used to reindex old hkl as
    //   h' = h C   (row vectors h, h')
    std::vector<double> best_sg_reindex_op() const;
    // Return "best" unit cell
    std::vector<double> best_unit_cell() const;
    // Return original unit cell with changed basis
    std::vector<double> orig_unit_cell() const;
    // Return maximum angular devation from constructor cell
    double Delta() const {return delta;}

  private:
    // symmetry contains unit_cell and spacegroup
    cctbx::crystal::symmetry input_symmetry_;
    cctbx::crystal::symmetry lattice_symmetry_;
    cctbx::crystal::symmetry best_symmetry_;
    // operator to change input setting to "best" setting
    cctbx::sgtbx::change_of_basis_op cb_op_inp_best;

    char lattice_type_;
    double max_delta_;
    double delta;
  };

}

#endif
