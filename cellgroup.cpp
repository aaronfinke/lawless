// cellgroup.cpp

#include "cellgroup.hh"
#include "hkl_symmetry.hh"
#include "cctbx_utils.hh"

namespace CCtbxSym
{
  //--------------------------------------------------------------
  CellGroup::CellGroup(const scala::SpaceGroup& spgp)
  // Construct from SpaceGroup object
  {
    //  initialise space-group
    group.reset();
    
    int Nsymp = spgp.num_primops();
    for (int k=0;k<Nsymp;k++) {
      // Make rt_mx matrices
      sgtbx::rt_mx R1 = MVutil::SetRtMx(MVutil::SetVMat33(spgp.symop(k).rot()));
      group.expand_smx(R1);
    }
    group.make_tidy();
  }
  //--------------------------------------------------------------
  scala::Scell CellGroup::constrain(const scala::Scell& cell) const
  // Return cell constrained by symmetry
  {
    scitbx::af::double6 dcell;
    for (int i=0;i<6;i++) dcell[i] = cell[i];
    uctbx::unit_cell uccell(dcell);
    uctbx::unit_cell uccell_constrained = group.average_unit_cell(uccell);
    std::vector<double> dbcell(6);
    for (int i=0;i<6;i++) {dbcell[i] = uccell_constrained.parameters()[i];}
    return scala::Scell(dbcell);
  }
}
