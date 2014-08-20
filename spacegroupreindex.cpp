// spacegroupreindex.cpp

#include <cctbx/crystal/symmetry.h>
#include <cctbx/sgtbx/lattice_symmetry.h>
#include <cctbx/uctbx/fast_minimum_reduction.h>

#include "spacegroupreindex.hh"
#include "pointgroup.hh"
#include "latsym.hh"

using namespace cctbx;
using clipper::Message_warn;

namespace scala {
  //--------------------------------------------------------------
  ReindexOp SpacegroupReindexOp(const std::string& from_SGname,
                                const std::string& to_SGname)
  // return reindex operator to convert alternative space groups
  // in same crystal system
  // Mainly (only?) for C2<->I2 and H3<->R3
  // Throws a clipper::Message_warn exception if the space groups are
  // not just alternative settings (ie have different reference groups)
  {
    // from space group
    std::string frname = CCtbxSym::CCTBX_SGsymbol_HorR(from_SGname);
    sgtbx::space_group from_SG =
      sgtbx::space_group(sgtbx::space_group_symbols(frname).hall());
    // Change of basis to reference setting
    sgtbx::change_of_basis_op ChB_ref_from = from_SG.type().cb_op();
    sgtbx::space_group from_SG_ref = from_SG.change_basis(ChB_ref_from);
    //^    CCtbxSym::PrintChBOp(ChB_ref_from); //^

    // to space group
    sgtbx::space_group to_SG =
      sgtbx::space_group(sgtbx::space_group_symbols
                         (CCtbxSym::CCTBX_SGsymbol_HorR(to_SGname)).hall());
    // Change of basis to reference setting
    sgtbx::change_of_basis_op ChB_ref_to = to_SG.type().cb_op();
    sgtbx::space_group to_SG_ref = to_SG.change_basis(ChB_ref_to);
    //^    CCtbxSym::PrintChBOp(ChB_ref_to); //^

    // Reference groups should be the same
    if (from_SG_ref != to_SG_ref) {
      std::string message =
        CCtbxSym::SpaceGroupName(from_SG_ref.type(), 'H')+
        " has different reference setting from "+
        CCtbxSym::SpaceGroupName(to_SG_ref.type(), 'H')+"\n";
      Message::message(Message_warn(message));
      throw Message_warn(message);
    }

    // We want the transformation from HKLIN to input
    ReindexOp reindex = CCtbxSym::SetReindexOp(ChB_ref_to.inverse() * ChB_ref_from);
    //^
    //    ReindexOp reindex_from = CCtbxSym::SetReindexOp(ChB_ref_from);
    //    ReindexOp reindex_to = CCtbxSym::SetReindexOp(ChB_ref_to);
    //    std::cout << "[H]from " << reindex_from.as_hkl() << "\n";
    //    std::cout << "[H]to   " << reindex_to.as_hkl() << "\n";
    //^

    return reindex;
  }
  //--------------------------------------------------------------
  bool SpacegroupReindex(const GlobalControls& GC,
                         const hkl_symmetry& HKLINsymm, const Scell& cell,
                         ReindexOp& Reindex, phaser_io::Output& output)
  // If SPACEGROUP is specified but no REINDEX operator, generate appropriate reindexing
  // to convert from input HKLIN file HKLINsymm to desired spacegroup
  // Probably really only useful (or indeed valid) for C2 <-> I2 & H3<->R3, or P222 groups
  //
  // input cell corresponds to HKLINsymm
  //
  // Returns true if Reindex is set
  // fails if the symmetries do not belong to same lattice group
{
  if (GC.Spacegroup() == "" ||  GC.Spacegroup() == "HKLIN" || GC.IsReindexSet()) return false;

  std::string HKLIN_SGname = HKLINsymm.symbol_xHM();
  hkl_symmetry NewSymm(GC.Spacegroup());
  std::string Input_SGname = NewSymm.symbol_xHM();
  if (NewSymm.CrysSys() != HKLINsymm.CrysSys()) {
      std::string message =
        "Specified SPACEGROUP "+GC.Spacegroup()+
        " must belong to same crystal system as input symmetry "+HKLIN_SGname+
        " unless REINDEX is explicitly given.\n";
      Message::message(Message_fatal(message));
  }

  try {
    // Get reindex operator
    Reindex = SpacegroupReindexOp(HKLIN_SGname, Input_SGname);
  }
  catch (Message_warn& warn) {
    std::string message =
      "Specified SPACEGROUP "+GC.Spacegroup()+
      " must have the same 'reference' space group as the input file symmetry "+HKLIN_SGname+
        " unless REINDEX is explicitly given.\n";
      Message::message(Message_fatal(message));
  }

  int AllowI2 =  GC.AllowI2();
  if (NewSymm.lattice_type() != 'I') {
    AllowI2 = 0;  // don't allow I lattice if we've asked for C2
  }

  if (NewSymm.lattice_type() != 'R') {  // Don't reduce to reference group for R lattice
    // Get reindex operator from required group Input_SGname to "reference" seting
    CCtbxSym::PointGroup PG(Input_SGname);
    Scell scell = cell.change_basis(Reindex);
    PG.SetCell(scell.UnitCell(), ReindexOp(), AllowI2);
    ReindexOp newreindex = PG.RefSGreindex(); // reindex cell -> best
    //^     std::cout << "Reindex " << Reindex.as_hkl() <<"\n"; //^
    //^     std::cout << "newreindex " << newreindex.as_hkl() <<"\n"; //^
    Reindex = Reindex * newreindex;
  }

  output.logTab(0,LOGFILE,
                "\nReindexing data with operator "+Reindex.as_hkl()+
                " from space group "+HKLIN_SGname+" to "+Input_SGname);

  return true;
}
}
