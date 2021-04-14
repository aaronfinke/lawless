// spacegroupreindex.cpp

#include <cctbx/crystal/symmetry.h>
#include <cctbx/sgtbx/lattice_symmetry.h>
#include <cctbx/uctbx/fast_minimum_reduction.h>

#include "spacegroupreindex.hh"
#include "pointgroup.hh"
#include "latsym.hh"
#include "string_util.hh"
#include "report_errors.hh"

using phaser_io::LOGFILE;
using phaser_io::LXML;

using namespace cctbx;

namespace scala {
  //--------------------------------------------------------------
  SpacegroupReindexOp::SpacegroupReindexOp(const std::string& from_SGname,
                                           const std::string& to_SGname)
  {
    sameintensitygroup = false;
    samereferencegroup = false;
    samepointgroup = false;
    reindex = ReindexOp(); // h,k,l

    // from space group
    std::string frname = CCtbxSym::CCTBX_SGsymbol_HorR(from_SGname);
    sgtbx::space_group from_SG =
      sgtbx::space_group(sgtbx::space_group_symbols(frname).hall());
    // Change of basis to reference setting
    sgtbx::change_of_basis_op ChB_ref_from = from_SG.type().cb_op();
    sgtbx::space_group from_SG_ref = from_SG.change_basis(ChB_ref_from);
    //CCtbxSym::PrintChBOp(ChB_ref_from); //^

    // to space group
    sgtbx::space_group to_SG =
      sgtbx::space_group(sgtbx::space_group_symbols
                         (CCtbxSym::CCTBX_SGsymbol_HorR(to_SGname)).hall());
    // Change of basis to reference setting
    sgtbx::change_of_basis_op ChB_ref_to = to_SG.type().cb_op();
    sgtbx::space_group to_SG_ref = to_SG.change_basis(ChB_ref_to);
    //CCtbxSym::PrintChBOp(ChB_ref_to); //^
    //^
    //    ReindexOp reindex_from = CCtbxSym::SetReindexOp(ChB_ref_from);
    //    ReindexOp reindex_to = CCtbxSym::SetReindexOp(ChB_ref_to);
    //    std::cout << "[H]from " << reindex_from.as_hkl() << "\n";
    //    std::cout << "[H]to   " << reindex_to.as_hkl() << "\n";
    //^

    //^^
    //    std::cout << "from_SG_ref "<<from_SGname<<"\n";
    //    CCtbxSym::PrintCctbxSymops(from_SG_ref);
    //    std::cout << "to_SG_ref "<<to_SGname<<"\n";
    //    CCtbxSym::PrintCctbxSymops(to_SG_ref);

    // Reference groups may be the same
    if (from_SG_ref == to_SG_ref) {
      samereferencegroup = true;
    }
    //      std::string message =
    //        CCtbxSym::SpaceGroupName(from_SG_ref.type(), 'H')+
    //        " has different reference setting from "+
    //        CCtbxSym::SpaceGroupName(to_SG_ref.type(), 'H')+"\n";
    //^std::cout <<"Message: "<<message<<"\n";
    //      Message::message(Message_warn(message));
    //      throw Message_warn(message);

    // We want the transformation from HKLIN to input
    reindex = CCtbxSym::SetReindexOp(ChB_ref_to.inverse() * ChB_ref_from);

    ////    } // same reference group

    // Test for same intensity group
    bool anom = false;  // ignore anomalous
    sgtbx::space_group from_IG = from_SG.build_derived_reflection_intensity_group(anom);
    sgtbx::space_group to_IG = to_SG.build_derived_reflection_intensity_group(anom);
    // Point groups
    sgtbx::space_group from_PG = from_SG.build_derived_point_group();
    sgtbx::space_group to_PG = to_SG.build_derived_point_group();

    if (from_IG == to_IG) {
      sameintensitygroup = true;
    }

    if (from_PG == to_PG) {
      samepointgroup = true;
    }
  }
  //--------------------------------------------------------------
  bool SpacegroupReindex(const GlobalControls& GC,
                         const hkl_symmetry& HKLINsymm, const Scell& cell,
                         ReindexOp& Reindex, const bool& failHere,
                         phaser_io::Output& output)
  // If SPACEGROUP is specified but no REINDEX operator, generate appropriate reindexing
  // to convert from input HKLIN file HKLINsymm to desired spacegroup
  // Probably really only useful (or indeed valid) for C2 <-> I2 & H3<->R3, or P222 groups
  //
  // input cell corresponds to HKLINsymm
  //
  // Returns false if Reindex is set
  // if failHere is true, fails if the symmetries do not belong to same lattice group
  {
    if (GC.Spacegroup() == "" ||  GC.Spacegroup() == "HKLIN" || GC.IsReindexSet()) return false;

    std::string HKLIN_SGname = HKLINsymm.symbol_xHM();
    hkl_symmetry NewSymm(GC.Spacegroup());
    std::string Input_SGname = NewSymm.symbol_xHM();
    if (failHere && (NewSymm.CrysSys() != HKLINsymm.CrysSys())) {
      std::string message =
        "Specified SPACEGROUP "+GC.Spacegroup()+
        " must belong to same crystal system and point group\nas the input space group "
        +HKLIN_SGname+
        " unless REINDEX is explicitly given";
      ReportErrors::printFatalError(message);
    }
    // Get reindex operator if needed
    SpacegroupReindexOp sgreindex(HKLIN_SGname, Input_SGname);

    Reindex = sgreindex.Reindex();
    if (sgreindex.sameIntensityGroup()) {
      if (Reindex.IsIdentity()) {
        output.logTab(0,LOGFILE,
                      "\nNo reindexing needed to convert from space group "
                      +HKLIN_SGname+" to "+Input_SGname);
      } else {
        output.logTab(0,LOGFILE,
                      "\nReindexing data with operator "+Reindex.as_hkl()+
                      " from space group "+HKLIN_SGname+" to "+Input_SGname);
      }
      /*
      // same intensity group, no reindexing needed
      output.logTab(0,LOGFILE,
                    "\nNo reindexing needed to convert from space group "
                    +HKLIN_SGname+" to "+Input_SGname);
      return false;
      */
    } else if (sgreindex.sameReferenceGroup()) {
      // eg C2 <-> I2, H3 <-> R3
      int AllowI2 =  GC.AllowI2();
      if (NewSymm.lattice_type() != 'I') {
        AllowI2 = 0;  // don't allow I lattice if we've asked for C2
      } else if (NewSymm.lattice_type() == 'I') {
        AllowI2 = -1;  // force I lattice if we've asked for I2
      }

      if (NewSymm.lattice_type() != 'R') {  // Don't reduce to reference group for R lattice
        // Get reindex operator from required group Input_SGname to "reference" setting
        CCtbxSym::PointGroup PG(Input_SGname);
        Scell scell = cell.change_basis(Reindex);
        PG.SetCell(scell.UnitCell(), ReindexOp(), AllowI2);
        ReindexOp newreindex = PG.RefSGreindex(); // reindex cell -> best

        //      std::cout << "\nSG: " << HKLIN_SGname<<", "<< Input_SGname<<"\n";
        //      std::cout <<"Cell: "<<cell.format()<<"\n";
        //      std::cout <<"SCell: "<<scell.format()<<"\n";
        //      std::cout <<"RefCell: "<< Scell(PG.TransformedCell()).format()<<"\n";
        //        std::cout << "Reindex " << Reindex.as_hkl() <<"\n"; //^
        //        std::cout << "newreindex " << newreindex.as_hkl() <<"\n"; //^

        Reindex = Reindex * newreindex;
      }

      output.logTab(0,LOGFILE,
                    "\nReindexing data with operator "+Reindex.as_hkl()+
                    " from space group "+HKLIN_SGname+" to "+Input_SGname);

      return true;
    } else {
      if (failHere) {
        std::string message =
          "Specified SPACEGROUP "+GC.Spacegroup()+
          " must have the same 'reference' space group as the input file symmetry "+HKLIN_SGname+
          " unless REINDEX is explicitly given";
        ReportErrors::printFatalError(message);
      }
    }
    return false;
  }
}
