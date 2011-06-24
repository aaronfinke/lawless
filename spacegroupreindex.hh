// spacegroupreindex.hh

#ifndef SPACEGROUPREINDEX_HEADER
#define SPACEGROUPREINDEX_HEADER

#include "hkl_datatypes.hh"
#include "globalcontrols.hh"
#include "hkl_symmetry.hh"
#include "Output.hh"

using phaser_io::LOGFILE;

namespace scala {
  //--------------------------------------------------------------
  //! return reindex operator to convert alternative space groups in same crystal system
  /*! Mainly (only?) for C2<->I2 and H3<->R3
    Throws a clipper::Message_warn exception if the space groups are
    not just alternative settings (ie have different reference groups) */
  ReindexOp SpacegroupReindexOp(const std::string& from_SGname,
				const std::string& to_SGname);

  // If SPACEGROUP is specified but no REINDEX operator, generate appropriate reindexing
  // to convert from input HKLIN file HKLINsymm to desired spacegroup
  // Probably really only useful (or indeed valid) for C2 <-> I2 & R3 <-> H3
  //
  // Returns true if Reindex is set
  // fails if the symmetries do not belong to same lattice group
  bool SpacegroupReindex(const GlobalControls& GC,
			 const hkl_symmetry& HKLINsymm, const Scell& cell,
			 ReindexOp& Reindex, phaser_io::Output& output);
}
#endif
