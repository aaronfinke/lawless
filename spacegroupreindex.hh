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

  class SpacegroupReindexOp {
    public:
    SpacegroupReindexOp(){}

    SpacegroupReindexOp(const std::string& from_SGname,
			const std::string& to_SGname);

    ReindexOp Reindex() const {return reindex;}

    // true if groups have same intensity group,
    // ie same point group, same lattice absences, and ignoring anomalous
    bool sameIntensityGroup() const {return sameintensitygroup;}

    // true if different settings of the same space group,
    //  eg C2 & I2, or R3 & H3
    bool sameReferenceGroup() const {return samereferencegroup;}

  private:
    bool status;  // true if set
    // true if groups have same intensity group,
    // ie same point group, same lattice absences, and ignoring anomalous
    bool sameintensitygroup;
    // true if different settings of the same space group,
    //  eg C2 & I2, or R3 & H3
    bool samereferencegroup;

    ReindexOp reindex;  // reindex from -> to, == h,k,l if !samereferencegroup
  };


  // If SPACEGROUP is specified but no REINDEX operator, generate appropriate reindexing
  // to convert from input HKLIN file HKLINsymm to desired spacegroup
  // Probably really only useful (or indeed valid) for C2 <-> I2 & R3 <-> H3
  //
  // Returns true if Reindex is set
  // fails if the symmetries do not belong to same lattice group
  bool SpacegroupReindex(const GlobalControls& GC,
			 const hkl_symmetry& HKLINsymm, const Scell& cell,
			 ReindexOp& Reindex, const bool& failHere,
			 phaser_io::Output& output);
}
#endif
