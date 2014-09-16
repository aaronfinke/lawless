//  GetSubGroups.hh

#ifndef GETSUBGROUPS_HEADER
#define GETSUBGROUPS_HEADER

#include "pointgroup.hh"
#include "hkl_symmetry.hh"

namespace scala
{
  //--------------------------------------------------------------
  // Make list of subgroups
  // if forceI2 true, convert any C2/m setting to I 
  std::vector<CCtbxSym::PointGroup> GetSubGroups(const hkl_symmetry& symm);
  //--------------------------------------------------------------
  // If file is merged, remove all subgroups which are subgroups of original group
  //  since these are implicit
  void RemoveImplicitSubgroups(std::vector<CCtbxSym::PointGroup>& subgroups,
			       const hkl_symmetry& OrigSymm, const ReindexOp& reindex);

} // scala

#endif
