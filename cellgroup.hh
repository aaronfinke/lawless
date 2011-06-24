// No CCtbx types are used in the external interface

#ifndef CELLGROUP_HEADER
#define CELLGROUP_HEADER

#include <cctbx/crystal/symmetry.h>
#include <cctbx/sgtbx/lattice_symmetry.h>
#include <cctbx/uctbx/fast_minimum_reduction.h>

using namespace cctbx;

#include "hkl_symmetry.hh"
// CCP4
#include "csymlib.h"    // CCP4 symmetry stuff

// Forward declarations
namespace scala{
  class hkl_symmetry;
  class Scell;
}
//namespace CSym{class CCP4SPG;}

namespace CCtbxSym
{

  //--------------------------------------------------------------
  class CellGroup
  {
    // A spacegroup object used to impose cell constraints
  public:
    CellGroup(){}
    // Construct from SpaceGroup object
    CellGroup(const scala::SpaceGroup& spgp);

    // Return cell constrained by symmetry
    scala::Scell constrain(const scala::Scell& cell) const;

  private:

    sgtbx::space_group group;
  };
}
#endif
