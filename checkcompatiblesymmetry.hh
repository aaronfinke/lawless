//
// checkcompatiblesymmetry.hh
//

#ifndef CHECKCOMPATIBLESYMMETRY_HEADER
#define CHECKCOMPATIBLESYMMETRY_HEADER

#include "hkl_symmetry.hh"

namespace scala {
  //--------------------------------------------------------------
  // Check that test set and reference sets are compatible
  // (1) they should have the same crystal system even if different Laue group
  // (2) if the test set is merged, then they should have the same Laue group
  // Fails (fatal) if this is not so
  // Return true if same Laue group
  bool CheckCompatibleSymmetry(const hkl_symmetry& RefSym,
			       const hkl_symmetry& TestSym,
			       const bool& TestDataMerged);
}

#endif
