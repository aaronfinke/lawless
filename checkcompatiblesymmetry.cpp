//
//  checkcompatiblesymmetry.cpp
//

#include "checkcompatiblesymmetry.hh"
#include "scala_util.hh"
#include "report_errors.hh"


namespace scala{
  //--------------------------------------------------------------
  bool CheckCompatibleSymmetry(const hkl_symmetry& RefSym,
                               const hkl_symmetry& TestSym,
                               const bool& TestDataMerged)
  // Check that test set and reference sets are compatible
  // (1) they should have the same crystal system even if different Laue group
  //     return true if same crystal system
  // (2) if the test set is merged, then they should have the same Laue group
  //     Fails (fatal) if this is not so
  //     Return true if same Laue group
  {
    bool OK = true;
    if (! (TestSym.CrysSys() == RefSym.CrysSys())) {
      OK = false;
    }
    if (!TestDataMerged) {return OK;}

    if (OK) {
      SpaceGroup pattSGref = RefSym.GetSpaceGroup().PattersonGroup();
      SpaceGroup pattSGtest = TestSym.GetSpaceGroup().PattersonGroup();

      // Only test Laue group if merged
      bool SameLaueGroup = (pattSGref == pattSGtest);

      if (!SameLaueGroup) {
        // Special for I2 / C2, allowed
        std::string nameref  = pattSGref.Symbol_hm();
        std::string nametest = pattSGtest.Symbol_hm();
        if (((nameref == "C 1 2/m 1") && (nametest == "I 1 2/m 1")) |   \
            ((nameref == "I 1 2/m 1") && (nametest == "C 1 2/m 1"))) {
          //OK
          SameLaueGroup = true;
        } else {
          OK = false;
        }
      }
    }

    if (!OK) {
      // Merged test data must have same Laue group as reference set
      std::string
        error("Merged test dataset (HKLIN) has different Laue symmetry to reference set");
      ReportErrors::printFatalError
        (error+"\n**** Incompatible symmetries ****");
    }
    return true;
  }
}
