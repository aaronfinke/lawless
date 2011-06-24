//
//  checkcompatiblesymmetry.cpp
//

#include "checkcompatiblesymmetry.hh"
#include "scala_util.hh"


namespace scala{
  //--------------------------------------------------------------
  bool CheckCompatibleSymmetry(const hkl_symmetry& RefSym,
			       const hkl_symmetry& TestSym,
			       const bool& TestDataMerged)
  // Check that test set and reference sets are compatible
  // (1) they should have the same crystal system even if different Laue group
  // (2) if the test set is merged, then they should have the same Laue group
  // Fails (fatal) if this is not so
  // Return true if same Laue group
  {  
    if (! TestSym.CrysSys() == RefSym.CrysSys()) {
      std::string
	error("Test dataset (HKLIN) has different lattice symmetry to reference set");
      Message::message(Message_fatal
		       (error+"\n**** Incompatible symmetries ****"));
    }
    
    bool SameLaueGroup = RefSym.GetSpaceGroup().PattersonGroup() ==
      TestSym.GetSpaceGroup().PattersonGroup();

    if (!SameLaueGroup && TestDataMerged) {
      // Merged test data must have same Laue group as reference set
	std::string
	  error("Merged test dataset (HKLIN) has different Laue symmetry to reference set");
	Message::message(Message_fatal
		       (error+"\n**** Incompatible symmetries ****"));
    }
    return SameLaueGroup;
  }
}
