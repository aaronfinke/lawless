//
// numbercomplete.cpp
//
// Number of reflections in each resolution bin in complete sphere
// 

#include "numbercomplete.hh"

// Clipper
#include <clipper/clipper.h>


namespace scala {
void NumberComplete(const ResoRange& ResRange,
		    const hkl_symmetry& symmetry,
		    const Scell& cell,
		    std::vector<int>& Nrefres,
		    std::vector<int>& Nrefacen)
// returns:
//   Nrefres    number in sphere for each resolution bin
//   Nrefacen   number acentric in sphere for each resolution bin 
{
  // Create clipper reflection list in P1 (hemisphere)
  double resmax = ResRange.ResHigh();
  clipper::HKL_info hkl_info_list(clipper::Spacegroup(clipper::Spacegroup::P1),
				  cell.ClipperCell(),
				  clipper::Resolution(resmax), true);
  
  clipper::HKL_info::HKL_reference_index ih;
  for (ih = hkl_info_list.first(); !ih.last(); ih.next()) {
    Hkl hkl(ih.hkl());
    // test for systematic absences (including lattice absences)
    if (symmetry.epsilon(hkl) > 0) { 
      clipper::ftype invresolsq = ih.invresolsq();
      int mres = ResRange.bin(invresolsq);
      Nrefres[mres] += 2;  // for sphere
      if (!symmetry.is_centric(hkl)) {
	Nrefacen[mres] += 2;
      }
    }
  }
}
}
