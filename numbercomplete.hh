//
// numbercomplete.hh
//
// Number of reflections in each resolution bin in complete sphere
// 


#ifndef NUMBERCOMPLETE_HEADER
#define NUMBERCOMPLETE_HEADER

#include <vector>

#include "range.hh"
#include "hkl_symmetry.hh"
#include "hkl_datatypes.hh"

namespace scala {
void NumberComplete(const ResoRange& ResRange,
		    const hkl_symmetry& symmetry,
		    const Scell& cell,
		    std::vector<int>& Nrefres,
		    std::vector<int>& Nrefacen);
}
#endif
