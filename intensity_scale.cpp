//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#include "intensity_scale.h"

namespace clipper {
	// compile templates
	
	template class Iscale_aniso<ftype32>;
	
	template class Iscale_aniso<ftype64>;
	
	
}


