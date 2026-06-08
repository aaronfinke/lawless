// writeoutputfiles.hh
//
// Write output reflection files as specified by outputcontrols
//

#ifndef WRITEOUTPUTFILES_HEADER
#define WRITEOUTPUTFILES_HEADER

#include <string>

#include "hkl_unmerge.hh"
#include "sdmodel.hh"
#include "globalcontrols.hh"
#include "Output.hh"
#include "mergedlist.hh"

namespace scala {
// ---------------------------------------------------------
  void WriteMergedOutputFiles(const MergedList & mergedlist,
			      const OutputControls& outputcontrols,
			      phaser_io::Output& output);
// ---------------------------------------------------------
  void WriteUnmergedOutputFiles(const std::string& title,
				const hkl_unmerge_list& hkl_list,
				const SDmodel& SDM,
				const float& Imax,
				const OutputControls& outputcontrols,
				phaser_io::Output& output);
// ---------------------------------------------------------
  // Create file containing columns
  //  h,k,l
  //  IMEAN, SIGIMEAN
  //  I(+), SIGI(+), I(-), SIGI(-)
  //
  //  outputformat  = +1  MTZ, +2 Scalepack
  void WriteMergedOutput(const int& outputformat,
			 const MergedList & mergedlist,
			 const OutputControls& outputcontrols,
			 phaser_io::Output& output);
}

#endif
