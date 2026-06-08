// rejectbatches.hh
//
// Reject batches if their scales are way out of line with the others
//

#ifndef REJECTBATCHES_HEADER
#define REJECTBATCHES_HEADER

#include "hkl_unmerge.hh"
#include "scalemodel.hh"
#include "controls.hh"
#include "Output.hh"

namespace scala {

  class RejectBatches {
  public:
    RejectBatches(){}
    
    RejectBatches(hkl_unmerge_list& hkl_list, ScaleModel& AllScales,
		  const all_controls& controls,
		  phaser_io::Output& output);
    
  private:

  };
}
#endif
