//
// scalerefine.cpp
//
// Main scaling routine
//


#include "scalerefine.hh"
#include "refinescale.hh"
#include <assert.h>
#define ASSERT assert

namespace scala {
// ---------------------------------------------------------
  void ScaleRefine(hkl_unmerge_list& hkl_list, ScaleModel& AllScales,
		   const SDmodel& SDM,
		   const all_controls& controls, const int& Ncycles,
		   const bool& print, phaser_io::Output& output)
  //  Main scaling
  // print = true to print scales
  // 
  {
    // Set up up refinement object:
    //  store addresses of reflection & scale objects 
    RefineScale refscl(hkl_list, AllScales, SDM,
		       controls.refinecontrol.Nprocs());

    // default protocols
    phaser::protocolPtr cPtr(new phaser::ProtocolScale(Ncycles));

    phaser::Minimizer Min;
    Min.run(refscl, cPtr, output);    // run minimiser

    AllScales.NormaliseParameters();  // Normalise result
    if (print) AllScales.PrintScales(output);

  }  // ScaleRefine

} // namespace scala
