// scala.h

#ifndef SCALA_HEADER
#define SCALA_HEADER

#define ASSERT assert
#include <assert.h>


#include <iostream>
#include <stdio.h>
#include <string>

// Clipper
#include <clipper/clipper.h>
#include "clipper/clipper-ccp4.h"
using clipper::Message;
using clipper::Message_fatal;



// CCP4
#include "csymlib.h"    // CCP4 symmetry stuff
#include "ccp4_parser.h"
#include "ccp4_general.h"
#include "ccp4_program.h"
#include "cmtzlib.h"    // CCP4 MTZlib headers (namespace CMtz)
#include "cvecmat.h"


// Scala bits
#include "util.hh" 
#include "hkl_datatypes.hh"
#include "score_datatypes.hh"
#include "spline.hh"
#include "scala_util.hh"
#include "range.hh"
#include "controls.hh"
#include "hkl_controls.hh"
#include "hash.hh"
#include "mtz_unmerge_io.hh"
#include "openinputfile.hh"
#include "hkl_symmetry.hh"
#include "hkl_unmerge.hh"
#include "hkl_merged_list.hh"
#include "normalise.hh"
#include "pgscore.hh"
#include "setscores.hh"
#include "scsignificance.hh"
#include "lattice.hh"
#include "average.hh"
#include "getsubgroups.hh"
#include "reindexscore.hh"
#include "copymergedmtz.hh"
#include "writeunmergedmtz.hh"
#include "intensitystatistics.hh"
#include "sysabstest.hh"
#include "centre.hh"
#include "interpretcommandline.hh"
#include "keywords.hh"
#include "InputAll.hh"
#include "Output.hh"
#include "datatypes.hh"

namespace scala
{
  void error_exit(int level, std::string routine, std::string message);

  // Utility things
  // Debug printing
  void print_hkl(std::string message, const Hkl& hkl);
  void print_sym(std::string message, const CSym::CCP4SPG * spacegroup);

}

#endif
