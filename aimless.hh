// aimless.hh

#ifndef AIMLESS_HEADER
#define AIMLESS_HEADER

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
#include "ccp4/csymlib.h"    // CCP4 symmetry stuff
#include "ccp4/ccp4_parser.h"
#include "ccp4/ccp4_general.h"
#include "ccp4/cmtzlib.h"    // CCP4 MTZlib headers (namespace CMtz)
#include "ccp4/cvecmat.h"


#include "interpretcommandline.hh"
#include "keywords_aimless.hh"
#include "InputAll.hh"
#include "Output.hh"
#include "printing.hh"

#include "analysesd.hh"
#include "sdmodel.hh"
#include "globalcontrols.hh"
#include "hkl_unmerge.hh"
#include "mtz_unmerge_io.hh"
#include "openinputfile.hh"

#endif
