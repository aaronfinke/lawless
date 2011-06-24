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
#include "csymlib.h"    // CCP4 symmetry stuff
#include "ccp4_parser.h"
#include "ccp4_general.h"
#include "cmtzlib.h"    // CCP4 MTZlib headers (namespace CMtz)
#include "cvecmat.h"


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
