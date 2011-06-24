// mtz_utils.hh

#ifndef MTZ_UTILS_HEADER
#define MTZ_UTILS_HEADER

// CCP4
//#include "csymlib.h"    // CCP4 symmetry stuff
#include "cmtzlib.h"    // CCP4 MTZlib headers (namespace CMtz)

#include "hkl_symmetry.hh"

namespace MtzIO 
{
  //--------------------------------------------------------------
  // Make clipper symop string from MTZ operators
  std::vector<clipper::Symop> ClipperSymopsFromMtzSYMGRP(const CMtz::SYMGRP& mtzsym);
  //--------------------------------------------------------------
  // Mtz symmetry from SpaceGroup
  CMtz::SYMGRP  spg_to_mtz(const scala::SpaceGroup& cspgp, const char& HorR);
  //--------------------------------------------------------------
  //&&&  CSym::CCP4SPG * spg_mtz_to_csym(const CMtz::SYMGRP& mtzsym);
 //--------------------------------------------------------------
  //&&&  CMtz::SYMGRP  spg_csym_to_mtz(const CSym::CCP4SPG * csym);
  //--------------------------------------------------------------
  // true if two MTZ-style symmetry structures are equal
  bool CmtzSymgrpEqual(const CMtz::SYMGRP& sg1, const CMtz::SYMGRP& sg2);

}

#endif
