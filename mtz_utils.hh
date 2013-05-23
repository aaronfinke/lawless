// mtz_utils.hh

#ifndef MTZ_UTILS_HEADER
#define MTZ_UTILS_HEADER

// CCP4
//#include "ccp4/csymlib.h"    // CCP4 symmetry stuff
#include "ccp4/cmtzlib.h"    // CCP4 MTZlib headers (namespace CMtz)

#include "hkl_symmetry.hh"

namespace MtzIO 
{
  //--------------------------------------------------------------
  // Make clipper symop string from MTZ operators
  std::vector<clipper::Symop> ClipperSymopsFromMtzSYMGRP(const CMtz::SYMGRP& mtzsym);
  //--------------------------------------------------------------
  // Mtz symmetry from SpaceGroup
  //  HorR   H or R for rhombohedral lattice
  //  spg_status  spg_confidence   P => pointgroup correct
  //                               E => spacegroup or enantiomorph
  //                               S => spacegroup is correct
  //                               X => flag not set
  CMtz::SYMGRP  spg_to_mtz(const scala::SpaceGroup& cspgp, const char& HorR,
			   const char& spg_status);
  //--------------------------------------------------------------
  //&&&  CSym::CCP4SPG * spg_mtz_to_csym(const CMtz::SYMGRP& mtzsym);
 //--------------------------------------------------------------
  //&&&  CMtz::SYMGRP  spg_csym_to_mtz(const CSym::CCP4SPG * csym);
  //--------------------------------------------------------------
  // true if two MTZ-style symmetry structures are equal
  bool CmtzSymgrpEqual(const CMtz::SYMGRP& sg1, const CMtz::SYMGRP& sg2);
  //--------------------------------------------------------------
  /* Write spacegroup info to mtzout */
  // Copied from Clipper code, ccp4_mtz_io.cpp::write_spacegroup
  //   written by Kevin Cowtan, copied with his permission 2013/05/20
  //  spg_status  spg_confidence   P => pointgroup correct
  //                               E => spacegroup or enantiomorph
  //                               S => spacegroup is correct
  //                               X => flag not set
  void ccp4_write_spacegroup(CMtz::MTZ* mtzout, const clipper::Spacegroup& sg,
			     const char& spg_status);

}

#endif
