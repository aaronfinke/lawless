// mtz_utils.cpp

#include <pthread.h>
#include "mtz_utils.hh"
#include "hkl_symmetry.hh"

// from CLIBS/cmtzlib.h
//  typedef struct { int spcgrp;           /**< spacegroup number */
//		 char spcgrpname[MAXSPGNAMELENGTH+1];  /**< spacegroup name */
//		 int nsym;             /**< number of symmetry operations */
//		 float sym[192][4][4]; /**< symmetry operations 
//                                          (translations in [*][3]) */
//		 int nsymp;            /**< number of primitive symmetry ops. */
//		 char symtyp;          /**< lattice type (P,A,B,C,I,F,R) */
//		 char pgname[11];      /**< pointgroup name */
//               } SYMGRP;

namespace MtzIO 
{
  //--------------------------------------------------------------
  std::vector<clipper::Symop> ClipperSymopsFromMtzSYMGRP(const CMtz::SYMGRP& mtzsym)
  // Make clipper symops from MTZ operators
  {
    std::vector<clipper::Symop> symops;
    clipper::Mat33<double> rot;
    clipper::Vec3<double>  trn;
    for ( int i = 0; i < mtzsym.nsym; i++ ) {
      for (int k = 0; k < 3; ++k) {
	for (int l = 0; l < 3; ++l) {
	  rot(k,l) = mtzsym.sym[i][k][l];
	}
	trn[k] = mtzsym.sym[i][k][3];
      }
      symops.push_back(clipper::Symop(RTop<>(rot,trn)));
    }
    return symops;
  }
  //--------------------------------------------------------------
  // Mtz symmetry from SpaceGroup
  CMtz::SYMGRP spg_to_mtz(const scala::SpaceGroup& cspgp, const char& HorR)
  {
    CMtz::SYMGRP mtzsym;
    mtzsym.spcgrp = cspgp.spacegroup_number();
    strcpy(mtzsym.spcgrpname, scala::SGnameHtoR(cspgp.symbol_hm(),HorR).c_str());
    mtzsym.nsym = cspgp.num_symops();
    mtzsym.nsymp = cspgp.num_primops();
    strcpy(mtzsym.pgname, scala::SGnameHtoR(cspgp.symbol_laue(),HorR).c_str());
    mtzsym.symtyp = cspgp.LatType();

    for (int i = 0; i < mtzsym.nsym; ++i) {
      for (int k = 0; k < 3; ++k) {
	for (int l = 0; l < 3; ++l) {
	  mtzsym.sym[i][k][l] = cspgp.symop(i).rot()(k,l);
	}
	mtzsym.sym[i][k][3] = cspgp.symop(i).trn()[k];
	for (int l = 0; l < 3; ++l) 
	  mtzsym.sym[i][3][l] = 0.0;
	mtzsym.sym[i][3][3] = 1.0;
      }
    }
    return mtzsym;
  }
//--------------------------------------------------------------
bool CmtzSymgrpEqual(const CMtz::SYMGRP& sg1, const CMtz::SYMGRP& sg2)
// true if two MTZ-style symmetry structures are equal
{
  //  typedef struct { int spcgrp;           /**< spacegroup number */
  //		 char spcgrpname[MAXSPGNAMELENGTH+1];  /**< spacegroup name */
  //		 int nsym;             /**< number of symmetry operations */
  //		 float sym[192][4][4]; /**< symmetry operations 
  //                                          (translations in [*][3]) */
  //		 int nsymp;            /**< number of primitive symmetry ops. */
  //		 char symtyp;          /**< lattice type (P,A,B,C,I,F,R) */
  //		 char pgname[11];      /**< pointgroup name */
  //               } SYMGRP;
  // Don't worry about names
  if (sg1.spcgrp != sg2.spcgrp) return false;
  if (sg1.nsym != sg2.nsym) return false;
  if (sg1.nsymp != sg2.nsymp) return false;
  for (int k=0;k<sg1.nsym;++k) {
    for (int j=0;j<4;++j) {
      for (int i=0;i<4;++i) {
	if (sg1.sym[k][j][i] != sg2.sym[k][j][i]) return false;
      }}}
  if (sg1.symtyp != sg2.symtyp) return false;
  return true;
}

}
