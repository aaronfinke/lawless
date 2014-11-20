// mtz_utils.cpp

#include "mtz_utils.hh"
#include "hkl_symmetry.hh"

// from CLIBS/mtzdata.h  extra item spg_confidence added from earlier versions
/** MTZ symmetry struct. */
// typedef struct { int spcgrp;           /**< spacegroup number */
//               char spcgrpname[MAXSPGNAMELENGTH+1];  /**< spacegroup name */
//               int nsym;             /**< number of symmetry operations */
//               float sym[192][4][4]; /**< symmetry operations
//                                          (translations in [*][3]) */
//               int nsymp;            /**< number of primitive symmetry ops. */
//               char symtyp;          /**< lattice type (P,A,B,C,I,F,R) */
//               char pgname[11];      /**< pointgroup name */
//                 char spg_confidence;  /**< L => Bravais lattice correct
//                                            P => pointgroup correct
//                                            E => spacegroup or enantiomorph
//                                            S => spacegroup is correct
//                                            X => flag not set */
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
  CMtz::SYMGRP spg_to_mtz(const scala::SpaceGroup& cspgp, const char& HorR,
                          const char& spg_status)
  {
    CMtz::SYMGRP mtzsym;
    //    mtzsym.spcgrp = cspgp.spacegroup_number();
    mtzsym.spcgrp = cspgp.CCP4_Spacegroup_number();
    strcpy(mtzsym.spcgrpname, scala::SGnameHtoR(cspgp.symbol_hm(),HorR).c_str());
    mtzsym.nsym = cspgp.num_symops();
    mtzsym.nsymp = cspgp.num_primops();
    strcpy(mtzsym.pgname, scala::SGnameHtoR(cspgp.symbol_laue(),HorR).c_str());
    mtzsym.symtyp = cspgp.LatType();
    mtzsym.spg_confidence = spg_status;

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
    //           char spcgrpname[MAXSPGNAMELENGTH+1];  /**< spacegroup name */
    //           int nsym;             /**< number of symmetry operations */
    //           float sym[192][4][4]; /**< symmetry operations
    //                                          (translations in [*][3]) */
    //           int nsymp;            /**< number of primitive symmetry ops. */
    //           char symtyp;          /**< lattice type (P,A,B,C,I,F,R) */
    //           char pgname[11];      /**< pointgroup name */
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
  //--------------------------------------------------------------
  using namespace clipper;
  /* Write spacegroup info to mtzout */
  // Modified from Clipper code, ccp4_mtz_io.cpp::write_spacegroup
  //   written by Kevin Cowtan, copied with his permission 2013/05/20
  void ccp4_write_spacegroup(CMtz::MTZ* mtzout, const scala::SpaceGroup& sg,
                          const char& spg_status)
  // Note clipper::Spacegroup scrambles symops
  {
    // tables of MTZ symbols
    char mtzlauetab[231][8]={"?","1","-1","2","2","2","m","m","m","m","2/m","2/m","2/m","2/m","2/m","2/m","222","222","222","222","222","222","222","222","222","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","4","4","4","4","4","4","-4","-4","4/m","4/m","4/m","4/m","4/m","4/m","422","422","422","422","422","422","422","422","422","422","4mm","4mm","4mm","4mm","4mm","4mm","4mm","4mm","4mm","4mm","4mm","4mm","-4m2","-4m2","-4m2","-4m2","-42m","-42m","-42m","-42m","-42m","-42m","-4m2","-4m2","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","3","3","3","3","-3","-3","312","321","312","321","312","321","32","3m1","31m","3m1","31m","3m","3m","-31m","-31m","-3m1","-3m1","-3m","-3m","6","6","6","6","6","6","-6","6/m","6/m","622","622","622","622","622","622","6mm","6mm","6mm","6mm","-62m","-62m","-6m2","-6m2","6/mmm","6/mmm","6/mmm","6/mmm","23","23","23","23","23","m-3","m-3","m-3","m-3","m-3","m-3","m-3","432","432","432","432","432","432","432","432","-43m","-43m","-43m","-43m","-43m","-43m","m-3m","m-3m","m-3m","m-3m","m-3m","m-3m","m-3m","m-3m","m-3m","m-3m"};
    unsigned int ccp4_code[][2] = { {0x90a34743,1003},{0x7851c0d1,1004},{0x6042abb6,4005},{0x458c0b8c,1005},{0x939daf66,1006},{0xd9a7d85a,1007},{0x5b2a5e1b,1009},{0x9b92779f,1010},{0x31194672,1011},{0x97e84b5c,1012},{0xc9e50c74,1013},{0xdc48ad2f,1014},{0xb0d2ac85,1015},{0xd756eb1b,1017},{0x10e6701,2017},{0xee8326bb,3018},{0x82570fbb,2018},{0xadb55f32,1059},{0xd5a0aa2d,1146},{0xd9a29bac,1148},{0xa20b8591,1155},{0xb951b4f7,1160},{0x219be015,1161},{0x1c80e47a,1166},{0xbb691c91,1167},{0x9b92779f,1010},{0xd9a29bac,1148},{0x1c80e47a,1166} };
    String mtzlaue = String("PG")+mtzlauetab[sg.spacegroup_number()];
    String mtzsymb = sg.symbol_hm();
    if ( sg.symbol_hm_ext() == "H" || sg.symbol_hm_ext() == "R" )
      mtzsymb = sg.symbol_hm_ext() + mtzsymb.substr(1);
    int mtzspgn = sg.spacegroup_number();
    for ( int s = 0; s < sizeof(ccp4_code)/sizeof(ccp4_code[0]); s++ )
      if ( ccp4_code[s][0] == sg.hash() ) mtzspgn = ccp4_code[s][1];
    // now write the records
    mtzout->mtzsymm.spcgrp = mtzspgn;
    mtzout->mtzsymm.nsym   = sg.num_symops();
    mtzout->mtzsymm.nsymp  = sg.num_primops();
    mtzout->mtzsymm.symtyp = mtzsymb[0];
    mtzout->mtzsymm.spg_confidence = spg_status;
    strncpy( mtzout->mtzsymm.spcgrpname, mtzsymb.c_str(), 11 );
    strncpy( mtzout->mtzsymm.pgname, mtzlaue.c_str(), 11 );
    for ( int i = 0; i < sg.num_symops(); i++ ) {
      for ( int j = 0; j < 3; j++ )
        for ( int k = 0; k < 3; k++ )
          mtzout->mtzsymm.sym[i][j][k] = sg.Symop(i).rot()(j,k);
      for ( int j = 0; j < 3; j++ )
        mtzout->mtzsymm.sym[i][j][3] = sg.Symop(i).trn()[j];
    }
  }
  //--------------------------------------------------------------
  //! Append to oldhistory
  //  uses version.hh to get program information
  #include "version.hh"
  std::vector<clipper::String> addToHistory
  (const std::vector<std::string> oldhistory)
  {
    char date[11];
    char time[9];
    CCP4::ccp4_utils_date(date);
    CCP4::ccp4_utils_time(time);
    clipper::String text = "From "+
      clipper::String(PROGRAM_NAME)+", version "+
      clipper::String(PROGRAM_VERSION)+
      ", run on "+clipper::String(date)+" at "+
      clipper::String(time);
    // History so far: NB new line goes at the beginning
    std::vector<clipper::String> historylines(1,text);
    for (size_t i=0; i<oldhistory.size(); i++) {
      historylines.push_back(oldhistory[i]);
    }
    return historylines;
  }
  //--------------------------------------------------------------
  //! Append to oldhistory and write to MTZ
  void MTZaddHistory(const std::vector<std::string> oldhistory,
                     CMtz::MTZ* mtzout)
  {
    // Make new history
    std::vector<clipper::String> historylines = addToHistory(oldhistory);
    // Add to MTZ
    int nlines = historylines.size();
    int Nhist;
    char line[MTZRECORDLENGTH];
    // Add in reverse order as MtzAddHistory reverses them
    for (int i=nlines-1;i>=0;--i) {
      if (historylines[i] != "") {
        strncpy(line, historylines[i].c_str(), MTZRECORDLENGTH-1);
        Nhist = MtzAddHistory(mtzout, &line, 1);
      }
    }
    Nhist = Nhist;
  }
  //--------------------------------------------------------------
}
