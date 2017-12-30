//  hkl_symmetry.cpp

#include <algorithm>

#include "hkl_symmetry.hh"
#include "matvec_utils.hh"
#include "string_util.hh"
#include "report_errors.hh"
#include "spacegroupreindex.hh"

namespace scala {
  //--------------------------------------------------------------
  LatticeCenteringList::LatticeCenteringList()
  // build library list
  {
    // Omit 'R' as it would just come up as primitive
    char AllowedTypes[] = {'P','A','B','C','I','F', 'H'};
    nlattices = sizeof(AllowedTypes);
    latticetypes.assign(AllowedTypes, AllowedTypes+nlattices);
    latcenops.resize(nlattices);

    for (int j=0;j<nlattices;++j) {
      // centred triclinic
      std::string sgname = std::string(1, latticetypes[j])+" 1";
      SpaceGroup SG1(sgname);
      int ncent = SG1.num_centering_symops();
      for (int i=0;i<ncent;++i) {
        latcenops[j].push_back(SG1.centering_symop(i));
      }
    }
  }
  //--------------------------------------------------------------
  // test whether cenops corresponds to a valid lattice centering,
  // in any order. Returns lattice type, or ' ' if not found
  char LatticeCenteringList::findCenteringOps
  (std::vector<clipper::Symop> cenops)
  {
    char lattyp = ' ';
    const double TOLERANCE = 0.001;
    if (cenops.size() == 0) {return lattyp;}
    for (int j=0;j<nlattices;++j) { // include P lattice
      int found = 0;
      for (size_t i=0;i<latcenops[j].size();++i) {
        for (size_t k=0;k<cenops.size();++k) {
          if (latcenops[j][i].equals(cenops[k], TOLERANCE)) {
            // op k == op j,i
            found++;
          }
        }
      }
      if (found == cenops.size()) {
        // all cenops found, so it is lattice type j
        lattyp = latticetypes[j];
        break;
      }
    }
    return lattyp;
  }
  //--------------------------------------------------------------
  // return centering operators for lattice type, "standard" (clipper) order
  std::vector<clipper::Symop> LatticeCenteringList::standardCenteringOps
  (const char& lattype) const
  {
    std::vector<clipper::Symop> cenops;
    std::vector<char>::const_iterator it =
      std::find(latticetypes.begin(), latticetypes.end(), lattype);
    if (it != latticetypes.end()) {
      int j = it - latticetypes.begin();
      cenops = latcenops.at(j);
    }
    return cenops;
  }
  //--------------------------------------------------------------
  std::string LatticeCenteringList::format() const
  {
    std::string s = "Possible centred lattices\n";
    for (int j=1;j<nlattices;++j) { // skip P lattice
      s += "\nLattice type : " + std::string(1, latticetypes[j]) + "\n";
      for (size_t i=0;i<latcenops[j].size();++i) {
        s += latcenops[j][i].format() + "\n";
      }
    }
    return s;
  }
  //--------------------------------------------------------------
  // for sorting
  bool operator < (const SymElement& a, const SymElement& b)
  // Sort on rotation order and axis direction (a,b,c)
  // 1. rotation order
  // 2. principle axis, c < b < a
  // 3. Sum(a+b+c) = 0
  //
  // This has the effect of sorting together trigonal dyads belonging to either
  // 321 or 312 point groups
  //  ( this isn't very clever)
  {
    if (a.Nfold != b.Nfold) return a.Nfold < b.Nfold;  // different rotation order
    // Along principle axis?
    int asum = abs(a.iaxis[0])+abs(a.iaxis[1])+abs(a.iaxis[2]);
    int bsum = abs(b.iaxis[0])+abs(b.iaxis[1])+abs(b.iaxis[2]);
    bool ap = (asum == 1);
    bool bp = (bsum == 1);
    if (ap && !bp) return true;  // a principle, b not
    if (!ap && bp) return false; // b principle, a not
    if (ap && bp) {
      // both principle axes, c < b < a
      int adx = MVutil::LargestVectorElement<int>(a.iaxis);
      int bdx = MVutil::LargestVectorElement<int>(b.iaxis);
      return adx > bdx;
    }
    // neither
    asum = a.iaxis[0]+a.iaxis[1]+a.iaxis[2];
    bsum = b.iaxis[0]+b.iaxis[1]+b.iaxis[2];
    return asum < bsum;
  }
  //--------------------------------------------------------------
  bool IsSymopIdent(const clipper::Symop& op)
  // Local function
  {
    clipper::RTop<double>  Iop =  clipper::RTop<double>::identity();
    return clipper::RTop<double>(op).equals(Iop,0.0001);
  }
  //--------------------------------------------------------------
  std::string SGnameHtoR(const std::string& sgname, const char& HorR)
  //! Change H to R or vv in space group name
  //! If 1st character of name is H or R, change to :
  //!   1. character after ":" if H or R
  //!   2. HorR if 'H' or 'R' [default 'H']
  //!   3. else leave
  //! eg  R 3 2 -> H 3 2
  //!     R 3 :H -> H 3
  {
    // remove leading and trailing spaces
    std::string sname = StringUtil::Trim(sgname);
    // First character should be a lattice type
    char lattype = toupper(sname[0]);
    char horr = toupper(HorR);
    if (lattype == 'H' || lattype == 'R') { // only do anything if H or R
      std::vector<std::string> parts = StringUtil::split(sname, ":");
      if (parts.size() > 1) {
        std::string ext = StringUtil::Trim(parts[1]);
        if (ext.size() == 1) ext = std::string(1,toupper(ext[0]));
        if (ext == "H" || ext == "R") {
          lattype = ext[0];
          sname = StringUtil::Trim(parts[0]);  // strip off ":X"
        }
      } else if (horr == 'H' || horr == 'R') {
        lattype = horr;
      }
      if (!scala::AllowedLatticeType(lattype)) {
        ReportErrors::printFatalError(clipper::String("Illegal lattice type ")+lattype);
      }
      sname = StringUtil::Trim(sname.substr(1));
      std::string lt = std::string(1,lattype);
      if (sgname.find(" ") != std::string::npos) {
        // if spaces in original name, add a space after lattice type
        lt += " ";
      }
      sname = lt+sname;
    }
    return sname;
  }
  // =============================================================
  //  class SpaceGroup
  //    constructors
  //--------------------------------------------------------------
  SpaceGroup::SpaceGroup(std::vector<clipper::Symop>& symops)
    : Nsymp(0), status(0), statusmessage("")
  {init(symops);}
  SpaceGroup::SpaceGroup(const int& SpgNumber)
    : Nsymp(0), status(0), statusmessage("")
  {init(SpgNumber);}
  SpaceGroup::SpaceGroup(const clipper::Spacegroup& ClpSG)
    : Nsymp(0), status(0), statusmessage("")
  {init(ClpSG);}
  //--------------------------------------------------------------
  SpaceGroup::SpaceGroup(const std::string& spgname)
    : Nsymp(0), status(0), statusmessage("")
  {SpaceGroup::init(spgname);}
  //--------------------------------------------------------------
  bool SpaceGroup::init(std::vector<clipper::Symop>& symops,
                        const bool& softfail)
  // if softfail = true, then return false on error, else fail here
  {
    csymops = symops;
    clipper::String symopString;
    for (size_t i=0;i<csymops.size();++i) {
      //      std::cout <<"Symop "<<i<<" "<<csymops[i].format()<<std::endl;
      symopString += csymops[i].format() + ";";
    }
    // Symop string ("x,y,z; ...")
    //    std::cout <<"\nclipper::Spacegroup constructed from operators: "
    //        << symopString<<std::endl;

    clipper::Spgr_descr spdescr(symopString, clipper::Spacegroup::Spgr_descr::Symops);
    //    std::cout << spdescr.spacegroup_number() <<"  ; "
    //        << spdescr.symbol_xhm() << std::endl;
    if (spdescr.spacegroup_number() <= 0) {
      if (softfail) {
        return false;
      }
      ReportErrors::printFatalError
        ("SpaceGroup initialisation from illegal operators (maybe from incorrect reindex):\n"+
         symopString);
    }

    clipper::Spacegroup::init(clipper::Spgr_descr(symopString,
                 clipper::Spacegroup::Spgr_descr::Symops));
    init();
    return true;
  }
  //--------------------------------------------------------------
  void SpaceGroup::init(const std::string& spgname)
  // Create from name
  {
    spacegroupname = spgname;
    // Change names with lattice type 'H' to 'R' for clipper construction
    std::string name = SGnameHtoR(spgname, 'R');
    char lt1 = StringUtil::Trim(spgname)[0];  // original lattice type
    clipper::Spgr_descr spgd;
    if (lt1 == 'R') {
      // special for rhombohedral setting of R lattice
      // Generate symbol for Clipper eg R 3 :R
      if (name.find(" ") != std::string::npos) {name += " ";}
      name = name+":R"; // add trailing :R
    }
    if (name.find(" ") == std::string::npos) {
      // if no spaces, add them by translating from the library file
      CSym::CCP4SPG* CCP4sg = CSym::ccp4spg_load_by_spgname(spgname.c_str());
      ////      if (CCP4sg == NULL) return;
      if (CCP4sg == NULL) {
        name = spgname;
      } else {
        name = CCP4sg->symbol_xHM;
      }
    }
    try {
      if (name.find(":") == std::string::npos) {
        // Special for centred triclinic
        if (isNameCentredTriclinic(name)) {
          // name modified to "X 1"
          spgd = clipper::Spgr_descr(name, Spacegroup::Hall);
        } else {
          spgd = clipper::Spgr_descr(name, Spacegroup::HM);
        }
      } else {
        spgd = clipper::Spgr_descr(name, Spacegroup::XHM);  // name contains ":"
      }
    }
    catch (Message_fatal) {
      // Failed to create clipper Spacegroup
      // try to create from symops via CCP4 library
      CSym::CCP4SPG* CCP4sg = CSym::ccp4spg_load_by_spgname(spgname.c_str());
      if (CCP4sg == NULL) return;
      std::string symopstring =  AllSymopsfromCCP4(CCP4sg);
      spgd = clipper::Spgr_descr(symopstring, Spacegroup::Symops);
    }
    clipper::Spacegroup::init(spgd);
    init();
  }
  //--------------------------------------------------------------
  bool SpaceGroup::isNameCentredTriclinic(std::string& name)
  // returns modified name if true: static
  {
    //    std::cout <<"isNameCentredTriclinic "<<name<<"\n"; //^
    std::string sname = StringUtil::Strip(name); // remove spaces
    if (sname.size() != 2) return false;  // test for "X1"
    if (sname[1] != '1') return false;  // not triclinic
    char lattices[] = {'P','A','B','C','I','F','R','H'};
    int nltypes = sizeof(lattices);;
    bool centred = false;
    for (int i=1;i<nltypes;++i) { // not P
      if (sname[0] == lattices[i]) {
        centred = true;
        name = std::string(1,lattices[i])+" 1";
      }
    }
    //    if (centred) std::cout << "centred\n"; //^
    return centred;
  }
  //--------------------------------------------------------------
  void SpaceGroup::init(const int& SpgNumber)
  // Create from number
  {
    clipper::Spacegroup::init(clipper::Spgr_descr(SpgNumber));
    init();
  }
  //--------------------------------------------------------------
  // Construct from clipper object
  void SpaceGroup::init(const clipper::Spacegroup& ClpSG)
  {
    clipper::Spacegroup::init(ClpSG);
    init();
  }
  //--------------------------------------------------------------
  void SpaceGroup::init()
  {
    if (csymops.size() == 0) {
      // Copy symops
      for (int i=0;i<num_symops();++i) {
        csymops.push_back(clipper::Spacegroup::symop(i));
      }
    }
    Nsymp = num_primops();
    // Generate list of rotation operators & inverse, discarding translations
    rotsymops.resize(csymops.size());
    invrotsymops.resize(csymops.size());
    for (size_t i=0;i<csymops.size();++i) {
      rotsymops[i] = clipper::Symop(clipper::RTop<>(csymops[i].rot()));
      invrotsymops[i] = clipper::Symop(clipper::RTop<>(csymops[i].rot()).inverse());
    }
    std::string spgname = spacegroupname;
    if (spgname == "Unknown") {
      spgname = CCP4spaceGroupName();  // construct from csymops
    } else if (spgname == "") {
      spgname = "Unknown";
    }
    spacegroupname = spgname;
    spacegroupnumber = spacegroup_number();
    if (spacegroupnumber == 0) {
      CCP4spacegroupnumber = 0;
    } else {
      spacegroupname = symbol_hm();
      CCP4spacegroupnumber = spacegroupnumber;
    }
    // Sort out name if unknown
    //^    std::cout << "SpaceGroup: HM: " << symbol_hm()
    //^       << "  Hall: " << symbol_hall() << "\n";
    // Always consult CCP4 libraries for space group numbers
    if (isNameCentredTriclinic(spacegroupname)) {
      spacegroupnumber = 0;
      CCP4spacegroupnumber = 0;
    } else {
      CCP4spaceGroupNumber();
    }
    SetLatType();
    // Change name to match lattice type if H or R
    spacegroupname = SGnameHtoR(spacegroupname, lattype);
  }
  //--------------------------------------------------------------
  //! Space group name
  std::string SpaceGroup::symbol_Hall() const
  {
    std::string hallSymbol = symbol_hall();
    if (hallSymbol == "Unknown") {
      if (CCP4spacegroupnumber == 0) {
        hallSymbol = spacegroupname;
      } else {
        hallSymbol =
          std::string(CSym::ccp4spg_load_by_ccp4_num(CCP4spacegroupnumber)->symbol_Hall);
      }
    }
    return hallSymbol;
  }
  //--------------------------------------------------------------
  clipper::HKL SpaceGroup::put_in_asu(const clipper::HKL& hkl, int& isym) const
  {
    for (int i=0;i<Nsymp;++i) {
      clipper::HKL h = hkl.transform(rotsymops[i]);
      if (recip_asu(h)) {
        isym = 2*i +1;
        return h;
      }
      if (recip_asu(-h)) {
        isym = 2*i +2;
        return -h;
      }
    }
    isym = 0;
    return clipper::HKL(0,0,0);
  }
  //--------------------------------------------------------------
  clipper::HKL SpaceGroup::get_from_asu(const clipper::HKL& hkl, const int& isym) const
  {
    // Check valid isym
    if (isym < 1 || isym > 2*Nsymp)
      ReportErrors::printFatalError("get_from_asu - ISYM out of range");
    clipper::HKL h = hkl.transform(invrotsymops[(isym-1)/2]);
    if (isym%2 == 0) {return -h;}
    else {return h;}
  }
  //--------------------------------------------------------------
  //! return the symop corresponding to Isym from put_in_asu
  clipper::Symop SpaceGroup::SymopFromIsym(const int& isym) const
  {
    if (isym < 1 || isym > 2*Nsymp)
      ReportErrors::printFatalError("SymopFromIsym - ISYM out of range");
    return csymops[(isym-1)/2];
  }
  //--------------------------------------------------------------
  SpaceGroup  SpaceGroup::NewLatticePointGroup(const char& Lattype) const
  //! Return group with new lattice, keeping point group operators
  {
    clipper::Spgr_descr::Symop_codes symcodes = generator_ops().pgrp_ops();
    //^
    //    for (int i=0;i<symcodes.size();++i) {
    //      std::cout << "SymCode " << i << " " << symcodes[i].symop().format() << "\n";
    //    } //^-
    clipper::String lat1 = std::string(1,Lattype)+" 1";  // P 1, C 1, I 1 etc
    clipper::Spgr_descr latops(lat1);
    clipper::Spgr_descr::Symop_codes latcodes = latops.generator_ops();
    for (size_t i=0;i<latcodes.size();++i) {
      //^      std::cout << "LatCode " << i << " " << latcodes[i].symop().format() << "\n"; //^-
      symcodes.push_back(latcodes[i]);
    }
    symcodes.expand();
    SpaceGroup SG;
    try {
      clipper::Spacegroup CSG = clipper::Spacegroup(clipper::Spgr_descr(symcodes));
      SG = SpaceGroup(CSG);
    }
    catch (Message_fatal) {
      /// return null group
    }
    return SG;
  }
  //--------------------------------------------------------------
  // All codes for ordered operators from clipper group
  std::vector<clipper::Symop_code> SpaceGroup::SymopCodes() const
  {
    std::vector<clipper::Symop_code> codes;
    for (int i=0;i<num_symops();++i) {
      codes.push_back(clipper::Symop_code(symop(i)));
    }
    return codes;
  }
  //--------------------------------------------------------------
  int SpaceGroup::ChangeBasis(const scala::ReindexOp& H)
  // Change basis by reindex operator
  // This is (I think) only valid if the reindexing operation preserves the point group
  // Also lattice centering operations need care
  // Probably only useful for permuting axes where appropriate, and I2<->C2
  // Returns status: = 0, OK, no lattice change
  //                 = -1 not OK, fail
  //                 = -2 not OK,
  //                 = -3 OK, origin shifted (translations not reindexed)
  //                 = +1 OK, invalid lattice type changed to P
  //                 = +2 OK, lattice type changed
  {
    reindexused = H;
    // List centering ops
    int ncent = num_centering_symops();
    std::vector<clipper::Symop> censymops(ncent);
    for (int i=0;i<ncent;++i) {
      censymops[i] = clipper::Spacegroup::centering_symop(i);
      //      std::cout <<"Centering Op: "<<censymops[i].format()<<"\n";
    }
    // Now deal with lattice centering
    std::vector<clipper::Symop> rcensymops(ncent);
    for (int i=0;i<ncent;++i) {
      rcensymops[i] = H.Symop(censymops[i]);
      //      std::cout << "Transformed cenop from " <<censymops[i].format()
      //                << " to " <<rcensymops[i].format()<<"\n";
    }
    LatticeCenteringList latcenlist;
    newlattype = latcenlist.findCenteringOps(rcensymops);  // reindexed lattice type
    // Get new centering operators in "standard" order
    std::vector<clipper::Symop> newcensymops =
      latcenlist.standardCenteringOps(newlattype);

    // Transform primitive operators (first in list)
    std::vector<clipper::Symop> newsymops(Nsymp);
    for (size_t js=0;js<Nsymp;js++) {
      //      std::cout <<"ChangeBasis Symop in:  " << csymops[js].format() <<"\n";
      // If [H] is the reindex operator,
      //  [S'] = [H]^-1 [S] [H]
      newsymops[js] = H.Symop(csymops[js]);  // now done in ReindexOp class
      //      clipper::Mat33<double> S = csymops[js].rot();
      //      clipper::Vec3<double>  t = csymops[js].trn();
      //      clipper::RTop<double> Sp(HR.inverse()*S*HR, HR.inverse() * t);
      //      csymops[js] = clipper::Symop(Sp);
      //      std::cout <<"*** ChangeBasis Symop out: " << newsymops[js].format() <<"\n";
    }

    // --- Messages
    //std::cout <<"newlattype: "<<newlattype<<"|\n";
    // newlattype = ' ' if reindexed centering operators are invalid
    //                  so reset to 'P'
    if (newlattype == ' ') {
      status = +1;
      statusmessage += "Invalid reindexed lattice centering\n";
      newlattype = 'P';
    } else if (newlattype != lattype) {
      status = +2;
      //      statusmessage += "Lattice type changed from "+std::string(1,lattype)+
      //        " to "+std::string(1,newlattype)+"\n";
      statusmessage = "Centering operators changed by reindexing "+H.as_hkl();
      if (newlattype == 'P') {
        statusmessage += std::string(" to a primitive lattice P\n")+
          "New operators would have been ";
      } else {
        statusmessage += " to lattice of type "+std::string(1,newlattype)+
          "\nNew operators: ";
      }
      for (int i=0;i<ncent;++i) {
        statusmessage += "("+rcensymops[i].format()+")";
        if (i < ncent-1) {statusmessage += ", ";}
      }
    }
    //    std::cout << "Status:\n"<<statusmessage<<"\n";
    // ---
    std::string spgname = spacegroupname;  // save space group name
    newSG(newcensymops, Nsymp, newsymops);
    if (status < 0) {
      // failed, check that pointgroup is unchanged
      SpaceGroup PG = PointGroup();
      std::vector<clipper::Symop> newpgsymops(PG.num_primops());
      for (size_t js=0;js<PG.num_primops();js++) {
        newpgsymops[js] = PG.Symop(js);
        //std::cout <<"ChangeBasis PG Symop in:  " << newpgsymops[js].format() <<"\n";
        // If [H] is the reindex operator,
        //  [S'] = [H]^-1 [S] [H]
        newpgsymops[js] = H.Symop(newpgsymops[js], true);
        //      std::cout <<"*** ChangeBasis Symop out: "
        //                << newpgsymops[js].format() <<"\n";
      }
      SpaceGroup RPG(newpgsymops);
      if (PG == RPG) {
        // OK, just use existing space group
        init(spgname);
        // success, set status to indicate possible origin shift
        status = -3;
        statusmessage = "Reindex operator "+H.as_hkl()+
          " converts space group "+spgname+
          " to an invalid space group in the same point group,"+
          " possibly with an origin shift";
      } else {
        // still wrong
        statusmessage = "Reindex operator "+H.as_hkl()+
          " converts space group "+spgname+" to an invalid space group";
        status = -1;
      }
    }
    return status;
  }
  //--------------------------------------------------------------
  void SpaceGroup::newSG(const std::vector<clipper::Symop>& newcensymops,
                         const int& Nsymp,
                         std::vector<clipper::Symop>& newsymops)
  // try to build new space group from reindexed operators
  // alters newsymops, adding lattice centering
  {
    if (newcensymops.size() > 1) { // non-primitive
      for (size_t i=1;i<newcensymops.size();++i) { // skip 1st identity
        //      std::cout <<"New Centering Op: "<<newcensymops[i].format()<<"\n";
        // Apply them to primitive operators
        newsymops.resize(Nsymp);
        for (size_t js=0;js<Nsymp;js++) {
          clipper::Symop cs(newsymops[js]*newcensymops[i]);
          newsymops.push_back(cs);
        }
      }
    }
    //    for (size_t i=0;i<newsymops.size();++i) {
    //      std::cout <<"New symop: " <<newsymops[i].format() <<"\n";
    //    }

    spacegroupname = "Unknown";  // force review of name
    bool validSG =  init(newsymops, true);
    if (!validSG) {
      status = -1;
    }
  }
  //--------------------------------------------------------------
  bool SpaceGroup::IsSymopIdentity(const int& symN) const
  {
    return IsSymopIdent(csymops.at(symN));
  }
  //--------------------------------------------------------------
  std::string SpaceGroup::formatAllSymops_as_xyz() const
  //!< real-space symops
  {
    std::string s;
    for (int i=0;i<num_symops();++i) {
      s += symop(i).format()+"\n";
    }
    return s;
  }
  //--------------------------------------------------------------
  std::string SpaceGroup::formatAllSymops_as_hkl() const
  //!< reciprocal symops
  {
    std::string s;
    for (int i=0;i<num_symops();++i) {
      s += MVutil::FormatSymop_as_hkl(symop(i), "[]")+"\n";
    }
    return s;
  }
  //--------------------------------------------------------------
  std::string SpaceGroup::formatISYM_as_hkl() const
  //!< reciprocal symops as ISYM table
  {
    clipper::RTop_frac rtf(clipper::String("-x,-y,-z"));
    clipper::Symop friedel(rtf);

    std::string s = "   Reciprocal space symmetry operators:\n";
    s += "    original indices for unique reflection hkl with symmetry number ISYM\n";
    for (int j=0;j<2;++j) { // loop +, -
      if (j == 0) {
        s +=
          "\n                              Bijvoet positive \n";
      } else {
        s +=
          "\n                              Bijvoet negative \n";
      }
      int nh = std::min(Nsymp,4);
      for (int k=0;k<nh;++k) {
        s += "   ISYM";
        if (k < nh-1) {s += "              ";}
      }
      s += "\n";
      for (int i=0;i<Nsymp;++i) {    // loop primitive operators
        int isym = 2*i + 1;
        clipper::Symop sop = invrotsymops[i];
        if (j == 1) { // Bijvoet negative
          sop = clipper::Symop(friedel * sop);
          isym += 1;
        }
        s += "  " + StringUtil::itos(isym, 4) + " " +
          StringUtil::LeftString(MVutil::FormatSymop_as_hkl(sop, "[]"), 14);
        if (i < Nsymp-1 && (i+1)%4 == 0) {
          s += "\n";
        }
      }
      s += "\n";
    }
    return s;
  }
  //--------------------------------------------------------------
  void SpaceGroup::SetLatType()
  {
    //  xHM symbol, colons removed (ie R 3 :H  converted to H 3)
    //  but allow for spacegroup names with ":1" etc in them
    // Take either the first character or the one after a colon
    lattype = 'P';
    char aftercolon = ' ';
    bool first = true;
    for (size_t i=0; i<spacegroupname.size(); i++) {
      if (first) {
        if (spacegroupname[i] != ' ') {
          lattype = spacegroupname[i];
          first = false;
        }
      }
      if (spacegroupname[i] == ':') {
        aftercolon = spacegroupname[i+1];
      }
    }
    if (AllowedLatticeType(aftercolon)) {
      lattype = aftercolon;
    }
    if (lattype == 'R' || lattype == 'H') {
      // Rhombohedral, check symmetry operators
      clipper::Symop_code zxy(clipper::Symop(RTop<>
                 (Mat33<double>(0,0,1.0, 1.0,0,0,  0,1.0,0)))); // operator z,x,y
      bool found = false;
      for (size_t i=0;i<csymops.size();++i) {
        if (clipper::Symop_code(csymops[i]) == zxy) {
          found = true;
          break;
        }
      }
      if (found) lattype = 'R';
      else lattype = 'H';
    }
  }
  //--------------------------------------------------------------
  SpaceGroup SpaceGroup::PattersonGroup() const
  {
    return SpaceGroup(clipper::Spacegroup
      (clipper::Spgr_descr(generator_ops().patterson_ops())));
  }
  //--------------------------------------------------------------
  SpaceGroup SpaceGroup::PointGroup() const
  {
    return SpaceGroup(clipper::Spacegroup
      (clipper::Spgr_descr(generator_ops().pgrp_ops())));
  }
  //--------------------------------------------------------------
  CSym::ccp4_symop SpaceGroup::MakeCCP4symop(const clipper::Symop& S) const
  // CCP4 symop from clipper symop
  {
    CSym::ccp4_symop op;
    for (int k = 0; k < 3; ++k) {
      for (int l = 0; l < 3; ++l) {
        op.rot[k][l] = S.rot()(k,l);
      }
      op.trn[k] = S.trn()[k];
    }
    return op;
  }
  //--------------------------------------------------------------
  clipper::Symop SpaceGroup::MakeClippersymop(const CSym::ccp4_symop& op) const
  // Clipper symop from CCP4 symop
  {
    clipper::Mat33<double> rot;
    clipper::Vec3<double>  trn;

    for (int k = 0; k < 3; ++k) {
      for (int l = 0; l < 3; ++l) {
        rot(k,l) = op.rot[k][l];
      }
      trn[k] = op.trn[k];
    }
    return clipper::Symop(RTop<>(rot,trn));
  }
  //--------------------------------------------------------------
  std::string SpaceGroup::CCP4spaceGroupName() const
  // Use CCP4 library routine to get name from operators csymops
  {
    int nsym = num_symops();
    std::vector<CSym::ccp4_symop> ccp4ops(nsym);
    for (int i=0;i<nsym;++i) {
      //^      std::cout << "Symop " << csymops[i].format() <<"\n"; //^-
      ccp4ops[i] = MakeCCP4symop(csymops[i]);
    }
    CSym::CCP4SPG* ccp4spg = CSym::ccp4_spgrp_reverse_lookup(nsym, &*ccp4ops.begin());
    if (ccp4spg == NULL) return "";
    //^    std::cout << "CCP4 spgname: " << std::string(ccp4spg->symbol_xHM) <<"\n"; //^
    return std::string(ccp4spg->symbol_xHM);
  }
  //--------------------------------------------------------------
  void SpaceGroup::CCP4spaceGroupNumber()
  // Use CCP4 library routine to get numbers from operators
  {
    int nsym = num_symops();
    std::vector<CSym::ccp4_symop> ccp4ops(nsym);
    for (int i=0;i<nsym;++i) {
      ccp4ops[i] = MakeCCP4symop(symop(i));
    }
    CSym::CCP4SPG* ccp4spg = CSym::ccp4_spgrp_reverse_lookup(nsym, &*ccp4ops.begin());
    if (ccp4spg == NULL) {
      return;
    }
    spacegroupnumber = ccp4spg->spg_num;
    CCP4spacegroupnumber = ccp4spg->spg_ccp4_num;
  }
  //--------------------------------------------------------------
  std::string SpaceGroup::AllSymopsfromCCP4(const CSym::CCP4SPG* ccp4sg) const
  // Make symop string from CCP4 space group
  {
    std::string s;
    int nsym = ccp4sg->nsymop;
    for (int i=0;i<nsym;++i) {
      s += MakeClippersymop((ccp4sg->symop)[i]).format();
      if (i < nsym-1) s += ";";
    }
    //^    std::cout << "SpaceGroup::AllSymopsfromCCP4:" << s << "\n"; //^
    return s;
  }
  //--------------------------------------------------------------
  bool operator == (const SpaceGroup& a,const SpaceGroup& b)
  {
    if (a != b) return false;
    return true;
  }
  //--------------------------------------------------------------
  bool operator != (const SpaceGroup& a,const SpaceGroup& b)
  {
    if (a.num_symops() != b.num_symops()) return true;
    if (a.spacegroup_number() != b.spacegroup_number()) return true;
    if (a.SymopCodes() != b.SymopCodes()) return true;
    return false;
  }
  // =============================================================
  //--------------------------------------------------------------
  void SymElement::print() const
  {
    std::cout << "Symmetry element " << Nfold
              << "-fold, comprises symops : ";
    for (size_t i=0;i<symops.size();i++) std::cout << symops[i]+1 << " ";
    std::cout << "\n";
  }
  //--------------------------------------------------------------
  // Constructors
  hkl_symmetry::hkl_symmetry()  {
    Nsymp = 0;  // Null initialisation
    spgname = "";
  }
  //--------------------------------------------------------------
  hkl_symmetry::hkl_symmetry(const std::string SpgName)
  //            ***********
  {
    spgname = SpgName;
    spaceGroup.init(spgname);
    set_symmetry();
  }
  //--------------------------------------------------------------
  hkl_symmetry::hkl_symmetry(const clipper::Spacegroup& ClpSG)
  //              ***********
  //    Construct from clipper spacegroup
  {
    spgname = ClpSG.symbol_hm();
    spaceGroup.init(ClpSG);
    set_symmetry();
  }
  //--------------------------------------------------------------
  hkl_symmetry::hkl_symmetry(const SpaceGroup& SG)
  //              ***********
  //    Construct from scala spacegroup
  {
    spgname = SG.symbol_xHM();
    spaceGroup = SG;
    set_symmetry();
  }
  //--------------------------------------------------------------
  hkl_symmetry::hkl_symmetry(const int& SpgNumber)
  //              ***********
  {
    spgname = "";
    spaceGroup.init(SpgNumber);
    set_symmetry();
  }
  //--------------------------------------------------------------
  void hkl_symmetry::set_symmetry()
  // Set up object
  {
    clipper::Symop S1S2inv, rot1, rot2, rot;

    Nsymp = spaceGroup.num_primops();
    symelmt.clear();
    elements.clear();

    // Extract lattice centre symbol
    LatType = spaceGroup.LatType();
    //^    std::cout << "spacegroup_number " << spaceGroup.spacegroup_number() << "\n";
    cryssys = CrysSysfromSGnumber(spaceGroup.Spacegroup_number());  // Crystal system

    centro = (spaceGroup.num_inversion_symops() > 1);
    chiral = ChiralTest();

    //    std::cout << "LatType " << LatType
    //        << " system " << CrystalType(cryssys,LatType).format() <<"\n";
    //    std::cout << "num_inversion_symops " << spaceGroup.num_inversion_symops() <<"\n";
    //    for (int i=0;i<spaceGroup.num_inversion_symops();++i){
    //      std::cout << i+1 <<"\n"<<spaceGroup.inversion_symop(i).format() <<"\n";
    //    }
    //    std::cout <<"Chiral: "<< chiral<<"\n";

    // Generate symcodes for all operators & inverses without translations
    inv_rot_symcodes.resize(Nsymp);
    rot_symcodes.resize(Nsymp);
    for (int i = 0; i < Nsymp; i++) {
      inv_rot_symcodes[i] = spaceGroup.InvRotSymopCode(i);
      rot_symcodes[i] = spaceGroup.RotSymopCode(i);
    }

    // Generate list of symmetry elements and associated symops
    // Generate list of elements by repeatedly applying each operator to
    // itself
    for (int i = 0; i < Nsymp; i++) {
      rot1 = spaceGroup.InvRotSymop(i);
      rot = rot1;
      if (i == 0) {
        if (!spaceGroup.IsSymopIdentity(0)) {
          // First symmetry operator is not identity, fatal
          ReportErrors::printFatalError(
            "hkl_symmetry: Illegal spacegroup, identity not first");
        }
      }
      // Start new group
      int k = 1;
      bool found = false;
      SymElement elmt;
      // First symop of element, store
      elmt.symops.push_back(i);

      if (i == 0) {
        // 1st operator = I
        found = true;
      } else {
        while (true) {
          k++;
          rot = clipper::Symop(rot * rot1);
          if (IsSymopIdent(rot)) {
            //
            //^std::cout << k <<"-fold rotation\n"
            //^ << rot1.format() << "\n";;

            // Found k-fold rotation
            found = true;
            break;  // from while(true) loop
          }
          //  find which symmetry operator this corresponds to
          found = false;
          clipper::Symop_code jcode(rot);
          for (int j = 0; j < Nsymp; j++) {
            if (j != i && spaceGroup.InvRotSymopCode(j) == jcode) {
              elmt.symops.push_back(j);
              found = true;
              break;
            }
          }
          if (!found) {
            // Incomplete symmetry element, already got this one
            break;}
        }
      }
      if (found) {
        // rot**k is identity, ie k-fold rotation axis
        // Store symmetry element incuding identity (k=1)
        elmt.Nfold = k;
        // Have we already go this one?
        bool known = false;
        if (elements.size() > 0) {
          for (size_t j=0;j<elements.size();j++)
            if (ElementEqual(elmt, elements[j])) {
              known = true;
              break;
            }
        }
        if (!known) {
          elmt.iaxis = AxisDirection(elmt);
          elements.push_back(elmt);
          // std::cout << "Element:  " <<  elmt.Nfold << " -fold axis. Symops ";
          // for (int jj=0; jj< elmt.symops.size();jj++)
          //  std::cout << elmt.symops[jj] << "  ";
          // std::cout << "\n";
        }
      }
    }

    // Sort symmetry elements by rotation order
    std::sort(elements.begin(), elements.end());

    Nelement_ = elements.size();

    // Make list of which symmetry element each symop belongs to
    // Accept only the first occurance, so that operators which
    // belong to more than one element belong only to the lower
    // rotation order
    // ie from these indices eg a 4-fold element will not contain
    // the 2-fold operator
    // Delete occurances after first from element list
    element_index.resize(Nsymp);
    for (int i=0;i<Nsymp;i++) element_index[i] = -1;
    // Loop elements
    for (size_t j=0;j<elements.size();j++) {
      std::vector<int> ops; // for copying symops with some perhaps deleted
      // Loop symops for element j
      for (size_t i=0;i<elements[j].symops.size();i++) {
        int k = elements[j].symops[i];  // symop number
        if (element_index[k] < 0) {
          element_index[k] = j;
          // save this symop
          ops.push_back(elements[j].symops[i]);
        }
      }
      // Replace symops for element
      elements[j].symops = ops;
    }
    // Check that all symops are assigned to an element
    for (int i=0;i<Nsymp;i++) {
      if (element_index[i] < 0)
ReportErrors::printFatalError("hkl_symmetry: symmetry operator not in element");
    }
    // Make lookup-table of ISYM/2 pairs pointing to corresponding
    // symmetry operator
    // Symmetry operator Jsym (= (Isym-1)/2, range 0:Nsymp-1)
    //  is the operator used to reduced original indices h
    //  to asymmetric unit
    // h(asu) = h(orig) S(Jsym)   where S is spacegroup.Symop
    //
    // Generate complete matrix (even though it is symmetric)
    // to simplify indexing
    for (int i1 = 0; i1 < Nsymp; i1++) {
      rot1 = spaceGroup.RotSymop(i1);
      //^
      //        std::cout << "\n****\n Op1 " << i1 << " "
      //                  << rot.format()
      //                  << "\n\n";
      //^-
      for (int i2 = 0; i2 < Nsymp; i2++) {
        //  rot = S(i1) S(i2)**-1
        rot2 = spaceGroup.InvRotSymop(i2);   // inverse operator
        S1S2inv =  clipper::Symop(rot1 * rot2);
        // S1S2inv should be an reciprocal-space operator in symop
        // Test all of them using symcodes
        clipper::Symop_code code(S1S2inv);
        int jsym = -1;
        for (int i = 0; i < Nsymp; i++) {
          if (code == spaceGroup.RotSymopCode(i)) {
            jsym = i;
            break;
          }
        }
        if (jsym < 0) {
          ReportErrors::printFatalError("hkl_symmetry - operator not in group");
        }
        //^+
        //          std::cout << "  Op2 " << i2 << " " << InvSymOpFormat(spacegroup_.symop[i2])
        //                    << " result " << jsym << " " << InvSymOpFormat(spacegroup_.symop[jsym])
        //                    << " element " << element_index[jsym]
        //                    << "\n";
        //^-
        // Entry for i1, i2, symmetry element number corresponding
        // to operator S(i1) S(i2)**-1
        symelmt.push_back(element_index[jsym]);
      }
    }
    if (spgname == "") spgname = symbol_xHM();
  }
  //--------------------------------------------------------------
  void hkl_symmetry::print_element(const int& idx) const
  //                 ^^^^^^^^^^^^
  {
    std::cout << "Element number " << idx+1 << " : ";
    elements[idx].print();
    //    std::cout << "\n";
    for (int i=0;i<Nsymp;i++) {
      std::cout << "Symop " << i+1 << "  "
                << spaceGroup.InvRotSymop(i).format()
                << "  belongs to element "
                << element_index[i]+1 << "\n";
    }
  }
  //--------------------------------------------------------------
  void hkl_symmetry::print_elements() const
  //                  ^^^^^^^^^^^^
  {
    std::cout << "\nSymmetry elements: number = " <<  elements.size()
              << "\n";
    for (size_t j=0;j<elements.size();j++) {
      print_element(j);
    }
  }
  //--------------------------------------------------------------
  clipper::Vec3<int>  hkl_symmetry::AxisDirection(const SymElement& elmt) const
  // Axis direction for kelement'th symmetry element
  {
    clipper::Mat33<double> RecSymOpM = spaceGroup.InvRotSymop(elmt.symops[0]).rot();
    clipper::Vec3<double> axis = MVutil::AxisDirection(RecSymOpM);
    // Make integral version of vector
    return MVutil::IntVec(axis);
  }
  //--------------------------------------------------------------
  clipper::Vec3<int>  hkl_symmetry::AxisDirection(const clipper::Symop& csymop)
  // Axis direction for symop
  {
    clipper::Mat33<double> RecSymOpM = csymop.rot();
    clipper::Vec3<double> axis = MVutil::AxisDirection(RecSymOpM);
    // Make integral version of vector
    return MVutil::IntVec(axis);
  }
  //--------------------------------------------------------------
  std::string hkl_symmetry::AxisString(const clipper::Symop& csymop)
  {
    std::string s = " ";
    clipper::Vec3<int> Iaxis = AxisDirection(csymop);
    // Along principle axis?
    if ((abs(Iaxis[0])+abs(Iaxis[1])+abs(Iaxis[2])) == 1) {
      if (Iaxis[0] == 1) s = "h";
      if (Iaxis[1] == 1) s = "k";
      if (Iaxis[2] == 1) s = "l";
    }
    return s;
  }
  //--------------------------------------------------------------
    std::string hkl_symmetry::AxisString(const int& kelement,
                                         clipper::Vec3<int>& Iaxis) const
  {
    std::string s = " ";
    // Get axis direction from sample (first) symop belonging
    // to this element
    Iaxis = AxisDirection(elements[kelement]);
    // Along principle axis?
    if ((abs(Iaxis[0])+abs(Iaxis[1])+abs(Iaxis[2])) == 1) {
      if (Iaxis[0] == 1) s = "h";
      if (Iaxis[1] == 1) s = "k";
      if (Iaxis[2] == 1) s = "l";
    }
    return s;
  }
  //--------------------------------------------------------------
  std::string hkl_symmetry::format_element(const int& kelement) const
    // make formatted version of symmetry element
    //  kelement (0 - Nelement-1)
    //  "N-fold axis along x|y|z|(lmn)"
  {
    if (kelement < 0 || kelement > Nelement_) {
        return "";
    }
    if (elements[kelement].Nfold == 1)
      return "identity";
    // Axis order
    clipper::String  s = clipper::String(elements[kelement].Nfold,1)
                       +"-fold ";

    // Get axis direction from sample (first) symop belonging
    // to this element
    clipper::Vec3<int> Iaxis;
    s += AxisString(kelement, Iaxis)+" ";
    s = s+"("+clipper::String(Iaxis[0], 2)+clipper::String(Iaxis[1], 2)
      +clipper::String(Iaxis[2], 2)+") ";

    // Symmetry operators (only those unique to this element)
    for (size_t i=0;i<elements[kelement].symops.size();i++) {
      clipper::Symop RecSymOp = spaceGroup.InvRotSymop(elements[kelement].symops[i]);
      s += MVutil::FormatSymop_as_hkl(RecSymOp,"{}");
    }

    return std::string(s);
  }
  //--------------------------------------------------------------
  std::string hkl_symmetry::XML_element(const int& kelement) const
    // make XMLormatted version of symmetry element
    //  kelement (0 - Nelement-1)
    //  <RotationOrder> Nfold </RotationOrder>
    //  <Axis> a b c </Axis>
  {
    if (kelement < 0 || kelement > Nelement_) {
      return "";
    }

    // Axis order
    clipper::String  s = "<RotationOrder> "+
      clipper::String(elements[kelement].Nfold,2)+" </RotationOrder>";

    clipper::Vec3<int> Iaxis(0,0,0);
    if (elements[kelement].Nfold > 1) {
      // Get axis direction from sample (first) symop belonging
      // to this element
      clipper::Mat33<double> RecSymOp = spaceGroup.InvRotSymop(elements[kelement].symops[0]).rot();
      clipper::Vec3<double> axis = MVutil::AxisDirection(RecSymOp);
      // Make integral version of vector
      Iaxis = MVutil::IntVec(axis);
    }

    s = s+" <Axis>"+clipper::String(Iaxis[0], 2)+clipper::String(Iaxis[1], 2)
      +clipper::String(Iaxis[2], 2)+" </Axis>\n";

    // String format
    s += StringUtil::MakeXMLtag("SymmetryElementString", format_element(kelement));

    return std::string(s);
  }
  //--------------------------------------------------------------
  Hkl hkl_symmetry::put_in_asu(const Hkl& hkl, int& isym) const
    //            ^^^^^^^^^
  {
    return Hkl(spaceGroup.put_in_asu(hkl.HKL(), isym));
  }

  //--------------------------------------------------------------
  Hkl hkl_symmetry::get_from_asu(const Hkl& hkl, const int& isym) const
    //            ^^^^^^^^^^^
  {
    return Hkl(spaceGroup.get_from_asu(hkl.HKL(), isym));
  }
  //--------------------------------------------------------------
  bool hkl_symmetry::is_centric(const Hkl& hkl) const
  {
    return spaceGroup.hkl_class(hkl.HKL()).centric();
  }
  //--------------------------------------------------------------
  float hkl_symmetry::epsilon(const Hkl& hkl) const
  {
    return spaceGroup.hkl_class(hkl.HKL()).epsilon();
  }
  //--------------------------------------------------------------
  bool hkl_symmetry::equals_r(const hkl_symmetry& other) const
    // true if spacegroup other is same ignoring translations
  {
    if (Nsymp != other.Nsymp) return false;
    if (LatType != other.LatType) return false;
    // compare symcodes corresponding to rotation-only parts
    // for all symops
    for (int i=0; i < Nsymp; i++) {
      bool found = false;
      for (int j=0; j<Nsymp; j++) {
        if (inv_rot_symcodes[i] == other.inv_rot_symcodes[j]) {
          found = true;
          break;
        }
      }
      if (! found) return false;
    }
    return true;
  }
  //--------------------------------------------------------------
  bool hkl_symmetry::equals_r_order(const hkl_symmetry& other) const
    // true if spacegroup other is same ignoring translations
    // with operators in the same order
  {
    if (Nsymp != other.Nsymp) return false;
    if (LatType != other.LatType) return false;
    // compare symcodes corresponding to rotation-only parts
    // for all symops in same order
    bool found = true;
    for (int i=0; i < Nsymp; i++) {
      if (inv_rot_symcodes[i] != other.inv_rot_symcodes[i]) {
        found = false;
      }
    }
    return found;
  }
  //--------------------------------------------------------------
  bool hkl_symmetry::equals_rt(const hkl_symmetry& other) const
  {
    if (Nsymp != other.Nsymp) return false;
    if (LatType != other.LatType) return false;
    if (spaceGroup != other.spaceGroup) return false;
    return true;
  }
  //--------------------------------------------------------------
  bool hkl_symmetry::IsElementIdent(const int& kelement) const
  {
    return (elements.at(kelement).Nfold == 1) ? true : false;
  }
  //--------------------------------------------------------------
  int hkl_symmetry::NopInElement(const int& kelement) const
    // Return number of operators in kelement'th symmetry element
  {
    if (kelement < 0 || kelement > Nelement_) {
      return -1;
    }
    return elements[kelement].symops.size();
  }
  //--------------------------------------------------------------
  std::vector<double> hkl_symmetry::SymopInElement(const int& lsym,
                                              const int& kelement) const
    // Return rotation part of lsym'th inverse symmetry operator of
    // kelement'th symmetry element
  {
    return MVutil::SetVMat33(spaceGroup.InvRotSymop(elements[kelement].symops[lsym]).rot());
  }
  //--------------------------------------------------------------
  clipper::Symop hkl_symmetry::ClipperSymopInElement(const int& lsym,
                                              const int& kelement) const
  // Return lsym'th inverse symmetry operator of  kelement'th symmetry element
  {
    return spaceGroup.InvRotSymop(elements[kelement].symops[lsym]);
  }
  //--------------------------------------------------------------
  std::vector<clipper::Symop> hkl_symmetry::RotSymopsInElement(const int& kelement) const
  {
    std::vector<clipper::Symop> ops(elements.at(kelement).symops.size());
    for (size_t i=0;i<elements[kelement].symops.size();++i) {
      ops[i] = spaceGroup.RotSymop(elements[kelement].symops[i]);
    }
    return ops;
  }
  //--------------------------------------------------------------
  std::vector<clipper::Symop> hkl_symmetry::PrimRotSymops() const
  {
    std::vector<clipper::Symop> ops(spaceGroup.num_primops());
    for (size_t i=0;i<ops.size();++i) {
      ops[i] = spaceGroup.RotSymop(i);
    }
    return ops;
  }
  //--------------------------------------------------------------
  int hkl_symmetry::get_symelmt(const int& isym1, const int& isym2) const
  // Return symmetry element (from 0) corresponding to the relation
  // between observations with isym1 and isym2
  // Isym = 2 * SymNumber + 1  for h+
  // Isym = 2 * SymNumber + 2  for h-
  // Ignore difference between h+ and h- for this purpose
  {
    int i1 = (isym1-1)/2;
    int i2 = (isym2-1)/2;
    return symelmt[i2*Nsymp + i1];
  }
  //--------------------------------------------------------------
  std::vector<int> hkl_symmetry::get_symelmt(const Hkl& h1, const Hkl& h2) const
  // Return list of symmetry elements (from 0)
  // corresponding to the relation
  // between observations with original indices h1 & h2
  // This call is used for centric reflections where more than one
  // symmetry operator relates two observations
  // Ignore difference between h+ and h- for this purpose
  {
    std::vector<int> listelmt;
    clipper::HKL h = h1.HKL();
    Hkl hh;

    // Loop Isym for 2*Nsymp operators
    for (int isym=1;isym<=Nsymp*2;isym++) {
      if (isym%2 != 0) {
        hh = Hkl(h.transform(spaceGroup.InvRotSymop((isym-1)/2)));
      } else {
        hh = Hkl((-h).transform(spaceGroup.InvRotSymop((isym-1)/2)));
      }
      if (hh == h2) {
        // Found operator which relates h1 to h2
        // store corresponding symmetry element
        listelmt.push_back(element_index[(isym-1)/2]);
      }
    }
    return listelmt;
  }
  //--------------------------------------------------------------
  std::string hkl_symmetry::symbol_Hall() const
  {
    return spaceGroup.symbol_Hall();
  }
  //--------------------------------------------------------------
  std::string hkl_symmetry::symbol_xHM(const char& HorR) const
  {
    return SGnameHtoR(spaceGroup.symbol_xHM(), HorR);
  }
  //--------------------------------------------------------------
  //! return formatted crystal system
  std::string hkl_symmetry::formatCrysSys() const
  {
    return CrystalType(cryssys, LatType).format(false);
  }
  //--------------------------------------------------------------
  std::vector<int> hkl_symmetry::CellConstraint() const
  {
    std::vector<int> lbcell(6);
    // cell flags
    //   [1]-[7] triclinic->cubic
    //   [8] rhombohedral
    //   [9] default
    int lb[9][6] =
      {
        {-1,-1,-1,-1,-1,-1},
        {-1,-1,-1, 0,-1, 0},
        {-1,-1,-1, 0, 0, 0},
        {-1, 1,-1, 0, 0, 0},
        {-1, 1,-1, 0, 0, 0},
        {-1, 1,-1, 0, 0, 0},
        {-1, 1, 1, 0, 0, 0},
        {-1, 1, 1,-1, 4, 4},
        {-1,-1,-1,-1,-1,-1}
      };

    int ic = -1;
    if (cryssys == TRICLINIC)
      {ic = 0;}
    else if (cryssys == MONOCLINIC)
      {ic = 1;}
    else if (cryssys == ORTHORHOMBIC)
      {ic = 2;}
    else if (cryssys == TETRAGONAL)
      {ic = 3;}
    else if (cryssys == TRIGONAL)
      {ic = 4;
        if (LatType == 'R') ic = 7;
      }
    else if (cryssys == HEXAGONAL)
      {ic = 5;}
    else if (cryssys == CUBIC)
      {ic = 6;}

    for (int i=0;i<6;i++) {lbcell[i] = lb[ic][i];}
    return lbcell;
  }
  //--------------------------------------------------------------
  bool hkl_symmetry::LatticePresent(const Hkl& hkl) const
  // true if present in lattice,
  // false if systematically absent from lattice
  {
    if (LatType == 'P' || LatType == 'R') return true;
    if (LatType == 'I') {
      return ((hkl.h()+hkl.k()+hkl.l())%2 == 0);
    } else if (LatType == 'C') {
      return ((hkl.h()+hkl.k())%2 == 0);
    } else if (LatType == 'A') {
      return ((hkl.k()+hkl.l())%2 == 0);
    } else if (LatType == 'B') {
      return ((hkl.h()+hkl.l())%2 == 0);
    } else if (LatType == 'F') {
      int test = abs(hkl.h()%2) + abs(hkl.k()%2) + abs(hkl.l()%2);
      return (test == 0) || (test == 3);
    } else if (LatType == 'H') {
      return ((-hkl.h()+hkl.k()+hkl.l())%3 == 0);
    } else {
      ReportErrors::printFatalError
        (clipper::String("hkl_symmetry::LatticePresent: unrecognised lattice ")
         +LatType);
    }
    return false;  // dummy, never gets here
  }
  //--------------------------------------------------------------
  int hkl_symmetry::ChangeBasis(const scala::ReindexOp& reindex)
  {
    int OK = spaceGroup.ChangeBasis(reindex);
    if (OK < 0) {return OK;}
    set_symmetry();
    return OK;
  }
  //--------------------------------------------------------------
  // true if jel'th element of this symmetry object == jel2'th element of other
  // for any of the symops in either element
  bool hkl_symmetry::equalElement(const int& jel,
          const hkl_symmetry& other, const int& jel2) const
  {
    int n1 = NopInElement(jel);
    int n2 = other.NopInElement(jel2);
    if (n1 != n2) {
      return false;  // not same number of operators
    }
    bool same = true;
    for (int i1=0;i1<n1;++i1) {
      clipper::Symop_code c1 = inv_rot_symcodes[elements[jel].symops[i1]];
      for (int i2=0;i2<n2;++i2) {
        clipper::Symop_code c2 = inv_rot_symcodes[elements[jel2].symops[i2]];
        //^
        //      std::cout <<"Op1:\n"<<
        //        spaceGroup.InvRotSymop(elements[jel].symops[i1]).format() <<"\n"
        //                <<" code " << c1 <<"\nOp2:\n" <<
        //        spaceGroup.InvRotSymop(elements[jel2].symops[i2]).format() <<"\n"
        //                <<" code " << c2 <<"\n";
        //^-
        if (c1 != c2) {same = false; break;}
      }
    }
    return same;
  }
  //--------------------------------------------------------------
  bool hkl_symmetry::ElementEqual(const SymElement& a, const SymElement& b) const
  {
    if (a.Nfold != b.Nfold) return false;
    if (a.symops.size() != b.symops.size()) return false;
    for (size_t i=0;i<a.symops.size();i++) {
      bool found = false;
      for (size_t j=0;j<b.symops.size();j++)
        if (a.symops[i] == b.symops[j]) {
          found = true;
          break;
        }
      if (!found) return false;
    }
    return true;
  }
  //--------------------------------------------------------------
  CrystalSystem hkl_symmetry::CrysSysfromSGnumber(const int& SpaceGroupNumber)
  {
    // based on true space group number
    // International Tables spacegroup number
    //   195-230 cubic
    //   168-194 hexagonal
    //   143-167 trigonal
    //   75-142  tetragonal
    //   16-74   orthorhombic
    //   3-15    monoclinic
    //   1-2     triclinic (and 0)
    CrystalSystem CrysSys;
    if (SpaceGroupNumber < 0 || SpaceGroupNumber > 230)
      {ReportErrors::printFatalError
                        ("Illegal space group number "+
                         clipper::String(SpaceGroupNumber));}
    if (SpaceGroupNumber > 194)
      {CrysSys = CUBIC;}
    else if (SpaceGroupNumber > 167)
      {CrysSys = HEXAGONAL;}
    else if (SpaceGroupNumber > 142)
      {CrysSys = TRIGONAL;}
    else if (SpaceGroupNumber > 74)
      {CrysSys = TETRAGONAL;}
    else if (SpaceGroupNumber > 15)
      {CrysSys = ORTHORHOMBIC;}
    else if (SpaceGroupNumber > 2)
      {CrysSys = MONOCLINIC;}
    else
      {CrysSys = TRICLINIC;}
    return CrysSys;
  }
  //--------------------------------------------------------------
  bool hkl_symmetry::ChiralTest()
  // return true if space group is chiral
  // ie no symops invert, tested by negative determinant
  {
    for (int l=0;l<Nsymp;l++) {
      if(spaceGroup.symop(l).rot().det() < 0.0) return false;
    }
    return true;
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  CheckReindexSymmetry::CheckReindexSymmetry(const ReindexOp& reindex,
                                             const bool& sgnamegiven,
                                             const std::string& spaceGroupName,
                                             const std::string& filespacegroupname,
                                             const hkl_symmetry& symmetry)
  // Check whether "reindex" on spacegroup filespacegroupname
  //
  //  reindex            reindex operator to be applied
  //  sgnamegiven        true if desired SG name comes from user input or reference
  //  spaceGroupName     desired SG name
  //  filespacegroupname SG name from input file
  //  symmetry           current symmetry
  //
  // For merged files: it is only legitimate to reindex within
  // the same point group
  //
  // Cases:
  // 1) Just explicit reindex (no spacegroup from input or reference)
  //    check that reindexed filespacegroupname is a legimate space group
  //    (status = -1 if not)
  // 2) Just explicit space group name: mainly relevant for C2<->I2, R3:H<->R3:R (may not work).
  //    and primitive orthorhombics
  //    A reindex operator will be generated as appropriate
  // 3) Both reindex and space group given: check that they are compatible and fail if not
  //    (status = -2 if intensity group (point group) is changed
  //
  {
    status = 0;  // default OK flag
    newsymmetry_ = symmetry;
    reindex_ = reindex;
    sgnamegiven_ = sgnamegiven;
    reindexfrominput_ = reindex.fromInput();
    spacegroupname_ = spaceGroupName;
    samereferencegroup_ = false;
    det_ = reindex_.det();

    // Special case R3 <-> H3
    if (checkRhombohedral(filespacegroupname)) {
      return;
    }

    if (!sgnamegiven_) { // space group name not given, so derive from MTZ name + reindex
      if (newsymmetry_.ChangeBasis(reindex_) != 0) {
        // An issue or failure in change basis with this operator,
        //   set status, = -1 fail, > 0 warning
        status = newsymmetry_.GetSpaceGroup().Status();
      }
      SG_ = newsymmetry_.GetSpaceGroup();
      spacegroupname_ = newsymmetry_.symbol_xHM();
      // check reindexed SG is same pointgroup
      SpacegroupReindexOp sgreindex(filespacegroupname, spacegroupname_);
      if (sgreindex.sameReferenceGroup()) {
        samereferencegroup_ = true;   // same reference group
      }
    } else {
      // Space group name given
      if (reindexfrominput_) {
        // reindex has come from reference, so let's trust it and SG name
        status = 0;
      } else {
        // if reindex is given as well as name, check compatibility
        SG_.init(filespacegroupname);
        if (SG_.ChangeBasis(reindex_) < 0) {
          // fail
          status = SG_.Status();
        } else {
          // check reindexed SG is same pointgroup
          SpacegroupReindexOp sgreindex(SG_.symbol_xHM(), spacegroupname_);
          if (sgreindex.sameReferenceGroup()) {
            samereferencegroup_ = true;   // same reference group
          } else if (!sgreindex.sameIntensityGroup()) {
            // not same
            status = -2;
          }
        }
      }
    }
    // SG changed if SG(filespacegroupname) != SG_
    spacegroupchanged_ = false;
    if (SpaceGroup(filespacegroupname) != SG_) {
      spacegroupchanged_ = true;
    }
    sameorigin_ = true;
    if (status == -3) {
      //  the origin might have changed, space group undetermined
      sameorigin_ = false;
    }
  }
  //--------------------------------------------------------------
  void CheckReindexSymmetry::report(const bool& verbose,
                                    phaser_io::Output& output) const
  {
    output.logTab(0,LXML, reindex_.as_hkl_XML()+"\n");
    output.logTab(0,LXML, reindex_.as_XML());

    std::string tf = "False";
    if (spacegroupchanged_) {tf = "True";}
    output.logTab(0,LXML,
                    StringUtil::MakeXMLtag("SpacegroupChanged",tf));

    tf = "True";
    if (sameorigin_) {tf = "False";}
    output.logTab(0,LXML,
                    StringUtil::MakeXMLtag("OriginChanged",tf));

    std::string statusmessage = SG_.statusMessage();
    std::string newspacegroupname = SG_.symbol_xHM();
    if (newspacegroupname == "Unknown") {newspacegroupname = spacegroupname_;}

    output.logTab(0,LXML,
                  StringUtil::MakeXMLtag("NewSpacegroupName",newspacegroupname));

    // Normally reindex determinant = 1.0, else changes cell volume
    const double TOL = 0.0001;
    if (!Close(det_, 1.0, TOL)) {
      std::string message = "Reindex operator changes the lattice volume";
      if (det_ < 0.99) {
        message += ", some reflections will be discarded";
      } else {
        message += ", some reflections will be missing in new lattice";
      }
      message += ", det(H) = "+StringUtil::ftos(det_,6,2);
      ReportErrors::printWarning(message, "ReindexDeterminantMessage", false);
    }
    if ((status == -3) && (sgnamegiven_)) {
      statusmessage =
        "Specified space group " + spacegroupname_ +
        " used, no check has been done on its validity";
      ReportErrors::printWarning(statusmessage, "ChangedSpaceGroupMessage", false);
    } else if (status < 0) {
      // Fatal, maybe die here
      if (status == -2) {
        statusmessage =
          "Specified SPACEGROUP "+spacegroupname_+
          " must belong to same crystal system and point group\n"+
          "as the reindexed input space group "+newspacegroupname+
          ", reindex operator "+reindex_.as_hkl();
      } else if (status == -3) {
        statusmessage +=
          "\nYou need to give an explicit space group to override the input one";
      }
      ReportErrors::printFatalError(statusmessage);
    } else if (status != 0) {
      ReportErrors::printWarning(statusmessage, "ReindexWarning", false);
    }
    if (!sgnamegiven_) {
      if (verbose) {
        output.logTabPrintf(1,LOGFILE,
                            "Reindexed space group : %s\n\n",
                            newsymmetry_.symbol_xHM().c_str());
      }
    } else {
      if (verbose) {
        output.logTabPrintf(1,LOGFILE,
                            "Space group : %s\n\n",
                            newsymmetry_.symbol_xHM().c_str());
      }
    }

    std::string s =
      "Copying merged data to output file in space group " + spacegroupname_;
    output.logTab(0,LXML,
                  StringUtil::MakeXMLtag("CopyMessage", s));
  }
  //--------------------------------------------------------------
  bool CheckReindexSymmetry::HorR(const char& X) const
  // true if 'X' == 'H' or 'R'
  {
    if (X == 'H') {return true;}
    if (X == 'R') {return true;}
    return false;
  }
  //--------------------------------------------------------------
  bool CheckReindexSymmetry::checkRhombohedral
  (const std::string filespacegroupname)
  // if filespacegroupname and spacegroupname_ indicate a transformation
  // between rhombohedral and hexagonal settings of rhombohedral lattice
  // reindex operator should already have been determined
  //
  // Returns true if so, after processing
  //         false if not so
  {
    char RH1 = filespacegroupname[0];
    char RH2 = spacegroupname_[0];
    if (!HorR(RH1)) {return false;}
    if (!HorR(RH2)) {return false;}
    if (RH1 == RH2) {return false;}

    // set things for report
    status = 0;
    SG_.init(spacegroupname_);
    newsymmetry_ = hkl_symmetry(spacegroupname_);
    samereferencegroup_ = true;
    spacegroupchanged_ = true;
    sameorigin_ = true;

    return true;
  }
  //--------------------------------------------------------------
}
