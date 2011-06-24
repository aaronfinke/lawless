
#include <string>
#include "lattice.hh"
#include "hkl_symmetry.hh"
#include "util.hh"
#include "string_util.hh"

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;

namespace scala
{
  //--------------------------------------------------------------
  std::string CheckInputRhombohedralSGname(const std::string& tag,
					   const std::string& sgname,
					   std::string& outstring)
  //! If sgname corresponds to a rhombohedral lattice, return in
  //! hexagonal setting unless specified as R xxx :R
  //! Print warning message with leading tag
  {
    char initialLatType = sgname[0];
    std::string revisedname = SGnameHtoR(sgname, 'H');  // H unless :R
    char newLatType = revisedname[0];
    if (initialLatType == 'R' && newLatType != 'R') {
      outstring = FormatOutput::logTab(0,
		    std::string("\n$TEXT:Warning:$$ $$\nWARNING: ")+
		    tag+"\n   Rhombohedral group name "+sgname+
		    " will be used in the hexagonal (H) setting "+revisedname+
		    "\n   To get the rhombohedral lattice setting, give name as "+
		    SGnameHtoR(sgname,'R')+" :R\n$$\n");
    } else if (initialLatType == 'R' && newLatType == 'R') {
      outstring = FormatOutput::logTab(0,
		    std::string("\n$TEXT:Warning:$$ $$\nWARNING: ")+
		    tag+"\n   Rhombohedral group name "+sgname+
		    " will be used in the rhombohedral setting "+revisedname+"\n$$\n");
    } else if (newLatType == 'H') {
      outstring = FormatOutput::logTab(0,
		    std::string("\n$TEXT:Warning:$$ $$\nWARNING: ")+
		    tag+"\n   Rhombohedral group name "+sgname+
		    " will be used in the hexagonal [H] setting "+revisedname+"\n$$\n");
    }
    return revisedname;
  }
  //--------------------------------------------------------------
  bool AllowedLatticeType(const char& LatType)
  //!< Returns true if LatType is allowed
  {
    char AllowedTypes[] = {'P','A','B','C','I','F','R','H'};
    for (size_t j=0;j<sizeof(AllowedTypes);j++) {
      if (LatType == AllowedTypes[j]) return true;
    }
    return false;
  }
  //--------------------------------------------------------------
  //! return type = RH ('R' or 'H') if type == 'R' or 'H'
  char RhombohedralLatType(const char& type, const char& RH)
  {
    if (type == 'R' || type == 'H') {
      if (!(RH == 'R' || RH == 'H')) {
	// invalid rhombohedral lattice type
	Message::message(Message_fatal("RhombohedralLatType must be R or H not "+RH) );
      }
      return RH;
    }
    return type;
  }
  //--------------------------------------------------------------
  //! replace 1st character of name with L if both are 'R' or 'H'
  std::string SetRlatticetype(const std::string& name, const char& L)
  {
    std::string s = name;
    if (s[0] == 'R' || s[0] == 'H') {
      if (L == 'R' || L == 'H') {
	s[0] = L;
      }
    }
    return s;
  }
  //--------------------------------------------------------------
  CrystalType::CrystalType(const hkl_symmetry& Symmetry)
  {
    crystalsystem = Symmetry.CrysSys();
    latticetype = Symmetry.lattice_type();
  }
  //--------------------------------------------------------------
  std::string CrystalType::format(const bool shortform) const
  {
    std::string CrysCode;
    if (shortform) {
      // Short form 
      if (crystalsystem == NOSYSTEM) return "00";
      if (crystalsystem == TRICLINIC)
	{CrysCode = 'a';} 
      else if (crystalsystem == MONOCLINIC)
	{CrysCode = 'm';}
      else if (crystalsystem == ORTHORHOMBIC)
	{CrysCode = 'o';}
      else if (crystalsystem == TETRAGONAL)
	{CrysCode = 't';}
      else if (crystalsystem == TRIGONAL || crystalsystem == HEXAGONAL)
	{CrysCode = 'h';}
      else if (crystalsystem == CUBIC)
	{CrysCode = 'c';}
      return CrysCode+latticetype;
    } else {
      // Long form
      if (crystalsystem == NOSYSTEM) return "undefined";
      if (crystalsystem == TRICLINIC)
 	{CrysCode = " triclinic";} 
      else if (crystalsystem == MONOCLINIC)
	{CrysCode = " monoclinic";}
      else if (crystalsystem == ORTHORHOMBIC)
	{CrysCode = " orthorhombic";}
      else if (crystalsystem == TETRAGONAL)
	{CrysCode = " tetragonal";}
      else if (crystalsystem == TRIGONAL)
	{CrysCode = " trigonal";}
      else if (crystalsystem == HEXAGONAL)
	{CrysCode = " hexagonal";}
      else if (crystalsystem == CUBIC)
	{CrysCode = " cubic";}
      std::string latcent;
      // Lattice types: 'P','A','B','C','I','F','R','H'
      if (latticetype == 'P') {latcent = "primitive";}
      else if (latticetype == 'A') {latcent = "A-centred";}
      else if (latticetype == 'B') {latcent = "B-centred";}
      else if (latticetype == 'C') {latcent = "C-centred";}
      else if (latticetype == 'I') {latcent = "body-centred";}
      else if (latticetype == 'F') {latcent = "face-centred";}
      else if (latticetype == 'R') {latcent = "rhombohedral (R)";}
      else if (latticetype == 'H') {latcent = "rhombohedral (H)";}
      return latcent+CrysCode;;
    }
  }
  //--------------------------------------------------------------
  bool operator == (const CrystalType& a,const CrystalType& b)
    // Crystal types are equal if their short formatted codes are
    // equal: this equates trigonal & heaxagonal
  {
    return (a.format() == b.format());
  }
  //--------------------------------------------------------------
  bool operator != (const CrystalType& a,const CrystalType& b)
    // Crystal types are equal if their short formatted codes are
    // equal: this equates trigonal & heaxagonal
  {
    return !(a.format() == b.format());
  }
  //--------------------------------------------------------------
  bool EquivalentCentredMonoclinic(CrystalType& a,const CrystalType& b)
  // True if types are either equal or both centred monoclinic (mC or mI)
  {
    if (a == b) return true;
    if (a.crystalsystem == MONOCLINIC && b.crystalsystem == MONOCLINIC) {
      if ((a.latticetype == 'C' || a.latticetype == 'I') &&
	  (b.latticetype == 'C' || b.latticetype == 'I')) {
	return true;
      }
    }
    return false;
  }
  //--------------------------------------------------------------
  bool EquivalentRhombohedral(CrystalType& a,const CrystalType& b)
  // True if types are either equal or both rhombohedral (R or H)
  {
    if (a == b) return true;
    if (a.crystalsystem == TRIGONAL && b.crystalsystem == TRIGONAL) {
      if ((a.latticetype == 'R' || a.latticetype == 'H') &&
	  (b.latticetype == 'R' || b.latticetype == 'H')) {
	return true;
      }
    }
    return false;
  }
  //--------------------------------------------------------------
  bool RhombohedralAxes(const std::vector<double>& unit_cell_dimensions)
    // True if not hexagonal setting (angles 90,90,120)
  {
    double tol =  0.2;  // tolerance on angles
    if (Close<double,double>(unit_cell_dimensions.at(3), 90.0, tol))
      if (Close<double,double>(unit_cell_dimensions.at(4), 90.0, tol))
	if (Close<double,double>(unit_cell_dimensions.at(5), 120.0, tol))
	  return false;
    return true;
  }
  //--------------------------------------------------------------
}
