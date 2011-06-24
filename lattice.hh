// lattice.hh

#ifndef LATTICE_HEADER
#define LATTICE_HEADER

#include <string>
#include <vector>

#include "hkl_datatypes.hh"

enum CrystalSystem {NOSYSTEM, TRICLINIC, MONOCLINIC, ORTHORHOMBIC, TETRAGONAL,
		    TRIGONAL, HEXAGONAL, CUBIC};

namespace scala
{
  class hkl_symmetry;
  //--------------------------------------------------------------
  //! If sgname corresponds to a rhombohedral lattice, return in
  //! hexagonal setting unless specified as R xxx :R
  //! Print warning message with leading tag
  std::string CheckInputRhombohedralSGname(const std::string& tag,
					   const std::string& sgname,
					   std::string& outstring);
  //--------------------------------------------------------------
  //! Return true if LatType is allowed lattice type (P, C, I etc)
  bool AllowedLatticeType(const char& LatType);
  //--------------------------------------------------------------
  //! return type = RH ('R' or 'H') if type == 'R' or 'H'
  char RhombohedralLatType(const char& type, const char& RH);
  //--------------------------------------------------------------
  //! True if not hexagonal setting (angles 90,90,120)
  bool RhombohedralAxes(const std::vector<double>& unit_cell_dimensions);
  //--------------------------------------------------------------
  //! replace 1st character of name with L if both are 'R' or 'H'
  std::string SetRlatticetype(const std::string& name, const char& L);
  //--------------------------------------------------------------
  //================================================================
  class CrystalType
  //! Crystal system + lattice centring
  {
  public:
    CrystalType()
    : crystalsystem(NOSYSTEM), latticetype(' ') {}
    
    CrystalType(const CrystalSystem& CrysSys, const char& LatType)
      : crystalsystem(CrysSys), latticetype(LatType) {}
    
    CrystalType(const hkl_symmetry& Symmetry);  // construct from symmetry object
    
    CrystalSystem crystalSystem() const {return crystalsystem;}
    
    bool isSet() const {return crystalsystem != NOSYSTEM;}
    
  // Return single letter type for crystal system + lattice type (short=true)
  // or longer version
    std::string format(const bool shortform=true) const;
    
    friend bool operator == (const CrystalType& a,const CrystalType& b);
    friend bool operator != (const CrystalType& a,const CrystalType& b);
    friend bool EquivalentCentredMonoclinic(CrystalType& a,const CrystalType& b);
    friend bool EquivalentRhombohedral(CrystalType& a,const CrystalType& b);
    
  private:
    CrystalSystem crystalsystem;
    char latticetype;
  };

}


#endif
