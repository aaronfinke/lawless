// writeunmerged.hh

#ifndef WRITEUNMERGED_HEADER
#define WRITEUNMERGED_HEADER

#include <string>

#include "hkl_unmerge.hh"
#include "Output.hh"
#include "hkl_symmetry.hh"
#include "sdmodel.hh"

// CCP4
#include "ccp4/csymlib.h"    // CCP4 symmetry stuff
#include "ccp4/ccp4_general.h"
#include "ccp4/cmtzlib.h"    // CCP4 MTZlib headers (namespace CMtz)
using namespace CMtz;

namespace MtzIO
{
  class WriteUnmerged {
  public:
    WriteUnmerged();

    //--------------------------------------------------------------
    // Write unmerged MTZ file from hkl_unmerge_list object
    // returns number written
    int writeUnmergedMTZ(const scala::hkl_unmerge_list& hkl_list,
			 const scala::SDmodel& SDM,
			 const bool& summedpartials,
			 const int& datasetIndex,
			 const std::string& filename_out,
			 const std::string& title);
    //--------------------------------------------------------------
    // Write unmerged scalepack file from hkl_unmerge_list object
    // returns number written
    int writeUnmergedSCA(const scala::hkl_unmerge_list& hkl_list,
			 const scala::SDmodel& SDM,
			 const int& datasetIndex,
			 const std::string& filename_out,
			 const float& maxintensity);

    //--------------------------------------------------------------
    int Nreflections() const {return nref;} // number written
    int Nmultiple() const {return nmultiple;} // number of multiples

  private:
    int nref; // total number of refecltions written
    int nmultiple; // number of multiples

    //--------------------------------------------------------------
    // write out all parts, all datasets
    std::vector<int>  writeParts(const scala::hkl_unmerge_list& hkl_list,
				 const int& NumCol,
				 CMtz::MTZ* mtzout,
				 CMtz::MTZCOL* col[]);
    //--------------------------------------------------------------
    std::vector<int>  writeObservations(const scala::hkl_unmerge_list& hkl_list,
					const scala::SDmodel& SDM,
					const int& NumCol,
					const int& datasetIndex,
					CMtz::MTZ* mtzout,
					CMtz::MTZCOL* col[]);
    // write out summed observations, for selected dataset(s)
    // omitting rejections

    void OptAddCol(const bool& coln,
		   MTZCOL* col[], int& ic, MTZ* mtzout, MTZSET* baseset,
		   const char* label, const char* type);
    // Conditional column addition, only if coln > 0

  };
}
#endif
