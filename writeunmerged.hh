// writeunmerged.hh

#ifndef WRITEUNMERGED_HEADER
#define WRITEUNMERGED_HEADER

#include <string>

#include "hkl_unmerge.hh"
#include "Output.hh"
#include "hkl_symmetry.hh"
#include "sdmodel.hh"

#include "ccp4/cmtzlib.h"    // CCP4 MTZlib headers (namespace CMtz)

namespace MtzIO
{
  //--------------------------------------------------------------
  int WriteUnmergedMTZ(const scala::hkl_unmerge_list& hkl_list,
			const scala::SDmodel& SDM,
			const bool& summedpartials,
			const int& datasetIndex,
			const std::string& filename_out,
			const std::string& title);
  //--------------------------------------------------------------
  int WriteUnmergedSCA(const scala::hkl_unmerge_list& hkl_list,
		       const scala::SDmodel& SDM,
		       const int& datasetIndex,
		       const std::string& filename_out,
		       const float& maxintensity);
  // Write unmerged scalepack file from hkl_unmerge_list object
  // returns number written
  //--------------------------------------------------------------
  // write out all parts, all datasets
  std::vector<int>  WriteParts(const scala::hkl_unmerge_list& hkl_list,
			       const int& NumCol,
			       CMtz::MTZ* mtzout,
			       CMtz::MTZCOL* col[]);
  //--------------------------------------------------------------
  std::vector<int>  WriteObservations(const scala::hkl_unmerge_list& hkl_list,
				      const scala::SDmodel& SDM,
				      const int& NumCol,
				      const int& datasetIndex,
				      CMtz::MTZ* mtzout,
				      CMtz::MTZCOL* col[]);
  // write out summed observations, for selected dataset(s)
  // omitting rejections
}
#endif
