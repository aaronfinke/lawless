// mergedlist.hh
//
// Merged data for one or more datasets, in clipper classes


#ifndef MERGEDLIST_HEADER
#define MERGEDLIST_HEADER

#include <vector>
#include <string>

// Clipper
#include <clipper/clipper.h>
#include "ccp4_mtz_io.h"
#include "clipper/clipper-ccp4.h"

#include "hkl_unmerge.hh"
#include "sdmodel.hh"
#include "alt_hkl_datatypes.h"

namespace scala {

  class MergedDatasetIntensities
  //! merged intensity data for one dataset
  {
  public:
    MergedDatasetIntensities(){}
    void init(const clipper::HKL_info& hkl_info_list,
	      const clipper::Cell& ccell);

    clipper::HKL_data<clipper::data32::I_sigI> Imean;   // Imean
    clipper::HKL_data<clipper::data32::J_sigJ_ano> Ipm; // I+, I-
    clipper::MTZdataset cset; 
    clipper::MTZcrystal cxtl; 
    clipper::String mtzpath;
    clipper::String mtzpathImean; // xname/dname/[IMEAN,SIGIMEAN]
    clipper::String mtzpathIpm;  // xname/dname/[I(+), SIGI(+), I(-), SIGI(-)]
  };
  //================================================================
  class MergedList {
    //! Merged data for one or more datasets
 
  public:
    MergedList(){}
    //! Fill from unmerged list & SD_model
    MergedList(const hkl_unmerge_list& hkl_list, const SDmodel& SDM);
    void init(const hkl_unmerge_list& hkl_list, const SDmodel& SDM);

    // Write data for datasetIndex to MTZ file
    // Return number of reflections written
    int WriteDatasetToMTZ(const std::string& outfilename,
			  const int& datasetIndex) const;
    int WriteDatasetToSCA(const std::string& outfilename,
			  const int& datasetIndex) const;

    int NumberDatasets() const {return ndatasets;}

    std::vector<Xdataset> Datasets() const {return xdatasets;};

    // return max(1/d^2) for each dataset
    std::vector<float> InvResMax() const {return resmaxdts;}
    // max(1/d^2) for given dataset, = 0 if unset
    float InvResMax(const int& datasetIndex) const;

    // Data access
    //! return reference to reflection list 
    clipper::HKL_info& HKLinfo() {return hkl_info_list;}
    //! return reference to Imean data for given dataset
    clipper::HKL_data<clipper::data32::I_sigI>&
    ImeanForDataset(const int& datasetindex);
   

  private:
    int ndatasets; // number of datasets
    std::vector<Xdataset> xdatasets;
    clipper::HKL_info hkl_info_list;  // hkl list
    // Merged intensities for each dataset
    std::vector<MergedDatasetIntensities> datasetdata;
    std::vector<int> nrefdts;  // number of reflections in each dataset
    std::vector<float> resmaxdts;  // maximum resolution (1/d^2) in each dataset    
    std::string title;
    double maxintensity;

    char spg_status; // aka spg_confidence in MTZ

    // returns false if I or sigI are Nan or sig = 0
    bool CheckNullImean(const clipper::data32::I_sigI& MIsig) const;
    // returns 0 if OK, -1 if both missing, +1 if I+ missing +2 if I- missing
    int CheckNullIano(const clipper::data32::J_sigJ_ano& MIsig) const;


  }; // MergedList

} // namespace scala 

#endif

