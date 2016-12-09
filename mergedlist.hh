// mergedlist.hh
//
// Merged data for one or more datasets, in clipper classes


#ifndef MERGEDLIST_HEADER
#define MERGEDLIST_HEADER

#include <vector>
#include <string>

// Clipper
#include <clipper/clipper.h>
// FIXME clipper update    #include "ccp4_mtz_io.h"
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
    clipper::HKL_data<clipper::data32::Flag> Npl; // N+
    clipper::HKL_data<clipper::data32::Flag> Nmn; // N-
    clipper::MTZdataset cset; 
    clipper::MTZcrystal cxtl; 
    clipper::String mtzpath;
    clipper::String mtzpathImean; // xname/dname/[IMEAN,SIGIMEAN]
    clipper::String mtzpathIpm;  // xname/dname/[I(+), SIGI(+), I(-), SIGI(-)]
    clipper::String mtzpathNpl;  // xname/dname/[N(+)]
    clipper::String mtzpathNmn;  // xname/dname/[N(-)]
  };
  //================================================================
  class MergedList {
    //! Merged data for one or more datasets
 
  public:
    MergedList(){}
    //! Fill from unmerged list & SD_model, for given dataset if index >=0
    //  default all datasets stored separately
    // Datasetindex = -2 to combine datasets together
    MergedList(const hkl_unmerge_list& hkl_list, const SDmodel& SDM,
	       const std::string& Title, const int& Datasetindex=-1);
    void init(const hkl_unmerge_list& hkl_list, const SDmodel& SDM,
	      const std::string& Title, const int& Datasetindex=-1);

    // Write data for datasetIndex to MTZ file
    // Return number of reflections written
    int WriteDatasetToMTZ(const std::string& outfilename,
			  const int& datasetIndex) const;
    int WriteDatasetToSCA(const std::string& outfilename,
			  const int& datasetIndex,
			  phaser_io::Output& output) const;

    int NumberDatasets() const {return ndatasets;}

    std::vector<Dataset> Datasets() const {return datasets;};

    // return max(1/d^2) for each dataset
    std::vector<float> InvResMax() const {return resmaxdts;}
    // max(1/d^2) for given dataset, = 0 if unset
    float InvResMax(const int& datasetIndex) const;
    double resHigh() const;  // for all datasets

    // Data access
    // Reset current reflection pointer to first reflection, for given dataset
    void start(const int& datasetIndex) const;
    // Next IsigI, returns false if end of list
    bool next(IsigI& Is) const;
    // get hkl for current reflection
    Hkl  hkl() const;
    // get clipper hkl for current reflection
    clipper::HKL HKL() const;

    //! return reference to reflection list 
    clipper::HKL_info& HKLinfo() {return hkl_info_list;}
    //! return reference to Imean data for given dataset
    const clipper::HKL_data<clipper::data32::I_sigI>&
    ImeanForDataset(const int& datasetindex) const;
    //! return reference to Imean data for given dataset
    clipper::HKL_data<clipper::data32::I_sigI>&
    ImeanForDataset(const int& datasetindex);

    // resolution of current reflection
    double invresolsq() const;

    clipper::Spacegroup spacegroup() const;
    clipper::Cell Cell() const;

    int num_reflections() const
    {return hkl_info_list.num_reflections();}

    double meanIntensity() const {return meanintensity;}
    double meanIovermeansigI() const {
      if (meansigintensity > 0.0) {
	return meanintensity/meansigintensity;
      }
      return 0.0;
    }

  private:
    int ndatasets; // number of datasets
    int dataset_index;  // dataset index, = -1 if all datasets stored
    std::vector<Dataset> datasets;
    clipper::HKL_info hkl_info_list;  // hkl list
    // Merged intensities for each dataset
    std::vector<MergedDatasetIntensities> datasetdata;
    std::vector<int> nrefdts;  // number of reflections in each dataset
    std::vector<float> resmaxdts;  // maximum resolution (1/d^2) in each dataset    
    std::string title;
    double maxintensity;
    double meanintensity;
    double meansigintensity;
    std::vector<std::string> historylines;

    char spg_status; // aka spg_confidence in MTZ

    // for access
    mutable int current_dataset_index;  // internal dataset index
    mutable clipper::HKL_info::HKL_reference_index hkl_index;
    mutable bool at_start;

    // returns false if I or sigI are Nan or sig = 0
    bool CheckNullImean(const clipper::data32::I_sigI& MIsig) const;
    // returns 0 if OK, -1 if both missing, +1 if I+ missing +2 if I- missing
    int CheckNullIano(const clipper::data32::J_sigJ_ano& MIsig) const;

    // return internal index to dataset datasetIndex
    // If dataset_index >=0, then only this dataset has been stored, so return 0
    // If dataset_index == -2, then all datasets have been stored together, return 0
    // If dataset_index <0, then all datasets have been stored, so return datasetIndex
    int InternalDTSindex(const int& datasetIndex) const;
  }; // MergedList

} // namespace scala 

#endif

