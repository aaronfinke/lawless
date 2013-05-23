// dataset.hh
//


#ifndef DATASET_HEADER
#define DATASET_HEADER

/*
// Clipper
#include <clipper/clipper.h>
#include "clipper/core/clipper_precision.h"
using clipper::ftype;
using clipper::Vec3;
using clipper::Mat33;
using clipper::String;
using clipper::Metric_tensor;
typedef clipper::Vec3<double> DVect3;
typedef clipper::Mat33<double> DMat33;
typedef clipper::Vec3<float> FVect3;
typedef clipper::Mat33<float> FMat33;
typedef clipper::Vec3<int> IVect3;

#include "ccp4/cmtzlib.h"    // CCP4 MTZlib headers (namespace CMtz)
#include "ccp4/csymlib.h"    // CCP4 symmetry stuff
#include "matvec_utils.hh"  // Matrix & vector utilities
#include "util.hh"
#include "range.hh"
*/

#include "hkl_datatypes.hh"

typedef float  Rtype;
typedef double Dtype;

typedef std::pair<int,int> IndexPair;
typedef std::pair<float,float> RPair;
typedef std::pair<double,double> DPair;


namespace scala
{
  //======================================================================
  class Dataset
  //! A set of Xdataset objects, since a dataset may come from multiple crystals
  {
  public:
    Dataset() {}
    //! construct from one Xdataset
    Dataset(const Xdataset& xdataset);

    void Clear() {xdatasets.clear();} //!< clear xdataset list

    //! return number of Xdatasets
    int Number() const {return xdatasets.size();}

    //! Add in a new xdataset if 1st or Dname is the same, return true if added
    bool AddXdataset(const Xdataset& xdataset);

    //! Add batch number to list for this dataset ID
    void add_batch(const int& setid, const int& batch_num);

    //! Add batch number to list for this PxdName
    void add_batch(const PxdName& pxdname, const int& batch_num);

    //! return list of PXDnames
    std::vector<PxdName> pxdnames() const; //!< return all PXD names

    //! return PXDname for given set ID
    PxdName pxdname(const int& setid) const;

    //! return consensus PXDname (just set Dname to "MultiCrystal")
    PxdName pxdname() const;

    //! return project name (all Xdatasets have same project name)
    std::string Pname() const;

    //! return dataset name (all Xdatasets have same dataset name)
    std::string Dname() const;

    //! return average (or sole) cell
    Scell cell() const;

    //! return average (or sole) wavelength
    float wavelength() const;

    //! return wavelength for named crystal
    float wavelength(const std::string& xname) const;

    //! return wavelength range
    Range wavelengthRange() const;

    //!< return average mosaicity
    float Mosaicity() const;

    void SetResRange(const ResoRange& resrange); //!< set resolution range
    ResoRange ResRange() const;   //!< return resolution range

    //! change basis: reindex to get new cell
    void change_basis(const ReindexOp& reindex_op);

    //! Set unit cells for all Xdatasets
    ///    void SetCell(const Scell& cell);
    void SetCellWavelength(const Scell& cell, const float& wavel);

    //! Set mosaicity for all Xdatasets
    void SetMosaicity(const float& mosaicity);

    //! add in another cell and wavelength, put into 1st Xdataset
    void AddCellWavelength(const Scell& newcell, const float& wavel);

    //! add to run index list for given Xdataset
    void AddRunIndex(const PxdName& pxdname, const int& RunIndex);

    //! clear all runs
    void ClearRunList();

    //!< return run index list
    std::vector<int> RunIndexList() const;

    //! return list of datasetIDs
    std::vector<int> SetIDs() const;

    //! return datasetID for given PXDname
    int GetID(const PxdName& pxdname) const;

    //! set datasetIDs to ID, ID+n etc, return last value used
    int StoreSetID(const int& ID);

    //! return true if xdataset with given PXD name is present here
    bool IsPxdPresent(const PxdName& pxdname) const;

    //! return true if xdataset with given ID is present here
    bool IsSetidPresent(const int& setid) const;

    void ClearBatchList();   //!< clear batch lists
    int num_batches() const; //!< number of batches in all Xdatasets

    //! return maximum ID of datasets here or given ID (if >0)
    int MaxID(const int& ID=-1) const;

    //! format all PxdNames
    std::string formatNames() const;

    std::string formatPrint() const; //!< format
    std::string format() const; //!< format all xdatasets

    //! equality, just tests all pxdnames
    friend bool operator == (const Dataset& a,const Dataset& b);

    UnitCellSet AllCellSet() const;
    std::vector<Scell> AllCells() const; //!< all cells

    std::vector<float> AllWavelengths() const;

    //! return number of cells/wavelengths
    int NumberofCells() const;

    std::string formatAllCells() const; //!< format cell & wavelength list if more than one

    //! return worst deviation (A), = 0 if only one
    double WorstDeviation() const;

  private:
    std::vector<Xdataset> xdatasets;

    // Get index for xdataset setid, = -1 if absent
    int XdatasetIndex(const int& setid) const;

    // Get index for xdataset PxdName, = -1 if absent
    int XdatasetIndex(const PxdName& pxdname) const;

    // Get index for Xname, = -1 if absent
    int XdatasetIndex(const std::string& xname) const;

    // Sanity check, all xdatasets should have same Dname
    void check() const;     // Dies here if not

    // Die here if nothing in list
    void dieIfEmpty(const std::string& where) const;

  };
  //======================================================================
  // Return true if dataset setid is in datasets list
  //  & return dataset index idataset (-1 if not)
  bool in_datasets(const int& setid,
		   const std::vector<Dataset>& datasets,
		   int& idataset);
  //======================================================================
  // Return true if dataset pxdname is in datasets list
  //  & return dataset index idataset (-1 if not)
  bool in_datasets(const PxdName& pxdname,
		   const std::vector<Dataset>& datasets,
		   int& idataset);
}  // namespace scala
#endif
// end  dataset.hh
