// referencelist.hh
//

#ifndef REFERENCELIST_HEADER
#define REFERENCELIST_HEADER

#include <cstring>

// Clipper
#include <clipper/clipper.h>
#include "clipper/clipper-ccp4.h"
using clipper::Message;
using clipper::Message_fatal;

#include "hkl_merge.hh"
#include "sdmodel.hh"
#include "referencescalemodel.hh"

namespace scala {
//======================================================================
  class ReferenceList {
    // container for a reference list of intensities, store in an
    // hkl_merge object
  public:
    ReferenceList() : hklrefname(""), inputtype(""), bulksolvent(false) {}

    // resoLimit < 0 to read whole file
    ReferenceList(const std::string& filename,
		  const std::string& labI, const std::string& labsigI,
		  const double& resoLimit,
		  const bool& verbose,
		  phaser_io::Output& output);

    void init(const std::string& filename,
	      const std::string& labI, const std::string& labsigI,
	      const double& resoLimit,
	      const bool& verbose,
	      phaser_io::Output& output);

    // Initialise from coordinate list
    void init(const std::string& xyzin,
	      const double& resoLimit,
	      const bool& verbose,
	      phaser_io::Output& output);

    // true if object is empty (no filename)
    bool IsEmpty() const {return (hklrefname == "");}

    // toleranceratio = 1.0 for difference > maximum resolution,
    //    larger tolerance is more lax
    // set status = -1 if the two lists have different symmetry (point group),
    // or +1 if cell is too different
    // return false if status != 0
    // see also private function
    bool checkCompatible(const hkl_unmerge_list& hkl_list,
			 const double& toleranceratio);


    bool scaleToObserved(const hkl_unmerge_list& hkl_list,
			 const int& datasetindex,
			 const SDmodel& SDM,
			 const double& toleranceratio,
			 phaser_io::Output& output);

    bool SFcalcScaleToObserved(const std::string& xyzin,
			       const hkl_unmerge_list& hkl_list,
			       const int& datasetindex,
			       const SDmodel& SDM,
			       const double& toleranceratio,
			       const bool& verbose,
			       phaser_io::Output& output);

    ReferenceScaleModel referenceScaleModel() const
    {return referencescalemodel;}

    hkl_symmetry symmetry() const {return hklmergelist.symmetry();}
    std::string SpaceGroupSymbol() const {return hklmergelist.SpaceGroupSymbol();}
    int num_obs() const {return hklmergelist.num_obs();}
    Scell Cell() const {return hklmergelist.Cell();}

    // column labels used
    std::vector<std::string> columnLabels() const;

    double resHigh() const {return resolimit;}

    // Return scaled I sigI for given hkl
    IsigI Isig(const Hkl& h) const;

    // Return unscaled I sigI for given hkl
    IsigI Isig0(const Hkl& h) const;

    //! format reason for any error
    std::string formatError() const
    {return referencescalemodel.formatError();}

    void recordScaleReference(phaser_io::Output& output) const;


  private:
    std::string hklrefname;
    std::string labI_;
    std::string labsigI_;
    double resolimit;

    std::string inputtype; // HKLREF or XYZIN
    bool bulksolvent;      // true if XYZIN and bulk solvent calculation

    hkl_merge hklmergelist;

    // Scales to test dataset
    ReferenceScaleModel referencescalemodel;

  };

} // namespace scala

#endif
