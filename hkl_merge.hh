// hkl_merge.hh
// a cut-down version of hkl_merged_list, just merged MTZ input, but including
// SF calculation
//

#ifndef HKL_MERGE_HEADER
#define HKL_MERGE_HEADER

#include <cstring>
#include <string>
#include <vector>

// Clipper
#include <clipper/clipper.h>
#include "clipper/clipper-ccp4.h"
using clipper::Message;
using clipper::Message_fatal;
using clipper::data32::F_phi;

#include "hkl_datatypes.hh"
#include "hkl_symmetry.hh"
#include "range.hh"
#include "mtz_merge_io.hh"
#include "columnlabels.hh"
#include "alt_hkl_datatypes.h"
#include "mergedlist.hh"

// phaser_io
#include "Output.hh"
using phaser_io::LOGFILE;

namespace MLIST
{
  enum mlist_status{EMPTY, HEADER, DATA, UNMERGED};
}

namespace scala{

  std::vector<clipper::String>
  ProcessLabels(const std::vector<clipper::String>& ColLab,
		const std::string& ColumnLabel,
		bool& IorF);
//======================================================================
  class hkl_merge
  //
  // Clipper object containing a set of merged data (intensities)
  // constructed from a file of hkl I (or F)
  //
  // This class must not be copied
  {
  public:
    hkl_merge() : status(MLIST::EMPTY), mtzread(0) {}
    //Construct from mtz file: just header stuff
    hkl_merge(const std::string& hklinname,
	      const bool& verbose,
	      phaser_io::Output& output);

    // Copy constructor throws exception
    hkl_merge(const hkl_merge& List);
    // Copy operator throws exception
    hkl_merge& operator= (const hkl_merge& List);

    // Initialise from MTZ file
    void init(const std::string& hklinname,
	      const bool& verbose,
	      phaser_io::Output& output);

    void read(const MtzIO::column_labels& column_list,
	      const double& ResoLimit,
	      const bool& verbose,
	      phaser_io::Output& output);

    // Initialise list from structure factor calculation from XYZIN coordinate file
    //  xyzin          file name for coordinate file
    //  spacegroup     spacegroup given on input to override that in xyzin file
    //                 blank if not given
    //  input_cell    cell if given on input
    //  FCresolution   maximum resolution to generate
    void CreateFromAtoms(const std::string& xyzin, 
			 const std::string& spacegroup,
			 const Scell& input_cell,
			 const double& FCresolution,
			 const bool& verbose,
			 phaser_io::Output& output);
    
    // Initialise list from structure factor calculation from XYZIN coordinate file
    //  with bulk solvent scaled to Fobs list
    //  xyzin          file name for coordinate file
    //  
    void createFromAtomsBulk(const std::string& xyzin,
			     const MergedList& mergedobslist,
			     const bool& verbose,
			     phaser_io::Output& output);


    // Access:
    ResoRange ResRange() const {return ResoRange(10000., ResMax);}
    double resHigh() const {return ResMax;}
    hkl_symmetry symmetry() const;
    std::string SpaceGroupSymbol() const;

    int num_obs() const;

    Scell Cell() const {return fcell;}
    // column labels used
    std::vector<std::string> columnLabels() const
    {return mtzfilein.columnLabelsUsed();}
    // true if generated from MTZ file, else from coordinates
    bool fromReflectionData() const;
    // true if from amplitudes F^2, else from intensities
    bool Amplitudes() const;

    // Return list of alternative indexing schemes
    //  max_delta tolerance on cell
    //    std::vector<scala::ReindexOp> OtherIndexing(const float& max_delta) const;

    // Return I sigI for given hkl
    IsigI Isig(const Hkl& h) const;

    // Reset current reflection pointer to first reflection
    void start() const;
    // Next IsigI, returns false if end of list
    bool next(IsigI& Is) const;
    // get hkl for current reflection
    Hkl  hkl() const;
    // resolution of current reflection
    double invresolsq() const;

    void PrintHeaderStuff(const bool& isreference,
			  phaser_io::Output& output) const;

    clipper::Spacegroup spacegroup() const;
    clipper::Cell cell() const;

  private:
    MLIST::mlist_status status;
    clipper::HKL_info hkl_info_list;
    clipper::HKL_data<clipper::data32::I_sigI> IsigData; // Imean
    clipper::HKL_data<clipper::data32::J_sigJ_ano> Ipm; // I+, I-
    mutable clipper::HKL_info::HKL_reference_index hkl_index;
    mutable bool at_start;

    int mtzread;  // +1 if mtzfile was read, -1 coordinates, 0 nothing
    MtzIO::MtzMrgFile mtzfilein;
    Scell fcell;
    hkl_symmetry fsymmetry;

    std::string filename;
    double ResMax;

    // Private function: Calculate structure factors
    // On input:
    //  xyzin         file name for coordinate file
    //  spacegroup     spacegroup given on input to override that in xyzin file
    //                 blank if not given
    //  input_cell    cell if given on input
    //  resolution    maximum resolution
    //
    // On exit:
    //  hkls          hkl list
    //  fc            Fcalc list
    void SFcalc(const std::string& xyzin,
		const std::string& spacegroup,
		const Scell& input_cell,
		const clipper::ftype64& resolution,
		clipper::HKL_info& hkls,
		clipper::HKL_data<F_phi>& fc);

    void SFcalcBulk(const std::string& xyzin,
		    const clipper::HKL_data<clipper::data32::F_sigF> fobs,
		    clipper::HKL_data<F_phi>& fc,
		    phaser_io::Output& output);

  };
}
#endif
