// mtz_merge_io.hh
//
// Phil Evans 2009
//
// MTZ io for merged files (read only)
// 


#ifndef MTZ_MERGE_IO_HEADER
#define MTZ_MERGE_IO_HEADER

#include <vector>

#include "hkl_controls.hh"

// Clipper
#include "clipper/clipper.h"
#include "clipper/clipper-ccp4.h"
using clipper::Message;
using clipper::Message_fatal;


#include "hkl_datatypes.hh"
#include "controls.hh"
#include "mtz_merge_io.hh"
#include "range.hh"
#include "columnlabels.hh"
#include "hkl_symmetry.hh"

namespace MtzIO {

  class MtzMrgFile {
  public:
    //! Constructor: does nothing
    MtzMrgFile();
    //! Destructor: close any file that was left open
    ~MtzMrgFile();

    // Open file for reading, read header & close again
    bool open_read(const std::string& filename);

    // Read all selected data from MTZ file into clipper objects
    // hkl_info_list, IsigData
    ClipperLabelPair ReadData(clipper::CCP4MTZfile& mtzin,
		  const double& ResoLimit,
		  const MtzIO::column_labels& column_list,
		  clipper::HKL_info& hkl_info_list,
		  clipper::HKL_data<clipper::data32::I_sigI>& IsigData,
		  clipper::MTZdataset& mtzdataset,
		  const bool& verbose,
		  std::string& output);

    // Fill an unmerged hkl_list from a merged file
    // On entry:
    //  mtzname            name of MTZ file (or logical name)
    //  file_sel           flags for general selection
    //                     - resolution limits
    //  column_list        list of column names wanted by the program
    //  InputPxdName       PXD name to override dataset information from
    //                     MTZ file
    //  output             output string for printing
    //  verbose            set verbosity level
    //                      = 0 silent, = +1 usual summary
    //                      >= +2 debug
    // 
    // On exit:
    //  hkl_list  has been filled and closed
    FileRead MakeHklList(const std::string& mtzname,
			 file_select& file_sel, 
			 col_controls& column_selection,
			 MtzIO::column_labels& column_list,
			 const scala::PxdName& InputPxdName,
			 const scala::Scell& cell,
			 std::string& output,
			 const int& verbose,
			 hkl_unmerge_list& hkl_list);

    // Access to information
    // Filename
    std::string Filename() const {return filenamein;}
    // Returns true if merged, false if unmerged
    bool Merged() const {return merged;}
    // Maximum resolution
    clipper::Resolution MtzResolution() const {return mtzfile_resolution;}
    // Cell
    Scell Cell() const {return mcell;}
    std::string Spacegroupsymbol() const {return spacegroupsymbol;}

  private:
    std::string filenamein; // input file name
    bool fileopen;  // open_read has been called
    bool merged;

    clipper::Resolution mtzfile_resolution;
    clipper::MTZdataset mtzdataset;

    double ResMax;
    Scell mcell;
    std::string spacegroupsymbol;
    scala::SpaceGroup spacegroup;

    bool at_start;

    // private functions
    // Next index, returns false if end of list
    bool next(clipper::HKL_info::HKL_reference_index& hkl_index);
  };


} // namespace MtzIO
#endif
