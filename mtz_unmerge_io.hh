// mtz_unmerge_io.hh
// Phil Evans August 2003
//
// MTZ io for  Unmerged files
//
// Written for Pointless


#ifndef MTZ_IO_HEAD
#define MTZ_IO_HEAD

#include <string>
#include <map>

// Clipper
#include "clipper/clipper.h"
#include "clipper/clipper-ccp4.h"
using clipper::Message;
using clipper::Message_fatal;


// CCP4
#include "csymlib.h"    // CCP4 symmetry stuff
#include "ccp4_parser.h"
#include "ccp4_general.h"
#include "cmtzlib.h"    // CCP4 MTZlib headers (namespace CMtz)

#include "hkl_datatypes.hh"
#include "controls.hh"
#include "columnlabels.hh"
#include "hkl_unmerge.hh"
#include "range.hh"
#include "hash.hh"


using namespace scala;

namespace scala{class FileRead;}

namespace MtzIO {  
  //--------------------------------------------------------------
  class MtzUnmrgFile
  //! Class to handle reading unmerged MTZ files: this class does not write MTZ files
    {
    public:
      //! Constructor: does nothing
      MtzUnmrgFile();
      //! Destructor: close any file that was left open
      ~MtzUnmrgFile();
      
      //! Open a file for read access, returns true if OK
      bool open_read( const std::string filename_in );
      //! Close a file after reading
      void close_read();

      // NB this code does not deal with writing files!
      
      //! Fill hkl_list with default options, select everything
      FileRead FillHklList(const std::string& mtzname,
			   std::string& output,
			   const int& verbose,
			   hkl_unmerge_list& hkl_list);

      //! Add this file to an hkl_unmerge_list object

      //! 1. If the hkl_list object is empty, put all data from this MTZ file
      //!   into hkl_list, but leave it open
      //! 2. If the hkl_list object is NOT empty, check if the new file is compatible
      //!   with the current list
      //!   (a) If it is compatible, add this file into hkl_list, but leave it open
      //!   (b) If it is not compatible, close MTZ file & return
      //!
      //! Returns true if the MTZ file has been read, false if it is incompatible
      FileRead AddHklList(const int& fileSeries,
			  const std::string& mtzname,
			  file_select& file_sel, 
			  col_controls& column_selection,
			  MtzIO::column_labels& column_list,
			  const all_controls& controls,
			  const scala::PxdName& InputPxdName,
			  const scala::Scell& cell,
			  const double& cellTolerance,
			  std::string& output,
			  const int& verbose,
			  hkl_unmerge_list& hkl_list);
      //!< On entry:
      //!< \param fileSeries    index number for file or file-series (from 1)
      //!<                      for batch exclusion in file_sel
      //!< \param mtzname       name of MTZ file (or logical name)
      //!< \param file_sel      flags for general selection
      //!<                       - dataset selection
      //!<                       - batch selection
      //!<                       - resolution limits
      //!<                       - detector coordinate rejection ranges
      //!< \param column_selection   flags for column selection
      //!<                       - PROFILE or INTEGRATED
      //!< \param column_list   list of column names wanted by the program
      //!< \param controls      run controls, partial controls
      //!< \param InputPxdName  PXD name to override dataset information from
      //!<                          MTZ file (forces one dataset)
      //!< \param output        output string for printing
      //!< \param verbose       set verbosity level
      //!<                        = 0 silent, = +1 usual summary
      //!<                       >= +2 debug
      //!< 
      //!< On exit:
      //!< \param hkl_list  has been filled or appended to, but not closed
      //!<              Subsequent calls are needed to "CloseDatasetBatch" &
      //!<              "close_part_list"
      
      //! Return number of reflections in file
      int Nrec() const {return Nrecl_file;}
      //! Return number of columns in file
      int Ncols() const {return Ncolumns;}
      //! Return true if file is sorted
      bool Sorted() const {return sorted;}
      //! True if file is merged
      bool Merged() const {return merged;}
      //! return spacegroup
      scala::SpaceGroup Spacegroup() const {return spacegroup_;}
      //! return batch list for all batches in file
      std::vector<Batch> BatchList();
      //! return cell
      scala::Scell Cell() const {return averagecell;}

      //! Return mtz structure
      CMtz::MTZ* mtz_file() const {return mtzin;}
      //! Return resolution range in file
      ResoRange ResRangeFile() const
      {return ResoRange(Range(minres, maxres));}
      //! Accepted resolution range
      ResoRange ResRange() const {return resrange;}

      //! reset to beginning of file for RRefl
      void Rewind();
      //! Read one record from MTZ file (raw data return)
      bool Rrefl(std::vector<float>& cols, std::vector<bool>& col_mnf);

      /*! Get kdataset'th dataset in crystal/dataset list
	Returns false if non-existent */
      bool get_dataset(const int& kdataset, Xdataset& this_dataset) const;
      /*! Get kbatch'th batch in batch list
	Returns false if non-existent */
      bool get_batch(const int& kbatch,  CMtz::MTZBAT& this_batch) const;
      //! Return title
      std::string get_title() const {return title;}
      //! Return filename
      std::string FilenameIn() const {return filename_in_;}
      //! Return column numbers and flags
      column_select ColumnSelect() const {return col_select;}

      column_labels ColumnLabels()const {return column_list;}

    private:
      enum MTZmode { NONE, READ, WRITE, APPEND };
      //! mtz object
      CMtz::MTZ* mtzin;
      MTZmode mode;
      std::string filename_in_, filename_out_;
      scala::SpaceGroup spacegroup_;
      CMtz::SYMGRP mtzsym;  // symmetry from MTZ file
      std::string title;

      std::vector<Xdataset> xdatasets;  // list of crystal/datasets in file
      std::vector<CMtz::MTZBAT*> mtzbatches;       // list of batches in file

      std::vector<Xdataset> datasets;  // Selected list
      std::vector<Batch> batches;      // List of batches with acceptances
      std::vector<Run> runs;
      std::vector<int> offsets;        // batch offsets for each run

      Scell averagecell;
      bool sorted;
      bool merged;

      int IrefCnt;
      int Ncolumns;

      int Nrecl_file;
      float minres, maxres;  // Resolution range in file
      ResoRange resrange;    // Resolution range accepted
      int Nrej_batch;  // Number of observation parts rejected by batch

      column_select col_select; // column selection flags
      MtzIO::column_labels column_list; // column list

      // Private member functions
      int read_datasets(const CMtz::MTZ* mtzin,
			 std::vector<Xdataset>& Xdatasets);
      bool read_batches(const CMtz::MTZ* mtzin, 
			std::vector<CMtz::MTZBAT*>& mtzbatches);

      Scell get_dset_batch_info(const int& fileSeries,
				const file_select& file_sel,
				const scala::PxdName& InputPxdName,
				const scala::Scell& cell,
				const column_select& col_sel,
				bool& DifferentCell);

      int get_refs(hkl_unmerge_list& hkl_list,
		   file_select& file_sel, 
		   const column_select& col_sel,
		   const Scell& averagecell,
		   bool& ChangeIndex);

      bool IsCompatible(const hkl_unmerge_list& hkl_list,
			const double& cellTolerance) const;
	
      // Make a list of runs from batches (for rebatching)
      void MakeRuns();

      // Return true if file is sorted and the index was not changed
      bool FileSorted() const;

      void OffsetBatches();

      // Get file column numbers for labels in column list
      // On exit:
      //   ColumnLabels  contains actual column numbers for labels found
      // Fails if compulsory column not found
      void get_col_lookup(column_labels& ColumnLabels);


  }; //   class MtzUnmrgFile;



}

#endif
