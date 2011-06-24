// mtz_unmerge_io.cpp
//
// mtz_unmerge_io   io package for unmerged mtz files
//
//   Some of this code is copied from Kevin Cowtan's
//   ccp4_mtz_io class
//
//   Uses cmtzlib functions
//
// All in MtzIO namespace
//

#include <algorithm>
#define ASSERT assert
#include <assert.h>

#include "mtz_unmerge_io.hh"
#include "mtz_utils.hh"
#include "scala_util.hh"
#include "checkcompatiblesymmetry.hh"
#include "columnlabels.hh"
#include "string_util.hh"
#include "openinputfile.hh"

namespace MtzIO 
{
  //--------------------------------------------------------------
  /*! Constructing an MtzUnmrgFile does nothing except flag the object as not
    attached to any file for either input or output */
  MtzUnmrgFile::MtzUnmrgFile()
  {
    mode = NONE;
    mtzsym.spcgrp = -1;  // set to null
  }

  //--------------------------------------------------------------
  /*! Close any files which were left open. */
  MtzUnmrgFile::~MtzUnmrgFile()
  {
    close_read();
  }
  //--------------------------------------------------------------
  /*! The file is opened for reading. This MtzUnmrgFile object will
    remain attached to this file until it is closed. Until that occurs,
    no other file may be opened with this object, however another
    MtzUnmrgFile object could be used to access another file. 
    \param  filename_in The input filename or pathname. */
  bool MtzUnmrgFile::open_read(const std::string filename_in)
  // returns false if fails
  {
    if ( mode != NONE )
      Message::message( Message_fatal( "MtzUnmrgFile: open_read - File already open" ) );
    if ( filename_in == "") 
      Message::message( Message_fatal( "MtzUnmrgFile: open_read - no filename given" ) );

    // store filename
    if (getenv(filename_in.c_str()) != NULL)
      filename_in_ = std::string(getenv(filename_in.c_str()));
    else
      filename_in_ = filename_in;

    // open file
    mtzin = CMtz::MtzGet( filename_in_.c_str(), 0 );
    if ( mtzin == NULL) 
      return false;
      //      Message::message( Message_fatal( "MtzUnmrgFile: open_read - failed assignment" ) );

    if (!CMtz::MtzAssignHKLtoBase( mtzin )) {return false;}
    // get the list of datasets and crystals
    // (returns 0 if no datasets in file & one was created)
    int Ndatasets_file = read_datasets(mtzin, xdatasets);
    if (Ndatasets_file <= 0) {
      {Message::message(
			Message_info( "MtzUnmrgFile: Warning, no datasets in file" ) );}
    }

    // get list of batches
    if (!read_batches(mtzin, mtzbatches)) {
      merged = true;
      return false;
    }
    merged = false;

    // Columns
    Ncolumns = CMtz::MtzNcol(mtzin);
    // Make list of required columns in unmerged HKLIN file
    column_list = MtzIO::setup_columns();
    // Set into column_list the actual column numbers for columns
    //    requested by program (in column_list)
    get_col_lookup(column_list);

    // Pick up global file parameters
    // get spacegroup by decoding symops
    bool no_symm = true;
    if (mtzsym.spcgrp > 0) {
      // We already have MTZ-style symmetry from an hkl_list object, check that it is
      // the same as this one
      if (CmtzSymgrpEqual(mtzsym, mtzin->mtzsymm)) {
	// Yes it is the same
	no_symm = false;
      }
    }
    if (no_symm) {
      mtzsym = mtzin->mtzsymm;  // Store MTZ symmetry
      std::vector<clipper::Symop> ops = ClipperSymopsFromMtzSYMGRP(mtzsym);
      spacegroup_.init(ops);
    }

    // get resolution:
    CMtz::MtzResLimits( mtzin, &minres, &maxres );

    Nrecl_file = CMtz::MtzNref(mtzin);

    // Title
    title = std::string(mtzin->title);

    // Is the file sorted?
    sorted = FileSorted();

    IrefCnt = 0;

    mode = READ;
    return true;
  }
  //--------------------------------------------------------------
  bool MtzUnmrgFile::get_dataset(const int& kdataset,
				 Xdataset& this_dataset) const
  // Get kdataset'th dataset in crystal/dataset list
  // kdataset from 0
  // Returns false if non-existent
  {
    // Find requested dataset
    if (unsigned(kdataset) < xdatasets.size())  {
      this_dataset = xdatasets[kdataset];
      return true;
    }
    // Not found
    return false;	  
  }
  //--------------------------------------------------------------
  bool MtzUnmrgFile::get_batch(const int& kbatch,
			       CMtz::MTZBAT& this_batch) const
  // Get kbatch'th batch in batch list
  // kbatch from 0
  // Returns false if non-existent
  {
    // Find requested batch
    if (unsigned(kbatch) < mtzbatches.size()) {
      this_batch = *mtzbatches[kbatch];
      return true;
    }
    // Not found
    return false;	  
  }
  //--------------------------------------------------------------
  bool MtzUnmrgFile::FileSorted() const
  // true if file is sorted
  // only true if labels are set and correct
  {
    std::vector<std::string> keys(5);

    keys[0] = "H";
    keys[1] = "K";
    keys[2] = "L";
    keys[3] = "M_ISYM";   // "/" changed to "_" in mtz library
    keys[4] = "BATCH";

    bool ok = true;

    for (int i=0;i<5;i++)
      {
	if (mtzin->order[i] != NULL)
	  if (mtzin->order[i]->label != NULL)
	    if (mtzin->order[i]->label == keys[i]) 
	      continue;
	ok = false;
	break;
      }
    return ok;
  }
  //--------------------------------------------------------------
  // return batch list for all batches in file
  std::vector<Batch> MtzUnmrgFile::BatchList()
  {
    if (batches.size() > 0) return batches; // list already filled (in AddHklList)
    // store all batches unconditionally
    batches.clear();
    CMtz::MTZBAT this_batch;
    bool accept = true;
    int idataset;
    int j = 0;
    while (get_batch(j, this_batch)) {
      idataset = this_batch.nbsetid;
      Batch batch(this_batch, accept, idataset);
      // store batch with "accept" flag
      batches.push_back(batch);
      ++j;
    }
    return batches;
  }
  //--------------------------------------------------------------
  void MtzUnmrgFile::Rewind()
  // reset to beginning of file for RRefl
  {
    CCP4::ccp4_file_seek(mtzin->filein, SIZE1, SEEK_SET); 
    IrefCnt = 0;
  }
  //--------------------------------------------------------------
  bool MtzUnmrgFile::Rrefl(std::vector<float>& cols, std::vector<bool>& col_mnf)
  // Read one reflection record into vector cols
  // vector col_mnf is true if data item missing (MNF), false if OK
  // return False on beyond last reflection
  {
    if (++IrefCnt <= Nrecl_file) {
      // NB the construction &*vector.begin() gives a pointer to
      //    the data part of a vector, which must have already been
      //    set up with the correct length
      CMtz::MtzRrefl(mtzin->filein, Ncolumns, &*cols.begin());
      for (int i = 0; i < Ncolumns; i++)
	col_mnf[i] = CMtz::ccp4_ismnf(mtzin, cols[i]);
      return true;
    } else {
      return false;
    }
  }
  //--------------------------------------------------------------
  FileRead MtzUnmrgFile::FillHklList(const std::string& mtzname,
				    std::string& output,
				    const int& verbose,
				    hkl_unmerge_list& hkl_list)
  // Fill hkl_list with default options
  {
    int fileSeries = 1;
    // File selection flags (resolution, datasets, batches etc)
    file_select file_sel;
    // Set Profile-fitted [default] or integrated intensity
    col_controls column_selection; 
    // Scala control classes (default settings)
    //  run controls
    //  partials controls
    all_controls controls;
    PxdName InputPxdName;
    Scell cell;
    double cellTolerance(2.0);
    return AddHklList(fileSeries, mtzname, file_sel, column_selection,
		      column_list, controls,  InputPxdName, cell, cellTolerance,
		      output, verbose, hkl_list);
  }
  //--------------------------------------------------------------
  FileRead MtzUnmrgFile::AddHklList(const int& fileSeries,
				    const std::string& mtzname,
				    file_select& file_sel, 
				    col_controls& column_selection,
				    MtzIO::column_labels& Column_list,
				    const all_controls& controls,
				    const scala::PxdName& InputPxdName,
				    const scala::Scell& cell,
				    const double& cellTolerance,
				    std::string& output,
				    const int& verbose,
				    hkl_unmerge_list& hkl_list)
  //
  // Add this file to an hkl_unmerge_list object
  //
  // 1. If the hkl_list object is empty, put all data from this MTZ file
  //   into hkl_list, but leave it open
  // 2. If the hkl_list object is NOT empty, check if the new file is compatible
  //   with the current list
  //   (a) If it is compatible, add this file into hkl_list, but leave it open
  //   (b) If it is not compatible, close MTZ file & return
  //
  // Returns FileRead.Opened() true if the MTZ file has been successfully opened
  // Returns FileRead.Read() true if the MTZ file has been read, false if it is incompatible
  //  
  // On entry:
  //  fileSeries         index number for file or file-series (from 1)
  //                     for batch exclusion in file_sel (MTZ files only)
  //  mtzname            name of MTZ file (or logical name)
  //  file_sel           flags for general selection
  //                     - dataset selection
  //                     - batch selection
  //                     - resolution limits
  //                     - detector coordinate rejection ranges
  //                     - input scale factor (MULTIPLY)
  //  column_selection   flags for column selection
  //                      - PROFILE or INTEGRATED
  //  column_list        list of column names wanted by the program
  //  controls           run controls, partial controls
  //  InputPxdName       PXD name to override dataset information from
  //                     MTZ file (forces one dataset)
  //  cell               if non-null, replace all cells from MTZ file
  //  output             output string for printing
  //  verbose            set verbosity level
  //                      = 0 silent, = +1 usual summary
  //                      >= +2 debug
  // 
  // On exit:
  //  hkl_list  has been filled, but not organised & partials
  //            assigned: this needs a call to "prepare" or
  //            "change_symmetry"
  //
  {
    column_list = Column_list;
    output = "";
    // If hkl_list has something in it already, pick up its MTZ-style symmetry,
    // if present
    if (!hkl_list.IsEmpty()) {
      mtzsym = hkl_list.MtzSym();
      std::vector<clipper::Symop> ops = ClipperSymopsFromMtzSYMGRP(mtzsym);
      spacegroup_.init(ops);
    } else {
      mtzsym.spcgrp = -1;  // set to null  //^!!
    }

    // Open MTZ file for reading if necessary
    bool opened = false;
    //^    std::cout <<"Mode: " << mode<<"\n"; //^
    if (mode == NONE)      {
      // Read headers etc
      opened = open_read(mtzname);
    } else if (mode == READ)      {
      // File already opened for reading, check it is the same file
      //^      std::cout << "Filenames:"<<filename_in_<<":"<<mtzname<<"\n"; //^
      if (filename_in_ != mtzname)	{
	// if not, close that file & open this one
	close_read();
	opened = open_read(mtzname);
      }
      opened = true;
    } else {
      Message::message
	(Message_fatal("MtzUnmrgFile::AddHklList: file not opened READ"));
    }
    if (!opened) {
      // Failed to open file
      return FileRead(false, false, false, 0);
    }
    if (merged) {  // file must be unmerged for this function
      Message::message
	(Message_fatal("MtzUnmrgFile::AddHklList: not an unmerged file"));
    }

    //   transfer the actual column numbers  into col_select for file reading
    col_select = column_select(column_list, column_selection);
    // Fudge for blank ROT column: if so flag to replace by batch number
    if (column_list.CNL("ROT").valuerange.AbsRange() < 0.001) {
      col_select.col_Rot = -1;
    }
    if (col_select.col_Rot < 0) {
      output += FormatOutput::logTab(0, 
   "**** WARNING: missing or empty ROT column in input file, BATCH number will be used instead");
    }

    // Extract selected dataset & batch information from mtz file,
    //  including unit cell things ready for resolution calculations
    bool DifferentCell;
    averagecell = get_dset_batch_info(fileSeries, file_sel, InputPxdName, cell,
				      col_select, DifferentCell);
    Scell accepted_cell = averagecell;

    if (DifferentCell) {
      output += FormatOutput::logTab(0, 
        "**** WARNING: input CELL is significantly different from cell from HKLIN file");
      output += FormatOutput::logTabPrintf(1,"Average HKLIN cell: ");
      for (int i=0;i<6;i++) output += FormatOutput::logTabPrintf(0,"%6.1f",averagecell[i]);
      output += FormatOutput::logTab(0,"\n");
      output += FormatOutput::logTabPrintf(1,"Input cell:         ");
      for (int i=0;i<6;i++) output += FormatOutput::logTabPrintf(0,"%6.1f",cell[i]);
      output += FormatOutput::logTab(0,"\n\n");
      accepted_cell = cell;
    }
    // MTZ headers read

    MakeRuns();   // Runs for these batches
    bool first = true;


    // Is hkl_list empty?
    if (hkl_list.IsEmpty()) {
      // Empty list, initialise
      // Initialise hkl_unmerge_list object
      hkl_list.init(title, Nrecl_file,
		    hkl_symmetry(spacegroup_), controls);
      offsets.assign(runs.size(),0);  // clear offsets
      hkl_list.SetMtzSym(mtzsym);
    } else {
      // Not empty, check for compatibility
      first = false;
      if (!IsCompatible(hkl_list, cellTolerance)) {
	// Not compatible, exit
	// Close mtz file
	close_read();
	return FileRead(true, false, false, Nrej_batch);
      }
      // otherwise, carry on
      // Apply batch offsets if any
      offsets = CompareRunRanges(hkl_list.RunList(), runs);
      ASSERT (offsets.size() == runs.size());
    }

    // Read all observations into hkl_list, subject to selection flags
    bool ChangeIndex;
    int Nread = get_refs(hkl_list, file_sel, col_select, accepted_cell, ChangeIndex);
    if (Nread <= 0)
      Message::message(Message_fatal
		       ("hkl_unmerge_list:: No reflections read") );

    // list is sorted if file was & no index is changed, for 1st file only
    sorted = sorted && !ChangeIndex && first;
    // Close observation part list: this doesn't stop more additions later
    hkl_list.close_part_list(resrange, sorted);

    OffsetBatches();  // offset batch numbers
    // Add in dataset list & batch list, but don't close list
    hkl_list.AddDatasetBatch(datasets, batches);
    // Store data presence flags (won't hurt to do this again)
    hkl_list.StoreDataFlags(col_select.DataFlags());
    // Store file name
    hkl_list.AppendFileName(mtzname);


    // Close mtz file
    close_read();

    if (verbose > 1)  {
      //        output += FormatOutput::logTabPrintf(0,
      //			    "\n---------------------------------------------------------------\n");
      output += FormatOutput::logTabPrintf(0,
			  "\nReflection list generated from file: %s\n",filename_in_.c_str());
      output += FormatOutput::logTabPrintf(0,
			  "\nTitle: %s\n\n", title.c_str());
      output += FormatOutput::logTabPrintf(0,
			  "   Space group from HKLIN file : %s\n",
			  spacegroup_.Symbol_hm().c_str());
      output += FormatOutput::logTabPrintf(0, "   Cell: ");
      for (int i=0;i<6;i++) output += FormatOutput::logTabPrintf(0,"%7.2f",
						   accepted_cell[i]);
      output += FormatOutput::logTab(0,"\n");
      output += FormatOutput::logTabPrintf(0,
			  "   Resolution range in file:  %8.2f    %8.2f\n",
			  ResRangeFile().ResLow(),
			  ResRangeFile().ResHigh());
      if (file_sel.Nrej_reso() > 0)
	output += FormatOutput::logTabPrintf(0,
			    "   Number of observation parts outside resolution limits = %d\n",
			    file_sel.Nrej_reso());
      if (file_sel.Nrej_mflag() > 0)
	output += FormatOutput::logTabPrintf(0,
			    "   Number rejected with M > 1 = %d\n",
			    file_sel.Nrej_mflag());
      bool offset = false;
      for (size_t i=0;i<offsets.size();++i) {
	if (offsets[i] != 0) offset = true;
      }
      if (offset) {
	if (runs.size() > 1) {
	  output += FormatOutput::logTabPrintf(0,
			      "   Batch numbers incremented by:");
	  for (size_t i=0;i<offsets.size();++i) {
	    output += FormatOutput::logTabPrintf(0," %5d", offsets[i]);
	  }
	  output += FormatOutput::logTabPrintf(0,"\n");
	}
      }
      //        output += FormatOutput::logTabPrintf(0,
      //			    "\n---------------------------------------------------------------\n");
    }
    return FileRead(true, true, true, Nrej_batch);
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  Scell MtzUnmrgFile::get_dset_batch_info(const int& fileSeries,
					  const file_select& file_sel,
					  const scala::PxdName& InputPxdName,
					  const scala::Scell& cell,
					  const column_select& col_sel,
					  bool& DifferentCell)
  // Select dataset & batch information from MTZ object into
  // hkl_list object for wanted datasets & batches
  //
  // Only wanted datasets are stored, but all batches in the
  // file are stored, so that automatic run assignment will
  // work properly: however, for rejected batches, only the batch header will be stored.
  // The actual observations will not be.
  // Rejected batches are flagged, including those from rejected datasets. 
  //
  // On entry:
  //  fileSeries     index number for file or file-series (from 1)
  //                 for batch exclusion in file_sel
  //  file_sel       selection flags for datasets & batches
  //  InputPxdName   PXD name from input
  //                 if not blank, force all datasets into one and rename
  //  cell           if not null, replace cells with this one
  //  col_sel        column selection info
  //                     - column numbers for each required item
  // 
  // On exit:
  //  DifferentCell        true if cell is "different" from averagecell
  //
  //  datasets, batches    lists of datasets & batches
  //
  // Returns average cell from file
  {
    const double TOLERANCE = 3.0;   // angular difference limit for warning
    DifferentCell = false;
    bool NewCell = false;
    if (cell[0] != 0.0) {NewCell = true;}   // use input cell

    datasets.clear();
    // If InputPxdName not blank, then we are going to put all input datasets
    // into one, but after rejecting datasets (if needed)
    bool ForceOneDataset = false;
    if (!InputPxdName.is_blank()) ForceOneDataset = true;

    // Select datasets from this MTZ file
    Xdataset this_dataset;
    int k = 0;
    while (get_dataset(k, this_dataset)) {
      // Do we want this dataset?
      if (file_sel.accept_dataset(this_dataset.pxdname())) {
	// yes, wanted
	datasets.push_back(this_dataset);
      }       
      ++k;
    }

    // Average unit cells over all datasets & store average
    int ndatasets = datasets.size();
    Scell averagecell = AverageDsetCell(datasets);
    float averagewvl = AverageDsetWavelength(datasets);

    Scell accepted_cell = averagecell;

    // If cell given, check agreement
    if (NewCell) {
      if (!averagecell.equalsTol(cell, TOLERANCE)) {
	DifferentCell = true;
      }
      accepted_cell = cell;
      for (size_t id=0;id<datasets.size();id++) {
      	datasets[id].cell() = accepted_cell;   // reset all dataset cells to 
      }
    }

    Xdataset OneDataset;
    int setid = 0;
    if (ForceOneDataset) {
      // Make new datasets to include everything
      OneDataset = Xdataset(InputPxdName, accepted_cell, averagewvl, setid);
    }

    // Now batches: store all batches even if not needed
    batches.clear();
    CMtz::MTZBAT this_batch;
    bool accept = false;
    int idataset;
    int j = 0;
    while (get_batch(j, this_batch)) {
      if (in_datasets(this_batch.nbsetid, datasets, idataset)) {
	// This batch is in accepted dataset
	// Do we want this batch? (batch exclusions etc)
	accept = file_sel.accept_batch(this_batch.num, fileSeries);
	// Store list of batch numbers for this dataset
	if (ForceOneDataset) {
	  OneDataset.add_batch(this_batch.num);
	  idataset = setid;
	} else {
	  datasets[idataset].add_batch(this_batch.num);
	}
      } else {
	accept = false;
	if (ForceOneDataset) {
	  idataset = setid;
	}
      }
      // Check if there is any valid time information
      if (col_sel.col_time < 0) {
	// No time column
	this_batch.time1 = 0.0;
	this_batch.time2 = 0.0;
      }
      Batch batch(this_batch, accept, idataset);
      if (NewCell) {
	batch.SetCell(accepted_cell);
      }
      // store batch with "accept" flag
      // idataset = -1 for rejected datasets
      batches.push_back(batch);
      ++j;
    }
    if (ForceOneDataset) {
      // Reset to one dataset
      datasets.clear();
      datasets.push_back(OneDataset);
      ndatasets = datasets.size();
    }	
    std::sort (batches.begin(), batches.end());
	    
    return averagecell;
  }
  //--------------------------------------------------------------
  float check_column(const std::vector<float>& cols,
		     const std::vector<bool>& col_mnf,
                     const int& mcol, bool& status)
  // extract column from reflection record, checking for MNFs
  //
  // On entry:
  //  cols         reflection record
  //  col_mnf      true for any column which is MNF, false for OK
  //  mcol         column index required (from 0), =-1 if absent from file
  //
  // On exit:
  //  returns value extracted (or default = 0.0)
  //  status       true if MNF found, else false
  //
  {
    float col_value = 0.0;  // set default = 0.0
    status = false;
    
    if (mcol >= 0)
      {
        // check column for MNF
        if (col_mnf[mcol]) 
          {
            status = true;
            return col_value;
          }
        col_value = cols[mcol];
      }
    return col_value;
  }
  //---------------------------------------------------------------------------
  int MtzUnmrgFile::get_refs(hkl_unmerge_list& hkl_list,
			     file_select& file_sel, 
			     const column_select& col_sel,
			     const Scell& averagecell,
			     bool& ChangeIndex)
  //
  // Read all (selected) reflections from MTZ file into hkl_unmerge object
  //
  // On entry:
  //  file_sel           flags for general selection
  //                     - dataset selection
  //                     - batch selection
  //                     - resolution limits
  //                     - detector coordinate rejection ranges
  //                     - input scale factor (MULTIPLY)
  //  col_sel     column selection info
  //                     - column numbers for each required item
  //
  // On exit:
  //  hkl_list        filled list  
  //  ChangeIndex     true if index changed
  //
  // Returns number of observation parts read
  {
    Hkl hkl;
    int misym, batch;
    Rtype I, sigI, Ipr, sigIpr,
      fraction_calc, Xdet, Ydet, phi,
      time, width, LP,  BgPkRatio;
    int Mpart, IObsFlag, Npart, Ipart;

    // Column buffers of correct length (as in file)
    std::vector<float>  cols(Ncols());
    std::vector<bool>   col_mnf(Ncols());

    Rtype scale = 1.0;
    Rtype sigscale = 0.0;
    bool StatusFlag;
    Nrej_batch = 0;

    int nbatches = batches.size();
    hash_table batch_lookup;
    // Make lookup table (hash table)
    // Setup up hash lookup table a bit larger than required
    batch_lookup.set_size( int(1.2 * nbatches));
    for (int i = 0; i < nbatches; i++)  {
      batch_lookup.add(batches[i].num(), i);
    }

    // For resolution range found in file
    Range InvResRange;

    // MTZ symmetry
    ChangeIndex = false;
    bool changeSymmetry = false;
    hkl_symmetry FileSym(spacegroup_);

    if (FileSym != hkl_list.symmetry()) {
      ChangeIndex = true;
      changeSymmetry = true;
    }

    int nread = 0;

    // <<<< Loop all reflection records
    while (Rrefl(cols, col_mnf)) {
      // check and copy all compulsory columns, MNF not allowed (fatal)
      bool flag = false;
      hkl.h() = Nint(check_column(cols, col_mnf, col_sel.col_h, StatusFlag));
      flag = flag || StatusFlag;
      hkl.k() = Nint(check_column(cols, col_mnf, col_sel.col_k, StatusFlag));
      flag = flag || StatusFlag;
      hkl.l() = Nint(check_column(cols, col_mnf, col_sel.col_l, StatusFlag));
      flag = flag || StatusFlag;
      misym = Nint(check_column(cols, col_mnf, col_sel.col_misym, StatusFlag));
      flag = flag || StatusFlag;
      int Mflag = misym/256;
      int isym = misym - Mflag*256;
      batch = Nint(check_column(cols, col_mnf, col_sel.col_batch, StatusFlag));
      flag = flag || StatusFlag;
      
      I = check_column(cols, col_mnf, col_sel.col_I, StatusFlag);
      flag = flag || StatusFlag;
      sigI = check_column(cols, col_mnf, col_sel.col_sigI, StatusFlag);
      flag = flag || StatusFlag;
      
      if (flag)
	Message::message(
			 Message_fatal( "get_refs: MNF in compulsory column near hkl "+hkl.format() ) );
      
      // >>> Rejection tests
      // Rejected batch (or dataset)
      if (!batches[batch_lookup.lookup(batch)].Accepted()) {
	Nrej_batch++; // count excluded records
	continue;
      }
      
      // Resolution range
      double s =  hkl.invresolsq(averagecell);
      if (! file_sel.in_reslimits(s)) {
	file_sel.incr_rej_reso();
	continue;
      }
      InvResRange.update(s);  //smin, smax
      
      // Mflag
      if (Mflag > 1 ) {
	file_sel.incr_rej_mflag();
	continue;
      }
      // <<<
      
      // check_column(const std::vector<float>& cols, std::vector<bool> col_mnf,
      //                int& mcol, float& col_default, float&  col_value)
      
      // Optional columns, set defaults if absent or MNF
      Ipr = check_column(cols, col_mnf, col_sel.col_Ipr, StatusFlag);
      sigIpr = check_column(cols, col_mnf, col_sel.col_sigIpr, StatusFlag);
      fraction_calc = check_column(cols, col_mnf, col_sel.col_fractioncalc, StatusFlag);
      Xdet = check_column(cols, col_mnf, col_sel.col_Xdet, StatusFlag);
      Ydet = check_column(cols, col_mnf, col_sel.col_Ydet, StatusFlag);
      phi = check_column(cols, col_mnf, col_sel.col_Rot, StatusFlag);
      width = check_column(cols, col_mnf, col_sel.col_Width, StatusFlag);
      LP = check_column(cols, col_mnf, col_sel.col_LP, StatusFlag);
      IObsFlag = Nint(check_column(cols, col_mnf, col_sel.col_ObsFlag, StatusFlag));
      BgPkRatio = check_column(cols, col_mnf, col_sel.col_BgPkRatio, StatusFlag);
      
      // phi default = batch number
      if (col_sel.col_Rot < 0) {
	phi = batch;
      }
      // time defaults = phi (Rot)
      if (col_sel.col_time < 0) {
	time = phi;
      } else {
	time = check_column(cols, col_mnf, col_sel.col_time, StatusFlag);
      }
      // Possible input scale
      sigscale = check_column(cols, col_mnf, col_sel.col_sigscale, StatusFlag);
      if (col_sel.col_scale >= 0) 
	{
	  scale = check_column(cols, col_mnf, col_sel.col_scale, StatusFlag);
	  // If the scale column is present but there is no valid
	  // scale then skip this observation
	  if (StatusFlag || scale == 0.0) 
	    continue;
	  // Apply input scale immediately
	  I *= scale;
	  sigI = sqrt(scale*sigI*scale*sigI + sigscale*I*sigscale*I);
	  if (col_sel.col_Ipr > 0) {
	    Ipr *= scale;
	    sigIpr = sqrt(scale*sigIpr*scale*sigIpr + sigscale*Ipr*sigscale*Ipr);
	  }
	}

      if (changeSymmetry) {
	//  reduce hkl to asymmetric unit
	int new_isym;
	Hkl hkl_new = hkl_list.symmetry().put_in_asu(FileSym.get_from_asu(hkl,isym), new_isym);
	if (new_isym != isym) {
	  // changed from input
	  ChangeIndex = true;
	}
	isym = new_isym;
	hkl = hkl_new;
      }

      // Process partial flags Mflag and Mpart
      Npart = 1;  // Default full, one part
      Ipart = 1;
      if (Mflag == 1) {
	// Partial
	Npart = -1;  // Number of parts unknown
	if (col_sel.col_Mpart > 0) {
	  Mpart = Nint(cols[col_sel.col_Mpart]);
	  if (Mpart == 10) 
	    // previously summed partial, treat as full
	    Npart = 1;
	  // Unpack predicted number of parts and serial
	  else if (Mpart > 200) {
	    Npart = Mpart/100;
	    Ipart = Mpart%100;
	  } else if (Mpart > 20) {
	    Npart = Mpart/10;
	    Ipart = Mpart%10;
	  }
	}
      }
        
      // Apply input scale (MULTIPLY)
      I *= file_sel.InputScale();
      sigI *= file_sel.InputScale();
      Ipr *= file_sel.InputScale();
      sigIpr *= file_sel.InputScale();
      
      // Offset batch number
      int irun = batches[batch_lookup.lookup(batch)].RunIndex();
      batch += offsets[irun];
      
      // Store this observation
      hkl_list.store_part(hkl, isym, batch, I, sigI, Ipr, sigIpr,
			  Xdet, Ydet, phi, time,
			  fraction_calc, width, LP,
			  Npart, Ipart, ObservationFlag(IObsFlag, BgPkRatio));
      nread++;
    }

    // Store accepted resolution range for this file
    resrange = ResoRange(InvResRange);
    return nread;
  }
  //--------------------------------------------------------------
  void MtzUnmrgFile::MakeRuns()
  // Make a list of runs from batches, just based on batch number
  {
    runs.clear();
    Run run;
    int lb = -1;
    int runindex = 0;
    for (size_t ib=0;ib<batches.size();++ib) {
      if ((lb < 0) || (batches[ib].num() == lb+1)) {
	run.AddBatch(batches[ib].num());
      } else {
	run.SortList();
	runs.push_back(run);
	runindex++;
	run.clear();
	run.AddBatch(batches[ib].num());
	lb = batches[ib].num();
      }
      batches[ib].SetRunIndex(runindex);
    }
    if (run.Nbatches() > 0) {
      runs.push_back(run);
    }
  }
  //--------------------------------------------------------------
  bool MtzUnmrgFile::IsCompatible(const hkl_unmerge_list& hkl_list,
				  const double& cellTolerance) const
  // Is the new MTZ file (header read) compatible with the previous list?
  {
    // Symmetry
    if (!CheckCompatibleSymmetry(hkl_list.symmetry(),
				hkl_symmetry(spacegroup_),
				false)) {
      // Symmetry fail
      return false;
    }
    //  Symmetry is OK, check cell
    if (!hkl_list.Cell().equalsTol(averagecell, cellTolerance)) {
      return false;
    }
    // Compatible
    return true;
  }
  //--------------------------------------------------------------
  void  MtzUnmrgFile::OffsetBatches()
  // Apply offsets to batches
  {
    for (size_t ib=0;ib<batches.size();++ib) {
      int irun = batches[ib].RunIndex();
      batches[ib].OffsetNum(offsets[irun]);
    }
  }
  //--------------------------------------------------------------
  //^ Close a file after reading
  void MtzUnmrgFile::close_read() {
    if (mode == READ) {
      MtzFree(mtzin);}
    mtzin = 0;
    mode = NONE;
  }
  //--------------------------------------------------------------
  void MtzUnmrgFile::get_col_lookup(column_labels& ColumnLabels)
  // Get file column numbers for labels in column list
  // Updates ColumnLabels
  //
  // On exit:
  //   ColumnLabels  contains actual column numbers for labels found
  //
  // Fails if compulsory column not found
  {
    if (ColumnLabels.size() == 0) {
      Message::message(Message_fatal("MtzUnrgFile::get_col_lookup - no columns in list")); 
    }
    ColumnLabels.start();  // start loop on column data
    ColumnNumberLabel CNL;

    while (ColumnLabels.next(CNL)) {
      CMtz::MTZCOL * col_data =
	CMtz::MtzColLookup(mtzin, CNL.label.c_str());  // lookup label in MTZ structure
      if (col_data) {
	// Column found, store index (from 1)
	CNL.number = col_data->source;
	CNL.type = col_data->type;
	CNL.valuerange = Range(col_data->min, col_data->max);
	ColumnLabels.Store(CNL); // store back it current position
      } else {
	if (CNL.number < 0) {  // compulsory column not found
	  Message::message( Message_fatal(
		  "Compulsory column not in input file - " + CNL.label));
	}
      }
    }    
  }
  //--------------------------------------------------------------
  /* Read crystals and datasets from mtzin */
  int MtzUnmrgFile::read_datasets(const CMtz::MTZ* mtzin,
				  std::vector<Xdataset>& xdatasets)
  {
    xdatasets.clear();
    Scell OverallCell;

    // Loop crystals
    for (int x=0; x < CMtz::MtzNxtal(mtzin); x++) {
      CMtz::MTZXTAL* xtl = CMtz::MtzIxtal(mtzin,x);

      // Loop datasets within crystal
      for (int s=0; s < CMtz::MtzNsetsInXtal(xtl); s++) {
	CMtz::MTZSET* set = CMtz::MtzIsetInXtal(xtl,s);
	// Don't store HKL_base, except the overall cell in case it is needed
	if (std::string(set->dname) == "HKL_base")
	  {
	    OverallCell = Scell(xtl->cell);
	  }
	else
	  {
	    Xdataset xdts(PxdName(xtl->pname, xtl->xname, set->dname),
				Scell(xtl->cell), set->wavelength, set->setid);
	    xdatasets.push_back(xdts);
	  }
      }
    }
    // If no datasets, create one
    int Ndatasets = xdatasets.size();
    if (Ndatasets == 0) {
      float wavelength = 0.0;
      int setid = 1;
      //      std::cout << "\n>>> WARNING: no datasets in file, creating one <<<\n";
      Xdataset xdts(PxdName("UnspecifiedProject",
			    "UnspecifiedCrystal",
			    "UnspecifiedDataset"),
		    OverallCell, wavelength, setid);
      xdatasets.push_back(xdts);
    }
    return Ndatasets;
  }
  //--------------------------------------------------------------
  bool MtzUnmrgFile::read_batches(const CMtz::MTZ* mtzin, 
				  std::vector<CMtz::MTZBAT*>& mtzbatches)
  // read list of batches from mtzin
  // Returns false if error
  {
    if (CMtz::MtzNbat(mtzin) == 0) {
      //*      {Message::message(
      //*			Message_info( "MtzUnmrgFile: no batches in file" ) );}
      return false;
    }
    mtzbatches.clear();
    int nbat = 0; 
    int idataset;
    // Batches are stored as linked list
    CMtz::MTZBAT *batch = mtzin->batch;

    while (batch != NULL) {
      mtzbatches.push_back(batch);
      // find dataset and add this batch to its list
      //  if nbsetid = 0 assign to first dataset
      if (in_datasets(batch->nbsetid, xdatasets, idataset)) {
	xdatasets[idataset].add_batch(batch->num);
      }
      batch = batch->next;  // pointer to next batch
      nbat++;
    }

    if (nbat != CMtz::MtzNbat(mtzin)) {
      Message::message(
		       Message_info( "MtzUnmrgFile: read_batches - wrong number of batch headers in file:" ) );
      return false;
    }
    return true;
  }
  //--------------------------------------------------------------
} // namespace MtzIO 


