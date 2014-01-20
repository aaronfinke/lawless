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
// If there are multiple lattices, then two possible schemes are used to represent
// the data in an MTZ file
// 
//  Scheme 1)    for each lattice
//   there will be 3 addition columns of the form Xn where "X" = H,K, or L
//   and n is the lattice number
//   Also a LATTNUM column, this is the lattice number for the main hkl
//   
//  Scheme 2) 
//   there will be 4 or 5 addition columns of the form
//    a) LATTNUMn where n is a sequential number (not the lattice number)
//    b) 3 columns Xn where "X" = H,K, or L
//    c) optionally a SCALEn column, giving the scale for this overlapped
//       reflection relative to the main HKL (LATTNUM)
//   Also a LATTNUM column, this is the lattice number for the main hkl
//   The number of extra column sets will depend on the maximum number of overlapped
//   spots in the dataset (nlatticecolumns)
//
// Scheme 2 may be distinguished by the presence of a LATTNUM1 column
// Note that in scheme 1 there is an entry Hn,Kn,Ln for the observation's main lattice
//   (ie that indicated by the LATTNUM column), but in scheme 2 the main lattice HKL is
//   not given in the extra index list, just in columns 1,2,3 (reduced hkl)

#include <algorithm>
#include <assert.h>
#define ASSERT assert

#include "mtz_unmerge_io.hh"
#include "mtz_utils.hh"
#include "scala_util.hh"
#include "checkcompatiblesymmetry.hh"
#include "columnlabels.hh"
#include "string_util.hh"
#include "openinputfile.hh"
//#include "timer.hh"

namespace MtzIO 
{
  //--------------------------------------------------------------
  /*! Constructing an MtzUnmrgFile does nothing except flag the object as not
    attached to any file for either input or output */
  MtzUnmrgFile::MtzUnmrgFile()
  {
    clear();
  }

  //--------------------------------------------------------------
  /*! Close any files which were left open. */
  MtzUnmrgFile::~MtzUnmrgFile()
  {
    close_read();
  }
  //--------------------------------------------------------------
  void MtzUnmrgFile::clear() {
    mode = NONE;
    mtzsym.spcgrp = -1;  // set to null
    mtzin = NULL;
    nlatticecolumns = 0;
    nlattices = 0;
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
    if ( mode != NONE ) {
      Message::message( Message_fatal( "MtzUnmrgFile: open_read - File already open" ) );}
    if ( filename_in == "") {
      Message::message( Message_fatal( "MtzUnmrgFile: open_read - no filename given" ) );}

    // store filename
    if (getenv(filename_in.c_str()) != NULL) {
      filename_in_ = std::string(getenv(filename_in.c_str()));
    } else{
      filename_in_ = filename_in;
    }
    // open file
    mtzin = CMtz::MtzGet( filename_in_.c_str(), 0 );
    if ( mtzin == NULL) {
      return false;
      //      Message::message( Message_fatal( "MtzUnmrgFile: open_read - failed assignment" ) );
    }
    if (!CMtz::MtzAssignHKLtoBase( mtzin )) {return false;}
    // get the list of datasets (fdatasets) from the file
    // (returns 0 if no datasets in file & one was created)
    int Ndatasets_file = read_datasets();
    if (Ndatasets_file <= 0) {
      {Message::message(
			Message_info( "MtzUnmrgFile: Warning, no datasets in file" ) );}
    }

    // get list of MTZ batches (vector<MTZBAT*> mtzbatches)
    // updates fdatasets with batch list
    if (!read_batches(mtzin, mtzbatches)) {
      merged = true;
      return false;
    }
    merged = false;

    // Do we have any multi-lattice entries in this file?
    //  If so count the columns
    CheckMultipleLattices();

    // Columns
    Ncolumns = CMtz::MtzNcol(mtzin);
    // Make list of required columns in unmerged HKLIN file (basic columns)
    column_label_list = MtzIO::setup_columns();
    column_label_list.SetNlatticeColumns(0);
    //
    if (nlatticecolumns > 0) {
      add_extra_columns(); // add extra lattice columns to column_label_list if required
    }
    // Set into column_label_list the actual column numbers for columns
    //    requested by program (in column_label_list)
    get_col_lookup(column_label_list);

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

    // History
    int nhistlines = mtzin->histlines;
    historylines.clear();
    char * phist = mtzin->hist;
    for (int i=0;i<nhistlines;++i) {
      std::string hline(phist, MTZRECORDLENGTH);
      historylines.push_back(hline);
      phist += MTZRECORDLENGTH;
    }

    // Is the file sorted?
    sorted = FileSorted();

    IrefCnt = 0;

    mode = READ;
    return true;
  }
  //--------------------------------------------------------------
  bool MtzUnmrgFile::get_dataset(const int& kdataset,
				 Dataset& this_dataset) const
  // Get kdataset'th dataset in dataset list
  // kdataset from 0
  // Returns false if non-existent
  {
    if (unsigned(kdataset) < fdatasets.size())  {
      this_dataset = fdatasets[kdataset];
      return true;
    }
    // Not found
    return false;	  
  }
  //--------------------------------------------------------------
  bool MtzUnmrgFile::get_batch(const int& kbatch,
			       CMtz::MTZBAT& this_batch) const
  // Get kbatch'th MTZ batch in batch list into this_batch
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
		      column_label_list, controls,  InputPxdName, cell, cellTolerance,
		      output, verbose, hkl_list);
  }
  //--------------------------------------------------------------
  FileRead MtzUnmrgFile::AddHklList(const int& fileSeries,
				    const std::string& mtzname,
				    file_select& file_sel, 
				    col_controls& column_selection,
				    const MtzIO::column_labels& column_label,
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
  //  column_label       list of column names wanted by the program, ignored if invalid
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
    if (column_label.Setup()) {
      column_label_list = column_label;  // only use if set
    }
    output = "";
    // If hkl_list has something in it already, pick up its MTZ-style symmetry,
    // if present
    if (!hkl_list.IsEmpty()) {
      mtzsym = hkl_list.MtzSym();
      std::vector<clipper::Symop> ops = ClipperSymopsFromMtzSYMGRP(mtzsym);
      spacegroup_.init(ops);
    } else {
      // Empty list, maybe containing symmetry
      mtzsym.spcgrp = -1;  // set to null
    }

    // Open MTZ file for reading if necessary
    bool opened = false;
    if (mode == NONE)      {
      // Read headers etc (MTZ datasets and batches, create array of Datasets, fdatasets)
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
    col_select = column_select(column_label_list, column_selection);
    // Fudge for blank ROT column: if so flag to replace by batch number
    if (column_label_list.CNL("ROT").valuerange.AbsRange() < 0.001) {
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
      // Is there symmetry in the "empty" list?
      hkl_symmetry symmset = hkl_list.symmetry();
      if (symmset.IsNull()) {
	symmset = hkl_symmetry(spacegroup_);
      } else {
	// Check for compatible symmetry
	if (! (symmset.CrysSys() == hkl_symmetry(spacegroup_).CrysSys())) {
	  std::string errormsg = FormatOutput::logTab(0, 
	     "**** ERROR: cannot combine files belonging to different crystal systems");
	  errormsg += "\n   Systems: "+symmset.formatCrysSys()+" : "+
	    hkl_symmetry(spacegroup_).formatCrysSys();
	  Message::message(Message_fatal
			   (errormsg+"\n**** Incompatible symmetries ****"));
	}
	if (spacegroup_.Symbol_hm() != symmset.symbol_xHM()) {
	  // Changing symmetry for this file
	  output += "\nChanging spacegroup on input from "+
	    spacegroup_.Symbol_hm()+" to "+symmset.symbol_xHM()+
	    +" to match first file\n";
	}
      }
      hkl_list.init(title, Nrecl_file,
		    symmset, controls);
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

      if (hkl_symmetry(spacegroup_) != hkl_list.symmetry()) {
	// Changing symmetry for this file
	output += "\nFor file "+mtzname+
	  "\n   change spacegroup on input from "+
	  spacegroup_.Symbol_hm()+" to "+hkl_list.symmetry().symbol_xHM()+
	  +" to match first file\n";
      }
    }

    // Append to history
    hkl_list.addHistory(historylines);

    // Read all observations into hkl_list, subject to selection flags
    bool ChangeIndex;
    //^    Timer timer;
    int Nread = get_refs(hkl_list, file_sel, col_select, accepted_cell, ChangeIndex);
    //^
    //    std::cout << "XDS::ReadObservations time " << timer.Dtime() << " elapsed " << timer.Etime() << "\n";
    //^-
    if (Nread <= 0)
      Message::message(Message_fatal
		       ("hkl_unmerge_list:: No reflections read") );

    // list is sorted if file was & no index is changed, for 1st file only
    sorted = sorted && !ChangeIndex && first;
    // Close observation part list: this doesn't stop more additions later
    hkl_list.close_part_list(resrange, sorted);

    OffsetBatches();  // offset batch numbers
    // Add in dataset list & batch list, but don't close list
    hkl_list.AddDatasetBatch(datasets, batches);  //  "datasets" from get_dset_batch_info
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
    Dataset this_dataset;
    int k = 0;
    while (get_dataset(k, this_dataset)) {
      // Do we want this dataset?
      // FIXME (see hkl_controls.hh) this always returns true for now, ie accept
      //   if we want to implement eg "EXCLUDE <pxdname>" we must remove it from the
      //   Dataset object and do something about it
      if (file_sel.accept_dataset(this_dataset.pxdnames())) {
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
	float wvl = datasets[id].wavelength();
	if (wvl < 0.001) {
	  wvl = averagewvl;
	}
      	datasets[id].SetCellWavelength(accepted_cell, wvl);   // reset all dataset cells to average
      }
    }

    PxdName pxdname;
    Dataset OneDataset;
    int setid = 1;
    if (ForceOneDataset) {
      // Make new dataset to include everything
      OneDataset.AddXdataset(Xdataset(InputPxdName, accepted_cell, averagewvl, setid));
      pxdname = InputPxdName; // for storing in batches
    }

    // Now batches: store all batches even if not needed
    batches.clear();
    CMtz::MTZBAT this_batch;
    bool accept = false;
    int idataset;
    int j = 0;
    while (get_batch(j, this_batch)) {
      if (!ForceOneDataset) {
	setid = this_batch.nbsetid; // SetID from file, if not OneDataset
      }	
      // in_datasets returns idataset as index into datasets array for found file setid
      if (in_datasets(this_batch.nbsetid, datasets, idataset)) {
	// This batch is in accepted dataset
	// Do we want this batch? (batch exclusions etc)
	accept = file_sel.accept_batch(this_batch.num, fileSeries);
	// Store list of batch numbers for this dataset
	if (ForceOneDataset) {
	  OneDataset.add_batch(setid, this_batch.num);
	  idataset = 0;
	} else {
	  datasets[idataset].add_batch(setid, this_batch.num);
	  pxdname = datasets[idataset].pxdname(setid); // name for this SetID
	}
      } else { // dataset not accepted, so reject batch
	accept = false;
	if (ForceOneDataset) {
	  idataset = 0; // one dataset
	} else {
	  pxdname = datasets[idataset].pxdname(setid); // name for this SetID
	}
      }
      // Check if there is any valid time information
      if (col_sel.col_time < 0) {
	// No time column
	this_batch.time1 = 0.0;
	this_batch.time2 = 0.0;
      }
      Batch batch(this_batch, accept, idataset);  // create Batch object
      batch.PXDname() = pxdname; 
      batch.DatasetID() = setid;
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
    
    if (mcol >= 0) {
      // check column for MNF
      if (col_mnf[mcol]) {
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

    // MAXNLATTICES is maximum number of lattices allowed
    numberinlattice.assign(MAXNLATTICES+1,0); // +1 as lattices are numbered from 1

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
      
      if (flag) {
	Message::message(
			 Message_fatal( "get_refs: MNF in compulsory column near hkl "+hkl.format() ) );
      }
      if (sigI <= 0.0) { // reject negative or zero sigma
	continue;
      }

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

      // Multiple lattice options
      std::vector<LatticeIndexInfo> lathkl;

      int latnum = 0;
      if (col_select.col_latnum > 0) {
	latnum = Nint(check_column(cols, col_mnf, col_select.col_latnum, StatusFlag));
	flag = flag || StatusFlag;
	// read extra hkl into lathkl, and lattnum, scale if scheme2
	Hkl hkln;
	for (int ih=0;ih<nlatticecolumns;++ih) {
	  bool flag = false;
	  int colnum = col_select.col_lathkl[ih];  // 1st column of group
	  // if scheme 2, read LATTNUMn
	  int latn = ih+1; // lattice number for this group if scheme 1
	  if (multilatscheme == 2) { // LATTNUMn
	    latn = Nint(check_column(cols, col_mnf, colnum, StatusFlag));
	    colnum++;
	  }
	  hkln.h() = Nint(check_column(cols, col_mnf, colnum, StatusFlag));
	  flag = flag || StatusFlag;
	  hkln.k() = Nint(check_column(cols, col_mnf, colnum+1, StatusFlag));
	  flag = flag || StatusFlag;
	  hkln.l() = Nint(check_column(cols, col_mnf, colnum+2, StatusFlag));
	  flag = flag || StatusFlag;
	  Rtype scale = 1.0;
	  if (col_select.col_latscale) {
	    scale = check_column(cols, col_mnf, colnum+3, StatusFlag);
	  }
	  if (hkln != Hkl(0,0,0)) {
	    // Store non-null extra indices with lattice number (from 1),
	    // but not if it belongs to the main lattice with the same hkl
	    if ((latn != latnum) ||
		(hkln != FileSym.get_from_asu(hkl,isym))) {
	      lathkl.push_back(LatticeIndexInfo(latn, hkln, scale));
	    }
	  }
	} // end loop lattices
 	if (flag) {
	  Message::message(
		 Message_fatal
		 ("get_refs: MNF in compulsory multilattice column near hkl "+hkl.format()));
	}
	// count entries for each lattice
	numberinlattice.at(latnum)++;
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
			  Npart, Ipart, ObservationFlag(IObsFlag, BgPkRatio),
			  latnum, lathkl);
      nread++;
    } // end loop read reflections

    // count lattices with non-zero entries
    nlattices = 0;
    if (nlatticecolumns > 0) {
      for (size_t j=1; j<numberinlattice.size(); j++) { // loop from 1
	if (numberinlattice[j] > 0) {
	  nlattices++;
	}
      }
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
      }
      lb = batches[ib].num();
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
  // Close a file after reading
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
  int MtzUnmrgFile::read_datasets()
  // sets fdatasets
  {
    fdatasets.clear();
    Scell OverallCell;
    columnlabels.clear();
    columntypes.clear();
    //^    std::cout << "\nMtzUnmrgFile::read_datasets\n"; //^
    // Loop crystals
    for (int x=0; x < CMtz::MtzNxtal(mtzin); x++) {
      CMtz::MTZXTAL* xtl = CMtz::MtzIxtal(mtzin,x);
      // Loop datasets within crystal
      for (int s=0; s < CMtz::MtzNsetsInXtal(xtl); s++) {
	CMtz::MTZSET* set = CMtz::MtzIsetInXtal(xtl,s);
	// Don't store HKL_base, except the overall cell in case it is needed
	if (std::string(set->dname) == "HKL_base") {
	  OverallCell = Scell(xtl->cell);
	} else {
	  Xdataset xdts(PxdName(xtl->pname, xtl->xname, set->dname),
			Scell(xtl->cell), set->wavelength, set->setid);
	  bool added = false;
	  if (fdatasets.size() > 0) {
	    // try to add this Xdataset to existing datasets
	    for (size_t idts=0;idts<fdatasets.size();++idts) {
	      added = fdatasets[idts].AddXdataset(xdts);
	      if (added) break;
	    }
	  }
	  if (!added) {
	    fdatasets.push_back(Dataset(xdts));
	  }
	}
	for (int c=0; c < CMtz::MtzNcolsInSet(set); c++) {
	  CMtz::MTZCOL* mc = CMtz::MtzIcolInSet(set,c);
	  std::string label(mc->label);
	  columnlabels.push_back(label);
	  std::string ctype(mc->type);
	  columntypes.push_back(ctype);
	  //^
	  //^	  std::cout << "Label, type " << c << " " << label
	  //^		    << " " << ctype << "\n";
	}
      }
    }
    // If no datasets, create a dummy one
    int Ndatasets = fdatasets.size();
    if (Ndatasets == 0) {
      float wavelength = 0.0;
      int setid = 1;
      //      std::cout << "\n>>> WARNING: no datasets in file, creating one <<<\n";
      Xdataset xdts(PxdName("UnspecifiedProject",
			    "UnspecifiedCrystal",
			    "UnspecifiedDataset"),
		    OverallCell, wavelength, setid);
      fdatasets.push_back(Dataset(xdts));
    }
    return Ndatasets;
  }
  //--------------------------------------------------------------
  bool MtzUnmrgFile::read_batches(const CMtz::MTZ* mtzin, 
				  std::vector<CMtz::MTZBAT*>& mtzbatches)
  // read list of batches from mtzin
  // returns batch list in mtzbatches
  // updates fdatasets with batch list
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
      //  if nbsetid = 0 assign to first dataset and reset setid to 1
      int setidb = batch->nbsetid;
      if (setidb <= 0) {
	setidb = 1;
	batch->nbsetid = setidb;
      }
      if (in_datasets(setidb, fdatasets, idataset)) {
	fdatasets[idataset].add_batch(setidb, batch->num);
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
  std::pair<int, int> 
  MtzUnmrgFile::CheckColumnlabelHKL(const std::string& label) const
  // Interpret a column label of the form H1, L3 etc
  // ie 1st character of column label is H, K or L & second character
  // a digit, then return
  //  first   1,2,3 for H, K, L
  //  second  digit value (lattice column number)
  // else 0,0
  {
    int jhkl = 0;
    int jlat = 0;
    char c1 = label[0];
    char c2 = label[1];
    if (c1 == 'H') {jhkl = 1;}
    else if (c1 == 'K') {jhkl = 2;}
    else if (c1 == 'L') {jhkl = 3;}
    if (jhkl > 0 && isdigit(c2)) {
      jlat = std::atoi(std::string(1,c2).c_str());
    } else {
      jhkl = 0;
    }
    return   std::pair<int, int>(jhkl, jlat);
  }
  //--------------------------------------------------------------
  int MtzUnmrgFile::CheckColumnlabelN(const std::string& label,
				      const std::string& basestring) const
  // Interpret a column label of the form <basestring>n, eg <basestring>1 etc
  // ie 1st part of column label == base & last character
  // a digit, then return last character as a digit, else 0 if fails
  {
    int jlat = 0;
    std::string::size_type i = label.find(basestring);
    if (i == std::string::npos) { // not found
      //      std::cout << "CheckColumnlabelN not found, " << label << " : " << basestring<<"\n"; //^-
      return jlat;
    }
    char c2 = label[label.size()-1];
    if (isdigit(c2)) {
      jlat = std::atoi(std::string(1,c2).c_str());
    }
    //^
    //    std::cout << "CheckColumnlabelN, " << label << " : " << basestring
    //	      <<" " << c2 <<" " <<jlat<<"\n"; //^-
    return jlat;
  }
  //--------------------------------------------------------------
  // Do we have any multi-lattice entries in this file?
  void MtzUnmrgFile::CheckMultipleLattices()
  // Multilattice if LATTNUM column present
  // sets nlatticecolumns extracted from columnlabels
  //
  // If there are multiple lattices, then
  // Scheme 1)    for each lattice
  // there will be 3 addition columns of the form Xn where "X" = H,K, or L
  // and n is the lattice number
  // Also a LATTNUM column, this is the lattice number for the main hkl
  // Scheme 2) 
  // there will be 4 or 5 addition columns of the form
  //   a) LATTNUMn where n is a sequential number
  //   b) 3 columns Xn where "X" = H,K, or L
  //   c) optionally a SCALEn column
  //  Also a LATTNUM column, this is the lattice number for the main hkl
  {
    nlatticecolumns = 0;  // usual case, no multiple lattices
    extracolumnlabels.clear();
    bool latnumcolumn = false;

    size_t ic=3;
    multilatscheme = 1;   // scheme 1

    while (ic<columnlabels.size()) { // skip 1st 3 columns
      if (columnlabels[ic] == "LATTNUM") {
	latnumcolumn = true;
	extracolumnlabels.push_back(columnlabels[ic]);
	ic++;
	continue;
      }  
      int nlatncol = CheckColumnlabelN(columnlabels[ic], "LATTNUM");
      if (nlatncol > 0) {
	multilatscheme = 2;   // scheme 2
	extracolumnlabels.push_back(columnlabels[ic]);
      }
      int nscalencol = CheckColumnlabelN(columnlabels[ic], "SCALE");
      if (nscalencol > 0) {
	extracolumnlabels.push_back(columnlabels[ic]);
      }
      std::pair<int, int> hkln = 
	CheckColumnlabelHKL(columnlabels.at(ic));
      if (hkln.first == 1) { // column"Hn" found
	// Sanity check: following columns should be Kn, Ln with same n
	//  and type = "H"
	int hn = hkln.second;  // n
	bool OK = true;
	if (columntypes.at(ic) != "H") {OK = false;} // fail
	hkln = CheckColumnlabelHKL(columnlabels.at(ic+1)); // Kn
	if (hkln.first != 2 || hkln.second != hn) {OK = false;} // fail
	if (columntypes.at(ic+1) != "H") {OK = false;} // fail
	hkln = CheckColumnlabelHKL(columnlabels.at(ic+2)); // Ln
	if (hkln.first != 3 || hkln.second != hn) {OK = false;} // fail
	if (columntypes.at(ic+2) != "H") {OK = false;} // fail
	if (!OK) {
	  Message::message
	    (Message_fatal
	     ("MtzUnmrgFile: CheckMultipleLattices - inconsistent HKL column labels or types" ) );
	}
	// Store extra column labels
	for (int i=0;i<3;++i) {
	  extracolumnlabels.push_back(columnlabels.at(ic+i));
	}
	nlatticecolumns++;
	ic += 2; // extra increment over 3 hkl labels
      } // end if Hn column
      ic++;
    } // end loop columns
    //^
    //    std::cout <<"CheckMultipleLattices, nlatticecolumns "
    //    	      << nlatticecolumns << "\n";;
    //    for (size_t j=0; j<extracolumnlabels.size(); j++) { 
    //      std::cout << " "<< extracolumnlabels[j];
    //    }
    //    std::cout <<"\n";
    //^-
  }
  //--------------------------------------------------------------
  // add extra lattice columns to column_label_list if required
  void  MtzUnmrgFile::add_extra_columns()
  {
    if (nlatticecolumns > 0) {
      for (size_t l=0;l<extracolumnlabels.size();++l) {
	// compulsory now we know it's there
	column_label_list.add(extracolumnlabels[l],OF_COMPULSORY);
      }
    }
    column_label_list.SetNlatticeColumns(nlatticecolumns);
  }
  //--------------------------------------------------------------
  //! return lattice number, = 0 for single lattice, = -1 for mixed lattices
  int MtzUnmrgFile::LatticeNumber() const
  {
    if (nlatticecolumns <= 0) return 0;
    int count = 0;
    int lat = 0;
    for (size_t i=0;i<numberinlattice.size();++i) {
      // for single lattice, only one of these should be > 0
      if (numberinlattice[i] > 0) {
	count++;
	lat = i;
      }
    }
    if (count == 1) return lat;
    return -1;
  }
  //--------------------------------------------------------------
  // Copy constructor throws exception unless object is EMPTY
  MtzUnmrgFile::MtzUnmrgFile(const MtzUnmrgFile& MUfile)
  {
    if (MUfile.mtzin != NULL)      {
      Message::message(Message_fatal
		       ("MtzUnmrgFile: illegal copy constructor"));
    }
    clear();
  }
  //--------------------------------------------------------------
  // Copy operator throws exception unless object is EMPTY
  MtzUnmrgFile& MtzUnmrgFile::operator= (const MtzUnmrgFile& MUfile)
  {
    if (MUfile.mtzin != NULL)      {
      Message::message(Message_fatal
		       ("MtzUnmrgFile: illegal copy operation"));
    }
    clear();
    return *this; 
  }
//--------------------------------------------------------------
} // namespace MtzIO 
