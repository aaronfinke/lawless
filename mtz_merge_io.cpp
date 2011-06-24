// mtz_merge_io.cpp
//
// Phil Evans 2009
//
// MTZ io for merged files (read only)
// 

#include "mtz_merge_io.hh"
#include "hkl_datatypes.hh"
#include "string_util.hh"

namespace MtzIO {
  //--------------------------------------------------------------
  std::string FixRhombohedralSymbolMess(const std::string& SGsymbol,
					const Scell& cell)
  {
    std::string Symbol = SGsymbol;
    int kR = Symbol.find("R");
    if (kR < int(Symbol.size())) {
      // Symbol contains "R"
      if (cell.AngleTest(90.,90.,120.))	{
	// Cell is hexagonal (angles 90,90,120), reset lattice symbol to "H"
	Symbol[kR] = 'H';
      }
    }
    return Symbol;
  }
  //--------------------------------------------------------------
  //! Constructor: does nothing
  MtzMrgFile::MtzMrgFile()
    : fileopen(false)
  {}
  //--------------------------------------------------------------
  MtzMrgFile::~MtzMrgFile()
  {}
  //--------------------------------------------------------------
  bool MtzMrgFile::open_read(const std::string& filename_in)
  // Open file for reading
  // returns false if fails
  {
    if (fileopen)
      Message::message( Message_fatal( "MtzMrgFile: open_read - File already open" ) );
    if ( filename_in == "") 
      Message::message( Message_fatal( "MtzMrgFile: open_read - no filename given" ) );

    // store filename
    if (getenv(filename_in.c_str()) != NULL) {
      filenamein = std::string(getenv(filename_in.c_str()));
    } else {
      filenamein = filename_in;
    }

    // open file
    clipper::CCP4MTZfile mtzin;
    try {
      mtzin.open_read(filenamein);
    }
    catch (Message_fatal) {
      // Failed to open file, missing or incomplete
      return false;
    }
    fileopen = true;

    mtzfile_resolution = mtzin.resolution();
    ResMax = mtzfile_resolution.limit();

    mcell = Scell(mtzin.cell());

    spacegroup.init(mtzin.spacegroup());
    spacegroupsymbol = spacegroup.Symbol_hm();

    // Clipper seems to return spacegroup R3 as "R3" even on hexagonal axes
    //  so change "R" to "H" if cell is hexagonal
    ///    spacegroupsymbol = FixRhombohedralSymbolMess
      ///      (mtzin.spacegroup().symbol_hm(), mcell);

    // Check for column label "M_ISYM" as marker for unmerged file
    merged = true;
    std::vector<clipper::String> ColLab = mtzin.column_labels();
    std::vector<std::vector<clipper::String> > LabelTypes(ColLab.size());
    bool dummy;
    ClipperLabelPair DummyColumns =
      ProcessLabels(ColLab, column_labels(), dummy);
    if (DummyColumns.path == "") {
      // File is unmerged
      merged = false;
    }
    mtzin.close_read();
    bool status = merged;
    return status;
  }
  //--------------------------------------------------------------
  FileRead MtzMrgFile::MakeHklList(const std::string& mtzname,
				   file_select& file_sel, 
				   col_controls& column_selection,
				   MtzIO::column_labels& column_list,
				   const scala::PxdName& InputPxdName,
				   const scala::Scell& cell,
				   std::string& output,
				   const int& verbose,
				   hkl_unmerge_list& hkl_list)
  // Fill an unmerged hkl_list from a merged file
  // On entry:
  //  mtzname            name of MTZ file (or logical name)
  //  file_sel           flags for general selection
  //                     - dataset selection
  //                     - batch selection
  //                     - resolution limits
  //                     - detector coordinate rejection ranges
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
  //
  // Returns FileRead.Opened() true if the MTZ file has been successfully opened
  // Returns FileRead.Read() true if the MTZ file has been read, false if has failed
  {
    if (!merged) {  // file must be merged for this function
      Message::message
	(Message_fatal("MtzMrgFile::MakeHklList: not a merged file"));
    }
    if (!hkl_list.IsEmpty()) {
      Message::message
	(Message_fatal("MtzMrgFile::MakeHklList: hkl_list is not empty"));
    }

    clipper::CCP4MTZfile mtzin;
    try {
      mtzin.open_read(filenamein);
    }
    catch (Message_fatal) {
      // Failed to open file, missing or incomplete
      return FileRead(false, false, false, 0);
    }
    output = "";
    if (verbose > 0) {
      output += FormatOutput::logTabPrintf(0,
	  "\nReflection list generated from merged file: %s\n",filenamein.c_str());
      output += FormatOutput::logTabPrintf(0,
		   "\nTitle: %s\n\n", mtzin.title().c_str());
      output += FormatOutput::logTabPrintf(0,
	      "   Space group from HKLIN file : %s\n",
		   mtzin.spacegroup().symbol_hm().c_str());
      output += FormatOutput::logTabPrintf(0, "   Cell: ");
      for (int i=0;i<6;i++) output += FormatOutput::logTabPrintf(0,"%7.2f",
				 mcell[i]);
      output += FormatOutput::logTab(0,"\n");
      output += FormatOutput::logTabPrintf(0,
			  "   Maximum resolution in file:  %8.2f\n",
			   mtzin.resolution().limit());
    }
    mtzfile_resolution = mtzin.resolution();
    // Read all data into Clipper objects
    clipper::ftype ResoLimit = file_sel.reslimits().ResHigh();
    clipper::HKL_info hkl_info_list;
    clipper::HKL_data<clipper::data32::I_sigI> IsigData;
    ClipperLabelPair labelthings = ReadData(mtzin, ResoLimit, column_list,
					    hkl_info_list, IsigData,
					    mtzdataset, (verbose>0), output);

    bool NoSigI = (labelthings.label2 == "");  // true if no sigma column

    // Now construct hkl_list
    std::string title;
    int Nref = hkl_info_list.num_reflections();

    // Construct dataset (clipper doesn't give us a project)
    std::vector<scala::Xdataset> DataSets;
    DataSets.push_back(scala::Xdataset(PxdName("", labelthings.xname,
					       labelthings.dname),
				       Scell(mtzin.cell()), mtzdataset.wavelength(), 1));
    // One batch
    std::vector<Batch> Batches(1);

    hkl_list.init(title, Nref,
		  hkl_symmetry(spacegroup), all_controls(),
		  ///		  hkl_symmetry(mtzin.spacegroup()), all_controls(),
		  DataSets, Batches);

    int isym = 1;
    int batch = 1;
    Rtype Xdet = 0.0;  Rtype Ydet = 0.0;
    Rtype phi = 0.0; Rtype time = 0.0;
    Rtype fraction_calc = 1.0;
    Rtype width = 0.0;
    Rtype LP = 0.0;
    int Npart = 1;
    int Ipart = 1;
    ObservationFlag ObsFlag;

    // Store all reflections
    // For resolution range found in file
    Range InvResRange;

    // Set start
    clipper::HKL_info::HKL_reference_index hkl_index = IsigData.first();
    at_start = true;
    IsigI Is;
    Rtype I, Ipr;
    Rtype sigI = 1.0;   // Dummy sigma = 1
    Rtype sigIpr = 1.0;

    while (next(hkl_index)) {
      Is = IsigData[hkl_index];
      if (!Is.missingI()) {  // I column OK
	if (Is.I() > 0.0) {
	  scala::Hkl hkl(hkl_index.hkl());  // hkl of current reflection
	  scala::Hkl hred = hkl_list.symmetry().put_in_asu(hkl, isym);
	  I = Is.I();
	  Ipr = Is.I();
	  if (!NoSigI) {
	    sigI = Is.sigI();
	    sigIpr = Is.sigI();
	  }
	  // Store this observation
	  hkl_list.store_part(hred, isym, batch, I, sigI, Ipr, sigIpr,
			      Xdet, Ydet, phi, time,
			      fraction_calc, width, LP,
			      Npart, Ipart, ObsFlag);
	  InvResRange.update( hkl_index.invresolsq());  //smin, smax
	}
      }
    }
    // 
    bool sorted = true;
    hkl_list.close_part_list(ResoRange(InvResRange), sorted);

    // Set into column_list the actual column numbers for columns requested by program
    //    (in column_list)
    //^!    column_list.get_col_lookup(*this);
    //^!    column_select col_select(column_list);

    // Store data presence flags
    //^!    hkl_list.StoreDataFlags(col_select.DataFlags());
    // Store file name
    hkl_list.AppendFileName(filenamein);
    return FileRead(true, true, true, 0);
  }
  //--------------------------------------------------------------
  ClipperLabelPair MtzMrgFile::ReadData(clipper::CCP4MTZfile& mtzin,
			    const double& ResoLimit,
			    const MtzIO::column_labels& column_list,
			    clipper::HKL_info& hkl_info_list,
		    clipper::HKL_data<clipper::data32::I_sigI>& IsigData,
			    clipper::MTZdataset& mtzdataset,
			    const bool& verbose,
			    std::string& output)
  // Read all selected data from MTZ file into clipper objects
  // hkl_info_list, IsigData, mtzdataset
  // Returns label things
  {
    ///    output = "";
    ResMax = mtzin.resolution().limit();
    if (ResoLimit > 0.0)
      ResMax = Max(ResoLimit, ResMax);
    
    // Set hkl list to desired resolution
    // reflections outside limits will be discarded
    hkl_info_list =
      ///      clipper::HKL_info(mtzin.spacegroup(), mtzin.cell(),
      clipper::HKL_info(spacegroup, mtzin.cell(),
			clipper::Resolution(ResMax));

    if (verbose) {
      if (ResoLimit > 0.0 && ResMax > mtzfile_resolution.limit()+0.001) {
	output += FormatOutput::logTabPrintf(0, "Maximum resolution in file %s: %8.3f",
					     filenamein.c_str(), mtzfile_resolution.limit());
	output += FormatOutput::logTabPrintf(0,"  restricted to %8.3f", ResMax);
	output += "\n";
      }
    }

    bool IorF;
    ClipperLabelPair labelthings
      = ProcessLabels(mtzin.column_labels(), column_list, IorF);

    clipper::String LabColI  = labelthings.label1;
    clipper::String LabColsigI  = labelthings.label2;
    bool NoSigI = (LabColsigI == "");  // true if there is no sigI column 

    if (verbose) {
      if (IorF) {
	// column found is F
	output += FormatOutput::logTab(1, "Columns for F, sigF (squared to I): "+
		      LabColI+"  "+LabColsigI+"\n");
      } else {
	// column found is I
	output += FormatOutput::logTab(1, "Columns for I, sigI: "+
		      LabColI+"  "+LabColsigI+"\n");
      }
    }

    // Clipper MTZ dataset
    mtzin.import_dataset(mtzdataset, labelthings.path);

    //......................................................
    // Read header info & hkl list
    mtzin.import_hkl_info(hkl_info_list);
    IsigData.init(hkl_info_list, hkl_info_list.cell());
    if (IorF) {
      // File contains F, square all values
      // I = F^2
      // sigI = 2 F sigF  + sigF^2
      clipper::HKL_data<clipper::data32::F_sigF> FsigData(hkl_info_list); // temporary for F
      // Read in data
      mtzin.import_hkl_data(FsigData, labelthings.path);
      mtzin.close_read();

      clipper::HKL_info::HKL_reference_index ih;
      clipper::data32::I_sigI Isig;
      double F;
      double sigF = 1.0;
      const double iscale = 0.01;  // scale down F^2

      for (ih = hkl_info_list.first(); !ih.last(); ih.next()) {
	// OK if
	// 1. NoSigI && F OK, or
	// 2. F OK
	if (!clipper::Util::is_null(FsigData[ih].f())) { // F not null
	  sigF = 1.0;
	  if (!NoSigI && !clipper::Util::is_null(FsigData[ih].sigf())) {
	    sigF = FsigData[ih].sigf();  // sigF if present and OK
	  }
	  F = FsigData[ih].f();
	  Isig.I() = F * F;
	  Isig.sigI() = 2.*F*sigF + sigF*sigF;
	  Isig.scale(iscale);	  
	  IsigData[ih] = Isig;
	}
      }
    } else {
      // Read in data
      mtzin.import_hkl_data(IsigData, labelthings.path);
      mtzin.close_read();
    }
    return labelthings;
  }
  //--------------------------------------------------------------
//--------------------------------------------------------------
  //--------------------------------------------------------------
  // Next IsigI, returns false if end of list
  bool MtzMrgFile::next(clipper::HKL_info::HKL_reference_index& hkl_index)
  {
    if (!at_start)
      // increment index
      hkl_index.next();
    at_start = false;
    if (hkl_index.last()) return false;
    return true;
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
} // namespace MtzIO

