// hkl_merged_list.cpp

#include "hkl_merged_list.hh"
#include "pointgroup.hh"

#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>

#include <assert.h>
#define ASSERT assert

namespace scala{
  //===================================================================
  hkl_merged_list::hkl_merged_list(const std::string& mtzname,
					 const bool& verbose,
					 phaser_io::Output& output)
    // Construct just the header information from the MTZ file
    // Do not create clipper objects yet (done in "read" method)
  {
    init(mtzname,verbose,output);
  }
  //------------------------------------------------------
  void hkl_merged_list::init(const std::string& mtzname,
			     const bool& verbose,
			     phaser_io::Output& output)
    // Initialise just the header information from the MTZ file
    // Do not create clipper objects yet (done in "read" method)
  {
    mtzfilein.open_read(mtzname); // file closed again after reading header

    // Retrieve information
    ResMax = mtzfilein.MtzResolution().limit();
    filename = mtzfilein.Filename();

    if (verbose) {
      output.logTab(0,LOGFILE, "\nHKL file "+filename+"\n");
      output.logTab(0,LOGFILE, "Spacegroup: "+
		    mtzfilein.Spacegroupsymbol()+"\n");
      output.logTabPrintf(0,LOGFILE, "Cell: ");
      for (int i=0;i<6;i++) {
	output.logTabPrintf(0,LOGFILE,"%7.2f",
			    mtzfilein.Cell().UnitCell()[i]);
      }
      output.logTab(0,LOGFILE,"\n");
    }

    status = MLIST::HEADER;
    if (!mtzfilein.Merged()) {
      status = MLIST::UNMERGED;
    }
  }
//--------------------------------------------------------------
  void hkl_merged_list::read(const MtzIO::column_labels& column_list,
			     const double& ResoLimit,
			     const bool& verbose,
			     phaser_io::Output& output)
  // Read data into arrays, construct clipper objects
  {
    if (status != MLIST::HEADER)
      Message::message(Message_fatal
		       ("hkl_merged_list: object not constructed"));

    clipper::CCP4MTZfile mtzin;
    mtzin.open_read(filename);

    clipper::MTZdataset mtzdataset;
    std::string outputstring;
    mtzfilein.ReadData(mtzin, ResoLimit, column_list, hkl_info_list,
		       IsigData, mtzdataset, verbose, outputstring);
    output.logTab(0,LOGFILE,outputstring);

    fcell = mtzfilein.Cell();
    // FIXME change R3 to H3
    fsymmetry = hkl_symmetry(mtzfilein.Spacegroupsymbol());

    status = MLIST::DATA;
    hkl_index = IsigData.first();
    at_start = true;
  }
//--------------------------------------------------------------
  IsigI hkl_merged_list::Isig(const Hkl& h) const
    // Return I sigI for given hkl
  {
    if (status != MLIST::DATA)
      Message::message(Message_fatal
		       ("hkl_merged_list: no data read"));

    IsigI IsI = IsigData[clipper::HKL(h.h(),h.k(),h.l())];
    if (IsI.missing()) {
      // no data, return 0,0
      return IsigI(0.0,0.0);
    } else {
      return IsigI(IsI);
    }
  }
//--------------------------------------------------------------
  hkl_symmetry hkl_merged_list::symmetry() const
  {
    // FIXME change R3 to H3
    return fsymmetry;
  }
//--------------------------------------------------------------
  std::string hkl_merged_list::SpaceGroupSymbol() const
  {
    return fsymmetry.symbol_xHM();
  }
//--------------------------------------------------------------
  int hkl_merged_list::num_obs() const
  {
    if (status == MLIST::DATA)
      return IsigData.num_obs();
    return 0;
  }
//--------------------------------------------------------------
  // Reset current reflection pointer to first reflection
  void hkl_merged_list::start() const
  {
    hkl_index = IsigData.first();
    at_start = true;
  }
  //--------------------------------------------------------------
  // Next IsigI, returns false if end of list
  bool hkl_merged_list::next(IsigI& Is) const
  {
    if (status != MLIST::DATA)
      Message::message(Message_fatal
		       ("hkl_merged_list: no data read"));
    if (!at_start)
      // increment index
      hkl_index.next();
    at_start = false;
    if (hkl_index.last()) return false;
    Is = IsigData[hkl_index];
    return true;
  }
  
  //--------------------------------------------------------------
  // get hkl for current reflection
  Hkl  hkl_merged_list::hkl() const
  {
    return Hkl(hkl_index.hkl());
  }
  //--------------------------------------------------------------
  // resolution of current reflection
  double hkl_merged_list::invresolsq() const
  {
    return hkl_index.invresolsq();
  }
  //--------------------------------------------------------------
// Copy constructor throws exception
  hkl_merged_list::hkl_merged_list(const hkl_merged_list& List)
  {
    Message::message(Message_fatal
		     ("hkl_merged_list: illegal copy constructor"));
  }
//--------------------------------------------------------------
  // Copy operator throws exception
  hkl_merged_list& hkl_merged_list::operator= (const hkl_merged_list& List)
  {
    Message::message(Message_fatal
		     ("hkl_merged_list: illegal copy operation"));
    return *this; // dummy
  }
//--------------------------------------------------------------
  void hkl_merged_list::PrintHeaderStuff(phaser_io::Output& output) const
  // Optional summary printing
  {
    //        output.logTabPrintf(0,LOGFILE,
    //			    "\nSummary of reflection list\n");
    output.logTabPrintf(0,LOGFILE,
			"   Highest resolution: %8.2f\n", ResMax);
    output.logTabPrintf(1,LOGFILE, "Unit cell: ");
    output.logTab(0,LOGFILE,fcell.formatPrint());
    output.logTabPrintf(1,LOGFILE,"");
    output.logTab(1,LOGFILE,"Space group: "+SpaceGroupSymbol());
  }
  //--------------------------------------------------------------
  // Initialise list from structure factor calculation from XYZIN coordinate file
  //  xyzin          file name for coordinate file
  //  spacegroup     spacegroup given on input to override that in xyzin file
  //                 blank if not given
  //  input_cell     cell if given on input
  //  FCresolution   maximum resolution to generate
  void hkl_merged_list::CreateFromAtoms(const std::string& xyzin,
					const std::string& spacegroup,
					const Scell& input_cell,
					const double& FCresolution,
					const bool& verbose,
					phaser_io::Output& output)
  {
    clipper::HKL_data<F_phi> fc;
    // fill hkl_info_list and fc
    SFcalc(xyzin, spacegroup, input_cell, FCresolution, hkl_info_list, fc);
    // Extract Fs and store as intensities
    //    clipper::HKL_data<clipper::data32::I_sigI> IsigData;
    IsigData.init(hkl_info_list, hkl_info_list.cell());
    const clipper::ftype32 SIGI = 1.0;
    int n = 0;
    clipper::HKL_info::HKL_reference_index ih;

    for (ih = hkl_info_list.first(); !ih.last(); ih.next()) {
      double f2 = fc[ih].f()*fc[ih].f();  // square F to I
      IsigData[ih] = clipper::data32::I_sigI(f2, SIGI);
      n++;
    }

    if (verbose) {
      output.logTab(0,LOGFILE,
		    "\nReference list generated by structure factor calculation\n  from coordinate file "+xyzin+"\n");
      output.logTab(0,LOGFILE, "Spacegroup: "+
		    hkl_info_list.spacegroup().descr().symbol_hm()+"\n");
      output.logTabPrintf(0,LOGFILE, "Cell: ");
      output.logTabPrintf(0,LOGFILE,"%7.2f", hkl_info_list.cell().descr().a());
      output.logTabPrintf(0,LOGFILE,"%7.2f", hkl_info_list.cell().descr().b());
      output.logTabPrintf(0,LOGFILE,"%7.2f", hkl_info_list.cell().descr().c());
      output.logTabPrintf(0,LOGFILE,"%7.2f", hkl_info_list.cell().descr().alpha_deg());
      output.logTabPrintf(0,LOGFILE,"%7.2f", hkl_info_list.cell().descr().beta_deg());
      output.logTabPrintf(0,LOGFILE,"%7.2f\n", hkl_info_list.cell().descr().gamma_deg());
      output.logTabPrintf(0,LOGFILE, "Maximum resolution used: %7.3f\n",
			  FCresolution);
      output.logTabPrintf(0,LOGFILE, "Number of reflections: %9d\n\n", n);
    }

    fcell = hkl_info_list.cell();
    fsymmetry = hkl_symmetry(hkl_info_list.spacegroup());
    ResMax = FCresolution;

    status = MLIST::DATA;
    hkl_index = IsigData.first();
    at_start = true;
  }
  //--------------------------------------------------------------
  void hkl_merged_list::SFcalc(const std::string& xyzin,
			       const std::string& spacegroup,
			       const Scell& input_cell,
			       const clipper::ftype64& resolution,
			       clipper::HKL_info& hkls,
			       clipper::HKL_data<F_phi>& fc)
  // Calculate structure factors
  //
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
  {
    // atomic model
    // Read file
    clipper::MMDBfile mfile;
    clipper::MiniMol  mmol;
    mfile.read_file(xyzin);
    mfile.import_minimol(mmol);
    
    // resolution
    clipper::Resolution Reso(resolution);
    
    // Get space group
    clipper::Spacegroup SG;
    if (spacegroup != "") {
      SG = scala::SpaceGroup(spacegroup);
      if (SG.is_null()) {
	Message::message(Message_fatal
			 ("Invalid SPACEGROUP given: "+spacegroup));
      }
    } else {
	SG.init(mmol.spacegroup());
	if (SG.is_null()) {
	  Message::message(Message_fatal
			   ("No valid space group in XYZIN file"));
	}
    }
    // Get cell
    clipper::Cell cell;
    if (input_cell.null()) {
      cell = mmol.cell();
      if (cell.is_null()) {
	Message::message(Message_fatal
			 ("No valid unit cell in XYZIN file"));
      }
    } else {
      cell = input_cell.ClipperCell();
    }

    // Generate hkl list
    hkls.init(SG, cell, Reso, true);
    
    // calculate structure factors
    fc.init( hkls, cell);
    clipper::SFcalc_iso_fft<float> sfc(fc, mmol.model().atom_list());
  }
}



