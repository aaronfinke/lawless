// hkl_merge.cpp
// a cut-down version of hkl_merged_list, just merged MTZ input

#include "hkl_merge.hh"
#include "pointgroup.hh"
#include "string_util.hh"
#include "hkl_datatypes.hh"
#include "mtz_utils.hh"
#include "cellgroup.hh"

#include <clipper/clipper.h>
#include "clipper/clipper-ccp4.h"
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>

#include <assert.h>
#define ASSERT assert

using phaser_io::LXML;

namespace scala {
  //===================================================================
  hkl_merge::hkl_merge(const std::string& hklinname,
		       const bool& verbose,
		       phaser_io::Output& output)
    // Construct just the header information from the MTZ file
    // Do not create clipper objects yet (done in "read" method)
  {
    init(hklinname,verbose,output);
  }
  //------------------------------------------------------
  void hkl_merge::init(const std::string& hklinname,
		       const bool& verbose,
		       phaser_io::Output& output)
  // Initialise just the header information from the MTZ file
  // Do not create clipper objects yet (done in "read" method)
  {
    // MTZ file, file closed again after reading header
    bool fileopen = mtzfilein.open_read(hklinname);
    if (!fileopen) {
      Message::message(Message_fatal
		       ("hkl_merge: open_read - not valid MTZ file"));
    }

    // Retrieve information
    ResMax = mtzfilein.MtzResolution().limit();
    filename = mtzfilein.Filename();
      
    if (verbose) {
      output.logTab(1,LOGFILE, "File name:"+filename+"\n");
      output.logTab(1,LOGFILE, "Spacegroup: "+
		    mtzfilein.Spacegroupsymbol()+"\n");
	output.logTabPrintf(1,LOGFILE, "Cell: ");
	for (int i=0;i<6;i++) {
	  output.logTabPrintf(0,LOGFILE,"%7.2f",
			      mtzfilein.Cell().UnitCell()[i]);
	}
	output.logTab(0,LOGFILE,"\n");
    }
    
    status = MLIST::HEADER;
    if (!mtzfilein.Merged()) {
      Message::message(Message_fatal
		       ("hkl_merge: input file must be merged"));
    }
  }
  //--------------------------------------------------------------
  void hkl_merge::read(const MtzIO::column_labels& column_list,
			     const double& ResoLimit,
			     const bool& verbose,
			     phaser_io::Output& output)
  // Read data into arrays, construct clipper objects
  {
    if (status != MLIST::HEADER)
      Message::message(Message_fatal
		       ("hkl_merge: object not constructed"));

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

    // reset resolution maximum
    ResMax = mtzfilein.resHigh();

    status = MLIST::DATA;
    hkl_index = IsigData.first();
    at_start = true;
  }
  //--------------------------------------------------------------
  IsigI hkl_merge::Isig(const Hkl& h) const
    // Return I sigI for given hkl
  {
    if (status != MLIST::DATA)
      Message::message(Message_fatal
		       ("hkl_merge: no data read"));
    if (hkl_info_list.index_of(h.HKL()) >= 0) {
      IsigI IsI = IsigData[h.HKL()];
      if (IsI.missing()) {
	// no data, return 0,-1
	return IsigI(0.0,-1.0);
      } else {
	return IsigI(IsI);
      }
    } else {
      // no data, return 0,-1
      return IsigI(0.0,-1.0);
    }
  }
//--------------------------------------------------------------
  hkl_symmetry hkl_merge::symmetry() const
  {
    // FIXME change R3 to H3
    return fsymmetry;
  }
//--------------------------------------------------------------
  std::string hkl_merge::SpaceGroupSymbol() const
  {
    return fsymmetry.symbol_xHM();
  }
//--------------------------------------------------------------
  int hkl_merge::num_obs() const
  {
    if (status == MLIST::DATA)
      return IsigData.num_obs();
    return 0;
  }
//--------------------------------------------------------------
  // Reset current reflection pointer to first reflection
  void hkl_merge::start() const
  {
    hkl_index = IsigData.first();
    at_start = true;
  }
  //--------------------------------------------------------------
  // Next IsigI, returns false if end of list
  bool hkl_merge::next(IsigI& Is) const
  {
    if (status != MLIST::DATA)
      Message::message(Message_fatal
		       ("hkl_merge: no data read"));
    if (!at_start) {
      // increment index
      hkl_index.next();
    }
    at_start = false;
    if (hkl_index.last()) return false;
    Is = IsigData[hkl_index];
    return true;
  }
  //--------------------------------------------------------------
  // get hkl for current reflection
  Hkl  hkl_merge::hkl() const
  {
    return Hkl(hkl_index.hkl());
  }
  //--------------------------------------------------------------
  // resolution of current reflection
  double hkl_merge::invresolsq() const
  {
    return hkl_index.invresolsq();
  }
  //--------------------------------------------------------------
// Copy constructor throws exception
  hkl_merge::hkl_merge(const hkl_merge& List)
  {
    Message::message(Message_fatal
		     ("hkl_merge: illegal copy constructor"));
  }
//--------------------------------------------------------------
  // Copy operator throws exception
  hkl_merge& hkl_merge::operator= (const hkl_merge& List)
  {
    Message::message(Message_fatal
		     ("hkl_merge: illegal copy operation"));
    return *this; // dummy
  }
//--------------------------------------------------------------
  void hkl_merge::PrintHeaderStuff(const bool& isreference,
					 phaser_io::Output& output) const
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

    if (! isreference) {
      // XML stuff
      std::string datatag = "ReflectionData";
      output.logTab(0,LXML,"<"+datatag+">");
      output.logTab(1,LXML,StringUtil::MakeXMLtag("MergedData","true"));
      output.logTab(1,LXML,
		    StringUtil::MakeXMLtag("ResolutionHigh",
					   StringUtil::ftos(ResMax,8,2)));
      output.logTab(1,LXML,
		    StringUtil::MakeXMLtag("NumberReflections",
					   StringUtil::itos(IsigData.num_obs(),8)));
      output.logTab(0,LXML,"</"+datatag+">");
    }
  }
  // ---------------------------------------------------------
  clipper::Spacegroup hkl_merge::spacegroup() const
  {
    return hkl_info_list.spacegroup();
  }
  // ---------------------------------------------------------
  clipper::Cell hkl_merge::cell() const
  {
    return hkl_info_list.cell();    
  }
  // ---------------------------------------------------------
  //--------------------------------------------------------------
  // Initialise list from structure factor calculation from XYZIN coordinate file
  //  xyzin          file name for coordinate file
  //  spacegroup     spacegroup given on input to override that in xyzin file
  //                 blank if not given
  //  input_cell     cell if given on input
  //  FCresolution   maximum resolution to generate
  void hkl_merge::CreateFromAtoms(const std::string& xyzin,
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
      if (!fc[ih].missing()) {
	double f2 = fc[ih].f()*fc[ih].f();  // square F to I
	IsigData[ih] = clipper::data32::I_sigI(f2, SIGI);
	n++;
      }
    }

    if (verbose) {
      output.logTab(0,LOGFILE,
		    "\nReference list generated by structure factor calculation\n  from coordinate file "+xyzin+"\n");
      output.logTab(0,LXML,
		    StringUtil::MakeXMLtag("XYZREF", xyzin));
      output.logTab(0,LOGFILE, "Spacegroup: "+
		    hkl_info_list.spacegroup().descr().symbol_hm()+"\n");
      output.logTab(0,LXML,
		    StringUtil::MakeXMLtag("XYZREFspacegroup",
					   hkl_info_list.spacegroup().descr().symbol_hm()));
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
  void hkl_merge::SFcalc(const std::string& xyzin,
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

    std::vector<clipper::Atom> atoms;
    // Strip zero occupancy atoms
    for (size_t i=0;i<mmol.model().atom_list().size();++i) {
      clipper::Atom atom = mmol.model().atom_list()[i];
      if (atom.occupancy() > 0.0) {
	atoms.push_back(atom);
      }
    }

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
	SG = mmol.spacegroup();
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
    clipper::SFcalc_iso_fft<float> sfc(fc, atoms);
  }
  //--------------------------------------------------------------
  // Initialise list from structure factor calculation from XYZIN coordinate file
  //  with bulk solvent scaled to Fobs list
  //  xyzin          file name for coordinate file
  //  
  void hkl_merge::createFromAtomsBulk(const std::string& xyzin,
				      const MergedList& mergedobslist,
				      const bool& verbose,
				      phaser_io::Output& output)
  {
    const clipper::HKL_data<clipper::data32::I_sigI>& isigi =
      mergedobslist.ImeanForDataset(0);

    hkl_info_list = isigi.base_hkl_info();

    clipper::HKL_data<clipper::data32::F_sigF> fobs(hkl_info_list);
    typedef clipper::HKL_data_base::HKL_reference_index HRI;
    for ( HRI ih = isigi.first(); !ih.last(); ih.next() ) {
      if (isigi[ih].sigI() > 0.0) {
	float I = isigi[ih].I();
	float sigI = isigi[ih].sigI();
	if (I >= 0.0) { // omit negatives
	  fobs[ih].f() = sqrt(I);
	  fobs[ih].sigf() = sqrt(sigI + I) - fobs[ih].f();
	}
      }
    }

    clipper::HKL_data<F_phi> fc(hkl_info_list);
    // fill hkl_info_list and fc
    SFcalcBulk(xyzin, fobs, fc, output);
    // Extract Fs and store as intensities
    //    clipper::HKL_data<clipper::data32::I_sigI> IsigData;
    IsigData.init(hkl_info_list, hkl_info_list.cell());
    const clipper::ftype32 SIGI = 1.0;
    int n = 0;
    clipper::HKL_info::HKL_reference_index ih;

    for (ih = hkl_info_list.first(); !ih.last(); ih.next()) {
      if (!fc[ih].missing()) {
	double f2 = fc[ih].f()*fc[ih].f();  // square F to I
	IsigData[ih] = clipper::data32::I_sigI(f2, SIGI);
	n++;
      }
    }
    
    double FCresolution = IsigData.base_hkl_info().resolution().limit();

    if (verbose) {
      output.logTab(0,LOGFILE,
		    "\nReference list generated by structure factor calculation\n  from coordinate file "+xyzin+"\n");
      output.logTab(0,LXML,
		    StringUtil::MakeXMLtag("XYZREF", xyzin));
      output.logTab(0,LOGFILE, "Spacegroup: "+
		    hkl_info_list.spacegroup().descr().symbol_hm()+"\n");
      output.logTab(0,LXML,
		    StringUtil::MakeXMLtag("XYZREFspacegroup",
					   hkl_info_list.spacegroup().descr().symbol_hm()));
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
  void hkl_merge::SFcalcBulk(const std::string& xyzin,
	     const clipper::HKL_data<clipper::data32::F_sigF> fobs,
			     clipper::HKL_data<F_phi>& fc,
			     phaser_io::Output& output)
  // Calculate structure factors
  //
  // On input:
  //  xyzin         file name for coordinate file
  //  fobs          fobs list
  //
  // On exit:
  //  fc            Fcalc list
  {
    // atomic model
    // Read file
    clipper::MMDBfile mfile;
    clipper::MiniMol  mmol;
    mfile.read_file(xyzin);
    mfile.import_minimol(mmol);

    std::vector<clipper::Atom> atoms;
    // Strip zero occupancy atoms
    for (size_t i=0;i<mmol.model().atom_list().size();++i) {
      clipper::Atom atom = mmol.model().atom_list()[i];
      if (atom.occupancy() > 0.0) {
	atoms.push_back(atom);
      }
    }

    // calculate structure factors
    clipper::SFcalc_obs_bulk<float> sfc(fc, fobs, atoms);

    double bulkfrc = sfc.bulk_frac();
    double bulkscl = sfc.bulk_scale();
  }
}   // namespace scala




