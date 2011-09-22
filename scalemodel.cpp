// scalemodel.cpp


#include "scalemodel.hh"
#include "util.hh"
#include "keywords_aimless.hh"
#include "jiffy.hh"
#include "tie.hh"
#include "string_util.hh"
#include "file_util.hh"
#include "fileread.hh"
#include "scala_util.hh"

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;

#include <assert.h>
#define ASSERT assert

using phaser_io::itos;
using phaser_io::ftos;
using phaser_io::LOGFILE;

namespace scala {
  //--------------------------------------------------------------
  ScaleModel::ScaleModel(const phaser_io::InputAll& input,
			 hkl_unmerge_list& hkl_list,
			 phaser_io::Output& output)
  // Construct from input commands and reflection list
  {
    init (input, hkl_list, output);
  }
  //--------------------------------------------------------------
  void ScaleModel::init(const phaser_io::InputAll& input,
			hkl_unmerge_list& hkl_list,
			phaser_io::Output& output)
  // initialise from input commands and reflection list
  // Scale specification(s) from input
  {
    // Setup scale model
    pole = 0;
    setup(input.getScaleSpecifications(), hkl_list, output);
    // Scale normalisation  ..............................
    // FIXME set up scale normalisation flags from input if necessary
    scalenormrun = -1;
    bfacnormrun = -1;
    std::vector<Run> runlist = hkl_list.RunList();

    for (int irun=0;irun<nruns;irun++) {
      if (scalenormrun < 0) {
	if (scalenormbatch >= 0) {
	  // Normalisation batch specified, is it in this run?
	  if (runlist[irun].IsInList(scalenormbatch)) {
	    // Yes, locate it
	    std::vector<int> batchnums = runlist[irun].BatchList(true); // accepted batches
	    for (size_t i=0;i<batchnums.size();++i) {
	      if (batchnums[i] == scalenormbatch) {
		scalenormbatch = i;  // index
		scalenormrun = irun; // run
		break;
	      }
	    }
	  }
	} else {
	  //  1st one
	  scalenormbatch = 0;
	  scalenormrun = irun;
	  // If smoothed scales && more than one, use second one
	  if (!primary_scales[scalenormrun].IsBatchScale() &&
	      primary_scales[scalenormrun].Number() > 2) {
	    scalenormbatch = 1;
	  }
	}
      }
      // Bfactors
      if (bfacnormbatch >= 0) {
	if (bfacnormrun < 0) {
	  // Normalisation batch specified, is it in this run?
	  if (runlist[irun].IsInList(bfacnormbatch)) {
	    // Yes, locate it (index j)
	    std::vector<int> batchnums = runlist[irun].BatchList(true);
	    for (size_t i=0;i<batchnums.size();++i) {
	      if (batchnums[i] == bfacnormbatch) {
		bfacnormbatch = i;
		bfacnormrun = irun;
		break;
	      }
	    }
	  }
	} 
	// else leave bfacnormbatch = -1, normalise on "best" batch later
      }
    }    // normalisation .............................

    // Ties
    SetupTies(input, hkl_list);
    // Always calculate all secondary beams & diffraction vectors
    hkl_list.CalcSecondaryBeams(pole);
  }
  //--------------------------------------------------------------
  void ScaleModel:: SetConstant(hkl_unmerge_list& hkl_list, phaser_io::Output& output)

  // set SCALE CONSTANT for all runs
  {
    scala::ScaleSpecification spec;
    spec.SetConstant();
    std::vector<scala::ScaleSpecification> scaleSpecs(1, spec);
    setup(scaleSpecs, hkl_list, output);
  }
  //--------------------------------------------------------------
  void ScaleModel::setup(const std::vector<scala::ScaleSpecification>& scaleSpecs,
			 hkl_unmerge_list& hkl_list,
			 phaser_io::Output& output)
  // Setup from scale specifications and reflection list
  // Sets pole, for ABSORPTION, = 1,2,3 for h,k,l, = -1 unspecified, = 0 SECONDARY
  {
    scalenormbatch= -1; // batch number for scale  normalisation, -1 for 1st
    bfacnormbatch = -1;   // batch number for B-factor  normalisation, -1 for best
    nruns = hkl_list.num_runs();
    std::vector<Run> runlist = hkl_list.RunList();
    // One primary scale object for each run
    primary_scales.resize(nruns);
    //  ... and one B-factor object per run
    //       (which may be null, ie nbfac == 0)
    relative_bfactors.resize(nruns);

    //^+
    //    std::cout << "\nScaleSpecs\n";
    //    int nspecs = scaleSpecs.size();
    //    for (int i=0;i<nspecs;i++) {
    //      scaleSpecs[i].dump();
    //    }
    //^-
    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    // Set up LINKs
    // By default, all runs belonging to the same dataset use the same
    // secondary scale
    sec_scale_index_run.assign(nruns,-1);
    int kscidx = -1;
    int j0 = -1;  // secondary
    std::vector<int> runDset(nruns,-1);
    int pole = 0; // for ABSORPTION, = 1,2,3 for h,k,l, = -1 unspecified, = 0 SECONDARY
    // Tile stuff
    int ktlidx = -1;
    int k0 = -1;  // tile
    int ndet0 = -1; // number of detectors
    detectortypes.resize(nruns);
    detector_scale_index_run.assign(nruns,-1);
    runnumbers.resize(nruns);

    // Check reflection list (batches for each run) to see what sort of scales would be valid
    ValidScaleModel validscalemodel(hkl_list);

    for (int irun=0;irun<nruns;++irun) { // loop runs for initial checks
      // Store run number
      runnumbers[irun] = runlist[irun].RunNumber();
      // which spec belongs to this run?
      int isp = scaleSpecIndex(irun, scaleSpecs, runlist);
      if (isp < 0) {
	// shouldn't happen
	Message::message(Message_fatal
			 ("ScaleModel: no scale specification for run "+
			  clipper::String(irun)));
      }

      // Secondary checks
      if (scaleSpecs[isp].sec_abs != scala::SecondaryScale::NONE &&
	  validscalemodel.ValidSecondary(irun)) {
	if (j0 < 0) {
	  // first run with a secondary correction
	  j0 = irun;  // 
	  sec_scale_index_run[j0] = ++kscidx;  // index for 1st run = 0
	  runDset[j0] = runlist[j0].DatasetIndex();
	  if (scaleSpecs[isp].sec_abs == SecondaryScale::ABSORPTION) {
	    //  ABSORPTION, store pole
	    pole = scaleSpecs[isp].pole;
	  } else {
	    pole = 0;  // SECONDARY
	  }
	} else {
	  // Check for all runs ABSORPTION or SECONDARY, not mixed
	  bool OK = true;
	  if (scaleSpecs[isp].sec_abs == SecondaryScale::ABSORPTION) {
	    if (pole == 0) {OK = false;}
	  } else {
	    // SECONDARY
	    if (pole != 0) {OK = false;}
	  }
	  if (!OK) {
	    Message::message(Message_fatal
			     ("You cannot mix SECONDARY and ABSORPTION"));
	  }
	  // check each run beyond 1st for same dataset
	  int idset = runlist[irun].DatasetIndex();
	  runDset[irun] =  idset;
	  // search previous runs for same dataset
	  bool found = false;
	  for (int jrun=0;jrun<irun;jrun++) {
	    if (runDset[jrun] == idset) {
	      // irun is same dataset as jrun
	      sec_scale_index_run[irun] = sec_scale_index_run[jrun];
	      found = true;
	      break;
	    }
	  }
	  if (!found) {
	    // irun is new dataset
	    sec_scale_index_run[irun] = ++kscidx;
	  }
	}
      } // end secondary checks
      // Detector (tile) checks
      if (scaleSpecs[isp].ntilex >= 0) {
	// Pick up detector size for 1st batch in this run
	int b0 = runlist[irun].BatchSerial0();  // batch serial 
	Batch bat0 = hkl_list.Batches()[b0];    // first batch in run
	detectortypes[irun] = DetectorType(bat0);

	if (k0 < 0) { // first run with tile correction
	  k0 = irun;
	  detector_scale_index_run[k0] = ++ktlidx;
	} else { // not first run
	  // check each run beyond 1st for same detector
	  // search previous runs for same detector
	  bool found = false;
	  for (int jrun=0;jrun<irun;jrun++) {
	    if (detectortypes[jrun] == detectortypes[k0]) {
	      // irun is same detector as jrun
	      detector_scale_index_run[irun] = detector_scale_index_run[jrun];
	      found = true;
	      break;
	    }
	  }
	  if (!found) {
	    // irun is new dataset
	    detector_scale_index_run[irun] = ++ktlidx;
	  }
	}
      } // end detector check
    }  // end run loop for initial checks
    nsecscales  = ++kscidx;
    ndetscales = ++ktlidx;
    // Process [UN]LINK commands
    //  ... when I've written them  FIXME
    //    (re)define sec_scale_index_run & nsecscales

    // make list of secondary scales, initialised to null 
    if (nsecscales > 0) {secondary_scales.assign(nsecscales, SecondaryScale());}
    // make list of detector scales, initialised to null 
    if (ndetscales > 0) {detector_scales.assign(ndetscales, DetectorScale());}
    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    // Set up runs
    // First set up runs which are explicitly specified
    std::vector<bool> runSetup(nruns, false);
    //  1st specification is always the general one, so search
    //  the others, if any
    for (size_t isp=1;isp<scaleSpecs.size();isp++) {
      int runNumber = scaleSpecs[isp].run;
      if (runNumber > 0) {
	int irun = FindRunIndex(runNumber, runlist);
	if (irun >= 0) {
	  // Specification for this run
	  std::string s = SetupScale(irun, scaleSpecs[isp], runlist[irun], validscalemodel);
	  if (s != "") {output.logTab(0,LOGFILE,s);}
	  runSetup[irun] = true;
	}
      }
    }
    // Then fill in the rest from the default specs
    for (int irun=0;irun<nruns;irun++) {
      if (!runSetup[irun]) {
	// Specification for this run
	std::string s = SetupScale(irun, scaleSpecs[0], runlist[irun], validscalemodel);
	if (s != "") {output.logTab(0,LOGFILE,s);}
	runSetup[irun] = true;
      }
    }
    // Count parameters, set nparameters & idxrun indices
    CountParameters();
  }
  //--------------------------------------------------------------
  void ScaleModel::SetupTies(const phaser_io::InputAll& input,
			     const hkl_unmerge_list& hkl_list)
  // Set up all ties from input list and defaults
  {
    sd_rotation = input.TIE_sd_rotation();
    sd_bfactor = input.TIE_sd_bfactor();
    sd_zerob = input.TIE_sd_zerob();
    sd_surface = input.TIE_sd_surface();
    sd_tile = input.TIE_sd_tile();

    ties.clear();  // clear tie list
    nties_rot = nties_bfac = nties_zerob = nties_surf = nties_tiles = 0;

    // Ties on primary scales
    if (sd_rotation > 0.0) {
      //  ... within each run
      for (int irun=0;irun<nruns;irun++) {
	// Set up ties for this run, using 1st global parameter index
	//  idxrun_primary_scales, append these to total list
	std::vector<Tie> these_ties = primary_scales[irun].
	  Ties(sd_rotation, idxrun_primary_scales[irun]);
	ties.insert(ties.end(), these_ties.begin(), these_ties.end());
	nties_rot += these_ties.size();
      }
    }
    // Ties on B-factors
    if (sd_bfactor > 0.0) {
      //  ... within each run
      for (int irun=0;irun<nruns;irun++) {
	// Set up ties for this run, using 1st global parameter index
	//  idxrun_bfactors, append these to total list
	std::vector<Tie> these_ties = relative_bfactors[irun].
	  Ties(sd_bfactor, idxrun_bfactors[irun]);
	ties.insert(ties.end(), these_ties.begin(), these_ties.end());
	nties_bfac += these_ties.size();
      }
    }
    // Ties on B-factors to zero
    if (sd_zerob > 0.0) {
      //  ... within each run
      for (int irun=0;irun<nruns;irun++) {
	// Set up ties for this run, using 1st global parameter index
	//  idxrun_bfactors, append these to total list
	std::vector<Tie> these_ties = relative_bfactors[irun].
	  ZeroTies(sd_zerob, idxrun_bfactors[irun]);
	ties.insert(ties.end(), these_ties.begin(), these_ties.end());
	nties_zerob += these_ties.size();
      }
    }
    // Ties on secondary beam scales
    if (sd_surface > 0.0) {
      //  ... for each surface
      if (nsecscales > 0) {
	for (int i=0;i<nsecscales;++i) {
	  std::vector<Tie> these_ties = secondary_scales[i].
	    Ties(sd_surface, idxrun_secondary[i]);
	  ties.insert(ties.end(), these_ties.begin(), these_ties.end());
	  nties_surf += these_ties.size();
	}
      }
    }
    // Ties on detector (tiles)
    if (ndetscales > 0) {
      if (sd_tile.size() > 0 && sd_tile[0] > 0.0) {
	for (int i=0;i<ndetscales;++i) {
	  std::vector<Tie> these_ties = detector_scales[i].
	    Ties(sd_tile, idxrun_detector[i]);
	  ties.insert(ties.end(), these_ties.begin(), these_ties.end());
	  nties_tiles += these_ties.size();
	}
      }
    }

    nties = ties.size();
  }
  //--------------------------------------------------------------
  std::string ScaleModel::ScaleParameterTypeString(const ScaleParameterType& type)
  // types of parameter
  //  enum ScaleParameterType {NONE, SCALE, BFACTOR, SECONDARY, TILE};
  // String representation
  {
    switch (type) {
    case NONE:
      return "NONE";
    case SCALE:
      return "SCALE";
    case BFACTOR:
      return "BFACTOR";
    case SECONDARY:
      return "SECONDARY";
    case TILE:
      return "TILE";
    default:
      return "UNKNOWN";
    }
  } 
  //--------------------------------------------------------------
  std::string ScaleModel::SetupScale(const int& irun,
			      const ScaleSpecification& scaleSpec,
			      const Run& run,
			      const ValidScaleModel&  validscalemodel)
  // Setup scales & B-factors for run irun
  // returns any warning messages
  {
    std::string s = ""; // warning messages
    bool batchscale = scaleSpec.batch; // true if BATCH mode
    // Check for valid smooth scaling if requested, fail is not
    if (scaleSpec.nscales > 1 ||
	scaleSpec.spacing > 0.0 ||
	scaleSpec.nbfac > 1 ||
	scaleSpec.bspacing > 0.0) {
      if (!batchscale) { // smooth mode
	if (!validscalemodel.ValidPrimary(irun)) {
	  Message::message(Message_fatal
			   ("\nERROR in ScaleModel: run "+clipper::String(run.RunNumber())+
			    " has insuffient information for smooth scaling\n"));
	}
      } else if (batchscale) {
	if (!validscalemodel.ValidBatch(irun)) { 
	  Message::message(Message_fatal
			   ("\nERROR in ScaleModel: run "+clipper::String(run.RunNumber())+
			    " has no useful batch information for BATCH scaling\n"));
	}
      }
    } // end check

    // Scales
    if (scaleSpec.nscales != 0) {
      if (batchscale) {
	// BATCH mode
	primary_scales[irun] =
	  PrimaryScale(run.BatchList(true));
      } else {
	// ROTATION mode
	if (scaleSpec.nscales > 0) {
	  // number given, number of intervals is one less
	  int NscaleIntervals = scaleSpec.nscales - 1;
	  primary_scales[irun] =
	    PrimaryScale(NscaleIntervals, run.PhiRange());
	} else if (scaleSpec.nscales < 0) {
	  // spacing given
	  primary_scales[irun] =
	    PrimaryScale(scaleSpec.spacing,
			 run.PhiRange());
	}
      }
    }
    // B-factors
    if (scaleSpec.nbfac != 0) {
      if (batchscale) {
	// BATCH mode
	relative_bfactors[irun] =
	  RelativeBfactor(run.BatchList(true));
      } else {
	// BROTATION mode
	if (scaleSpec.nbfac > 0) {
	  // number given, number of intervals is one less
	  int NbfacIntervals = scaleSpec.nbfac - 1;
	  relative_bfactors[irun] =
	    RelativeBfactor(NbfacIntervals, run.TimeRange());
	} else if (scaleSpec.nbfac < 0) {
	  // spacing given
	  relative_bfactors[irun] =
	    RelativeBfactor(scaleSpec.bspacing,
			    run.TimeRange());
	}
      }
    }
    // Secondary scales
    if (scaleSpec.sec_abs != SecondaryScale::NONE) {
      // Check for valid data: note that if primary data is missing, then so must be secondary
      if (!validscalemodel.ValidSecondary(irun)) {
	  s += FormatOutput::logTabPrintf(0,
     "WARNING: Run %3d has insuffient information for secondary beam scaling\n",
					  run.RunNumber());
      } else {
	// Which secondary scale object corresponds to this run?
	int jsc = sec_scale_index_run.at(irun);
	// Is this one already set up?
	if (secondary_scales[jsc].Type() == SecondaryScale::NONE) {
	  // No, so do it
	  secondary_scales[jsc] = SecondaryScale(scaleSpec.sec_abs,
						 scaleSpec.lmax, scaleSpec.lmaxodd,
						 scaleSpec.pole);
	}
      }
    }
    // Tiles, detector
    if (scaleSpec.detectorscaletype != DetectorScale::NONE) {
      if (!validscalemodel.ValidTile(irun)) {
	  s += FormatOutput::logTabPrintf(0,
     "WARNING: Run %3d has insuffient information for detector (tile) scaling\n",
					  run.RunNumber());
      } else {
	int jsc = detector_scale_index_run.at(irun);
	// Is this one set up?
	if (detector_scales[jsc].Type() ==
	    DetectorScale::NONE) { // if not, do it
	  DetectorScale detscale(scaleSpec.detectorscaletype,
				 scaleSpec.ntilex, scaleSpec.ntiley,
				 detectortypes[irun]);
	  if (detscale.Valid()) {
	    detector_scales[jsc].init(scaleSpec.detectorscaletype,
				      scaleSpec.ntilex, scaleSpec.ntiley,
				      detectortypes[irun]);
	  } else {
	    s += FormatOutput::logTabPrintf(0,
     "WARNING: Run %3d has insuffient information for detector (tile) scaling\n",
					  run.RunNumber());
	  }
	}
      }
    } // end detector
    return s;
  }
  //--------------------------------------------------------------
  int ScaleModel::scaleSpecIndex(const int& irun,
				 const std::vector<ScaleSpecification>& scaleSpecs,
				 const std::vector<Run>& runList) const
  // Returns index into scale specification list for run index irun
  //  returns -1 if not found
  {
    // Look for explicitly specified run matching irun, & find the general spec
    int isp0 = -1;
    for (size_t isp=0;isp<scaleSpecs.size();isp++) {
      int runNumber = scaleSpecs[isp].run;
      if (runNumber > 0) {
	int ir = FindRunIndex(runNumber, runList);
	if (ir == irun) {return isp;}
      } else {
	isp0 = isp;  // this is the general spec
      }
    }
    return isp0;
  }
  //--------------------------------------------------------------
  void ScaleModel::PrintLayout(phaser_io::Output& output)
  {
    output.logTab(0,LOGFILE,"\n>>>> Layout of scale factors: <<<<\n");
    for (size_t irun=0;irun<runnumbers.size();++irun) {
      output.logTabPrintf(0,LOGFILE,
			  "\nRun %5d\n", runnumbers[irun]);
      output.logTab(0,LOGFILE, primary_scales.at(irun).format());
      if (nbfactors == 0) {
	output.logTab(0,LOGFILE, "No relative B-factors");
      } else {
	output.logTab(0,LOGFILE, relative_bfactors.at(irun).format());
      }
      if (sec_scale_index_run.at(irun) >= 0) {
	output.logTab(0,LOGFILE,
		      secondary_scales[sec_scale_index_run.at(irun)].format());
      }
      if (detector_scale_index_run.at(irun) >= 0) {
	output.logTab(0,LOGFILE,
		      detector_scales[detector_scale_index_run.at(irun)].format());
      }
    }
    if (nties > 0) {
      output.logTab(0,LOGFILE,"\n");
      if (nties_rot > 0) {
	output.logTabPrintf(0,LOGFILE,
			    (std::string("Primary scales will be TIED together along ")+
			     "the rotation axis with a standard deviation of %6.3f, "+
			     "number of ties %5d\n").c_str(), sd_rotation, nties_rot);
      }
      if (nties_bfac > 0) {
	output.logTabPrintf(0,LOGFILE,
			    (std::string("B-factors will be TIED together along ")+
			     "the rotation axis with a standard deviation of %6.3f, "+
			     "number of ties %5d\n").c_str(), sd_bfactor, nties_bfac);
      }
      if (nties_zerob > 0) {
	output.logTabPrintf(0,LOGFILE,
			    (std::string("B-factors will be TIED to zero ")+
			     "with a standard deviation of %6.3f, "+
			     "number of ties %5d\n").c_str(), sd_zerob, nties_zerob);
      }
      if (nties_surf > 0) {
	output.logTabPrintf(0,LOGFILE,
			    (std::string("Secondary beam parameters will be TIED to zero, ")+
			     "ie restrained to a sphere, \n"+
			     "      with a standard deviation of %6.3f, "+
			     "number of ties %5d\n").c_str(), sd_surface, nties_surf);
      }
      if (nties_tiles > 0) {
	output.logTabPrintf(0,LOGFILE,
			    (std::string("Detector parameters r,w,A will be TIED across the tiles,")+
			     " with standard deviations %6.3f,%6.3f,%6.3f,\n"+
			     "  and tile centre positions (x0,y0) will be tied to the true centre with standard deviation %6.3f\n").c_str(), sd_tile[0], sd_tile[1], sd_tile[2], sd_tile[3]);
      }
    }
    output.logTabPrintf(0,LOGFILE,"\n");
  }
  //--------------------------------------------------------------
  std::string ScaleModel::PrintTwoWrappingLines(const std::vector<double>& v,
						const std::string& t1,
						const std::vector<int>& n,
						const std::string& t2,
						const int& fw,
						const int& fd)
  // Print wrapping lines:
  //   line 1, values v (double), label t1
  //   line 2, values n (int),    label t2
  //   fw  field width
  //   fd  number of decimal points for v
  {
    std::string ss;
    int tlen = Max(t1.size(), t2.size()); // maximum label length
    const int PAGEWIDTH = 100; // characters across
    int nperline = (PAGEWIDTH - tlen - 1)/fw; // number/line
    ASSERT (v.size() == n.size());

    int nitem = v.size();
    int i1 = 0; // 1st item in line
    while (i1 < nitem) {
      ss += "\n"; // end of line
      int i2 = Min(i1 + nperline, nitem) - 1; // last item in line
      // v
      ss += StringUtil::LeftString(t1+":", tlen+1);
      for (int i=i1;i<=i2;++i) {
	ss += StringUtil::ftos(float(v[i]), fw, fd);
      }
      ss += "\n"; // end of line
      // n
      ss += StringUtil::LeftString(t2+":", tlen+1);
      for (int i=i1;i<=i2;++i) {
	ss += clipper::String(n[i], fw);
      }
      ss += "\n"; // end of line
      i1 = i2+1;
    }
    return ss;
  }
  //--------------------------------------------------------------
  void ScaleModel::PrintScales(phaser_io::Output& output)
  // Print all scale parameters
  {
    const int nperline = 10;
    output.logTab(0,LOGFILE,"\nScale parameters:\n");
    for (size_t irun=0;irun<runnumbers.size();++irun) {
      output.logTabPrintf(0,LOGFILE,
			  "\nRun %5d\n", runnumbers[irun]);
      // Primary
      output.logTab(0,LOGFILE,"Primary scales and number of observations");
      std::vector<double> pscales = primary_scales[irun].Scales();
      for (size_t i=0;i<pscales.size();++i) {pscales[i]=1.0/pscales[i];}
      std::vector<int> nobsPar = primary_scales[irun].Nobservations();
      output.logTab(0,LOGFILE,
		    PrintTwoWrappingLines(pscales, "Scales", nobsPar, "Nobs",9,3));
      // B-factors
      std::vector<double> bfacs = relative_bfactors[irun].Bfactors();
      nobsPar = relative_bfactors[irun].Nobservations();
      if (bfacs.size() > 0) {
	output.logTab(0,LOGFILE,"\nRelative B-factors and number of observations");
	output.logTab(0,LOGFILE,
		      PrintTwoWrappingLines(bfacs, "B-factors", nobsPar, "Nobs",9,3));
	output.logTabPrintf(0,LOGFILE,"\n");
      }
    } // end run loop
    if (nsecscales > 0) {   // Secondary
      output.logTab(0,LOGFILE,"Secondary scales");
      for (int j=0;j<nsecscales;++j) {
	std::vector<double> secsclpar = secondary_scales[j].Coefficients();
	for (size_t i=0;i<secsclpar.size();++i) {
	  if (i%nperline == 0) {output.logTabPrintf(0,LOGFILE,"\n");}
	  output.logTabPrintf(0,LOGFILE," %9.3f", secsclpar[i]);
	}
      }
    } // End Secondary
    if (ndetscales > 0) { // Tiles
      output.logTabPrintf(0,LOGFILE,"\n\n");
      for (int i=0;i<ndetscales;++i) {
	output.logTab(0,LOGFILE,detector_scales[i].formatparameters());
      }
    }
    output.logTabPrintf(0,LOGFILE,"\n\n");
  }
  //--------------------------------------------------------------
  void ScaleModel::CountParameters()
  // Set all parameter counts
  // Order of parameters:
  //   1. all primary scale parameters (nprimaryscale)
  //   2. all B-factors (nbfactors)
  //   3. secondary parameters (nsecondaryscale)
  //   4. detector (ntilescale)
  //   5. ..
  {
    nparameters = 0;
    nprimaryscale = 0;
    nbfactors = 0;
    // Primary scales & B-factors, one per run
    ASSERT (int(primary_scales.size()) == nruns);
    ASSERT (int(relative_bfactors.size()) == nruns);
    idxrun_primary_scales.resize(nruns);
    idxrun_bfactors.resize(nruns);
    idxrun_secondary.resize(nsecscales);
    idxrun_detector.resize(ndetscales);

    // Primary
    for (int irun=0;irun<nruns;irun++) {
      idxrun_primary_scales[irun] = nparameters;   // index to 1st primary parameter
      nparameters += primary_scales[irun].Number();
      nprimaryscale += primary_scales[irun].Number();
    }
    // B-factors
    for (int irun=0;irun<nruns;irun++) {
      idxrun_bfactors[irun] = nparameters;   // index to 1st bfactor parameter
      nparameters += relative_bfactors[irun].Number();
      nbfactors += relative_bfactors[irun].Number();
    }
    // Secondary
    ASSERT (nsecscales <= nruns);
    nsecondaryscale = 0;
    if (nsecscales > 0) {
      for (int i=0;i<nsecscales;++i) {
	idxrun_secondary[i] = nparameters;   // index to 1st secondary parameter
	nsecondaryscale += secondary_scales[i].Number();
	nparameters += secondary_scales[i].Number();
      }
    }
    // Detector
    ASSERT (ndetscales <= nruns);
    ntilescale = 0;
    if (ndetscales > 0) {
      for (int i=0;i<ndetscales;++i) {
	idxrun_detector[i] = nparameters;   // index to 1st secondary parameter
	ntilescale += detector_scales[i].Number();
	nparameters += detector_scales[i].Number();
      }
    }
    // Other things ...
  }
  //--------------------------------------------------------------
  std::vector<double> ScaleModel::GetParameters() const
  // Get vector of parameters
  // Order of parameters:
  //   1. all primary scale parameters (nprimaryscale)
  //   2. all B-factors (nbfactors)
  //   3. secondary parameters (nsecondaryscale)
  //   4. detector (ntilescale)
  //   5. ..
  {
    std::vector<double> params;

    for (int irun=0;irun<nruns;irun++) {
      // Append group of scales
      std::vector<double> pscales = primary_scales[irun].Scales();
      params.insert(params.end(), pscales.begin(), pscales.end());
    }
    ASSERT (int(params.size()) == nprimaryscale);
    // B-factors
    for (int irun=0;irun<nruns;irun++) {
      // Append group of B-factors
      std::vector<double> bfacs = relative_bfactors[irun].Bfactors();
      params.insert(params.end(), bfacs.begin(), bfacs.end());
    }
    ASSERT (int(params.size()) == nprimaryscale+nbfactors);
    // Secondary
    for (int i=0;i<nsecscales;++i) {
      std::vector<double> secsclpar = secondary_scales[i].Coefficients();
      params.insert(params.end(), secsclpar.begin(), secsclpar.end());
      //      std::cout <<  secondary_scales[i].Number() << " " << params.size() << "\n";
    }
    ASSERT (int(params.size()) == nprimaryscale+nbfactors+nsecondaryscale);
    // Detector
    for (int i=0;i<ndetscales;++i) {
      std::vector<double> detsclpar = detector_scales[i].Parameters();
      params.insert(params.end(), detsclpar.begin(), detsclpar.end());
      //      std::cout <<  detector_scales[i].Number() << " " << params.size() << "\n";
    }
    ASSERT (int(params.size()) == nprimaryscale+nbfactors+nsecondaryscale+ntilescale);
    // >>>>
    // Just check numbers
    ASSERT (int(params.size()) == nparameters);
    return params;
  }      
  //--------------------------------------------------------------
  std::vector<ScaleModel::ScaleParameterType> ScaleModel::GetParameterType() const
  // get type for all parameters
  {
    std::vector<ScaleModel::ScaleParameterType> partype(nparameters, ScaleModel::NONE);
    // Order of parameters:
    //   1. all primary scale parameters (nprimaryscale)
    //   2. all B-factors (nbfactors)
    //   3. secondary parameters (nsecondaryscale)
    //   4. detector (ntilescale)
    //   5. ..
    int k=-1;
    for (int i=0;i<nprimaryscale;++i) {partype[++k] = ScaleModel::SCALE;}
    for (int i=0;i<nbfactors;++i) {partype[++k] = ScaleModel::BFACTOR;}
    for (int i=0;i<nsecondaryscale;++i) {partype[++k] = ScaleModel::SECONDARY;}
    for (int i=0;i<ntilescale;++i) {partype[++k] = ScaleModel::TILE;}
    ASSERT (++k == nparameters);
    return partype;
  }
  //--------------------------------------------------------------
  ScaleModel::ScaleParameterType ScaleModel::GetParameterType(const int& Ipar) const
  // get type for a parameter
  {
    if (Ipar < 0 || Ipar >= nparameters) {
      clipper::Message::message(Message_fatal
				("GetLowerBound: parameter number out of range"+
				 clipper::String(Ipar)));
    }
    // Order of parameters:
    //   1. all primary scale parameters (nprimaryscale)
    //   2. all B-factors (nbfactors)
    //   3. secondary parameters (nsecondaryscale)
    //   4. ..
    if (Ipar < nprimaryscale) {return ScaleModel::SCALE;}
    if (Ipar < nprimaryscale+nbfactors) {return ScaleModel::BFACTOR;}
    if (Ipar < nprimaryscale+nbfactors+nsecondaryscale) {return ScaleModel::SECONDARY;}
    if (Ipar < nprimaryscale+nbfactors+nsecondaryscale+ntilescale) {
      return ScaleModel::TILE;
    }
    return ScaleModel::NONE;
  }
  //--------------------------------------------------------------
  //! return detector scale number & scale index into detector parameter list
  std::pair<int,int> ScaleModel::DetectorParameterNumber(const int& Ipar) const
  {
    // Find detector scale
    int idetsc;
    bool found = false;
    for (idetsc=0;idetsc<idxrun_detector.size();++idetsc) {
      if (Ipar >= idxrun_detector[idetsc]) {
	found = true; break;
      }
    }
    if (!found) {idetsc = idxrun_detector.size()-1;}
    // Index into detector parameter list
    ASSERT (idetsc >= 0);
    return std::pair<int,int>(idetsc, Ipar - idxrun_detector[idetsc]);
  }
  //--------------------------------------------------------------
  // Set all parameters from vector
  void ScaleModel::SetParameters(const std::vector<float>& params,
				 const std::vector<int>& Nobs)
  {
    std::vector<double> pars(params.size());
    for (size_t i=0;i<params.size();++i) {pars[i] = params[i];}
    SetParameters(pars, Nobs);
  }
  //--------------------------------------------------------------
  // Set all parameters from vector
  void ScaleModel::SetParameters(const std::vector<double>& params,
				 const std::vector<int>& Nobs)
  {
    ASSERT (int(params.size()) == nparameters);
    std::vector<double>::const_iterator pos1 = params.begin();  // start of range
    std::vector<double>::const_iterator pos2;                   // end of range
    std::vector<int>::const_iterator posn1 = Nobs.begin();   // start of range
    std::vector<int>::const_iterator posn2;                  // end of range
    for (int irun=0;irun<nruns;irun++) {
      // Extract group of scales
      pos2 = pos1 + primary_scales[irun].Number();
      primary_scales[irun].StoreScales(std::vector<double>(pos1, pos2));
      pos1 = pos2;
      posn2 = posn1 + primary_scales[irun].Number();
      primary_scales[irun].StoreNobservations(std::vector<int>(posn1, posn2));
      posn1 = posn2;
    }
    for (int irun=0;irun<nruns;irun++) {
      // Extract group of B-factors
      pos2 = pos1 + relative_bfactors[irun].Number();
      relative_bfactors[irun].StoreBfactors(std::vector<double>(pos1, pos2));
      pos1 = pos2;
      posn2 = posn1 + relative_bfactors[irun].Number();
      relative_bfactors[irun].StoreNobservations(std::vector<int>(posn1, posn2));
      posn1 = posn2;
    }
    for (int i=0;i<nsecscales;++i) {
      // Secondary
      pos2 = pos1 + secondary_scales[i].Number();
      secondary_scales[i].StoreCoefficients(std::vector<double>(pos1, pos2));
      pos1 = pos2;
      posn2 = posn1 + secondary_scales[i].Number();
      secondary_scales[i].StoreNobservations(std::vector<int>(posn1, posn2));
      posn1 = posn2;
    }
    for (int i=0;i<ndetscales;++i) {
      // Detector
      pos2 = pos1 + detector_scales[i].Number();
      detector_scales[i].StoreParameters(std::vector<double>(pos1, pos2));
      pos1 = pos2;
      posn2 = posn1 + detector_scales[i].Number();
      detector_scales[i].StoreNobservations(std::vector<int>(posn1, posn2));
      posn1 = posn2;
    }
    NormaliseParameters();
  }
  //--------------------------------------------------------------
  void ScaleModel::SetInitialScales(const std::vector<double>& gscales,
				    const std::vector<int>& numobsrotrange)
  // Set initial primary scales, eg from InitialScales
  // also store count of number of observations
  {
    // For each run, number of initial scales = number of intervals, depends on type
    //  1) for batch scales, number = number of scales
    //  2) for smooth scales, Nintervals = Nscales - 2 
    //  3) ... unless Nscales = 1 or 2, in which case Nintervals = 1
    int niscl = 0;  // count number expected
    for (int ir=0;ir<nruns;++ir) {
      niscl += primary_scales[ir].Nintervals();
    }
    ASSERT (niscl == int(gscales.size()));
    ASSERT (niscl == int(numobsrotrange.size()));

    int jsr = 0; // index into gscales
    int j;       // index into scale array for each run
    for (int ir=0;ir<nruns;++ir) { // loop runs
      std::vector<double> gsclrun(primary_scales[ir].Number()); // scale parameters for this run
      std::vector<int>    nobsrun(primary_scales[ir].Number()); // number of observations for this run
      for (int i=0;i<primary_scales[ir].Nintervals();++i) {
	if (primary_scales[ir].IsBatchScale() || primary_scales[ir].Nintervals() == 1) {
	  // No leading scale
	  j = i;
	} else {
	  j = i+1;
	  // extra scale at start, duplicate of 1st interval
	  if (i == 0) {
	    gsclrun[i] = gscales[jsr];
	    nobsrun[i] = numobsrotrange[jsr];
	  }
	}
	gsclrun[j] = gscales[jsr];
	nobsrun[j] = numobsrotrange[jsr];
	jsr++;
      } // end loop scales intervals
      if (!primary_scales[ir].IsBatchScale() && primary_scales[ir].Number() >= 2) {
	// extra scale at end: j is one beyond intitial scales array
	gsclrun[j+1] = gsclrun[j];
	nobsrun[j+1] = numobsrotrange[j-1];
	j++;
      }
      ASSERT (j+1 == primary_scales[ir].Number());
      primary_scales[ir].StoreScales(gsclrun);
      primary_scales[ir].StoreNobservations(nobsrun);
    } // end loop runs
  }
  //--------------------------------------------------------------
  // Return true if model is refinable, ie not just one scale and one B-factor
  bool ScaleModel::IsRefinable() const
  {
    return (nprimaryscale > 1) || (nbfactors > 1) ||
      (nsecondaryscale > 0) || (ntilescale > 0);
  }
  //--------------------------------------------------------------
  double ScaleModel::ScaleObs(observation& obs, const Rtype& invresolsq) const
  // Scale observation, returns scale applied
  {
    // Run
    int irun = obs.run();
    double g = 1.0;
    double ps;    // Primary scale
    if (primary_scales[irun].IsBatchScale()) {
      ps = primary_scales[irun].Scale(obs.Batch());
    } else {
      ps = primary_scales[irun].Scale(obs.phi());
    }
    // B-factor scale
    double bs = 1.0;

    if (relative_bfactors[irun].Number() > 0) {
      if (relative_bfactors[irun].IsBatchBfactor()) {
	bs = relative_bfactors[irun].BfactorScale(obs.Batch(), invresolsq);
      } else {
	bs = relative_bfactors[irun].BfactorScale(obs.time(), invresolsq);
      }
    }
    // Secondary   FIXME (why?)
    double ss = 1.0;
    if (nsecscales > 0) {
      if (sec_scale_index_run[irun] >= 0) {
	double thetap, phip;
	obs.GetS2(thetap, phip);
	ss = secondary_scales[sec_scale_index_run[irun]].Scale(thetap, phip);
      }
    }
    // Detector
    double ds = 1.0;
    if (ndetscales > 0) {
      if (detector_scale_index_run[irun] >= 0) {
	ds = detector_scales[detector_scale_index_run[irun]].Scale(obs.XYdet());
      }
    }
    g = ps*bs*ss*ds;
    obs.SetGscale(g);
    return g;
  }
  //--------------------------------------------------------------
  double ScaleModel::ScaleObs(observation& obs, const Rtype& invresolsq,
			      std::vector<double>& dghldp) const
  // Scale observation, returns scale applied and
  // partial derivative vector d(ghl)/dp
  {
    // Run
    int irun = obs.run();
    double g = 1.0;

    dghldp.assign(nparameters, 0.0);


    double ps;    // Primary scale
    std::vector<double> dgdpm;  // derivatives for this run only
    if (primary_scales[irun].IsBatchScale()) {
      ps = primary_scales[irun].ScaleDeriv(obs.Batch(), dgdpm);
    } else {
      ps = primary_scales[irun].ScaleDeriv(obs.phi(), dgdpm);
    }
    
    // B-factor scale
    std::vector<double> dgdB;  // derivatives for this run only
    double bs = 1.0;
    if (relative_bfactors[irun].Number() > 0) {
      if (relative_bfactors[irun].IsBatchBfactor()) {
	bs = relative_bfactors[irun].BfactorScaleDeriv(obs.Batch(), invresolsq, dgdB);
      } else {
	bs = relative_bfactors[irun].BfactorScaleDeriv(obs.time(), invresolsq, dgdB);
      }
    }
    
    // Secondary 
    double ss = 1.0;
    std::vector<double> dgds;  // derivatives for this run only
    if (nsecscales > 0) {
      if (sec_scale_index_run[irun] >= 0) {
	// Return scale & derivatives for secondary beam directions
	// polar angles thetap, phip (calculated in ScaleModel constructor)
	double thetap, phip;
	obs.GetS2(thetap, phip);
	int k = sec_scale_index_run[irun];
	ss = secondary_scales[k].ScaleDeriv(thetap, phip, dgds);
      }
    }

    // Detector
    double ds = 1.0;
    std::vector<double> dgdd;  // derivatives for this run only
    if (ndetscales > 0) {
      if (detector_scale_index_run[irun] >= 0) {
	ds = detector_scales[detector_scale_index_run[irun]].
	  ScaleDeriv(obs.XYdet(), dgdd);
      }
    }

    // dghl/dp = dg(primary)/dp * bs * ss * ds
    for (size_t i=0;i<dgdpm.size();++i) {
      dgdpm[i] *= bs * ss *ds;
    }
    std::copy(dgdpm.begin(), dgdpm.end(),
	      dghldp.begin() + idxrun_primary_scales[irun]);

    // dghl/dp = dg(B)/dp * ps * ss * ds
    for (size_t i=0;i<dgdB.size();++i) {
      dgdB[i] *= ps * ss * ds;
    }
    std::copy(dgdB.begin(), dgdB.end(),
	      dghldp.begin()+idxrun_bfactors[irun]);

    if (nsecscales > 0) {
      // dghl/dp = dg(sec)/dp * ps * bs
      for (size_t i=0;i<dgds.size();++i) {
	dgds[i] *= ps * bs * ds;
      }
      int k = sec_scale_index_run[irun];
      std::copy(dgds.begin(), dgds.end(), dghldp.begin()+idxrun_secondary[k]);
    }

    if (ndetscales > 0) {
      // dghl/dp = dg(det)/dp * ps * bs *ss
      for (size_t i=0;i<dgdd.size();++i) {
	dgdd[i] *= ps * bs * ss;
      }
      int k = detector_scale_index_run[irun];
      std::copy(dgdd.begin(), dgdd.end(), dghldp.begin()+idxrun_detector[k]);
    }

    g = ps*bs*ss*ds;
    obs.SetGscale(g);

    return g;
  }
  //--------------------------------------------------------------
  void ScaleModel::NormaliseParameters()
  // Normalise scales & B-factors
  // Order of parameters:
  //   1. all primary scale parameters (nprimaryscale)
  //   2. all B-factors (nbfactors)
  //   3. secondary parameters (nsecondaryscale)
  //   4. detector (ntilescale)
  //   5. ..
  {
    // Scales for normalisation run
    ASSERT (scalenormbatch >= 0);
    // normalisation factor for scales
    double scnorm = primary_scales[scalenormrun].Scales()[scalenormbatch];
    if (scnorm <= 0.0) {
      std::vector<double> pscales = primary_scales[scalenormrun].Scales();
      clipper::Message::message(Message_fatal
				("Normalisation scale < 0"));
    }
    for (int irun=0;irun<nruns;irun++) {
      // Scales for this run
      std::vector<double> pscales = primary_scales[irun].Scales();
      for (size_t i=0;i<pscales.size();++i) {
	pscales[i] /= scnorm;
      }
      primary_scales[irun].StoreScales(pscales);
    }
    // B-factors
    double bfnorm = -1000000.;  
    if (bfacnormbatch >= 0) {
      // Normalisation batch specified
      bfnorm = relative_bfactors[bfacnormrun].Bfactors()[bfacnormbatch];
    } else {
      // Find largest Bfactor
      for (int irun=0;irun<nruns;irun++) {
	//  B-factors for this run
	std::vector<double> bfacs = relative_bfactors[irun].Bfactors();
	// Normalisation on "best" batch: if smoothed Bfactors && > 2, omit first & last
	int i1 = 0;
	int i2 = bfacs.size();
	if (!relative_bfactors[irun].IsBatchBfactor() &&
	    i2 > 2) {
	  i1 = 1;
	  i2--;
	}
	for (int i=i1;i<i2;++i) {
	  bfnorm = Max(bfnorm, bfacs[i]);
	}
      }
      if (bfnorm > -999999.) {
	for (int irun=0;irun<nruns;irun++) {
	  //  B-factors for this run
	  std::vector<double> bfacs = relative_bfactors[irun].Bfactors();
	  for (size_t i=0;i<bfacs.size();++i) {
	    bfacs[i] -= bfnorm;
	  }
	  relative_bfactors[irun].StoreBfactors(bfacs);
	}
      }
    }
  }      
  //--------------------------------------------------------------
  double ScaleModel::TieValues
  (const bool& DoGradient, const bool& DoHessian,
   const std::vector<double>& params,
   std::vector<double>& dRdpi,
   std::vector<TieHessian>& Htie)
  // Return restraint target, and optionally gradient & Hessian contributions
  //
  // On entry:
  //  DoGradient  true to calculate gradient dRdpi
  //  DoHessian   true to calculate Hessian terms Htie
  //
  // On exit:
  //  dRdpi(nparameters)  gradient contribution for each parameter
  //  Htie                list of indexed Hessian contributions
  //  
  {
    if (DoGradient) {
      dRdpi.assign(nparameters, 0.0); // clear derivatives
      if (DoHessian) {
	Htie.clear();                   // and Hessian
      }}

    double R = 0.0;
    for (int itie=0;itie<nties;++itie) { // loop ties
      R += ties[itie].R(params);
      if (DoGradient) {
	std::vector<TieGradient> grad = ties[itie].Gradient(params);
	for (size_t i=0;i<grad.size();++i) {
	  dRdpi[grad[i].index] = grad[i].Grad;
	}
	if (DoHessian) {
	  std::vector<TieHessian> hessian = ties[itie].Hessian(params);
	  Htie.insert(Htie.end(), hessian.begin(), hessian.end());
	}
      }
    }
    return R;
  }
  //--------------------------------------------------------------
  bool ScaleModel::GetLowerBound(const int& Ipar, double& Lower) const
  // get lower bound for parameter, depending on type:
  // return false if unbounded
  {
    switch (GetParameterType(Ipar)) {
    case ScaleModel::SCALE:
      // A primary scale factor, cannot be negative
      Lower = 0.01;
      return true;
    case ScaleModel::BFACTOR:
      // A relative B-factor, leave unbounded
      return false;
    case ScaleModel::SECONDARY:
      // A secondary beam parameter, leave unbounded
      return false;
    case ScaleModel::TILE:
      {    
	// A detector beam parameter, but which one?
	std::pair<int,int> idxpar = DetectorParameterNumber(Ipar);
	return detector_scales[idxpar.first].LowerBound(idxpar.second, Lower);
      }
    default:
      return false;
    }
    return false;
  }
  //--------------------------------------------------------------
  bool  ScaleModel::GetUpperBound(const int& Ipar, double& Upper) const
  // get upper bound for parameter, depending on type:
  // return false if unbounded
  {
    switch (GetParameterType(Ipar)) {
    case ScaleModel::SCALE:
      // A primary scale factor, no upper bound
      return false;
    case ScaleModel::BFACTOR:
      // A relative B-factor, leave unbounded
      return false;
    case ScaleModel::SECONDARY:
      // A secondary beam parameter, leave unbounded
      return false;
    case ScaleModel::TILE:
      {
	// A detector beam parameter, but which one?
	std::pair<int,int> idxpar = DetectorParameterNumber(Ipar);
	return detector_scales[idxpar.first].UpperBound(idxpar.second, Upper);
      }
    default:
      return false;
    }
    return false;
  }
  //--------------------------------------------------------------
  double ScaleModel::GetLargeShift(const int& Ipar) const
  // get "large shift" value for parameter, depending on type
  {
    switch (GetParameterType(Ipar)) {
    case ScaleModel::SCALE:
      // A primary scale factor
      return 1.0;
    case ScaleModel::BFACTOR:
      // A relative B-factor
      return 2.0;
    case ScaleModel::SECONDARY:
      // A secondary beam parameter
      return 0.1;
    case ScaleModel::TILE:
      {    
	// A detector beam parameter, but which one?
	std::pair<int,int> idxpar = DetectorParameterNumber(Ipar);
	return detector_scales[idxpar.first].LargeShift(idxpar.second);
      }
    default:
      return 0.0;
    }
    return false;
  }
  //--------------------------------------------------------------
  void ScaleModel::Save(const std::string& dumpfilename,
			const std::vector<Run>& runlist) const
  // Dump scale model to file
  {
    FILE* dumpfile = OpenFile(dumpfilename, true);
    if (dumpfile == NULL) {
      clipper::Message::message(Message_fatal
				("Failed to open ScaleModelDumpFile "+dumpfilename));
    }
    fprintf(dumpfile, "%s", FormatSave(runlist).c_str());
    fclose(dumpfile);
  }
  //--------------------------------------------------------------
  std::string ScaleModel::FormatSave(const std::vector<Run>& runlist) const
  // Format scalemodel for dump/restore
  //  runlist  runs to be saved
  // NB ties are not saved
  {
    ASSERT (nruns == int(runlist.size()));
    const std::string SCALEMODELVERSION = "V1.1";
    std::string ds = "ScaleModel "+SCALEMODELVERSION+" {\n";
    ds += "Nruns "+ clipper::String(nruns)+"\n";;

    // Save essential run stuff to be matched on restore
    for (int i=0;i<nruns;++i) {
      ds += "RunNumber "+clipper::String(runnumbers[i])+"\n";
      ds += runlist[i].FormatSave();
    }
    // Primary beam things
    // number of runs == number of primary models
    ASSERT (int(primary_scales.size()) == nruns);
    for (int j=0;j<nruns;++j) { // primary scales
      ds += primary_scales[j].FormatSave();
    }
    // Relative B-factor
    int nbfacs = relative_bfactors.size();
    ds += "Nbfactors "+ clipper::String(nbfacs)+"\n";
    if (nbfacs > 0) { // B-factors
      ASSERT (nruns == nbfacs);
      for (int j=0;j<nbfacs;++j) {
	ds += relative_bfactors[j].FormatSave();
      }
    }
    // Secondary beam
    ds += "Nsecscales "+ clipper::String(nsecscales)+"\n";
    if (nsecscales > 0) {
      for (int i=0;i<nsecscales;++i) {
	ds += secondary_scales[i].FormatSave();
      }
      ASSERT (int(sec_scale_index_run.size()) == nruns);
      ds += "Sec_scale_index_run\n"+
	StringUtil::FormatSaveVector(sec_scale_index_run);
    }
    //t  // save tile parameters

    // sds
    ds += "sd_rotation "+ clipper::String(sd_rotation)+"\n";
    ds += "sd_bfactor "+ clipper::String(sd_bfactor)+"\n";
    ds += "sd_zerob "+ clipper::String(sd_zerob)+"\n";
    ds += "sd_surface "+ clipper::String(sd_surface)+"\n";
    ds += "sd_tile_number "+clipper::String(int(sd_tile.size()))+"\n";
    for (size_t i=0;i<sd_tile.size();++i) {
      ds += "sd_tile "+clipper::String(sd_tile[i])+"\n";
    }

    // Normalisation
    ds += "Scalenormrun "+ clipper::String(scalenormrun)+"\n";
    ds += "Scalenormbatch "+ clipper::String(scalenormbatch)+"\n";
    ds += "Bfacnormrun "+ clipper::String(bfacnormrun)+"\n";
    ds += "Bfacnormbatch "+ clipper::String(bfacnormbatch)+"\n";

    return ds+"}\n";
  }
  //--------------------------------------------------------------
  void ScaleModel::Restore(const std::string& restorefilename,
			   const std::vector<Run>& runlist)
  // Restore from dump file, for those runs which are defined now
  // Note that ties are not restored
  {
    std::ifstream scalesin(restorefilename.c_str());
    Fileread FR(scalesin, restorefilename, "RESTORE");

    FR.ReadTag("ScaleModel"); // Note ReadTag fails if tag is wrong
    if (FR.GetTag() != "V1.1") {  // version check
      clipper::Message::message(Message_fatal
				("RESTORE incompatible version in "+restorefilename));
    }
    FR.Skip(); // skip "{"
    FR.ReadTag("Nruns");
    int svnruns = FR.Int();  // number of runs in save file
    nruns = runlist.size();  // number of runs in runlist = number to be used
    runnumbers.resize(nruns);
    // list of runs from save file with corresponding run serials in runlist
    std::vector<int> runsfromsavefile(svnruns, -1);
    int nrfound = 0;
    

    // Loop runs in save file
    for (int ir=0;ir<svnruns;++ir) {
      FR.ReadTag("RunNumber");
      int runnum = FR.Int();
      runnum = runnum;
      FR.ReadTag("Run");
      if (FR.GetTag() != "V1") {  // version check
	clipper::Message::message(Message_fatal
				  ("RESTORE incompatible run version in "+restorefilename));
      }
      FR.Skip(); // skip "{"
      FR.ReadTag("Batch_number_list");
      int nbat = FR.Int();
      std::vector<int> batchnumbers = FR.IntVec(nbat);
      if (!FR.CheckEnd()) {
	clipper::Message::message(Message_warn
				  ("ScaleModel::Restore unexpected tag "+FR.Tag()));
      }
      // Does this run match any in runlist?
      // return index in runlist, -1 if not found
      int irun = RunNotFound(runlist, batchnumbers);
      if (irun >= 0) {
	// Build list of runs from save file which should be kept
	runsfromsavefile[ir] = irun;
	nrfound++;
	runnumbers[irun] = runlist[irun].RunNumber();
      }
    } // end loop runs in save file
      // runsfromsavefile now contains the index in runlist for each run in save file
    if (nrfound < int(runlist.size())) {
      clipper::Message::message(Message_fatal
				("RESTORE not all runs found in save file"));
    }
    // number of primary scales == number of runs
    // number to use = nruns
    int jpr = 0; // index to accepted runs/primary scales
    for (int ipr = 0;ipr<svnruns;++ipr) { // loop primary scales in file
      // Do we want this one?
      if (runsfromsavefile[ipr] < 0) {
	// No, skip it
	FR.SkipSection(0); // skip section
      } else {
	// Yes, read it
	primary_scales[jpr++].Restore(FR);
      }
    } // end loop primary scales

    FR.ReadTag("Nbfactors");
    int nbf = FR.Int();
    if (nbf != svnruns) {
      clipper::Message::message(Message_fatal
				("RESTORE number of B-factors != number of runs"));
    }
    jpr = 0;
    for (int ipr = 0;ipr<svnruns;++ipr) { // loop B-factors in file
      // Do we want this one?
      if (runsfromsavefile[ipr] < 0) {
	// No, skip it
	FR.SkipSection(0); // skip section
      } else {
	// Yes, read it
	relative_bfactors[jpr++].Restore(FR);
      }
    } // end loop B-factors

    FR.ReadTag("Nsecscales");
    int nssc = FR.Int();
    for (int i = 0;i<nssc;++i) { // loop secondary scales in file
      secondary_scales[i].Restore(FR);
    } // end loop secondary scales
    FR.ReadTag("Sec_scale_index_run");
    // Index for each run
    sec_scale_index_run.clear();
    for (int ipr = 0;ipr<svnruns;++ipr) { // loop sec scale indices
      int js = FR.Int();
      // Do we want this one?
      if (runsfromsavefile[ipr] >= 0) {
	sec_scale_index_run.push_back(js);
      }
    }
    if (int(sec_scale_index_run.size()) != nruns) {
      clipper::Message::message(Message_fatal
				("RESTORE number of Sec_scale_index_runs != number of runs"));
    }
    //t  tile

    // sds
    FR.ReadTag("sd_rotation"); sd_rotation = FR.Double();
    FR.ReadTag("sd_bfactor"); sd_bfactor = FR.Double();
    FR.ReadTag("sd_zerob"); sd_zerob = FR.Double();
    FR.ReadTag("sd_surface"); sd_surface = FR.Double();
    int sdtn;
    FR.ReadTag("sd_tile_number"); sdtn = FR.Int();
    sd_tile.resize(sdtn);
    for (size_t i=0;i<sdtn;++i) {
      FR.ReadTag("sd_tile"); sd_tile[i] = FR.Double();
    }
    // Normalisation
    FR.ReadTag("Scalenormrun"); scalenormrun = FR.Int();
    FR.ReadTag("Scalenormbatch"); scalenormbatch = FR.Int();
    FR.ReadTag("Bfacnormrun"); bfacnormrun = FR.Int();
    FR.ReadTag("Bfacnormbatch"); bfacnormbatch = FR.Int();

    CountParameters(); // set parameter counts etc
  }
  //--------------------------------------------------------------
  // Does this run match any in runlist?
  // return index in runlist, -1 if not found
  int ScaleModel::RunNotFound(const std::vector<Run>& runlist,
			      const std::vector<int>& batchnumbers) const
  {
    int nbn = batchnumbers.size();
    for (size_t ir=0;ir<runlist.size();++ir) { // loop runs
      std::vector<int> bl = runlist[ir].BatchList();
      if (int(bl.size()) == nbn) { // same size
	if (bl == batchnumbers) { // all same
	  return ir;
	}}
    }
    return -1;
  }
  //--------------------------------------------------------------
  void ScaleModel::WriteImage(const std::string fname) const
  //! Write image[s] for each detector scale
  {
    std::string imagefilename = fname;
    if (getenv(imagefilename.c_str()) != NULL) { // it's an environment variable
      imagefilename = std::string(getenv(imagefilename.c_str()));
    }

    for (int idsc=0;idsc<ndetscales;++idsc) {
      std::string basename = FileNameNoExtension(imagefilename);
      std::string ext = FileNameExtension(imagefilename);
      if (ext == "") {ext = "img";}
      std::string name = imagefilename;
      name = basename+"_"+StringUtil::Strip(StringUtil::itos(idsc+1,4))+"."+ext;
      detector_scales[idsc].WriteImage(name);
    }
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  ValidScaleModel::ValidScaleModel(hkl_unmerge_list& hkl_list)
  {
    init(hkl_list);
  }
  //--------------------------------------------------------------
  void ValidScaleModel::init(hkl_unmerge_list& hkl_list)
  // Interrogate reflection list batches for each run to determine if there is sufficient
  // information to support each type of scale type
  {
    nruns = hkl_list.num_runs();

    // Set all flags to start with, depending on which columns are present
    // smooth primary scale (& B-factor)
    data_flags colflags = hkl_list.DataFlags();
    validprimary.assign(nruns,colflags.is_Rot);
    // batch number
    validbatch.assign(nruns,colflags.is_batch);
    // secondary beam, ie full geometry
    validsecondary.assign(nruns,true);
    // tile information
    validtile.assign(nruns,(colflags.is_Xdet && colflags.is_Ydet));

    // For valid batch numbers, we need more than one batch in the run
    std::vector<Run> runs = hkl_list.RunList();
    for (int irun=0;irun<nruns;++irun) {
      if (runs[irun].Nbatches() <= 1) {
	validbatch.at(irun) = false;
      }}

    std::vector<Batch> batches = hkl_list.Batches();  // all batch data
    int nbatches = batches.size();

    for (int ib=0;ib<nbatches;++ib) { // loop batches
      int irun = batches[ib].RunIndex();  // run serial number
      if (irun >= 0) {  // only for batches assigned to a run
	// For smooth primary beam corrections (scale & B-factor) we need valid Phi values
	if (!batches[ib].ValidPhi() || std::abs(batches[ib].PhiRange()) < 0.0001) {
	  validprimary.at(irun) = false;
	  // Primary data missing implies secondary as well
	  validsecondary.at(irun) = false;
	}
	// For valid secondary beam calculation, we need full geometry
	if (!batches[ib].ValidOrientation()) {
	  validsecondary.at(irun) = false;
	}
	// For valid tile information, we need detector [pixel] coordinates
	DetectorType dettype(batches[ib]);
	if (dettype.XdetRange().AbsRange() < 0.001) {
	  validtile.at(irun) = false;
	}
	if (dettype.YdetRange().AbsRange() < 0.001) {
	  validtile.at(irun) = false;
	}
      }   // if valid run
    } // end loop batches
  }
  //--------------------------------------------------------------
}
