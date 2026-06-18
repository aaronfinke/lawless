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
#include "restore.hh"
#include "report_errors.hh"
#include "median.hh"

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
using phaser_io::LXML;

namespace scala {
  //--------------------------------------------------------------
  ScaleModel::ScaleModel(const phaser_io::InputAll& input,
                         hkl_unmerge_list& hkl_list,
			 all_controls& controls,
                         phaser_io::Output& output)
  // Construct from input commands and reflection list
  {
    nwavscale = 0;
    idxwavscale = 0;
    wavelength_only_mode_ = false;
    gpr_requested = false;
    init (input, hkl_list, controls, output);
  }
  //--------------------------------------------------------------
  void ScaleModel::init(const phaser_io::InputAll& input,
                        hkl_unmerge_list& hkl_list,
			all_controls& controls,
                        phaser_io::Output& output)
  // initialise from input commands and reflection list
  // Scale specification(s) from input
  {
    status = 0;
    // Refine reference, allow on scale parameter
    refinereference = controls.refinecontrol.Reference();

    // Setup scale model
    pole = 0;
    std::vector<scala::ScaleSpecification> scaleSpecs =
      input.getScaleSpecifications();
    LinkSpecs linkspecs = input.getLINKs();
    bool notile = input.noTile();  // true if explicit NOTILE command
    if (!notile) {
      // set automatic TILE settings if appropriate
      autoTiles(scaleSpecs, hkl_list, output);
    }
    nwavscale = 0;
    idxwavscale = 0;
    wavelength_only_mode_ = false;
    gpr_requested = input.HasGPR();
    gpr_lambda_ref = input.getLambdaRef();
    if (gpr_requested) gpr_control = input.getGPRControl();
    setup(scaleSpecs, linkspecs, hkl_list, output);
    bool haveWavelength = input.IsLaue() &&
      (!input.getWavelengthRanges().empty() || gpr_requested);
    if (status < 0) {
      if (!haveWavelength) {return;}  // truly insufficient information
      // Wavelength normalization requested but primary scaling has
      // insufficient information (e.g. LAMBDAONLY without SCALES CONSTANT).
      // Fall back to constant primary scaling so only the wavelength
      // Chebyshev coefficients are refined.
      output.logTab(0, phaser_io::LOGFILE,
        "\nPrimary scaling has insufficient information; "
        "using constant primary scale with wavelength normalization only\n");
      SetConstant(hkl_list, output);
    }
    status = +1;

    // Wavelength (Chebyshev) normalization if LAUE keyword given
    if (haveWavelength) {
      wavelength_scale = WavelengthChebyshevScale(input.getWavelengthRanges(),
                                                   input.getLambdaRef());
      CountParameters();   // recount to include wavelength parameters
      VC.resize(nparameters, nparameters, 0.0);
      varpar.assign(nparameters, 0.0);
    }

    std::vector<Run> runlist = hkl_list.RunList();
    // run number for each lattice number-1 (lattices are numbered from 1)
    // only set to >= 0 for active main lattices corresponding to a run
    idxrunlattice.assign(hkl_list.NumberofLattices(), -1);
    for (int irun=0;irun<nruns;irun++) {
      int latnum = runlist[irun].LatticeNumber();
      if (latnum > 0) {
        idxrunlattice.at(latnum-1) = irun;
      }
    }

    // Scale normalisation  ..............................
    // set up scale normalisation flags from input if necessary
    scalenormrun = -1;
    bfacnormrun = -1;
    bfacnormbatch = input.getBfacNormBatchNumber();
    if (bfacnormbatch == -2) {
      // Specified Run number only if FIRST batch
      if (input.getBfacNormRunNumber() >= 0) {
	bfacnormrun = input.getBfacNormRunNumber();
	bool ok = false;
	for (int irun=0;irun<nruns;irun++) {
	  if (runlist[irun].RunNumber() == bfacnormrun) {
	    ok = true;
	    bfacnormrun = irun;  // run serial number
	  }
	}
	if (!ok) {
	  std::string runmessage = "BFACTOR RUN number "+
	    StringUtil::itos(bfacnormrun)+" is not known";
	  ReportErrors::printFatalError(runmessage);
	}
      }
    }
    normalisebfac = true;
    nfreedom = -1; // no use of variance data
    parametersdusage = input.getUSESDPARAMETER(); // set use flag from input
    if (parametersdusage != scala::ScaleSpecification::NONE) {
      nfreedom = 0; // will be set later
    }

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
      if (relative_bfactors[irun].Number() == 0) {
        // no Bfactor for this run, so no normalisation needed for any run
        normalisebfac = false;
      } else {
	// Bfactors
        // if bfacnormbatch == -2, then use first one in specified run
        // bfacnormrun may be -1 to use the first one
        if (bfacnormbatch > 0) {  // specified batch
	  // NB this does npt work and is disabled
	  //  need to work out which rotation range contains the specified batch
          if (bfacnormrun < 0) {  // should always be -1 if bfacnormbatch > 0
            // Normalisation batch specified, is it in this run?
            if (runlist[irun].IsInList(bfacnormbatch)) {
              // Yes, locate it (index j)
              std::vector<int> batchnums = runlist[irun].BatchList(true);
              for (size_t i=0;i<batchnums.size();++i) {
                if (batchnums[i] == bfacnormbatch) {
                  bfacnormbatchserial = i;   // batch serial
                  bfacnormrun = irun;  // run serial
                  break;
                }
              }
            }
          }
        }
        // else leave bfacnormbatch = -1, normalise on "best" batch later
      }
    }    // normalisation .............................

    // Ties
    SetupTies(input);
    // Always calculate all secondary beams & diffraction vectors
    bool secbeamsOK = hkl_list.CalcSecondaryBeams(pole);
    negativeSecScale = 0;
    negativeSecScaleLast = 0;
    negativeSecScaleOccurred = 0;

    if (nsecscales > 0 && !secbeamsOK) {
      // if we have secondary scales we must have secondary beams calculated
      //  but we shouldn't get here anyway!
      Message::message(Message_fatal
         ("ScaleModel: can't use Secondary Scale unless we have valid orientation information"));
    }
  }
  //--------------------------------------------------------------
  // return > 0 if a negative secondary scale occurred, and clear the flag
  int ScaleModel::NegativeSecScaleOccurred() {
    int count = negativeSecScaleOccurred;
    negativeSecScaleOccurred = 0;
    return count;
  }
  //--------------------------------------------------------------
  void ScaleModel::SetConstant(hkl_unmerge_list& hkl_list, phaser_io::Output& output)

  // set SCALE CONSTANT for all runs
  {
    scala::ScaleSpecification spec;
    spec.SetConstant();
    std::vector<scala::ScaleSpecification> scaleSpecs(1, spec);
    setup(scaleSpecs, LinkSpecs(), hkl_list, output);
  }
  //--------------------------------------------------------------
  void ScaleModel::setup(const std::vector<scala::ScaleSpecification>& scaleSpecs,
                         const LinkSpecs& linkspecs,
                         hkl_unmerge_list& hkl_list,
                         phaser_io::Output& output)
  // Setup from scale specifications and reflection list
  // Sets pole, for ABSORPTION, = 1,2,3 for h,k,l, = -1 unspecified, = 0 SECONDARY
  {
    scalenormbatch= -1; // batch number for scale  normalisation, -1 for 1st
    bfacnormbatch = -1;   // batch number for B-factor  normalisation, -1 for best
    bfacnormbatchserial = 0; // for safety
    nruns = hkl_list.num_runs();
    std::vector<Run> runlist = hkl_list.RunList();
    std::vector<Dataset> datasets = hkl_list.AllDatasets();

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
    std::vector<std::string> runXname(nruns,"");
    std::vector<double> runwavelengths(nruns,0.0);
    pole = 0; // for ABSORPTION, = 1,2,3 for h,k,l, = -1 unspecified, = 0 SECONDARY
    // Tile stuff
    int ktlidx = -1;
    int k0 = -1;  // tile
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
      //  Assign secondary scales to runs here, but these assignments will be
      //  replaced later if there are any LINK SURFACE commands
      int datasetidx;
      if (scaleSpecs[isp].sec_abs != scala::SecondaryScale::NONE &&
          validscalemodel.ValidSecondary(irun)) {
        if (j0 < 0) {
          // first run with a secondary correction
          j0 = irun;  //
          sec_scale_index_run[j0] = ++kscidx;  // index for 1st run = 0
          runXname[j0] = runlist[j0].PXDname().xname();
          datasetidx = runlist[j0].DatasetIndex();
          runwavelengths[j0] = datasets[datasetidx].wavelength(runXname[j0]);
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
          runXname[irun] = runlist[irun].PXDname().xname();
          datasetidx = runlist[irun].DatasetIndex();
          runwavelengths[irun] = datasets[datasetidx].wavelength(runXname[irun]);
          // search previous runs for same dataset
          bool found = false;
          int jrun;
          for (jrun=0;jrun<irun;jrun++) {
            if (runXname[jrun] == runXname[irun]) {
              // irun is same dataset as jrun
              found = true;
              break;
            }
          }
          if (found) {
            // irun should have same scale as jrun
            // but check for different wavelengths
            if (testwavelengths(runwavelengths[irun], runwavelengths[jrun])) {
              // similar wavelengths
              sec_scale_index_run[irun] = sec_scale_index_run[jrun];
            } else {
              // different wavelengths, irun is new dataset
              sec_scale_index_run[irun] = ++kscidx;
            }
          } else {
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
    //    (re)define sec_scale_index_run & nsecscales
    processLinks(linkspecs, runlist, output);

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
          std::string s = SetupScale(irun, scaleSpecs[isp], runlist[irun],
                                     validscalemodel, output);
          if (s != "") {output.logTab(0,LOGFILE,s);}
          if (status < 0) {return;}
          runSetup[irun] = true;
        }
      }
    }
    // Then fill in the rest from the default specs
    for (int irun=0;irun<nruns;irun++) {
      if (!runSetup[irun]) {
        // Specification for this run
        std::string s = SetupScale(irun, scaleSpecs[0], runlist[irun],
                                   validscalemodel, output);
        if (s != "") {output.logTab(0,LOGFILE,s);}
        if (status < 0) {return;}
        runSetup[irun] = true;
      }
    }
    // Count parameters, set nparameters & idxrun indices
    CountParameters();
    // initialise variances
    VC.resize(nparameters,nparameters,0.0);
    varpar.assign(nparameters, 0.0);
  }
  //--------------------------------------------------------------
  bool ScaleModel::testwavelengths(const double& wavelength1,
                                   const double& wavelength2) const
  // return true if wavelengths are similar
  {
    const double TEST = 0.1;
    bool tested = Close<double, double>(wavelength1, wavelength2, TEST);
    return tested;
  }
  //--------------------------------------------------------------
  void ScaleModel::processLinks(const LinkSpecs& linkspecs,
                                const std::vector<Run>& runlist,
                                phaser_io::Output& output)
  // (re)define sec_scale_index_run & nsecscales
  {
    // sec_scale_index_run will contain the scale index for each run
    if (!linkspecs.isSet()) {return;}  // no links, nothing to do
    if (nsecscales <= 0) {return;}

    ASSERT (sec_scale_index_run.size() == runlist.size());
    int numruns = runlist.size();

    if (linkspecs.linkAll()) {
      // LINK ALL
      for (int irun=0;irun<numruns;++irun) {
        sec_scale_index_run[irun] = 0;
      }
      nsecscales = 1;
    } else if (linkspecs.unlinkAll()) {
      // UNLINK ALL
      for (int irun=0;irun<numruns;++irun) {
        sec_scale_index_run[irun] = irun;
      }
      nsecscales = numruns;
    }
    // LINK r2 TO r1
    if (linkspecs.Nlinks() > 0) {
      std::vector<std::pair<int,int> > links = linkspecs.Links();
      for (size_t k=0; k<links.size(); k++) {
        std::pair<int,int> runidx = checkLink(links[k], runlist);
        if (runidx.first >= 0) {
          // scale index to remove
          int sclidx = sec_scale_index_run[runidx.first];
          // list runs with this index
          std::vector<int> runswithscaleindex = runsWithScaleIndex(sclidx);
          for (size_t k=0; k<runswithscaleindex.size(); k++) {
            sec_scale_index_run[runswithscaleindex[k]] =
              sec_scale_index_run[runidx.second];
          }
          // we have now eliminated scale index sclidx, decrement all
          // indices greater than this
          for (size_t irun=0; irun<sec_scale_index_run.size(); irun++) {
            if (sec_scale_index_run[irun] > sclidx) {
              sec_scale_index_run[irun]--;
            }
          }
        }
      }
      // reset and check nsecscales
      nsecscales = numberSecondaryScale();
    } // end explicit LINKs
    // UNLINK r2 TO r1
    if (linkspecs.Nunlinks() > 0) {
      std::vector<std::pair<int,int> > unlinks = linkspecs.Unlinks();
      for (size_t k=0; k<unlinks.size(); k++) {
        std::pair<int,int> runidx = checkLink(unlinks[k], runlist);
        if (runidx.first >= 0) {
          // check that these runs are linked, ignre if not
          if (sec_scale_index_run[runidx.first] ==
            sec_scale_index_run[runidx.second]) {
            // reassign runidx.first to next available slot
            sec_scale_index_run[runidx.first] = nsecscales;
            nsecscales++; // increment number
          }
        }
      }
      // reset and check nsecscales
      nsecscales = numberSecondaryScale();
    }  // end explicit UNLINKs
  }
  //--------------------------------------------------------------
  int ScaleModel::numberSecondaryScale() const
  // get number of secondary scales from sec_scale_index_run, check that all
  // are present
  {
    // can't be more scales than number of runs
    std::vector<int> sclidx(sec_scale_index_run.size(), -1);
    int maxsclidx = -1;
    for (size_t irun=0; irun<sec_scale_index_run.size(); irun++) {
      int idxscl = sec_scale_index_run[irun];
      maxsclidx = std::max(idxscl, maxsclidx);
      sclidx[idxscl]++;
    }
    for (int i=0;i<=maxsclidx;++i) {
      if (sclidx[i] < 0) {
        // this scale i has not been allocated to a run
        Message::message(Message_fatal
                         ("ScaleModel: secondary scale number "+
                          clipper::String(i)+" is not allocated to a run"));
      }
    }
    return maxsclidx+1;
  }
  //--------------------------------------------------------------
  std::vector<int> ScaleModel::runsWithScaleIndex(const int& sclidx) const
  // list of run indices which share scale index sclidx
  {
    std::vector<int> runswithscaleindex;
    for (size_t irun=0; irun<sec_scale_index_run.size(); irun++) {
      if (sec_scale_index_run[irun] == sclidx) {
        runswithscaleindex.push_back(irun);
      }
    }
    return runswithscaleindex;
  }
  //--------------------------------------------------------------
  std::pair<int,int> ScaleModel::checkLink(const std::pair<int,int>& link,
                                           const std::vector<Run>& runlist) const
  // returns run indices for both ends of the link, first = -1 if not found
  {
    int numruns = runlist.size();
    int r1 = FindRunIndex(link.first, runlist);  //run index
    int r2 = FindRunIndex(link.second, runlist);
    if ((r1 < 0) || (r2 < 0)) {return std::pair<int,int>(-1,-1);} // not found
    return std::pair<int,int>(r1, r2);
  }
  //--------------------------------------------------------------
  void ScaleModel::autoTiles(std::vector<scala::ScaleSpecification>& scaleSpecs,
                         hkl_unmerge_list& hkl_list,
                         phaser_io::Output& output)
  // set automatic TILE settings if appropriate
  // modifies scaleSpecs
  {
    nruns = hkl_list.num_runs();
    std::vector<Run> runlist = hkl_list.RunList();
    for (size_t isp=0; isp<scaleSpecs.size(); isp++) { // loop specifications
      if ((scaleSpecs[isp].detectorscaletype == DetectorScale::AUTOMATIC) ||
          (scaleSpecs[isp].ntilex < 0)) { // no TILE specified for this run
        int irun = scaleSpecs[isp].run;
        bool generalspec = false;
        if (irun < 0) { // the general specification, always the 1st slot
          generalspec = true;
          if (nruns > 1) {
            // find a run which does not have a specific spec
            int krun = -1;
            for (size_t ksp=0; ksp<scaleSpecs.size(); ksp++) {
              if (scaleSpecs[ksp].run < 0) {
                krun = ksp;
                break;
              }
            }
            if (krun >= 0) {
              irun = krun;
            }
          } else {
            irun = 0;
          }
          if (irun < 0) {irun = 0;}
        }
        // Pick up detector size for 1st batch in this run
        int b0 = runlist[irun].BatchSerial0();  // batch serial
        Batch bat0 = hkl_list.Batches()[b0];    // first batch in run
        DetectorType detectortype(bat0);
        //      std::cout << "\nScaleModel::autoTiles DetectorType = "
        //                << detectortype.TypeLabel() << "\n"; //^
        if (detectortype.type() == DetectorType::CCD3x3) {
          scaleSpecs[isp].ntilex = detectortype.NtileX();
          scaleSpecs[isp].ntiley = detectortype.NtileY();
          scaleSpecs[isp].detectorscaletype = scala::DetectorScale::CCD2;
          std::string rn = StringUtil::itos(runlist[irun].RunNumber(),3);
          std::string s = "\nNB The detector type for run "+rn+
            " appears to be a 3x3 tiled CCD detector";
          output.logTab(0,LOGFILE,s);
          if (generalspec) {
            output.logTab(0,LOGFILE,
  " A TILE scale model has therefore been added to the general SCALES specification");
          } else {
            output.logTab(0,LOGFILE,
  " A TILE scale model has therefore been added to the SCALES specification for run "+rn);
          }
          output.logTab(0,LOGFILE,
        " This may be switched off using the command SCALES NOTILE\n");
        } else if (scaleSpecs[isp].detectorscaletype == DetectorScale::AUTOMATIC) {
          // Undefined detector type
          Message::message(Message_fatal
                           ("\nERROR in ScaleModel: undefined detector type for TILE"));
        }
      } // end if !TILE

    } // end loop specs

  }
  //--------------------------------------------------------------
  void ScaleModel::SetupTies(const phaser_io::InputAll& input)
  // Set up all ties from input list and defaults
  {
    sd_rotation = input.TIE_sd_rotation();
    sd_bfactor = input.TIE_sd_bfactor();
    sd_zerob = input.TIE_sd_zerob();
    sd_surface = input.TIE_sd_surface();
    tie_tile = input.TIE_tile();
    SetupTies();
  }
  //--------------------------------------------------------------
  void ScaleModel::SetupTies()
  // Set up all ties from input list and defaults

  {
    ties.clear();  // clear tie list
    nties_rot = nties_bfac = nties_zerob = nties_surf = nties_tiles = 0;

    // Ties on primary scales
    if ((sd_rotation > 0.0) && (nprimaryscale > 0)) {
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
      if (tie_tile.size() > 0 && tie_tile[0] > 0.0) {
        for (int i=0;i<ndetscales;++i) {
          std::vector<Tie> these_ties =
            detector_scales[i].Ties(tie_tile, idxrun_detector[i]);
          ties.insert(ties.end(), these_ties.begin(), these_ties.end());
          nties_tiles += these_ties.size();
        }
      }
    }

    nties = ties.size();
  }
  //--------------------------------------------------------------
  // If true, allow tile corrections to vary azimuthally
  // if false, force to be radially symmetric
  // returns true if anything has changed
  bool ScaleModel::symmetricTiles(const bool& symmetric)
  {
    if (ndetscales <= 0) {return false;} // no tiles
    for (size_t i=0; i<detector_scales.size(); i++) {
      detector_scales[i].setSymmetric(symmetric);
    }
    // Recalculate number of parameters and indices etc
    int npar0 = Nparameters();
    CountParameters();
    SetupTies();  // reset tie list
    if (npar0 != Nparameters()) {
      return true; // Nparameters has changed
    }
    return false;
  }
  //--------------------------------------------------------------
  bool ScaleModel::switchSecondaryScales(const bool& on)
  // switch secondary scales On (true) or Off (false)
  {
    if (on) {
      nsecscales = std::abs(nsecscales);
    } else {
      nsecscales = -std::abs(nsecscales);
    }
    // Recalculate number of parameters and indices etc
    int npar0 = Nparameters();
    CountParameters();
    SetupTies();  // reset tie list
    if (npar0 != Nparameters()) {
      nfreedom = 0; // can't use variances
      return true; // Nparameters has changed
    }
    return false;
  }
  //--------------------------------------------------------------
  void ScaleModel::setBatchReject(const std::vector<bool>& usebatch,
                                  const std::vector<int>& batchnumbers)
  // Set reject list for batches, relevant for BATCH scale mode only (fail if not)
  {
    if (!isAllBatch()) {
      Message::message(Message_fatal
       ("Cannot use REJECT BATCH unless BATCH scaling is used"));
    }
    ASSERT (usebatch.size() == batchnumbers.size());
    for (int irun=0;irun<nruns;irun++) {
      primary_scales[irun].setBatchReject(usebatch, batchnumbers);
      if (relative_bfactors[irun].Number() > 0) {
        relative_bfactors[irun].setBatchReject(usebatch, batchnumbers);
      }
    }
    // Fix up things that may have changed
    CountParameters();
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
    case WAVELENGTH:
      return "WAVELENGTH";
    default:
      return "UNKNOWN";
    }
  }
  //--------------------------------------------------------------
  std::string ScaleModel::SetupScale(const int& irun,
                                     const ScaleSpecification& scaleSpec,
                                     const Run& run,
                                     const ValidScaleModel&  validscalemodel,
                                     phaser_io::Output& output)
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
          s = "\nERROR in ScaleModel: run "+clipper::String(run.RunNumber())+
            " has insufficient information for smooth scaling\n";
          ReportErrors::printWarning(s, "ScaleModelError");
          status = -1;
          return s;
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
     "WARNING: Run %3d has insufficient information for secondary beam scaling\n",
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
     "WARNING: Run %3d has insufficient information for detector (tile) scaling\n",
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
     "WARNING: Run %3d has insufficient information for detector (tile) scaling\n",
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
  void ScaleModel::PrintLayout(phaser_io::Output& output) const
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
      if (nsecscales > 0) {
        if (sec_scale_index_run.at(irun) >= 0) {
          output.logTab(0,LOGFILE,
                        secondary_scales[sec_scale_index_run.at(irun)].format());
        }
      }
      if (detector_scale_index_run.at(irun) >= 0) {
        output.logTab(0,LOGFILE,
                      detector_scales[detector_scale_index_run.at(irun)].format());
      }
    }

    if (bfacnormbatch == -1) {
      // auto
      output.logTab(0,LOGFILE,
		    "\nB-factors will be 'normalised' to the 'best' (largest) relative B");
    } else if (bfacnormbatch == -2) {
      int irun = Max(0, bfacnormrun);
      output.logTab(0,LOGFILE,
		    "\nB-factors will be 'normalised' to the first range in run "+\
		    StringUtil::itos(runnumbers[irun]));
    } else {
      // bfacnormbatch is batch number for normalisation batch
      // NB this option is disabled
      output.logTab(0,LOGFILE,
		    "\nB-factors will be 'normalised' to batch "+ \
		    StringUtil::itos(bfacnormbatch));
    }

    // Assignment of secondary scales to runs
    if (nsecscales > 0 && runnumbers.size() > 1) {
      output.logTab(0,LOGFILE,
    "\nAllocation of secondary (SURFACE) scales to runs (automatic or from LINKs)");
      for (int k=0;k<nsecscales;++k) {
        output.logTab(1,LOGFILE,formatSecondaryrunset(k));
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
        output.logTab(0,LOGFILE, detector_scales[0].formatTies());
      }
    }

    // Parameter SD usage
    std::string s;
    if (parametersdusage == scala::ScaleSpecification::NONE) {
        s = "Parameter variances will not be used for sigma(I) estimates";
    } else if (parametersdusage == scala::ScaleSpecification::DIAGONAL) {
        s = "Parameter variances (DIAGONAL) will be used for sigma(I) estimates";
    } else if (parametersdusage == scala::ScaleSpecification::COVARIANCE) {
        s = "Parameter covariances (COVARIANCE) will be used for sigma(I) estimates";
    }
    output.logTab(0,LOGFILE, "\n"+s);

    output.logTabPrintf(0,LOGFILE,"\n");
  }
  //--------------------------------------------------------------
  std::string ScaleModel::PrintWrappingLines(const std::vector<double>& v,
                                             const std::string& t1,
                                             const std::vector<double>& v2,
                                             const std::string& t2,
                                             const std::vector<int>& n,
                                             const std::string& t3,
                                             const int& fw,
                                             const int& fd)
  // Print wrapping lines:
  //   line 1, values v (double), label t1
  //   line 2, values v2 (int),   label t2  (optional, if size > 0)
  //   line 3, values n (int),    label t3  (optional)
  //   fw  field width
  //   fd  number of decimal points for v
  {
    std::string ss;
    int tlen = Max(Max(t1.size(), t2.size()), t3.size()); // maximum label length
    const int PAGEWIDTH = 100; // characters across
    int nperline = (PAGEWIDTH - tlen - 1)/fw; // number/line
    if (n.size() > 0) {
      ASSERT (v.size() == n.size());
    }
    if (v2.size() > 0) {
      ASSERT (v2.size() == v.size());
    }
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
      // v2
      if (v2.size() > 0) {
        ss += StringUtil::LeftString(t2+":", tlen+1);
        for (int i=i1;i<=i2;++i) {
          ss += StringUtil::ftos(float(v2[i]), fw, fd);
        }
        ss += "\n"; // end of line
      }
      // n
      if (n.size() > 0) {
        ss += StringUtil::LeftString(t3+":", tlen+1);
        for (int i=i1;i<=i2;++i) {
          ss += clipper::String(n[i], fw);
        }
        ss += "\n"; // end of line
      }
      i1 = i2+1;
    }
    return ss;
  }
  //--------------------------------------------------------------
  std::string ScaleModel::PrintWrappingLinesWithSD(const std::vector<double>& v,
                                                   const std::string& t1,
                                                   const std::vector<double>& sds,
                                                   const int& fw,
                                                   const int& fd)
  // Print wrapping lines as val(SD):
  //   values v (double), label t1, sd(v) in sds
  //   fw  field width
  //   fd  number of decimal points for v
  // Assume sd << v
  {
    ASSERT (sds.size() == v.size());
    int nitem = v.size();
    std::string ss;
    int tlen = t1.size(); // label length
    const int PAGEWIDTH = 100; // characters across
    // each item will be of form xxx.ddd(sd)

    // guess total field width to get number per line
    int totfw = fw + 2 + fd; // allowing for brackets
    int nperline = (PAGEWIDTH - tlen - 1)/totfw; // number/line
    int i1 = 0; // 1st item in line
    while (i1 < nitem) {
      int i2 = Min(i1 + nperline, nitem) - 1; // last item in line
      ss += StringUtil::LeftString(t1+":", tlen+1);
      for (int i=i1;i<=i2;++i) {
        ss += StringUtil::valueSD(v[i], sds[i], totfw, fw, fd);
      }
      ss += "\n"; // end of line
      i1 = i2+1;
    }
    return ss;
  }
  //--------------------------------------------------------------
  std::vector<double> ScaleModel::extractSDs(const int& idxsd, const size_t& nsd) const
  // extract nsd SDs from variance, beginning at index idxsd
  {
    std::vector<double> sds(nsd, 0.0);
    for (size_t i=0;i<nsd;++i) {
      if (varpar[idxsd + i] > 0.0) {
        sds[i] = sqrt(varpar[idxsd + i]);
      }
    }
    return sds;
  }
  //--------------------------------------------------------------
  std::string ScaleModel::formatSecondaryrunset(const int& scaleset) const
  // format as "Set <scaleset>, run[s]: <runnumbers>"
  {
    // Get list of runs corresponding to this secondary scale set
    std::vector<int> runswithscaleindex = runsWithScaleIndex(scaleset);
    std::string s = "";
    bool first = true;
    std::string runtext = "run";
    for (size_t k=0; k<runswithscaleindex.size(); k++) {
      if (!first) {
        s += ", ";
        runtext = "runs";
      }
      first = false;
      s +=
        StringUtil::Strip(StringUtil::itos(runnumbers[runswithscaleindex[k]],3));
    }
    s = "Scale set "+StringUtil::itos(scaleset+1,2)+" used for "+runtext+" "+s;
    return s;
  }
  //--------------------------------------------------------------
  void ScaleModel::PrintScales(phaser_io::Output& output) const
  // Print all scale parameters, with SDs if available
  {
    const int nperline = 10;
    std::vector<double> sds;

    output.logTab(0,LOGFILE,"\nScale parameters:\n");
    for (size_t irun=0;irun<runnumbers.size();++irun) {
      output.logTabPrintf(0,LOGFILE,
                          "\nRun %5d\n", runnumbers[irun]);
      // Primary
      output.logTab(0,LOGFILE,"Primary scales and number of observations");
      std::vector<double> pscales = primary_scales[irun].Scales();
      for (size_t i=0;i<pscales.size();++i) {pscales[i]=1.0/pscales[i];}
      std::vector<int> nobsPar = primary_scales[irun].Nobservations();
      if (nfreedom <= 0) {
        sds.clear();
        output.logTab(0,LOGFILE,
                      PrintWrappingLines(pscales, "Scales", sds, "",
                                         nobsPar, "Nobs",9,3));
      } else {
        if (nprimaryscale == 0) {
          ASSERT (pscales.size() == 1);
          sds.assign(1, 0.0);
        } else {
          sds = extractSDs(idxrun_primary_scales[irun], pscales.size());
          for (size_t k=0; k<sds.size(); k++) {
            // sd(k) = sd(1/g) = sd(g) * k^2
            sds[k] *= pscales[k]*pscales[k];
          }
        }
        output.logTab(0,LOGFILE,
                      PrintWrappingLines(pscales, "Scales", sds, "Sd",
                                         nobsPar, "Nobs",9,3));
      }
      // B-factors
      std::vector<double> bfacs = relative_bfactors[irun].Bfactors();
      nobsPar = relative_bfactors[irun].Nobservations();
      if (bfacs.size() > 0) {
        output.logTab(0,LOGFILE,"\nRelative B-factors and number of observations");
      if (nfreedom <= 0) {
        sds.clear();
        output.logTab(0,LOGFILE,
                      PrintWrappingLines(bfacs, "B-factors", sds, "",
                                         nobsPar, "Nobs",9,3));
      } else {
        sds = extractSDs(idxrun_bfactors[irun], bfacs.size());
        output.logTab(0,LOGFILE,
                      PrintWrappingLines(bfacs, "B-factors", sds, "Sd",
                                         nobsPar, "Nobs",9,3));
      }
        output.logTabPrintf(0,LOGFILE,"\n");
      }
    } // end run loop
    if (nsecscales > 0) {   // Secondary
      output.logTab(0,LOGFILE,"Secondary scales");
      for (int j=0;j<nsecscales;++j) {
        std::vector<double> secsclpar = secondary_scales[j].Coefficients();
        sds.clear();
        if (nfreedom > 0) {
          sds = extractSDs(idxrun_secondary[j], secsclpar.size());
        };
        if (runnumbers.size() > 1) {
          output.logTab(0,LOGFILE,"\n"+formatSecondaryrunset(j));
        } else {
          output.logTabPrintf(0,LOGFILE,"\nScale set %3d\n",j+1);
        }
        std::vector<int> ndummy;
        if (nfreedom <= 0) {
          output.logTab(0,LOGFILE,
                        PrintWrappingLines(secsclpar, "Coefficient", sds, "",
                                           ndummy, "",8,4));
        } else {
          output.logTab(0,LOGFILE,
                        PrintWrappingLinesWithSD(secsclpar,
                                                 "Coefficient(Sd)", sds, 8,4));
        }
      }
    } // End Secondary
    if (ndetscales > 0) { // Tiles
      output.logTabPrintf(0,LOGFILE,"\n\n");
      for (int i=0;i<ndetscales;++i) {
        sds.clear();
        if (nfreedom > 0) {
          sds = extractSDs(idxrun_detector[i], detector_scales[i].Number());
        };
        output.logTab(0,LOGFILE,detector_scales[i].formatparameters(sds));
      }
    }
    output.logTabPrintf(0,LOGFILE,"\n\n");
  }
  //--------------------------------------------------------------
  void ScaleModel::PrintSecondaryCorrections(phaser_io::Output& output) const
  // Print secondary corrections as 2D array
  {
    if (nsecscales > 0) {   // Secondary
      output.logTab(0,LOGFILE,
                    "\nSecondary scale corrections");
      output.logTab(0,LOGFILE,
                    "\nCalculated for polar angles of theta (colatitude from 0 at N pole) and phi (longitude)\n");



      double angleinterval = 10.0;
      int nskip = 2;  // for log file print
      // Phi values
      int nphi = Nint(360.0/angleinterval) + 1;
      std::vector<double> phivalues(nphi);
      for (size_t k=0; k<phivalues.size(); k++) {
        phivalues[k] = k * angleinterval;
      }
      // Theta values
      int ntheta = Nint(180.0/angleinterval) + 1;
      std::vector<double> thetavalues(ntheta);
      for (size_t k=0; k<thetavalues.size(); k++) {
        thetavalues[k] = k * angleinterval;
      }

      std::string line;
      for (int j=0;j<nsecscales;++j) {
        output.logTabPrintf(0,LOGFILE,"\nSecondary scale number %3d\n", j+1);
        line = " Phi ";
        for (size_t kp=0; kp<phivalues.size(); kp++) {
          if (kp%nskip == 0) {
            line += StringUtil::ftos(phivalues[kp], 5, 0);
          }
        }
        output.logTab(0,LOGFILE, line);
        output.logTab(0,LOGFILE, " Theta");

        output.logTab(0,LXML,"<SecondaryCorrection>");
        output.logTab(0,LXML,
                      StringUtil::MakeXMLtag("ScaleNumber", j+1));
        output.logTab(0,LXML,
                      StringUtil::MakeXMLtag("PhiValues",
                      StringUtil::FormatSaveVector(phivalues)));

        for (size_t kt=0; kt<thetavalues.size(); kt++) {
          line = "";
          if (kt%nskip == 0) {
            line = StringUtil::ftos(thetavalues[kt], 5, 0)+" ";
          }
          output.logTab(0,LXML,"<secscales>");
          output.logTab(0,LXML,
                      StringUtil::MakeXMLtag("Theta", thetavalues[kt]));
          std::vector<double> corrections(phivalues.size());

          for (size_t kp=0; kp<phivalues.size(); kp++) {
            corrections[kp] =
              secondary_scales[j].Scale(clipper::Util::d2rad(thetavalues[kt]),
                                        clipper::Util::d2rad(phivalues[kp]));
            if (corrections[kp] != 0.0) {
              corrections[kp] = 1.0/corrections[kp];
            }

            if ((kp%nskip == 0) && (kt%nskip == 0)) {
              line += StringUtil::ftos(corrections[kp],5,2);
            }
          }
          output.logTab(0,LOGFILE, line);
          output.logTab(0,LXML,
                        StringUtil::MakeXMLtag("corrections",
                       StringUtil::FormatSaveVector(corrections)));

          output.logTab(0,LXML,"</secscales>");

        }
        output.logTab(0,LXML,"</SecondaryCorrection>");
      }
    }
  }
  //--------------------------------------------------------------
  void ScaleModel::CountParameters()
  // private
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
    if (nsecscales >= 0) {
      idxrun_secondary.resize(nsecscales);
    } else {
      idxrun_secondary.clear();
    }
    idxrun_detector.resize(ndetscales);

    // Primary
    // if SCALES CONSTANT and one run, then set nparameters = 0, unles refine reference
    if ((nruns == 1) && (primary_scales[0].Number() <= 1) && !refinereference) {
      nparameters = 0;
      nprimaryscale = 0;
      idxrun_primary_scales[0] = 0;
    } else {
      for (int irun=0;irun<nruns;irun++) {
        idxrun_primary_scales[irun] = nparameters;   // index to 1st primary parameter
        nparameters += primary_scales[irun].Number();
        nprimaryscale += primary_scales[irun].Number();
      }
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
        idxrun_detector[i] = nparameters;   // index to 1st detector parameter
        ntilescale += detector_scales[i].Number();
        nparameters += detector_scales[i].Number();
      }
    }
    // Wavelength (Chebyshev) normalization
    nwavscale = wavelength_scale.Number();
    if (nwavscale > 0) {
      idxwavscale = nparameters;
      nparameters += nwavscale;
    }
  }
  //--------------------------------------------------------------
  std::vector<double> ScaleModel::GetParameters() const
  // Get vector of parameters
  // Order of parameters:
  //   1. all primary scale parameters (nprimaryscale) (this may == 0)
  //   2. all B-factors (nbfactors)
  //   3. secondary parameters (nsecondaryscale)
  //   4. detector (ntilescale)
  //   5. ..
  {
    std::vector<double> params;

    if (nprimaryscale > 0) {
      for (int irun=0;irun<nruns;irun++) {
        // Append group of scales
        std::vector<double> pscales = primary_scales[irun].Scales();
        params.insert(params.end(), pscales.begin(), pscales.end());
      }
      ASSERT (int(params.size()) == nprimaryscale);
    }
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
      //std::cout <<  secondary_scales[i].Number() << " " << params.size() << "\n"; //^^
    }
    ASSERT (int(params.size()) == nprimaryscale+nbfactors+nsecondaryscale);
    // Detector
    for (int i=0;i<ndetscales;++i) {
      std::vector<double> detsclpar = detector_scales[i].Parameters();
      params.insert(params.end(), detsclpar.begin(), detsclpar.end());
      //      std::cout <<  detector_scales[i].Number() << " " << params.size() << "\n";
    }
    ASSERT (int(params.size()) == nprimaryscale+nbfactors+nsecondaryscale+ntilescale);
    // Wavelength (Chebyshev) normalization
    if (nwavscale > 0) {
      std::vector<double> wavpar = wavelength_scale.Coefficients();
      params.insert(params.end(), wavpar.begin(), wavpar.end());
    }
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
    for (int i=0;i<nwavscale;++i) {partype[++k] = ScaleModel::WAVELENGTH;}
    ASSERT (++k == nparameters);
    return partype;
  }
  //--------------------------------------------------------------
  ScaleModel::ScaleParameterType ScaleModel::GetParameterType(const int& Ipar) const
  // get type for a parameter
  {
    if (Ipar < 0 || Ipar >= nparameters) {
      clipper::Message::message(Message_fatal
                                ("GetParameterType: parameter number out of range"+
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
    if (Ipar < nprimaryscale+nbfactors+nsecondaryscale+ntilescale+nwavscale) {
      return ScaleModel::WAVELENGTH;
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
    for (idetsc=0;idetsc<int(idxrun_detector.size());++idetsc) {
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
                                 const std::vector<int>& Nobs,
                                 const bool& donormalise)
  {
    std::vector<double> pars(params.size());
    for (size_t i=0;i<params.size();++i) {pars[i] = params[i];}
    SetParameters(pars, Nobs, donormalise);
  }
  //--------------------------------------------------------------
  // Set all parameters from vector
  void ScaleModel::SetParameters(const std::vector<double>& params,
                                 const std::vector<int>& Nobs,
                                 const bool& donormalise)
  // optional normalisation
  {
    ASSERT (int(params.size()) == nparameters);
    std::vector<double>::const_iterator pos1 = params.begin();  // start of range
    std::vector<double>::const_iterator pos2;                   // end of range
    std::vector<int>::const_iterator posn1 = Nobs.begin();   // start of range
    std::vector<int>::const_iterator posn2;                  // end of range
    if (nprimaryscale > 0) {
      for (int irun=0;irun<nruns;irun++) {
        // Extract group of scales
        pos2 = pos1 + primary_scales[irun].Number();
        primary_scales[irun].StoreScales(std::vector<double>(pos1, pos2));
        pos1 = pos2;
        posn2 = posn1 + primary_scales[irun].Number();
        primary_scales[irun].StoreNobservations(std::vector<int>(posn1, posn2));
        posn1 = posn2;
      }
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
    if (nwavscale > 0) {
      // Wavelength (Chebyshev) normalization
      pos2 = pos1 + nwavscale;
      wavelength_scale.StoreCoefficients(std::vector<double>(pos1, pos2));
      pos1 = pos2;
      posn1 += nwavscale;  // no per-observation count for wavelength params
    }
    if (donormalise) {NormaliseParameters();}
  }
  //--------------------------------------------------------------
  //
  void ScaleModel::clearCounts()
    // Clear observation counts at beginning of cycle, mainly relevant for tiles
  {
    for (int irun=0;irun<nruns;irun++) {
      primary_scales[irun].StoreNobservations(std::vector<int>(primary_scales[irun].Number(), 0));
    }
    for (int irun=0;irun<nruns;irun++) {
      relative_bfactors[irun].StoreNobservations(std::vector<int>(relative_bfactors[irun].Number(), 0));
    }
    for (int i=0;i<nsecscales;++i) {
      // Secondary
      secondary_scales[i].StoreNobservations(std::vector<int>(secondary_scales[i].Number(), 0));
    }
    for (int i=0;i<ndetscales;++i) {
      // Detector
      detector_scales[i].StoreNobservations(std::vector<int>(detector_scales[i].Number(), 0));
      detector_scales[i].clearCounts();
    }
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
    //  4) ... or Nscales = 3, in which case Nintervals = 2
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

      if (primary_scales[ir].IsBatchScale()) {  // Batch scale  -------
        ASSERT (primary_scales[ir].Number() == primary_scales[ir].Nintervals());
        for (int i=0;i<primary_scales[ir].Nintervals();++i) {
          gsclrun[i] = gscales[jsr];
          nobsrun[i] = numobsrotrange[jsr];
          jsr++;
        }
      } else { // Smooth scale  -------
        if (primary_scales[ir].Number() <= 2) { // 1 or 2 scales, 1 interval
          ASSERT (primary_scales[ir].Nintervals() == 1);
          for (int i=0;i<primary_scales[ir].Number();++i) {
            gsclrun[i] = gscales[jsr];
            nobsrun[i] = numobsrotrange[jsr];
          }
          jsr++;
        } else if (primary_scales[ir].Number() == 3) { // 3 scales, 2 intervals
          ASSERT (primary_scales[ir].Nintervals() == 2);
          // Scales at 0,1,2 average for middle one
          gsclrun[0] = gscales[jsr];
          nobsrun[0] = numobsrotrange[jsr];
          gsclrun[2] = gscales.at(jsr+1);
          nobsrun[2] = numobsrotrange.at(jsr+1);
          gsclrun[1] = 0.5*(gscales[jsr] + gscales.at(jsr+1)); // average scale
          nobsrun[1] = (numobsrotrange[jsr] + numobsrotrange.at(jsr+1))/2; // average number
          jsr += 2;
        } else { // >3 scales n-2 intervals
          ASSERT (primary_scales[ir].Number() == primary_scales[ir].Nintervals()+2);
          for (int i=0;i<primary_scales[ir].Nintervals();++i) {
            j = i+1;
            // extra scale at start, duplicate of 1st interval
            if (i == 0) {
              gsclrun[i] = gscales[jsr];
              nobsrun[i] = numobsrotrange[jsr];
            }
            gsclrun[j] = gscales[jsr];
            nobsrun[j] = numobsrotrange[jsr];
            jsr++;
          } // end loop scales intervals
          // extra scale at end: j is one beyond intitial scales array
          gsclrun[j+1] = gsclrun[j];
          nobsrun[j+1] = numobsrotrange[j-1];
          j++;
          ASSERT (j+1 == primary_scales[ir].Number());
        }
      }
      primary_scales[ir].StoreScales(gsclrun);
      primary_scales[ir].StoreNobservations(nobsrun);
    } // end loop runs
  }
  //--------------------------------------------------------------
  // Return true if model is refinable, ie not just one scale and one B-factor
  bool ScaleModel::IsRefinable() const
  {
    if (status < 0) {return false;}  // insufficient information
    return (nprimaryscale > 1) || (nbfactors > 1) ||
      (nsecondaryscale > 0) || (ntilescale > 0) || (nwavscale > 0);
  }
  //--------------------------------------------------------------
  // return reason for being not refinable:
  //   insufficient information or only one parameter, or blank if it is
    std::string ScaleModel::whyNotRefineable() const
  {
    if (IsRefinable()) {return "";}
    if (status < 0) {
      return "insufficient information in file";
    }  // insufficient information
    if (nparameters <= 1) {
      return "only one parameter";
    }
    return "";
  }
  //--------------------------------------------------------------
  // Number of normalisation parameters: usually 2 (k,B);
  //   1 if fixed Bfactors for any run; = 0 if not refineable
  int ScaleModel::Nfilter() const
  {
    int nfilter = 0;
    if (IsRefinable()) {
      if (normalisebfac) {
        nfilter = 2;
      } else {
        nfilter = 1;
      }
    }
    return nfilter;
  }
  //--------------------------------------------------------------
  bool ScaleModel::isAllBatch() const
  // Return true if BATCH scales for all runs
  {
    bool allbatch = true;
    for (int irun=0;irun<nruns;irun++) {
      if (!primary_scales[irun].IsBatchScale()) {
        allbatch = false;
        break;
      }
      if (relative_bfactors[irun].Number() > 0) {
        if (!relative_bfactors[irun].IsBatchBfactor()) {
          allbatch = false;
          break;
        }
      }
    }
    return allbatch;
  }
  //--------------------------------------------------------------
  double ScaleModel::ScaleObs(observation& obs, const Rtype& invresolsq,
                              const bool& onlyUseSingletons) const
  // Scale observation, returns scale applied
  // if onlyUseSingletons true, do not attempt to apply scales to overlaps
  {
    int irun = obs.run();
    double g, varg;
    if (nfreedom <= 0) {
      g = ScaleFactor(irun, obs, invresolsq); // main scale
      obs.SetGscale(g);
    } else {
      std::vector<double> dghldp;
      g =  ScaleFactorDeriv(irun, obs, invresolsq, dghldp);
      varg = VarScale(dghldp);
      obs.SetGscaleVar(g, varg);
    }
    if (g <= 0.0) {
      std::string message = "ScaleObs: non-positive scale "+StringUtil::ftos(g)+
        " "+obs.hkl_original().format();
      ReportErrors::printWarning(message,"NegativeScale",false);
    }

    // For multiple lattice observations, get appropriate inverse scale factors
    if (!obs.IsSingleton() && !onlyUseSingletons) {
      int mainlatnum = obs.MainLatticeNumber();  // for main hkl
      std::vector<LatticeIndexInfo> lathkl = obs.lathkl();
      //^
      //      if (lathkl.size() > 0) {
      //        std::cout << "ScaleObs " << obs.hkl_original().format() <<
      //          " Mainlat "<< mainlatnum <<
      //          " Batch " << obs.Batch() << " gm " << g <<"\n";
      //      }
      //^-

      for (size_t l=0; l<lathkl.size(); l++) {
        int latnum = lathkl[l].latnum;
        if (latnum > 0) {
          double grel = 0.0;   // relative scale, g(latticeN)/g(latticeMain)
          double sdgrel = 0.0;
          int jscale = idxrunlattice[latnum-1];
          if (jscale >= 0) {
            if (nfreedom <= 0) {
              grel = ScaleFactor(jscale, obs, invresolsq) / g;
            } else {
              std::vector<double> dghldp;
              double glat = ScaleFactorDeriv(jscale, obs, invresolsq, dghldp);
              double varglat = VarScale(dghldp);
              grel = glat / g;
              // Var(grel)/grel = Var(glat)/glat + Var(g)/g
              sdgrel = sqrt(grel*(varglat/glat + varg/g));
            }
          }
          lathkl[l].gscale = grel;      // relative lattice fraction from scales
          lathkl[l].sdgscale = sdgrel;  // sd(relative lattice fraction)

          //^
          //      std::cout << "  lattice " << latnum << " " <<
          //        lathkl[l].hkl.format() << " g " << glat <<
          //        " relative g " << lathkl[l].gscale <<"\n";
          //^-
        }
      }
      obs.StoreLathkl(lathkl);
    }
    return g;
  }
  //--------------------------------------------------------------
  ////  double ScaleModel::VarScale(const std::vector<double>& dghldp) const
  double ScaleModel::VarScale(const std::vector<double>& dghldp) const
  // Variance(gscale)
  {
    double varg = 0.0;
    if (parametersdusage == scala::ScaleSpecification::DIAGONAL) {
      ASSERT (nfreedom > 0);
      ASSERT (dghldp.size() == varpar.size());
      // Diagonal approximation: Var(g) = Sum (dghldp[i]^2 * Var(p[i]))
      for (size_t i=0; i<dghldp.size(); i++) {
        if (dghldp[i] != 0.0) {
          varg += dghldp[i] * dghldp[i] * varpar[i];
          //      if (DEBUG) { //^^
          //        std::string ptype = ScaleParameterTypeString(GetParameterType()[i]);
          //        std::cout << "Par " << i <<", type "<<ptype
          //                  <<", varpar "<<varpar[i]<<", dghldp "<<dghldp[i]
          //                  << ", Delvarg "<<dghldp[i] * dghldp[i] * varpar[i]
          //                  <<", varg "<<varg<<"\n";
          //      }
        }
      }
    } else if (parametersdusage == scala::ScaleSpecification::COVARIANCE) {
      ASSERT (nfreedom > 0);
      ASSERT (dghldp.size() == varpar.size());
      // Full covariance matrix:
      //   Var(g) = (dghldp)T [VC] (dghldp)
      double v;
      for (size_t k=0; k<dghldp.size(); k++) {
        if (dghldp[k] != 0.0) {
          v = 0.0;
          for (size_t j=0; j<dghldp.size(); j++) {
            if (dghldp[j] != 0.0) {
              v += VC(k,j) * dghldp[j];  // {[VC](dghldp)}[k]
            }
          }
          varg += dghldp[k] * v;
        }
      }
    }
    return varg;
  }
  //--------------------------------------------------------------
  double ScaleModel::ScaleFactor(const int& jscale,
                                 const observation& obs,
                                 const Rtype& invresolsq) const
  // Returns scale for observation, using scale set jscale (== irun for main observation
  {
    double g = 1.0;
    double ps;    // Primary scale
    if (primary_scales[jscale].IsBatchScale()) {
      ps = primary_scales[jscale].Scale(obs.Batch());
    } else {
      ps = primary_scales[jscale].Scale(obs.phi());
    }

    // B-factor scale
    double bs = 1.0;

    if (relative_bfactors[jscale].Number() > 0) {
      if (relative_bfactors[jscale].IsBatchBfactor()) {
        bs = relative_bfactors[jscale].BfactorScale(obs.Batch(), invresolsq);
      } else {
        bs = relative_bfactors[jscale].BfactorScale(obs.time(), invresolsq);
      }
    }
    // Secondary   FIXME (why?)
    double ss = 1.0;
    if (nsecscales > 0) {
      if (sec_scale_index_run[jscale] >= 0) {
        double thetap, phip;
        obs.GetS2(thetap, phip);
        ss = secondary_scales[sec_scale_index_run[jscale]].Scale(thetap, phip);
        if (ss <= 0.0) {
          std::cout <<"neg scale " << ss << " "<<thetap<<" "<<phip<<"\n"; //^^
          negativeSecScale++; // record occurance of negative scale
          ss = 0.5;
        }
      }
    }
    // Detector
    double ds = 1.0;
    if (ndetscales > 0) {
      if (detector_scale_index_run[jscale] >= 0) {
        ds = detector_scales[detector_scale_index_run[jscale]].Scale(obs.XYdet());
      }
    }
    // Wavelength normalization (Chebyshev refinable * GPR fixed correction)
    double ws = 1.0;
    if (nwavscale > 0) {
      ws = wavelength_scale.Scale(obs.lambda());
    }
    if (gpr_scale.IsActive()) {
      ws *= gpr_scale.Scale(obs.lambda());
    }
    g = ps*bs*ss*ds*ws;
    //^^^
    ////    if (g <= 0.0) {
    if (g <= 1.0e-8) {
      std::cout <<"g too small "
                << g <<" "<< ps<<" "<<bs<<" "<<ss<<" "<<ds<<" "<<ws<<"\n";
    }
    return g;
  }
  //--------------------------------------------------------------
  double ScaleModel::ScaleObs(observation& obs, const Rtype& invresolsq,
                              std::vector<double>& dghldp) const
  // Scale observation, returns scale applied and
  // partial derivative vector d(ghl)/dp
  {
    ASSERT (obs.IsSingleton()); // otherwise trouble!
                                // called here from refinement, which can't handle overlaps
    // Run
    int irun = obs.run();
    double g = ScaleFactorDeriv(irun, obs, invresolsq, dghldp);
    obs.SetGscale(g);
    return g;
  }
  //--------------------------------------------------------------
  double ScaleModel::ScaleFactorDeriv(const int& jscale,
                                      observation& obs,
                                      const Rtype& invresolsq,
                                      std::vector<double>& dghldp) const
  // Returns scale for observation, using scale set jscale (== irun for main observation
  // and partial derivative vector d(ghl)/dp
  {
    dghldp.assign(nparameters, 0.0);
    double ps;    // Primary scale
    std::vector<double> dgdpm;  // derivatives for this run only
    // Special for a single scale for one run only
    if (nprimaryscale <= 0) {
      ps = 1.0;
      dgdpm.clear();
    } else {
      if (primary_scales[jscale].IsBatchScale()) {
        ps = primary_scales[jscale].ScaleDeriv(obs.Batch(), dgdpm);
      } else {
        ps = primary_scales[jscale].ScaleDeriv(obs.phi(), dgdpm);
      }
    }
    //^^^
    if (ps <= 0.0) {
      std::cout <<"ps <= 0 " << ps <<" "<<obs.phi()<<"\n";
    } //-^^

    // B-factor scale
    std::vector<double> dgdB;  // derivatives for this run only
    double bs = 1.0;
    if (relative_bfactors[jscale].Number() > 0) {
      if (relative_bfactors[jscale].IsBatchBfactor()) {
        bs = relative_bfactors[jscale].BfactorScaleDeriv(obs.Batch(), invresolsq, dgdB);
      } else {
        bs = relative_bfactors[jscale].BfactorScaleDeriv(obs.time(), invresolsq, dgdB);
      }
    }

    // Secondary
    double ss = 1.0;
    std::vector<double> dgds;  // derivatives for this run only
    if (nsecscales > 0) {
      if (sec_scale_index_run[jscale] >= 0) {
        // Return scale & derivatives for secondary beam directions
        // polar angles thetap, phip (calculated in ScaleModel constructor)
        double thetap, phip;
        obs.GetS2(thetap, phip);
        int k = sec_scale_index_run[jscale];
        ss = secondary_scales[k].ScaleDeriv(thetap, phip, dgds);
        if (ss <= 0.0) {
          negativeSecScale++; // record occurance of negative scale
          ss = 0.5;
        }
      }
    }

    // Detector
    double ds = 1.0;
    std::vector<double> dgdd;  // derivatives for this run only
    if (ndetscales > 0) {
      if (detector_scale_index_run[jscale] >= 0) {
        ds = detector_scales[detector_scale_index_run[jscale]].
          ScaleDeriv(obs.XYdet(), dgdd);
      }
    }

    // Wavelength normalization.  ws = (Chebyshev refinable) * (GPR fixed).
    // The GPR factor has no refinable parameters; it multiplies every scale
    // component (and the Chebyshev derivatives) just like a constant.
    double ws = 1.0;
    std::vector<double> dgdw;  // derivatives for wavelength params
    if (nwavscale > 0) {
      ws = wavelength_scale.ScaleDeriv(obs.lambda(), dgdw);
    }
    double wsgpr = gpr_scale.IsActive() ? gpr_scale.Scale(obs.lambda()) : 1.0;
    ws *= wsgpr;

    // dghl/dp = dg(primary)/dp * bs * ss * ds * ws
    if (dgdpm.size() > 0) {
      for (size_t i=0;i<dgdpm.size();++i) {
        dgdpm[i] *= bs * ss * ds * ws;
      }
      std::copy(dgdpm.begin(), dgdpm.end(),
                dghldp.begin() + idxrun_primary_scales[jscale]);
    }

    // dghl/dp = dg(B)/dp * ps * ss * ds * ws
    for (size_t i=0;i<dgdB.size();++i) {
      dgdB[i] *= ps * ss * ds * ws;
    }
    std::copy(dgdB.begin(), dgdB.end(),
              dghldp.begin()+idxrun_bfactors[jscale]);

    if (nsecscales > 0) {
      // dghl/dp = dg(sec)/dp * ps * bs * ds * ws
      for (size_t i=0;i<dgds.size();++i) {
        dgds[i] *= ps * bs * ds * ws;
      }
      int k = sec_scale_index_run[jscale];
      std::copy(dgds.begin(), dgds.end(), dghldp.begin()+idxrun_secondary[k]);
    }

    if (ndetscales > 0) {
      // dghl/dp = dg(det)/dp * ps * bs * ss * ws
      for (size_t i=0;i<dgdd.size();++i) {
        dgdd[i] *= ps * bs * ss * ws;
      }
      int k = detector_scale_index_run[jscale];
      std::copy(dgdd.begin(), dgdd.end(), dghldp.begin()+idxrun_detector[k]);
    }

    if (nwavscale > 0) {
      // dghl/dp = dg(wav)/dp * ps * bs * ss * ds * wsgpr
      // (wsgpr is the fixed GPR factor; 1.0 when GPR inactive)
      for (size_t i=0;i<dgdw.size();++i) {
        dgdw[i] *= ps * bs * ss * ds * wsgpr;
      }
      std::copy(dgdw.begin(), dgdw.end(), dghldp.begin()+idxwavscale);
    }

    // In wavelength-only mode, zero all non-wavelength derivatives so the
    // least-squares engine moves only the Chebyshev coefficients
    if (wavelength_only_mode_) {
      for (int i = 0; i < nparameters; ++i) {
        if (GetParameterType(i) != WAVELENGTH) dghldp[i] = 0.0;
      }
    }

    double g = ps*bs*ss*ds*ws;
    //^^^
    ////    if (g <= 0.0) {
    if (g <= 1.0e-8) {
      std::cout <<"g too small "
                << g <<" "<< ps<<" "<<bs<<" "<<ss<<" "<<ds<<" "<<ws<<"\n";
    }
    obs.SetGscale(g);
    return g;
  }
  //--------------------------------------------------------------
  void ScaleModel::PrintWavelengthNormalization(phaser_io::Output& output) const
  {
    if (!HasWavelengthScale()) return;
    output.logTab(0, LOGFILE, wavelength_scale.PrintNormalization(12));
  }
  //--------------------------------------------------------------
  void ScaleModel::FitGPRWavelength(const std::vector<double>& lambdas,
                                    const std::vector<double>& logratios,
                                    const std::vector<double>& weights,
                                    phaser_io::Output& output)
  // Fit the GP wavelength normalization from per-observation samples and
  // store it as a fixed multiplicative correction (applied in ScaleFactor).
  {
    if (!gpr_requested) return;
    std::string fitlog;
    int nb = gpr_scale.Fit(lambdas, logratios, weights,
                           gpr_control, gpr_lambda_ref, fitlog);
    output.logTab(0, LOGFILE, fitlog);
    if (nb == 0) {
      output.logTab(0, LOGFILE,
        "GPR wavelength normalization fit failed; no correction applied\n");
    }
  }
  //--------------------------------------------------------------
  void ScaleModel::PrintGPRWavelengthNormalization(phaser_io::Output& output) const
  {
    if (!gpr_scale.IsActive()) return;
    output.logTab(0, LOGFILE, gpr_scale.PrintNormalization(12));
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
    if (nprimaryscale == 0) return;  // nothing to normalise
    // Scales for normalisation run
    ASSERT (scalenormbatch >= 0);
    if (nprimaryscale > 0) {
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
    }
    // B-factors
    if (normalisebfac) {  // only if all runs have variable B-factors
      double bfnorm = -1000000.;
      if (bfacnormbatch == -2) {
	// first batch in specified run bfacnormrun or first run
	int bfrun = bfacnormrun;
	if (bfacnormrun < 0) {
	  bfrun = 0;
	}
        bfnorm = relative_bfactors[bfrun].Bfactors()[0];
      } else if (bfacnormbatch >= 0) {
        // Normalisation batch specified
        bfnorm = relative_bfactors[bfacnormrun].Bfactors()[bfacnormbatchserial];
      } else {
        // Find largest Bfactor
        for (int irun=0;irun<nruns;irun++) {
          //  B-factors for this run
          std::vector<double> bfacs = relative_bfactors[irun].Bfactors();
	  std::vector<int> nobspar = relative_bfactors[irun].Nobservations();
	  ASSERT (bfacs.size() == nobspar.size());
          // Normalisation on "best" batch: if smoothed Bfactors && > 2, omit first & last
          int i1 = 0;
          int i2 = bfacs.size();
	  int nthreshold = 0;
          if (!relative_bfactors[irun].IsBatchBfactor() &&
              i2 > 2) {
	    // Also omit parameters with relatively few observations
	    // This may help to limit the problems with big gaps in the data
	    Median<int> medianNobs(nobspar);
	    const double FRACTIONOFMEDIAN = 0.6;
	    nthreshold = int(double(medianNobs.median())*FRACTIONOFMEDIAN);
            // ... but for smoothed values add in the mean of first & last pairs
	    // Always include the first pair
            double bf = 0.5 * (bfacs[0] + bfacs[1]);
	    bfnorm = Max(bfnorm, bf);
            bf = 0.5 * (bfacs[i2-2] + bfacs[i2-1]);
	    if (nobspar[1] > nthreshold) {bfnorm = Max(bfnorm, bf);}
            i1 = 1;  // now omit 1st & last
            i2--;
          }
          for (int i=i1;i<i2;++i) {
	    if (nobspar[i] > nthreshold) {
	      bfnorm = Max(bfnorm, bfacs[i]);
	    }
          }
        } // end loop runs
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
    } // bfactors
    // Fix up secondary scales if there is a negative scale for any observation
    if (negativeSecScale > 0) {
      fixupSecondaryScales();
    }
  }
  //--------------------------------------------------------------
  void ScaleModel::fixupSecondaryScales()
  // Fix up secondary scales if there is a negative scale for any observation
  {
    if (nsecscales == 0) {return;}
    //std::cout << "ScaleModel::SetParameters negative sec scale\n";
    // scale down all coefficients
    const double FACTOR = 0.5;
    for (int i=0;i<nsecscales;++i) {
      std::vector<double> pars = secondary_scales[i].Coefficients();
      for (size_t k=0; k<pars.size(); k++) {
        pars[k] *= FACTOR;
      }
      secondary_scales[i].StoreCoefficients(pars); // store back
    }
    negativeSecScaleOccurred = Max(negativeSecScale, negativeSecScaleOccurred);
    negativeSecScale = 0;  // clear flag
  }
  //--------------------------------------------------------------
  void ScaleModel::clearNegativeSecScaleFlag()
  // call at beginning of each refinement cycle to clear negative sec scale flag
  // and store previous value for end
  {
    negativeSecScaleLast = negativeSecScale; // record previous value
    negativeSecScale = 0;
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
      // A relative B-factor unbounded
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
    case ScaleModel::WAVELENGTH:
      // Log parameterisation: a_k are unbounded (exp ensures positivity)
      return false;
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
      return 0.5;
    case ScaleModel::BFACTOR:
      // A relative B-factor
      return 2.0;
    case ScaleModel::SECONDARY:
      // A secondary beam parameter
      return 0.05;
    case ScaleModel::TILE:
      {
        // A detector beam parameter, but which one?
        std::pair<int,int> idxpar = DetectorParameterNumber(Ipar);
        return detector_scales[idxpar.first].LargeShift(idxpar.second);
      }
    case ScaleModel::WAVELENGTH:
      return 0.5;
    default:
      return 0.0;
    }
    return false;
  }
  //--------------------------------------------------------------
  //! store parameter variance information
  void ScaleModel::storeParameterVariances
  (const TNT::Fortran_Matrix<floatType>& H,
   const double& wd2in,
   const int& nminusm)
  {
    ASSERT (H.dim(1) == H.dim(2));
    int Np = H.dim(1);  // dimension = number of parameters
    ASSERT (Np == nparameters);
    VC.resize(Np,Np);
    varpar.resize(Np);
    wd2 = wd2in;
    nfreedom = nminusm;

    double scale = wd2/double(nfreedom);
    for (int i=0;i<Np;++i) {
      varpar[i] = scale * H(i+1,i+1);  // H indexed from 1
      for (int j=0;j<Np;++j) {
        VC(i,j) = scale * H(i+1,j+1);
      }
    }
  }
  //--------------------------------------------------------------
  // return false if number of variance parameters is not same as nparameters
  // OK (true) if no variance used
  bool ScaleModel::checkVarianceNumbers() const
  {
    if (parametersdusage == scala::ScaleSpecification::NONE) {
      return true;
    }
    if (nparameters != int(varpar.size())) {
      return false;
    }
    return true;
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
    // save tile parameters
    ds += "Ndetscales "+ clipper::String(ndetscales)+"\n";
    if (ndetscales > 0) {
      for (int i=0;i<ndetscales;++i) {
        ds += detector_scales[i].FormatSave();
      }
      ASSERT (int(detector_scale_index_run.size()) == nruns);
      ds += "Detector_scale_index_run\n"+
        StringUtil::FormatSaveVector(detector_scale_index_run);
    }

    // sds
    ds += "sd_rotation "+ clipper::String(sd_rotation)+"\n";
    ds += "sd_bfactor "+ clipper::String(sd_bfactor)+"\n";
    ds += "sd_zerob "+ clipper::String(sd_zerob)+"\n";
    ds += "sd_surface "+ clipper::String(sd_surface)+"\n";
    ds += "tie_tile_number "+clipper::String(int(tie_tile.size()))+"\n";
    for (size_t i=0;i<tie_tile.size();++i) {
      ds += "tie_tile "+clipper::String(tie_tile[i])+"\n";
    }

    // Normalisation
    ds += "Scalenormrun "+ clipper::String(scalenormrun)+"\n";
    ds += "Scalenormbatch "+ clipper::String(scalenormbatch)+"\n";
    ds += "Bfacnormrun "+ clipper::String(bfacnormrun)+"\n";
    ds += "Bfacnormbatch "+ clipper::String(bfacnormbatch)+"\n";

    // SD estimates
    // Variance/covariance information
    ds += "Variances{\n";
    ds += "Nparameters "+clipper::String(nparameters)+"\n";
    ds += "nfreedom "+clipper::String(nfreedom)+"\n";
    ds += "wD2 "+StringUtil::ftos(wd2)+"\n";
    ds += "ParameterVariance\n"+StringUtil::FormatSaveVector(varpar)+"\n";
    ds += "VarianceCovariance {\n";
    ds += StringUtil::FormatSaveArray(VC)+"\n";
    ds += "}\n}\n";
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

    RunsFromSavefile savefileruns(FR, runlist);
    if (savefileruns.NumberRunsFound() < int(runlist.size())) {
      clipper::Message::message(Message_fatal
                        ("RESTORE not all runs found in save file"));
    }

    int svnruns = savefileruns.NumberRunsInSaveFile();
    std::vector<int> runsfromsavefile = savefileruns.Runsfromsavefile();
    // runsfromsavefile now contains the index in runlist
    // for each run in save file, or -1 if not wanted

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
    if (nssc != nsecscales) {
      ReportErrors::printWarning("RESTORE incompatible secondary scale models",
                                 "RestoreChangeSecScaleModel",false);
    }
    nsecscales = nssc;
    secondary_scales.resize(nsecscales);
    sec_scale_index_run.clear();
    if (nssc > 0) {
      for (int i = 0;i<nssc;++i) { // loop secondary scales in file
        secondary_scales[i].Restore(FR);
      } // end loop secondary scales
      FR.ReadTag("Sec_scale_index_run");
      // Index for each run
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
    }

    FR.ReadTag("Ndetscales");
    int ndsc = FR.Int();
    if (ndsc != ndetscales) {
        clipper::Message::message(Message_fatal
          ("RESTORE incompatible detector models"));
    }
    if (ndetscales > 0) {
      for (int i=0;i<ndetscales;++i) {
        detector_scales[i].Restore(FR);
      }
      FR.ReadTag("Detector_scale_index_run");
      ASSERT (int(detector_scale_index_run.size()) == nruns);
      detector_scale_index_run = FR.IntVec(nruns);
    }

    // sds
    FR.ReadTag("sd_rotation"); sd_rotation = FR.Double();
    FR.ReadTag("sd_bfactor"); sd_bfactor = FR.Double();
    FR.ReadTag("sd_zerob"); sd_zerob = FR.Double();
    FR.ReadTag("sd_surface"); sd_surface = FR.Double();
    int sdtn;
    FR.ReadTag("tie_tile_number"); sdtn = FR.Int();
    tie_tile.resize(sdtn);
    for (int i=0;i<sdtn;++i) {
      FR.ReadTag("tie_tile"); tie_tile[i] = FR.Double();
    }
    // Normalisation
    FR.ReadTag("Scalenormrun"); scalenormrun = FR.Int();
    FR.ReadTag("Scalenormbatch"); scalenormbatch = FR.Int();
    FR.ReadTag("Bfacnormrun"); bfacnormrun = FR.Int();
    FR.ReadTag("Bfacnormbatch"); bfacnormbatch = FR.Int();

    CountParameters(); // set parameter counts etc

    // Read (optional) variances for nparameters
    if (FR.GetTag() != "Variances") { // variances present
      FR.ReadTag("Nparameters"); int npar = FR.Int();
      if (npar != nparameters) {
        clipper::Message::message(Message_fatal
          ("RESTORE Variances: wrong number of parameters "+
           clipper::String(npar)+", "+clipper::String(nparameters)));
      }
      FR.ReadTag("nfreedom"); nfreedom = FR.Int();
      FR.ReadTag("wD2"); wd2 = FR.Double();
      FR.ReadTag("ParameterVariance"); varpar = FR.DoubleVec(nparameters);
      FR.ReadTag("VarianceCovariance"); FR.Skip();
      VC = FR.Array2d(nparameters,nparameters);
      if (!FR.CheckEnd()) {
        clipper::Message::message(Message_warn
                                  ("Restore error missing '}'"));
      }
    }

    scalesin.close();
  }
  //--------------------------------------------------------------
  void ScaleModel::WriteImage(const std::string fname,
                              phaser_io::Output& output) const
  //! Write image[s] for each detector scale
  {
    std::string imagefilename = fname;
    if (getenv(imagefilename.c_str()) != NULL) { // it's an environment variable
      imagefilename = std::string(getenv(phaser_io::stoup(imagefilename).c_str()));
    }

    for (int idsc=0;idsc<ndetscales;++idsc) {
      std::string basename = FileNameNoExtension(imagefilename);
      std::string ext = FileNameExtension(imagefilename);
      if (ext == "") {ext = "img";}
      ext = "."+ext;
      std::string name = imagefilename;
      if (ndetscales > 1) {
        name = basename+"_"+StringUtil::Strip(StringUtil::itos(idsc+1,4))+ext;
      } else {
        name = basename+ext;
      }
      detector_scales[idsc].WriteImage(name);
      output.logTab(0,LOGFILE,
                    "\nDetector scale image (x1000) written to file "+name);
    }
  }
  //--------------------------------------------------------------
  // Get secondary factor for ksecscale'th secscale, theta, phi (radians)
  double ScaleModel::secscale(const size_t& ksecscale,
                              const double& theta, const double& phi) const
  {
    return secondary_scales[ksecscale].Scale(theta, phi);
  }
  //--------------------------------------------------------------
  std::vector<std::string> ScaleModel::secondaryscaletypes() const
  // for each secondary scale (if any),
  //  return "SECONDARY" or "ABSORPTION"+pole
  {
    std::vector<std::string> sectypes;
    if (nsecscales == 0) {return sectypes;}

    for (size_t k=0; k<secondary_scales.size(); k++) {
      if (secondary_scales[k].Type() == SecondaryScale::NONE) {
        sectypes.push_back("");
      } else if (secondary_scales[k].Type() == SecondaryScale::SECONDARY) {
        sectypes.push_back("SECONDARY");
      } else if (secondary_scales[k].Type() == SecondaryScale::ABSORPTION) {
        sectypes.push_back("ABSORPTION " + secondary_scales[k].formatPole());
      }
    }
    return sectypes;
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
