// aimless.cpp

// [Scale &] analyse unmerged diffraction data
// (Scala replacement)
//
// Phil Evans, MRC Laboratory of Molecular Biology
//             Hills Road, Cambridge CB2 0QH
// 2006-?
//

#include "aimless.hh"
#include "scalemodel.hh"
#include "applyscales.hh"
#include "initialscales.hh"
#include "scalerefine.hh"
#include "scalerefinefh.hh"
#include "reject.hh"
#include "statistics.hh"
#include "summarystatistics.hh"
#include "icering.hh"
#include "analyseanom.hh"
#include "selectscalingreflections.hh"
#include "sdmodel.hh"
#include "selectsdcorrreflections.hh"
#include "writeoutputfiles.hh"
#include "timer.hh"
#include "string_util.hh"
#include "optimisecombine.hh"
#include "file_util.hh"
#include "observationstatuscontrol.hh"
#include "analyseoverlaps.hh"
#include "referencelist.hh"
#include "rejectbatches.hh"
#include "runcorrelations.hh"
#include "report_errors.hh"
#include "secondaryscalestats.hh"
#include "comparetoreference.hh"

#ifdef _MSC_VER
#include <io.h>
#define STDIN_FILENO 0
#endif

#if _OPENMP
#include <omp.h>
#endif

#include "version.hh"
#include "ccp4/ccp4_program.h"

using namespace scala;
using phaser_io::LOGFILE;
using phaser_io::LXML;

//--------------------------------------------------------------
bool IsOnline()
// Returns true if stdin is not connected to a file
{
  if (isatty(STDIN_FILENO)== 0)
    return false;
  else
    return true;
}
//--------------------------------------------------------------
int main(int argc, char* argv[])
{
  Timer overalltime;

  // Initialise output object, CCP4 mode, write header
  phaser_io::Output output;
  output.setPackageCCP4();
  output.SetMaxLineWidth(600);
  //  output.openOutputStreams("DEBUG");
  //  output.openOutputStreams("VERBOSE");
  //  output.setVerbose(true, true);

  ReportErrors reportErrors(output);

  phaser_io::InterpretCommandLine CL(argc, argv, output);
  if (!CL.Run()) {
    return 0;
  }

  // Initialise CCP4 command line parser
  CCP4::ccp4fyp(argc, argv);
  CCP4::ccp4ProgramName (PROGRAM_NAME.c_str());
  std::string rcsdate = "$Date: "+std::string(PROGRAM_DATE2)+"$";
  CCP4::ccp4RCSDate     (rcsdate.c_str());
  CCP4::ccp4_prog_vers(PROGRAM_VERSION.c_str());
  CCP4::ccp4_banner();

  CL.printCommandLine(output);  // echo command line arguments

  GlobalControls GC;
  std::string hklin_filename = "";
  std::string hklref_filename = "";
  std::string xyzref_filename = "";

  Timer timer;
  int errorstatus = 0;

  try {
    // Input from command line: optional HKLIN filename
    //
    ///    phaser_io::InterpretCommandLine CL(argc, argv, output);
    hklin_filename = CL.getHKLIN1();
    hklref_filename = CL.getHKLREF();
    xyzref_filename = CL.getXYZIN();

    bool logHeaderdone = false;
    if (CL.getXMLOUT() != "") {
      output.setXmlout(CL.getXMLOUT());
      // Print header and optionly initialise XML output
      output.logHeader(LOGFILE);
      logHeaderdone = true;
    }

    // Read input & store unless "-n" or "--no-input" switches given
    phaser_io::InputAll input(CL.noInput(), output);
    ////    phaser_io::InputAll input(IsOnline(), output);
    input.Analyse();

    // Optional HKLIN
    if (input.getHKLIN() != "") {
      hklin_filename = input.getHKLIN();
    }

    // Optional HKLREF
    if (input.getHKLREF() != "") {
      hklref_filename = input.getHKLREF();
    }
    if (input.getXYZIN()  != "") {
      xyzref_filename = input.getXYZIN();
    }
    if (hklref_filename != "" && xyzref_filename != "") {
      ReportErrors::printFatalError
        ("Cannot have both HKLREF and XYZIN filenames given");
    }

    // XMLOUT command
    if (input.getXMLOUT() != "")
      {output.setXmlout(input.getXMLOUT());}

    // ROGUES command
    std::string roguesfilename = "ROGUES";
    if (input.getROGUES() != "")
      {roguesfilename = input.getROGUES();}

    if (!logHeaderdone) {
      output.logHeader(LOGFILE);
    }
    PrintTitle(output);
    if (hklin_filename == "")
      ReportErrors::printFatalError("HKLIN filename not given");

   // TITLE command, defaults to title from HKLIN file (see below)
    std::string runTitle = input.Title();

    // Set Bitflag control, to reject any flagged observations
    ObservationFlagControl ObsFlagControlRejectall;

    std::string SpaceGroup = "";

    // *************************************************************
    //  Setup up controls
    //    these are set from command input (or default)

    // Control of flow through program
    FlowControl FC;

    // Scala control classes
    //  run controls
    //  partials controls
    //  outlier controls
    all_controls controls;

    // Scaling refinement options: parallel stuff set later
    controls.refinecontrol = input.RefineControl();

    if (controls.refinecontrol.Ncyc1() <= 0) FC.SetRoughScale(false);
    if (controls.refinecontrol.Ncycles() <= 0) FC.SetMainScale(false);
    if (input.Onlymerge()) FC.SetOnlyMerge();  // No scaling option ONLYMERGE

    if (input.InitialUnity()) {
      FC.SetInitialScale(false);  // INITIAL UNITY, no initial scaling
    }
    // Output settings
    //  Merge/unmerge options
    OutputControls outputcontrols = input.Outputcontrols();

    // set filenames from command line or environment
    std::string hklout_filename = CL.getHKLOUT();
    if (input.getHKLOUT() != "") {
      hklout_filename = input.getHKLOUT();
    }
    std::string hkloutunmerged_filename = CL.getHKLOUTUNMERGED();
    if (input.getUNMERGEDOUT() != "") {
      hkloutunmerged_filename = input.getUNMERGEDOUT();
    }

    outputcontrols.SetFilenames(hklout_filename, hkloutunmerged_filename,
                                CL.getSCAOUT(), CL.getSCAOUTUNMERGED());

    //  Setup up controls for reflection & column selection etc
    // Set Profile-fitted [default] or integrated intensity
    col_controls column_selection;
    // store column selection from input (keyword INTENSITIES)
    // for later storage into static class SelectI
    // this may be overridden if no IPR column is present in the file
    column_selection.SetIcolFlag(input.GetIcolFlag(), input.GetCombineImid(),
                                 input.GetCombinePower());

    // File selection flags (resolution, datasets, batches etc)
    file_select file_sel(input,0);

    // Explicit run definition, if present
    controls.runs.StoreRunBatchSelection(input.Runsetselection());
    // Store input resolution ranges by run, if any set
    controls.runs.StoreResoByRun(input.GetResoByRun());

    // ANOMALOUS ON|OFF
    controls.anomalouscontrol.Anomalous = input.getANOMALOUS();
    // True if "Anomalous" command given
    controls.anomalouscontrol.FlagInput = input.AnomalousFlagInput();

    controls.partials = partial_controls(input.getFracLimMin(),
                                         input.getFracLimMax(),
                                         input.getSclMinLim(),
                                         input.getCheck(),
                                         input.getMaxGap());

    controls.analysis = AnalysisControls(input.getResoBins(),
                                         input.getIntBins(),
                                         input.ConeAngle(),
                                         input.MinimumHalfdatasetCC(),
                                         input.MinimumHalfdatasetAnomCC(),
                                         input.MinimumIoverSigma(),
                                         input.MinimumBatchIoverSigma(),
                                         input.SmoothStatisticsRange(),
                                         input.BatchGroupRange(),
                                         input.DetectorAnalysis());

    // OutlierControl from input or defaults
    // Set number of datasets later
    if (FC.OnlyMerge() && !input.isrejectset()) {
      // Onlymerge with no explicit setting, turn off outlier tests
      controls.outlierScale.SetNoreject();
      controls.outlierMerge.SetNoreject();
    } else {
      controls.outlierScale = input.GetOutlierControlsScale();
      controls.outlierMerge = input.GetOutlierControlsMerge();
    }

    // Controls on acceptable observation flags, for merging only
    // All flagged observations will be omitted from scaling
    controls.observationflagcontrol = input.Observationflagcontrol();

    // Make list of required columns in unmerged HKLIN file
    MtzIO::column_labels column_list = MtzIO::setup_columns();

    controls.plotcontrol.xmgraceoutput = input.XMGRoutput();

    //<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
    // Construct reflection list by reading MTZ file
    hkl_unmerge_list hkl_list;
    output.logTab(0,LOGFILE,
                   "\n---------------------------------------------------------------\n");
    output.logTab(0,LOGFILE, "\nReading data from HKLIN filename: " + hklin_filename + "\n");
    const int verbose = +3;
    MtzIO::MtzUnmrgFile mtzin;
    int fileSeries = 0;
    PxdName InputPxdName = input.getNAME();
    double Tolerance = 2.0;
    std::string outputstring;
    timer.Start();
    Scell inputcell = input.getCELL();
    mtzin.AddHklList(fileSeries,
                     hklin_filename, file_sel, column_selection, column_list,
                     controls, InputPxdName, inputcell,
                     Tolerance, outputstring, verbose,
                     hkl_list);
    output.logTab(0,LOGFILE,outputstring);
    output.logTab(0,LOGFILE,
                  "\nTime for reading HKLIN: "+timer.format(true));

    hkl_list.ResetObsAccept(ObsFlagControlRejectall);

    // Do we need to change symmetry?
    if (SpaceGroup != "" || GC.IsReindexSet()) {
      output.logTab(0,LOGFILE, "\nReindexing or changing symmetry\n");
      hkl_symmetry new_symm(SpaceGroup);
      if (SpaceGroup == "") {
        new_symm = hkl_list.symmetry();
      } else {
        output.logTab(0,LOGFILE,"  New space group: "+SpaceGroup);
      }
      output.logTab(0,LOGFILE, "  Reindex operator from input: "
                    + GC.Reindex().as_hkl() + "\n");

      // reindex, sort & reorganise hkl list
      hkl_list.change_symmetry(new_symm, GC.Reindex());
    }
    // sort & organise reflection list as required
    int Nrefl = hkl_list.prepare();
    Nrefl = Nrefl;
    // Sum partials, return number of partials
    int Npart = hkl_list.sum_partials();
    int Nobs = hkl_list.num_observations();
    Nobs = Nobs;
    std::vector<std::string> dummycolumnlabels;  // dummy for unmerged file
    PrintUnmergedHeaderStuff(hkl_list, output, verbose);
    PrintFileInfoToXML("HKLIN",hklin_filename,
                       hkl_list.cell(),
                       hkl_list.symmetry().symbol_xHM(),
                       dummycolumnlabels,
                       output);

    bool multilattice = hkl_list.MultiLattice();
    // Set to use singletons only for scaling etc, if multiple lattices present
    SetOverlapFlags(true, hkl_list);

    bool onlyUseSingletons = false;
    if (hkl_list.NumberofLattices() != hkl_list.NumberofMainLattices()){
      // Not all (multi)lattices present in file, set global flag to use singletons only
      onlyUseSingletons = true;
    }

    // Read optional reference file for statistics, check for compatibility
    ReferenceList hklreflist;
    bool refOK;
    double toleranceratio = 1.0;
    // use SF calculation with bulk solvent, do it later so that
    //  it can be scaled
    bool SF_BULK_SOLVENT = true;
    bool referencedata = (hklref_filename != "" || xyzref_filename != "");
    bool readRefFirst = referencedata;
    if (referencedata) {
      // reference data given
      // true to read reference list here, else later
      readRefFirst = true;
      // read/generate reference first if HKLREF
      if (xyzref_filename != "" && SF_BULK_SOLVENT) {
        if (controls.refinecontrol.Reference()) {
          // refining against Fcalc^2 doesn't work well
          ReportErrors::printFatalError
            ("Cannot refine against coordinate reference XYZIN");
        } else {
          readRefFirst = false;
        }
      }
    }

    if (readRefFirst) {  // reference from HKLREF or XYZIN with no bulk solvent
      bool verbose = true;
      if (hklref_filename != "") {
        std::string refmessage = "\nReference file for analysis (HKLREF)";
        if (controls.refinecontrol.Reference()) {
          refmessage = "\nReference file for scaling and analysis (HKLREF)";
        }
        output.logTab(0,LOGFILE, refmessage);
        hklreflist.init(hklref_filename,
                        input.getLABREF_I(), input.getLABREF_sigI(),
                        hkl_list.ResRange().ResHigh(), verbose,
                        output);
        PrintFileInfoToXML("HKLREF",hklref_filename,
                           hklreflist.Cell(),
                           hklreflist.SpaceGroupSymbol(),
                           hklreflist.columnLabels(),
                           output);
      } else {
        // xyzref (XYZIN) coordinates given
        hklreflist.init(xyzref_filename,
                        hkl_list.ResRange().ResHigh(), verbose,
                        output);
        PrintFileInfoToXML("XYZIN",xyzref_filename,
                           hklreflist.Cell(),
                           hklreflist.SpaceGroupSymbol(),
                           dummycolumnlabels,
                           output);
      }
      refOK = hklreflist.checkCompatible(hkl_list, toleranceratio);
      if (!refOK) {
        std::string s = "HKLREF file is incompatible with HKLIN file\n";
        s += hklreflist.formatError();
        output.logTab(0,LXML, StringUtil::MakeXMLtag("FatalErrorMessage",s));
        Message::message(Message_fatal(s));
      }
    }

    // Set number of datasets for anomalous outliers
    controls.outlierScale.SetNdatasets(0);  // no anomalous rejections in scaling
    controls.outlierMerge.SetNdatasets(hkl_list.num_datasets());

    bool optimiseCombine = true;
    if (column_selection.IcolFlag() > 0) {
      // INTENSITIES COMBINE option
      if (column_selection.IsImidSet()) {
        // Imid set explicitly, so it won't be changed
        output.logTab(0,LOGFILE, "\n"+SelectI::format());
        optimiseCombine = false;
      } else {
        output.logTab(0,LOGFILE, "\nSelection of intensity type (Isum or Ipr) will be optimised");
        output.logTab(0,LOGFILE, "Profile fitted value Ipr will be used for 1st scaling");
        SelectI::SetIcolFlag(-1, -1.0);
      }
    } else { // INTENSITIES PROFILE or INTEGRATED
      output.logTab(0,LOGFILE, "\n"+SelectI::format());
      optimiseCombine = false;
    }
    if (Npart > 0) {
      PrintPartialCounts(hkl_list, controls, output);
    }

    // Set up SD correction model for all runs, fulls & partials for each run
    // from input or by default
    // If Onlymerge and SDcorrection not set explciitly, set null correction
    bool setnull = (FC.OnlyMerge() && !input.isSDcorrectionSet());
    SDmodel SD_model = CreateSDmodel(input, hkl_list.RunList(), setnull);
    if (SD_model.SampleSD()) {
      SelectedObservations::SetSampleSD(SD_model.MinimumSample());
    }
    FC.sdoptimise = true;
    if (FC.OnlyMerge() && !input.SDC_RefineSet()) {
      FC.sdoptimise = false;  // normally optimise SD correction unless onlymerge
    }
    // Print outlier information
    PrintOutlierSettings(controls, output);

    // Default title
    if (runTitle.size() == 0) runTitle = hkl_list.Title();

    // Set up scale model
    scala::ScaleModel AllScales;
    SecondaryScaleStats secondaryscalestats;
    bool initialscale = FC.initialScale;
    double overallmeankI = -1000000.0;  // mean I

    #if _OPENMP
    // Parallel stuff
    output.logTab(0,LOGFILE,
                  "\nParallisation of refinement:\n");
    if (controls.refinecontrol.Nprocs() < 0) {
      // NPROC AUTO, set number of processors from number of observations Nobs
      // maximum number of processors to use
      const int MAXUSEDPROCS = 8;
      // number of observations/processor: what is the "best" value?
      const double NUMOBSPERPROC = 200000;
      // Number to use
      int nproc = Max(1,Nint(double(Nobs)/NUMOBSPERPROC));
      nproc = Min(MAXUSEDPROCS, nproc);
      controls.refinecontrol.SetNprocs(double(nproc));
      output.logTab(0,LOGFILE,
                    std::string("Number of processors determined automatically\n")+
                    "  from number of observations "+
                    StringUtil::Strip(clipper::String(Nobs))+
                    " and number/processor "+
                    StringUtil::Strip(clipper::String(NUMOBSPERPROC)));
    }
    output.logTab(0,LOGFILE,
                  controls.refinecontrol.format());

    omp_set_num_threads(controls.refinecontrol.Nprocs());
#endif

    Normalise NormRes;
    // Put Ice ring object into hkl_list
    Rings icerings;
    int ringlisttype = input.ringListType();
    bool iceringreject = input.iceRingReject();

    icerings.setIceRings(ringlisttype);
    if (ringlisttype == 0) {
      // Keep icerings for everything
      icerings.SetRejectAll(false);
    } else {
      // Reject icerings for scaling etc
      icerings.SetRejectAll(true);
    }
    hkl_list.SetIceRings(icerings);

    Rings normicerings;  // for Normalisation
    normicerings.setIceRings(ringlisttype);
    normicerings.SetRejectAll(false);

    ApplyScales applyscales;

    // Restoring scales from file?
    FC.restore = input.Restore();
    if (FC.restore) {
      AllScales.init(input, hkl_list, controls, output);
      AllScales.Restore(input.RestoreFileName(),
                        hkl_list.RunList());
      initialscale = false; // no initial scales
      if (!AllScales.IsRefinable()) {
        FC.SetOnlyMerge();
        AllScales.SetConstant(hkl_list, output);
        std::string s = "No scaling done, "+AllScales.whyNotRefineable();
        output.logTab(0,LOGFILE, s);
        output.logTab(0,LXML,
                      StringUtil::MakeXMLwithclass("ScaleModelFail",
                                                   s, false,
                                                   "warningmessage"));
      }
      AllScales.PrintLayout(output);
      AllScales.PrintScales(output);

      secondaryscalestats.init(AllScales);
      applyscales.scale(AllScales, hkl_list, onlyUseSingletons);
      secondaryscalestats.getScaleStats(hkl_list, onlyUseSingletons);
      secondaryscalestats.PrintSecondaryCorrections(output);
      overallmeankI = applyscales.meanI();

      // Restore SD correction
      if (input.SDC_NumberInput() != 0) {
        FC.sdcorrectionsinput = true;
      }
      if (!input.SDC_RefineSet() && !FC.sdcorrectionsinput) {
        SD_model.Restore(input.RestoreFileName(),
                         hkl_list.RunList());
      }
      if (FC.OnlyMerge()) {
        if (!input.SDC_RefineSet()) {
          // no sdoptimisation if restore and onlymerge and SDCORR REFINE not set
          FC.sdoptimise = false;
        }
        // ----- Optimise Combine settings
        if (optimiseCombine) {
          OptimiseCombine OptCombine(hkl_list, output);
          if (OptCombine.IsOptimised()) {
            output.logTab(0,LOGFILE,
                          "\nTime for optimisation of intensity type selection: "+
                          timer.format(true));
            // Revaluate summed partials for scaling
            hkl_list.sum_partials(true);
          }
        }
      }
    } else {
      // not RESTORE
      if (FC.OnlyMerge()) {
        // Onlymerge, set scales CONSTANT
        AllScales.SetConstant(hkl_list, output);
        if(!input.SDC_RefineSet()) {
          // turn off sd optimisation unless requested
          FC.sdoptimise = false;
        }
      } else {
        // Set up scale model, from input commands & reflection list
        AllScales.init(input, hkl_list, controls, output);
        // If no refinable parameters, set OnlyMerge
        if (!AllScales.IsRefinable()) {
	  if (!controls.refinecontrol.Reference()) {
	    // but not refine reference
	    FC.SetOnlyMerge();
	    AllScales.SetConstant(hkl_list, output);
	    // turn off sd optimisation unless explicit
	    if (!input.SDC_RefineSet()) {
	      FC.sdoptimise = false;
	    }
	    std::string s = "No scaling done, "+AllScales.whyNotRefineable();
	    output.logTab(0,LOGFILE, s);
	    output.logTab(0,LXML,
			  StringUtil::MakeXMLwithclass("ScaleModelFail",
                                                     s, false,
                                                     "warningmessage"));
	  }
        } else {
          AllScales.PrintLayout(output);
        }
      }
    }

    // control SD analyses:
    //  = 0  two analyses expected, after rough scaling and main scaling
    //  = +1    same, after first analysis
    //  = -1 only one analysis expected (onlymerge) ie first & last
    int firstSDanalysis = 0;
    if (FC.OnlyMerge()) {
      firstSDanalysis = -1;
    }

    // If space group is non-chiral, turn off anomalous, unless explicit
    bool chiralsg = hkl_list.symmetry().IsChiral();
    bool lowmultiplicity = false;
    if (!chiralsg) {
      if (controls.anomalouscontrol.FlagInput &&
          controls.anomalouscontrol.Anomalous) {
        // Explicit anomalous ON for non-chiral space group, print warning
        ReportErrors::printWarning
          ("Explicit anomalous ON for non-chiral space group",
           "NonChiralSpacegroupWarning");
      } else {
        // switch off Anomalous
        controls.anomalouscontrol.Anomalous = false;
        controls.anomalouscontrol.AnomalousSDcorr = false;
      }
    } else {
      // By default do SD correction optimisation only within I+/I- sets
      // in case there is anomalous, unless multiplicity is low
      controls.anomalouscontrol.AnomalousSDcorr = true;

      lowmultiplicity = false;
      double multiplicity = double(hkl_list.num_observations())/
        double(hkl_list.num_reflections_valid());
      if (controls.anomalouscontrol.AnomalousSDcorr) {
        // If multiplicity low, combine I+ & I- for SD correction
        // * tried this but didn't always work on bad data
        // Try again with lower threshold
        const double MINMULTFORSDCORR = 1.5;
        // Separate I+ & I- for SD correction, unless multiplicity is low
        if (multiplicity < MINMULTFORSDCORR) {
          output.logTab(0,LOGFILE,
                        std::string("WARNING: multiplicity low, ")+
                        StringUtil::Strip(StringUtil::ftos(multiplicity,8,1))+
                        " (below threshold "+
                        StringUtil::Strip(StringUtil::ftos(MINMULTFORSDCORR,8,1))+
                        "), so combine I+ and I- for SD correction\n");
          controls.anomalouscontrol.AnomalousSDcorr = false;
          lowmultiplicity = true;
        }
      }
      output.logFlush();
    }

    // Check valid reference scaling
    bool scaletoreference = false;
    if (!FC.OnlyMerge()) {  // no check if only merge
      if (controls.refinecontrol.Reference()) {
        // Must have an appropriate reference set
        if (hklreflist.num_obs() <= 0) {
          // No reference set read, why not?
          std::string refmessage = "Cannot refine against reference data:";
          if (!referencedata) {
            refmessage += " no reference file assigned";
          } else if (!readRefFirst && xyzref_filename != "") {
            refmessage += " coordinates XYZIN not allowed as reference";
          }
          ReportErrors::printFatalError(refmessage);
        }
        FC.initialScale = false; //  no initial scale for now
        initialscale = false; // no initial scales
        scaletoreference = true;
        hklreflist.recordScaleReference(output); // write info to log and XML
      }
    }

    bool suppressScaling = false;  // maybe suppress scaling (onlymerge)
    double minimum_overlap = input.Minimum_overlap();
    int allowed_gap = input.Maximum_gap();

    // ----- Initial scales
    bool allowgap = false;
    if (initialscale) {
      timer.Start();
      InitialScales initialscales(hkl_list, AllScales, initialscale, controls, output);
      output.logTab(0,LOGFILE,
		    "\nTime for initial scaling: "+timer.format(true));
      output.logFlush();

      // test for enough data for scaling
      if (! initialscales.enoughData(minimum_overlap, allowed_gap)) {
	//  there is a gap in a run, should scaling be suppressed?
	if (controls.runs.Explicit()) {
	  if (initialscales.numberofrotationranges() == 1) {
	    // Just one range, suppress scaling
	    suppressScaling = true;
	  } else {
	    // OK if there are explicit runs
	    //   unless minimum_overlap is set explicitly to > 0
	    std::string s = "There is a gap in a run, with overlap < "+
	      StringUtil::ftos(std::abs(minimum_overlap));
	    ReportErrors::printWarning(s, "OverlapWarning");
	    allowgap = true;
	    if (minimum_overlap > 0.0) {
	      suppressScaling = true;
	    }
	  }
	} else {
	  suppressScaling = true;
	}
      }

      if (minimum_overlap == 0.0 && !allowgap) {
	output.logTab(0,LOGFILE,
		      "\nNo test for minimum fractional overlap between rotation ranges (INITIAL MINIMUM_OVERLAP)");
      }
      initialscales.reportOverlapXML(output, allowgap);

      // Option to reject batches based on extreme scale factors
      // relevant for eg XFEL data
      if (controls.outlierScale.Reject(ALL).batchrejectfactor > 0.0) {
        RejectBatches rejectBatches(hkl_list, AllScales, controls, output);
      }
    } else {
      if (input.InitialUnity()) {
        output.logTab(0,LOGFILE,
              "\n========= Initial scales all set to 1.0 =========\n");
      }
    }

    if (suppressScaling) {
      // set onlymerge
      FC.SetOnlyMerge();
      output.logTab(0,LOGFILE,
         "\n**** NB No scaling will be done as there seems to be insufficient data\n");
      output.logTabPrintf(0,LOGFILE,
                "        Minimum threshold for fractional overlap between rotation ranges= %7.2f\n",
                          std::abs(minimum_overlap));
      AllScales.SetConstant(hkl_list, output);
      SD_model.SetRefine(0);  // SDCORRECTION NOREFINE
    } else {
      if (!FC.restore || !FC.OnlyMerge()) {
        if (allowgap) {
          output.logTabPrintf(0,LOGFILE,
               "\nThere is a gap in rotation ranges are above the minimum threshold for fractional overlap\n between rotation ranges = %5.2f, but this is allowed with explicit run definition\n",
                              std::abs(minimum_overlap));
        } else {
        output.logTabPrintf(0,LOGFILE,
                            "\nSufficient rotation ranges are above the minimum threshold for fractional overlap between rotation ranges = %5.2f\n",
                            std::abs(minimum_overlap));
        }
      }
    }

    if (FC.OnlyMerge()) {
      output.logTab(0,LXML,
                    "<OnlyMerge/>");
    }

    // Set weighting for SD model
    //   (doesn't make a huge difference at least in some tests)
    SD_model.SetWeight(input.SDCweightType());
    //    SD_model.SetVarianceWeights();
    //    SD_model.SetSqrtScaleWeights();
    //    SD_model.SetUnitWeights();


    bool anomOn = false;  // no anomalous for scaling
    // ----- first rough scaling
    WriteRogues DummyRogues;
    if (FC.roughScale) {
      timer.Start();
      double IovSDmin = controls.refinecontrol.IovSDmin();
      double E2min = -1.0;   // no |E^2| selection here
      double E2max = -1.0;
      Normalise NormResDummy;    // dummy, not used
      std::vector<int> selrej =
        SelectScalingReflections(hkl_list, SD_model, AllScales, IovSDmin,
                                 NormResDummy, E2min, E2max);
      output.logTabPrintf(0,LOGFILE,
                          "\n========= First round scaling =========\n");
      output.logTabPrintf(0,LOGFILE,
          "\nFirst scaling: %7d reflections selected from %8d with I/sd > %6.2f",
			  selrej[0], selrej[1], IovSDmin);
      if (selrej[2] > 1) {
        output.logTabPrintf(0,LOGFILE,
                            ", using every %3d'th reflection above that limit",
                          selrej[2]);
      }
      output.logTabPrintf(0,LOGFILE,"\n");

      // For the 1st round, force any tile corrections to be radially symmetric
      bool nparchanged = AllScales.symmetricTiles(true);
      // and switch off secondary scaling
      AllScales.switchSecondaryScales(false);

      if (controls.refinecontrol.BFGS() ||
          controls.refinecontrol.Reference()) {
        if (nparchanged && AllScales.haveParameterVariances()) {
          // don't use incorrect parameter variances (from restore) if number has changed
          AllScales.ignoreParameterVariances();
        }
        ScaleRefine(hkl_list, hklreflist, AllScales, SD_model, controls,
                    controls.refinecontrol.Ncyc1(), false, output);
        if (nparchanged) {
          // don't use incorrect parameter variances if number has changed
          AllScales.ignoreParameterVariances();
        }
      } else {
        ScaleRefineFH(hkl_list, AllScales, controls,
                      controls.refinecontrol.Ncycles(), output);
      }

      // Apply all scales (ie store g for each observation, the original I is unchanged)
      // All observations are scaled, including rejected ones
      //AllScales.PrintScales(output); //^^
      applyscales.scale(AllScales, hkl_list, onlyUseSingletons);
      overallmeankI = applyscales.meanI();

      output.logFlush();

      // Overall Normalisation
      double MinIsigRatio = 0.6;  // resolution cutoff for Emaxtest
      // Clear all outlier & other status flags (except ObsFlags)
      ClearObsStatus(hkl_list);
      NormRes.init(hkl_list, MinIsigRatio, normicerings, 0);

      // -- 1st outlier rejection
      // Use outlier flags appropriate for scaling
      std::vector<int> nrejs;
      if (controls.outlierScale.GetOutlierPolicy() != OutlierControl::NOREJECT) {
        RejectOutlier(hkl_list, SD_model, NormRes, anomOn,
                      controls.outlierScale, DummyRogues);
        nrejs = CountOutliers(hkl_list);
        output.logTabPrintf(0,LOGFILE,
                            "\nNumber of outliers within I+ || I- sets: %6d,  between I+ & I- %6d, on |E|max %6d\n",
                            nrejs[0], nrejs[1], nrejs[2]);
      } else {
        output.logTab(0,LOGFILE,"\nNo outlier rejection in scaling");
      }
      output.logFlush();
      // -- End 1st outlier rejection

      // For the 2nd round, allow tile corrections to vary azimuthally
      AllScales.symmetricTiles(false);
      // and switch on secondary scaling
      AllScales.switchSecondaryScales(true);

      hkl_list.ResetObsAccept(ObsFlagControlRejectall);  // count observation flag rejects
      output.logTab(0,LOGFILE,
                    "\nTime for 1st scaling: "+timer.format(true));

      // ----- Optimise Combine settings
      if (optimiseCombine) {
        OptimiseCombine OptCombine(hkl_list, output);
        if (OptCombine.IsOptimised()) {
          output.logTab(0,LOGFILE,
                        "\nTime for optimisation of intensity type selection: "+
                        timer.format(true));
          // Revaluate summed partials for scaling
          hkl_list.sum_partials(true);
        }
      }
      // -----

      // First rough SD analysis before main scaling: assumes scales have been applied
      // SD corrections are not applied, but SD_model is updated
      //  hkl_list is const
      if (SD_model.Refine()) {
        // NormRes just used for intensity binning
        AnalyseSD(SD_model, hkl_list, controls, NormRes,
                  firstSDanalysis, output);
        firstSDanalysis = +1;
        output.logFlush();
      }

      // -- 2nd outlier rejection
      // Use outlier flags appropriate for scaling
      if (controls.outlierScale.GetOutlierPolicy() != OutlierControl::NOREJECT) {
        RejectOutlier(hkl_list, SD_model, NormRes, anomOn,
                      controls.outlierScale, DummyRogues);
        nrejs = CountOutliers(hkl_list);
        output.logTabPrintf(0,LOGFILE,
                            "\nNumber of outliers within I+ || I- sets: %6d,  between I+ & I- %6d, on |E|max %6d\n",
                            nrejs[0], nrejs[1], nrejs[2]);
      } else {
        output.logTab(0,LOGFILE,"\nNo outlier rejection in scaling");
      }
      output.logFlush();
      // -- End 2nd outlier rejection

    }   // end 1st scaling

    // ----- Main scaling
    if (FC.mainScale) {
      timer.Start();
      double IovSDmin = 0.0;
      double E2min = controls.refinecontrol.E2min();
      double E2max = controls.refinecontrol.E2max();
      std::vector<int> selrej =
        SelectScalingReflections(hkl_list, SD_model, AllScales, IovSDmin,
                                 NormRes, E2min, E2max);
      output.logTabPrintf(0,LOGFILE,
                          "\n========= Main scaling =========\n");
      output.logTabPrintf(0,LOGFILE,
         "\nMain scaling: %7d reflections selected from %8d with |E^2| > %6.2f and |E^2| < %6.2f\n\n",
			  selrej[0], selrej[1], E2min, E2max);
      int Ncyc = controls.refinecontrol.Ncycles();
      if (controls.refinecontrol.BFGS() ||
          controls.refinecontrol.Reference()) {
        ScaleRefine(hkl_list, hklreflist,
                    AllScales, SD_model, controls, Ncyc, true, output);
      } else {
        ScaleRefineFH(hkl_list, AllScales, controls, Ncyc, output);
      }
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // Apply all scales
      // All observations are scaled, including rejected ones
      hkl_list.ResetReflAccept();  // set to accept everything
      secondaryscalestats.init(AllScales);
      applyscales.scale(AllScales, hkl_list, onlyUseSingletons);
      output.logTab(0,LOGFILE,
                    "\nTime for main scaling: "+timer.format(true));

      secondaryscalestats.getScaleStats(hkl_list, onlyUseSingletons);
      secondaryscalestats.PrintSecondaryCorrections(output);
      AllScales.WriteImage("TILEIMAGE", output);

      output.logFlush();
    }

    // Resolution ranges
    ResoRange ResRange = hkl_list.ResLimRange();
    // set number of bins to override default if required
    int nresbin = controls.analysis.NresoBins();
    if (nresbin > 0) {
      ResRange.SetNbins(nresbin);
    } else {
      nresbin =  ResRange.Nbins();
    }
    double resrangewidth = ResRange.Width(); // bin width for maximum resolution

    // Overall Normalisation
    double MinIsigRatio = 0.6;  // resolution cutoff for Emaxtest
    int printlevel = 0;  // 1 to dump to norm.plot
    NormRes.init(hkl_list, MinIsigRatio, normicerings, printlevel);

    hkl_list.ResetObsAccept(ObsFlagControlRejectall);  // count observation flag rejects

    if (FC.sdoptimise) {
      // Clear all outlier & other status flags (except ObsFlags)
      //  temporary outlier rejection is done in AnalyseSD
      //    (actually in SumsforSDcorrection in refinesdcorrection.cpp)
      ClearObsStatus(hkl_list);
      // SD analysis: assumes scales have been applied
      // SD corrections are not applied, but SD_model is updated
      //  hkl_list is const
      // NormRes just used for intensity binning
      // firstSDanalysis = 0  two analyses expected, after rough scaling and main scaling
      //                 = +1    same, after first analysis
      //                 = -1 only one analysis expected (onlymerge) ie first & last
      AnalyseSD(SD_model, hkl_list, controls, NormRes, firstSDanalysis, output);
      output.logFlush();
    } else {
      // No optimisation
      if (FC.restore) {
        if (FC.sdcorrectionsinput) {
          output.logTab(0,LOGFILE,
                        "\nSD correction parameters input\n"+
                        SD_model.format());
        } else {
          output.logTab(0,LOGFILE,
                "\nSD correction parameters restored from SCALES file\n"+
                      SD_model.format());
        }
        output.logTab(0,LXML,SD_model.asXML());
      } else {
        output.logTab(0,LOGFILE,
        "\nSD correction parameters\n"+SD_model.format());
      }
      AnalyseNormalProbability(SD_model, hkl_list, controls, true, output);
    }

    // Always dump scale model and SDmodel
    WriteToFile(input.DumpFileName(),
                AllScales.FormatSave(hkl_list.RunList())+SD_model.FormatSave());

    firstSDanalysis = +2;

    AllSummaryStatistics allsummarystatistics;

    //   slopes of anomalous normal probability plots for each dataset
    //   Plot ANOMPLOT normal probability plot
    //   Inflate anomalous rejection criterion according to analysis on DelAnom (in controls)
    //  hkl_list is const
    anomOn = controls.anomalouscontrol.Anomalous;  // from input
    ResoRange resrangeanom = ResRange;
    // For statistics, reset range to go from same "infinite" resolution,
    //  to the maximum resolution over all datasets
    double lowres = 10000.;
    resrangeanom.SetRange(lowres, ResRange.ResHigh());
    resrangeanom.SetNbins(nresbin);
    AnalyseAnom analysanom(hkl_list, SD_model, controls, resrangeanom, true, output);
    std::vector<double> anomProbSlopes = analysanom.Slopes();
    output.logFlush();

    // Analyse distribution anomalous differences to get estimate of
    //  maximum likely values for final statistics
    //  hkl_list is const
    AllAnomDistributions allAnomDistributions(hkl_list, SD_model, controls,
                                              analysanom, resrangeanom);
    allAnomDistributions.SetSlope(anomProbSlopes);
    if (hkl_list.num_accepted_datasets() > 1) {
      allAnomDistributions.Print(output);
    }

    // Do we really have anomalous?
    //  -1 no data, 0 weak anomalous, +1 significant
    int isanom =
      allAnomDistributions.IsAnomalous(controls); // true if anomalous
    bool anomfound = (isanom > 0);
    if (isanom < 0) {
      // no information
      allsummarystatistics.SetAnomStatus(AnomalousStatus::NO_ANOMALOUS_DATA);
    } else {
      // Should we change the options?
      if (controls.anomalouscontrol.FlagInput || !chiralsg) {
	// explicit anomalous on or off from input
	if (controls.anomalouscontrol.Anomalous) { // On
	  if (anomfound) {
	    allsummarystatistics.SetAnomStatus(AnomalousStatus::ANOMALOUS_ON_FOUND);
	  } else {
	    allsummarystatistics.SetAnomStatus(AnomalousStatus::ANOMALOUS_ON_ABSENT);
	  }
	} else { // Off
	  if (anomfound) {
	    allsummarystatistics.SetAnomStatus(AnomalousStatus::ANOMALOUS_OFF_FOUND);
	  } else {
	    allsummarystatistics.SetAnomStatus(AnomalousStatus::ANOMALOUS_OFF_ABSENT);
	  }
	}
      } else { // No explicit flag given, set appropriately
	if (anomfound) {
	  controls.anomalouscontrol.Anomalous = true;
	  controls.anomalouscontrol.AnomalousSDcorr = true;
	  if (lowmultiplicity) {controls.anomalouscontrol.AnomalousSDcorr = false;}
	  allsummarystatistics.SetAnomStatus(AnomalousStatus::ANOMALOUS_FOUND);
	} else {
	  controls.anomalouscontrol.Anomalous = false;
	  controls.anomalouscontrol.AnomalousSDcorr = false;
	  allsummarystatistics.SetAnomStatus(AnomalousStatus::ANOMALOUS_ABSENT);
	}
      }
    }
    output.logTab(0,LOGFILE,"\n"+
                  AnomDistribution::formatStatus(allsummarystatistics.AnomStatus())+
                  "\n");

    output.logTab(0,LOGFILE,"\nOutlier analysis\n================\n");

    // Reset flag controls to possibly accept some observations flagged as rejected,
    // eg profile-fitted overloads etc
    hkl_list.ResetObsAccept(controls.observationflagcontrol);
    // Check for outliers & reject them
    // Start rogues output, ROGUES file & ROGUEPLOT
    //  for plotting ice rings, use shortest wavelength all datasets
    double wavelength = 100000000.;
    for (int i=0;i<hkl_list.num_datasets();++i) {
      wavelength = Min(wavelength, hkl_list.dataset(i).wavelength());
    }

    // ----   Get anisotropy to use for normaliastion and Emax test,
    //        assume same for all datasets
    // and Normalisation 
    Timer anisotime;
    AnisotropicAnalysis anisoanal;
    int idts = -1;  // all together
    // Get principal axes of anisotropy depending on symmetry and data
    anisoanal.init(hkl_list, idts, SD_model);
    anisoanal.SetConeAngle(controls.analysis.ConeAngle());
    output.logTab(0, LOGFILE,
         "\nTime for determination of anisotropic axes: "+anisotime.format(true));

    // Update normalisation
    bool aniso_normalisation = input.anisotropicNormalisation();
    if (anisoanal.IsCubic()) {
      aniso_normalisation = false;
    }
    std::string normmsg = "";
    if (aniso_normalisation) {
      clipper::U_aniso_frac uanisofrac = anisoanal.U_aniso_frac();
      NormRes.setAniso(uanisofrac);
      output.logTab(0, LOGFILE,
		    "Normalisation with anisotropy correction");
      normmsg = "True";
    } else {
      NormRes.noAniso();
      output.logTab(0, LOGFILE,
		    "Normalisation NO anisotropy correction");
      normmsg = "False";
    }
    output.logTab(0,LXML,
	     StringUtil::MakeXMLtag("anisotropicNormalisation",normmsg));
    // Overall Normalisation
    //printlevel = 1;  // 1 to dump to norm.plot
    printlevel = 0;  // 1 to dump to norm.plot
    NormRes.init(hkl_list, MinIsigRatio, normicerings, printlevel,true);

    // ----

    if (controls.outlierMerge.GetOutlierPolicy() != OutlierControl::NOREJECT) {
      // doRoguePlot true as long as we have geometric data for all batches
      // to calculate detector position
      bool doRoguePlot = hkl_list.validOrientation();  // false if no orientation
      if (doRoguePlot && FC.OnlyMerge()) {
        // we need secondary beam directions
        int pole = 0;  // whatever
        hkl_list.CalcSecondaryBeams(pole);
      }

      Rings ploticerings;
      int rltype = ringlisttype;
      if (ringlisttype == 0) {
	rltype = 2;
      }
      ploticerings.setIceRings(rltype);

      WriteRogues RoguesList(roguesfilename, true, doRoguePlot, multilattice,
                             runTitle, hkl_list.Srange().max(), wavelength,
			     ploticerings,
                             controls.outlierMerge,
                             controls.plotcontrol.xmgraceoutput);
      //  hkl_list is updated for status, but SDs are not changed
      RejectOutlier(hkl_list, SD_model, NormRes,
                    controls.anomalouscontrol.Anomalous,
                    controls.outlierMerge, RoguesList);
      RoguesList.End();
      std::vector<int> nrejs = CountOutliers(hkl_list);

      if (controls.outlierMerge.isEmaxTest()) {
        // yes there is an Emax test
        output.logTab(0,LOGFILE, "Test for Emax ");
      } else {
        output.logTab(0,LOGFILE, "No Emax test");
      }
      output.logTabPrintf(0,LOGFILE,
                          "Number of rejected outliers within I+ || I- sets: %6d,  between I+ & I- %6d, on |E|max %6d\n",
                          nrejs[0], nrejs[1], nrejs[2]);
      output.logTab(0,LXML,CountOutliersXML(nrejs));
      output.logFlush();
      if (doRoguePlot) {
        // ROGUEPLOT to XML
        output.logTab(0,LXML, RoguesList.formatXML());
      }
    } else {
      output.logTab(0,LOGFILE,"\nNo outlier rejection in merging");
    }

    output.logTab(0,LXML, controls.outlierMerge.formatXML());

    // for each dataset
    std::vector<Analyseoverlaps>  analyseoverlaps(hkl_list.num_datasets());
    if (multilattice && !onlyUseSingletons) {
      // Check all overlapped observations for self overlaps, ie overlaps with the same
      // or symmetry-related spots (but not Friedel-related).
      // Now the data are scaled, these overlap set can be combined into a pseudo-singleton
      // Also accumulate overlap statistics etc
      SetOverlapFlags(false, hkl_list); // include overlaps
      //      bool verbose = true;
      bool verbose = false;
      for (int idts=0;idts<hkl_list.num_datasets();++idts) {
        // Project/Crystal/Dataset for this dataset
        PxdName dataset_pxd = hkl_list.dataset(idts).pxdname();
        output.logTab(0,LOGFILE,"\nCheck for self-overlaps, dataset "+
                      dataset_pxd.dname());
        analyseoverlaps[idts].init(hkl_list, idts, verbose, output);
        int nmerged = analyseoverlaps[idts].NumberMerged();
        if (nmerged > 0) {
            output.logTabPrintf(0,LOGFILE,
     "  Number of self-overlapped observations reclassified as singletons = %5d\n",
                            nmerged);
        }
      } // end loop datasets
      SetOverlapFlags(true, hkl_list); // exclude overlaps
    }

    output.logTab(0,LOGFILE,
                  "\n********************\n* Final statistics *\n********************\n");
    output.logTab(0,LOGFILE,controls.observationflagcontrol.PrintCounts());
    output.logTab(0,LXML,controls.observationflagcontrol.asXML());

    // Scale optional reference file for statistics to observed data
    if (referencedata) {
      int datasetindex = -2;  // combine all datasets together
      if (xyzref_filename != "" && SF_BULK_SOLVENT) {
        // calculate SF from atoms with bulk solvent
        refOK =
          hklreflist.SFcalcScaleToObserved(xyzref_filename,
                                           hkl_list, datasetindex,
                                           SD_model, toleranceratio,
                                           true, output);
        ASSERT (refOK); // checked earlier
        PrintFileInfoToXML("XYZIN",xyzref_filename,
                           hklreflist.Cell(),
                           hklreflist.SpaceGroupSymbol(),
                           dummycolumnlabels,
                           output);
      }
      // (re)scale to observed, wherever the F list has come from
      refOK =
        hklreflist.scaleToObserved(hkl_list, datasetindex,
                                   SD_model, toleranceratio, output);
        ASSERT (refOK); // checked earlier

	// Comparison of multiple datasets to reference
	if (!hklreflist.IsEmpty() && (hkl_list.num_datasets() > 1)) {
	  CompareToReference comparetoreference(hkl_list, hklreflist, ResRange);
	  comparetoreference.printTable(output);
	}
    }

    // get maximum batch width
    std::vector<Batch> batches = hkl_list.Batches();

    // Gather & print all statistics for each dataset ------------------------------------
    for (int idts=0;idts<hkl_list.num_datasets();++idts) {
      if (hkl_list.dataset(idts).accepted()) {
        //  hkl_list is const
        // NormRes just used for intensity binning
        // Store summary statistics for this dataset
        // Resolution range for this dataset
        ResoRange resrangedataset = hkl_list.dataset(idts).ResRange();

        // For statistics, reset range to go from same "infinite" resolution
        double lowres = 10000.;
        resrangedataset.SetRange(lowres, resrangedataset.ResHigh());
        // Use same resolution bin width for all datasets
        resrangedataset.SetWidth(resrangewidth);
        // and check the number of bins is not > nresbin
        if (resrangedataset.Nbins() > nresbin) {
          resrangedataset.SetNbins(nresbin);
        }

        AnomDistribution anomds = allAnomDistributions.Anomdistribution(idts);
        double aslope = anomProbSlopes.at(idts);

	// Put icering definition into hkl_list copying reject status from input
	// Default no rejections
 	Rings finalicerings = hkl_list.getIceRings();
	finalicerings.SetRejectAll(iceringreject);
	finalicerings.CopyRejRings(finalicerings);
	hkl_list.SetIceRings(finalicerings);

        SummaryStatistics sumstat = Statistics(AllScales, hkl_list, SD_model,
                                               controls, idts, resrangedataset,
                                               NormRes, anisoanal, anomds, aslope,
                                               hklreflist,
                                               output);

        allsummarystatistics.AddSummaryStatistics(sumstat);

        if (multilattice && !onlyUseSingletons) {
          analyseoverlaps[idts].PrintOverlapTable(output);
        }

        bool Result = true;
        if (hkl_list.num_datasets() != 1) {
          output.logTab(0,LOGFILE,
                        "==============================================================\n");
          Result = false;
        } else {
          applyscales.print(batches, output);
        }

        // Print summary as a Results table if one dataset, otherwise just to logfile
        allsummarystatistics.PrintOneSummaryTable(idts, Result, output);
        output.logFlush();
      }
    } // end loop datasets -----------------------------------------

    // run-run correlations: don't do them if there are too many
    const int MAXRUNCORRELATION = 20;
    if ((hkl_list.num_runs() > 1) &&
        (hkl_list.num_runs() < MAXRUNCORRELATION)) {
      RunCorrelations runcorrelations(hkl_list, SD_model, NormRes, nresbin);
      runcorrelations.formatTable(output);
    }

    if (hkl_list.num_datasets() > 1) { // summary for multiple datasets
      applyscales.print(batches, output);
      bool Result = true;
      allsummarystatistics.PrintSummaryTable(Result,
                                     controls.anomalouscontrol.Anomalous, output);
    }
    output.logTab(0,LOGFILE,
     "\n==============================================================\n");
    output.WriteResult();
    output.logTab(0,LOGFILE,
     "==============================================================\n");

    // Output merged data
    MergedList mergedlist;
    if (outputcontrols.Merged()) {
      // Put all data into clipper classes
      mergedlist.init(hkl_list, SD_model, runTitle);
      // Output merged reflections file(s)
      scala::WriteMergedOutputFiles(mergedlist,
                                    outputcontrols, output);
    }

    if (multilattice && !onlyUseSingletons) {
      SetOverlapFlags(false, hkl_list); // include overlaps
    }
    if (outputcontrols.UnMerged()) {
      // Output unmerged reflections file(s), including overlaps
      scala::WriteUnmergedOutputFiles(runTitle, hkl_list, SD_model,
                                      NormRes.Imax(), outputcontrols, output);
    }
    // Record output options to XML
    outputcontrols.recordToXML(output);


  }  // end try


  catch (phaser_io::PreprocessorError& capErr) {
    output.logWarning(LOGFILE, capErr.partialEcho()+"PREPROCESSOR ERROR: " + capErr.Message()
                      + "\n");
    errorstatus = 1;
  }
  catch (phaser_io::SyntaxError& ccp4Err) {
    output.logWarning(LOGFILE, ccp4Err.Echo()+"\nSYNTAX ERROR: " + ccp4Err.Message()
                      + "\n");
    errorstatus = 2;
  }
  catch (phaser_io::InputError& inpErr) {
    output.logWarning(LOGFILE, inpErr.Echo()+"\nINPUT ERROR: " + inpErr.Message()
                      + "\n");
    errorstatus = 3;
  }

  catch (Message_fatal& message) {
    output.logWarning(LOGFILE, "\nFATAL ERROR message: \n"
                      + message.text() + "\n");
    errorstatus = 4;
  }

  catch (std::bad_alloc const& err) {
    output.logWarning(LOGFILE,
                      std::string("\nERROR: ")
                      + std::string(err.what())+"\n"+
                      " You have run out of memory to store the data\n"+
                      "  you may need more memory or a 64-bit machine.\n");
    errorstatus = 5;
  }

  catch (std::exception const& err) {
    output.logWarning(LOGFILE, "\nUNHANDLED EXCEPTION: " + std::string(err.what())+"\n");
    errorstatus = 6;
  }
  catch (...) {
    output.logWarning(LOGFILE, "\nUNKNOWN EXCEPTION TYPE\n");
    errorstatus = 7;
  }

  if (output.doXmlout()) {
    output.logTab(0, LXML,"</AIMLESS>");
    output.unsetXmlout();
  }

  Citation citation
    ("P.R.Evans and G.N.Murshudov, 'How good are my data and what is the resolution?'"+
                    std::string(" Acta Cryst. D69, 1204-1214  (2013)."),
     "http://journals.iucr.org/d/issues/2013/07/00/ba5190/index.html");
  output.logTab(0,LOGFILE, "\n"+citation.MakeLogCitation());


  output.logTab(0, LOGFILE,
                "\nEnd of aimless job, total time: "+overalltime.format(true)+"\n\n");

  return errorstatus;
}
