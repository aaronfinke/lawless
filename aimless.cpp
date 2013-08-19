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

  // Initialise CCP4 command line parser
  CCP4::ccp4fyp(argc, argv);

  CCP4::ccp4ProgramName (PROGRAM_NAME.c_str());
  std::string rcsdate = "$Date: "+std::string(PROGRAM_DATE2)+"$";
  CCP4::ccp4RCSDate     (rcsdate.c_str());
  CCP4::ccp4_prog_vers(PROGRAM_VERSION.c_str());
  CCP4::ccp4_banner();

  // Initialise output object, CCP4 mode, write header
  phaser_io::Output output;
  output.setPackageCCP4();
  output.SetMaxLineWidth(600);
  //  output.openOutputStreams("DEBUG");
  //  output.openOutputStreams("VERBOSE");
  //  output.setVerbose(true, true);

  GlobalControls GC;
  std::string hklin_filename = "";

  Timer timer;

  try {
    // Input from command line: optional HKLIN filename
    // 
    phaser_io::InterpretCommandLine CL(argc, argv);
    hklin_filename = CL.getHKLIN1();
    if (CL.getXMLOUT() != "")
      {output.setXmlout(CL.getXMLOUT());}

    // Read input & store if not "online"
    phaser_io::InputAll input(IsOnline(), output);
    input.Analyse();

    // XMLOUT command
    if (input.getXMLOUT() != "")
      {output.setXmlout(input.getXMLOUT());}

    output.logHeader(LOGFILE);
    PrintTitle(output);
    if (hklin_filename == "")
      Message::message(Message_fatal("HKLIN filename not given"));


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
    outputcontrols.SetFilenames(CL.getHKLOUT(), CL.getHKLOUTUNMERGED(),
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
    controls.runs.StoreRunBatchSelection(input.RunBatches());
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
    int NbatchsmoothDefault = -1;
    controls.analysis = AnalysisControls(input.getResoBins(),
					 input.getIntBins(),
					 input.ConeAngle(),
					 input.MinimumHalfdatasetCC(),
					 input.MinimumIoverSigma(),
					 input.MinimumBatchIoverSigma(),
					 NbatchsmoothDefault,
					 input.DetectorAnalysis());

    // OutlierControl from input or defaults
    // Set number of datasets later
    controls.outlierScale = input.GetOutlierControlsScale();
    controls.outlierMerge = input.GetOutlierControlsMerge();

    // Controls on acceptable observation flags, for merging only
    // All flagged observations will be omitted from scaling
    controls.observationflagcontrol = input.Observationflagcontrol();

    // Make list of required columns in unmerged HKLIN file
    MtzIO::column_labels column_list = MtzIO::setup_columns();

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
    mtzin.AddHklList(fileSeries,
		     hklin_filename, file_sel, column_selection, column_list,
		     controls, InputPxdName, scala::Scell(),
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
    PrintUnmergedHeaderStuff(hkl_list, output, verbose);
    PrintFileInfoToXML("HKLIN",hklin_filename,
		       hkl_list.cell(),
		       hkl_list.symmetry().symbol_xHM(),
		       output);

    bool multilattice = hkl_list.MultiLattice();
    // Set to use singletons only for scaling etc, if multiple lattices present
    SetOverlapFlags(true, hkl_list);

    bool onlyUseSingletons = false;
    if (hkl_list.NumberofLattices() != hkl_list.NumberofMainLattices()){
      // Not all (multi)lattices present in file, set gloabl flag to use singletons only
      onlyUseSingletons = true;
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
    SDmodel SD_model = CreateSDmodel(input, hkl_list.RunList());
    FC.sdoptimise = true;  // normally optimise SD correction unless onlymerge && restore 

    // Print outlier information
    PrintOutlierSettings(controls, output);

    // Default title
    if (runTitle.size() == 0) runTitle = hkl_list.Title();

    // Set up scale model
    scala::ScaleModel AllScales;
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
      const float NUMOBSPERPROC = 200000;
      // Number to use
      int nproc = Max(1,Nint(float(Nobs)/NUMOBSPERPROC));
      nproc = Min(MAXUSEDPROCS, nproc);
      controls.refinecontrol.SetNprocs(float(nproc));
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

    // Restoring scales from file?
    FC.restore = input.Restore();
    if (FC.restore) {
      AllScales.init(input, hkl_list, output);
      AllScales.Restore(input.RestoreFileName(),
			hkl_list.RunList());
      initialscale = false; // no initial scales
      AllScales.PrintLayout(output);
      AllScales.PrintScales(output);
      overallmeankI = ApplyScales(AllScales, hkl_list, onlyUseSingletons);

      // Restore SD correction
      if (!input.SDC_RefineSet()) {
	SD_model.Restore(input.RestoreFileName(),
			 hkl_list.RunList());
      }
      if (FC.OnlyMerge()) {
	if (!input.SDC_RefineSet()) {
	  // no sdoptimisation if restore and onlymerge and SDCORR REFINE not set
	  FC.sdoptimise = false;
	}
      }
    } else {
      if (FC.OnlyMerge()) {
	// Onlymerge, set scales CONSTANT
	AllScales.SetConstant(hkl_list, output);
      } else {
	// Set up scale model, from input commands & reflection list
	AllScales.init(input, hkl_list, output);
	AllScales.PrintLayout(output);
	// If no refinable parameters, set OnlyMerge
	if (!AllScales.IsRefinable()) {
	  FC.SetOnlyMerge();
	  output.logTab(0,LOGFILE,
			"No refinable parameters");
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
  
    // By default do SD correction optimisation only within I+/I- sets
    // in case there is anomalous, unless multiplicity is low
    controls.anomalouscontrol.AnomalousSDcorr = true;

    bool lowmultiplicity = false;
    if (controls.anomalouscontrol.AnomalousSDcorr) {
      // If multiplicity low, combine I+ & I- for SD correction
      // * tried this but didn't always work on bad data
      // Try again with lower threshold
      float multiplicity = float(hkl_list.num_observations())/
	float(hkl_list.num_reflections_valid());
      const float MINMULTFORSDCORR = 1.5;
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

    // ----- Initial scales
    if (initialscale) {
      timer.Start();
      InitialScales(hkl_list, AllScales, controls, output);
      output.logTab(0,LOGFILE,
		    "\nTime for initial scaling: "+timer.format(true));
      output.logFlush();
    } else {
      if (input.InitialUnity()) {
	output.logTab(0,LOGFILE,
	      "\n========= Initial scales all set to 1.0 =========\n");
      }
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
      float IovSDmin = controls.refinecontrol.IovSDmin();
      float E2min = -1.0;   // no |E^2| selection here
      float E2max = -1.0;
      std::pair<int,int> selrej = 
	SelectScalingReflections(hkl_list, SD_model, AllScales, IovSDmin, E2min, E2max);
      output.logTabPrintf(0,LOGFILE,
			  "\n========= First round scaling =========\n");
      output.logTabPrintf(0,LOGFILE,
	  "\nFirst scaling: %7d reflections selected from %8d with I/sd > %6.2f",
			  hkl_list.num_reflections()-selrej.first,
			  hkl_list.num_reflections(), IovSDmin);
      if (selrej.second > 1) {
	output.logTabPrintf(0,LOGFILE,
			    ", using every %3d'th reflection above that limit",
			  selrej.second);
      }
      output.logTabPrintf(0,LOGFILE,"\n");

      // For the 1st round, force any tile corrections to be radially symmetric
      AllScales.symmetricTiles(true);

      if (controls.refinecontrol.BFGS()) {
	ScaleRefine(hkl_list, AllScales, controls,
		    controls.refinecontrol.Ncyc1(), false, output);
      } else {
	ScaleRefineFH(hkl_list, AllScales, controls,
		      controls.refinecontrol.Ncycles(), output);
      }
      // Apply all scales (ie store g for each observation, the original I is unchanged)
      // All observations are scaled, including rejected ones
      overallmeankI = ApplyScales(AllScales, hkl_list, onlyUseSingletons);

      output.logFlush();

      // Overall Normalisation 
      double MinIsigRatio = -1.0;  // no resolution cutoff
      bool Overall = true;  // no run|time variation, just one curve
      Rings NoRings;        // no omission of ice rings
      // Set up resolution bins for normalisation, allowing for number of observations
      ResoRange ResRangeN(hkl_list.RRange().min(), hkl_list.RRange().max(),
			  hkl_list.num_observations());  // shouldn't be changed in Normalise
      Normalise NormRes = SetNormalise(hkl_list, MinIsigRatio, Overall,
				       ResRangeN, NoRings, 0);

      // -- 1st outlier rejection
      // Use outlier flags appropriate for scaling
      RejectOutlier(hkl_list, SD_model, NormRes, anomOn,
		    controls.outlierScale, DummyRogues);
      std::vector<int> nrejs = CountOutliers(hkl_list);
      output.logTabPrintf(0,LOGFILE,
       	  "\nNumber of outliers within I+ || I- sets: %6d,  between I+ & I- %6d, on |E|max %6d\n",
			  nrejs[0], nrejs[1], nrejs[2]);
      output.logFlush();
      // -- End 1st outlier rejection

      // For the 2nd round, allow tile corrections to vary azimuthally
      AllScales.symmetricTiles(false);
      
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
      RejectOutlier(hkl_list, SD_model, NormRes, anomOn,
		    controls.outlierScale, DummyRogues);
      nrejs = CountOutliers(hkl_list);
      output.logTabPrintf(0,LOGFILE,
       	  "\nNumber of outliers within I+ || I- sets: %6d,  between I+ & I- %6d, on |E|max %6d\n",
			  nrejs[0], nrejs[1], nrejs[2]);
      output.logFlush();
      // -- End 2nd outlier rejection

    }   // end 1st scaling

    // ----- Main scaling
    if (FC.mainScale) {
      timer.Start();
      float IovSDmin = 0.0;
      float E2min = controls.refinecontrol.E2min();
      float E2max = controls.refinecontrol.E2max();
      std::pair<int,int> selrej = 
	SelectScalingReflections(hkl_list, SD_model, AllScales, IovSDmin, E2min, E2max);
      output.logTabPrintf(0,LOGFILE,
			  "\n========= Main scaling =========\n");
      output.logTabPrintf(0,LOGFILE,
	 "\nMain scaling: %7d reflections selected from %8d with |E^2| > %6.2f and |E^2| < %6.2f\n\n",
			  hkl_list.num_reflections()-selrej.first,
			  hkl_list.num_reflections(), E2min, E2max);
      int Ncyc = controls.refinecontrol.Ncycles();
      if (controls.refinecontrol.BFGS()) {
	ScaleRefine(hkl_list, AllScales, controls, Ncyc, true, output);
      } else {
	ScaleRefineFH(hkl_list, AllScales, controls, Ncyc, output);
      }
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // Apply all scales
      // All observations are scaled, including rejected ones
      hkl_list.ResetReflAccept();  // set to accept everything
      ApplyScales(AllScales, hkl_list, onlyUseSingletons);
      output.logTab(0,LOGFILE,
		    "\nTime for main scaling: "+timer.format(true));

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

    // Overall Normalisation 
    double MinIsigRatio = -1.0;  // no resolution cutoff
    bool Overall = true;  // no run|time variation, just one curve
    Rings NoRings;        // no omission of ice rings
    // Set up resolution bins for normalisation, allowing for number of observations
    ResoRange ResRangeN(hkl_list.RRange().min(), hkl_list.RRange().max(),
			  hkl_list.num_observations());  // shouldn't be changed in Normalise
    Normalise NormRes = SetNormalise(hkl_list, MinIsigRatio, Overall,
				     ResRangeN, NoRings, 0);

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
	output.logTab(0,LOGFILE,
      "\nSD correction parameters restored from SCALES file\n"+
		      SD_model.format());
	output.logTab(0,LXML,SD_model.asXML());
      } else {
	output.logTab(0,LOGFILE,
	"\nSD correction parameters\n"+SD_model.format());
      }
      AnalyseNormalProbability(SD_model, hkl_list, controls, true, output);
    }

    // if scaling done, dump scale model and SDmodel
    if (FC.mainScale) {
      WriteToFile(input.DumpFileName(),
		  AllScales.FormatSave(hkl_list.RunList())+SD_model.FormatSave());
    }

    firstSDanalysis = +2;

    AllSummaryStatistics allsummarystatistics;

    //   slopes of anomalous normal probability plots for each dataset
    //   Plot ANOMPLOT normal probability plot
    //   Inflate anomalous rejection criterion according to analysis on DelAnom (in controls)
    //  hkl_list is const
    anomOn = controls.anomalouscontrol.Anomalous;  // from input
    ResoRange resrangeanom = ResRange;
    // For statistics, reset range to go from same "infinite" resolution
    float lowres = 10000.;
    resrangeanom.SetRange(lowres, ResRange.ResHigh());
    resrangeanom.SetNbins(nresbin);
    AnalyseAnom analysanom(hkl_list, SD_model, controls, resrangeanom, true, output);
    std::vector<float> anomProbSlopes = analysanom.Slopes();
    output.logFlush();

    // Analyse distribution anomalous differences to get estimate of
    //  maximum likely values for final statistics
    //  hkl_list is const
    AllAnomDistributions allAnomDistributions(hkl_list, SD_model, controls,
					      analysanom,
    					      resrangeanom, NormRes);
    allAnomDistributions.SetSlope(anomProbSlopes);
    allAnomDistributions.Print(output);

    // Do we really have anomalous?
    bool anomfound =
      allAnomDistributions.IsAnomalous(controls); // true if anomalous
    // Should we change the options?
    if (input.AnomalousFlagInput()) {
      // explicit anomalous on or off from input
      if (controls.anomalouscontrol.Anomalous) { // On
	if (anomfound) { 
	  allsummarystatistics.SetAnomStatus(AnomDistribution::ANOMALOUS_ON_FOUND);
	} else {
	  allsummarystatistics.SetAnomStatus(AnomDistribution::ANOMALOUS_ON_ABSENT);
	}
      } else { // Off
	if (anomfound) { 
	  allsummarystatistics.SetAnomStatus(AnomDistribution::ANOMALOUS_OFF_FOUND);
	} else {
	  allsummarystatistics.SetAnomStatus(AnomDistribution::ANOMALOUS_OFF_ABSENT);
	}
      }
    } else { // No explicit flag given, set appropriately
      if (anomfound) { 
	controls.anomalouscontrol.Anomalous = true;
	controls.anomalouscontrol.AnomalousSDcorr = true;
	if (lowmultiplicity) {controls.anomalouscontrol.AnomalousSDcorr = false;}
	allsummarystatistics.SetAnomStatus(AnomDistribution::ANOMALOUS_FOUND);
      } else {
	controls.anomalouscontrol.Anomalous = false;
	controls.anomalouscontrol.AnomalousSDcorr = false;
	allsummarystatistics.SetAnomStatus(AnomDistribution::ANOMALOUS_ABSENT);
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
    float wavelength = 100000000.;
    for (int i=0;i<hkl_list.num_datasets();++i) {
      wavelength = Min(wavelength, hkl_list.dataset(i).wavelength());
    }
    // doRouguePlot true as long as we have geometric data for all batches
    // to calculate detector position
    bool doRouguePlot = hkl_list.validOrientation();  // false if no orientation
    WriteRogues RoguesList(true, doRouguePlot, multilattice,
			   runTitle, hkl_list.Srange().max(), wavelength,
			   controls.outlierMerge);
    //  hkl_list is updated for status, but SDs are not changed
    RejectOutlier(hkl_list, SD_model, NormRes, controls.anomalouscontrol.Anomalous,
		  controls.outlierMerge, RoguesList);
    RoguesList.End();
    std::vector<int> nrejs = CountOutliers(hkl_list);
    output.logTabPrintf(0,LOGFILE,
	"Number of rejected outliers within I+ || I- sets: %6d,  between I+ & I- %6d, on |E|max %6d\n",
			nrejs[0], nrejs[1], nrejs[2]);
    output.logTab(0,LXML,CountOutliersXML(nrejs));
    output.logFlush();
    // ROGUEPLOT to XML
    output.logTab(0,LXML, RoguesList.formatXML());


    // for each dataset
    std::vector<Analyseoverlaps>  analyseoverlaps(hkl_list.num_datasets());
    if (multilattice && !onlyUseSingletons) {
      // Check all overlapped observations for self overlaps, ie overlaps with the same
      // or symmetry-related spots (but not Friedel-related).
      // Now the data are scaled, these overlap set can be combined into a pseudo-singleton
      // Also accumulate overlap statistics etc
      SetOverlapFlags(false, hkl_list); // include overlaps
      bool verbose = true;
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

    // Smoothing of batch statistics
    double smoothwidth = input.SmoothStatisticsRange(); // angular range, -1 if unset
    if (smoothwidth <= 0.0) {
      // Set from scale rotation range if known
      double spacing = AllScales.primary_scale(0).Spacing();
      if (spacing > 0.0) {
	smoothwidth = spacing;  // = spacing
      }
    }
    // get maximum batch width
    std::vector<Batch> batches = hkl_list.Batches();
    float width = 0.0;
    for (size_t ib=0;ib<batches.size();++ib) {
      width = Max(width, batches[ib].PhiRange());
    }

    int nbatchsmooth = 5;
    if (width > 0.0 && smoothwidth > 0.0) {
      nbatchsmooth = (Nint(smoothwidth/width)/2)*2 + 1; // force to be odd
    }
    controls.analysis.SetNbatchSmooth(nbatchsmooth);
    double resrangewidth;

    // Gather & print all statistics
    for (int idts=0;idts<hkl_list.num_datasets();++idts) {
      //  hkl_list is const
      // NormRes just used for intensity binning
      // Store summary statistics for this dataset
      // Resolution range for this dataset
      ResoRange resrangedataset = hkl_list.dataset(idts).ResRange();

      // For statistics, reset range to go from same "infinite" resolution
      float lowres = 10000.;
      resrangedataset.SetRange(lowres, resrangedataset.ResHigh());
      if (idts == 0) { // 1st dataset
	resrangedataset.SetNbins(nresbin);
	resrangewidth = resrangedataset.Width();
      } else if (idts > 0) {
	// Use same resolution bin width for all datasets
	resrangedataset.SetWidth(resrangewidth);
      }

      AnomDistribution anomds = allAnomDistributions.Anomdistribution(idts);
      float aslope = anomProbSlopes[idts];
      SummaryStatistics sumstat = Statistics(AllScales, hkl_list, SD_model, controls, idts,
					     resrangedataset, NormRes, anomds, aslope, output);
      allsummarystatistics.AddSummaryStatistics(sumstat);

      if (multilattice && !onlyUseSingletons) {
	analyseoverlaps[idts].PrintOverlapTable(output);
      }

      bool Result = true;
      if (hkl_list.num_datasets() != 1) {
	output.logTab(0,LOGFILE,
		      "==============================================================\n");
	Result = false;
      }
      // Print summary as a Results table if one dataset
      allsummarystatistics.PrintOneSummaryTable(idts, Result, output);
      output.logFlush();
    } // end loop datasets
    if (hkl_list.num_datasets() > 1) { // summary for multiple datasets
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
      // Analyse overlaps
      SetOverlapFlags(true, hkl_list); // include overlaps
      //# AnalyseOverlaps(hkl_list)
    }
    if (outputcontrols.UnMerged()) {
      // Output unmerged reflections file(s), including overlaps
      scala::WriteUnmergedOutputFiles(runTitle, hkl_list, SD_model,
				      NormRes.Imax(), outputcontrols, output);
    }
  }  // end try

  catch (phaser_io::PreprocessorError& capErr) {
    output.logWarning(LOGFILE, capErr.partialEcho()+"PREPROCESSOR ERROR: " + capErr.Message()
                      + "\n");
  }
  catch (phaser_io::SyntaxError& ccp4Err) {
    output.logWarning(LOGFILE, ccp4Err.Echo()+"\nSYNTAX ERROR: " + ccp4Err.Message()
                      + "\n");
  }
  catch (phaser_io::InputError& inpErr) {
    output.logWarning(LOGFILE, inpErr.Echo()+"\nINPUT ERROR: " + inpErr.Message()
                      + "\n");
  }

  catch (Message_fatal& message)
    {
      output.logWarning(LOGFILE, "\nFATAL ERROR message: \n"
                        + message.text() + "\n");
    } 

  catch (std::exception const& err) {
    output.logWarning(LOGFILE, "\nUNHANDLED EXCEPTION: " + std::string(err.what())+"\n");
  }
  catch (...) {
    output.logWarning(LOGFILE, "\nUNKNOWN EXCEPTION TYPE\n");
  }

  if (output.doXmlout()) output.logTab(0, LXML,"</AIMLESS>");

  Citation citation
    ("P.R.Evans and G.N.Murshudov, 'How good are my data and what is the resolution?'"+
		    std::string(" Acta Cryst. D69, 1204-1214  (2013)."),
     "http://journals.iucr.org/d/issues/2013/07/00/ba5190/index.html");
  output.logTab(0,LOGFILE, "\n"+citation.MakeLogCitation());


  output.logTab(0, LOGFILE,
		"\nEnd of aimless job, total time: "+overalltime.format(true)+"\n\n");

  return 0;

}
