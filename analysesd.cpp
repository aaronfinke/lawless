// analysesd.cpp

// // #include "testgradient.hh"

#include "analysesd.hh"
#include "normalprobanal.hh"
#include "plotfiles.hh"
#include "selectedobservations.hh"
#include "jiffy.hh"
#include "selectsdcorrreflections.hh"
#include "tablegraph.hh"
#include "optimisesdcorr.hh"
#include "refinesdcorrection.hh"
#include "file_util.hh"
#include "timer.hh"
#include "reject.hh"
using phaser_io::LOGFILE;

namespace scala
{
  //^--------------------------------------------------------------
  class DumpDeltaPairs {
  public:
    DumpDeltaPairs(); // open file

    void dump(SelectedObservations& SelObs);

    void end();
		
  private:      
    FILE* file;
    correl_coeff cc;
    LinearFit linefit;
  };
  //--------------------------------------------------------------
  DumpDeltaPairs::DumpDeltaPairs()
  {
    file = OpenFile("dumpdeltapairs.dat", true);
  }
  //--------------------------------------------------------------
  void DumpDeltaPairs::dump(SelectedObservations& SelObs)
  {
    if (SelObs.Number() <= 1) return;
    std::vector<float> delta = SelObs.Deviations(); // delta
    std::vector<float> delI = SelObs.DelI();  // I - <I>
    std::vector<float> sigI = SelObs.sigmaI(); // sigma(I)
    float an = SelObs.Number();
    float fac = sqrt((an-1.0)/an);
    //    float fac = 1.0;
    for (size_t i=0;i<delta.size();++i) {
      if (delta[i] != 0.0) {
	// delta2 = delI/[fac * sigma^2)
	float d2 = delI.at(i)/(fac * sigI.at(i));
	fprintf(file, "%8.4f %8.4f\n", delta[i], d2);
	cc.add(delta[i], d2);
	linefit.add(delta[i], d2, 1.0);
      }
    }
  }
  //--------------------------------------------------------------
  void DumpDeltaPairs::end() 
  {
    std::cout << "\nDumpDeltaPairs CC = " << cc.result().val;
    std::cout << "  slope = " << linefit.slope() << "\n";
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  void AnalyseSD (SDmodel& SDM, hkl_unmerge_list& hkl_list,
		  all_controls& controls,
		  const Normalise& NormRes,
		  const int& firstAnalysis,
		  phaser_io::Output& output)
  // On entry:
  //  SDM        SD correction model
  //  hkl_list   reflection data
  //  controls   may update controls.anomalouscontrol.AnomalousSDcorr
  //  NormRes    normalisation
  //  firstAnalysis
  //     = 0  two analyses expected, after rough scaling and main scaling
  //     = +1    same, after first analysis
  //     = -1 only one analysis expected ie first & last
  {
    const bool REFINESD = true; // true to refine, false use Simplex

    // Intensity bins etc
    int NintBin = controls.analysis.NiBins();
    //   Number of bins, number of "reference" bin,
    //   intensity at "reference" bin, maximum intensity
    float Iav = NormRes.Imean();
    float Jmax = NormRes.Imax();
    IntensityBin Irange(NintBin, NintBin/2, Iav, Jmax);
    bool plot = false;
    int  max_cycles = 3;
    double tolerance = 0.001;  // relative tolerance in parameter shift
    double rtolerance = 0.005; // tolerance on change of residual

    if (firstAnalysis != 0) {
      output.logTab(0,LOGFILE,
		    std::string("\nOptimisation and analysis of standard deviations\n")+
		    "================================================\n");
    } else {
      output.logTab(0,LOGFILE,
		    std::string("\nFirst rough optimisation and analysis of standard deviations\n")+
		                  "============================================================\n");
    }


    if (SDM.Refine()) {
      Timer timer;
      int npass = 1; // normally one pass through refinement at each stage
      if (firstAnalysis < 0) { // two passes if following restore; sdcorrection refine
	npass = 2;
      }
      output.logTab(0,LOGFILE,"\nWeighting scheme for averages: "+
		    SDM.formatWeightType());
      

      output.logTab(0,LOGFILE,"\n"+SDM.formatFullPartialInfo());

      if (firstAnalysis <= 0) {
	// Initial correction from normal probability analysis
	bool updated = true;
	while (true) {  // maybe more than one shot at this
	  SDMdataNumbers sdmnum =
	    UpdateSDMfromNPlot(SDM, hkl_list, controls, false, output);
	  //^^
	  //	  for (int i=0;i<sdmnum.nfnp.size();++i) {
	  //	    std::cout << "Set " << i << " " << sdmnum.nfnp[i].first
	  //		      <<" " << sdmnum.nfnp[i].second <<"\n";
	  //	  }
	  //^-
	  
	  // Are there some "sets" with no data?
	  if (!sdmnum.enoughdata) { // at least one set with insufficientdata
	    if (controls.anomalouscontrol.AnomalousSDcorr) {
	      // keeping I+ & I- separate, try combining them
	      controls.anomalouscontrol.AnomalousSDcorr = false;
	      output.logTab(0,LOGFILE,
			    "Try combining I+ and I-");
	      continue;
	    } else {
	      if (SDM.Nsets()) {
		// One runs or all runs same, and I+ and I- together, bail out
		output.logTab(0,LOGFILE,
		      "\n!!!! Insufficient data to determine SD correction factors, no refinement");
		SDM.SetRefine(false);
		updated = false;
		break;
	      } else {
		// try combining runs
		SDM.SetAllRunsSame(true);
		output.logTab(0,LOGFILE,
			    "Try combining runs");
		continue;
	      }
	    }
	  } else { // sufficent data for everything
	    break;
	  }
	} // end infinite while
	if (updated) {output.logTab(0,LOGFILE,
				    "\nSD correction parameters after normal probability correction\n"+SDM.format());
	}
      } else {  // not first
	output.logTab(0,LOGFILE,
		      "\nCurrent SD correction parameters\n"+SDM.format());
      }

      if (SDM.Refine()) { // may have changed!
	if (controls.anomalouscontrol.AnomalousSDcorr) {
	  output.logTab(0,LOGFILE,
			"\nI+ and I- will be kept separate in SD optimisation");
	} else {
	  output.logTab(0,LOGFILE,
			"\nSeparation of anomalous I+ and I- suppressed in SD optimisation due to low multiplicity");
	}


	// Save tie settings from SDM to restore later
	SDties savedties = SDM.Ties();
	
	for (int ipass=0;ipass<npass;++ipass) { // loop one or two passes
	  bool initialpass = (firstAnalysis == 0) || (npass > 1 && ipass == 0); // initial pass
	  // Set parameters for this pass
	  //   plot       true to write normal probability plots, for each run, and partials & fulls
	  //   NintensBinTarget  target number/intensity bin
	  //   FixSdB     true to fix SdB even if it is variable in SDM
	  //   tolerance  convergence
	  //   max_cycles maximum number of cycles
	  // Values for initial rough analysis
	  int NintensBinTarget = 100;
	  bool FixSdB = true; 
	  tolerance = 0.001;
	  rtolerance = 0.01;
	  max_cycles = 6;
	  if (initialpass) { // values for first analysis
	    // Put a tie on SDadd if there isn't one already, to help stabilise refinement
	    int ksdadd = savedties.targets.size()-1; // SDadd parameter is the last one
	    if (std::abs(savedties.sdtargets.at(ksdadd)) < 0.000001) {
	      double sdaddTarget = 0.02;    // target value for SDadd
	      double sdaddSDTarget = 0.1;   // and its SD
	      SDM.ResetTie(ksdadd, sdaddTarget, sdaddSDTarget);
	    }
	  } else if (!initialpass) { // values for final analysis
	    plot = true;
	    NintensBinTarget = 400;
	    FixSdB = false; 
	    tolerance = 0.0004;
	    rtolerance = 0.001;
	    max_cycles = 20;  // 5
	  }
	  
	  // ==== Outlier rejection, no check between I+ & I- (outlier.ndatasets = 0)
	  // Use special SD model for outlier rejection, with inflated SDadd since
	  // we don't know this, and we don't want to reject all the strong reflections
	  SDmodel SDMoutlier = SDM;
	  const double SDADDOUTLIER = 0.05;
	  SDMoutlier.SetSDadd(SDADDOUTLIER);
	  
	  OutlierControl outliercontrol(0);
	  outliercontrol.Combine() = false; // don't check between datasets
	  float sdrej = 30.0;
	  float sdrej2 = sdrej;
	  //./	  scala::RejectFlags::Reject2Policy Rej2policy = scala::RejectFlags::KEEP;
	  scala::RejectFlags::Reject2Policy Rej2policy = scala::RejectFlags::REJECT;
	  if (!initialpass) { // values for finalanalysis
	    sdrej = 30.0;
	    sdrej2 = sdrej;
	    //./	    Rej2policy = scala::RejectFlags::KEEP;
	    Rej2policy = scala::RejectFlags::REJECT;
	  }
	  RejectFlags rejflags(sdrej, sdrej2, Rej2policy);
	  outliercontrol.SetReject(rejflags, scala::ALL);
	  WriteRogues DummyRogues;
	  // Check for outliers & reject them
	  //  hkl_list is updated for status, but SDs are not changed
	  RejectOutlier(hkl_list, SDMoutlier,
			NormRes, controls.anomalouscontrol.AnomalousSDcorr,
			outliercontrol, DummyRogues);
	  std::vector<int> nrejs = CountOutliers(hkl_list);
	  output.logTabPrintf(0,LOGFILE,
			      "\nFor SD optimisation, number of outliers within I+ || I- sets: %6d,  between I+ & I- %6d, on |E|max %6d\n",
			      nrejs[0], nrejs[1], nrejs[2]);
	  // Select subset of reflections to speed up optimisation
	  int nacc = SelectSDcorrReflections(hkl_list, controls, NintensBinTarget);
	  output.logTabPrintf(0,LOGFILE,
			      "\n%7d reflections selected for SD optimisation out of %8d in file\n",
			      nacc, hkl_list.num_reflections());
	
	
	  bool saveSdBfix = SDM.NoSDb();
	  if (FixSdB) {
	    SDM.SetNoSDb(true);
	  }
	      
	  // Optimise SD model
	  if (REFINESD) {
	    SDanalysis sdanal =
	      RefineSDcorrectionFactors(SDM, hkl_list, controls, Irange,
					tolerance, rtolerance, max_cycles, output);
	    //^
	    //	  PrintSDanalysis(sdanal, SDanalysis(), RejectFlags(), Irange, hkl_list.RunList(),
	    //		       SDM, -1, PxdName(), false, output);
	    //^-
	  } else {
	    OptimiseSDcorr(SDM, hkl_list, controls, Irange,
			   tolerance, max_cycles, output);
	  }
	  output.logTab(0,LOGFILE,
			"\nSD correction parameters after optimisation\n"+SDM.format());

	  if (FixSdB) {
	    SDM.SetNoSDb(saveSdBfix);  // restore saved fixSdB flag
	  }
	
	  hkl_list.ResetReflAccept();  // set to accept (ie cancel SelectSDcorrReflections)
	  // Reset SDM ties
	  SDM.ResetTies(savedties);
	} // end loop one or two passes
	// Clear all outlier & other status flags (except ObsFlags)
	ClearObsStatus(hkl_list);

	// I did try to do a final normal probability correction, but this may make it worse
	//      if (firstAnalysis != 0) {
	//	// Final correction of SDfac
	//	/////*/	UpdateSDMfromNPlot(SDM, hkl_list, controls, true, output);
	//	output.logTab(0,LOGFILE,
	//		      "\nSD correction parameters after 2nd normal probability correction\n"+
	//		      SDM.format());
	//      }
	output.logTab(0,LOGFILE,
		      "\nTime for SD optimisation = "+timer.format(true));
      } // norefine because of no data 
    } else { // norefine explicit
	output.logTab(0,LOGFILE,
      "\nNo refinement of SD correction parameters\n"+SDM.format());      
	//^
	//	hkl_list.ResetReflAccept();  // set to accept everything
	//	// Select subset of reflections to speed up optimisation
	//	//	int nacc = SelectSDcorrReflections(hkl_list, controls, NintensBinTarget);
	//	//	output.logTabPrintf(0,LOGFILE,
	//	//			    "\n%7d reflections selected for SD residual out of %8d in file\n",
	//	//			    nacc, hkl_list.num_reflections());
	//	max_cycles = 0; // no refinement
	//	SDanalysis sdanal =
	//	  RefineSDcorrectionFactors(SDM, hkl_list, controls, Irange,
	//				    tolerance, rtolerance, max_cycles, output);
	//	hkl_list.ResetReflAccept();  // set to accept everything
	//^-
    } // end refine/norefine

    //^
    //    MakeSDplot(SDM, hkl_list, controls, NormRes, output);
    //^-

    // Normal probability analysis for each run
    std::vector<Run> Runs = hkl_list.RunList();
    int Nruns = hkl_list.num_runs();
    std::vector<NormalProbAnal> normalprobanal(2*Nruns);
    NormalProbPlot NPPlot;
    if (plot) {
      NPPlot.init("NORMPLOT", true, "Normal probability plot","");
    }

    // Replot with corrections
    AccumulateNormProb(SDM, hkl_list, controls.anomalouscontrol.AnomalousSDcorr, normalprobanal);
    if (plot) {
      for (int irun=0;irun<Nruns;irun++) {
	// Dataset name
	std::string dname = hkl_list.xdataset(Runs[irun].DatasetIndex()).pxdname().dname();
	std::string fp;
	for (int ip=2*irun;ip<=2*irun+1;ip++) {
	  if (ip == 2*irun) {
	    fp = " fulls";
	  } else {
	    fp = " partials";
	  }
	  std::string legend = "Run "+phaser_io::itos(hkl_list.RunList()[irun].RunNumber());
	  normalprobanal[ip].Plot(NPPlot, legend+": "+dname+fp);
	}
      }
      
      NPPlot.ClosePlot();
    }
  }  // AnalyseSD
  //--------------------------------------------------------------
  void AddInDelta(SelectedObservations& SelObs, const bool& allsame,
		  std::vector<NormalProbAnal>& normalprobanal)
  // if allsame true, treat all runs the same
  {
    std::vector<float> delta;
    int NumUsed = SelObs.Number();
    if (NumUsed > 1){
      int jnp;
      delta = SelObs.Deviations();
      for (size_t i=0;i<delta.size();i++)	{
	// use only delta != 0.0
	if (delta[i] != 0.0) {
	  int irun = SelObs.Run(i);  // run index for i'th observation (from 0)
	  if (allsame) irun = 0;
	  if (SelObs.Full(i)) {
	    jnp = 2*irun; // i'th observation is Full
	  } else {
	    jnp = 2*irun + 1; // i'th observation is partial
	  }
	  normalprobanal[jnp].AddDelta(delta[i]);
	}
      }
    }
  }	    
  //--------------------------------------------------------------
  int AccumulateNormProb(const SDmodel& SDM, const hkl_unmerge_list& hkl_list,
			  const bool& Anomalous,
			  std::vector<NormalProbAnal>& normalprobanal)
  // Go through reflection list, apply current SD correction model, &
  // accumulate deltas into normalprobanal (vector for full, partial for each run)
  // Returns number rejected
  {
    float sdrej = 5.0;     // for now, FIXME
    float sdrej2 = sdrej;
    scala::RejectFlags::Reject2Policy Rej2policy = scala::RejectFlags::KEEP;

    int Ndatasets = hkl_list.num_datasets();

    reflection this_refl;
    int Nrej = 0;
    int nref = 0;
    hkl_list.rewind();
    //^    DumpDeltaPairs ddp;

    bool allsame = SDM.AllRunsSame();  // true to combine parameters for all runs

    //   analyse deviations within each run & full/partial
    while (hkl_list.next_reflection(this_refl) >= 0)  {
      bool Centric = hkl_list.symmetry().is_centric(this_refl.hkl());
      //  Apply current SD correction to reflection (all observations)
      SDM.CorrectReflection(this_refl);  // using average intensity
      nref++;
      // loop datasets
      for (int id=0;id<Ndatasets;id++) {
	if (Centric || !Anomalous) {
	  // No anomalous, treat all observations together
	  SelectedObservations Selobs(this_refl, id, ALL);
	  // Reject outliers
	  Nrej += Selobs.Outliers(RejectFlags(sdrej, sdrej2, Rej2policy));
	  AddInDelta(Selobs, allsame, normalprobanal);
	  //^	  ddp.dump(Selobs);
	} else {
	  // Anomalous, treat I+ & I- separately
	  // FIXME outlier rejection between I+ & I- not done yet
	  SelectedObservations Selobs(this_refl, id, IPLUS);
	  // Reject outliers
	  Nrej += Selobs.Outliers(RejectFlags(sdrej, sdrej2, Rej2policy));
	  AddInDelta(Selobs, allsame, normalprobanal);
	  //^	  ddp.dump(Selobs);
	  Selobs.init(this_refl, id, IMINUS);
	  // Reject outliers
	  Nrej += Selobs.Outliers(RejectFlags(sdrej, sdrej2, Rej2policy));
	  AddInDelta(Selobs, allsame, normalprobanal);
	  //^	  ddp.dump(Selobs);
	}
      } // end loop datasets
    } // end loop reflections
    //^    ddp.end();

    return Nrej;
  }
  //--------------------------------------------------------------
  SDMdataNumbers  UpdateSDMfromNPlot(SDmodel& SDM, const hkl_unmerge_list& hkl_list,
				     const all_controls& controls, const bool& fixup,
				     phaser_io::Output& output)
  // if fixup is true, then we are doing a final fix of the full||partial
  // values which were not optimised due to too few data, but still check
  // that there are enough data to do this
  // Returns for each set count of number of fulls & number of partials used
  {
    // for identifying the dataset for each observation
    std::vector<Run> Runs = hkl_list.RunList();
    int Nsets = SDM.Nsets();
    
    SDMdataNumbers sdmnums(Nsets);

    std::vector<NormalProbAnal> normalprobanal(2*Nsets);
    std::vector<std::pair<float,float> > slopes(Nsets);

    // Get all data into normal probability plots
    AccumulateNormProb(SDM, hkl_list, controls.anomalouscontrol.AnomalousSDcorr, normalprobanal);

    const int TOOFEW = 20;  // minimum number for NP plot
    // Loop sets for analyses
    int leastdata = TOOFEW*2;
    bool enoughdata = true;
    for (int iset=0;iset<Nsets;iset++) {
      // use slope of selected data (central part of distribution)
      // to update SD correction models for each run
      slopes[iset].first  = normalprobanal[2*iset].Slope();
      slopes[iset].second = normalprobanal[2*iset+1].Slope();
      
      sdmnums.nfnp[iset].first  = normalprobanal[2*iset].Number(); // counts
      sdmnums.nfnp[iset].second = normalprobanal[2*iset+1].Number();
      if (sdmnums.Number(iset) <= TOOFEW) {
	// insufficient data for sensible fitting
	enoughdata = false;
	leastdata = Min(leastdata, sdmnums.Number(iset));
	sdmnums.enoughdata = false;
      }
    }
    if (!enoughdata) {
	output.logTabPrintf(0,LOGFILE,
	    "\n!!!! Insufficient data pairs! smallest number in run = %3d\n",
			    leastdata);
	return sdmnums;
    }

    bool copyflag;
    for (int iset=0;iset<Nsets;iset++) {
      float slopef = slopes[iset].first;
      float slopep = slopes[iset].second;
      int nf = sdmnums.nfnp[iset].first;
      int np = sdmnums.nfnp[iset].second;
      if (fixup) {  // never true at present!
	// We are fixing up values
	if (SDM.UseFlag(iset) == +2) {
	   // few partials, test number
	  if (np > TOOFEW) {
	    copyflag = false; // OK to accept value for partials
	    // don't change fulls
	    slopef = 1.0;
	  }
	} else if (SDM.UseFlag(iset) == -2) {
	   // few fulls, test number
	  if (nf > TOOFEW) {
	    copyflag = false; // OK to accept value for fulls
	    // don't change partials
	    slopep = 1.0;
	  }
	}
      } else { // not fixup
	// initial values
	copyflag = true;
      }
      if (SDM.AllRunsSame()) {
	output.logTab(0,LOGFILE,"\nFor all runs, ");
      } else {
	output.logTabPrintf(0,LOGFILE,
			    "\nFor run %d, ", Runs[iset].RunNumber());
      }

      SDM.UpdateFactor(iset, copyflag, slopef, slopep);
      if (nf > 0 && np > 0) {
	// both fulls and partials
	output.logTabPrintf(0,LOGFILE,
			    "slopes (full, partial) of central part of normal probability plot = %6.2f, %6.2f\n",
			    slopef, slopep);
	if (SDM.UseFlag(iset) == 0) {
	  output.logTab(0,LOGFILE,
			"  Correction applied to parameters for fulls and partials");
	} else if (SDM.UseFlag(iset) > 0) {
	  output.logTab(0,LOGFILE,
			"  Correction from fulls applied to parameters for fulls and partials");
	} else {
	  output.logTab(0,LOGFILE,
			"  Correction from partials applied to parameters for fulls and partials");
	}
      } else if (nf > 0) {
	// only fulls
	output.logTabPrintf(0,LOGFILE,
			    "slope of central part of normal probability plot = %6.2f\n",
			    slopef);
	output.logTab(0,LOGFILE,
		      "  Correction applied to parameters for fulls");
      } else if (np > 0) {
	// only partials
	output.logTabPrintf(0,LOGFILE,
			    "slope of central part of normal probability plot = %6.2f\n",
			    slopep);
	output.logTab(0,LOGFILE,
		      "  Correction applied to parameters for partials");
      }
      for (int ip=2*iset;ip<=2*iset+1;ip++) {
	normalprobanal[ip].Clear();  // reset normal probability plot
      }	
    } // end loop sets
  return sdmnums;
  }
  //---------------------------------------------------------------
  void MakeSDplot(SDmodel& SDM,
		  hkl_unmerge_list& hkl_list,
		  const all_controls& controls,
		  const Normalise& NormRes,
		  phaser_io::Output& output)
  // for debuggery
  {
    // Intensity bins etc
    int NNintBin = controls.analysis.NiBins();
    //   Number of bins, number of "reference" bin,
    //   intensity at "reference" bin, maximum intensity
    float IIav = NormRes.Imean();
    float JJmax = NormRes.Imax();
    IntensityBin IIrange(NNintBin, NNintBin/2, IIav, JJmax);

    SDcorrResidual sdr(SDM, hkl_list, IIrange, controls.anomalouscontrol.AnomalousSDcorr);
    sdr(SDM.GetParameters());
    std::cout <<"\n****** All data SD plot ******\n\n";
    std::cout << "\n" << SDM.format() <<"\n";
    sdr.PrintTable(output);
  }
  //---------------------------------------------------------------
  //---------------------------------------------------------------
}  // end namespace scala

