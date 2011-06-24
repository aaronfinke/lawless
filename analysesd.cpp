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
		  const all_controls& controls,
		  const Normalise& NormRes,
		  const int& firstAnalysis,
		  phaser_io::Output& output)
  // On entry:
  //  SDM        SD correction model
  //  hkl_list   reflection data
  //  controls   
  //  NormRes    normalisation
  //  firstAnalysis
  //     = 0  two analyses expected, after rough scaling and main scaling
  //     = +1    same, after first analysis
  //     = -1 only one analysis expected ie first & last
  {
    //  plot       true to write normal probability plots, for each run, and partials & fulls
    //  NintensBinTarget  target number/intensity bin
    //  FixSdB     true to fix SdB even if it is variable in SDM
    //  tolerance  convergence
    //  max_cycles maximum number of cycles
    // Values for initial rough analysis
    bool plot = false;
    int NintensBinTarget = 100;
    bool FixSdB = true; 
    double tolerance = 0.001;
    int  max_cycles = 3;
    if (firstAnalysis != 0) { // values for finalanalysis
      plot = true;
      NintensBinTarget = 300;
      FixSdB = false; 
      tolerance = 0.0004;
      max_cycles = 10;  // 5
    }

    const bool REFINESD = true; // true to refine, false use Simplex
    if (firstAnalysis != 0) {
      output.logTab(0,LOGFILE,
		    std::string("\nOptimisation and analysis of standard deviations\n")+
		    "================================================\n");
    } else {
      output.logTab(0,LOGFILE,
		    std::string("\nFirst rough optimisation and analysis of standard deviations\n")+
		                  "============================================================\n");
    }

    if (controls.Anomalous && (controls.AnomalousSDcorr != controls.Anomalous)) {
      output.logTab(0,LOGFILE,
    "\nSeparation of anomalous I+ and I- suppressed in SD optimisation due to low multiplicity");
    }

    output.logTab(0,LOGFILE,"\n"+SDM.formatFullPartialInfo());

    // Intensity bins etc
    int NintBin = controls.analysis.NiBins();
    //   Number of bins, number of "reference" bin,
    //   intensity at "reference" bin, maximum intensity
    float Iav = NormRes.Imean();
    float Jmax = NormRes.Imax();
    IntensityBin Irange(NintBin, NintBin/2, Iav, Jmax);

    if (SDM.Refine()) {
      if (firstAnalysis <= 0) {
	// Initial correction from normal probability analysis
	UpdateSDMfromNPlot(SDM, hkl_list, controls, false, output);
	output.logTab(0,LOGFILE,
		      "\nSD correction parameters after normal probability correction\n"+SDM.format());
      } else {
	output.logTab(0,LOGFILE,
		      "\nCurrent SD correction parameters\n"+SDM.format());
      }
      // Select subset of reflections to speed up optimisation
      int nacc = SelectSDcorrReflections(hkl_list, controls, NintensBinTarget);
      ///      int nacc = hkl_list.num_reflections();
      output.logTabPrintf(0,LOGFILE,
	  "\n%7d reflections selected for SD optimisation out of %8d in file\n",
	  nacc, hkl_list.num_reflections());

      Timer timer;

      bool saveSdBfix = SDM.NoSDb();
      if (FixSdB) {
	SDM.SetNoSDb(true);
      }
	      
      // Optimise SD model
      if (REFINESD) {
	SDanalysis sdanal =
	  RefineSDcorrectionFactors(SDM, hkl_list, controls, Irange,
				    tolerance, max_cycles, output);
      } else {
	OptimiseSDcorr(SDM, hkl_list, controls, Irange,
		       tolerance, max_cycles, output);
      }
      output.logTab(0,LOGFILE,
	    "\nSD correction parameters after optimisation\n"+SDM.format());

      if (FixSdB) {
	SDM.SetNoSDb(saveSdBfix);  // restore saved fixSdB flag
      }
      
      hkl_list.ResetReflAccept();  // set to accept everything

      // I did try to do a final normal probability correction, but this may make it worse
      //      if (firstAnalysis != 0) {
      //	// Final correction of SDfac
      //	/////*/	UpdateSDMfromNPlot(SDM, hkl_list, controls, true, output);
      //	output.logTab(0,LOGFILE,
      //		      "\nSD correction parameters after 2nd normal probability correction\n"+
      //		      SDM.format());
      //      }
      double dtime = timer.Stop();
      output.logTabPrintf(0,LOGFILE,
			  "\nTime for SD optimisation = %10.1f secs\n",
			  dtime);
    } else { // norefine
	output.logTab(0,LOGFILE,
      "\nNo refinement of SD correction parameters\n"+SDM.format());      
	//^
	hkl_list.ResetReflAccept();  // set to accept everything
	// Select subset of reflections to speed up optimisation
	//	int nacc = SelectSDcorrReflections(hkl_list, controls, NintensBinTarget);
	//	output.logTabPrintf(0,LOGFILE,
	//			    "\n%7d reflections selected for SD residual out of %8d in file\n",
	//			    nacc, hkl_list.num_reflections());
	max_cycles = 0;
	SDanalysis sdanal =
	  RefineSDcorrectionFactors(SDM, hkl_list, controls, Irange,
					   tolerance, max_cycles, output);
	hkl_list.ResetReflAccept();  // set to accept everything
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
    AccumulateNormProb(SDM, hkl_list, controls.AnomalousSDcorr, normalprobanal);
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
	  int irun = SelObs.Run(i);  // run index for i'th observation (from 1)
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
    float sdrej2 = 5.0;
    scala::RejectFlags::Reject2Policy Rej2policy = scala::RejectFlags::REJECT;

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
  void UpdateSDMfromNPlot(SDmodel& SDM, const hkl_unmerge_list& hkl_list,
			  const all_controls& controls, const bool& fixup,
			  phaser_io::Output& output)
  // if fixup is true, then we are doing a final fix of the full||partial
  // values which were not optimised due to too few data, but still check
  // that there are enough data to do this
  {
    // for identifying the dataset for each observation
    std::vector<Run> Runs = hkl_list.RunList();
    int Nsets = SDM.Nsets();
    
    std::vector<NormalProbAnal> normalprobanal(2*Nsets);

    // Get all data into normal probability plots
    AccumulateNormProb(SDM, hkl_list, controls.AnomalousSDcorr, normalprobanal);

    const int TOOFEW = 20;  // minimum number for NP plot
    // Loop sets for analyses
    bool copyflag;
    for (int iset=0;iset<Nsets;iset++) {
      // use slope of selected data (central part of distribution)
      // to update SD correction models for each run
      float slopef = normalprobanal[2*iset].Slope();
      float slopep = normalprobanal[2*iset+1].Slope();
      int nf = normalprobanal[2*iset].Number(); // counts
      int np = normalprobanal[2*iset+1].Number();
      if (fixup) {
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
      } else {
	// initial values
	copyflag = true;
      }
      SDM.UpdateFactor(iset, copyflag, slopef, slopep);
      if (SDM.AllRunsSame()) {
	output.logTab(0,LOGFILE,"\nFor all runs, ");
      } else {
	output.logTabPrintf(0,LOGFILE,
			    "\nFor run %d, ", Runs[iset].RunNumber());
      }
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
    }
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

    SDcorrResidual sdr(SDM, hkl_list, IIrange, controls.AnomalousSDcorr);
    sdr(SDM.GetParameters());
    std::cout <<"\n****** All data SD plot ******\n\n";
    std::cout << "\n" << SDM.format() <<"\n";
    sdr.PrintTable(output);
  }
  //---------------------------------------------------------------
  //---------------------------------------------------------------
}  // end namespace scala

