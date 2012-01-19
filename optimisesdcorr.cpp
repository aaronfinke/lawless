// optimisesdcorr.cpp
//
#include "optimisesdcorr.hh"
#include "analysesd.hh"

namespace scala {
  //---------------------------------------------------------------
  SDcorrResidual::SDcorrResidual(SDmodel& sdm,
				 const hkl_unmerge_list& Hkl_list,
				 IntensityBin& Irange,
				 const bool& Anomalous)
  {
    init(sdm,Hkl_list,Irange,Anomalous);
  }
  //---------------------------------------------------------------
  void SDcorrResidual::init(SDmodel& sdm,
			    const hkl_unmerge_list& Hkl_list,
			    IntensityBin& Irange,
			    const bool& Anomalous)
  {
    // Store addresses of data objects
    SDM = &sdm;
    hkl_list = &Hkl_list;
    irange = Irange;
    anomalous = Anomalous;
    nrunsused = hkl_list->num_runs();
    if (SDM->AllRunsSame()) {
      nrunsused = +1;
    }
    combine.resize(nrunsused);
    for (int ir=0;ir<nrunsused;++ir) {  // runs
      combine[ir] = (SDM->UseFlag(ir) != 0); // true if few or no fulls||partials
    }
  }
  //---------------------------------------------------------------
  void SDcorrResidual::init(const SDanalysis& Sdanal)
  {
    sdanal = Sdanal;
    //  addresses of null data objects
    SDM = NULL;
    hkl_list = NULL;

    irange = Sdanal.Irange();
    anomalous = false;
    nrunsused = Sdanal.NumberRunsUsed();
    combine.resize(nrunsused);
    for (int ir=0;ir<nrunsused;++ir) {  // runs
      // true if few or no fulls||partials
      combine[ir] = !(sdanal.FullsUsed(ir) && sdanal.PartialsUsed(ir));
    }
  }
  //---------------------------------------------------------------
  double SDcorrResidual::operator() (const std::vector<double>& args) const
  // Return function value given the parameter vector args
  {
    ASSERT (SDM != NULL);

    // Update SDmodel
    SDM->SetParameters(args);

    float sdrej = 5.0;     // for now, FIXME
    float sdrej2 = 5.0;
    scala::RejectFlags::Reject2Policy Rej2policy = scala::RejectFlags::REJECT;

    int Ndatasets = hkl_list->num_datasets();

    sdanal.init(irange, *SDM, SDM->AllRunsSame());

    reflection this_refl;
    int Nrej = 0;
    int nref = 0;
    hkl_list->rewind();
    std::vector<float> delta, delta2;

    //  deviations within each run & full/partial
    while (hkl_list->next_reflection(this_refl) >= 0)  {
      bool Centric = hkl_list->symmetry().is_centric(this_refl.hkl());
      SDM->CorrectReflection(this_refl);
      // Average I <I> over all observations
      SelectedObservations Selobs(this_refl, -1, ALL);
      float Iav = Selobs.AverageScaleWt().I();  // average intensity for SD correction 

      int mint = irange.bin(Iav);
      nref++;
      // loop datasets
      for (int id=0;id<Ndatasets;id++) {
	if (Centric || !anomalous) {
	  // No anomalous, treat all observations together
	  SelectedObservations Selobs(this_refl, id, ALL);
	  // Reject outliers
	  Nrej += Selobs.Outliers(RejectFlags(sdrej, sdrej2, Rej2policy));
	  sdanal.AddSelobsDelta(Selobs, mint);
	} else {
	  // Anomalous, treat I+ & I- separately
	  // FIXME outlier rejection between I+ & I- not done yet
	  SelectedObservations Selobs(this_refl, id, IPLUS);
	  // Reject outliers
	  Nrej += Selobs.Outliers(RejectFlags(sdrej, sdrej2, Rej2policy));
	  sdanal.AddSelobsDelta(Selobs, mint);
	  Selobs = SelectedObservations(this_refl, id, IMINUS);
	  // Reject outliers
	  Nrej += Selobs.Outliers(RejectFlags(sdrej, sdrej2, Rej2policy));
	  sdanal.AddSelobsDelta(Selobs, mint);
	}
      } // end loop datasets
    } // end loop reflections

    return Residual().overallR;
  }
  //---------------------------------------------------------------
  SDresiduals SDcorrResidual::Residual() const
  {
    ASSERT (!sdanal.Empty());

    // Residual = 1/nrunsused SumRuns [(SumIntensity/Full/Part w*(1-Sigma)^2)/Sum w]
    // w = Sqrt Number in intensity bin (all runs)  (for now anyway)
    double sumr = 0.0;
    double sumw = 0.0;

    // Weight array for each intensity bin
    // relative weight for main residual compared to restraint residual
    const double WTREL = 0.01;
    // "SD" of residual in each intensity bin, for weighting by 1/SD^2
    const double SDRESID = 0.04;

    // weight for each intensity bin, equal (unit) weights
    double w1 = WTREL/(SDRESID*SDRESID*(irange.NumberBins()));
    std::vector<double> wib(irange.NumberBins(), w1);

    //%/  sqrt(N) weights
    //%/    std::vector<double> wib(irange.NumberBins(), 0.0);
    //%/    for (int i=0;i<irange.NumberBins();++i) { // intensity bins
    //%/      for (int ir=0;ir<nrunsused;++ir) {  // runs
    //%/	for (int j=0;j<2;++j) {  // full/partial
    //%/	  wib[i] += double(sdanal.GetMeanSD(i,ir,j).Count()); // Sum N
    //%/	}
    //%/      }
    //%/      wib[i] = sqrt(wib[i]);
    //%/    }

    double r;
    int nr = 0;
    std::vector<double> classRs; // for each parameter group
    std::vector<int> nRs; // for each parameter group

    for (int ir=0;ir<nrunsused;++ir) {  // runs
      std::vector<double> sumrrun(2, 0.0); // full, partial
      std::vector<double> sumwrun(2, 0.0); // full, partial
      std::vector<int> ninrun(2,0);
      for (int i=0;i<irange.NumberBins();++i) { // intensity bins
	if (combine[ir]) {
	  // combine full & partial (or just one)
	  MeanSD mnsd = sdanal.GetMeanSD(i,ir,0) + sdanal.GetMeanSD(i,ir,1);
	  r = 1.0 - mnsd.SD();
	  sumrrun[0] += wib[i] * r * r;
	  sumwrun[0] += wib[i];
	  ninrun[0] += mnsd.Count();
	} else {
	  for (int j=0;j<2;++j) {  // full/partial
	    r = 1.0 - sdanal.GetMeanSD(i,ir,j).SD();
	    sumrrun[j] += wib[i] * r * r;
	    sumwrun[j] += wib[i];
	    ninrun[j] += sdanal.GetMeanSD(i,ir,j).Count();
	  }
	}
      } // end loop intensity bins
      for (int j=0;j<2;++j) {  // full/partial
	if (sumwrun[j] > 0.0) {
	  sumr += sumrrun[j];
	  sumw += sumwrun[j];
	  nr += ninrun[j];
	  double resid = 0.5*sumrrun[j]/sumwrun[j];
	  classRs.push_back(resid);
	  nRs.push_back(ninrun[j]);
	}
      }
    }  // end loop runs

    double resid = 0.0;
    if (sumw > 0.0) resid = 0.5*sumr/sumw;
    return SDresiduals(resid, nr, classRs, nRs);
  }
  //---------------------------------------------------------------
  void SDcorrResidual::PrintTable(phaser_io::Output& output)
  {
    PrintSDanalysis(sdanal, SDanalysis(), RejectFlags(),
	       irange, hkl_list->RunList(), *SDM, -1, PxdName(), output);
  }
  //---------------------------------------------------------------
  SDcorrRefine::SDcorrRefine(SDmodel& sdm,
			     const hkl_unmerge_list& Hkl_list,
			     IntensityBin& Irange,
			     const bool& anomalous)
  {
    SDCresid.init(sdm,  Hkl_list, Irange, anomalous);
    nparams = sdm.Nparams();
  }
  //---------------------------------------------------------------
  double SDcorrRefine::operator() (const std::vector<double>& args) const
  {
    return SDCresid(args);
  }
  //---------------------------------------------------------------
  //---------------------------------------------------------------
  //---------------------------------------------------------------
  std::vector<std::vector<double> > StartValues(const SDmodel& SDM,
						const double& scale)
  {
    std::vector<double> params = SDM.GetParameters();
    // Initial shifts for each parameter type, scaled by "scale"
    std::vector<double> shifts = SDM.GetShifts(scale);
    std::vector<std::vector<double> > starts;
    starts.push_back(params);  // initial parameters
    int np = params.size();
    for (int i=0;i<np;++i) {
      std::vector<double> pars = params;
      // Shift i'th parameter
      pars[i] += shifts[i];
       starts.push_back(pars);
    }
    return starts;
  }
  //---------------------------------------------------------------
  void OptimiseSDcorr(SDmodel& SDM,
		      const hkl_unmerge_list& hkl_list,
		      const all_controls& controls,
		      IntensityBin& Irange,
		      const double& tolerance,
		      const int&  max_cycles,
		      phaser_io::Output& output)
  // Optimise SD correction model
  {
    // Set up Simplex minimiser
    Optimiser_simplex::TYPE type = Optimiser_simplex::NORMAL;
    Optimiser_simplex optimiserSimplex(tolerance, max_cycles, type);
    optimiserSimplex.debug(0);
    // Object to calculate residuals
    SDcorrRefine SDCref(SDM, hkl_list, Irange, controls.anomalouscontrol.Anomalous);
    // Set up vector of vectors of starting values
    std::vector<std::vector<double> > start = StartValues(SDM, 1.0);
    std::vector<double> newparams = optimiserSimplex(SDCref, start);
    SDM.SetParameters(newparams);
    output.logTabPrintf(0,LOGFILE,
	"\nAfter %4d optimisation cycles, final residual = %7.4f\n",
			optimiserSimplex.ncycles(),
			optimiserSimplex.bestResidual());
  }
  //---------------------------------------------------------------
}
