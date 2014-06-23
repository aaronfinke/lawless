// sdanalysis.cpp

#include "sdanalysis.hh"
#include "plotfiles.hh"
#include "jiffy.hh"
#include "tablegraph.hh"
#include "string_util.hh"
#include "optimisesdcorr.hh"

namespace scala
{
  //---------------------------------------------------------------
  void SDparameterGroup::init(const int& Numintbins, const int& Nparam)
  {
    numintbins = Numintbins;
    nparams = Nparam;
    mndelta.assign(numintbins, MeanSD());
    if (nparams > 0) {
      sumddeltadp.resize(nparams);
      sumddelta2dp.resize(nparams);
      nsumddeltadp.resize(nparams);
      
      for (int i=0;i<nparams;++i) {
	sumddeltadp[i].assign(numintbins, 0.0);
	sumddelta2dp[i].assign(numintbins, 0.0);
	nsumddeltadp[i].assign(numintbins, 0);
      }
    }
  }
  //---------------------------------------------------------------
  // Add in delta for intensity bin
  void SDparameterGroup::AddDelta(const int& mint, const float& delta)
  {
    mndelta[mint].Add(delta);
  }
  //---------------------------------------------------------------
  // Add in derivative for intensity bin and full/partial
  // kpl is local parameter number
  void SDparameterGroup::AddDerivative(const int& mint, const int& kpl,
				 const float& dddp, const float& delta)
  {
    // sum( d(delta(i in j)) / d(p(k)))
    sumddeltadp[kpl][mint] += dddp;
    // sum( d(delta[i] in j^2) / d(p(k)))
    sumddelta2dp[kpl][mint] +=
      2. * delta * dddp;
    nsumddeltadp[kpl][mint]++;
  }
  //---------------------------------------------------------------
  //! return sd(delta) for each intensity bin
  std::vector<double> SDparameterGroup::SDdelta() const
  {
    std::vector<double> sdd(numintbins);
    for (int i=0;i<numintbins;++i) {
      sdd[i] = mndelta[i].SD();
    }
    return sdd;
  }
  //---------------------------------------------------------------
  //! return MeanSD(delta) for each intensity bin
  std::vector<MeanSD> SDparameterGroup::MeanDelta() const
  {
    return mndelta;
  }
  //---------------------------------------------------------------
  std::vector<std::vector<double> > SDparameterGroup::dSigmaDeltaDp() const
  // return d(sigma(delta[inb]))/dp[k] for intensity bins inb
  // outer vector length nparams (k), inner length numintbins (inb)
  {
    std::vector<std::vector<double> > dsddp;
    if (nparams == 0) return dsddp;
    dsddp.resize(numintbins);
    for (int inb=0;inb<numintbins;++inb) {
      std::vector<double>  dsddpk(nparams, 0.0); // for each parameter
      if (mndelta[inb].Count() > 2) {
	double meanDelta = mndelta[inb].Mean();
	double sigDel = mndelta[inb].SD();
	for (int k=0;k<nparams;++k) {
	  double ninb = nsumddeltadp[k][inb]; // N(inb) 
	  // d(Var(delta(inb)))/dp = (1/N)Sum(d(delta[inb]^2) / d(p(k))) -
	  //     (2/N) <delta[inb]> Sumd(delta[inb]) / d(p(k))
	  if (ninb > 0.0) {
	    double dvar = (1./ninb)*(sumddelta2dp[k][inb] -
				     2.*meanDelta*sumddeltadp[k][inb]);
	    // dsigma/dp = (1/2sigma)dVar/dp
	    dsddpk[k] = 0.5*dvar/sigDel;
	  }
	} // end loop parameters
      }
      dsddp[inb] = dsddpk;
    } // end loop intensity bins
    return dsddp;
  }
  //---------------------------------------------------------------
  SDparameterGroup& SDparameterGroup::operator +=
  (const SDparameterGroup& other)
  {
    ASSERT (numintbins == other.numintbins);
    ASSERT (numintbins == int(mndelta.size()));
    ASSERT (nparams == other.nparams);
    for (int i=0;i<numintbins;++i) {
      mndelta[i] += other.mndelta[i];
    }
    if (nparams > 0) {
      for (int j=0;j<nparams;++j) {
	for (int i=0;i<numintbins;++i) {
	  sumddeltadp[j][i] += other.sumddeltadp[j][i];
	  sumddelta2dp[j][i] += other.sumddelta2dp[j][i];
	  nsumddeltadp[j][i] += other.nsumddeltadp[j][i];
	}
      }
    }
    return *this;
  }
  //---------------------------------------------------------------
  SDparameterGroup& operator +
  (const SDparameterGroup& a, const SDparameterGroup& b)
  {
    SDparameterGroup c = a;
    return c += b;
  }
  //---------------------------------------------------------------
  //---------------------------------------------------------------
  // Construct and clear
  SDanalysis::SDanalysis(const IntensityBin& Irange,  const SDmodel& SDM,
			const bool& Allsamerun, const bool& Derivatives)
  {
    init(Irange, SDM, Allsamerun, Derivatives);
  }
  //---------------------------------------------------------------
  // Initialise and clear
  void SDanalysis::init(const IntensityBin& Irange, const SDmodel& SDM,
			const bool& Allsamerun, const bool& Derivatives)
  {
    irange = Irange;
    numintensitybins = Irange.NumberBins();
    numruns = SDM.Nruns();
    numrunsused = numruns;
    allsamerun = Allsamerun;
    if (allsamerun) numrunsused = 1;  // combine all runs into one
    nparams = 0;  // no derivatives
    if (Derivatives) {
      nparams = SDM.Nparams();
    }
    fulls  = false;
    partials = false;
    deltalimit = 20.; // ignore delta above this limit
    rungroupnumberfulls.resize(numrunsused);
    rungroupnumberpartials.resize(numrunsused);

    numparametergroups = 0;
    useflag.resize(numrunsused);
    for (int irun=0;irun<numrunsused;++irun) { // loop runs
      useflag[irun] = SDM.UseFlag(irun);
      // number of groups in run, 1 or 2
      int nparrun = (SDM.UseFlag(irun) == 0) ? 2 : 1;
      for (int i=0;i<nparrun;++i) {
	grouprunnumber.push_back(irun);  // run number for each group
	bool groupisfull = true;
	if (nparrun > 1 && i>0) groupisfull = false;
	if (nparrun == 1 && SDM.UseFlag(irun) < 0) groupisfull = false;
	groupfull.push_back(groupisfull);
      }
      if (nparrun == 1) {
	// group number for run
	if (SDM.UseFlag(irun) > 0) { // only fulls, point partials to fulls
	  rungroupnumberfulls[irun] = numparametergroups;
	  rungroupnumberpartials[irun] = numparametergroups;
	} else { // only partials, point fulls to partials
	  rungroupnumberfulls[irun] = numparametergroups;
	  rungroupnumberpartials[irun] = numparametergroups;
	}
      } else {
	// fulls and partials
	rungroupnumberfulls[irun] = numparametergroups;
	rungroupnumberpartials[irun] = numparametergroups+1;
      }
      numparametergroups  += nparrun;
    } // end loop runs
    ASSERT (int(grouprunnumber.size()) == numparametergroups);
    nbinclasses = numintensitybins * numparametergroups;

    sdparametergroup.resize(numparametergroups);
    for (int jpc=0;jpc<numparametergroups;++jpc) { // loop parameter groups
      if (Derivatives) {
	int irun = grouprunnumber[jpc];
	int npargrp =  SDM.NparamsFull(irun);
	if (!groupfull[jpc]) npargrp = SDM.NparamsPartial(irun);
	sdparametergroup[jpc].init(numintensitybins, npargrp);
      } else {
	sdparametergroup[jpc].init(numintensitybins, 0);
      }
    } // end loop parameter groups

    if (Derivatives) {
      // Make list of local parameter indices for each global parameter
      idxlocal.resize(nparams);
      int kpg = 0; // global
      int kpl = 0; // local
      idxparametergroup.resize(numparametergroups);
      for (int jpc=0;jpc<numparametergroups;++jpc) { // loop parameter groups
	int npargrp = sdparametergroup[jpc].Nparam();
	// index into global parameter list of 1st parameter for this group
	idxparametergroup[jpc] = kpg;
	kpl = 0; // local
	for (int i=0;i<npargrp;++i) {
	  // local parameter index for each global parameter
	  idxlocal[kpg++] = kpl++;
	}
      }
      ASSERT (nparams == kpg);
      //^
      //      std::cout << "SDanalysis Global/local indices ";
      //      for (int i=0;i<nparams;++i) {
      //	std::cout << " : " << i << " " << idxlocal[i];
      //      }
      //      std::cout <<"\n";
      //^-
    }
  }
  //---------------------------------------------------------------
  //! returns true if object is empty
  bool SDanalysis::Empty() const
  {
    if (numintensitybins == 0) return true;
    if (!(fulls || partials)) return true;
    return false;
  }
  //---------------------------------------------------------------
  //! return true if there are fulls and they were used (cf useflag)
  bool SDanalysis::FullsUsed(const int& irun) const
  {
    return (fulls && (useflag.at(irun) >= 0));
  }
  //---------------------------------------------------------------
  //! return true if there are partials and they were used (cf useflag)
  bool SDanalysis::PartialsUsed(const int& irun) const
  {
    return (partials && (useflag.at(irun) <= 0));
  }
  //---------------------------------------------------------------
  int SDanalysis::ParameterGroup(const int& irun, const bool& full) const
  // parameter group number for run and full/partial
  {
    int ir = irun;
    if (allsamerun) ir = 0;
    if (full) {
      return rungroupnumberfulls[ir];
    } else {
      return rungroupnumberpartials[ir];
    }
  }
  //---------------------------------------------------------------
  bool SDanalysis::AddDelta(const float& delta, const int& mint,
			    const int& irun, const bool& full)
  // add in delta, for intensity bin mint, run irun, full == false for partial
  {
    if (std::abs(delta) > deltalimit) return false;
    int ir = irun;
    if (allsamerun) ir = 0;
    if (full) {
      fulls = true;
    } else {
      partials = true;
    }
    //^
    //    if (ParameterGroup(irun, full) < 0 || ParameterGroup(irun, full) >= numparametergroups) {
    //      std::cout <<"help\n";
    //    } //^-
    sdparametergroup[ParameterGroup(irun, full)].AddDelta(mint, delta);
    return true;
  }
  // ------------------------------------------------------------
  void SDanalysis::AddSelobsDelta(SelectedObservations& selobs,
				  const int& mint)
  // Add in delta1 contributions for intensity bin mint
  //  delta1 = (Ihl - <Ih>others)/sqrt(SD(Ihl)^2 + SD(Ihothers)^2)
  //   where <Ih>others is the average over all observations of reflection h,
  //   excluding Ihl itself and SD(Ihothers) is its ESD
  //   cf AddSelobsDelta2 which uses delta2 definition
  {
    if (selobs.Number() > 1) {
      std::vector<float> delta = selobs.Deviations();
      for (size_t i=0;i<delta.size();++i) {
	if (delta[i] != 0.0) { // valid delta
	  AddDelta(delta[i], mint, selobs.Run(i), selobs.Full(i));
	  //^
	  //	  //	  if (mint == numintensitybins-1) {
	  //	  std::cout <<"AddSelobsDelta " << mint <<" "
	  //		    << selobs.hkl().format() << " " << delta[i] <<" " <<used;
	  //	  if (std::abs(delta[i]) > 3.0) std::cout << " ****";
	  //	  std::cout <<"\n";
	  //	  } //^-
	}
      }
    }
  }
  // ------------------------------------------------------------
  void SDanalysis::AddSelobsDelta2(SelectedObservations& selobs,
				  const int& mint)
  // Add in delta2 contributions for intensity bin mint
  //  delta2 = sqrt(n/n-1) (Ihl - <Ih>)/SD(Ihl)
  //   where <Ih> is the average over all observations of reflection h,
  //   including Ihl itself
  //   cf AddSelobsDelta which uses delta1 definition
  {
    if (selobs.Number() > 1) {
      std::vector<float> delta2 = selobs.Delta2();

      for (size_t i=0;i<delta2.size();++i) {
	if (delta2[i] != 0.0) { // valid delta
	  //^
	  //	  std::cout <<"AddSelobsDelta2 " << selobs.hkl().format()
	  //		    <<" " << mint << " " << delta2[i]
	  //		    <<" run "<<selobs.Run(i) <<" full "<<selobs.Full(i)
	  //		    <<" "<<selobs.Reflection().get_observation(i).I()
	  //		    <<" "<<selobs.Reflection().get_observation(i).sigI()
	  //		    <<" "<<i<<std::endl;
	  //	  }
	  AddDelta(delta2[i],
		   mint, selobs.Run(i), selobs.Full(i));
	  //^
	  //	  if (mint == numintensitybins-1) {
	  //	    std::cout <<"AddSelobsDelta2 " << selobs.hkl().format()
	  //		      <<" " << mint << " " << delta2[i]
	  //		      <<" "<<selobs.Reflection().get_observation(i).I()
	  //		      <<" "<<selobs.Reflection().get_observation(i).sigI()
	  //		      <<"\n";
	  //	  } //^-
	}
      }
    }
  }
  // ------------------------------------------------------------
  void SDanalysis::AddDerivatives(SelectedObservations& selobs,
				  const int& mint,
				  const std::vector <std::vector<double> >& ddeltadp)
  // Add in to sums for derivatives
  //  ddeltadp[iobs][k] is d(delta(iobs))/dp(k) for the iobs'th observation in selobs
  //  p(k) is the k'th parameter of nparams (here k is global parameter index)
  {
    if (selobs.Number() > 1) {
      std::vector<float> delta = selobs.Delta2();
      for (size_t iobs=0;iobs<delta.size();++iobs) {  // loop observations
	if (delta[iobs] != 0.0) { // valid delta
	  int irun = selobs.Run(iobs);
	  if (allsamerun) irun = 0;
	  int jpc = ParameterGroup(irun, selobs.Full(iobs));
	  ASSERT (int(ddeltadp[iobs].size()) == nparams);
	  for (int k=0;k<nparams;++k) { // loop global parameters k
	    if (ddeltadp[iobs][k] != 0.0) {
	      int kpl = idxlocal[k]; // local parameter number within run/full/partial
	      sdparametergroup[jpc].AddDerivative(mint, kpl,
						     ddeltadp[iobs][k], delta[iobs]);
	    }
	  } // end loop global parameters k
	}
      }  // end loop observations
    }
  }
  //---------------------------------------------------------------
  //! SDs for each "bin class" jpc, 
  std::vector<double> SDanalysis::SDdelta() const
  {
    std::vector<double> sdd;
    std::vector<double> sd;
    for (int jpc=0;jpc<numparametergroups;++jpc) { // loop parameter groups
      sd = sdparametergroup[jpc].SDdelta(); // for each intensity bin
      if (sd.size() > 0) sdd.insert(sdd.end(), sd.begin(), sd.end());
    }
    ASSERT (int(sdd.size()) == nbinclasses);
    return sdd;
  }
  //---------------------------------------------------------------
  std::vector<MeanSD> SDanalysis::MeanDelta() const
  //!< MeanSD for deltas for all bin classes
  {
    std::vector<MeanSD> mdd;
    std::vector<MeanSD> md;
    for (int jpc=0;jpc<numparametergroups;++jpc) { // loop parameter groups
      md = sdparametergroup[jpc].MeanDelta(); // for each intensity bin
      if (md.size() > 0) mdd.insert(mdd.end(), md.begin(), md.end());
    }
    ASSERT (int(mdd.size()) == nbinclasses);
    return mdd;
  }
  //---------------------------------------------------------------
  std::vector<std::vector<std::vector<double> > > SDanalysis::Derivatives() const
  // SDs and their derivatives w.r.t. parameters p(k) for
  // each "bin class" jpc
  // A bin class jpc is run/intensity/full|partial
  //   returns dsigDeldp[jpc][inb][k] partial derivatives d(sigma(delta(jpc)))/dp(k) for
  //       parameter class jpc, intensity bin inb, parameter k, k is index local to class jpc
  {
    // d(sigma(delta(jpc,inb))/dp(k)
    std::vector<std::vector<std::vector<double> > > dsigDeldp;

    std::vector<std::vector<double> > dsdp;   // for each intensity bin for each local parameter
    for (int jpc=0;jpc<numparametergroups;++jpc) { // loop parameter groups
      dsdp = sdparametergroup[jpc].dSigmaDeltaDp(); // for each local parameter
      if (dsdp.size() > 0) dsigDeldp.push_back(dsdp);
    }    
    ASSERT (int(dsigDeldp.size()) == numparametergroups);
    return dsigDeldp;
  }
  //---------------------------------------------------------------
  // Return one element
  MeanSD SDanalysis::GetMeanSD(const int& mint,
			       const int& irun, const bool& full) const
  {
    int jpc = ParameterGroup(irun, full);
    if (jpc >= 0) {
      return sdparametergroup[jpc].MeanDelta()[mint];
    }
    return MeanSD();  // no information
  }
  //---------------------------------------------------------------
  // Return one element
  MeanSD SDanalysis::GetMeanSD(const int& mint,
			       const int& irun, const int& fullpart) const
  {
    return GetMeanSD(mint, irun, (fullpart == 0));
  }
  //---------------------------------------------------------------
  // Return vector for all intensity bins
  std::vector<MeanSD> SDanalysis::GetMeanSD(const int& irun, const bool& full) const
  {
    int jpc = ParameterGroup(irun, full);
    if (jpc >= 0) {
      return sdparametergroup[jpc].MeanDelta();
    }
    return std::vector<MeanSD>();  // no information
  }
  //---------------------------------------------------------------
  std::vector<int> SDanalysis::NumberinIntbins() const
  //! return number of observations in each intensity bin, over all runs and full/partials
  {
    std::vector<int> ninbins(numintensitybins, 0);
    for (int irun=0;irun<numrunsused;++irun) {
      for (int mint=0;mint<numintensitybins;++mint) {
	if (fulls) ninbins[mint] +=
	  GetMeanSD(mint, irun, true).Count();
	if (partials) ninbins[mint] +=
	  GetMeanSD(mint, irun, false).Count();
      }
    }
    return ninbins;
  }
  //---------------------------------------------------------------
  std::vector<int> SDanalysis::NumberinClass() const
  //! return number of observations in each bin class
  {
    std::vector<int> ninbins(nbinclasses, 0);
    for (int jpc=0;jpc<numparametergroups;++jpc) { // loop parameter groups
      for (int mint=0;mint<numintensitybins;++mint) {
	ninbins[BinClass(jpc,mint)] =
	  sdparametergroup[jpc].MeanDelta()[mint].Count();
      }
    }
    return ninbins;
  }
  //---------------------------------------------------------------
  std::pair<int,int> SDanalysis::Number() const
  //! return total number of observations, fulls & partials
  {
    std::pair<int,int> num(0,0);
    for (int jpc=0;jpc<numparametergroups;++jpc) { // loop parameter groups
      for (int mint=0;mint<numintensitybins;++mint) {
	if (groupfull[jpc]) {
	  num.first +=
	    sdparametergroup[jpc].MeanDelta()[mint].Count();
	} else {
	  num.second +=
	    sdparametergroup[jpc].MeanDelta()[mint].Count();
	}
      }
    }
    return num;
  }
  //---------------------------------------------------------------
  //! return bin class for parameter class and intensity bin
  int SDanalysis::BinClass(const int& jpc, const int& mint) const
  {
    return jpc*numintensitybins + mint;
  }
  //---------------------------------------------------------------
  //! return number of parameters in parameter class
  int SDanalysis::Nparam(const int& paramgroup) const
  {
    return sdparametergroup[paramgroup].Nparam();
  }
  //---------------------------------------------------------------
  //! return index in complete parameter list to first parameter in class
  int SDanalysis::IdxParam(const int& paramgroup) const
  {
    return idxparametergroup[paramgroup];
  }
  //---------------------------------------------------------------
  //---------------------------------------------------------------
  void AddMsdData(std::vector<MeanSD>& a, const std::vector<MeanSD>& b)
  // add vector b to a
  {
    ASSERT (a.size() == b.size());
    for (size_t i=0;i<a.size();++i) {
      a[i] += b[i];
    }
  }
  //---------------------------------------------------------------
  void PrintSDanalysis(const SDanalysis& sdanal1, const SDanalysis& sdanal2,
		       const RejectFlags& rejflags,
		       const IntensityBin& Irange,
		       const std::vector<Run>& runlist,
		       const SDmodel& SDM,
		       const int& datasetIndex, const PxdName& dataset_pxd,
		       const bool& fullprint,
		       phaser_io::Output& output)
  // Print table from one or two SDanalysis objects
  //   sdanal1 if both are present this is for the "core" data
  //   sdanal2 if both are present this is for the "core" data, else null
  // fullprint == false for brief printing
  {
    if (fullprint) {
      output.logTab(0,LOGFILE,
		  std::string("\nAnalysis of standard deviations\n")+
		  "===============================\n"+
		  "This analyses the distribution of the normalised deviations\n"+
		  "Delta = (Ihl - Mn(Iothers) )/sqrt[sd(Ihl)**2 + sd(Mn(I))**2]\n"+
		  "If the SD is a true estimate of the error, this distribution should have\n"+
		  " Mean=0.0 and Sigma=1.0 for all ranges of intensity\n"		  
		  "\nThe analysis is repeated for ranges of increasing Imean\n"+
		  "The Mean is expected to increase with Imean since the latter\n"+
		  "is a weighted mean and sd(Ihl) & Ihl are correlated\n"+
		  "\nIf the Sigma increases with Imean, increase the value of SdAdd\n\n"); 	   
    }
    if (sdanal1.Empty()) {
      // no data
      output.logTab(0,LOGFILE,"\n!!!! No data !!!!\n\n");
      return;
    }
    // Is there a second SDanalysis object?
    bool secondsdanal = !sdanal2.Empty();

    output.logTab(0,LOGFILE,
		  "SD corrections:- SdFac * Sqrt[sd](I**2 + SdB I + (SdAdd I)**2)\n");

    //^ Residuals
    //    SDresiduals resids = SDcorrResidual(sdanal1).Residual();
    //    std::cout << "Overall residual: " << resids.overallR
    //	      << " number " << resids.nR <<"\n";
    //    for (int i=0;i<resids.classRs.size();++i) {
    //      std::cout << "Residual for class " << i+1
    //		<< " " << resids.classRs[i]
    //		<< " number " << resids.cnR[i] <<"\n";
    //    }
    //^-

    // Print SD correction factors
    output.logTab(0,LOGFILE,SDM.format());

    // numruns is total number
    int numruns = sdanal1.NumberRuns();
    bool allsamerun = sdanal1.AllSameRun();

    // If 2nd SDanalysis object, check for compatibility
    if (secondsdanal) {
      ASSERT (numruns == sdanal2.NumberRuns());
      ASSERT (allsamerun == sdanal2.AllSameRun());
    }

    int nsd = 1;  // number of full/partial/1st/2nd datasets

    // full/partial flag, +1 fulls, -1 partials, 0 both
    int fullpartialall = -1000;

    int kr=0;

    std::vector<bool> whichones(4, false);  // for full1, partial1, full2, partial2

    for (int ir=0;ir<numruns;++ir) {   // Loop all runs
      int irun = ir;
      if (allsamerun) irun = 0;
      if (datasetIndex < 0 || datasetIndex == runlist[ir].DatasetIndex() ||
	  (allsamerun && ir==0)) {
	kr++;  // count runs printed

	bool fulls = sdanal1.FullsUsed(irun);
	bool partials = sdanal1.PartialsUsed(irun);
	bool fulls2 = false;
	bool partials2 = false;
	if (secondsdanal) {
	  fulls2 = sdanal2.FullsUsed(irun);
	  partials2 = sdanal2.PartialsUsed(irun);
	}

	// full/partial flag, +1 fulls, -1 partials, 0 both
	int fullpartial = 0;
	if (!fulls) fullpartial = -1;
	if (!partials) fullpartial = +1;
	if (secondsdanal) {
	  if (fullpartial < 0 && fulls2) fullpartial = 0;
	  if (fullpartial > 0 && partials2) fullpartial = 0;
	}
	if (fullpartialall < -999) {
	  fullpartialall = fullpartial;
	} else if (fullpartialall != fullpartial) {
	  fullpartialall = 0;
	}

	std::vector<std::vector<MeanSD> > msddata;
	if (fulls) {	// fulls for 1st object
	  msddata.push_back(sdanal1.GetMeanSD(irun, true));
	  whichones[0] = true;
	}  
	if (partials) { // partials for 1st object
	  msddata.push_back(sdanal1.GetMeanSD(irun, false));
	  whichones[1] = true;
	}
	if (secondsdanal) {
	  if (fulls2) {	// fulls for 2nd object
	    msddata.push_back(sdanal2.GetMeanSD(irun, true));
	    whichones[2] = true;
	  }  
	  if (partials2) { // partials for 2nd object
	    msddata.push_back(sdanal2.GetMeanSD(irun, false));
	    whichones[3] = true;
	  }
	}
	std::string ttitle = " Run "+clipper::String(runlist[irun].RunNumber(),3)+
	  ", standard deviation v. Intensity, "+
	  dataset_pxd.dname();

	PrintSDanalysisTable(rejflags, Irange, ttitle, fullpartial, secondsdanal,
			     msddata, output);
      }
    }  // end loop runs

    nsd = 0; // count objects
    for (size_t i=0; i<whichones.size(); i++) { 
      if (whichones[i]) {nsd++;}
    }

    if (kr > 1 && fullprint) {
      // > 1 run in this dataset, print totals over all relevant runs
      std::vector<std::vector<MeanSD> > msddata(nsd);
      for (int i=0;i<nsd;++i) {
	msddata[i].assign(Irange.NumberBins(),MeanSD());
      }

      int kf1 = -1; // indices into msddata array for fulls 1
      int kp1 = -1; //  " partials 1
      int kf2 = -1; //  " fulls 2
      int kp2 = -1; //  " partials 2
      if (fullpartialall == 0) { // both full & partial
	kf1 = 0; kp1 = 1;
	if (nsd == 4) {kf2 = 2; kp2 = 3;}
      } else if (fullpartialall == -1) { // only partials
	kp1 = 0;
	if (nsd == 2) {kp2 = 1;}
      } else if (fullpartialall == +1) { // only fulls
	kf1 = 0;
	if (nsd == 2) {kf2 = 1;}
      }
      
      // full/partial flag, +1 fulls, -1 partials, 0 both
      //      int fullpartial = 0;

      for (int ir=0;ir<numruns;++ir) {   // Loop all runs
	int irun = ir;
	if (allsamerun) irun = 0;
	if (datasetIndex < 0 || datasetIndex == runlist[ir].DatasetIndex() ||
	    (allsamerun && ir==0)) {
	  if (sdanal1.FullsUsed(irun)) {	// fulls for 1st object
	    AddMsdData(msddata.at(kf1), sdanal1.GetMeanSD(irun, true));
	  }  
	  if (sdanal1.PartialsUsed(irun)) { // partials for 1st object
	    AddMsdData(msddata.at(kp1), sdanal1.GetMeanSD(irun, false));
	  }
	  if (secondsdanal) {
	    if (sdanal2.FullsUsed(irun)) {	// fulls for 2nd object
	      AddMsdData(msddata.at(kf2), sdanal2.GetMeanSD(irun, true));
	    }  
	    if (sdanal2.PartialsUsed(irun)) { // partials for 2nd object
	      AddMsdData(msddata.at(kp2), sdanal2.GetMeanSD(irun, false));
	    }
	  }
	}
      } // end loop all runs
      std::string ttitle = " All runs, standard deviation v. Intensity, "+
	dataset_pxd.dname();

      PrintSDanalysisTable(rejflags, Irange, ttitle, fullpartialall, secondsdanal,
			   msddata, output);
    }  // multiple runs
  }
  //---------------------------------------------------------------
  void PrintSDanalysisTable(const RejectFlags& rejflags,
			    const IntensityBin& Irange,
			    const std::string& ttitle,
			    const int& fullpartial,
			    const bool& outer,
			    const std::vector<std::vector<MeanSD> >& msdanal,
			    phaser_io::Output& output)
  // Print table for one run
  //  Irange       intensity binning
  //  gtitle       header title for table
  //  fullpartial  +1 fulls, -1 partials, 0 both
  //  outer        true if "outer" data present 
  //  msdanal    vector for each intensity for each full/partial/core/outer
  {
    const int MINNUMBER = 8; // minimum number to print statistic
    // number of analysis sets outer/core full/partial (1,2, or 4)
    int nanalsets;
    bool both = false;
    if (fullpartial == 0) both = true;
    if (both) {nanalsets = 2;}
    else {nanalsets = 1;}
    if (outer) {nanalsets *= 2;}
    ASSERT (nanalsets == int(msdanal.size()));
    
    int c[] =  {2,5,8,11,14}; // possible column numbers for graphs
    int c1[] =  {2,8,14};      // possible column numbers for graphs
    int cc[] =  {2,11,14};    // possible column numbers for graphs

    TableGraph table;
    table.init(ttitle);
    output.logTab(0,LOGFILE,table.formatTitle());
    table.StoreID("Graph-SDanalysis");

    std::vector<int> cln(c,c+nanalsets+1);
    std::string graphtitle = " Sigma(scatter/SD)";
    std::string corelimit;
    if (outer) {
      corelimit = clipper::String(rejflags.sdrej, 3, 2) + " sd";
      graphtitle += ", within "+corelimit;
      if (both) {
	cln.assign(cc,cc+nanalsets/2+1); // {2,11,14}
      } else {
	cln.assign(c1,c1+nanalsets/2+1); // {2,8,14}
      }
    };
    output.logTab(0,LOGFILE, table.Graph(graphtitle,"N",cln));
    if (outer) {
      corelimit = clipper::String(rejflags.sdrej, 3, 2) + " sd";
      graphtitle += ", all and within "+corelimit;
      cln.assign(c,c+nanalsets+1);
      output.logTab(0,LOGFILE, table.Graph(graphtitle,"N",cln));
    };

    std::string fpclabel;

    std::vector<std::string> collabels;
    collabels.push_back("Range");         // 1
    collabels.push_back("Mn(I)");           // 2
    for (int i=0;i<nanalsets;++i) {
      std::string ctyp;
      if (both) {
	ctyp = "F";  //full
	if (i%2!=0) { // partial
	  ctyp = "P";
	}
      } else {  // not both
	if (fullpartial > 0) { // full
	  ctyp = "F";  //full
	} else { 
	  ctyp = "P";
	}
      }
      std::string fpcl;
      if (ctyp == "F") {
	fpcl = "   Fulls";
      } else {
	fpcl = "Partials";
      }
      if (outer) {
	if (i >= nanalsets/2) {
	  ctyp += "c";  // core data only
	  fpcl += ", <"+corelimit;
	} else {
	  fpcl += ", all";
	}
      }
      collabels.push_back("N"+ctyp);      // 3, 6,  9, 12
      collabels.push_back("Mn"+ctyp);     // 4, 7, 10, 13
      collabels.push_back("Sd"+ctyp);    // 5, 8, 11, 14
      fpclabel += StringUtil::CentreString(fpcl, 21); // centre in 21char field
    }
    int nc = collabels.size(); // number of columns

    // flags for possible null entries
    bool z[] = {false, false, false, false, true, false, false, true,
		false, false, true, false, false, true};
    std::vector<bool> Zero(z, z+nc);
    // format omitting 1st 2
    std::string fmt =
      "%9d%6.2f%6.2f%9d%6.2f%6.2f%9d%6.2f%6.2f%9d%6.2f%6.2f\n";
    if (nanalsets == 1) {
      fmt = "%9d%6.2f%6.2f\n";
    } else if (nanalsets == 2) {
      fmt = "%9d%6.2f%6.2f%9d%6.2f%6.2f\n";
    }

    output.logTab(0,LOGFILE,"\n              "+fpclabel);
    output.logTab(0,LOGFILE,
		  table.ColumnFields(collabels, Zero, "%4d%8.0f"+fmt));

    std::vector<MeanSD> mnsdoverall(nanalsets);  // overall values

    for (int mint=0;mint<Irange.NumberBins();++mint) {
      // each set has a count, mean & SD
      std::vector<int> mcount(nanalsets);
      std::vector<double> mmean(nanalsets);
      std::vector<double> msd(nanalsets);
      for (int i=0;i<nanalsets;++i) {
	mcount[i] = msdanal[i][mint].Count();
	mmean[i] = msdanal[i][mint].Mean();
	msd[i] = msdanal[i][mint].SD();
	if (mcount[i]  < MINNUMBER) msd[i] = 0.0;
	mnsdoverall[i] += msdanal[i][mint];
      }
      if (nanalsets == 1) {
	output.logTab(0,LOGFILE,
		      table.Line(nc, mint+1, Irange.mean(mint),
				 mcount[0], mmean[0], msd[0]));
      } else if (nanalsets == 2) {
	output.logTab(0,LOGFILE,
		      table.Line(nc, mint+1, Irange.mean(mint),
				 mcount[0], mmean[0], msd[0],
				 mcount[1], mmean[1], msd[1]));
      } else if (nanalsets == 4) {
	output.logTab(0,LOGFILE,
		      table.Line(nc, mint+1, Irange.mean(mint),
				 mcount[0], mmean[0], msd[0],
				 mcount[1], mmean[1], msd[1],
				 mcount[2], mmean[2], msd[2],
				 mcount[3], mmean[3], msd[3]));
      }
    } // end loop intensity bins

    output.logTab(0,LOGFILE,
		  table.CloseTable());
	  
    fmt = "Overall:    "+fmt+"\n";
    if (nanalsets == 1) {
      output.logTabPrintf(0,LOGFILE,fmt.c_str(),
			  mnsdoverall[0].Count(),
			  mnsdoverall[0].Mean(),
			  mnsdoverall[0].SD());
    } else if (nanalsets == 2) {
      output.logTabPrintf(0,LOGFILE,fmt.c_str(),
			  mnsdoverall[0].Count(),
			  mnsdoverall[0].Mean(),
			  mnsdoverall[0].SD(),
			  mnsdoverall[1].Count(),
			  mnsdoverall[1].Mean(),
			  mnsdoverall[1].SD());
    } else if (nanalsets == 4) {
      output.logTabPrintf(0,LOGFILE,fmt.c_str(),
			  mnsdoverall[0].Count(),
			  mnsdoverall[0].Mean(),
			  mnsdoverall[0].SD(),
			  mnsdoverall[1].Count(),
			  mnsdoverall[1].Mean(),
			  mnsdoverall[1].SD(),
			  mnsdoverall[2].Count(),
			  mnsdoverall[2].Mean(),
			  mnsdoverall[2].SD(),
			  mnsdoverall[3].Count(),
			  mnsdoverall[3].Mean(),
			  mnsdoverall[3].SD());
    }
  }
  //---------------------------------------------------------------
  //---------------------------------------------------------------
}  // end namespace scala
