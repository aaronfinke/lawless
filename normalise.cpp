// normalise.cpp

// Conversion factor to convert I to E is done by run
// Within each run, the dividing factor is a function of s2 = sSqr = 1/d^2 and
// "time" (typically relative rotation as a proxy for time)
//
//    c(s2, t) = c0 exp[(B0 + B*t)*s2]
// ie the "B-factor" is a linear function of time t
//
// For determination, this is linearised by taking logs
//    ln(c) = ln(c0) + B0*s2 + B*s2*t
//
// Normalisation is then E = (1/c) * I

#include "normalise.hh"
#include "score_datatypes.hh"
#include "scala_util.hh"
#include "string_util.hh"

namespace scala {
  //--------------------------------------------------------------
  void BinSums::add(const Rtype& I, const Rtype& sSqr, const Rtype& Time) {
    // Add in I & s^2
    sum_I += I;
    sum_sSqr += sSqr;
    sum_time += Time-time0;
    max_time = Max(Time, max_time);
    N++;
  }
  //--------------------------------------------------------------
  // <I>
  double BinSums::MeanI() const {
    if (N == 0) return 0.0;
    return sum_I/double(N);
  }
  //--------------------------------------------------------------
  // <s^2>
  double BinSums::MeanSSqr() const {
    if (N == 0) return 0.0;
    return sum_sSqr/double(N);
  }
  //--------------------------------------------------------------
  // <time> (relative to time0)
  double BinSums::MeanTime() const {
    if (N == 0) return 0.0;
    return sum_time/double(N);
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  void BfactorModel::Setup(const bool& Constant,
			   const double& Time0) {
      if (Constant) {nparam = 2;}
      else {nparam = 3;}
      time0 = Time0;
    }
  //--------------------------------------------------------------
  // Store parameters
  void BfactorModel::Set(const std::vector<double>& params) {
    ASSERT (int(params.size()) == nparam);
    lnC0 = params[0];
    B0 = params[1];
    if (nparam > 2) {
      B = params[2];
    }
    else {B = 0.0;}
    valid = true;
    //^
    //    std::cout << "BfactorModel::Set " << params[0] << " " << params[1];
    //    if (nparam > 2) {std::cout << " " << params[2];}
    //    std::cout << "\n";
    //^-
  }
  //--------------------------------------------------------------
  // Return measurement vector = (1.0  sSqr,  sSqr*t)
  std::vector<double> BfactorModel::MeasurementVector(const double& sSqr, const double& t) {
    std::vector<double> mv(nparam);
    mv[0] = 1.0;
    mv[1] = sSqr;
    if (nparam > 2) {mv[2] = sSqr * t;}
    return mv;
  }
  //--------------------------------------------------------------
  // Return parameters
  std::vector<double> BfactorModel::Params() const {
    std::vector<double> params(nparam);
    params[0] = lnC0;
    params[1] = B0;
    if (nparam > 2) params[2] = B;
    return params;
  }
  //--------------------------------------------------------------
  // Return dividing correction factor for resolution sSqr &
  // "time" (actual time, not relative to start
  double BfactorModel::Factor(const double& sSqr, const double& t) const {
    double lc = lnC0 + (B0 + B*(t-time0))*sSqr;
    return exp(lc);
  }
  //--------------------------------------------------------------
  // time relative to start
  double BfactorModel::FactorT(const double& sSqr, const double& t) const {
    double lc = lnC0 + (B0 + B*t)*sSqr;
    return exp(lc);
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  Normalise::Normalise(const std::vector<clipper::Array2d<BinSums> >& sums_st)
    : E2max(30.0f)
  {
    type = +2;
    Nruns = sums_st.size();
    ASSERT (Nruns > 0);
    // set "Overall" true if there is no run or batch dependence
    // ie only one run & one batch
    // sums_st.cols are batches
    // sums_st.rows are resolution bins
    bool Overall = false;
    if (Nruns == 1 && sums_st[0].cols() == 1) {
      Overall = true;
    }

    RunFactor.resize(Nruns);   // separate factors for each run (LSQlinear)
    Bfactors.resize(Nruns);    // BfactorModel
    double w = 1.0;   // unit weights
    // Minimum number of batches for variation of B-factor with time
    const int MIN_NBATCH_TIME_VARIATION = 20;
    int max_nparam = -1;

    // Maximum Number of resolution bins in any run
    Nbins = -1;
    for (int ir=0;ir<Nruns;ir++) {        // loop runs
      Nbins = Max(Nbins, sums_st[ir].rows());
    }
    sSqrmnI.resize(Nbins);
    for (int is=0;is<Nbins;is++) {sSqrmnI[is].clear();}
    // Average "time" for each run
    std::vector<double> avgtime(Nruns, 0.0);

    // + + + + + + + + + + + + 
    for (int ir=0;ir<Nruns;ir++) {        // loop runs
      int nbatch = sums_st[ir].cols();    // Number of batches in run
      float Time0 = sums_st[ir](0,0).Time0();  // Time start point
      float DTime = 0.0;  // Time length

      // Number of parameters for each run = 2 or 3
      //  If the number of batches (proportional to time range) is "small", then it does not
      //  make sense to have a time variation in normalisation (ie radiation damage
      //  correction
      bool Constant = false;      
      if (nbatch < MIN_NBATCH_TIME_VARIATION) {Constant = true;}
      Bfactors[ir].Setup(Constant, Time0);  // store time0
      int Nparam = Bfactors[ir].Nparam();   // 2 or 3
      max_nparam = Max(max_nparam, Nparam);
      RunFactor[ir].clear(Nparam);
      std::vector<double> mv(Nparam);        // measurement vector

      for (int is=0;is<Nbins;is++) {         // loop resolution bins
	for (int ib=0;ib<nbatch;ib++) {   // loop batches
	  // Mean I & mean sSqr for this bin (is,ib)
	  double mnI    = sums_st[ir](is,ib).MeanI();
	  double mnSSqr = sums_st[ir](is,ib).MeanSSqr();
	  double mnTime = sums_st[ir](is,ib).MeanTime(); // from start
	  DTime = Max(DTime, sums_st[ir](is,ib).MaxTime());// maximum time from start of run
	  if (mnI > 0.001) {
	    // Measurement vector = (1 s2 [s2*t])
	    mv = Bfactors[ir].MeasurementVector(mnSSqr, mnTime);
	    RunFactor[ir].add(log(mnI), mv, w);  // add into LSQ sums
	  }
	}
      }
      // Solve for parameters for this run
      if (RunFactor[ir].Nobs() > 0) {
	Bfactors[ir].Set(RunFactor[ir].solve()); // solve & store parameters
	//  The number of parameters depends on "Constant" flag, 2 or 3
      } else {
	Bfactors[ir].Invalid();
      }
      DTime -= Time0;
      avgtime[ir] = 0.5*DTime;       // average time for run (relative)
      Bfactors[ir].DTime() = DTime;  // store maximum time for run
      //^
      //      std::cout << "run, time0, dtime " << ir << " "
      //      		<< Time0 << " " << DTime << "\n";
    }   // end loop runs
    // + + + + + + + + + + + + 

    // Get average correction term, over all runs & time,
    // ie just as function of s^2
    std::vector<double> params(max_nparam,0.0);
    std::vector<double> p(max_nparam);
    int n=0;
    for (int ir=0;ir<Nruns;ir++) {        // loop runs
      if (Bfactors[ir].Valid()) {
	p = Bfactors[ir].Params();
	// Reset B-factor parameters to correspond to mid-time point
	if (max_nparam > 2) {
	  p[1] += avgtime[ir] * p[2];
	  p[2] = 0.0;
	}
	for (int i=0;i<max_nparam;i++) {params[i] += p[i];}
	n++;
      }
    }
    if (n > 0) {
      for (int i=0;i<max_nparam;i++) {
	params[i] /= float(n);
      }
    }
    AvgFactor.Setup(true);
    params.resize(2);
    AvgFactor.Set(params);  // average correction factor just function of s^2
    validAll = true;
    validRun = true;
    if (Overall) validRun = false;

    if (type == +2) {
      MeanSD ImeanSD;  // for overall <I>
      imax = 0.0;
      
      for (int ir=0;ir<Nruns;ir++) {        // loop runs
	if (Bfactors[ir].Valid()) {
	  int nbatch = sums_st[ir].cols();    // Number of batches in run
	  for (int is=0;is<Nbins;is++) {         // loop resolution bins
	    for (int ib=0;ib<nbatch;ib++) {   // loop batches
	      if (sums_st[ir](is,ib).num() > 0) {
		// Mean I & mean sSqr for this bin (is,ib)
		double mnI    = sums_st[ir](is,ib).MeanI();
		double mnSSqr = sums_st[ir](is,ib).MeanSSqr();
		ImeanSD.Add(sums_st[ir](is,ib).MeanI());
		imax = Max(imax, sums_st[ir](is,ib).MeanI());
		// Apply correction so far (for time relative to start)
		double factor = Bfactors[ir].FactorT(mnSSqr, sums_st[ir](is,ib).MeanTime());
		if (factor != 0.0) {
		  // sums for overall resolution
		  sSqrmnI[is].add(mnI/factor, mnSSqr);
		}
	      }
	    }
	  }
	}
      }
      
      imean = ImeanSD.Mean();   // <I>
      sSqrmnIcorr.resize(Nbins+1);
      // First point at s = 0. at intercept of B-factor (= 1)
      sSqrmnIcorr[0].first = 0.0;
      sSqrmnIcorr[0].second = 1.0;
      // Minimum value to avoid negatives
      bcmin = 1.0;
      int ismax = -2;
      
      int k = 1;
      for (int i=0;i<Nbins;i++)  {
	if (sSqrmnI[i].num() > 0) {
	  sSqrmnIcorr[k].first = sSqrmnI[i].MeanSSqr();  // s^2
	  double mnI = sSqrmnI[i].MeanI();    // <I> corrected
	  sSqrmnIcorr[k].second = mnI;    // <I> corrected
	  if (mnI > 0.00001) {
	    bcmin = Min(bcmin, mnI);
	    // Maximum bin
	    ismax = k;
	  }
	  k++;	      
	}
      }
      if (ismax > 0) {sSqrmnIcorr.resize(ismax+1);}
      bincorr = Spline(sSqrmnIcorr);
      bcmin = Max(0.5*bcmin, 0.1);  // minimum correction factor (dividing) = 1/2 minimum
      // this avoids negative corrections
    }
  }
  //--------------------------------------------------------------
  float Normalise::apply(const float& I, const float& sSqr,
			 const int& irun, const float& time) const
  // "time" here is actual time or rotation, not relative to start of run:
  // offset to make relative to start of run is done internally
  {
    // Dividing scale
    float c = Bfactors.at(irun).Factor(sSqr, time);
    if (type == +1) 
      return I / c;
    if (type == +2)
      {
	float scorr = Max(bcmin, bincorr.Interpolate(sSqr));
	return I / (c * scorr);
      }
    return 0.0; // dummy
  }
  //--------------------------------------------------------------
  IsigI Normalise::apply(const IsigI& Is, const float& sSqr,
			 const int& irun, const float& time) const
  // "time" here is actual time or rotation, not relative to start of run:
  // offset to make relative to start of run is done internally
  {
    // Dividing scale
    float c = Bfactors.at(irun).Factor(sSqr, time);
    float scorr;
    if (type == +2)
      scorr = Max(bcmin, bincorr.Interpolate(sSqr));
    else
      scorr = 1.0;
    IsigI IsScl(Is);
    IsScl.scale(1./(c * scorr));
    return IsScl;
  }
  //--------------------------------------------------------------
  float Normalise::applyAvg(const float& I, const float& sSqr) const
  // Average correction
  {
    if (!validAll)
      Message::message(Message_fatal("Applying invalid B-Normalise"));
    // Dividing scale
    float c = AvgFactor.Factor(sSqr);
    if (type == +1) 
      return I / c;
    if (type == +2)
      {
	float scorr = Max(bcmin, bincorr.Interpolate(sSqr));
	return I / (c * scorr);
      }
    return 0.0; // dummy
  }
  //--------------------------------------------------------------
  IsigI Normalise::applyAvg(const IsigI& Is, const float& sSqr) const
  {
    if (!validAll)
      Message::message(Message_fatal("Applying invalid B-Normalise"));
    // Dividing scale
    float c = AvgFactor.Factor(sSqr);
    float scorr;
    if (type == +2)
      scorr = Max(bcmin, bincorr.Interpolate(sSqr));
    else
      scorr = 1.0;
    IsigI IsScl(Is);
    IsScl.scale(1./(c * scorr));
    return IsScl;
  }
  //--------------------------------------------------------------
  // Total correction factor (multiplying factor)
  float Normalise::Corr(const float& sSqr,
			const int& irun, const float& time) const
  {
    float c = apply(1.0, sSqr, irun, time);
    return c;
  }
  //--------------------------------------------------------------
  int normdumpnumber = 0;
  void NormDump(const Normalise& NormRes,
		const std::vector<int>& Batch0,
		const int& Nbin,
		const std::vector<clipper::Array2d<BinSums> >& sums_st)
  {
    // **** Open files for dumping  ****
    FILE* file;
    normdumpnumber++;
    std::string name = StringUtil::Strip("I"+StringUtil::itos(normdumpnumber,3))+".plot";
    file = fopen(name.c_str(), "w");
    if (file == NULL)
      Message::message(Message_fatal("Can't open file "+name));
    
    int Nruns = sums_st.size();
    // Constant term (intercept), dividing scale
    float SclC = NormRes.SclCor();
    std::vector<float> vsSqr;
    
    for (int ir=0;ir<Nruns;ir++) {        // loop runs
      fprintf(file, "# Run %6d\n", ir);
      int nbatch = sums_st[ir].cols();    // Number of batches in run
      for (int ib=0;ib<nbatch;ib++) {   // loop batches
	fprintf(file, "# Batch %6d\n", ib);
	for (int is=0;is<Nbin;is++) {         // loop resolution bins
	  if (sums_st[ir](is,ib).num() > 0) {
	    float mnI  = sums_st[ir](is,ib).MeanI();
	    if (mnI > 0.01) {
	      float mnsSqr = sums_st[ir](is,ib).MeanSSqr();
	      vsSqr.push_back(mnsSqr);
	      float Icorr = NormRes.apply(mnI, mnsSqr, ir, float(ib+Batch0[ir]));
	      fprintf(file, " %8.5f %8.4f  %10.6f %8d\n",
		      mnsSqr,mnI/SclC,Icorr,sums_st[ir](is,ib).num());
	    }
	  }
	}
	fprintf(file, "&\n");
      }
      // Interpolated values
      const int DIV = 5;
      if (vsSqr.size() > 0) {
	for (size_t i=0;i<vsSqr.size()-1;++i) {
	  for (int j=0;j<DIV;++j) {
	    float s = vsSqr[i] + j*(vsSqr[i+1]-vsSqr[i])/float(DIV);
	    float Icorr = NormRes.apply(1.0, s, ir, 0.0);
	    fprintf(file, " %8.5f %10.6f\n",
		    s, Icorr);
	  }
	}
      }
      fprintf(file, "&\n");
    }
    /*
      int nb = 50;
      double sbin = (ResRange.SResHigh() - ResRange.SResLow())/float(nb);
      double Idummy = 1.0;
      for (int i=0;i<nb;i++)
      {
      double s = ResRange.SResLow() + i * sbin;
      double corr = NormRes.apply(Idummy, s);
      fprintf(file, " %8.5f %10.6f\n", s, corr);
      }
      fprintf(file, "&\n");
    */
    
    int status = fclose(file);
    status = status;
  // **** end dump  ****
  }
  //--------------------------------------------------------------
  Normalise SetNormalise(const hkl_unmerge_list& hkl_list,
			 const double& MinIsigRatio,
			 const bool& Overall,
			 ResoRange& ResRange,
			 Rings& Icerings,
			 const int PrintLevel)
  // Set up intensity normalisation object NormRes
  //
  // Use binned <I> to get normalisation object
  // Binned of resolution (from ResRange) and batch within each run
  // Optionally output plot to file lnI.plot (if Printlevel > 0)
  // Return updated ResoRange, setting resolution cutoffs and bins
  //   (if MinIsigRatio < 0, don't change resolution range)
  // If Overall true, generate a single normalisation for all runs, ie
  // assume the data is already scaled
  //
  // MinIsigRatio   minimum I/sigI ratio on unaveraged data
  //                if < 0, do not reset resolution range

  {
    observation this_obs;
    int Nbin = ResRange.Nbins();
    
    // Sums for overall I/sigI etc by resolution for resolution cut-off
    std::vector<double> sum_I(Nbin);
    std::vector<double> sum_sigI(Nbin);
    std::vector<double> sum_sSqr(Nbin);
    std::vector<int> n_I(Nbin);
    
    // Ice rings
    Icerings.ClearSums();
    float IceTolerance = 4.0;
    
    // Binning by resolution and "time" (often rotation relative to start of run)
    // for each run
    int Nruns = hkl_list.num_runs();
    std::vector<Run> Runs = hkl_list.RunList();
    if (Overall) {Nruns = 1;}     // for overall normalisation, only one "run"
    std::vector<clipper::Array2d<BinSums> > sums_st(Nruns);
    std::vector<int> Batch0(Nruns);
    std::vector<int> NBatch(Nruns);
    
    for (int ir=0;ir<Nruns;ir++) {
      std::pair<int,int> batchrange = Runs[ir].BatchRange(); // min & max batch number
      int nbatches = batchrange.second - batchrange.first + 1; // not all may be present
      Batch0[ir] = batchrange.first;
      if (Overall) {   // for overall normalisation, only one "batch"
	nbatches = 1;
	Batch0[ir] = 1;
      }
      NBatch[ir] = nbatches;
      sums_st[ir] = clipper::Array2d<BinSums>(Nbin, nbatches);
      for (int i=0;i<Nbin;i++) {
	for (int j=0;j<nbatches;j++) {
	  sums_st[ir](i,j).clear(double(Runs[ir].TimeRange().min()));  // Set time0
	}}
    }
    
    int irun;
    int ib;

    reflection this_refl;
    hkl_list.rewind();  // Just in case
    int numobs = 0; // count contributions
    int zonalrefs = 0;
    std::vector<MeanSD> resmnI(Nbin); // <I>
    std::vector<MeanSD> resmnSdI(Nbin); // <sdI>

    //^
    //    FILE* ndump = OpenFile("NORMDUMP", true); //^-

    while (hkl_list.next_reflection(this_refl) >= 0) {
      // use only general reflections h!=k!=l!=0 to avoid
      // problems with unknown epsilon
      if (this_refl.hkl().IsGeneral()) {
	// Resolution 
	float sSqr = this_refl.invresolsq();
	int rbin = ResRange.tbin(sSqr);
	if (rbin >= 0) { // test that reflection is in range
	  // Is this in an ice ring? Omit these from averages
	  int Iring = Icerings.InRing(sSqr);
	  while (this_refl.next_observation(this_obs) >= 0) {
	    if (Iring >= 0) {
	      Icerings.AddObs(Iring, this_obs.kI_sigI(), sSqr);
	    } else {
	      sum_I[rbin] += this_obs.kI();
	      sum_sigI[rbin] += this_obs.ksigI();
	      sum_sSqr[rbin] += sSqr;
	      //^
	      //	      if (rbin == 0) { /// inner bin only
	      //		///		fprintf(ndump,"%8.5f %8.1f\n", sSqr, this_obs.ksigI());
	      //		fprintf(ndump,"%10.1f %8.1f %8.5f\n", this_obs.kI(), this_obs.ksigI(), sSqr);
	      //	      } //^-
	      resmnI[rbin].Add(this_obs.kI());
	      resmnSdI[rbin].Add(this_obs.ksigI());
	      n_I[rbin]++;
	      // By run things
	      if (Overall) {
		irun = 0;
		ib = 0;
	      } else {
		irun = this_obs.run();
		ib = this_obs.Batch() - Batch0[irun];
		ASSERT (ib < NBatch[irun]);
		//^		if (ib >= NBatch[irun]) {
		//^		  std::cout << ib << " relative batch\n";
		//^		}
	      }
	      // I, s^2, time
	      sums_st[irun](rbin,ib).add(this_obs.kI(), sSqr, this_obs.time());
	      numobs++;
	      //^
	      //	      std::cout << "sums_st[irun](rbin,ib) " << irun <<" "<< rbin <<" "<< ib <<" "
	      //			<< this_obs.kI() <<"\n"; //^-
	    }
	  }
	}
      } else {
	zonalrefs++; // count reflections in potential zones
      }
    }

    //^    fclose(ndump); //^

    if (numobs == 0) {
      clipper::String msg = "No general reflections accepted in normalisation:\n";
      msg += clipper::String(zonalrefs)+" reflections in potential zero levels (eg h,k,l = 0) not used";
      Message::message(Message_fatal(msg));
    }

    // (<I>) v. sSqr & time
    Normalise NormRes(sums_st);
     
    // Means in resolution bins
    std::vector<MeanIsdIsSqr> meanisdissqr(Nbin);
    //^
    //^    std::cout << "MeanIsdIsSqr " << Nbin << "\n";
    //^-
    for (int i=0;i<Nbin;++i) {
      if (n_I[i] > 0) {
	meanisdissqr[i].MnI = resmnI[i].Mean();
	meanisdissqr[i].SdMnI = resmnI[i].SD();
	meanisdissqr[i].MnSdI = resmnSdI[i].Mean();
	meanisdissqr[i].SdMnSdI = resmnSdI[i].SD();
	meanisdissqr[i].MnsSqr = sum_sSqr[i]/n_I[i];
	meanisdissqr[i].N = n_I[i];
	//^
	//^	std::cout << meanisdissqr[i].MnsSqr
	//^		  << " " << meanisdissqr[i].MnI << " " << meanisdissqr[i].SdMnI 		  <<" "<<meanisdissqr[i].MnSdI<<" "<<meanisdissqr[i].SdMnSdI
	//^		  <<" "<<meanisdissqr[i].N <<"\n";
	//^-
      }
    }
    NormRes.StoreMeanIsdIsSqr(meanisdissqr);

    // Reset high resolution cutoff to cut out weak high resolution bins
    // unless MinIsigRatio < 0
    float smax = ResRange.SResHigh();
    if (MinIsigRatio > 0.0) {
      // Skip 0'th bin in case of low resolution funnies
      int i = Nbin-1;
      for (i=1;i<Nbin;i++) {
	if (n_I[i] > 0) {
	  if (sum_I[i]/sum_sigI[i] < MinIsigRatio) break;
	}
      }
      if (i > 0) i--;
      smax = ResRange.bounds(i).second;
      float highreso = 1./sqrt(smax);
      // Reset low resolution limit
      float lowreso = Max(9.0, ResRange.ResLow());
      if (lowreso < highreso*3.0) lowreso = Min(highreso*3.0, lowreso);
      ResRange.SetRange(lowreso, highreso);
    }

    // Check if ice ring averages are out of line
    for (int Iring=0; Iring<Icerings.Nrings(); Iring++) {
      bool IceReject = false;
      float sSqr = Icerings.MeanSSqr(Iring);
      if (sSqr <= smax) {
	float Imean = Icerings.MeanI(Iring);
	float sigImean = Icerings.MeanSigI(Iring);
	float expectedI = 1./NormRes.applyAvg(1.0, sSqr);
	if ((Imean-expectedI)/sigImean > IceTolerance) {
	  IceReject = true;}
      }
      Icerings.SetReject(Iring, IceReject);
    }
    
    // **** Open files for dumping  ****
    if (PrintLevel > 0) {
      NormDump(NormRes, Batch0, Nbin, sums_st);
    }
  // **** end dump  ****

  return NormRes;
}
//--------------------------------------------------------------
Normalise SetNormaliseMerged(const hkl_merged_list& ref_list,
			     const double& MinIsigRatio,
			     ResoRange& ResRange)
// Set up intensity normalisation object NormRes from merged data
//
// Use binned <I> to get normalisation object
// Optionally, if MinIsigRatio > 0.0
//   return updated ResoRange, setting resolution cutoffs and bins
//
// MinIsigRatio   minimum I/sigI ratio on unaveraged data
{
  int Nbin = ResRange.Nbins();
  std::vector<double> sum_I(Nbin);
  std::vector<double> sum_sigI(Nbin);
  std::vector<double> sum_sSqr(Nbin);
  std::vector<int> n_I(Nbin);

  ref_list.start();

  IsigI Isig;

  int Nruns = 1;
  std::vector<clipper::Array2d<BinSums> > sums_st(Nruns);
  sums_st[0] = clipper::Array2d<BinSums>(Nbin, 1);
  int ntot=0;

  // Loop all reflections
  while (ref_list.next(Isig)) {
    // use only general reflections h!=k!=l!=0 to avoid
    // problems with unknown epsilon
    if ( ! Isig.missing() && ref_list.hkl().IsGeneral()) {
      // Resolution 
      float sSqr = ref_list.invresolsq();
      unsigned int rbin = ResRange.bin(sSqr);
      sum_I[rbin] += Isig.I();
      sum_sigI[rbin] += Isig.sigI();
      sum_sSqr[rbin] += sSqr;
      n_I[rbin]++;
      sums_st[0](rbin,0).add(Isig.I(), sSqr);
      ntot++;
    }
  }

  // (<I>) v. sSqr

  if (ntot == 0) {
    Message::message(Message_fatal("SetNormaliseMerged: no data"));}

    Normalise NormRes(sums_st);
    //^
    //    std::vector<int> Batch0(Nruns,0);
    //    NormDump(NormRes, Batch0, Nbin, sums_st); //^
    //^-

  // Reset high resolution cutoff to cut out weak high resolution bins
  // unless MinIsigRatio < 0
  int i = Nbin;
  if (MinIsigRatio > 0.0) {
    for (i=0;i<Nbin;i++) {
      if (n_I[i] > 0) {
	if (sum_I[i]/sum_sigI[i] < MinIsigRatio) break;
      }
    }
  }
  if (i > 0) i--;
  float smax = ResRange.bounds(i).second;
  float highreso = 1./sqrt(smax);
  // Reset low resolution limit
  float lowreso = Max(9.0, ResRange.ResLow());
  if (lowreso < highreso*3.0) lowreso = Min(highreso*3.0, lowreso);
  ResRange.SetRange(lowreso, highreso);
  return NormRes;
}
}
