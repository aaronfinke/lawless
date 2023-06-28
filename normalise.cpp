// normalise.cpp
//
//  Simplified version, all runs, all batches
// Conversion factor to convert I to E

#include "normalise.hh"
#include "score_datatypes.hh"
#include "scala_util.hh"
#include "string_util.hh"
#include "report_errors.hh"
#include "selectedobservations.hh"
//#include "file_util.hh"     // debug

namespace scala {
  //--------------------------------------------------------------
  Normalise::Normalise(const hkl_unmerge_list& hkl_list,
                       const double& MinIsigRatio,
                       Rings& Icerings,
		       const clipper::U_aniso_frac& u_aniso_frac,
		       const int PrintLevel,
		       bool final)
  {
    setAniso(u_aniso_frac);
    init( hkl_list, MinIsigRatio, Icerings, PrintLevel, final);
  }
  //--------------------------------------------------------------
  void Normalise::init(const hkl_unmerge_list& hkl_list,
                       const double& MinIsigRatio,
                       Rings& Icerings,
		       const int PrintLevel,
		       bool final)
  // Set up intensity normalisation object
  //
  // Use binned <I> to get normalisation object
  // Binned of resolution (from ResRange) over all batches and runs
  // Optionally output plot to file lnI.plot (if Printlevel > 0)
  //
  // MinIsigRatio   minimum I/sigI ratio on averaged data
  //                ranges beyond this threshold get reset
  //                if < 0, no reset
  {
    // We want more resolution bins than there are normally for analysis
    resorange = ResoRange(hkl_list.RRange().min(), hkl_list.RRange().max(),
                          hkl_list.num_observations());
    nrbin = resorange.Nbins();

    // Sums for overall I/sigI etc by resolution for resolution cut-off
    std::vector<MeanValue> mnSqrv(nrbin);

    std::vector<Median<float> > medI(nrbin);
    std::vector<Median<float> > anisoMedI(nrbin);
    // weightedMean used just for SD(mean) from weights
    std::vector<MeanVariance> weightedMean(nrbin);

    // Ice rings
    Icerings.ClearSums();
    float IceTolerance = 3.0;  // changed from 4 in 1.6.14

    reflection this_refl;
    observation this_obs;
    hkl_list.rewind();  // Just in case
    int numobs = 0; // count contributions
    int zonalrefs = 0;

    imax = -100000.;
    double w = 1.0;

    // Use medians of samples
    int numobsinfile = hkl_list.num_observations();

    while (hkl_list.next_reflection(this_refl) >= 0) {
      // use only general reflections h!=k!=l!=0 to avoid
      // problems with unknown epsilon
      // Can't use test for centric, in case space group is centrosymmetric
      if (this_refl.hkl().IsGeneral()) {
        // Resolution
        float sSqr = this_refl.invresolsq();
        int rbin = resorange.tbin(sSqr);
        if (rbin >= 0) { // test that reflection is in range
          // Is this in an ice ring? Omit these from averages
          int Iring = Icerings.InRing(sSqr);
          if (Iring < 0) {
            SelectedObservations sel(this_refl, -1, ALL);
            IsigI Isigav = sel.Average();
            double avI = Isigav.I();
            double w = 1.0/(Isigav.sigI()*Isigav.sigI());
            weightedMean[rbin].Add(double(Isigav.I()), w);
	    if (useAniso) {
	      // apply u_aniso
	      double anisoscale = exp(u_aniso_frac_scaled.quad_form(this_refl.hkl().real()));
	      double aI = avI * anisoscale;
	      anisoMedI[rbin].add(aI);
	    }
            medI[rbin].add(avI);
            mnSqrv[rbin].Add(sSqr);
            imax = std::max(imax, double(Isigav.I()));
            numobs++;
          }
        }
      } else {
        zonalrefs++; // count reflections in centric zones
      }
    }

    if (numobs == 0) {
      clipper::String msg = "No general reflections accepted in normalisation:\n";
      msg += clipper::String(zonalrefs)+" reflections in potential zero levels (eg h,k,l = 0) not used";
      ReportErrors::printFatalError(msg);
    }

    // Remove last bin if it is empty
    while (weightedMean[nrbin-1].Count() == 0) {
      nrbin--;
      resorange.SetNbins(nrbin);
    }

    // Store <I>, sd<I>, count for each bin,
    // weak ones may be replaced below
    setstores();
    int nneg = 0;  // count negative bins
    for (int i=0;i<nrbin;i++) {
      // sets medianI, mnI, mcount, anisomedI, sdI for each bin
      store(i, mnSqrv[i].Mean(), weightedMean[i], medI[i], anisoMedI[i]);
      if (mnI[i] <= 0.0) {
        nneg++;
      }
      //      std::cout << medianI[i] <<" "<< mnI[i] <<" "
      //		<< mcount[i] <<" "<< anisomedI[i] <<" "<< sdI[i] <<"\n";
    }

    if (nneg > 0) {
      // Fix up negative bins, replace by sdI THIS IS A FUDGE
         for (int i=0;i<nrbin;i++) {
	   if (mnI[i] <= 0.0) {
	     mnI[i] = sdI[i];
	   }
	 }
	 std::string message = "Normalisation: "+StringUtil::itos(nneg,3)+
	   " negative bins have been reset to sd(<I>";	 
	 if (useAniso) {
	   message += ", and anisotropy has been turned off";
	   // turn off anisotropic correction
	   useAniso = false;
	 }
	 if (final) {
	   ReportErrors::printWarning(message, "NegativeNormalisation");
	 }
    }
    imean = MeanValue(mnI).Mean();

    std::vector<double> mnI0 = mnI;

    // fit highresolution part to exponential curve, replace mnI
    if (nneg == 0) {mnI = fitLogCurve(MinIsigRatio);}

    // We want mnI, sdI, mcount for each resolution bin
    // Reset weak high resolution bins unless MinIsigRatio < 0
    // NOT USED NOW, replaced by fitLogCurve
    //if (MinIsigRatio > 0.0) {
      // Weak high resolution bins are unreliable, so (pending a better method)
      // replace <I> by a value extrapolated from the last accepted bin
      // Very crude!!
    //  resetweak(MinIsigRatio);
    //}

    // Set up spline
    std::vector<RPair> sSqrmnI; // extra slots at beginning and end
    sSqrmnI.push_back(RPair(resorange.min(), mnI[0])); // duplicate 1st
    for (int irbin=0;irbin<nrbin;++irbin) {
      if (mcount[irbin] > 0) {
        sSqrmnI.push_back(RPair(mnsSqr[irbin], mnI[irbin]));
      }
    }
    // Duplicate last point
    if (mcount[nrbin-1] > 0) {
      sSqrmnI.push_back(RPair(resorange.max(), mnI[nrbin-1]));
    }
    //    for (size_t i=0; i<sSqrmnI.size(); i++) { 
    //      std::cout << sSqrmnI[i].first << " "<< sSqrmnI[i].second <<" <<\n";
    //    }

    bincorr = Spline(sSqrmnI);  // make spline
    valid = true;

    // **** Open files for dumping  ****
    if (PrintLevel > 0) {
      dump();
    }
    // **** end dump  ****

    return;
  }
  //--------------------------------------------------------------
  void Normalise::setAniso(const clipper::U_aniso_frac& u_aniso_frac)
  {
    // Store u_aniso_frac tensor (scaled)
    // Negated and scaled by twoPi^2 for use in quadratic form,
    // u_aniso_frac may be null
    useAniso = true;
    if (u_aniso_frac.is_null()) {
      u_aniso_frac_scaled =
	clipper::U_aniso_frac(clipper::Mat33sym<double>().null());
      useAniso = false;
    } else {
      u_aniso_frac_scaled = -clipper::Util::twopi2()*u_aniso_frac;
      //std::cout<<"Normalise u_aniso_frac_scaled:\n"<<u_aniso_frac_scaled.format()<<"\n";
    }
  }
  //--------------------------------------------------------------
  void Normalise::setstores()
  {
    mnsSqr.resize(nrbin);
    mnI.resize(nrbin);
    medianI.resize(nrbin);
    sdI.resize(nrbin);
    mcount.resize(nrbin);
    anisomedI.resize(nrbin);
  }
  //--------------------------------------------------------------
  void Normalise::store(const int& ibin, const double& sSqr,
                        const MeanVariance& mnv,
                        Median<float>& medI,
                        Median<float>& anisoMedI)
  // Store:
  //   mnI        from median, mean of trimmed range,
  //               after anisotrpic correction if useAniso true
  //   sdI        from weighted mean
  //   medI       median
  {
    const float TRIMFRAC=0.01;
    // from simulations of exponential distribution, factor to correct for
    // 1% trimming top & bottom
    const float TRIMFACTOR=1.038;
    mnsSqr[ibin]  = sSqr;
    medianI[ibin] = medI.median();
    if (useAniso) {
      mnI[ibin] = anisoMedI.meanofrange(1.0f-TRIMFRAC, TRIMFRAC)*TRIMFACTOR;
      mcount[ibin] = anisoMedI.count();
    } else {
      mnI[ibin] = medI.meanofrange(1.0f-TRIMFRAC, TRIMFRAC)*TRIMFACTOR;
      mcount[ibin] = medI.count();
    }
    anisomedI[ibin] = anisoMedI.median();
    sdI[ibin] = mnv.SDofMeanfromWeights();
  }
  //--------------------------------------------------------------
  double Normalise::iovsig(const int& ibin) const
  // return I/sigI for resolution bin
  {
    double snratio = 0.0;
    if (sdI[ibin] > 0.0) {
      snratio = mnI[ibin]/sdI[ibin];
    }
    return snratio;
  }
  //--------------------------------------------------------------
  std::vector<double> Normalise::fitLogCurve(const double& minIsigRatio)
  // fit exponential curve to data beyond resolution of ~3A
  {

    // Find out where the really weak data start
    int i;
    int kfirst = -1;
    // Skip 0'th bin in case of low resolution funnies
    for (i=1;i<nrbin;i++) {
      if (mcount[i] > 0) {
        double mnIovsig = iovsig(i);
	//std::cout <<"***Iovsig, mnI "<<i<<" "<<mnIovsig<<" "<< mnI[i]<<"\n";
        if (kfirst < 0 && mnIovsig < minIsigRatio) {
          kfirst = i;  // 1st bin below threshold
        }
      }
    }
    // kfirst is first bin below threshold


    // fit data beyond this resolution if there are enough points
    const double RESBEYOND = 3.0;
    const int MINPOINTS = 5;
    double sSqrbeyond = 1.0/(RESBEYOND*RESBEYOND);
    int kbin = resorange.tbin(sSqrbeyond);
    std::vector<double> mnIb(nrbin, 0.0);

    if (nrbin-kbin < MINPOINTS) {    // too few points
      return mnI;  // return unchecged data
    }

    if (kfirst < 0) {
      kfirst = kbin;
    }
    LinearFit linfit;
    for (int i=kbin;i<nrbin;++i) {
      if (mnI[i] > 0.0) {  // omit negatives
	double v = log(mnI[i]);
	linfit.add(mnsSqr[i], v, 1.0);
      }
    }

    RPair sfit = linfit.result();
    for (int i=0;i<nrbin;++i) {
      if (i<kfirst) {
	mnIb[i] = mnI[i];
      } else {
	double lnI = sfit.first * mnsSqr[i] + sfit.second;
	mnIb[i] = exp(lnI);
      }
    }
    return mnIb;
  }
  //--------------------------------------------------------------
  void Normalise::resetweak(const double& minIsigRatio)
  {
    // Weak high resolution bins are unreliable, so (pending a better method)
    // replace <I> by a value extrapolated from the last accepted bin
    // Very crude!!   NOT USED NOW

    if (minIsigRatio < 0.0) return;

    // Skip 0'th bin in case of low resolution funnies
    int i;
    int k = -1;
    for (i=1;i<nrbin;i++) {
      if (mcount[i] > 0) {
        double mnIovsig = iovsig(i);
        if (k < 0 && mnIovsig < minIsigRatio) {
          k = i;  // 1st bin below threshold
        }
      }
    }

    if (k > 0) {
      k--; // last accepted bin
      int k2 = nrbin-1;  // last bin
      if (k < k2) {
        // extrapolate values
        double v1 = mnI[k];  // last "reliable" <I>
        double v2 = 0.6*v1;  // value at end, arbitrary
        double d = (v1-v2)/double(k2-k); // difference/bin
        for (int i=k+1;i<nrbin;++i) {
          double xI = v1 - double(i-k) * d;
          double xsd = sdI[k];  // just propagate sd
          mnI[i] = xI;   // store modified values
          sdI[i] = xsd;
        }
      }
    }
  }
  //--------------------------------------------------------------
  float Normalise::apply(const float& I, const float& sSqr,
			 const DVect3& rhkl) const
  // apply correction to generate E^2
  //  The relevant <I> for hkl is anisoCorr/<I>, so divide by that
  {
    if (!valid) {
      ReportErrors::printFatalError("Normalise not set");
    }
    // Dividing scale, = <I>/anisoscale
    float scorr = bincorr.Interpolate(sSqr);
    // apply u_aniso if useAniso true
    double anisoscale = anisoCorr(rhkl);
    return I * (anisoscale / scorr);
  }
  //--------------------------------------------------------------
  IsigI Normalise::apply(const IsigI& Is, const float& sSqr,
			 const DVect3& rhkl) const
  // apply correction to generate E^2
  //  The relevant <I> for hkl is anisoCorr/<I>, so divide by that
  {
    if (!valid) {
      ReportErrors::printFatalError("Normalise not set");
    }
    // Dividing scale, = <I>
    float scorr = bincorr.Interpolate(sSqr);
    IsigI IsScl(Is);
    // apply u_aniso  if useAniso true
    double anisoscale = anisoCorr(rhkl);
    IsScl.scale(anisoscale / scorr);
    return IsScl;
  }
  //--------------------------------------------------------------
  float Normalise::anisoCorr(const DVect3& rhkl) const
  // anisotropic part of correction
  // This is the scale needed to bring an individual I to match the average
  {
    double anisoscale = 1.0;
    if (useAniso) {
      // apply u_aniso
      anisoscale = exp(u_aniso_frac_scaled.quad_form(rhkl));
    }
    return anisoscale;
  }
  //--------------------------------------------------------------
  float Normalise::Corr(const float& sSqr, const DVect3& rhkl) const
  // correction, multiplying scale
  {
    if (!valid) {
      ReportErrors::printFatalError("Normalise not set");
    }
    float scorr = bincorr.Interpolate(sSqr);
    double anisoscale = 1.0;
    if (std::abs(rhkl*rhkl) > 1.0e-20) {
      anisoscale = anisoCorr(rhkl);
    }
    return anisoscale / scorr;
  }
  //--------------------------------------------------------------
  // mean/median ratio, averaged over some low resolution bins
  double Normalise::mmratio() const
  {
    MeanValue meanratio;
    int totalcount = 0;
    const int MINNUMBER = 1000;

    for (int is=0;is<nrbin;is++) {         // loop resolution bins
      if (mcount[is] > 0) {
        if (totalcount > MINNUMBER) {break;}
        float mnmedratio = mnI[is]/medianI[is];
        meanratio.Add(mnmedratio);
        totalcount += mcount[is];
      }
    }
    return meanratio.Mean();
  }
  //--------------------------------------------------------------
  void Normalise::dump(const std::string& filename)
  {
    // **** Open file for dumping  ****
    FILE* file;
    std::string name = filename;
    if (filename == "") {
      name = "norm.plot";
    }
    file = fopen(name.c_str(), "w");
    if (file == NULL) {
      ReportErrors::printFatalError("Can't open file "+name);
    }
    fprintf(file,
      "   sSqr      mnI    sdI  medianI AnisoMed   Corr    IovSd  mean/median    N\n");

    for (int is=0;is<nrbin;is++) {         // loop resolution bins
      if (mcount[is] > 0) {
        // mean from removing extremes
        float meanI  = mnI[is];
        float mnmedratio = meanI/medianI[is];
        fprintf(file,
                "%8.4f %7.1f %6.2f %7.1f %7.1f %8.4f %8.1f %8.3f  %7d\n",
                mnsSqr[is], mnI[is], sdI[is], 
		medianI[is], anisomedI[is],
                Corr(mnsSqr[is], DVect3(0,0,0)), iovsig(is),
                mnmedratio, mcount[is]);
      }
    }

    int status = fclose(file);
    status = status;
    // **** end dump  ****
  }
}
