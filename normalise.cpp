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
                       const int PrintLevel)
  {
    init( hkl_list, MinIsigRatio, Icerings, PrintLevel);
  }
  //--------------------------------------------------------------
  void Normalise::init(const hkl_unmerge_list& hkl_list,
                       const double& MinIsigRatio,
                       Rings& Icerings,
                       const int PrintLevel)
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
    std::vector<Median<float> > medI2(nrbin);
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
            medI[rbin].add(avI);
            medI2[rbin].add(avI*avI);
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

    // Store <I>, sd<I>, count for each bin,
    // weak ones may be replaced below
    setstores();
    int nneg = 0;  // count negative bins
    for (int i=0;i<nrbin;i++) {
      store(i, mnSqrv[i].Mean(), weightedMean[i], medI[i]);
      if (weightedMean[i].Mean() <= 0.0) {
        nneg++;
      }
    }

    imean = MeanValue(mnI).Mean();

    // We want mnI, sdI, mcount for each resolution bin
    // Reset weak high resolution bins unless MinIsigRatio < 0
    if (MinIsigRatio > 0.0) {
      // Weak high resolution bins are unreliable, so (pending a better method)
      // replace <I> by a value extrapolated from the last accepted bin
      // Very crude!!
      resetweak(MinIsigRatio);
    }

    // Set up spline
    std::vector<RPair> sSqrmnI; // extra slots at beginning and end
    sSqrmnI.push_back(RPair(resorange.min(), mnI[0])); // duplicate 1st
    for (int irbin=0;irbin<nrbin;++irbin) {
      if (mcount[irbin] > 0) {
        sSqrmnI.push_back(RPair(mnsSqr[irbin], mnI[irbin]));
      }
    }
    // Duplicate last point
    sSqrmnI.push_back(RPair(resorange.max(), mnI[nrbin-1]));
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
  void Normalise::setstores()
  {
    mnsSqr.resize(nrbin);
    mnI.resize(nrbin);
    medianI.resize(nrbin);
    sdI.resize(nrbin);
    mcount.resize(nrbin);
  }
  //--------------------------------------------------------------
  void Normalise::store(const int& ibin, const double& sSqr,
                        const MeanVariance& mnv,
                        Median<float>& medI)
  // Store:
  //   mnI      from median, mean of trimmed range
  //   sdI      from weighted mean
  //   medianI  median
  {
    const float TRIMFRAC=0.01;
    // from simulations of exponential distribution, factor to correct for
    // 1% trimming top & bottom
    const float TRIMFACTOR=1.038;
    mnsSqr[ibin]  = sSqr;
    medianI[ibin] = medI.median();
    mnI[ibin] = medI.meanofrange(1.0f-TRIMFRAC, TRIMFRAC)*TRIMFACTOR;
    sdI[ibin] = mnv.SDofMeanfromWeights();
    mcount[ibin] = medI.count();
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
  void Normalise::resetweak(const double& minIsigRatio)
  {
    // Weak high resolution bins are unreliable, so (pending a better method)
    // replace <I> by a value extrapolated from the last accepted bin
    // Very crude!!

    if (minIsigRatio < 0.0) return;

    // Skip 0'th bin in case of low resolution funnies
    int i;
    int k = -1;
    for (i=1;i<nrbin;i++) {
      if (mcount[i] > 0) {
        double mnIovsig = iovsig(i);
        //^     std::cout <<"***Iovsig "<<i<<" "<<mnIovsig<<"\n";
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
        //      std::cout << "Nresetweak "<<k<<" "<<v1<<" "<<sdI[k]<<std::endl; //^
        for (int i=k+1;i<nrbin;++i) {
          double xI = v1 - double(i-k) * d;
          double xsd = sdI[k];  // just propagate sd
          mnI[i] = xI;   // store modified values
          sdI[i] = xsd;
          //          std::cout << "Nresetweak "<<i<<" "<<xI<<" "<<xsd<<std::endl; //^
        }
      }
    }
  }
  //--------------------------------------------------------------
  float Normalise::apply(const float& I, const float& sSqr) const
  // apply correction
  {
    if (!valid) {
      ReportErrors::printFatalError("Normalise not set");
    }
    // Dividing scale, = <I>
    float scorr = bincorr.Interpolate(sSqr);
    return I / scorr;
  }
  //--------------------------------------------------------------
  IsigI Normalise::apply(const IsigI& Is, const float& sSqr) const
  {
    if (!valid) {
      ReportErrors::printFatalError("Normalise not set");
    }
    // Dividing scale, = <I>
    float scorr = bincorr.Interpolate(sSqr);
    IsigI IsScl(Is);
    IsScl.scale(1./scorr);
    return IsScl;
  }
  //--------------------------------------------------------------
  float Normalise::Corr(const float& sSqr) const
  // correction, multiplying scale
  {
    if (!valid) {
      ReportErrors::printFatalError("Normalise not set");
    }
    float scorr = bincorr.Interpolate(sSqr);
    return 1.0 / scorr;
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
      "   sSqr      mnI      sdI   medianI    Corr    IovSd  mean/median     N\n");

    for (int is=0;is<nrbin;is++) {         // loop resolution bins
      if (mcount[is] > 0) {
        // mean from removing extremes
        float meanI  = mnI[is];
        float mnmedratio = meanI/medianI[is];
        fprintf(file,
                "%8.4f %8.1f %8.2f %8.1f %8.4f %8.1f %8.3f  %8d\n",
                mnsSqr[is], mnI[is], sdI[is], medianI[is],
                Corr(mnsSqr[is]), iovsig(is),
                mnmedratio, mcount[is]);
      }
    }

    int status = fclose(file);
    status = status;
    // **** end dump  ****
  }
}
