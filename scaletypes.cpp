// scaletypes.cpp

#include "scaletypes.hh"
#include "util.hh"
#include "jiffy.hh"
#include "range.hh"
#include "string_util.hh"
#include "restore.hh"

#include <assert.h>
#define ASSERT assert

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;

using phaser_io::itos;
using phaser_io::dtos;
using phaser_io::ftos;


namespace scala {
  //--------------------------------------------------------------
  void ScaleSpecification::dump() const
  {
    if (run < 0) {
    std::cout << "\nScaleSpecification for all runs\n";
    } else {
      std::cout << "\nScaleSpecification for run " << run << "\n";
    }
    if (batch) {
      std::cout << "BATCH mode\n";
      std::cout
        << "nscales " << nscales << "\n"
        << "nbfac " << nbfac << "\n";
    } else {
      std::cout << "ROTATION mode\n";
      if (nscales >= 0) {
        std::cout  << "nscales " << nscales << "\n";
      } else {
        std::cout << "spacing " << spacing << "\n";
      }
      if (nbfac >= 0) {
        std::cout << "nbfac " << nbfac << "\n";
      } else {
        std::cout << "bspacing " << bspacing << "\n";
      }
    }
    if (sec_abs == SecondaryScale::NONE)
      {std::cout << "sec_abs NONE\n";}
    else if (sec_abs == SecondaryScale::SECONDARY) {
      std::cout << "sec_abs SECONDARY\n";
      std::cout    << "lmax " << lmax << " " << lmaxodd << "\n";
    } else if (sec_abs == SecondaryScale::ABSORPTION) {
      std::cout << "sec_abs ABSORPTION\n";
      std::cout    << "lmax " << lmax << " " << lmaxodd << "\n"
                   << "pole " << pole << "\n";
    }
  }
  //--------------------------------------------------------------
  SmoothedValue::SmoothedValue(const Range& Xrange, const int& Ns)
  // Construct from range of raw unnormalised coordinate & number of
  // sample intervals
  // Set smoothing values to defaults, Nav = 3
  //
  // Cases:
  //  1) Ns = 0  single value
  //  2) Ns = 1 or 2,  2 values at z = 0, 1
  //  2) Ns = 3, 3 values at z = 0, 1.0, 2.0
  //  3) Ns > 1, Ns+2 values at -0.5, 0.5, 1.5 ...
  {
    xmin = Xrange.first();
    xmax = Xrange.last();
    // Ok for negative range! xmin may be > xmax, spacing negative
    x0 = xmin;  // coordinate of z = 0
    nsample = Ns;  // number of intervals
    if (nsample == 0) {
      // Special for Ns = 0, single scale, no interpolation
      nvalues = 1;
      spacing = 0.0;
    } else {
      if (nsample == 1) { // one interval
        nvalues = 2;  // if Ns = 1, use 2 values
      } else if (nsample == 2) {
        nvalues = nsample + 1;  // Ns = 2 or 3, use 2 or 3 values
      } else {
        nvalues = nsample + 2;
      }
      //  smoothing spacing
      spacing = (xmax - xmin)/double(nsample);
    }

    values.resize(nvalues);
    positions.resize(nvalues);

    for (int i=0;i<nvalues;++i) {
      values[i] = 0.0;    // just in case
      if (nvalues <= 2) {
        positions[i] = i; // nvalues 1 or 2
      } else if (nvalues == 3) {
        positions[i] = double(i);  // 0.0, 1.0, 2.0
      } else {
        positions[i] = double(i)-0.5; // -0.5, 0.5, 1.5, ...
      }
    }

    //^
    //    std::cout << "Positions: ";
    //    for (int i=0;i<nvalues;++i) {std::cout << " " << positions[i];}
    //    std::cout << "\n";

    const int NAVGDEFAULT = 3;
    SetSmoothing(NAVGDEFAULT, -1.0);  // set default smoothing parameters
  }
  //--------------------------------------------------------------
  void SmoothedValue::SetSmoothing(const int& Naverage, const double& Sigma)
  // Set smoothing values: number of points, sigma
  // If sigma < 0, set to "optimum" (!) (or at least suitable) value from Naverage
  {
    naverage = Naverage;
    sigma = Sigma;

    if (naverage > nvalues) {
      naverage = nvalues;
    }
    if (naverage < 1 || naverage > 5) {
      Message::message(Message_fatal
                       ("SmoothedValue:: Naverage must be between 1 & 5"));
    }
    half_nav = double(naverage)/2.0;
    if (sigma < 0.0) {
      //  Default values 0.65, 0.7, 0.75, 0.8 for nav = 2,3,4,5
        sigma = 0.65 + 0.05 * (naverage-2);
    }
    //^
    //^    std::cout << "Nav, sigma " << naverage << " " << sigma << "\n";
  }
  //--------------------------------------------------------------
  // Store values (nvalues points)
  void SmoothedValue::StoreValues(const std::vector<double>& Values)
  {
    ASSERT (int(Values.size()) == nvalues);
    values = Values;
  }
  //--------------------------------------------------------------
  // Store values (nvalues points)
  void SmoothedValue::StoreValues(const std::vector<float>& Values)
  {
    ASSERT (int(Values.size()) == nvalues);
    for (int i=0;i<nvalues;++i) {
      values[i] = Values[i];
    }
  }
  //--------------------------------------------------------------
  // Store values (nvalues points) to same value
  void SmoothedValue::StoreValue(const double& Value)
  {
    for (int i=0;i<nvalues;++i) {
      values[i] = Value;
    }
  }
  //--------------------------------------------------------------
  double SmoothedValue::Value(const double& x) const
  // Return interpolated value at point (original unnormalised coordinate)
  {
    double sumw;
    std::vector<double> w(nvalues);
    return ValueWeight(x, w, sumw);  // discard weight vector
  }
  //--------------------------------------------------------------
  double SmoothedValue::ValueWeight(const double& x,
                                   std::vector<double>& weight, double& sumweight) const
  // Return interpolated value at point, plus weights at each point,
  // for original unnormalised coordinate
  {
    ASSERT (int(weight.size()) == nvalues);
    double value;
    if (nvalues == 1) {
      // Special for single value
      value = values[0];
      weight[0] = 1.0;
      sumweight = 1.0;
    } else {
      // clear weight vector
      for (int i=0;i<nvalues;++i) {weight[i] = 0.0;}
      // Normalised coordinate
      double z = (x-x0)/spacing;
      double sumwv = 0.0;
      sumweight = 0.0;

      int i1, i2;
      if (nvalues <= 3) { // 2 or 3 points, positions are at ends
        i1 = 0;
        i2 = nvalues;
      } else { // > 3 points
        // 1st point in array (index 0) is at position -0.5
        // Reduce number of points at ends, but not less than 2
        i1 = Nint(z - half_nav)+1;
        i2 = i1 + naverage;
        if (i1 < 0) {
          // beginning of range
          i1 = 0;
          i2 = Max(2,i2);
        }
        if (i2 > nvalues) {
          i2 = nvalues;
          i1 = Min(i1, nvalues-2);
        }
      }
      //^
      //^      std::cout << "z, i1, i2 sigma " << z << " " << i1 << " " << i2
      //^               << " " << sigma << "\n";
      for (int i=i1;i<i2;++i) {
        double ds = (z - positions[i])/sigma;
        weight[i] = exp(-ds*ds);
        sumwv += weight[i] * values[i];
        sumweight  += weight[i];
        //^
        //^     std::cout << "    i, ds, weight, values[i] "
        //^               << i << " " << ds << " " << weight[i] << " " << values[i] << "\n";
      }
      if (sumweight > 0.0) {
        value = sumwv/sumweight;
      } else {
        value = 0.0;
      }
      //^      std::cout << "  Value " << value << "\n";
    }
    return value;
  }
  //--------------------------------------------------------------
  std::string SmoothedValue::FormatSave() const
  // format for save/restore
  {
    std::string dump = "SmoothedValue V1 {\n"; // with version number
    dump += "Nvalues "+itos(nvalues)+"\n";
    dump += "Xmin "+clipper::String(xmin)+"\n";
    dump += "Xmax "+clipper::String(xmax)+"\n";
    dump += "X0 "+clipper::String(x0)+"\n";
    dump += "Spacing "+clipper::String(spacing)+"\n";
    dump += "Nsample "+itos(nsample)+"\n";
    ASSERT (int(values.size()) == nvalues);
    ASSERT (int(positions.size()) == nvalues);
    dump += "Values\n"+StringUtil::FormatSaveVector(values);
    dump += "Positions\n"+StringUtil::FormatSaveVector(positions);
    dump += "Naverage "+itos(naverage)+"\n";
    dump += "Sigma "+clipper::String(sigma)+"\n";;
    return dump+"}\n";
  }
  //--------------------------------------------------------------
  void SmoothedValue::Restore(Fileread& FR)
  // Restore values from file
  {
    FR.ReadTag("SmoothedValue"); // fails if tag does not match
    if (FR.GetTag() != "V1") {  // version check
      clipper::Message::message(Message_fatal
        ("SmoothedValue::Restore incompatible version in "+FR.Filename()));
    }
    FR.Skip();
    FR.ReadTag("Nvalues"); nvalues = FR.Int();
    FR.ReadTag("Xmin"); xmin = FR.Double();
    FR.ReadTag("Xmax"); xmax = FR.Double();
    FR.ReadTag("X0"); x0 = FR.Double();
    FR.ReadTag("Spacing"); spacing = FR.Double();
    FR.ReadTag("Nsample"); nsample = FR.Int();
    FR.ReadTag("Values"); values = FR.DoubleVec(nvalues);
    FR.ReadTag("Positions"); positions = FR.DoubleVec(nvalues);
    FR.ReadTag("Naverage"); naverage = FR.Int();
    FR.ReadTag("Sigma"); sigma = FR.Double();
    half_nav = double(naverage)/2.0;
    if (!FR.CheckEnd()) {
      clipper::Message::message(Message_warn
        ("SmoothedValue::Restore unexpected tag "+FR.Tag()));
    }
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  PrimaryScale::PrimaryScale(const int& NscaleIntervals,
                             const Range& phirange)
  // Construct smooth scaling from number of scale intervals
  // Note that actual number nscales will be NscalesIntervals+2
  //   unless NscalesIntervals = 0, in which nscales = 1
  {
    nscaleintervals = Max(0,NscaleIntervals);
    init(phirange);
  }
  //--------------------------------------------------------------
  PrimaryScale::PrimaryScale(const double& scaleSpacing, const Range& phirange)
  // Construct smooth scaling from scale spacing
  //  scalespacing will be adjusted to give an integral number of intervals
  // Always at least 2 scales
  // note that phirange may be descending
  {
    nscaleintervals = Max(1,Nint(phirange.AbsRange()/scaleSpacing));
    init(phirange);
  }
  //--------------------------------------------------------------
  void PrimaryScale::init(const Range& phirange)
  {
    batchscale = false;
    phi0 = phirange.first();  //starting phi
    // Set up smoothed scales object
    smoothscale = SmoothedValue(phirange, nscaleintervals);
    // Get actual number of scales parameters
    //  may be = 1 for single value
    nscales = smoothscale.Nvalues();
    // Get revised scale interval
    scalespacing = smoothscale.Spacing();
    // Set averaging parameters
    smoothscale.SetSmoothing(navgscale, sdwz);
    // Default to 1.0
    smoothscale.StoreValue(1.0);
    allbatches = true;    // use all batches by default
    nobsPar.assign(nscales, 0);
  }
  //--------------------------------------------------------------
  // Construct batch scaling from list of batch numbers (from run)
  PrimaryScale::PrimaryScale(const std::vector<int>& BatchNumbers)
  {
    batchscale = true;
    nscales = BatchNumbers.size();
    batchscales.assign(nscales, 1.0);
    nscaleintervals = nscales;
    allbatches = true;    // use all batches by default
    // Make lookup table (hash table)
    // Setup up hash lookup table a bit larger than required
    batch_lookup.set_size( int(1.2 * nscales));
    for (int i = 0; i < nscales; i++)  {
      batch_lookup.add(BatchNumbers[i], i);
    }
  }
  //--------------------------------------------------------------
  std::vector<Tie> PrimaryScale::Ties(const double& sdtie, const int& idx0)
  // Return list of ties: sdtie is sd for weight, idx0 is index to first global
  // parameter for setting ties, since they refer to the global parameter index
  {
    std::vector<Tie> ties;
    if (sdtie > 0.0) {
      double weight = 1./(sdtie*sdtie);
      // parameters are tied together in pairs
      // loop from 2nd scale values, batch or smoothed
      if (nscales > 1) {
        for (int i=1;i<nscales;++i) {
          ties.push_back(Tie(idx0+i-1 ,idx0+i, weight));
        }
      }
    }
    return ties;
  }
  //--------------------------------------------------------------
  // valid scale at batch number batchnum
  bool PrimaryScale::ValidScale(const int& batchnum) const
  {
    // Batch serial number
    int jb = batch_lookup.lookup(batchnum);
    if (!allbatches) {
      jb = batchscaleindex[jb];
    }
    if (jb >= 0) {return true;}
    return false;
  }
  //--------------------------------------------------------------
  int PrimaryScale::batchSerialIndex(const int& batchnum) const
  // index number in this batch list for batch number batchnum, omitting rejects
  {
    // Batch serial number
    int jb = batch_lookup.lookup(batchnum);
    if (!allbatches) {
      jb = batchscaleindex[jb];
    }
    return jb;
  }
  //--------------------------------------------------------------
  double PrimaryScale::Scale(const double& phi) const
  // Scale at position phi
  {
    std::vector<double> dgdp;
    double scale;
    ScaleDeriv(false, phi, scale, dgdp);
    return scale;
  }
  //--------------------------------------------------------------
  void PrimaryScale::ScaleDeriv(const double& phi,
                               double& scale, std::vector<double>& dgdp) const
  {
    ScaleDeriv(true, phi, scale, dgdp);
  }
  //--------------------------------------------------------------
  double PrimaryScale::ScaleDeriv(const double& phi,
                                 std::vector<double>& dgdp) const
  {
    double scale;
    ScaleDeriv(true, phi, scale, dgdp);
    return scale;
  }
  //--------------------------------------------------------------
  void PrimaryScale::ScaleDeriv(const bool& Deriv, const double& phi,
                               double& scale, std::vector<double>& dgdp) const
  // Return scale and derivative vector at position phi, smooth scaling
  {
    if (Deriv) {
      // get scale & weight vector
      dgdp.assign(nscales,0.0);
      std::vector<double> w(nscales);
      double sumw;
      scale = smoothscale.ValueWeight(phi, w, sumw);
      for (int i=0;i<nscales;i++) {
        if (sumw > 0.0) {
          dgdp[i] = w[i]/sumw;
        }
      }
    } else {
      // No derivative, just the scale
      scale = smoothscale.Value(phi);
    }
  }
  //--------------------------------------------------------------
  double PrimaryScale::Scale(const int& batch) const
  // Scale at batch number batch
  {
    std::vector<double> dgdp;
    double scale;
    ScaleDeriv(false, batch, scale, dgdp);
    return scale;
  }
  //--------------------------------------------------------------
  void PrimaryScale::ScaleDeriv(const int& batch,
                               double& scale, std::vector<double>& dgdp) const
  {
    ScaleDeriv(true, batch, scale, dgdp);
  }
  //--------------------------------------------------------------
  double PrimaryScale::ScaleDeriv(const int& batch,
                                 std::vector<double>& dgdp) const
  {
    double scale;
    ScaleDeriv(true, batch, scale, dgdp);
    return scale;
  }
  //--------------------------------------------------------------
  void PrimaryScale::ScaleDeriv(const bool& Deriv, const int& batch,
                               double& scale, std::vector<double>& dgdp) const
  // Return scale and derivative vector for batch number batch
  {
    // Batch serial number
    int jb = batch_lookup.lookup(batch);
    if (!allbatches) {
      //      std::cout <<"ScaleDeriv "<<batch<<" "<<jb<<" "<<batchscaleindex.at(jb)<<"\n";
      jb = batchscaleindex[jb];
    }
    ASSERT (jb < nscales);
    if (jb < 0) {
      if (Deriv) {
        dgdp.assign(nscales,0.0);
      }
      scale =  1.0;
    } else {
      scale = batchscales[jb];
      if (Deriv) {
        dgdp.assign(nscales,0.0);
        dgdp[jb] = 1.0;
      }
    }
  }
  //--------------------------------------------------------------
  void PrimaryScale::StoreScales(const std::vector<double>& Scales)
  // Store scales vector
  {
    ASSERT (int(Scales.size()) == nscales);
    if (batchscale) {
      if (allbatches) {
        batchscales = Scales;
      } else {
        for (size_t i=0; i<Scales.size(); i++) {
          batchscales.at(scalebatchindex[i]) = Scales[i];
          //      std::cout <<i<<" "<<scalebatchindex[i]<<" scdx bidx (StoreScales)\n";
        }
      }
    } else {
      smoothscale.StoreValues(Scales);
    }
  }
  //--------------------------------------------------------------
  void PrimaryScale::StoreNobservations(const std::vector<int>& Nobs)
  // Store Nobs vector (length nscales)
  {
    ASSERT (int(Nobs.size()) == nscales);
    if (allbatches) {
      nobsPar = Nobs;
    } else {
      for (size_t i=0; i<Nobs.size(); i++) {
          nobsPar.at(scalebatchindex[i]) = Nobs[i];
      }
    }
  }
  //--------------------------------------------------------------
  std::vector<double>  PrimaryScale::Scales() const
  // Retrieve scales
  {
    if (batchscale) {
      if (allbatches) {
        return batchscales;
      } else {
        std::vector<double> scales(nscales);
        for (size_t ib=0; ib<batchscaleindex.size(); ib++) {
          int k = batchscaleindex[ib];
          if (k >= 0) {
            scales.at(k) = batchscales[ib];
            k++;
          }
        }
        return scales;
      }
    } else {
      return smoothscale.Values();
    }
  }
  //--------------------------------------------------------------
  std::vector<int>  PrimaryScale::Nobservations() const
  // Retrieve number of observations
  {
    if (batchscale) {
      if (allbatches) {
        return nobsPar;
      } else {
        std::vector<int> nobs(nscales);
        for (size_t ib=0; ib<batchscaleindex.size(); ib++) {
          int k = batchscaleindex[ib];
          if (k >= 0) {
            nobs.at(k) = nobsPar[ib];
            k++;
          }
        }
        return nobs;
      }
    } else {
      return nobsPar;
    }
  }
  //--------------------------------------------------------------
  std::string PrimaryScale::format() const
  {
    std::string text;
    if (batchscale) {
      // Batch scaling
      text = "Batch scaling for batches   "+
        clipper::String(batch_lookup.number(0), 8)+
        " to "+clipper::String(batch_lookup.number(nscales-1), 8);
    } else {
      if (nscales == 1) {
        text = "Single scale factor";
      } else {
        text = "Smooth scaling:   "+clipper::String(nscales,3)+
        " scales at intervals of "+clipper::String(scalespacing,6,4)+
        " over range "+clipper::String(phi0,7,5)+" to "+
        clipper::String(phi0+nscaleintervals*scalespacing, 7,5)+
          " in "+clipper::String(nscaleintervals,3)+" parts";
      }
    }
    return text;
  }
  //--------------------------------------------------------------
  void PrimaryScale::SetNavgScale(const int& NavgScale)
  {
    if (NavgScale < 2 || NavgScale > 5) {
      Message::message(Message_fatal
                       ("PrimaryScale:: NavgScale must be 2 to 5"));
    }
    navgscale = NavgScale;
  }
  //--------------------------------------------------------------
  std::string PrimaryScale::FormatSave() const
  // Format all information into a labelled save format for later restoration
  {
    std::string dump = "PrimaryScale V2 {\n";
    dump += "NscaleIntervals "+itos(nscaleintervals)+"\n";
    dump += "Nscales "+itos(nscales)+"\n";
    dump += "NobsPar\n"+StringUtil::FormatSaveVector(nobsPar);
    if (batchscale) {
      dump += "Batch\n";
      dump += batch_lookup.FormatSave();
      if (allbatches) {
        ASSERT (int(batchscales.size()) == nscales);
        dump += "Allbatches\n";
        dump += "BatchScales\n"+StringUtil::FormatSaveVector(batchscales);
      } else {
        dump += "Somebatches\n";
        dump += "Nbatches "+itos(batchscaleindex.size())+"\n";
        dump += "BatchScales\n"+StringUtil::FormatSaveVector(batchscales);
        dump += "Batchscaleindex\n"+StringUtil::FormatSaveVector(batchscaleindex);
        dump += "Nscaleindex "+itos(scalebatchindex.size())+"\n";
        dump += "Scalebatchindex\n"+StringUtil::FormatSaveVector(scalebatchindex);
      }
    } else { // smooth
      dump += "Smooth\n";
      dump += "Scalespacing "+ ftos(scalespacing)+"\n";;
      dump += "Phi0 "+ftos(phi0)+"\n";
      dump += "Scales\n"+smoothscale.FormatSave();
    }
    return dump+"}\n";
  }
  //--------------------------------------------------------------
  // restore
  void PrimaryScale::Restore(Fileread& FR)
  {
    FR.ReadTag("PrimaryScale"); // fails if tag does not match
    std::string version = FR.GetTag();
    int versionnumber = atoi(version.substr(1).c_str());
    if (version != "V1" && version != "V2") {  // version check
      clipper::Message::message(Message_fatal
        ("PrimaryScale::Restore incompatible version in "+FR.Filename()));
    }
    FR.Skip();
    FR.ReadTag("NscaleIntervals"); nscaleintervals = FR.Int();
    FR.ReadTag("Nscales"); nscales = FR.Int();
    FR.ReadTag("NobsPar"); nobsPar = FR.IntVec(nscales);
    std::string tag = FR.GetTag();
    if (tag == "Batch") {
      batchscale = true;
      batch_lookup.Restore(FR);  // batch lookup table
      if (versionnumber >= 2) {
        std::string tag = FR.GetTag();
        if (tag == "Allbatches") {
          allbatches = true;
          FR.ReadTag("BatchScales"); batchscales = FR.DoubleVec(nscales);
        } else if (tag == "Somebatches") {
          FR.ReadTag("Nbatches"); int nbatches = FR.Int();
          FR.ReadTag("BatchScales"); batchscales = FR.DoubleVec(nscales);
          FR.ReadTag("Batchscaleindex"); batchscaleindex = FR.IntVec(nbatches);
          FR.ReadTag("Nscaleindex"); int nscaleindex = FR.Int();
          FR.ReadTag("Scalebatchindex"); scalebatchindex = FR.IntVec(nscaleindex);
        } else {
          Message::message(Message_fatal
                           ("PrimaryScale::Restore unrecognised tag "+tag+
                            " in "+FR.Filename()));
        }
      } else {
        FR.ReadTag("BatchScales"); batchscales = FR.DoubleVec(nscales);
      }
    } else if (tag == "Smooth") { // smooth
      batchscale = false;
      FR.ReadTag("Scalespacing"); scalespacing = FR.Double();
      FR.ReadTag("Phi0"); phi0 = FR.Double();
      FR.ReadTag("Scales"); smoothscale.Restore(FR);
    } else {
      Message::message(Message_fatal
        ("PrimaryScale::Restore unrecognised tag "+tag+
         " in "+FR.Filename()));
    }
    if (!FR.CheckEnd()) {
      clipper::Message::message(Message_warn
        ("PrimaryScale::Restore unexpected tag "+FR.Tag()));
    }
  }
  //--------------------------------------------------------------
  void PrimaryScale::setBatchReject(const std::vector<bool>& Usebatch,
                                    const std::vector<int>& batchnumbers)
  {
    ASSERT (batchscale);
    ASSERT (allbatches == true);
    bool anyreject = false;
    for (size_t ib=0; ib<Usebatch.size(); ib++) {
      int batchnum = batchnumbers[ib];
      int jb = batch_lookup.lookup(batchnum);
      if (jb >= 0) { // we have this batch in this run
        if (!Usebatch[ib]) {
          anyreject = true;
          break;
        }}
    }
    if (!anyreject) return;  // nothing to do

    ASSERT (Usebatch.size() == batchnumbers.size());
    allbatches = false;
    usebatch = Usebatch; // all batches even if not in this run

    // set index lists
    //    std::vector<int> batchscaleindex; // if !allbatches, index into scale list for this batch
    //    std::vector<int> scalebatchindex; // if !allbatches, index into batch list for this scale
    int nbatches = nscales;
    batchscaleindex.assign(nbatches, -1);  // for each batch in this run
    scalebatchindex.clear();
    int k = 0;
    for (size_t ib=0; ib<batchnumbers.size(); ib++) {
      int batchnum = batchnumbers[ib];
      int jb = batch_lookup.lookup(batchnum);
      if (jb >= 0) { // we have this batch in this run
        if (usebatch[ib]) {
          batchscaleindex[jb] = k;  // jb'th batch uses k'th scale
          scalebatchindex.push_back(jb);  // k'th scale corresponds to the jb'th batch
          //^
          //      std::cout << "Use Batch " <<batchnum<<" serial " <<jb
          //                <<" scale index "<<k<<"\n";
          //^-
          k++;
        } else {
          //      std::cout << "Reject batch " << batchnum <<"\n";
        }
      }
    }
    //    std::cout <<"nscales changed from "<<nscales<<" to "<<scalebatchindex.size()<<"\n"; //^
    nscales = scalebatchindex.size();
    nscaleintervals = nscales;
    batchscales.assign(nbatches,1.0);  //B batch scales
    nobsPar.assign(nbatches, 0);
  }
  //--------------------------------------------------------------
  //  Number of points in moving average, 3, 4 or 5
  int PrimaryScale::navgscale = 3;
  double PrimaryScale::sdwz = -1.0;    // "SD" for weighting
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  RelativeBfactor::RelativeBfactor() // null setting
    : batchbfac(false), nbfacintervals(0) ,nbfac(0), time0(0.0)
  {
    bfactors.clear();
  }
  //--------------------------------------------------------------
  RelativeBfactor::RelativeBfactor(const int& NbfacIntervals,
                                   const Range& timerange)
  // Construct smooth B-factors from number of Bfactor intervals
  // Note that actual number nbfac will be NbfacIntervals+2
  //   unless NbfacIntervals = 0, in which nbfac = 1
  {
    nbfacintervals = Max(0,NbfacIntervals);
    init(timerange);
  }
  //--------------------------------------------------------------
  RelativeBfactor::RelativeBfactor(const double& bfacSpacing,
                                   const Range& timerange)
  // Construct smooth B-factors from spacing
  //  spacing will be adjusted to give an integral number of intervals
  // Always at least 2 B-factors
  {
    nbfacintervals = Max(1, Nint(timerange.AbsRange()/bfacSpacing));
    init(timerange);
  }
  //--------------------------------------------------------------
  void RelativeBfactor::init(const Range& timerange)
  {
    batchbfac = false;
    time0 = timerange.min();  //starting time
    // Set up smoothed Bfactor object
    smoothB = SmoothedValue(timerange, nbfacintervals);
    // Get actual number of Bfactor parameters, may be 1 for single value
    nbfac = smoothB.Nvalues();
    // get revised interval
    bfacspacing = smoothB.Spacing();
    // Set averaging parameters
    smoothB.SetSmoothing(navgbfac, sdwt);
    // Default = 0.0
    smoothB.StoreValue(0.0);
    allbatches = true;    // use all batches by default
    nobsPar.assign(nbfac, 0);
  }
  //--------------------------------------------------------------
  RelativeBfactor::RelativeBfactor(const std::vector<int>& BatchNumbers)
  // Construct batch B-factors from  list of batch numbers (from run)
  {
    batchbfac = true;
    nbfac = BatchNumbers.size();
    bfactors = std::vector<double>(nbfac, 0.0);
    allbatches = true;    // use all batches by default
    // Make lookup table (hash table)
    // Setup up hash lookup table a bit larger than required
    batch_lookup.set_size( int(1.2 * nbfac));
    for (int i = 0; i < nbfac; i++)  {
      batch_lookup.add(BatchNumbers[i], i);
    }
  }
  //--------------------------------------------------------------
  std::vector<Tie> RelativeBfactor::Ties(const double& sdtie, const int& idx0) const
  // Return list of ties: sdtie is sd for weight, idx0 is index to first global
  // parameter for setting ties, since they refer to the global parameter index
  {
    std::vector<Tie> ties;
    if (sdtie > 0.0) {
      double weight = 1./(sdtie*sdtie);
      // parameters are tied togther in pairs
      // loop from 2nd value, batch or smoothed
      if (nbfac > 1) {
        for (int i=1;i<nbfac;++i) {
          ties.push_back(Tie(idx0+i-1 ,idx0+i, weight));
        }
      }
    }
    return ties;
  }
  //--------------------------------------------------------------
  std::vector<Tie> RelativeBfactor::ZeroTies(const double& sdtie, const int& idx0) const
  // Return list of ties: sdtie is sd for weight, idx0 is index to first global
  // parameter for setting ties, since they refer to the global parameter index
  // Tie B-factor to zero
  {
    const double TARGET = 0.0;
    std::vector<Tie> ties;
    if (sdtie > 0.0) {
      double weight = 1./(sdtie*sdtie);
      // parameters are tied to zero
      for (int i=0;i<nbfac;++i) {
        ties.push_back(Tie(idx0+i, TARGET, weight));
      }
    }
    return ties;
  }
  //--------------------------------------------------------------
  // valid B-factor at batch number batchnum
  bool RelativeBfactor::ValidBfactor(const int& batchnum) const
  {
    // Batch serial number
    int jb = batch_lookup.lookup(batchnum);
    if (!allbatches) {
      jb = batchbfacindex[jb];
    }
    if (jb >= 0) {return true;}
    return false;
  }
  //--------------------------------------------------------------
  int RelativeBfactor::batchSerialIndex(const int& batchnum) const
  // index number in this batch list for batch number batchnum, omitting rejects
  {
    // Batch serial number
    int jb = batch_lookup.lookup(batchnum);
    if (!allbatches) {
      jb = batchbfacindex[jb];
    }
    return jb;
  }
  //--------------------------------------------------------------
  double RelativeBfactor::BfactorScale(const double& time, const double& invresolsq) const
  // B-Factor at position time, 4(sin theta/lambda)^2 = invresolsq
  {
    std::vector<double> dgdp;
    double gbfac;
    BfactorScaleDeriv(false, time, invresolsq, gbfac, dgdp);
    return gbfac;
  }
  //--------------------------------------------------------------
  void RelativeBfactor::BfactorScaleDeriv(const double& time, const double& invresolsq,
                                     double& gbfac, std::vector<double>& dgdp) const
  {
    BfactorScaleDeriv(true, time, invresolsq, gbfac, dgdp);
  }
  //--------------------------------------------------------------
  double RelativeBfactor::BfactorScaleDeriv(const double& time, const double& invresolsq,
                                     std::vector<double>& dgdp) const
  {
    double gbfac;
    BfactorScaleDeriv(true, time, invresolsq, gbfac, dgdp);
    return gbfac;
  }
  //--------------------------------------------------------------
  void RelativeBfactor::BfactorScaleDeriv(const bool& Deriv, const double& time,
                                     const double& invresolsq,
                                     double& gbfac, std::vector<double>& dgdp) const
  // Return B-factor scale and derivative vector at position time &
  // 4(sin theta/lambda)^2 = invresolsq
  {
    if (Deriv) {
      // get B-factor & weight vector
      dgdp.assign(nbfac,0.0);
      std::vector<double> w(nbfac);
      double sumw;
      gbfac = exp(0.5 * invresolsq * smoothB.ValueWeight(time, w, sumw));
      for (int i=0;i<nbfac;i++) {
        if (sumw > 0.0) {
          dgdp[i] = 0.5 * invresolsq * gbfac * w[i]/sumw;
        }
      }
    } else {
      // No derivative, just the B-factor scale
      gbfac = exp(0.5 * invresolsq * smoothB.Value(time));
    }
  }
  //--------------------------------------------------------------
  double RelativeBfactor::BfactorValue(const double& time)
  // Bfactor value at time "time"
  {
    if (nbfac <= 0) return 0.0;
    return smoothB.Value(time);
  }
  //--------------------------------------------------------------
  double RelativeBfactor::BfactorScale(const int& batch,
                                 const double& invresolsq) const
  // B-Factor scale at batch number batch
  {
    std::vector<double> dgdp;
    double gbfac;
    BfactorScaleDeriv(false, batch, invresolsq, gbfac, dgdp);
    return gbfac;
  }
  //--------------------------------------------------------------
  void RelativeBfactor::BfactorScaleDeriv(const int& batch, const double& invresolsq,
                               double& gbfac, std::vector<double>& dgdp) const
  {
    BfactorScaleDeriv(true, batch, invresolsq, gbfac, dgdp);
  }
  //--------------------------------------------------------------
  double RelativeBfactor::BfactorScaleDeriv(const int& batch, const double& invresolsq,
                                           std::vector<double>& dgdp) const
  {
    double gbfac;
    BfactorScaleDeriv(true, batch, invresolsq, gbfac, dgdp);
    return gbfac;
  }
  //--------------------------------------------------------------
  void RelativeBfactor::BfactorScaleDeriv(const bool& Deriv, const int& batch,
                                     const double& invresolsq,
                                     double& gbfac, std::vector<double>& dgdp) const
  // Return B-factor and derivative vector for batch number batch
  {
    // Batch serial number
    int jb = batch_lookup.lookup(batch);
    if (!allbatches) {
      jb = batchbfacindex[jb];
    }
    ASSERT (jb >= 0 && jb < nbfac);
    gbfac = exp(0.5 * invresolsq * bfactors[jb]);
    if (Deriv) {
      dgdp = std::vector<double>(nbfac,0.0);
      dgdp[jb] = 0.5 * invresolsq * gbfac;
    }
  }
  //--------------------------------------------------------------
  double RelativeBfactor::BfactorValueB(const  int& batch) const
  {
    if (nbfac == 0) return 0.0;
    int jb = batch_lookup.lookup(batch);
    if (!allbatches) {
      jb = batchbfacindex[jb];
    }
    ASSERT (jb >= 0 && jb < nbfac);
    return bfactors[jb];
  }
  //--------------------------------------------------------------
  void RelativeBfactor::StoreBfactors(const std::vector<double>& Bfacs)
  // Store B-factor vector
  {
    ASSERT (int(Bfacs.size()) == nbfac);
    if (batchbfac) {
      if (allbatches) {
        bfactors = Bfacs;
      } else {
        for (size_t i=0; i<Bfacs.size(); i++) {
          bfactors.at(bfacbatchindex[i]) = Bfacs[i];
        }
      }
    } else {
      smoothB.StoreValues(Bfacs);
    }
  }
  //--------------------------------------------------------------
  void RelativeBfactor::StoreNobservations(const std::vector<int>& Nobs)
  // Store Nobs vector (length nscales)
  {
    ASSERT (int(Nobs.size()) == nbfac);
    if (allbatches) {
      nobsPar = Nobs;
    } else {
      for (size_t i=0; i<Nobs.size(); i++) {
          nobsPar.at(bfacbatchindex[i]) = Nobs[i];
      }
    }
  }
  //--------------------------------------------------------------
  std::vector<double> RelativeBfactor::Bfactors() const
  // Retrieve B-factors
  {
    if (nbfac <= 0) {
      return bfactors;  // should be empty
    } else if (batchbfac) {
      if (allbatches) {
        return bfactors;
      } else {
        std::vector<double> bfacs(nbfac);
        for (size_t ib=0; ib<batchbfacindex.size(); ib++) {
          int k = batchbfacindex[ib];
          if (k >= 0) {
            bfacs.at(k) = bfactors[ib];
            k++;
          }
        }
        return bfacs;
      }
    } else {
      return smoothB.Values();
    }
  }
  //--------------------------------------------------------------
  std::vector<int>  RelativeBfactor::Nobservations() const
  // Retrieve number of observations
  {
    if (batchbfac) {
      if (allbatches) {
        return nobsPar;
      } else {
        std::vector<int> nobs(nbfac);
        for (size_t ib=0; ib<batchbfacindex.size(); ib++) {
          int k = batchbfacindex[ib];
          if (k >= 0) {
            nobs.at(k) = nobsPar[ib];
            k++;
          }
        }
        return nobs;
      }
    } else {
      return nobsPar;
    }
  }
  //--------------------------------------------------------------
  std::string RelativeBfactor::format() const
  {
    std::string text;
    if (batchbfac) {
      // Batch scaling
      text = "Batch B-factors for batches "+
        clipper::String(batch_lookup.number(0), 8)+
        " to "+clipper::String(batch_lookup.number(nbfac-1), 8);
    } else {
      if (nbfac == 0) {
        text = "No B-factors";
      } else if (nbfac == 1) {
        text = "Single relative B-factor";
      } else {
        text = "Smooth B-factors: "+clipper::String(nbfac,3)+
          " scales at intervals of "+clipper::String(bfacspacing,6,4)+
          " over range "+clipper::String(time0,7,5)+" to "+
          clipper::String(time0+nbfacintervals*bfacspacing,7,5)+
          " in "+clipper::String(nbfacintervals,3)+" parts";
      }
    }
    return text;
  }
  //--------------------------------------------------------------
  std::string RelativeBfactor::FormatSave() const
  // Format all information into a labelled save format for later restoration
  {
    std::string dump = "RelativeBfactor V2 {\n";
    dump += "NbfacIntervals "+itos(nbfacintervals)+"\n";
    dump += "Nbfac "+itos(nbfac)+"\n";
    if (nbfac > 0) {
      dump += "NobsPar\n"+StringUtil::FormatSaveVector(nobsPar);
      if (batchbfac) {
        dump += "Batch\n";
        dump += batch_lookup.FormatSave();
        if (allbatches) {
          ASSERT (int(bfactors.size()) == nbfac);
          dump += "Allbatches\n";
          dump += "Bfactors\n"+StringUtil::FormatSaveVector(bfactors);
        } else {
          dump += "Somebatches\n";
          dump += "Nbatches "+itos(batchbfacindex.size())+"\n";
          dump += "Bfactors\n"+StringUtil::FormatSaveVector(bfactors);
          dump += "Batchbfacindex\n"+StringUtil::FormatSaveVector(batchbfacindex);
          dump += "Nbfacindex "+itos(bfacbatchindex.size())+"\n";
          dump += "Bfacbatchindex\n"+StringUtil::FormatSaveVector(bfacbatchindex);
        }
      } else { // smooth
        dump += "Smooth\n";
        dump += "Bfacspacing "+ clipper::String(bfacspacing)+"\n";;
        dump += "Time0 "+clipper::String(time0)+"\n";
        dump += "Bfactors\n"+smoothB.FormatSave();
      }
    }
    return dump+"}\n";
  }
  //--------------------------------------------------------------
  // restore
  void RelativeBfactor::Restore(Fileread& FR)
  {
    FR.ReadTag("RelativeBfactor"); // fails if tag does not match
    std::string version = FR.GetTag();
    int versionnumber = atoi(version.substr(1).c_str());
    if (version != "V1" && version != "V2") {  // version check
      clipper::Message::message(Message_fatal
        ("RelativeBfactor::Restore incompatible version in "+FR.Filename()));
    }
    FR.Skip();
    FR.ReadTag("NbfacIntervals"); nbfacintervals = FR.Int();
    FR.ReadTag("Nbfac"); nbfac = FR.Int();
    if (nbfac > 0) {
      FR.ReadTag("NobsPar"); nobsPar = FR.IntVec(nbfac);
      std::string tag = FR.GetTag();
      if (tag == "Batch") {
        batchbfac = true;
        batch_lookup.Restore(FR);
        if (versionnumber >= 2) {
          std::string tag = FR.GetTag();
          if (tag == "Allbatches") {
            allbatches = true;
            FR.ReadTag("Bfactors"); bfactors = FR.DoubleVec(nbfac);
          } else if (tag == "Somebatches") {
            FR.ReadTag("Nbatches"); int nbatches = FR.Int();
            FR.ReadTag("Bfactors"); bfactors = FR.DoubleVec(nbfac);
            FR.ReadTag("Batchbfacindex"); batchbfacindex = FR.IntVec(nbfac);
            FR.ReadTag("Nbfacindex"); int nbfacindex = FR.Int();
            FR.ReadTag("Bfacbatchindex"); bfacbatchindex = FR.IntVec(nbfacindex);
          } else {
            Message::message(Message_fatal
                             ("RelativeBfactor::Restore unrecognised tag "+tag+
                              " in "+FR.Filename()));
          }
        } else {
          FR.ReadTag("Bfactors"); bfactors = FR.DoubleVec(nbfac);
        }
      } else if (tag == "Smooth") { // smooth
        batchbfac = false;
        FR.ReadTag("Bfacspacing"); bfacspacing = FR.Double();
        FR.ReadTag("Time0"); time0 = FR.Double();
        FR.ReadTag("Bfactors"); smoothB.Restore(FR);
      } else {
        clipper::Message::message(Message_fatal
                                  ("RelativeBfactor::Restore unrecognised tag "+tag+
                                   " in "+FR.Filename()));
      }
    }
    if (!FR.CheckEnd()) {
      clipper::Message::message(Message_warn
        ("RelativeBfactor::Restore unexpected tag "+FR.Tag()));
    }
  }
  //--------------------------------------------------------------
  void RelativeBfactor::setBatchReject(const std::vector<bool>& Usebatch,
                                       const std::vector<int>& batchnumbers)
  {
    ASSERT (batchbfac);
    ASSERT (allbatches == true);
    bool anyreject = false;
    for (size_t ib=0; ib<Usebatch.size(); ib++) {
      int batchnum = batchnumbers[ib];
      int jb = batch_lookup.lookup(batchnum);
      if (jb >= 0) { // we have this batch in this run
        if (!Usebatch[ib]) {
          anyreject = true;
          break;
        }}
    }
    if (!anyreject) return;  // nothing to do

    ASSERT (Usebatch.size() == batchnumbers.size());
    allbatches = false;
    usebatch = Usebatch; // all batches even if not in this run

    // set index lists
    //    std::vector<int> batchbfacindex; // if !allbatches, index into bfac list for this batch
    //    std::vector<int> bfacbatchindex; // if !allbatches, index into batch list for this bfac
    int nbatches = nbfac;
    batchbfacindex.assign(nbatches, -1);  // for each batch in this run
    bfacbatchindex.clear();
    int k = 0;
    for (size_t ib=0; ib<batchnumbers.size(); ib++) {
      int batchnum = batchnumbers[ib];
      int jb = batch_lookup.lookup(batchnum);
      if (jb >= 0) { // we have this batch in this run
        if (usebatch[ib]) {
          batchbfacindex[jb] = k;  // jb'th batch uses k'th bfac
          bfacbatchindex.push_back(jb);  // k'th bfac corresponds to the jb'th batch
          //^
          //      std::cout << "Use Batch " <<batchnum<<" serial " <<jb
          //                <<" bfac index "<<k<<"\n";
          //^-
          k++;
        } else {
          //      std::cout << "Reject batch " << batchnum <<"\n";
        }
      }
    }
    //    std::cout <<"nbfacs changed from "<<nbfacs<<" to "<<bfacbatchindex.size()<<"\n"; //^
    nbfac = bfacbatchindex.size();
    nbfacintervals = nbfac;
    bfactors.assign(nbatches,1.0);  //B batch bfacs
    nobsPar.assign(nbatches, 0);
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  //  Number of points in moving average, 3, 4 or 5
  int RelativeBfactor::navgbfac = 3;
  double RelativeBfactor::sdwt = -1;    // "SD" for weighting
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  SecondaryScale::SecondaryScale(const SecondaryScaleType& secSclType,
                                 const int& lMax, const int& lMaxOdd,
                                 const int& Pole)
    : secscltype(secSclType), lmax(lMax), lmaxodd(lMaxOdd), pole(Pole)
  {
    sphHarmonic = SphericalHarmonic(lMax, lMaxOdd);
    ncoeffs = sphHarmonic.Nterms();  // excluding 00 term
    sphcoefficients.assign(ncoeffs, 0.0);  // clear coefficients list
  }
  //--------------------------------------------------------------
  std::vector<Tie> SecondaryScale::Ties(const double& sdtie, const int& idx0)
  // Return list of ties: sdtie is sd for weight, idx0 is index to first global
  // parameter for setting ties, since they refer to the global parameter index
  {
    std::vector<Tie> ties;
    if (sdtie > 0.0) {
      double weight = 1./(sdtie*sdtie);
      double target = 0.0;
      // parameters are tied to target value of 0.0
      // loop all values
      if (ncoeffs > 0) {
        for (int i=0;i<ncoeffs;++i) {
          ties.push_back(Tie(idx0+i, target, weight));
        }
      }
    }
    return ties;
  }
  //--------------------------------------------------------------
  void SecondaryScale::StoreCoefficients(const std::vector<double>& Sphcoefficients)
  // Store coefficient vector (length ncoeffs)
  {
    ASSERT (int(Sphcoefficients.size()) == ncoeffs);
    sphcoefficients = Sphcoefficients;
  }
  //--------------------------------------------------------------
  void SecondaryScale::StoreNobservations(const std::vector<int>& Nobs)
  // Store Nobs vector (length ncoeffs)
  {
    ASSERT (int(Nobs.size()) == ncoeffs);
    nobsPar = Nobs;
  }
  //--------------------------------------------------------------
  double SecondaryScale::Scale(const double& thetap, const double& phip) const
  // Return scale for secondary beam directions
  // polar angles thetap, phip
  {
    std::vector<double> dgdp;
    return ScaleDeriv(thetap, phip, dgdp);
  }
  //--------------------------------------------------------------
  void SecondaryScale::ScaleDeriv(const double& thetap, const double& phip,
                                  double& scale, std::vector<double>& dgdp) const
  // Return scale & derivatives for secondary beam directions
  // polar angles thetap, phip
  {
    // Spherical harmonic values == dgdp  (no 00 term)
    dgdp = sphHarmonic.Ylm(thetap, phip);
    double sc = 1.0;  // constant term
    // omit Y00 = constant
    for (int i=0;i<ncoeffs;++i) {
      sc += sphcoefficients[i] * dgdp[i];
    }
    scale = sc;
  }
  //--------------------------------------------------------------
  double SecondaryScale::ScaleDeriv(const double & thetap, const double& phip,
                                   std::vector<double>& dgdp) const
  // Return scale & derivatives for secondary beam directions
  //polar angles thetap, phip
  {
    double scale;
    ScaleDeriv(thetap, phip, scale, dgdp);
    return scale;
  }
  //--------------------------------------------------------------
  std::string SecondaryScale::format() const
  {
    std::string text;
    if (secscltype == NONE) {
      text += "No secondary beam correction";
    } else {
      if (secscltype == SECONDARY) {
        text+= "Secondary beam correction in camera frame, lmax = "+
          clipper::String(lmax)+", "+clipper::String(lmaxodd);
      } else if (secscltype == ABSORPTION) {
        text+= "Secondary beam correction in crystal frame, lmax = "+
          clipper::String(lmax)+", "+clipper::String(lmaxodd);
        text += ", pole = "+formatPole(pole);
      }
    }
    return text;
  }
  //--------------------------------------------------------------
  std::string SecondaryScale::formatPole(const int& pole) const
  // Return h, k, l for pole = 1,2,3, else "none"
  {
    std::string s = "automatic";
    if (pole == 1) {s = "h";}
    else if (pole == 2) {s = "k";}
    else if (pole == 3) {s = "l";}
    return s;
  }
  //--------------------------------------------------------------
  std::string SecondaryScale::FormatSave() const
  // Format all information into a labelled save format for later restoration
  {
    std::string dump = "SecondaryScale V1 {\n";
    dump += "Lmax "+itos(lmax)+"\n";
    dump += "LmaxOdd "+itos(lmaxodd)+"\n";
    dump += "Pole "+itos(pole)+"\n";
    dump += "Ncoeffs "+itos(ncoeffs)+"\n";
    dump += "NobsPar\n"+StringUtil::FormatSaveVector(nobsPar);
    dump += "SphCoefficients\n"+StringUtil::FormatSaveVector(sphcoefficients);
    return dump+"}\n";
  }
  //--------------------------------------------------------------
  // restore
  void SecondaryScale::Restore(Fileread& FR)
  {
    FR.ReadTag("SecondaryScale"); // fails if tag does not match
    if (FR.GetTag() != "V1") {  // version check
      clipper::Message::message(Message_fatal
        ("SecondaryScale::Restore incompatible version in "+FR.Filename()));
    }
    FR.Skip();
    FR.ReadTag("Lmax"); int lmax_in = FR.Int();
    if (lmax_in != lmax) {ScaleRestoreError::error("lmax");}
    lmax = lmax_in;
    FR.ReadTag("LmaxOdd"); lmaxodd = FR.Int();
    FR.ReadTag("Pole"); pole = FR.Int();
    FR.ReadTag("Ncoeffs"); int ncoeffs_in = FR.Int();
    if (ncoeffs_in != ncoeffs) {ScaleRestoreError::error("ncoeffs");}
    ncoeffs = ncoeffs_in;
    FR.ReadTag("NobsPar"); nobsPar = FR.IntVec(ncoeffs);
    FR.ReadTag("SphCoefficients"); sphcoefficients = FR.DoubleVec(ncoeffs);
    if (pole == 0) {
      secscltype = SECONDARY;
    } else {
      secscltype = ABSORPTION;
    }
    if (!FR.CheckEnd()) {
      clipper::Message::message(Message_warn
        ("SecondaryScale::Restore unexpected tag "+FR.Tag()));
    }
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  void ScaleSpecification::SetConstant(const int& irun)
  // SCALES CONSTANT
  {
    run = irun;
    batch = false;
    nscales = 1;
    spacing = 0.0;
    nbfac = 0;
    bspacing = 0.0;
    sec_abs = SecondaryScale::NONE;
  }
  //--------------------------------------------------------------
} // namespace scala
