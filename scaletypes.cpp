// scaletypes.cpp

#include "scaletypes.hh"
#include "util.hh"
#include "jiffy.hh"
#include "range.hh"
#include "string_util.hh"
#include "restore.hh"

#include <assert.h>
#define ASSERT assert

#include <cmath>
#include <limits>
#include <algorithm>
// Eigen3 (header-only) for the GPR Cholesky solve
#include <Eigen/Dense>

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
        text += ", pole = "+formatPole();
      }
    }
    return text;
  }
  //--------------------------------------------------------------
  std::string SecondaryScale::formatPole() const
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
  // =================== WavelengthChebyshevScale ====================
  //--------------------------------------------------------------
  WavelengthChebyshevScale::WavelengthChebyshevScale(
    const std::vector<WavelengthRange>& Ranges, const double& LambdaRef)
    : lambda_ref(LambdaRef), f_ref(0.0)
  {
    ranges = Ranges;
    ncoeffs = 0;
    for (size_t i=0; i<ranges.size(); ++i) {
      ranges[i].offset = ncoeffs;
      ncoeffs += ranges[i].degree + 1;
    }
    // Log parameterisation: h(lambda) = sum_k a_k T_k(z)
    // ws = exp(h(lambda) - h(lambda_ref))
    // Initial: all a_k = 0 => h=0 everywhere => ws=1 (flat spectrum)
    coeffs.assign(ncoeffs, 0.0);
    ref_range = rangeIndex(lambda_ref);
    f_ref = 0.0;  // h(lambda_ref) = 0 initially
  }
  //--------------------------------------------------------------
  void WavelengthChebyshevScale::StoreCoefficients(const std::vector<double>& Coeffs)
  {
    if (int(Coeffs.size()) != ncoeffs) {
      Message::message(Message_fatal
        ("WavelengthChebyshevScale::StoreCoefficients: wrong number of coefficients"));
    }
    coeffs = Coeffs;
    // recompute log_ref = h(lambda_ref) for log parameterisation
    if (ref_range >= 0) {
      f_ref = evalRange(lambda_ref, ref_range);  // f_ref stores h(lambda_ref) = log-poly at ref
    } else {
      f_ref = 0.0;
    }
  }
  //--------------------------------------------------------------
  double WavelengthChebyshevScale::chebeval(const std::vector<double>& c,
                                             const double& z) const
  // Clenshaw's algorithm: sum_{k=0}^{n} c[k] T_k(z)
  {
    int n = int(c.size()) - 1;
    if (n < 0) return 0.0;
    if (n == 0) return c[0];
    double D = 0.0, DD = 0.0;
    for (int j = n; j >= 1; --j) {
      double sv = D;
      D = 2.0*z*D - DD + c[j];
      DD = sv;
    }
    return z*D - DD + c[0];
  }
  //--------------------------------------------------------------
  void WavelengthChebyshevScale::chebbasis(const int& degree, const double& z,
                                            std::vector<double>& T) const
  // Chebyshev basis: T[k] = T_k(z) for k=0..degree
  {
    T.resize(degree+1);
    T[0] = 1.0;
    if (degree >= 1) T[1] = z;
    for (int k=2; k<=degree; ++k) {
      T[k] = 2.0*z*T[k-1] - T[k-2];
    }
  }
  //--------------------------------------------------------------
  int WavelengthChebyshevScale::rangeIndex(const double& lambda) const
  // Return range index for lambda, -1 if outside all ranges
  {
    for (size_t i=0; i<ranges.size(); ++i) {
      if (lambda >= ranges[i].lam_min && lambda <= ranges[i].lam_max) {
        return int(i);
      }
    }
    return -1;
  }
  //--------------------------------------------------------------
  double WavelengthChebyshevScale::mapToZ(const double& lambda,
                                           const int& irange) const
  {
    return (2.0*lambda - ranges[irange].lam_min - ranges[irange].lam_max)
           / (ranges[irange].lam_max - ranges[irange].lam_min);
  }
  //--------------------------------------------------------------
  double WavelengthChebyshevScale::evalRange(const double& lambda,
                                              const int& irange) const
  {
    double z = mapToZ(lambda, irange);
    int off = ranges[irange].offset;
    int deg = ranges[irange].degree;
    std::vector<double> c(coeffs.begin()+off, coeffs.begin()+off+deg+1);
    return chebeval(c, z);
  }
  //--------------------------------------------------------------
  double WavelengthChebyshevScale::Scale(const double& lambda) const
  // Log parameterisation: coefficients represent h(lambda) = log(f(lambda))
  // ws(lambda) = exp(h(lambda) - h(lambda_ref))
  {
    if (ncoeffs == 0) return 1.0;
    int ir = rangeIndex(lambda);
    if (ir < 0) return 1.0;  // out of range: no correction
    double h = evalRange(lambda, ir);   // h(lambda) = log-polynomial value
    return std::exp(h - f_ref);         // f_ref holds h(lambda_ref)
  }
  //--------------------------------------------------------------
  double WavelengthChebyshevScale::ScaleDeriv(const double& lambda,
                                               std::vector<double>& dgdp) const
  // Log parameterisation: h(lambda) = sum_k a_k T_k(z)
  // ws = exp(h(lambda) - h(lambda_ref))
  // d(ws)/d(a_k) = ws * (T_k(z) - T_k(z_ref))
  {
    dgdp.assign(ncoeffs, 0.0);
    if (ncoeffs == 0) return 1.0;

    int ir = rangeIndex(lambda);
    if (ir < 0) return 1.0;

    double z = mapToZ(lambda, ir);
    int off = ranges[ir].offset;
    int deg = ranges[ir].degree;

    double h = evalRange(lambda, ir);   // h(lambda)
    double w = std::exp(h - f_ref);     // ws = exp(h - h_ref)

    // Chebyshev basis at z
    std::vector<double> Tz;
    chebbasis(deg, z, Tz);

    // Chebyshev basis at z_ref (same range)
    std::vector<double> Tzref;
    if (ref_range == ir) {
      double zref = mapToZ(lambda_ref, ir);
      chebbasis(deg, zref, Tzref);
    } else {
      Tzref.assign(deg+1, 0.0);
    }

    for (int k=0; k<=deg; ++k) {
      dgdp[off+k] = w * (Tz[k] - Tzref[k]);
    }
    return w;
  }
  //--------------------------------------------------------------
  std::string WavelengthChebyshevScale::PrintNormalization(const int& npoints) const
  {
    if (!IsActive()) return "";
    std::string s = "\n Wavelength normalization (Chebyshev)\n";
    s += FormatOutput::logTabPrintf(1,
           " Reference wavelength: %6.4f A\n", lambda_ref);
    for (int ir = 0; ir < int(ranges.size()); ++ir) {
      const WavelengthRange& r = ranges[ir];
      s += FormatOutput::logTabPrintf(1,
             "\n Range %d: %6.4f - %6.4f A   degree %d\n",
             ir+1, r.lam_min, r.lam_max, r.degree);
      s += FormatOutput::logTabPrintf(1, "   Log-coefficients:");
      for (int k = 0; k <= r.degree; ++k) {
        s += FormatOutput::logTabPrintf(0, " %9.5f", coeffs[r.offset + k]);
      }
      s += "\n";
      s += FormatOutput::logTabPrintf(1,
             "   %8s  %10s\n", "lambda", "w(lambda)");
      double step = (r.lam_max - r.lam_min) / (npoints - 1);
      for (int ip = 0; ip < npoints; ++ip) {
        double lam = r.lam_min + ip * step;
        double w = Scale(lam);
        std::string marker = (std::abs(lam - lambda_ref) < 0.5 * step) ? " <- ref" : "";
        s += FormatOutput::logTabPrintf(1,
               "   %8.4f  %10.5f%s\n", lam, w, marker.c_str());
      }
      s += AsciiPlot(ir);
    }
    return s;
  }
  //--------------------------------------------------------------
  std::string WavelengthChebyshevScale::AsciiPlot(const int& irange,
                                                  const int& width,
                                                  const int& height) const
  // ASCII line plot of w(lambda) over one range.  The '*' trace shows the
  // normalization curve; the reference wavelength is marked with a ':' column.
  {
    const WavelengthRange& r = ranges[irange];
    if (r.lam_max <= r.lam_min || width < 2 || height < 2) return "";

    // Sample w(lambda) across the range
    std::vector<double> w(width);
    double wmin = 1.0e30, wmax = -1.0e30;
    for (int col = 0; col < width; ++col) {
      double lam = r.lam_min + (r.lam_max - r.lam_min) * col / double(width - 1);
      w[col] = Scale(lam);
      wmin = std::min(wmin, w[col]);
      wmax = std::max(wmax, w[col]);
    }
    if (wmax - wmin < 1.0e-6) {  // flat curve: pad so it sits mid-plot
      double mid = 0.5 * (wmin + wmax);
      wmin = mid - 0.5;
      wmax = mid + 0.5;
    }

    // Column of the reference wavelength, -1 if outside this range
    int refcol = -1;
    if (lambda_ref >= r.lam_min && lambda_ref <= r.lam_max) {
      refcol = int((lambda_ref - r.lam_min) / (r.lam_max - r.lam_min)
                   * (width - 1) + 0.5);
    }

    // Build grid (height rows top=wmax, bottom=wmin)
    std::vector<std::string> grid(height, std::string(width, ' '));
    if (refcol >= 0) {
      for (int row = 0; row < height; ++row) grid[row][refcol] = ':';
    }
    for (int col = 0; col < width; ++col) {
      int row = int((wmax - w[col]) / (wmax - wmin) * (height - 1) + 0.5);
      if (row < 0) row = 0;
      if (row >= height) row = height - 1;
      grid[row][col] = '*';
    }

    std::string s = FormatOutput::logTabPrintf(1, "   w(lambda):\n");
    for (int row = 0; row < height; ++row) {
      double wlabel = wmax - (wmax - wmin) * row / double(height - 1);
      s += FormatOutput::logTabPrintf(1, "   %8.4f |%s\n",
                                      wlabel, grid[row].c_str());
    }
    // x-axis
    s += FormatOutput::logTabPrintf(1, "            +%s\n",
                                    std::string(width, '-').c_str());
    std::string lo = clipper::String(r.lam_min, 6, 3).trim();
    std::string hi = clipper::String(r.lam_max, 6, 3).trim();
    std::string axis(width, ' ');
    for (size_t i = 0; i < lo.size() && i < axis.size(); ++i) axis[i] = lo[i];
    int hipos = width - int(hi.size());
    if (hipos < 0) hipos = 0;
    for (size_t i = 0; hipos + i < axis.size() && i < hi.size(); ++i)
      axis[hipos + i] = hi[i];
    s += FormatOutput::logTabPrintf(1, "            %s  lambda (A)\n",
                                    axis.c_str());
    if (refcol >= 0) {
      s += FormatOutput::logTabPrintf(1,
             "            (':' marks reference wavelength %6.4f A)\n",
             lambda_ref);
    }
    return s;
  }
  //--------------------------------------------------------------
  std::string WavelengthChebyshevScale::asXML() const
  {
    if (!IsActive()) return "";
    std::string s = "<WavelengthNormalisation>\n";
    s += "  " + StringUtil::MakeXMLtag("ReferenceWavelength", lambda_ref, 8, 4);
    s += "\n";
    const int npoints = 21;
    for (int ir = 0; ir < int(ranges.size()); ++ir) {
      const WavelengthRange& r = ranges[ir];
      s += "  <Range number=\"" + clipper::String(ir + 1).trim() + "\">\n";
      s += "    " + StringUtil::MakeXMLtag("LambdaMin", r.lam_min, 8, 4);
      s += StringUtil::MakeXMLtag("LambdaMax", r.lam_max, 8, 4);
      s += StringUtil::MakeXMLtag("Degree", r.degree, 4) + "\n";
      std::string coefstr;
      for (int k = 0; k <= r.degree; ++k) {
        coefstr += clipper::String(coeffs[r.offset + k], 12, 6);
      }
      s += "    " + StringUtil::MakeXMLtag("LogCoefficients", coefstr) + "\n";
      s += "    <Normalisation>\n";
      double step = (r.lam_max - r.lam_min) / (npoints - 1);
      for (int ip = 0; ip < npoints; ++ip) {
        double lam = r.lam_min + ip * step;
        s += "      <point>";
        s += StringUtil::MakeXMLtag("lambda", lam, 8, 4);
        s += StringUtil::MakeXMLtag("w", Scale(lam), 10, 5);
        s += " </point>\n";
      }
      s += "    </Normalisation>\n";
      s += "  </Range>\n";
    }
    s += "</WavelengthNormalisation>\n";
    return s;
  }
  //--------------------------------------------------------------
  std::string WavelengthChebyshevScale::GnuplotScript(const std::string& title,
                                                      const std::string& version) const
  // Self-contained gnuplot script (LAMBDANORM): w(lambda), one segment per range.
  {
    if (!IsActive()) return "";
    // sanitise title for use inside a double-quoted gnuplot string
    std::string t = title;
    for (size_t i = 0; i < t.size(); ++i) {
      if (t[i] == '"') t[i] = '\'';
      if (t[i] == '\n' || t[i] == '\r') t[i] = ' ';
    }
    t = clipper::String(t).trim();
    std::string lref = clipper::String(lambda_ref, 8, 4).trim();

    const int nr = int(ranges.size());
    const int npoints = 100;
    // distinct per-range colours (cycled if more than 5 ranges, but MAXRANGES=5)
    const char* col[5] = {"#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e"};

    // overall plotted wavelength domain = union of all ranges
    double gmin = ranges[0].lam_min, gmax = ranges[0].lam_max;
    for (int i = 0; i < nr; ++i) {
      gmin = std::min(gmin, ranges[i].lam_min);
      gmax = std::max(gmax, ranges[i].lam_max);
    }
    // per-range scale using THAT range's polynomial (extrapolated outside it):
    //   w_ir(lam) = exp(h_ir(lam) - h(lambda_ref))
    // y-range is taken from the USED (in-interval) portions only, so divergent
    // extrapolations simply run off-screen and are clipped.
    double ymin = 1.0e30, ymax = -1.0e30;
    for (int ir = 0; ir < nr; ++ir) {
      const WavelengthRange& r = ranges[ir];
      for (int ip = 0; ip < npoints; ++ip) {
        double lam = r.lam_min + (r.lam_max - r.lam_min) * ip / double(npoints - 1);
        double w = std::exp(evalRange(lam, ir) - f_ref);
        if (std::isfinite(w)) { ymin = std::min(ymin, w); ymax = std::max(ymax, w); }
      }
    }
    if (!(ymax > ymin)) { ymin = 0.0; ymax = 1.0; }
    double ymargin = 0.12 * (ymax - ymin);
    double ylo = ymin - ymargin; if (ylo < 0.0) ylo = 0.0;
    double yhi = ymax + ymargin;

    std::string s;
    s += "# LAMBDANORM - Laue wavelength normalization curve\n";
    s += "# Generated by " + version + "\n";
    if (!t.empty()) s += "# Title: " + t + "\n";
    s += "#\n";
    s += "# View this plot with gnuplot:\n";
    s += "#     gnuplot -p LAMBDANORM\n";
    s += "# (the -p flag keeps the plot window open after drawing)\n";
    s += "#\n";
    s += "# Method: Chebyshev wavelength normalization, fitted in log space.\n";
    s += "# w(lambda) = exp(h(lambda) - h(lambda_ref)); reference lambda = "
         + lref + " A.\n";
    s += "# Each range is drawn in full in its own colour: SOLID where the range\n";
    s += "# is used (within its own interval), DASHED where the range's polynomial\n";
    s += "# is extrapolated outside its interval (not used for scaling).\n";
    s += "# Data columns (per datablock):\n";
    s += "#   1 lambda (A)\n";
    s += "#   2 w(lambda)   normalization scale for that range\n";
    s += "\n";
    s += "set title \"Laue wavelength normalization (Chebyshev)";
    if (!t.empty()) s += "\\n" + t;
    s += "  -  " + version + "\"\n";
    s += "set xlabel \"Wavelength {/Symbol l} / Angstrom\"\n";
    s += "set ylabel \"Normalization w({/Symbol l})\"\n";
    s += "set grid\n";
    s += "set key top right\n";
    s += "set xrange [" + clipper::String(gmin, 8, 4).trim() + ":"
         + clipper::String(gmax, 8, 4).trim() + "]\n";
    s += "set yrange [" + clipper::String(ylo, 8, 4).trim() + ":"
         + clipper::String(yhi, 8, 4).trim() + "]\n";
    if (lambda_ref > 0.0 && rangeIndex(lambda_ref) >= 0) {
      s += "set arrow from " + lref + ", graph 0 to " + lref
           + ", graph 1 nohead dashtype 2 lc rgb \"#888888\"\n";
      s += "set label \"ref\" at " + lref + ", graph 0.96 center tc rgb \"#888888\"\n";
    }
    s += "\n";

    // one solid datablock (used interval) and one dashed datablock (extrapolation
    // outside the interval, on either side) per range.  Boundary lambdas are
    // shared so solid and dashed segments meet.  Non-finite extrapolated points
    // are skipped (the dashed line just stops).
    for (int ir = 0; ir < nr; ++ir) {
      const WavelengthRange& r = ranges[ir];
      // solid: used interval
      s += "$LNs" + clipper::String(ir).trim() + " << EOD\n";
      for (int ip = 0; ip < npoints; ++ip) {
        double lam = r.lam_min + (r.lam_max - r.lam_min) * ip / double(npoints - 1);
        s += clipper::String(lam, 10, 4) + " "
           + clipper::String(std::exp(evalRange(lam, ir) - f_ref), 12, 5) + "\n";
      }
      s += "EOD\n";
      // dashed: extrapolation left of and right of the used interval
      s += "$LNd" + clipper::String(ir).trim() + " << EOD\n";
      bool wrote = false;
      if (r.lam_min > gmin) {   // left extrapolation, ending at r.lam_min
        for (int ip = 0; ip < npoints; ++ip) {
          double lam = gmin + (r.lam_min - gmin) * ip / double(npoints - 1);
          double w = std::exp(evalRange(lam, ir) - f_ref);
          if (std::isfinite(w)) {
            s += clipper::String(lam, 10, 4) + " " + clipper::String(w, 12, 5) + "\n";
          }
        }
        wrote = true;
      }
      if (r.lam_max < gmax) {   // right extrapolation, starting at r.lam_max
        if (wrote) s += "\n";   // pen-up between the two extrapolated arms
        for (int ip = 0; ip < npoints; ++ip) {
          double lam = r.lam_max + (gmax - r.lam_max) * ip / double(npoints - 1);
          double w = std::exp(evalRange(lam, ir) - f_ref);
          if (std::isfinite(w)) {
            s += clipper::String(lam, 10, 4) + " " + clipper::String(w, 12, 5) + "\n";
          }
        }
      }
      s += "EOD\n";
    }
    s += "\n";

    // plot: solid used segment + dashed extrapolation, one colour per range
    s += "plot \\\n";
    for (int ir = 0; ir < nr; ++ir) {
      const WavelengthRange& r = ranges[ir];
      std::string c = col[ir % 5];
      std::string idx = clipper::String(ir).trim();
      std::string legend;
      if (nr == 1) {
        legend = "w({/Symbol l})";
      } else {
        legend = "Range " + clipper::String(ir + 1).trim() + " ("
               + clipper::String(r.lam_min, 6, 2).trim() + "-"
               + clipper::String(r.lam_max, 6, 2).trim() + " A)";
      }
      s += "  $LNs" + idx + " using 1:2 with lines lw 2.5 lc rgb \"" + c
         + "\" title \"" + legend + "\"";
      bool hasextrap = (r.lam_min > gmin) || (r.lam_max < gmax);
      if (hasextrap) {
        s += ", \\\n";
        s += "  $LNd" + idx + " using 1:2 with lines lw 1.2 dashtype 2 lc rgb \""
           + c + "\" notitle";
      }
      if (ir < nr - 1) s += ", \\\n";
    }
    s += "\n";
    return s;
  }
  //--------------------------------------------------------------
  std::string WavelengthChebyshevScale::FormatSave() const
  {
    std::string s = "WavelengthChebyshevScale\n";
    s += "Nranges " + clipper::String(int(ranges.size())) + "\n";
    s += "LambdaRef " + clipper::String(lambda_ref) + "\n";
    for (size_t i=0; i<ranges.size(); ++i) {
      s += "Range " + clipper::String(int(i)) +
           " LamMin " + clipper::String(ranges[i].lam_min) +
           " LamMax " + clipper::String(ranges[i].lam_max) +
           " Degree " + clipper::String(ranges[i].degree) + "\n";
    }
    s += "Ncoeffs " + clipper::String(ncoeffs) + "\n";
    s += "Coefficients";
    for (int i=0; i<ncoeffs; ++i) {
      s += " " + clipper::String(coeffs[i]);
    }
    s += "\n";
    s += "End\n";
    return s;
  }
  //--------------------------------------------------------------
  void WavelengthChebyshevScale::Restore(Fileread& FR)
  {
    FR.ReadTag("Nranges"); int nr = FR.Int();
    FR.ReadTag("LambdaRef"); lambda_ref = FR.Double();
    ranges.resize(nr);
    ncoeffs = 0;
    for (int i=0; i<nr; ++i) {
      FR.ReadTag("Range"); FR.Int(); // index
      FR.ReadTag("LamMin"); ranges[i].lam_min = FR.Double();
      FR.ReadTag("LamMax"); ranges[i].lam_max = FR.Double();
      FR.ReadTag("Degree"); ranges[i].degree = FR.Int();
      ranges[i].offset = ncoeffs;
      ncoeffs += ranges[i].degree + 1;
    }
    FR.ReadTag("Ncoeffs"); int nc = FR.Int();
    if (nc != ncoeffs) {
      Message::message(Message_fatal("WavelengthChebyshevScale::Restore: ncoeffs mismatch"));
    }
    FR.ReadTag("Coefficients"); coeffs = FR.DoubleVec(ncoeffs);
    ref_range = rangeIndex(lambda_ref);
    f_ref = (ref_range >= 0) ? evalRange(lambda_ref, ref_range) : 1.0;
    if (!FR.CheckEnd()) {
      Message::message(Message_warn
        ("WavelengthChebyshevScale::Restore unexpected tag "+FR.Tag()));
    }
  }
  //--------------------------------------------------------------
  // =================== WavelengthGPRScale =========================
  //--------------------------------------------------------------
  double WavelengthGPRScale::kernelValue(const double& xa,
                                          const double& xb) const
  // Covariance between two points in the warped coordinate x = log(lambda),
  // using the stored hyperparameters
  {
    double sigf2 = sigf * sigf;
    double r = (xa - xb) / lengthscale;
    if (kernel == MATERN32) {
      double d = std::sqrt(3.0) * std::fabs(r);
      return sigf2 * (1.0 + d) * std::exp(-d);
    }
    // squared-exponential (default)
    return sigf2 * std::exp(-0.5 * r * r);
  }
  //--------------------------------------------------------------
  namespace {
    // one per-observation sample for the GP training set
    struct GPRSample {
      double lam;   // wavelength
      double rho;   // I_obs / <I>_mates
      double u;     // weight for rho: <I>_mates / sigma^2
    };
    bool GPRSampleByLambda(const GPRSample& a, const GPRSample& b)
    {return a.lam < b.lam;}
  }
  //--------------------------------------------------------------
  int WavelengthGPRScale::Fit(const std::vector<double>& lambdas,
                               const std::vector<double>& ratios,
                               const std::vector<double>& weights,
                               const GPRControl& ctrl,
                               const double& LambdaRef,
                               std::string& fitlog)
  // Fit the GP and build the lookup table. Returns number of training bins.
  {
    active = false;
    grid_g.clear();
    grid_sd.clear();
    train_lam.clear();
    train_h.clear();
    train_hsd.clear();
    ntrain = 0;
    gmean = 0.0;
    noisescale = 1.0;
    kernel = ctrl.kernel;
    lambda_ref = LambdaRef;

    const int Nsamp = int(lambdas.size());
    if (Nsamp < 50) {
      fitlog += " GPR wavelength fit: too few samples ("
        + clipper::String(Nsamp) + "); no correction applied\n";
      return 0;
    }

    // ----- wavelength range
    lam_min = ctrl.lam_min;
    lam_max = ctrl.lam_max;
    if (lam_min <= 0.0 || lam_max <= 0.0 || lam_max <= lam_min) {
      lam_min =  std::numeric_limits<double>::max();
      lam_max = -std::numeric_limits<double>::max();
      for (int n = 0; n < Nsamp; ++n) {
        lam_min = std::min(lam_min, lambdas[n]);
        lam_max = std::max(lam_max, lambdas[n]);
      }
    }
    double range = lam_max - lam_min;
    if (range <= 0.0 || lam_min <= 0.0) {
      fitlog += " GPR wavelength fit: degenerate wavelength range; skipped\n";
      return 0;
    }

    // ----- usable samples, sorted by wavelength (for equal-count binning)
    std::vector<GPRSample> smp;
    smp.reserve(Nsamp);
    for (int n = 0; n < Nsamp; ++n) {
      double lam = lambdas[n];
      if (lam < lam_min || lam > lam_max) continue;
      double u = weights[n];
      double rho = ratios[n];
      if (!(u > 0.0) || !std::isfinite(u) || !std::isfinite(rho)) continue;
      GPRSample s;
      s.lam = lam; s.rho = rho; s.u = u;
      smp.push_back(s);
    }
    const int N = int(smp.size());
    if (N < 50) {
      fitlog += " GPR wavelength fit: too few usable samples ("
        + clipper::String(N) + "); no correction applied\n";
      return 0;
    }
    std::sort(smp.begin(), smp.end(), GPRSampleByLambda);

    // ----- Binning: EQUAL-COUNT, with a cap on the bin width in log(lambda).
    // Uniform-width bins give wildly unequal precision across the spectrum (a
    // Laue dataset has orders of magnitude more observations near the peak of
    // the spectrum than in its tails), and the resulting noisy tail bins are
    // what makes an unconstrained GP oscillate.  Purely equal-count bins fix
    // that but collapse the sparse tail into one very wide bin whose centroid
    // sits far from the end of the range, leaving the GP to extrapolate (and
    // run away) over the last stretch.  Closing a bin on whichever limit is
    // reached first keeps training points spread over the whole range while
    // still giving each one usable precision.
    const int MINCOUNT = 20;
    int nbins = (ctrl.nbins > 0) ? ctrl.nbins : 50;
    if (nbins > N / MINCOUNT) nbins = N / MINCOUNT;
    if (nbins < 4) {
      fitlog += " GPR wavelength fit: too few samples for binning ("
        + clipper::String(N) + "); no correction applied\n";
      return 0;
    }
    const double xspan = std::log(smp[N-1].lam) - std::log(smp[0].lam);
    const double target = double(N) / nbins;            // nominal bin population
    const double wcap = 4.0 * xspan / nbins;            // max bin width, log(lambda)
    const int mincnt = std::max(MINCOUNT, int(0.25 * target));
    std::vector<int> bstart, bend;
    {
      int i0 = 0;
      while (i0 < N) {
        double x0 = std::log(smp[i0].lam);
        int i1 = i0 + 1;
        while (i1 < N) {
          int cnt = i1 - i0;
          if (cnt >= mincnt &&
              (cnt >= target || std::log(smp[i1-1].lam) - x0 >= wcap)) break;
          ++i1;
        }
        if (N - i1 < mincnt) i1 = N;   // absorb a short trailing bin
        bstart.push_back(i0);
        bend.push_back(i1);
        i0 = i1;
      }
    }

    // ----- one training point per bin, by the linear-space ratio estimator
    //   r  = sum(u rho) / sum(u)              (unbiased, tolerates rho <= 0)
    //   var(r) = sum(u^2 (rho-r)^2) / (sum u)^2     (sandwich/ratio variance)
    //   target h = log(r),  var(h) = var(r)/r^2
    // A mean of per-observation log-ratios (the previous estimator) is badly
    // biased and heavy-tailed for weak observations, and its "SEM^2" used the
    // raw bin count even though 1/sigma^2 weights make the effective sample
    // size several times smaller -- so the noise was underestimated by a large
    // factor and the marginal likelihood drove the length scale to its floor.
    std::vector<double> tl, ty, tn;   // wavelength, target, noise variance
    for (size_t b = 0; b < bstart.size(); ++b) {
      int i0 = bstart[b];
      int i1 = bend[b];
      if (i1 - i0 < MINCOUNT) continue;
      double su = 0.0, sur = 0.0, sul = 0.0;
      for (int i = i0; i < i1; ++i) {
        su  += smp[i].u;
        sur += smp[i].u * smp[i].rho;
        sul += smp[i].u * smp[i].lam;
      }
      if (!(su > 0.0)) continue;
      double r = sur / su;
      if (!(r > 0.0)) continue;      // cannot take a log of a non-positive mean
      double sv = 0.0;
      for (int i = i0; i < i1; ++i) {
        double e = smp[i].u * (smp[i].rho - r);
        sv += e * e;
      }
      double varr = sv / (su * su);
      double varh = varr / (r * r);
      if (!std::isfinite(varh)) continue;
      tl.push_back(sul / su);        // weighted-mean wavelength in the bin
      ty.push_back(std::log(r));
      tn.push_back(varh);
    }
    int B = int(tl.size());
    if (B < 4) {
      fitlog += " GPR wavelength fit: too few usable bins ("
        + clipper::String(B) + "); no correction applied\n";
      return 0;
    }

    // ----- noise floor, relative to the median bin variance
    std::vector<double> tnsort = tn;
    std::sort(tnsort.begin(), tnsort.end());
    double nmedian = tnsort[tnsort.size()/2];
    double nfloor = std::max(1.0e-6, 0.01 * nmedian);
    for (int i = 0; i < B; ++i) tn[i] = std::max(tn[i], nfloor);

    // ----- WARPED abscissa x = log(lambda).  A white-beam spectrum varies
    // rapidly over a narrow interval at short wavelength and slowly over a wide
    // interval at long wavelength; in log(lambda) a single stationary length
    // scale describes both, so the fit no longer has to choose between
    // following the short-wavelength rise and smoothing the long-wavelength
    // tail (the compromise that produced the spurious oscillations).
    std::vector<double> tx(B);
    for (int i = 0; i < B; ++i) tx[i] = std::log(tl[i]);
    double xmin = tx[0], xmax = tx[0];
    for (int i = 1; i < B; ++i) {
      xmin = std::min(xmin, tx[i]);
      xmax = std::max(xmax, tx[i]);
    }
    double xrange = xmax - xmin;
    if (xrange <= 0.0) {
      fitlog += " GPR wavelength fit: degenerate training range; skipped\n";
      return 0;
    }
    std::vector<double> dxs;
    for (int i = 1; i < B; ++i) dxs.push_back(std::fabs(tx[i] - tx[i-1]));
    std::sort(dxs.begin(), dxs.end());
    double dxmed = dxs[dxs.size()/2];

    // ----- constant GP mean.  Away from the data the posterior reverts to this
    // level rather than to zero, so sparsely sampled wavelengths get the mean
    // response instead of an artificial dip towards ws = exp(-g_ref).
    double swm = 0.0, swmy = 0.0;
    for (int i = 0; i < B; ++i) {
      double wi = 1.0 / tn[i];
      swm += wi; swmy += wi * ty[i];
    }
    gmean = swmy / swm;
    Eigen::VectorXd Y(B);
    for (int i = 0; i < B; ++i) Y(i) = ty[i] - gmean;

    // ----- hyperparameter candidates.  sigma_f and a noise-inflation factor
    // are searched jointly with the length scale: fixing sigma_f at the raw
    // spread of the bin targets (which includes the noise) biased the fit
    // towards explaining noise as signal.
    double yvar = 0.0, nmean = 0.0;
    for (int i = 0; i < B; ++i) {
      yvar += Y(i) * Y(i);
      nmean += tn[i];
    }
    yvar /= B; nmean /= B;
    double sig0 = std::sqrt(std::max(yvar - nmean, 1.0e-4));

    std::vector<double> ellcand;
    double elllo = std::max(2.0 * dxmed, 0.02 * xrange);
    double ellhi = 1.5 * xrange;
    bool ellfixed = (ctrl.lengthscale > 0.0);
    if (ellfixed) {
      // user length scale is given in Angstrom; convert to log(lambda) units
      // at the reference wavelength (dx = dlambda/lambda)
      double lscale = (lambda_ref > 0.0) ? lambda_ref : std::sqrt(lam_min*lam_max);
      ellcand.push_back(ctrl.lengthscale / lscale);
    } else {
      const int NELL = 16;
      for (int k = 0; k < NELL; ++k) {
        double f = double(k) / (NELL - 1);
        ellcand.push_back(elllo * std::pow(ellhi / elllo, f));  // log-spaced
      }
    }
    std::vector<double> sigcand;
    const int NSIGF = 9;
    for (int k = 0; k < NSIGF; ++k) {
      double f = -0.7 + 1.4 * double(k) / (NSIGF - 1);
      sigcand.push_back(sig0 * std::pow(10.0, f));
    }
    std::vector<double> alfcand;
    const int NALF = 8;
    for (int k = 0; k < NALF; ++k) {
      double f = -0.3 + 1.5 * double(k) / (NALF - 1);
      alfcand.push_back(std::pow(10.0, f));
    }

    const double jitter = 1.0e-10;
    double bestlml = -std::numeric_limits<double>::max();
    double bestell = ellcand[0], bestsig = sigcand[0], bestalf = alfcand[0];
    Eigen::VectorXd bestalpha;
    Eigen::MatrixXd K0(B, B);
    for (size_t e = 0; e < ellcand.size(); ++e) {
      lengthscale = ellcand[e];
      sigf = 1.0;                     // unit-amplitude kernel, scaled below
      for (int i = 0; i < B; ++i)
        for (int j = 0; j < B; ++j) K0(i, j) = kernelValue(tx[i], tx[j]);
      for (size_t sgi = 0; sgi < sigcand.size(); ++sgi) {
        double sf2 = sigcand[sgi] * sigcand[sgi];
        for (size_t ai = 0; ai < alfcand.size(); ++ai) {
          Eigen::MatrixXd A = sf2 * K0;
          for (int i = 0; i < B; ++i) A(i, i) += alfcand[ai] * tn[i] + jitter;
          Eigen::LLT<Eigen::MatrixXd> llt(A);
          if (llt.info() != Eigen::Success) continue;
          Eigen::VectorXd alpha = llt.solve(Y);
          double logdet = 0.0;
          Eigen::MatrixXd L = llt.matrixL();
          for (int i = 0; i < B; ++i) logdet += std::log(L(i, i));
          double lml = -0.5 * Y.dot(alpha) - logdet
                       - 0.5 * B * std::log(2.0 * M_PI);
          if (lml > bestlml) {
            bestlml = lml;
            bestell = ellcand[e];
            bestsig = sigcand[sgi];
            bestalf = alfcand[ai];
            bestalpha = alpha;
          }
        }
      }
    }
    if (bestalpha.size() == 0) {
      fitlog += " GPR wavelength fit: Cholesky factorisation failed; skipped\n";
      return 0;
    }
    lengthscale = bestell;
    sigf = bestsig;
    noisescale = bestalf;

    // mean training noise (for reporting)
    double nsum = 0.0;
    for (int i = 0; i < B; ++i) nsum += std::sqrt(noisescale * tn[i]);
    signoise = nsum / B;

    // ----- refactorise at the chosen hyperparameters (for the posterior variance)
    Eigen::LLT<Eigen::MatrixXd> bestllt;
    {
      Eigen::MatrixXd A(B, B);
      for (int i = 0; i < B; ++i) {
        for (int j = 0; j < B; ++j) A(i, j) = kernelValue(tx[i], tx[j]);
        A(i, i) += noisescale * tn[i] + jitter;
      }
      bestllt.compute(A);
    }

    // ----- evaluate posterior mean and SD on a dense lookup grid, uniform in
    // lambda (so Scale() stays a cheap linear interpolation)
    //   mean g(lam)  = gmean + k_*^T alpha
    //   var  g(lam)  = k(x,x) - v^T v,  v = L^{-1} k_*
    // Beyond the outermost training points the mean is HELD CONSTANT at the
    // edge value.  There the posterior is driven by the local gradient of the
    // last bins, which is exactly the runaway extrapolation that makes a
    // polynomial fit diverge at the sparse ends of the spectrum.  The SD is
    // still evaluated at the true wavelength, so the widening band continues to
    // show that the extremes are unsupported by data.
    const double xlim_lo = tx[0], xlim_hi = tx[B-1];
    ngrid = std::max(200, 4 * nbins);
    grid_step = range / (ngrid - 1);
    grid_g.assign(ngrid, 0.0);
    grid_sd.assign(ngrid, 0.0);
    Eigen::VectorXd kstar(B);
    for (int m = 0; m < ngrid; ++m) {
      double x = std::log(lam_min + m * grid_step);
      double xg = std::min(std::max(x, xlim_lo), xlim_hi);
      double g = 0.0;
      for (int i = 0; i < B; ++i) g += kernelValue(xg, tx[i]) * bestalpha(i);
      grid_g[m] = gmean + g;
      for (int i = 0; i < B; ++i) kstar(i) = kernelValue(x, tx[i]);
      double kss = kernelValue(x, x);
      Eigen::VectorXd v = bestllt.matrixL().solve(kstar);
      double var = kss - v.dot(v);
      if (var < 0.0) var = 0.0;   // guard against round-off
      grid_sd[m] = std::sqrt(var);
    }

    // ----- reference value g(lambda_ref)
    if (lambda_ref >= lam_min && lambda_ref <= lam_max) {
      double xref = std::min(std::max(std::log(lambda_ref), xlim_lo), xlim_hi);
      double g = 0.0;
      for (int i = 0; i < B; ++i) g += kernelValue(xref, tx[i]) * bestalpha(i);
      g_ref = gmean + g;
    } else {
      g_ref = gmean;
    }

    train_lam = tl;
    train_h   = ty;
    train_hsd.resize(B);
    for (int i = 0; i < B; ++i) train_hsd[i] = std::sqrt(noisescale * tn[i]);

    ntrain = B;
    active = true;
    fitlog += " GPR wavelength fit: " + clipper::String(B)
      + " bins from " + clipper::String(N)
      + " observations (nominally " + clipper::String(int(target)) + " per bin)\n"
      + "   length scale " + clipper::String(lengthscale, 6, 4)
      + " in log(lambda) (" + clipper::String(lengthscale * lambda_ref, 6, 4)
      + " A at the reference wavelength), sigma_f "
      + clipper::String(sigf, 6, 4)
      + ", noise inflation " + clipper::String(noisescale, 6, 2) + "\n"
      + "   log marginal likelihood " + clipper::String(bestlml, 8, 3) + "\n";
    if (!ellfixed) {
      if (lengthscale < 1.01 * elllo) {
        fitlog += "   WARNING: length scale is at the lower search limit; the fit"
                  " may be following noise -- consider fewer NORMGPRBINS\n";
      } else if (lengthscale > 0.99 * ellhi) {
        fitlog += "   WARNING: length scale is at the upper search limit; the fit"
                  " is nearly featureless\n";
      }
    }
    return B;
  }
  //--------------------------------------------------------------
  double WavelengthGPRScale::Scale(const double& lambda) const
  // ws(lambda) = exp(g(lambda) - g(lambda_ref)); 1.0 if inactive or out of range
  {
    if (!active || ngrid < 2) return 1.0;
    if (lambda < lam_min || lambda > lam_max) return 1.0;
    double x = (lambda - lam_min) / grid_step;
    int i = int(x);
    if (i < 0) i = 0;
    if (i >= ngrid - 1) i = ngrid - 2;
    double f = x - i;
    double g = grid_g[i] * (1.0 - f) + grid_g[i + 1] * f;
    return std::exp(g - g_ref);
  }
  //--------------------------------------------------------------
  double WavelengthGPRScale::Uncertainty(const double& lambda) const
  // GP posterior SD in log space ~ relative (fractional) uncertainty of ws.
  {
    if (!active || ngrid < 2 || int(grid_sd.size()) != ngrid) return 0.0;
    if (lambda < lam_min || lambda > lam_max) return 0.0;
    double x = (lambda - lam_min) / grid_step;
    int i = int(x);
    if (i < 0) i = 0;
    if (i >= ngrid - 1) i = ngrid - 2;
    double f = x - i;
    return grid_sd[i] * (1.0 - f) + grid_sd[i + 1] * f;
  }
  //--------------------------------------------------------------
  std::string WavelengthGPRScale::PrintNormalization(const int& npoints) const
  {
    if (!IsActive()) return "";
    std::string s = "\n Wavelength normalization (Gaussian process)\n";
    s += FormatOutput::logTabPrintf(1,
           " Reference wavelength: %6.4f A\n", lambda_ref);
    s += FormatOutput::logTabPrintf(1,
           " Range: %6.4f - %6.4f A   kernel: %s\n",
           lam_min, lam_max, (kernel == MATERN32) ? "Matern-3/2" : "squared-exp");
    s += FormatOutput::logTabPrintf(1,
           " Length scale: %7.4f in log(lambda) (%7.4f A at lambda_ref)\n",
           lengthscale, lengthscale * lambda_ref);
    s += FormatOutput::logTabPrintf(1,
           " sigma_f: %7.4f   noise inflation: %6.2f   training bins: %d\n",
           sigf, noisescale, ntrain);
    s += FormatOutput::logTabPrintf(1, "   %8s  %10s\n", "lambda", "w(lambda)");
    double step = (lam_max - lam_min) / (npoints - 1);
    for (int ip = 0; ip < npoints; ++ip) {
      double lam = lam_min + ip * step;
      double w = Scale(lam);
      std::string marker = (std::fabs(lam - lambda_ref) < 0.5 * step) ? " <- ref" : "";
      s += FormatOutput::logTabPrintf(1,
             "   %8.4f  %10.5f%s\n", lam, w, marker.c_str());
    }
    s += AsciiPlot();
    return s;
  }
  //--------------------------------------------------------------
  std::string WavelengthGPRScale::AsciiPlot(const int& width,
                                            const int& height) const
  // ASCII line plot of w(lambda) over the fitted range.  The '*' trace shows
  // the normalization curve; the reference wavelength is marked with a ':' column.
  {
    if (!active || lam_max <= lam_min || width < 2 || height < 2) return "";

    // Sample w(lambda) across the range
    std::vector<double> w(width);
    double wmin = 1.0e30, wmax = -1.0e30;
    for (int col = 0; col < width; ++col) {
      double lam = lam_min + (lam_max - lam_min) * col / double(width - 1);
      w[col] = Scale(lam);
      wmin = std::min(wmin, w[col]);
      wmax = std::max(wmax, w[col]);
    }
    if (wmax - wmin < 1.0e-6) {  // flat curve: pad so it sits mid-plot
      double mid = 0.5 * (wmin + wmax);
      wmin = mid - 0.5;
      wmax = mid + 0.5;
    }

    // Column of the reference wavelength, -1 if outside the range
    int refcol = -1;
    if (lambda_ref >= lam_min && lambda_ref <= lam_max) {
      refcol = int((lambda_ref - lam_min) / (lam_max - lam_min)
                   * (width - 1) + 0.5);
    }

    // Build grid (height rows top=wmax, bottom=wmin)
    std::vector<std::string> grid(height, std::string(width, ' '));
    if (refcol >= 0) {
      for (int row = 0; row < height; ++row) grid[row][refcol] = ':';
    }
    for (int col = 0; col < width; ++col) {
      int row = int((wmax - w[col]) / (wmax - wmin) * (height - 1) + 0.5);
      if (row < 0) row = 0;
      if (row >= height) row = height - 1;
      grid[row][col] = '*';
    }

    std::string s = FormatOutput::logTabPrintf(1, "   w(lambda):\n");
    for (int row = 0; row < height; ++row) {
      double wlabel = wmax - (wmax - wmin) * row / double(height - 1);
      s += FormatOutput::logTabPrintf(1, "   %8.4f |%s\n",
                                      wlabel, grid[row].c_str());
    }
    // x-axis
    s += FormatOutput::logTabPrintf(1, "            +%s\n",
                                    std::string(width, '-').c_str());
    std::string lo = clipper::String(lam_min, 6, 3).trim();
    std::string hi = clipper::String(lam_max, 6, 3).trim();
    std::string axis(width, ' ');
    for (size_t i = 0; i < lo.size() && i < axis.size(); ++i) axis[i] = lo[i];
    int hipos = width - int(hi.size());
    if (hipos < 0) hipos = 0;
    for (size_t i = 0; hipos + i < axis.size() && i < hi.size(); ++i)
      axis[hipos + i] = hi[i];
    s += FormatOutput::logTabPrintf(1, "            %s  lambda (A)\n",
                                    axis.c_str());
    if (refcol >= 0) {
      s += FormatOutput::logTabPrintf(1,
             "            (':' marks reference wavelength %6.4f A)\n",
             lambda_ref);
    }
    return s;
  }
  //--------------------------------------------------------------
  std::string WavelengthGPRScale::asXML() const
  {
    if (!IsActive()) return "";
    std::string s = "<WavelengthNormalisationGPR>\n";
    s += "  " + StringUtil::MakeXMLtag("ReferenceWavelength", lambda_ref, 8, 4);
    s += "\n";
    s += "  " + StringUtil::MakeXMLtag("LambdaMin", lam_min, 8, 4);
    s += StringUtil::MakeXMLtag("LambdaMax", lam_max, 8, 4) + "\n";
    s += "  " + StringUtil::MakeXMLtag("Kernel",
           std::string((kernel == MATERN32) ? "Matern-3/2" : "squared-exp"));
    s += "\n";
    // length scale is in log(lambda) units (the kernel abscissa)
    s += "  " + StringUtil::MakeXMLtag("LengthScale", lengthscale, 8, 4);
    s += StringUtil::MakeXMLtag("SigmaF", sigf, 8, 4);
    s += StringUtil::MakeXMLtag("NoiseInflation", noisescale, 8, 3);
    s += StringUtil::MakeXMLtag("TrainingBins", ntrain, 4) + "\n";
    s += "  <Normalisation>\n";
    const int npoints = 21;
    double step = (lam_max - lam_min) / (npoints - 1);
    for (int ip = 0; ip < npoints; ++ip) {
      double lam = lam_min + ip * step;
      s += "    <point>";
      s += StringUtil::MakeXMLtag("lambda", lam, 8, 4);
      s += StringUtil::MakeXMLtag("w", Scale(lam), 10, 5);
      s += StringUtil::MakeXMLtag("uncertainty", Uncertainty(lam), 10, 5);
      s += " </point>\n";
    }
    s += "  </Normalisation>\n";
    s += "</WavelengthNormalisationGPR>\n";
    return s;
  }
  //--------------------------------------------------------------
  std::string WavelengthGPRScale::GnuplotScript(const std::string& title,
                                                const std::string& version) const
  // Self-contained gnuplot script (LAMBDANORM): w(lambda) with a 1-sigma band.
  {
    if (!IsActive()) return "";
    // sanitise title for use inside a double-quoted gnuplot string
    std::string t = title;
    for (size_t i = 0; i < t.size(); ++i) {
      if (t[i] == '"') t[i] = '\'';
      if (t[i] == '\n' || t[i] == '\r') t[i] = ' ';
    }
    t = clipper::String(t).trim();
    std::string lref = clipper::String(lambda_ref, 8, 4).trim();

    std::string s;
    s += "# LAMBDANORM - Laue wavelength normalization curve\n";
    s += "# Generated by " + version + "\n";
    if (!t.empty()) s += "# Title: " + t + "\n";
    s += "#\n";
    s += "# View this plot with gnuplot:\n";
    s += "#     gnuplot -p LAMBDANORM\n";
    s += "# (the -p flag keeps the plot window open after drawing)\n";
    s += "#\n";
    s += "# Method: Gaussian-process wavelength normalization, fitted in log space.\n";
    s += "# w(lambda) = exp(g(lambda) - g(lambda_ref)); reference lambda = "
         + lref + " A.\n";
    s += "# Data columns:\n";
    s += "#   1 lambda (A)\n";
    s += "#   2 w(lambda)         normalization scale\n";
    s += "#   3 rel_uncertainty   GP posterior SD in log space (fractional error on w)\n";
    s += "#   4 w_lo = w*(1-rel)  lower 1-sigma\n";
    s += "#   5 w_hi = w*(1+rel)  upper 1-sigma\n";
    s += "#\n";
    s += "# $LAMBDABINS holds the binned observations the GP was fitted to\n";
    s += "# (lambda, w_bin, sigma(w_bin)); the curve should follow them without\n";
    s += "# chasing individual points.\n";
    s += "\n";
    s += "set title \"Laue wavelength normalization (GP)";
    if (!t.empty()) s += "\\n" + t;
    s += "  -  " + version + "\"\n";
    s += "set xlabel \"Wavelength {/Symbol l} / Angstrom\"\n";
    s += "set ylabel \"Normalization w({/Symbol l})\"\n";
    s += "set grid\n";
    s += "set key top right\n";
    s += "set xrange [" + clipper::String(lam_min, 8, 4).trim() + ":"
         + clipper::String(lam_max, 8, 4).trim() + "]\n";
    s += "set style fill transparent solid 0.30 noborder\n";
    s += "set arrow from " + lref + ", graph 0 to " + lref
         + ", graph 1 nohead dashtype 2 lc rgb \"#888888\"\n";
    s += "set label \"ref\" at " + lref + ", graph 0.96 center tc rgb \"#888888\"\n";
    s += "\n";
    s += "$LAMBDANORM << EOD\n";
    const int npoints = 100;
    double step = (lam_max - lam_min) / (npoints - 1);
    for (int ip = 0; ip < npoints; ++ip) {
      double lam = lam_min + ip * step;
      double w = Scale(lam);
      double u = Uncertainty(lam);
      s += clipper::String(lam, 10, 4) + " " + clipper::String(w, 11, 5) + " "
         + clipper::String(u, 11, 5) + " "
         + clipper::String(w * (1.0 - u), 11, 5) + " "
         + clipper::String(w * (1.0 + u), 11, 5) + "\n";
    }
    s += "EOD\n";
    s += "\n";
    s += "$LAMBDABINS << EOD\n";
    for (size_t ib = 0; ib < train_lam.size(); ++ib) {
      double wb = std::exp(train_h[ib] - g_ref);
      s += clipper::String(train_lam[ib], 10, 4) + " "
         + clipper::String(wb, 11, 5) + " "
         + clipper::String(wb * train_hsd[ib], 11, 5) + "\n";
    }
    s += "EOD\n";
    s += "\n";
    s += "plot \\\n";
    s += "  $LAMBDANORM using 1:4:5 with filledcurves lc rgb \"#9ec3e8\" "
         "title \"1{/Symbol s} band\", \\\n";
    s += "  $LAMBDABINS using 1:2:3 with yerrorbars pt 7 ps 0.5 "
         "lc rgb \"#999999\" title \"binned data\", \\\n";
    s += "  $LAMBDANORM using 1:2 with lines lw 2 lc rgb \"#1f4e96\" "
         "title \"w({/Symbol l})\"\n";
    return s;
  }
  //--------------------------------------------------------------
  std::string WavelengthGPRScale::FormatSave() const
  {
    std::string s = "WavelengthGPRScale\n";
    s += "Active " + clipper::String(active ? 1 : 0) + "\n";
    s += "LamMin " + clipper::String(lam_min) +
         " LamMax " + clipper::String(lam_max) +
         " LambdaRef " + clipper::String(lambda_ref) +
         " Gref " + clipper::String(g_ref) + "\n";
    s += "Kernel " + clipper::String(int(kernel)) +
         " Lengthscale " + clipper::String(lengthscale) +
         " Sigf " + clipper::String(sigf) +
         " Gmean " + clipper::String(gmean) +
         " Noisescale " + clipper::String(noisescale) +
         " Ntrain " + clipper::String(ntrain) + "\n";
    s += "Ngrid " + clipper::String(ngrid) +
         " GridStep " + clipper::String(grid_step) + "\n";
    s += "GridG\n" + StringUtil::FormatSaveVector(grid_g);
    s += "GridSD\n" + StringUtil::FormatSaveVector(grid_sd);
    s += "End\n";
    return s;
  }
  //--------------------------------------------------------------
  void WavelengthGPRScale::Restore(Fileread& FR)
  {
    FR.ReadTag("Active"); active = (FR.Int() != 0);
    FR.ReadTag("LamMin"); lam_min = FR.Double();
    FR.ReadTag("LamMax"); lam_max = FR.Double();
    FR.ReadTag("LambdaRef"); lambda_ref = FR.Double();
    FR.ReadTag("Gref"); g_ref = FR.Double();
    FR.ReadTag("Kernel"); kernel = KernelType(FR.Int());
    FR.ReadTag("Lengthscale"); lengthscale = FR.Double();
    FR.ReadTag("Sigf"); sigf = FR.Double();
    FR.ReadTag("Gmean"); gmean = FR.Double();
    FR.ReadTag("Noisescale"); noisescale = FR.Double();
    FR.ReadTag("Ntrain"); ntrain = FR.Int();
    FR.ReadTag("Ngrid"); ngrid = FR.Int();
    FR.ReadTag("GridStep"); grid_step = FR.Double();
    FR.ReadTag("GridG"); grid_g = FR.DoubleVec(ngrid);
    FR.ReadTag("GridSD"); grid_sd = FR.DoubleVec(ngrid);
    if (!FR.CheckEnd()) {
      Message::message(Message_warn
        ("WavelengthGPRScale::Restore unexpected tag "+FR.Tag()));
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
