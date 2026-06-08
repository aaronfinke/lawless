// resolutionlimit.cpp
//
// Class to determine resolution "limit"
//

#include "resolutionlimit.hh"
#include "radialfunction.hh"
#include "fitresolutiondata.hh"
#include "radialfunction.hh"
#include "simpleminimise.hh"
#include "median.hh"
#include "string_util.hh"
#include "score_datatypes.hh"

#include <assert.h>
#define ASSERT assert

namespace scala {
  // ------------------------------------------------------------
  ResolutionLimit::ResolutionLimit(const std::vector<double> score,
                                   const ResoRange& ResRange,
                                   const double& Limit,
                                   const FitType& Fittype)
    : limit(0.0), highres(0.0), status(-3), sufficientdata(false),
      nrej(0), fittype(NONE)
  //! construct from score list, resolution range and minimum score
  {
    init(score, ResRange, Limit, Fittype);
  }
  // ------------------------------------------------------------
  void ResolutionLimit::init(const std::vector<double> score,
                             const ResoRange& ResRange,
                             const double& Limit,
                             const FitType& Fittype)
  //! initialise from score list, resolution range and minimum score
  {
    ASSERT (int(score.size()) == ResRange.Nbins());
    //! A high resolution limit is determined from a score list for
    //! each resolution bin, by finding the point at which the score
    //! drops below the limit. If the scores in all resolution
    //! bins are above the limit, then the maximum resolution is used.

    fittype = Fittype;
    fitted = false;
    limit = Limit;
    status = -2;
    nrej = 0;
    sufficientdata = true;

    if (testinsufficientdata(score)) {
      // not enough data
      highres = 0.0;
      sufficientdata = false;
    }

    highres = ResRange.ResHigh();

    int i1 = -1;
    int nbins = ResRange.Nbins();
    if (nbins <= 1) {
      highres = ResRange.ResHigh();
      status = +1;
    } else {
      int nsafe = 0;
      int n = 0;
      int ibinmax = -1;  // last occupied bin
      for (int i=nbins-1;i>0;--i) { // loop backwards but ignore 1st bin
        if (score[i] != 0.0) {  // ignore scores of zero
          if (ibinmax < 0) {ibinmax = i;}
          n++;
          if (score[i] >= limit) {
            if (i1 < 0) {
              i1 = i;  // last one above limit
            }
            nsafe++;
          }
        }
      }
      lastfilledbin = ibinmax;

      if ((nsafe == n) && (ibinmax == nbins-1)) { // all above limit and all present
        highres = ResRange.ResHigh();
        status = +1;
      }

      if (i1 < 0) {
        // All bins below limit, clear highres & set status
        highres = 0.0;
        status = -1;
      }
    }

    if (!sufficientdata) {
      highres = 0.0;
      return;
    } // no fit unless there is sufficient data

    double reshigh = highres;

    if (fittype != NONE) { // fitting function to scores
      // Do the fit even if all are above or all below limit
      //  (flagged in status)
      reshigh = fit(score, ResRange);
      //      std::cout << "status "<<status<<" highres " <<highres<<" reshigh "<<reshigh<<std::endl;
      //  Status usually = -2
      //   = +1 all above threshold, and fitted (highres = max res)
      //   = -1 all below threshold, or fit failed (reshigh <= 0.0)
      if (status == -2) {
        if (reshigh > 0.0) {
          highres = 1.0/sqrt(reshigh);
          //std::cout << " highres " <<highres<<std::endl;
          status = 0;
        } else {
          if (limit > 0.001) {
            // negative highres, fall back to LINEAR
            fittype = LINEAR;
            highres = 0.0;
            status = -1;
          }
        }
        fitted = true;
      } else if (status == -1) {
        highres = 0.0;
        fitted = true;
      } else if (status >= 0) {
        fitted = true;
      }
    }
    if (status != 0 && status != -1) { // simply interpolate unless fit was successful
      if (nbins > 1) {
        if (fitted && (fittype == LINEAR)) {
          // we have a fall-back linear fit,
          //  find the point at which this falls below limit
        } else {
          // linear interpolate on 1/d^2 between bins i1 and i1+1
          if (status == +1 || i1+1 >= int(score.size())) {
            highres = ResRange.ResHigh();
          } else {
            double den = score[i1]-score[i1+1];
            if (score[i1+1] == 0.0) {
              den = score[i1];
            }
            double f = 1.0;
            if (den > 0.0) { // trap for divide by 0
              f = (score[i1]-limit)/den;
            }
            highres = ResRange.middle(i1) +
              f * (ResRange.middle(i1+1) - ResRange.middle(i1));
            if (highres > 0.0) highres = 1.0/sqrt(highres);
            //std::cout << "RL interpolate " <<i1<<" "<<highres <<" "<<reshigh<<" Status " << status<<std::endl;
          }
        }
        if (status == -2) {
          status = 0;
        }
      }
    }
    if (status != -1) {
      highres = std::max(highres, ResRange.ResHigh());
    }

    //std::cout << "RL resolution " <<highres <<" "<<reshigh<<" Status " << status
    //        <<std::endl<<std::endl;
  }
  // ------------------------------------------------------------
  ResolutionLimit::ResolutionLimit(const std::vector<MeanSD> mnsd,
                                   const ResoRange& ResRange,
                                   const double& Limit,
                                   const FitType& Fittype)
  //! construct from MeanSD list (score = Mean), resolution range and minimum score
  {
    init (mnsd, ResRange, Limit, Fittype);
  }
  // ------------------------------------------------------------
  void ResolutionLimit::init(const std::vector<MeanSD> mnsd,
                             const ResoRange& ResRange,
                             const double& Limit,
                             const FitType& Fittype)
  //! initialise from MeanSD list (score = Mean), resolution range and minimum score
  {
    std::vector<double> score(mnsd.size());
    for (size_t i=0;i<mnsd.size();++i) {
      score[i] = mnsd[i].Mean();
    }
    init (score, ResRange, Limit, Fittype);
    //std::cout << "ResolutionLimit " << highres <<std::endl; //^^
  }
  // ------------------------------------------------------------
  double ResolutionLimit::fit(const std::vector<double>& score,
                            const ResoRange& ResRange)
  // return reshigh (1/d^2)
  //   if no valid fit,
  //     reshigh (1/d^2) <= 0
  //     also status = -1, fittype = LINEAR
  // otherwise status is unchanged
  {
    // Straight line fit option
    LinearFit linefit;

    // Don't ignore 1st point if small number of ranges,
    // or first range is not low resolution
    const int MINSIZE = 5;
    const double RESOLOW = 10.0;
    bool ignorefirstpoint = true;
    if ((score.size() <= MINSIZE) ||
        (ResRange.middleA(0) < RESOLOW)) {
      ignorefirstpoint = false;
    }

    std::vector<ResolutionData> data(score.size());
    int ndata = 0; // number of valid data
    int nnegs = 0; // number of negative scores
    for (size_t i=0; i<score.size(); i++) {
      data[i].s = ResRange.middle(i);
      data[i].v = score[i];
      data[i].w = 1.0;      // weight, <=0 to ignore
      if (i == 0 && ignorefirstpoint) {
        // ignore first point unless small number of ranges
        data[i].w = 0.0;
      }
      if (score[i] == 0.0) {
        data[i].w = 0.0;      // weight, <=0 to ignore
      } else {
        ndata++;
        linefit.add(data[i].s, score[i], data[i].w);
      }
      if (score[i] < 0.0) {
        nnegs++;
      }
      //^
      //      printf("%8.3f %8.3f %8.3f\n", data[i].s, data[i].v, data[i].w);
    }

    double smax = ResRange.SResHigh();  // maximum s = 1/d^2

    // check straight-line fit for values which cannot fit the TANH function
    slope = linefit.result().first;
    intercept = linefit.result().second;
    if (slope > 0.0) {
      // Use Line fit
      fittype = LINEAR;
      //      std::cout << "ResolutionLimit::fit type set to LINEAR\n";
      return 0.0;
    }

    // Should we use 2 or 3 parameters? Only use 3 if there are a
    //  "significant" number of negatives
    int npar = 2;
    const double MAXNEGS = 0.3;
    if (double(nnegs)/double(ndata) > MAXNEGS) {
      npar = 3;
    }

    std::vector<double> params(npar);
    params[0] = 0.5*smax;  // d0
    params[1] = smax/double(ndata);      //r
    //    std::cout << "Start params " << params[0] <<" "<<params[1]<<std::endl; //^-
    if (npar == 3) {
      params[2] = 1.0;       // dcc, offset to allow for CC < 0
    }

    radialfunction.init(params);  // initial parameters
    // set up data and radial function
    FitResolutionData fitresolutiondata(data, radialfunction);
    // Use LNCOSH target
    fitresolutiondata.setQuadratic(false);
    ////fitresolutiondata.setQuadratic(true);  // quadratic

    int ncycles = 10;
    double tolerance = 0.01;
    double damp = 0.2;
    DampedGaussNewton dampedgaussnewton(fitresolutiondata,
                                        ncycles, tolerance, damp);

    double reshigh = -1.0;  // current best estimate of 1/d^2
    if (dampedgaussnewton.valid()) {
      reshigh = radialfunction.inverse(limit);
    }

    // Check for outliers
    // multiple of "error" (inter-quartile range of differences)
    int nrej = 0;
    if (ndata > npar+2) {
      double reject = 4.0;
      nrej = rejectoutliers(data, reject);
      if ((ndata-nrej) <= npar+2) {
        // give up
        return reshigh;  // -1
      }
    }
    if (nrej > 0) {
      //^^
      //      std::cout << "ResolutionLimit::fit before reject\n"<<
      //        radialfunction.format() << "\n";
      //      reshigh = radialfunction.inverse(limit);
      //      if (reshigh > 0.0) {
      //        reshigh = 1.0/sqrt(reshigh);
      //      }
      //      std::cout <<"Rejected "<<nrej<<" High res before " << reshigh<<"\n";
      //^-
      dampedgaussnewton.run(fitresolutiondata,
                        ncycles, tolerance, damp);
      if (dampedgaussnewton.valid()) {
        reshigh = radialfunction.inverse(limit);
      }
    }

    //    std::cout << "ResolutionLimit::fit ResHigh "<<reshigh
    //        <<"\n"<< radialfunction.format() << "\n\n"; //^
    if (limit > 0.0) {  // valid limit given
      if (reshigh < 0.0) {
        // negative highres, fall back to LINEAR
        fittype = LINEAR;
        if (slope >= 0.0) {
          // no answer if slope not negative
          reshigh = 0.0;
          status = -1;
        } else {
          reshigh = linefit.xvalueaty(limit);  // may still be negative
          if (reshigh <= 0.0) {
            reshigh = 0.0;
            status = -1;
          }
        }
      }
    }
    return reshigh;
  }
  // ------------------------------------------------------------
  int ResolutionLimit::rejectoutliers(std::vector<ResolutionData>& data,
                                      const double& reject) const
  {
    std::vector<double> deltaabs;  // vector of obs-calc differences for median
    std::vector<double> delta(data.size(), 0.0);  // vector of obs-calc differences

    nrej = 0;

    int j = 0;
    for (size_t k=0; k<data.size(); k++) { // loop data
      if (data[k].w > 0.0) {
        double s = data[k].s;
        double v = data[k].v;
        double vcalc = radialfunction.value(s);
        deltaabs.push_back(std::abs(v - vcalc));
        delta[k] = (v - vcalc);
        //      std::cout << "j, v, vc, del " <<j<<" "
        //                <<v<<" "<<vcalc<<" "<<v-vcalc<<"\n";
        j++;
      }
    }

    Median<double> mdn(deltaabs);
    double median = mdn.median();
    double iqr = mdn.interquartilerange();

    //    std::cout << "Median difference: " << median<<", IQR "<<iqr<<std::endl;

    for (size_t k=0; k<data.size(); k++) { // loop data
      double w = data[k].w;
      if (w > 0) {
        double rdel = (std::abs(delta[k]) - median) / iqr;
        if (std::abs(rdel) > reject) {
          // reject this one
          data[k].w = 0.0;
          //      std::cout << "Reject item "<<k<<", rdel = "<<rdel<<"\n";
          nrej++;
        }
      }
    }
    return nrej;
  }
  // ------------------------------------------------------------
  std::string ResolutionLimit::format(const bool& anomalous) const
  // anomalous == true for assessment of anomalous signal (changes wording)
  {
    std::string s;
    if (status <= -3) {
      s += "Resolution limit not initialised\n";
      return s;
    }
    if (status <= -2) {
      s += "No resolution limit determined\n";
      return s;
    }
    s += "Threshold (see ANALYSIS keyword): "+
        StringUtil::ftos(limit, 6, 2)+"\n";

    if (sufficientdata) {  // sufficent data to determine limit
      if (fittype == NONE) { // no function fit
        s += "Resolution limit determined from the point at which the score drops below threshold\n";
      } else if (fittype == LINEAR) {
        s += "Resolution limit determined from a linear fit\n";
      } else {
        s += "Resolution limit determined from a curve fit to the function "+
          radialfunction.formatfunction();
        if (nrej > 0) {
          s += ", after rejecting "+StringUtil::itos(nrej,2)+" bins";
        }
        s += "\n";
      }
    } else { // insufficient data
      s += "Insufficient data to determine resolution limit\n";
    }

    if (status == -1) {
      if (fitted) {
        s += "All scores are below the threshold";
      } else {
        s += "Fitted scores are below the threshold";
      }
      if (anomalous) {
        s += ", ie there is no apparent anomalous signal from CCanom\n";
      } else {
        s += ", ie the data are very poor even at the lowest resolution\n";
      }
    } else if (status == +1) {
      s += "All scores are above the threshold,";
      if (anomalous) {
        s += " ie there is a significant anomalous signal for all resolutions to the maximum of "+
        StringUtil::ftos(highres, 6, 2)+"A\n";
      } else {
        s += " ie data extends to the maximum resolution of "+
          StringUtil::ftos(highres, 6, 2)+"A\n";
      }
    } else if (status == 0) {
      if (anomalous) {
        s += "Estimate of the resolution limit for a significant anomalous signal"+
        StringUtil::ftos(highres, 6, 2)+"A";
      } else {
        s += "Estimate of the overall resolution limit "+
        StringUtil::ftos(highres, 6, 2)+"A";
      }
      if (fittype != NONE) {
        s += "\n        estimated from the point where the fit drops below threshold\n";
      }
    }
    return s;
  }
  // ------------------------------------------------------------
  std::string ResolutionLimit::formatbrief(const bool& anomalous) const
  // anomalous == true for assessment of anomalous signal (changes wording)
  // brief version
  {
    std::string s;
    if (status <= -2) {
      s += "No resolution limit determined";
      return s;
    }
    std::string sthreshold =  StringUtil::Strip(StringUtil::ftos(limit, 6, 2));

    if (sufficientdata) {  // sufficent data to determine limit
      if (status == -1) {
        s += "All scores are below the threshold "+sthreshold;
        if (anomalous) {
          s += ", ie there is no significant anomalous signal";
        } else {
          s += ", ie the data are very poor even at the lowest resolution";
        }
      } else if (status == +1) {
      s += "All scores are above the threshold "+sthreshold;
      if (anomalous) {
        s += ", ie there is a significant anomalous signal to the edge at "+
        StringUtil::ftos(highres, 6, 2)+"A";
      } else {
        s += " ie data extends to the maximum resolution of "+
          StringUtil::ftos(highres, 6, 2)+"A";
      }
    } else if (status == 0) {
        if (anomalous) {
          s += "Estimate of the resolution limit for a significant anomalous signal"+
            StringUtil::ftos(highres, 6, 2)+"A";
        } else {
          s += "Estimate of the overall resolution limit: "+
            StringUtil::ftos(highres, 6, 2)+"A";
        }
        if (fittype == NONE) { // no function fit
          s += ",  from the point at which the score drops below threshold"+sthreshold;
        } else {
          s += ", from the point where the fit drops below threshold "+sthreshold;
        }
      }
    } else { // insufficient data
      s += "Insufficient data to determine resolution limit";
    }
    return s;
  }
  // ------------------------------------------------------------
  std::string ResolutionLimit::formatparameters() const {
    if (fittype == LINEAR) {
      return "Linear fit: slope = " + StringUtil::ftos(slope,9,4)+
        ", intercept = " + StringUtil::ftos(intercept,9,4);
    } else {
      return radialfunction.format();
    }
  }
  // ------------------------------------------------------------
  bool ResolutionLimit::testinsufficientdata(const std::vector<double> score) const
  // return true if there are not enough data points to determine
  // maximum resolution
  {
    int nbins = score.size();
    if (nbins == 0) {return true;}
    bool insufficientdata = false;

    // list of bin numbers that are not zero
    std::vector<int> filledbins;
    for (int i=0;i<nbins;++i) {
      if (score[i] != 0.0) {
        filledbins.push_back(i);
      }
    }

    //^^
    //    for (int i=0;i<nbins;++i) {
    //      std::cout <<" "<<score[i];
    //    }
    //    std::cout <<std::endl;
    //^-

    //^std::cout <<"ResolutionLimit: nbins " <<nbins <<" nfilled "<<filledbins.size()<<"\n"; //^^

    if (filledbins.size() == 0) {
      return true; // no data
    }

    const size_t MINNUMBER = 4;
    if (filledbins.size() < MINNUMBER) {
      // too few datapoints
      //      std::cout <<"insufficientdata, Nfilledbins "<<filledbins.size()<<"\n"; //^
      insufficientdata = true;
    }

    // Fraction of bins with data
    double filledfraction = double(filledbins.size())/double(nbins);
    const double MINFILLEDFRACTION = 0.2;
    if (filledfraction < MINFILLEDFRACTION) {
      // too few datapoints
      //      std::cout <<"insufficientdata, filledfraction "<<filledfraction<<"\n"; //^
      insufficientdata = true;
    }

    // Last filled bin should not be too far from end
    double lastfraction = double(filledbins.back()+1)/double(nbins);
    double minlastfraction = std::min(0.9, double(nbins-1)/double(nbins)-0.001);
    if (lastfraction < minlastfraction) {
      //^std::cout <<"insufficientdata, lastfraction "<<lastfraction
      //^               <<" "<<minlastfraction<<"\n"; //^
      // but allow it if score has already fallen below threshold
      if (score[filledbins.back()] < limit) {
        insufficientdata = false;
      } else {
        insufficientdata = true;
      }
    }

    if (insufficientdata) {
      // not sufficient but have another go
      // Count number of bins above a good high threshold
      const double GOODCC = 0.9;
      int ngood = 0;
      for (int i=0;i<nbins;++i) {
        if (score[i] > GOODCC) {ngood++;}
      }
      if (ngood >= MINNUMBER) { // accept
        insufficientdata = false;
      }
    }
    return insufficientdata;
  }
  // ------------------------------------------------------------
  //! return fitted value at s = 1/d^2
  double ResolutionLimit::fitvalue(const double& s) const
  {
    if (fitted) {
      if (fittype == LINEAR) {
        return intercept + s * slope;
      } else {
        return radialfunction.value(s);
      }
    }
    return 0.0;
  }
  // ------------------------------------------------------------
}
