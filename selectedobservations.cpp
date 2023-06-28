// selectedobservations.cpp
//
//  Class to store and work with a subset of observations belonging to a reflection
//  eg from one dataset, or I+, I-
//

#include <algorithm>

#include "selectedobservations.hh"
#include "controls.hh"
#include "scala_util.hh"

#define ASSERT assert
#include <assert.h>

namespace scala
{
  bool SelectedObservations::sampleSD = false;   // static
  int SelectedObservations::minimumsample = 6;

  // ------------------------------------------------------------
  // constructor for selecting datasets & anomalous class
  SelectedObservations::SelectedObservations(const reflection& Refl,
                                             const int& datasetIndex,
                                             const AnomalousClass& Anomclass)
  {
    init(Refl, datasetIndex, Anomclass, WeightType::VARIANCE);
  }
  // ------------------------------------------------------------
  // constructor for selecting datasets & anomalous class, weight type
  SelectedObservations::SelectedObservations(const reflection& Refl,
                                             const int& datasetIndex,
                                             const AnomalousClass& Anomclass,
                                             const WeightType::AverageWeightType& weightType)
  {
    init(Refl, datasetIndex, Anomclass, weightType);
  }
  // ------------------------------------------------------------
  // Initialise, selecting datasets & anomalous class
  void SelectedObservations::init(const reflection& Refl,
                                  const int& datasetIndex,
                                  const AnomalousClass& Anomclass)
  {
    init(Refl, datasetIndex, Anomclass, WeightType::VARIANCE);
  }
  // ------------------------------------------------------------
  // Initialise, selecting datasets & anomalous class
  void SelectedObservations::init(const reflection& Refl,
                                  const int& datasetIndex,
                                  const AnomalousClass& Anomclass,
                          const WeightType::AverageWeightType& weightType)
  // if datasetIndex < 0 select all data
  // Note that this will only select "accepted" observations
  // Npart is number of parts to split into, if required
  {
    this_ref = &Refl;
    observation this_obs;
    int index;
    int PlusMinus;
    nobs = Refl.num_observations();
    ASSERT (nobs > 0);
    // Set size of internal vectors, for all observations in this
    // reflection including unselected ones
    use.assign(nobs, false);
    outliers.assign(nobs, 0);
    wI.resize(nobs);
    wj.assign(nobs, 0.0);
    wv.resize(nobs);
    delta.assign(nobs,0.0);
    deltaall.assign(nobs,0.0);
    part.assign(nobs,0);
    Nused = 0;

    Refl.reset(); // reset obs count
    // loop accepted observations
    while ((index = Refl.next_observation(this_obs)) >= 0) {
      if (datasetIndex < 0 ||
          this_obs.datasetIndex() == datasetIndex) {
        if ((Anomclass == IPLUS) || (Anomclass == IMINUS)) {
          PlusMinus = this_obs.Isym()%2;
          if (Anomclass == IPLUS && PlusMinus == 0) continue;
          if (Anomclass == IMINUS && PlusMinus != 0) continue;
        }
        use[index] = true; // set used flag
        Nused++;
      }
    }
    npart = -1;
    State = 0;
    weighttype = weightType;
    samplemean();  // calculate weights and mv
    Average();  // calculate average
    nextobs = -1;
    discrepant = false;
  }
  // ------------------------------------------------------------
  void SelectedObservations::SetNpart(const int& Npart)
  // Divide into Npart parts: Npart = 2 or 4
  {
    ASSERT (Npart == 2 || Npart == 4);
    if (Npart != npart) {  // only do this if changed
      npart = Npart;
      part.assign(nobs,0);
      std::vector<std::pair<int,float> > idxrandom;
      for (int index=0;index<nobs;++index) {
        part[index] = 0;
        if (use[index]) {
          // index to this observation & a random number between 0 & 1
          idxrandom.push_back(std::pair<int,float>(index, FRandom(1.0)));
        }
      }
      // sort into order on random number
      std::sort(idxrandom.begin(), idxrandom.end(), CompareIFpair);

      int n2 = Nused/2;
      int n1 = n2/2;
      int n3 = n1+n2;

      for (size_t i=0;i<idxrandom.size();++i) {
        //  index for random part-datasets, 0 -> Npart-1
        int index = idxrandom[i].first;
        if (npart == 2) {
          if (int(i) >= n2) part[index] = 1;
        } else {
          // npart == 4
          if (int(i) >= n3) {part[index] = 3;}
          else if (int(i) >= n2) {part[index] = 2;}
          else if (int(i) >= n1) {part[index] = 1;}
        }
        ASSERT (part[index] < npart);  // check for now
      }
    }
  }
  // ------------------------------------------------------------
  // Next used & accepted observation, returns -1 if end
  int SelectedObservations::next_observation(observation& obs) const
  {
    while (++nextobs < nobs) {
      if (use[nextobs]) {
        obs = this_ref->get_observation(nextobs);
        return nextobs;
      }
    }
    nextobs = -1;
    return nextobs;
  }
  // ------------------------------------------------------------
  bool SelectedObservations::HalfAverages(float& I1, float& I2)
  // Get average I for each random part
  //  return false unless both are present
  {
    if (Nused < 2) return false;
    SetNpart(2);
    bool gotBoth = true;
    I1 = AveragePart(0).I();
    if (I1 == 0.0) gotBoth = false;
    I2 = AveragePart(1).I();
    if (I2 == 0.0) gotBoth = false;
    // Reset state to force average calculation next time
    State = 0;
    return gotBoth;
  }
  // ------------------------------------------------------------
  bool SelectedObservations::HalfAveragesSigI(IsigI& I1sig, IsigI& I2sig)
  // Get average I, sigI for each random part
  //  return false unless both are present
  {
    if (Nused < 2) return false;
    SetNpart(2);
    bool gotBoth = true;
    I1sig = AveragePart(0);
    if (I1sig.sigI() <= 0.0) gotBoth = false;
    I2sig = AveragePart(1);
    if (I2sig.sigI() <= 0.0) gotBoth = false;
    // Reset state to force average calculation next time
    State = 0;
    return gotBoth;
  }
  // ------------------------------------------------------------
  bool SelectedObservations::PartAverages(std::vector<float>& Is)
  // Get average I for each random part
  //  return false unless all are present
  {
    if (Nused < npart) return false;
    bool gotall = true;
    Is.resize(npart);
    for (int i =0;i<npart;++i) {
      Is[i] = AveragePart(i).I();
      if (Is[i] == 0.0) gotall = false;
    }
    // Reset state to force average calculation next time
    State = 0;
    return gotall;
  }
  // ------------------------------------------------------------
  //! Set weight
  void SelectedObservations::SetWeight(const WeightType::AverageWeightType& weightType)
  {
    weighttype = weightType;
    State = 0;
    Average();  // recalculate average with new weights
  }
  // ------------------------------------------------------------
  double SelectedObservations::Weight(const double& sd, const double& g) const
  // Return weight calculated from sd (scaled) or g according to weighttype
  {
    if (weighttype == WeightType::VARIANCE) {
      return 1.0f/(sd*sd);
    } else if (weighttype == WeightType::SQRTSCALE) {
      if (g <= 0.0f) return 0.0f;
      // a strong observation has small scale ie large g and high weight
      return sqrt(g);
    } else if (weighttype == WeightType::SCALE) {
      if (g <= 0.0f) return 0.0f;
      // a strong observation has small scale ie large g and high weight
      return g;
    }
    return 1.0f;
  }
  // ------------------------------------------------------------
  //! Set sample variance, minimum number of values (<0 to switch off), static
  void SelectedObservations::SetSampleSD (const int& minSample)
  {
    sampleSD = false;
    if (minSample > 0) {
      sampleSD = true;
      minimumsample = minSample;
    }
  }
  // ------------------------------------------------------------
  IsigI SelectedObservations::Average()
  // Weight depends on weighttype
  //   UNIT        unit weights
  //   VARIANCE    weight = 1/variance
  //   SQRTSCALE   weight = sqrt(g)  g = 1/scale
  //   SCALE       weight = g        g = 1/scale
  //  Note that a smaller scale = larger g = larger weight
  {
    //  <I> = Sum(w I') / Sum (w)
    if (State == +1) {return avIsigI;}
    if (Nused == 0) {return IsigI(0.0,0.0);}

    sumwI = 0.0;
    sumwj = 0.0;
    sumwv = 0.0;  // sum(1/v) or sum(w^2 v)
    sdIsample.assign(nobs, 0.0);

    Nused = 0;

    double w;
    IsigI Is;
    Rtype sd;
    double Vs; // sample variance
    sampleSDused = false;
    sdIs = 0.0;

    for (int i=0;i<nobs;i++)  {
      if (use[i]) {
        Nused++;
        // scaled I', sigI'
        Is = this_ref->get_observation(i).kI_sigI();
        sd = Is.sigI();
        if (sd <= 0.0) {
          std::cout << sd << " !!SD\n";
        }
        ASSERT (sd > 0.0);
        w = wj[i];
        wI[i] = w * Is.I(); // w I'
        sumwI += wI[i];
        sumwj += wj[i];
        if (weighttype == WeightType::VARIANCE) {
          wv[i] = w;  // 1/sd^2
        } else {
          wv[i] = w * w * sd *sd;
        }
        sumwv += wv[i];  // sum for var(Imean)
      }
    }
    if (Nused > 0) {
      // always get sdI
      if (weighttype == WeightType::VARIANCE) {
        sdI = sqrt(1.0/sumwv);
      } else {
        sdI = sqrt(sumwv)/sumwj;
      }
      if (sampleSD && (Nused > minimumsample)) {
        Vs = mv.SampleVariance();
        sdIs = sqrt(Vs/double(Nused)); //  SD(<I>)
        sampleSDused = true;
        //^
        //        std::cout << "Vs, SDmn, wSD, mnI, N, sumwj, sumwv " <<Vs<<" "<<sdI<<
        //          " "<<sqrt(1.0/sumwj)<<
        //          " "<<sumwI/sumwj<<" "<< Nused<<" "<<sumwj<<" "<<sumwv<<"\n";
        //      if (sdI > 10000.) {
        //        std::cout <<"Large SDs, sampleSD " << mv.SampleSD() <<"\n";
        //        for (int i=0;i<nobs;i++)  {
        //          if (use[i]) {
        //            std::cout << "I, sigI, gscale "
        //                      << this_ref->get_observation(i).kI()<<" "
        //                      << this_ref->get_observation(i).ksigI()<<" "
        //                      <<this_ref->get_observation(i).Gscale() <<"\n";
        //          }
        //        }
        //      }
        //^-
        avIsigI = IsigI(sumwI/sumwj, sdIs);
      } else {
        avIsigI = IsigI(sumwI/sumwj, sdI);
      }
      State = +1;

      if (sampleSDused && (sumwj > 0.0)) {
        // individual SDs from sample sdI
        double fac = Vs * sumwj / double(Nused);
        double totui = 0.0; //^ sanity check
        double totv = 0.0;
        for (int i=0;i<nobs;i++)  {
          if (use[i] && wj[i] > 0.0) {
            // Note wj[i] = w, sumwj = Sum(wj)
            // var(I[i]) = var(<I>)/u[i] where Sum(u[i]) = 1,
            // var(<I>) = Vs/N  Vs = sample variance
            // u[i] = wj[i]/sumwj
            // ie var(I[i]) = Vs sumwj / N wj[i]
            // fac = Vs sumwj / N
            //    sd(I[i]) = sqrt(fac / wj[i])
            sdIsample[i] = sqrt(fac/wj[i]);
            totv += fac/wj[i];
            totui += wj[i]/sumwj;
          }
        }
        //      std::cout << "Sanity checks "
        //                <<sdIs<<" "<<Vs/double(Nused)<<" "<< totv/double(Nused) <<" "<<totui<<"\n";
      }
    } else {
      avIsigI = IsigI(0.0,0.0);
      State = -1;
    }
    return avIsigI;
  }
  // ------------------------------------------------------------
  void SelectedObservations::samplemean()
  // sets mv, and wj
  {
    mv.clear();
    double w;
    double g;
    IsigI Is;
    Rtype sd;
    for (int i=0;i<nobs;i++)  {
      if (use[i]) {
        g = this_ref->get_observation(i).Gscale();
        // scaled I', sigI'
        Is = this_ref->get_observation(i).kI_sigI();
        sd = Is.sigI();
        ASSERT (sd > 0.0);
        w = Weight(sd, g);   // weight according to weighttype
        wj[i] = w;
        // scaled I', sigI'
        Is = this_ref->get_observation(i).kI_sigI();
        // always add into sums for sampleSD
        mv.Add(double(Is.I()), w);
      }
    }
  }
  // ------------------------------------------------------------
  //! Var(<I>) from sample variance, unconditional from Average()
  //  just samplevariance/N
  double SelectedObservations::Variance() const
  {
    return mv.VarianceofMean();
  }
  // ------------------------------------------------------------
  //! Sample variance, unconditional from Average()
  double SelectedObservations::sampleVariance() const
  {
    return mv.SampleVariance();
  }
  // ------------------------------------------------------------
  IsigI SelectedObservations::AveragePart(const int& WhichPart)
  // If WhichPart >= 0, use only selected random part
  //              < 0  use all accepted
  {
    //  <I> = Sum(w g I) / Sum (w g^2)
    if (Nused == 0) {return IsigI(0.0,0.0);}

    sumwI = 0.0;
    sumwj = 0.0;
    sumwv = 0.0;

    IsigI Is;
    Rtype sd;

    double w;
    double g;

    int Nu = 0;
    MeanVariance mvp;

    for (int i=0;i<nobs;i++)  {
      if (use[i]) {
        if (WhichPart < 0 || part[i] == WhichPart) {
          Nu++;
          // scaled I', sigI'
          Is = this_ref->get_observation(i).kI_sigI();
          sd = Is.sigI();
          ASSERT (sd > 0.0);
          w = wj[i];
          wI[i] = w * Is.I(); // w I'
          sumwI += wI[i];
          sumwj += wj[i];
          if (weighttype == WeightType::VARIANCE) {
            wv[i] = w;  // 1/sd^2
          } else {
            wv[i] = w * w * sd *sd;
          }
          sumwv += wv[i];  // sum for var(Imean)
          if (sampleSD) {
            mvp.Add(wI[i], w);
          }
        }
      }
    }
    IsigI avIsigIpart;
    if (Nu > 0) {
      // always get sdI
      if (weighttype == WeightType::VARIANCE) {
        sdI = sqrt(1.0/sumwv);
      } else {
        sdI = sqrt(sumwv)/sumwj;
      }
      if (sampleSD && (Nused > minimumsample)) {
        sdIs = mvp.SDofMean(); // SD of mean from sample variance
        avIsigIpart = IsigI(sumwI/sumwj, sdIs);
      } else {
        avIsigIpart = IsigI(sumwI/sumwj, sdI);
      }
    } else {
      avIsigIpart = IsigI(0.0,0.0);
    }
    return avIsigIpart;
  }
  // ------------------------------------------------------------
  // List of deviations delta (ie delI/sigma(I) ) where delI
  //  is difference from mean of other observations
  //   returns delta(NobsRefl), unused slots set = 0.0 ie not closed down
  //   Rejected observations are left with their original delta,
  //     accepted ones are reevaluated after rejection
  //   delta.size() = total number of observations in reflection
  //  If (fromSample && sampleSDused), then use sample variance, else
  //   use individual variance
  //
  std::vector<float> SelectedObservations::Deviations(const bool& fromSample)
  {
    if (Nused <= 0) {return delta;}
    if (State == 0) Average();
    if (State >= 2) {return delta;}  // already calculated
    if (Nused > 1) {
      // we need at least 2 observations
      float Iothers;
      float varothers;
      float vv;
      // true iff both sampleSDused && fromSample true
      bool usesamplesd = (sampleSDused && fromSample);
      //      bool DEBUG = true;
      bool DEBUG = false;
      //^^
      //      for (size_t k=0; k<outliers.size(); k++) {
      //        if (outliers[k]) {DEBUG = true;}
      //      }
      //      std::cout <<"Dev "<<this_ref->hkl().format()
      //                <<" "<<sampleSDused<<" "<<fromSample<<" "
      //                <<usesamplesd <<std::endl;

      if (DEBUG) {
        std::cout << "SelectedObservations::Deviations() "<<this_ref->hkl().format()
                  <<" "<<avIsigI.I()<<" "<<sdI<<" "<<sdIs<<"\n";
      }
      //^
      double sum = 0.0;
      int n = 0;
      //^-
      std::vector<IvarI> mnothers = MeanIothers(); // mean of other observations

      delta.assign(nobs, 0.0);


      for (int i=0;i<nobs;i++) {
        if (use[i]) {
          // <I>(others)  ie excluding this observation
          //  and its variance (scaled to this observation)
          varothers = mnothers[i].varI();
          Iothers = mnothers[i].I();
          if (usesamplesd) {
            vv = varothers * (sumwj - wj[i]) / (float(Nused-1) * wj[i]);
          } else {
            vv = this_ref->get_observation(i).ksigI()*
              this_ref->get_observation(i).ksigI();
            //?//vv += varothers;
          }
          //^
          //      if (vv <= 0.0) {
          //        std::cout << "varothers " << varothers<<" "<<vv <<std::endl;
          //      }
          if (vv <= 0.0) {
            std::cout << vv << " !!VV\n";
          }
          ASSERT (vv > 0.0);
          delta[i] = (this_ref->get_observation(i).kI() - Iothers)/sqrt(vv);
          if (DEBUG) {
            std::string s = "  ";
            if (std::abs(delta[i]) > 3.0) {s = " *";}
            printf("kI, ksigI, Ioth, varOth, varI, delI, sd, delta %8.1f %8.1f %8.1f %8.1f %8.1f %8.3f %8.3f %8.4f %s\n"
                   , this_ref->get_observation(i).kI()
                   , this_ref->get_observation(i).ksigI()
                   , Iothers ,varothers,vv-varothers,
                   (this_ref->get_observation(i).kI() - Iothers), sqrt(vv),
                   delta[i], s.c_str());
            sum += delta[i]*delta[i];
            n++;
          }
        }
      }
      if (DEBUG) {
        std::cout <<"MeanD^2 " << sum/double(n) <<"\n";
      }
    }

    State = +2;
    return delta;
  }
  // ------------------------------------------------------------
  // For each observation, return mean of other observations, scaled to each observation
  //   returns mnothers(NobsRefl), unused slots set = 0.0 ie not closed down
  //
  std::vector<IvarI> SelectedObservations::MeanIothers()
  {
    std::vector<IvarI> mnothers(nobs, IvarI(0.0,0.0));
    if (Nused <= 0) return mnothers;
    if (State == 0) Average();
    if (Nused > 1) {
      // we need at least 2 observations
      double wjothers, wvothers;

      for (int i=0;i<nobs;i++) {
        if (use[i]) {
          // <I>(others)  ie excluding this observation
          //  and its variance
          const double MINW = 1.0e-30;
          wjothers = Max(sumwj - wj[i], MINW); // trap very small wj for rounding errors
          wvothers = Max(sumwv - wv[i], MINW); // trap very small wv for rounding errors
          mnothers[i].I() = (sumwI - wI[i])/wjothers;
          if (sampleSDused) {
            mnothers[i].varI() =
              mv.VarianceofMeanOmit1(this_ref->get_observation(i).kI(), wj[i]);
          } else {
            if (weighttype == WeightType::VARIANCE) {
              mnothers[i].varI() = 1.0  / wvothers;
            } else {
              mnothers[i].varI() = wvothers/(wjothers*wjothers);
            }
          }
        }
      }
    }
    State = +1;
    return mnothers;
  }
  // ------------------------------------------------------------
  std::vector<float> SelectedObservations::DeltaAll() const
  {
    std::vector<float> delall = delta;
    for (size_t i=0;i<delta.size();++i) {
      if (delta[i] == 0.0) {
	delall[i] = deltaall[i];
      }
    }
    return delall;
  }
  // ------------------------------------------------------------
  std::vector<float> SelectedObservations::Delta2(const bool& fromSample)
  // List of deviations delta2 (ie delI/sigma(I) ) where delI
  //  is difference from mean of all observations
  //   returns delta2(NobsRefl), unused slots set = 0.0 ie not closed down
  //   delta2.size() = total number of observations in reflection
  //  If (fromSample && sampleSDused), then use sample variance, else
  //   use individual variance
  //
  {
    std::vector<float> delta2(nobs,0.0);
    if (Nused <= 0) return delta2;
    if (State == 0) Average();
    if (Nused > 1) {
      // true iff both sampleSDused && fromSample true
      bool usesamplesd = (sampleSDused && fromSample);
      float sigmai;
      float fac = sqrt(float(Nused)/(Nused-1));
      for (int i=0;i<nobs;i++) {
        if (use[i]) {
          float delI = (this_ref->get_observation(i).kI() - avIsigI.I());
          if (usesamplesd) {
            sigmai = sdIsample[i];
          } else {
            sigmai = this_ref->get_observation(i).ksigI();
          }
          //^^
          //      observation obs = this_ref->get_observation(i);
          //      std::cout << fac<<" "<<delI<<" "<<sigmai
          //                <<" "<<obs.Gscale()<<" "<<obs.sigI()<<" "<<obs.ksigI()
          //                <<" "<<obs.hkl_original().format()
          //                <<" fac, delI, sigmaI, g, sigI, ksigI\n"; //^^
          delta2[i] = fac * delI/sigmai;
          //      std::cout <<"Delta2, fac, delI, sigmai, I, avI "<<
          //        delta2[i]<<" "<<fac<<" "<<delI<<" "<<sigmai<<
          //        " "<<this_ref->get_observation(i).kI()<<" "<<
          //        " "<<avIsigI.I()<<"\n";
        }
      }
    }
    return delta2;
  }
  // ------------------------------------------------------------
  std::vector<float> SelectedObservations::Delta3()
  // List of deviations delta3 (ie delI/SDsample(I) ) where delI
  //  is difference from mean of all observations
  // Only relevant if SDsample is used
  //   returns delta3(NobsRefl), unused slots set = 0.0 ie not closed down
  //   delta3.size() = total number of observations in reflection
  {
    std::vector<float> delta3(nobs,0.0);
    if (Nused <= 0) {return delta3;}
    if (State == 0) Average();
    if (Nused > 1 && SampleSDused()) {
      float fac = sqrt(float(Nused)/(Nused-1));
      for (int i=0;i<nobs;i++) {
        if (use[i]) {
          float delI = (this_ref->get_observation(i).kI() - avIsigI.I());
          delta3[i] = fac * delI/sdIsample[i];
        }
      }
    }
    return delta3;
  }
  // ------------------------------------------------------------
  // List of delI (scaled)
  //   returns delI(nobs), unused slots set = 0.0 ie not closed down
  std::vector<float> SelectedObservations::DelI()
  {
    std::vector<float> delI = std::vector<float>(nobs,0.0);
    if (Nused <= 0) return delI;
    if (State == 0) Average();
    if (Nused > 1) {
      // we need at least 2 observations
      for (int i=0;i<nobs;i++) {
        if (use[i]) {
          delI[i] = (this_ref->get_observation(i).kI() - avIsigI.I());
        }
      }
    }
    return delI;
  }
  // ------------------------------------------------------------
  // List of sigma(I) (scaled)
  //   returns sigmaI(nobs), unused slots set = 0.0 ie not closed down
  std::vector<float> SelectedObservations::sigmaI() const
  {
    std::vector<float> sigmai = std::vector<float>(nobs,0.0);
    if (Nused <= 0) return sigmai;
    for (int i=0;i<nobs;i++) {
      if (use[i]) {
        sigmai[i] = this_ref->get_observation(i).ksigI();
      }
    }
    return sigmai;
  }
  // ------------------------------------------------------------
  // List of I, sigma(I) (scaled)
  //   returns I,sigI(nobs), unused slots set = 0.0 ie not closed down
  std::vector<IsigI> SelectedObservations::IsigIlist() const
  {
    std::vector<IsigI> isigi(nobs);
    if (Nused <= 0) return isigi;
    for (int i=0;i<nobs;i++) {
      if (use[i]) {
        isigi[i] = this_ref->get_observation(i).kI_sigI();
      }
    }
    return isigi;
  }
  // ------------------------------------------------------------
  int SelectedObservations::Outliers(const RejectFlags& rejflags)
  //  Reject outliers
  //
  // On entry:
  // from rejflags::
  //  sdrej   rejection limit (sds) for more than 2 observations
  //  sdrej2  rejection limit (sds) for 2 observations
  //  rej2policy    what to do with two observations
  //           = REJECT           reject both
  //           = KEEP             keep both
  //           = REJECTLARGER     reject larger
  //           = REJECTSMALLER    reject smaller
  /*
      Algorithm:

      (1) if there are 2 observations (left), then

      (a) for each observation Ihl, test deviation
      Delta(hl) = (Ihl - ghl Iother) / sqrt[sigIhl^2 + (ghl*sdIother)^2]
      against sdrej2, where Iother = the other observation

      (b) if either Delta(hl) > sdrej2, then reject reflection
      according to flag rej2policy (qv)

      (2) if there 3 or more observations left, then

      (a) for each observation Ihl,

      (i)   calculate weighted mean of all other
      observations <I>n-1 & its sd(<I>n-1)
      (ii)  deviation Delta(hl) =
      (Ihl - ghl <I>n-1>) / sqrt[sigIhl^2 + (ghl*sd(<I>n-1))^2]
      (iii) find largest deviation max|Delta(hl)|
      (iv)  count number of observations for which
      Delta(hl) .ge. 0 (ngt), & for which
      Delta(hl) .lt. 0 (nlt)

      (b) if max|Delta(hl)| > sdrej, then reject one observation, but
      which one?
      (i) if ngt == 1 .or. nlt == 1, then one observation is a long
      way from the others, and this one is rejected
      (ii) else reject the one with the worst deviation max|Delta(hl)

      (3) iterate from beginning
  */
  {
    bool DEBUG = false;
    float sdrej  = rejflags.sdrej;
    float sdrej2 = rejflags.sdrej2;
    RejectFlags::Reject2Policy rej2policy = rejflags.rej2policy;
    //  enum Reject2Policy {REJECT, KEEP, REJECTLARGER, REJECTSMALLER};

    int Nrej = 0;
    if (Nused <= 0) {
      Deviations();
      return Nrej;
    }
    if (State == 0) Average();

    float I;
    int lsmaller, llarger;
    bool first = true;

    // Loop until only one observation if necessary
    // ie iterative outlier rejection
    while (Nused > 1) {
      //recalculate deviations if required
      if (State < +2) {
	Deviations();
	if (first) {
	  // save initial deviations, for printing
	  deltaall = delta;
	  first = false;
	}
      }
      if (Nused == 2) {
        //  -  -  -  -  -  -  two observations -  -  -  -  -  -
        float smaller = 1.0e+20;
        llarger = -1;
        bool reject = false;
        for (int i=0;i<nobs;i++) {
          if (use[i]) {
            if (std::abs(delta[i]) > sdrej2) reject = true;
            ASSERT (this_ref->get_observation(i).Gscale() != 0.0);
            I = this_ref->get_observation(i).I()/
              this_ref->get_observation(i).Gscale();
            if (I < smaller) {
              smaller = I;
              lsmaller = i;
              if (llarger < 0) llarger = i;
            } else
              {llarger = i;}
          }
        }
        //---
        if (reject) {
          discrepant = true;
          if (rej2policy == scala::RejectFlags::REJECT) {
            use[lsmaller] = false;
            use[llarger] = false;
            outliers[lsmaller] = +1;
            outliers[llarger] = +1;
            Nrej += 2;
            Nused -= 2;
            //^+
            if (DEBUG) std::cout << "RejBoth " << delta[lsmaller]
                                 << " " << delta[llarger]
                                 << "  hkl " << this_ref->hkl().format() << "\n";

          } else if (rej2policy == scala::RejectFlags::REJECTSMALLER) {
            use[lsmaller] = false;
            outliers[lsmaller] = +1;
            Nrej++;
            Nused--;
            //^+
            if (DEBUG) std::cout << "RejSmaller " << delta[lsmaller]
                                 << "  hkl " << this_ref->hkl().format() << "\n";
          } else if (rej2policy == scala::RejectFlags::REJECTLARGER) {
            use[llarger] = false;
            outliers[llarger] = +1;
            Nrej++;
            Nused--;
            //^+
            if (DEBUG) std::cout << "RejLarger " << delta[llarger]
                                 << "  hkl " << this_ref->hkl().format() << "\n";
          } else if (rej2policy == scala::RejectFlags::KEEP) {
            // keep both but flag them
            outliers[lsmaller] = -1;
            outliers[llarger] = -1;
          }
          State = 0;  // flag to force recalculation of deviations
          if (rej2policy == scala::RejectFlags::KEEP) {break;}
        } else {
          break;}  // no rejection
      } else if (Nused > 2) {
        //  -  -  -  -  -  -  three or more observations -  -  -
        float worst = -1.0e+10;
        int ngt = 0;
        int nlt = 0;
        int lworst, lgt, llt;
        // Find worst
        for (int i=0;i<nobs;i++) {
          if (use[i]) {
            // use deviation from "others" as rejection criterion
            if (std::abs(delta[i]) > worst)
              {worst = std::abs(delta[i]); lworst=i;}

          }
        }
        if (worst > sdrej) {
          // Reject one but which one?
          //  If one is on the opposite side of the mean deviation to all others,
          //   reject that one
          //  else reject worst
          discrepant = true;
          // Get mean deviation
          double devmean = 0.0;
          for (int i=0;i<nobs;i++) {
            if (use[i]) {
              devmean += delta[i];
            }
          }
          devmean /= double(Nused);
          for (int i=0;i<nobs;i++) {
            if (use[i]) {
              // count positives & negative difference from mean deviation
              //  find largest
              if ((delta[i]-devmean) >= 0.0) {
                ngt++; lgt=i;
              } else {
                nlt++; llt=i;
              }
            }
          }
          if (nlt == 1) {lworst = llt;}
          else if (ngt == 1) {lworst = lgt;}
          // reject the unique one on one side of the mean deviation, or the worst
          use[lworst] = false;
          outliers[lworst] = true;
          Nrej++;
          Nused--;
          State = 0;
          //^+
          if (DEBUG) {
            std::cout << "Reject " << delta[lworst]
                      << "  hkl " << this_ref->hkl().format() << "\n";
          }
        } else {
          break;}  // no rejection
      } else {break;}   // no rejection
    }
    // recalculate deviations and sample mean if required
    if (State < +2) {
      samplemean();
      Deviations();
    }
    State = +3;
    return Nrej;
  }
  // ------------------------------------------------------------
  int SelectedObservations::Run(const int& kobs) const
  // Run number for kobs'th observation
  // note that kobs is the index into the whole list, including unused ones
  {
    return this_ref->get_observation(kobs).run();
  }
  // ------------------------------------------------------------
  bool SelectedObservations::Full(const int& kobs) const
  // True if kobs'th observation is Fully recorded
  {
    return this_ref->get_observation(kobs).IsFull();
  }
  // ------------------------------------------------------------
  bool SelectedObservations::AnyFull() const
  // True if any observations are Fully recorded
  {
    for (int i=0;i<nobs;++i) {
      if (this_ref->get_observation(i).IsFull()) return true;
    }
    return false;
  }
  // ------------------------------------------------------------
  std::vector<int> SelectedObservations::OutlierIndexList
  (const RejectFlags& rejflags)
  // Return list of (index number+1) for each outlier observation, if any,
  //   or negated for Deviant but not rejected
  // Calls "Outliers" first
  // only called from RejectList::Check
  {
    Outliers(rejflags);
    std::vector<int> idxlist;
    for (int i=0;i<nobs;++i) {
      if (outliers[i] > 0) {
        // reject
        idxlist.push_back(i+1);
      } else if (outliers[i] < 0) {
        // deviant
        idxlist.push_back(-i-1);
      }
    }
    return idxlist;
  }
  // ------------------------------------------------------------
  double SelectedObservations::chiSq(const bool& fromSample)
  // Mean Chi^2, goodness of fit, Mean((I-<I>)/sigma)
  //  if (fromSample && sampleSDused), then use sample variance
  //  if (fromSample && !sampleSDused), return 0.0
  // else  use individual variance
  {
    MeanValue chisq;
    std::vector<float> delta;
    State = 0; // force recalculation
    if (!fromSample) {
      //      State = std::min(State, +1);  // force recalculation
      delta = Delta2(fromSample);
    } else {
      if (sampleSDused) {
        delta = Delta2(fromSample);
      } else {
        return 0.0;
      }
    }
    if (delta.size() == 0.0) {return 0.0;}
    double sum = 0.0;
    int n = 0;
    for (size_t k=0; k<delta.size(); k++) {
      if (delta[k] != 0.0) {
        //^
        //        std::string s = " * ";
        //        if (fromSample) {s = " $ ";}
        //      std::cout << "delta " <<s<< delta[k]<<" "
        //                <<delta[k]*delta[k] <<std::endl; //^^
        sum += delta[k]*delta[k];
        n++;
        //^-
        chisq.Add(delta[k]*delta[k], 1.0f);
      }
    }
    //    std::cout << "Mean I, Mean D^2 " << avIsigI.I()<<" "
    //        << sum/double(n) << std::endl;
    return chisq.Mean();
  }
  // ------------------------------------------------------------
  std::vector<Iwgrp> SelectedObservations::iwgroup() const
  {
    std::vector<Iwgrp> iwb(Nused);
    size_t k = 0;
    for (int i=0;i<nobs;i++) {
      if (use[i]) {
        iwb[k++] = Iwgrp(this_ref->get_observation(i).kI_sigI().I(),
                         wj[i], this_ref->get_observation(i).Batch());
        //      std::cout << "iwgroup I, s "<<
        //        this_ref->get_observation(i).kI_sigI().I() <<" "<<
        //        this_ref->get_observation(i).kI_sigI().sigI() <<"\n";
      }
    }
    iwb.resize(k);
    return iwb;
  }
}
