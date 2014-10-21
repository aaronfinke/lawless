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
    // Set size of internal vectors, for all observations in this
    // reflection including unselected ones
    use.assign(nobs, false);
    outliers.assign(nobs, false);
    wgI.resize(nobs);
    wg2.resize(nobs);
    delta.assign(nobs,0.0);
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
    Average();  // calculate average
    nextobs = -1;
    discrepant = false;
  }
  // ------------------------------------------------------------
  bool CompareIFpair(const std::pair<int,float>& p1,const std::pair<int,float>& p2)
  {return (p1.second < p2.second);}
  // ------------------------------------------------------------
  void SelectedObservations::SetNpart(const int& Npart)
  // Divide into Npart parts: Npart = 2 or 4
  {
    ASSERT (Npart == 2 || Npart == 4);
    if (Npart != npart) {  // only do this if changed
      npart = Npart;
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
    Average();  // recalculate average with new weights
  }
  // ------------------------------------------------------------
  Rtype SelectedObservations::Weight(const Rtype& sd, const Rtype& g) const
  // Return weight calculated from val according to weighttype
  {
    if (weighttype == WeightType::VARIANCE) {
      return 1.0f/(sd*sd);
    } else if (weighttype == WeightType::SQRTSCALE) {
      if (g <= 0.0f) return 0.0f;
      // a strong observation has small scale ie large g and high weight not this      return 1.0f/sqrt(g);
      return sqrt(g);
    }
    return 1.0f;
  }
  // ------------------------------------------------------------
  //! Set sample variance, minimum number of values (<0 to switch off)
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
  //   SQRTSCALE   weight = 1/sqrt(g)  g = 1/scale
  //  Note that a smaller scale = larger g = larger weight
  {
    //  <I> = Sum(w g I) / Sum (w g^2)
    if (State == +1) {return avIsigI;}
    if (Nused == 0) {return IsigI(0.0,0.0);}


    sumwgI = 0.0;
    sumwg2 = 0.0;

    double w;
    double g;

    Nused = 0;
    MeanVariance mv;

    for (int i=0;i<nobs;i++)  {
      if (use[i]) {
        Nused++;
        g = this_ref->get_observation(i).Gscale();
        Rtype sd = this_ref->get_observation(i).sigI();
	//	if (sd <= 0.0) {
	//	  observation obs = this_ref->get_observation(i);
	//	  std::cout <<"SelObsAverage: "<<g<<" "<<sd
	//		    <<" "<<obs.hkl_original().format()<<" "
	//		    <<obs.Batch()<<std::endl;
	//	} //^^-
        ASSERT (sd > 0.0);
        w = Weight(sd, g);   // weight according to weighttype
        wgI[i] = w * g * this_ref->get_observation(i).I();
        sumwgI += wgI[i];
        wg2[i] = w * g * g;
        sumwg2 += wg2[i];
        if (sampleSD) {
          mv.Add(double(this_ref->get_observation(i).kI()), w/(g*g));
        }
      }
    }
    if (Nused > 0) {
      if (sampleSD && (Nused > minimumsample)) {
        sdI = mv.SD(); // sample SD
        //^
        //      std::cout << "SampleSD, wSD, I, N " <<sdI<<" "<<sqrt(1.0/sumwg2)<<
        //        " "<<sumwgI/sumwg2<<" "<< Nused<<"\n";
        //      if (sdI > 10000.) {
        //        std::cout <<"Large SDs, sampleSD " << mv.SampleSD() <<"\n";
        //        for (int i=0;i<nobs;i++)  {
        //          if (use[i]) {
        //            std::cout << "I, sigI, gscale "
        //              << this_ref->get_observation(i).kI()<<" "
        //              << this_ref->get_observation(i).ksigI()<<" "
        //              <<this_ref->get_observation(i).Gscale() <<"\n";
        //          }
        //        }
        //      }
        //^-
      } else {
        sdI = sqrt(1.0/sumwg2);
      }
      avIsigI = IsigI(sumwgI/sumwg2, sdI);
      State = +1;
    } else {
      avIsigI = IsigI(0.0,0.0);
      State = -1;
    }
    return avIsigI;
  }
  // ------------------------------------------------------------
  IsigI SelectedObservations::AveragePart(const int& WhichPart)
  // If WhichPart >= 0, use only selected random part
  //              < 0  use all accepted
  {
    //  <I> = Sum(w g I) / Sum (w g^2)
    if (Nused == 0) {return IsigI(0.0,0.0);}

    sumwgI = 0.0;
    sumwg2 = 0.0;

    double w;
    double g;

    int Nu = 0;
    MeanVariance mv;

    for (int i=0;i<nobs;i++)  {
      if (use[i]) {
        if (WhichPart < 0 || part[i] == WhichPart) {
          Nu++;
          g = this_ref->get_observation(i).Gscale();
          w = Weight(this_ref->get_observation(i).sigI(), g);
          wgI[i] = w * g * this_ref->get_observation(i).I();
          sumwgI += wgI[i];
          wg2[i] = w * g * g;
          sumwg2 += wg2[i];
          mv.Add(wgI[i], w);
        }
      }
    }
    IsigI avIsigIpart;
    if (Nu > 0) {
      if (sampleSD && (Nused > minimumsample)) {
        sdI = mv.SD(); // sample SD
      } else {
        sdI = sqrt(1.0/sumwg2);
      }
      avIsigIpart = IsigI(sumwgI/sumwg2, sdI);
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
  //
  std::vector<float> SelectedObservations::Deviations()
  {
    if (Nused <= 0) return delta;
    if (State == 0) Average();
    if (Nused > 1) {
      // we need at least 2 observations
      float Iothers;
      float varothers;

      std::vector<IsigI> mnothers = MeanIothers(); // mean of other observations

      for (int i=0;i<nobs;i++) {
        if (use[i]) {
          // <I>(others)  ie excluding this observation
          //  and its variance (scaled to this observation)
          varothers = mnothers[i].sigI() * mnothers[i].sigI();
          Iothers = mnothers[i].I();
          float vv = this_ref->get_observation(i).sigI()*
            this_ref->get_observation(i).sigI() +
            varothers;
          if (!(vv > 0.0)) {
            std::cout << "Aaargh " << vv << " "
                      << this_ref->get_observation(i).sigI()
                      << " " << varothers << "\n";
          }
          ASSERT (vv > 0.0);
          delta[i] = (this_ref->get_observation(i).I() - Iothers)/
            sqrt(this_ref->get_observation(i).sigI()*
                 this_ref->get_observation(i).sigI() +
                 varothers);
        }
      }
    }
    State = +2;
    return delta;
  }
  // ------------------------------------------------------------
  // For each observation, return mean of other observations, scaled to each observation
  //   returns mnothers(NobsRefl), unused slots set = 0.0 ie not closed down
  //
  std::vector<IsigI> SelectedObservations::MeanIothers()
  {
    std::vector<IsigI> mnothers(nobs, IsigI(0.0,0.0));
    if (Nused <= 0) return mnothers;
    if (State == 0) Average();
    if (Nused > 1) {
      // we need at least 2 observations
      float varothers;
      float g;
      double wg2others;

      for (int i=0;i<nobs;i++) {
        if (use[i]) {
          // <I>(others)  ie excluding this observation
          //  and its variance
          const double MINWG2 = 1.0e-30;
          wg2others = Max(sumwg2 - wg2[i], MINWG2); // trap very small wg2 for rounding errors
          varothers = 1./wg2others;
          g = this_ref->get_observation(i).Gscale();
          mnothers[i].I() = g * (sumwgI - wgI[i]) * varothers;
          mnothers[i].sigI() = g * sqrt(varothers);
        }
      }
    }
    State = +1;
    return mnothers;
  }
  // ------------------------------------------------------------
  std::vector<float> SelectedObservations::Delta2()
  // List of deviations delta2 (ie delI/sigma(I) ) where delI
  //  is difference from mean of all observations
  //   returns delta2(NobsRefl), unused slots set = 0.0 ie not closed down
  //   delta2.size() = total number of observations in reflection
  {
    std::vector<float> delta2(nobs,0.0);
    if (Nused <= 0) return delta2;
    if (State == 0) Average();
    if (Nused > 1) {
      float fac = sqrt(float(Nused)/(Nused-1));
      for (int i=0;i<nobs;i++) {
        if (use[i]) {
          float delI = (this_ref->get_observation(i).kI() - avIsigI.I());
          float sigmai = this_ref->get_observation(i).ksigI();
          delta2[i] = fac * delI/sigmai;
        }
      }
    }
    return delta2;
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
  std::vector<float> SelectedObservations::sigmaI()
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
    if (Nused <= 0) return Nrej;
    if (State == 0) Average();

    float I;
    int lsmaller, llarger;

    // Loop until only one observation if necessary
    // ie iterative outlier rejection
    while (Nused > 1) {
      //recalculate deviations if required
      if (State < +2) Deviations();
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
            outliers[lsmaller] = true;
            outliers[llarger] = true;
            Nrej += 2;
            Nused -= 2;
            //^+
            if (DEBUG) std::cout << "RejBoth " << delta[lsmaller]
                                 << " " << delta[llarger]
                                 << "  hkl " << this_ref->hkl().format() << "\n";

          } else if (rej2policy == scala::RejectFlags::REJECTSMALLER) {
            use[lsmaller] = false;
            outliers[lsmaller] = true;
            Nrej++;
            Nused--;
            //^+
            if (DEBUG) std::cout << "RejSmaller " << delta[lsmaller]
                                 << "  hkl " << this_ref->hkl().format() << "\n";
          } else if (rej2policy == scala::RejectFlags::REJECTLARGER) {
            use[llarger] = false;
            outliers[llarger] = true;
            Nrej++;
            Nused--;
            //^+
            if (DEBUG) std::cout << "RejLarger " << delta[llarger]
                                 << "  hkl " << this_ref->hkl().format() << "\n";
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
          if (DEBUG) std::cout << "Reject " << delta[lworst]
                               << "  hkl " << this_ref->hkl().format() << "\n";
        } else {
          break;}  // no rejection
      } else {break;}   // no rejection
    }
    //recalculate deviations if required
    if (State < +2) {
      Deviations();}
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
  // Return list of index numbers for each outlier observation, if any
  // Calls "Outliers" first
  {
    Outliers(rejflags);
    std::vector<int> idxlist;
    for (int i=0;i<nobs;++i) {
      if (outliers[i]) {
        idxlist.push_back(i);
      }
    }
    return idxlist;
  }
  // ------------------------------------------------------------
  //! return formatted version of weight
  std::string SelectedObservations::formatWeightType
    (const WeightType::AverageWeightType& weighttype )
  // Weight type for averaging
  //   UNIT        unit weights
  //   VARIANCE    weight = 1/variance
  //   SQRTSCALE   weight = 1/sqrt(g)  g = 1/scale
  {
    if (weighttype == WeightType::UNIT) {return "unit weights";}
    if (weighttype == WeightType::VARIANCE) {return "variance weights";}
    if (weighttype == WeightType::SQRTSCALE) {return "SquareRoot(scale) weights";}
    return "";
  }
  // ------------------------------------------------------------

}
