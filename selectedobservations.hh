// selectedobservations.hh
//
//  Class to store and work with a subset of observations belonging to a reflection
//  eg from one dataset, or I+, I-
//
#ifndef SELECTEDOBSERVATIONS_HEADER
#define SELECTEDOBSERVATIONS_HEADER

#include "hkl_datatypes.hh"
#include "hkl_unmerge.hh"
#include "weighttype.hh"

namespace scala
{
  //  class reflection;
  //--------------------------------------------------------------
  class IvarI {
    // Just a pair, I, var(I), to save some sqrt & square
  public:
    IvarI() {clipper::Util::set_null(I_); clipper::Util::set_null(varI_); }
    IvarI(const float& I, const float& varI) {I_=I; varI_=varI;}

    // Accessors
    const float&    I() const {return  I_;}
    const float& varI() const {return varI_;}
    // write access
    float&    I() {return    I_;}
    float& varI() {return varI_;}

  private:
    float I_, varI_;
  };
  // ------------------------------------------------------------
  class Iwgrp {
    // just for CumulativeCChalf
  public:
    Iwgrp(){}
    Iwgrp(const double& II, const double& ww, const int& ibatch)
      : I(II), w(ww), batch(ibatch) {}

    double I; // intensity
    double w; // weight
    int batch;
  };
  //--------------------------------------------------------------
  class SelectedObservations
  //! Observations selected from a reflection for:-
  //!   (a) a dataset or all datasets
  //!   (b) all I+ & I- (ALL); just I+ (IPLUS); or just I- (IMINUS)
  {
  public:
    SelectedObservations(){}
    
    //! constructor for selecting datasets & anomalous class
    //! if datasetIndex < 0 select everything
    SelectedObservations(const reflection& Refl, const int& datasetIndex,
			 const AnomalousClass& Anomclass);

    //! constructor for selecting datasets & anomalous class
    //! if datasetIndex < 0 select everything, set weight type
    SelectedObservations(const reflection& Refl, const int& datasetIndex,
			 const AnomalousClass& Anomclass,
			 const WeightType::AverageWeightType& weightType);

    //! Initialise, selecting datasets & anomalous class
    void init(const reflection& Refl,
	      const int& datasetIndex,
	      const AnomalousClass& Anomclass);

    //! Initialise, selecting datasets & anomalous class, weight type
    void init(const reflection& Refl,
	      const int& datasetIndex,
	      const AnomalousClass& Anomclass,
	      const WeightType::AverageWeightType& weightType);

    //! Reset State to force recalculation of average etc
    void reset() {State = 0;}

    //! Divide into Npart parts
    void SetNpart(const int& Npart);

    // Access
    int Number() const {return Nused;}  //!< number selected
    int Nobs() const {return nobs;}     //!< number of reflection

    //! Next used & accepted observation, returns -1 if end 
    int next_observation(observation& obs) const;  
    void reset_next() {nextobs = -1;}

    Hkl hkl() const {return this_ref->hkl();}  //!< reduced hkl

    //! return reflection
    reflection Reflection() const {return *this_ref;}

    //! Set weight
    void SetWeight(const WeightType::AverageWeightType& weightType);

    //! Set variance weights
    void SetVarianceWeights() {SetWeight(WeightType::VARIANCE);}

    //! Set SqrtScale weights
    void SetSqrtScaleWeights() {SetWeight(WeightType::SQRTSCALE);}

    //! Set Scale weights
    void SetScaleWeights() {SetWeight(WeightType::SCALE);}

    //! Set unitweights
    void SetUnitWeights() {SetWeight(WeightType::UNIT);}

    //! Average I, weight as specified
    IsigI Average();

    //! SD(<I>) from variances
    double SDvariance() const {return sdI;}

    //! SD(<I>) from sampleSD, if done (else 0.0)
    double SDsample() const {return sdIs;}

    //! Var(<I>) from sample variance, unconditional from Average()
    double Variance() const;

    //! Sample variance, unconditional from Average()
    double sampleVariance() const;

    //! Set sample variance, minimum number of values (<0 to switch off)
    static void SetSampleSD (const int& minSample=10);

    //! Average I, variance weight for part of data
    // If WhichPart >= 0 (0 -> npart-1), use only selected random part
    //              < 0  use all accepted
    IsigI AveragePart(const int& WhichPart = -1);

    //! Get average I for each random, return false unless both are present
    bool HalfAverages(float& I1, float& I2);
    //! Get average I for each random part, return false unless all are present
    bool PartAverages(std::vector<float>& Is);
    //! Get average IsigI for each random part, return false unless both are present
    bool HalfAveragesSigI(IsigI& I1sig, IsigI& I2sig);

    // individual SDs from sample variance, if used,= 0.0 for unused slots
    std::vector<double> sampleSDI() const {return sdIsample;}
    bool SampleSDused() const {return sampleSDused;}

    //! List of deviations delta (ie delI/sigma(I) ) where delI
    //!  is difference from mean of other observations
    //!   returns deltasize() = NobsRefl, unused slots set = 0.0 ie not closed down
    //  If (fromSample && sampleSDused), then use sample variance, else
    //   use individual variance
    std::vector<float> Deviations(const bool& fromSample=false);    
    //! Deviations including rejected outliers
    std::vector<float> DeltaAll() const {return delta;}
    //! List of deviations delta2 (ie fac * delI/sigma(I) ) where delI
    //!  is difference from mean of all observations and
    //!  fac = sqrt(n/n-1)
    //  If (fromSample && sampleSDused), then use sample variance, else
    //   use individual variance
    std::vector<float> Delta2(const bool& fromSample=false);

    //! List of deviations delta3 (ie delI/SDsample(I) ) where delI
    //!  is difference from mean of all observations
    std::vector<float> Delta3();

    // Mean Chi^2, goodness of fit, Mean((I-Iothers)/sigma)
    //  if (fromSample && sampleSDused), then use sample variance
    //  if (fromSample && !sampleSDused), return 0.0
    // else  use individual variance
    double chiSq(const bool& fromSample);

    //! weights for each observation, = 0.0 if not used
    std::vector<double> weights() const {return  wj;}

    // For each observation, return mean of other observations,
    //   scaled to each observation
    //   returns mnothers(NobsRefl), unused slots set = 0.0 ie not closed down
    // 
    std::vector<IvarI> MeanIothers();

    //! List of delI (scaled)
    //!  returns delI(Nobs), unused slots set = 0.0 ie not closed down
    std::vector<float> DelI();        // variance-weighted <I>

    //! List of sigma(I)
    //!   returns sigmaI(Nobs), unused slots set = 0.0 ie not closed down
    std::vector<float> sigmaI() const;

    // List of I, sigma(I) (scaled)
    //   returns I,sigI(nobs), unused slots set = 0.0 ie not closed down
    std::vector<IsigI> IsigIlist() const;

    // List of (I, weight, batch)
    std::vector<Iwgrp> iwgroup() const;

    //! Run number for kobs'th observation
    //! note that kobs is the index into the whole list, including unused ones
    int Run(const int& kobs) const;
    //! True if kobs'th observation is Fully recorded
    //! note that kobs is the index into the whole list, including unused ones
    bool Full(const int& kobs) const;
    //! True if any observations are Fully recorded
    bool AnyFull() const;

    //!  Reject outliers
    int Outliers(const RejectFlags& rejflags);
    //! Return list of index numbers for each outlier observation, if any
    //! Calls "Outliers" first
    std::vector<int> OutlierIndexList(const RejectFlags& rejflags);
    //! Return true if outliers found even if not rejected
    bool Discrepant() const {return discrepant;}

  // On entry:
  //  sdrej   rejection limit (sds) for more than 2 observations
  //  sdrej2  rejection limit (sds) for 2 observations
  //  rej2policy    what to do with two observations
  //           = REJECT           reject both
  //           = KEEP             keep both
  //           = REJECTLARGER     reject larger
  //           = REJECTSMALLER    reject smaller

  private:
    const reflection* this_ref;
    std::vector<bool> use;         // use flags, initially all true
    std::vector<bool> outliers;    // outlier flags, initially all false
    std::vector<int> part;         // randomly assigned to 0 -> npart-1
    bool discrepant;     // true if outliers found even if not rejected
    int npart;                     // number of part lists
    int nobs;     // number of observations in reflection
    int Nused;    // number used (without outliers)
    // individual SDs from sample variance, if used,= 0.0 for unused slots
    std::vector<double> sdIsample;
    bool sampleSDused;  // true if sample SD used for this reflection

    // Sums for deviations & outliers
    std::vector<double> wI;  // w I'
    std::vector<double> wj;  // w
    double sumwI;  // Sum(w I')
    double sumwj;  // Sum(w)

    // wv = 1/var(I) if variance weight
    //    else = w^2 var(I)
    std::vector<double> wv; 
    double sumwv;  // Sum(wv)
    // Deviations delI/sigma  from "others"
    std::vector<float> delta;     // for current list
    WeightType::AverageWeightType weighttype;   // type of weighting for average
    MeanVariance mv;  // for sample variance

    // true to calculate sample SD in avIsigI from values instead of error propagation,
    // provided that the number is > minimumsample
    static bool sampleSD;
    static int minimumsample;

    IsigI avIsigI;
    double sdI;     // SD <I> from weights
    double sdIs;    // SD <I> from sample (if done)
    // state:
    // =  0  observations stored
    // = +1  average calculated
    // = +2  deviations calculated
    // = +3  outliers calculated
    int State; 
    mutable int nextobs;  // index to next observation

    double Weight(const double& sd, const double& g) const;

    void samplemean();


  };
}
#endif
