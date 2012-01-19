// selectedobservations.hh
//
//  Class to store and work with a subset of observations belonging to a reflection
//  eg from one dataset, or I+, I-
//
#ifndef SELECTEDOBSERVATIONS_HEADER
#define SELECTEDOBSERVATIONS_HEADER

#include "hkl_datatypes.hh"
#include "hkl_unmerge.hh"

namespace scala
{
   
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
			 const AnomalousClass& Anomclass=ALL);

    //! Initialise, selecting datasets & anomalous class
    void init(const reflection& Refl,
	      const int& datasetIndex,
	      const AnomalousClass& Anomclass=ALL);

    //! Divide into Npart parts
    void SetNpart(const int& Npart);

    // Access
    int Number() const {return Nused;}  //!< number selected
    int Nobs() const {return nobs;}     //!< number of reflection

    //! Next used & accepted observation, returns -1 if end 
    int next_observation(observation& obs) const;  

    Hkl hkl() const {return this_ref->hkl();}  //!< reduced hkl

    //! return reflection
    reflection Reflection() const {return *this_ref;}

    //! Average I, variance weight
    IsigI Average();

    //! Average I, variance weight for part of data
    // If WhichPart >= 0 (0 -> npart-1), use only selected random part
    //              < 0  use all accepted
    IsigI AveragePart(const int& WhichPart = -1);

    //! Average I, sqrt(1/g) weight
    // Don't store average
    IsigI AverageScaleWt() const;

    //! Get average I for each random, return false unless both are present
    bool HalfAverages(float& I1, float& I2);
    //! Get average I for each random part, return false unless all are present
    bool PartAverages(std::vector<float>& Is);

    //! List of deviations delta (ie delI/sigma(I) ) where delI
    //!  is difference from mean of other observations
    //!   returns deltasize() = NobsRefl, unused slots set = 0.0 ie not closed down
    std::vector<float> Deviations();    
    //! Deviations including rejected outliers
    std::vector<float> DeltaAll() const {return delta;}
    //! List of deviations delta2 (ie fac * delI/sigma(I) ) where delI
    //!  is difference from mean of all observations and
    //!  fac = sqrt(n/n-1)
    std::vector<float> Delta2();        // variance-weighted <I>
    std::vector<float> Delta2scalewt(); // scale-weighted <I>

    // For each observation, return mean of other observations,
    //   scaled to each observation
    //   returns mnothers(NobsRefl), unused slots set = 0.0 ie not closed down
    // 
    std::vector<IsigI> MeanIothers();

    //! List of delI (scaled)
    //!  returns delI(Nobs), unused slots set = 0.0 ie not closed down
    std::vector<float> DelI();        // variance-weighted <I>
    std::vector<float> DelIscalewt(); // scale-weighted <I>

    //! List of sigma(I)
    //!   returns sigmaI(Nobs), unused slots set = 0.0 ie not closed down
    std::vector<float> sigmaI();

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
    // Sums for deviations & outliers
    std::vector<double> wgI;  // w g I
    std::vector<double> wg2;  // w g^2
    double sumwgI;  // Sum(w g I)
    double sumwg2;  // Sum(w g^2)
    // Deviations delI/sigma  from "others"
    std::vector<float> delta;     // for current list
    IsigI avIsigI;
    // =  0  observations stored
    // = +1  average calculated
    // = +2  deviations calculated
    // = +3  outliers calculated
    int State; 
    mutable int nextobs;  // index to next observation
  };
}
#endif
