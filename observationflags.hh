// observationflags.hh

// Classes for handling observation flags & status

//  1) Mosflm reflection flags ObsFlag
//     Information as read from file & stored as is
//
//     Classes:
//        ObservationFlag          stores an ObsFlag for an
//                                 observation or observation_part
//        ObservationFlagControl   controls whether flagged
//                                 observations are accepted
// 
//  2) volatile status flags (Status)
//
//      Classes
//        ObservationStatus          stores status for an observation
//        ObservationStatusSet       sets status flags depending on control
//
//
#ifndef OBSERVATIONFLAGS_HEADER
#define OBSERVATIONFLAGS_HEADER

#include <string>

namespace scala
{
  //--------------------------------------------------------------
  float PackFlagValues(const float& BGratio,
		       const float& PKratio, const float& Gradient);
  //--------------------------------------------------------------
  class ObservationFlag
  // Flags & status for an observation, or part
  // 
  // (1)  Observation flag ObsFlag
  // Bit lmflag       packed flag
  //  0     1         BGRATIO too large                      FLAG_BGRATIO
  //  1     2         PKRATIO too large                      FLAG_PKRATIO
  //  2     4         Negative > 5*sigma                     FLAG_TOONEGATIVE
  //  3     8         Gradient too high                      FLAG_GRADIENT
  //  4    16         Profile fitted overload                FLAG_OVERLOAD
  //  5    32         Profile fitted "edge" reflection       FLAG_EDGE
  //  
  //
  // (2) bgpk       packed BG/PK ratios
  //   packed as 
  //    bgpk = int(Float(Nint(PkRatio * 100.)*100) +
  //               Min(BgRatio*10.,99.)) + Min(Gradient, 0.99)
  {
  public:
    ObservationFlag() : bitflags(0), bgpk(0.0) {}
    ObservationFlag(const int& ObsFlag, const float& bgpkratio)
      : bitflags(ObsFlag), bgpk(bgpkratio){}
    ObservationFlag(const ObservationFlag& flag);

    // Combine with another flags
    void AddFlag(const ObservationFlag& other);
        
    // true if no flags set
    bool OK() const {return (bitflags == 0);}
    
    // Return complete packed flags
    unsigned int Flags() const {return bitflags;}
    float BgPk() const {return bgpk;}

    // Return unpacked values of BGratio, PKratio, gradient
    void BgPkValues(float& BGratio, float& PKratio, float& Gradient) const;

    // Return individual bit flags
    //  true if bit set
    bool TestBGratio() const {return (bitflags & 1) != 0;}
    bool TestPKratio() const {return (bitflags & 2) != 0;}
    bool TestTooNeg() const {return (bitflags & 4) != 0;}
    bool TestGradient() const {return (bitflags & 8) != 0;}
    bool TestOverload() const {return (bitflags & 16) != 0;}
    bool TestEdge() const {return (bitflags & 32) != 0;}

    //! return brief formatted version of which flags are set
    std::string format() const;

    static const int FLAG_BGRATIO;     //      1
    static const int FLAG_PKRATIO;     //      2
    static const int FLAG_TOONEGATIVE; //      4
    static const int FLAG_GRADIENT;    //      8
    static const int FLAG_OVERLOAD;    //     16
    static const int FLAG_EDGE;        //     32
    
  private:
    unsigned int bitflags;
    float bgpk;
  };
  //--------------------------------------------------------------
  class ObservationFlagControl
  // Control of observation flags for all observations
  // Tests for acceptance
  // Counts of numbers 
  {
  public:
    ObservationFlagControl();

    // Construct setting flags
    //  BGrlim, PKrlim, Gradlim < 0.0
    //       means always reject if bit flag set
    //       > 0.0 accept if bit flag set and value < limit
    //  AcceptOverload
    ObservationFlagControl(const float& BGrlim,  const float& PKrlim,
		   const float& Gradlim,
		   const bool& AcceptOverload, const bool& AcceptEdge);

    void SetBGRlimit(const float& BGrlim);
    void SetPKRlimit(const float& PKrlim);
    void SetGradlimit(const float& Gradlim);
    void SetAcceptOverload();
    void SetAcceptEdge();

    // Returns true is observation accepted, & count them
    bool IsAccepted(const ObservationFlag& flag);

    // Initialise to defaults (no acceptances & clear counts)
    void Init();
    // Clear counts
    void Clear();
    
    std::string PrintCounts() const;

    int NumAccOverload() const {return Naccoverload;}

  private:
    float bgrlim;  // maximum on BGratio   < 0.0 for no test
    float pkrlim;  // maximum on PKratio   < 0.0 for no test
    float grdlim;  // maximum on gradient  < 0.0 for no test
    bool acceptoverload; // true to accept overloads
    bool acceptedge;     // true to accept edge reflections
    
    int NBGratio;    // count of BGratio flags
    int NPKratio;    // count of PKratio flags
    int NTooNeg;      // count of too negative flags
    int NGradient;   // count of gradient flags
    int Noverload;   // count of overload flags
    int Nedge;       // count of edge flags
    
    int NaccBGratio;    // count of BGratio acceptances
    int NaccPKratio;    // count of PKratio acceptances
    int NaccTooNeg;     // count of Too Negative acceptances
    int NaccGradient;   // count of gradient acceptances
    int Naccoverload;   // count of overload acceptances
    int Naccedge;       // count of edge acceptances

    float MaxBGratio;    // Maximum BGratio
    float MaxPKratio;    // Maximum PKratio 
    float MaxGradient;   // Maximum gradient

    float MaxAccBGratio;    // Maximum accepted BGratio
    float MaxAccPKratio;    // Maximum accepted PKratio 
    float MaxAccGradient;   // Maximum accepted gradient
  };
  // --------------------------------------------------------------
  class ObservationStatus
  //  Volatile status
  // Bit  
  //  0     1     rejected based on ObsFlag
  //  1     2     outside (possibly volatile) resolution limits eg within run
  //  2     4     outlier (deviation too large): within I+ or I- set
  //  3     8     outlier (deviation too large): between I+ & I- sets
  //  4    16     > Emax limit
  //  5    32     too strong for scaling   (also Emax)
  //  6    64     too weak for scaling     (Emin)
  {
  public:
    enum ObsStatusFlag {OBSSTAT_FLAG=1, OBSSTAT_RESOLUTION=2, OBSSTAT_OUTLIER=4, OBSSTAT_OUTLIERANOM=8,
	       OBSSTAT_EMAX=16, OBSSTAT_STRONG=32, OBSSTAT_WEAK=64};

    ObservationStatus() :bitflags(0){}
    ObservationStatus(const unsigned int& flags) :bitflags(flags){}

    void Clear() {bitflags = 0;}
    // Status access
    //  accepted if no status bits are set
    bool IsAccepted() const {return (bitflags == 0);}
    // Clear all flags except ObsFlag & resolution flag
    void ResetStatus() {bitflags &= 3;}

    unsigned int Bitflags() const {return bitflags;}

    // ObsFlag acceptance
    void SetObsFlag() {bitflags |= OBSSTAT_FLAG;}  // set
    void UnsetObsFlag() {bitflags &= (wordmask-OBSSTAT_FLAG);}  // unset
    //  true if bit set
    bool TestObsFlag() const {return (bitflags & OBSSTAT_FLAG) != 0;} // test

    // Outside (run) resolution limits 
    void SetResolution() {bitflags |= OBSSTAT_RESOLUTION;}
    void UnSetResolution() {bitflags &= (wordmask-OBSSTAT_RESOLUTION);}
    bool TestResolution() const {return (bitflags & OBSSTAT_RESOLUTION) != 0;}

    // Outlier within I+ or I-
    void SetOutlier() {bitflags |= OBSSTAT_OUTLIER;}
    void UnsetOutlier() {bitflags &= (wordmask-OBSSTAT_OUTLIER);}
    bool TestOutlier() const {return (bitflags & OBSSTAT_OUTLIER) != 0;}

    // Outlier between I+ or I-
    void SetOutlierAnom() {bitflags |= OBSSTAT_OUTLIERANOM;}
    void UnsetOutlierAnom() {bitflags &= (wordmask-OBSSTAT_OUTLIERANOM);}
    bool TestOutlierAnom() const {return (bitflags & OBSSTAT_OUTLIERANOM) != 0;}

    // Too large, E > Emax
    void SetEmax() {bitflags |= OBSSTAT_EMAX;}
    void UnsetEmax() {bitflags &= (wordmask-OBSSTAT_EMAX);}
    bool TestEmax() const {return (bitflags & OBSSTAT_EMAX) != 0;}

    // Too strong for scaling
    void SetTooStrong() {bitflags |= OBSSTAT_STRONG;}
    void UnsetTooStrong() {bitflags &= (wordmask-OBSSTAT_STRONG);}
    bool TestTooStrong() const {return (bitflags & OBSSTAT_STRONG) != 0;}

    // Too weak for scaling
    void SetTooWeak() {bitflags |= OBSSTAT_WEAK;}
    void UnsetTooWeak() {bitflags &= (wordmask-OBSSTAT_WEAK);}
    bool TestTooWeak() const {return (bitflags & OBSSTAT_WEAK) != 0;}
    
  private:
    unsigned int bitflags;
    static const unsigned int wordmask = 0xFFFF;
  };
  //--------------------------------------------------------------
  class ObservationStatusSet
  // Control of volatile status for an observation
  {
  };
}

#endif
