// sdctypes.hh

#ifndef SDCTYPES_HEADER
#define SDCTYPES_HEADER

#include "hkl_datatypes.hh"

namespace scala
{

  class observation; // forward declaration

  class SDcorrection
  // sd' = SDfac * Sqrt(sd^2 + SdB * I + (SDadd * I)^2)
  {
  public:
    SDcorrection() : sdfac(1.0), sdadd(0.0), fixsdb(true) {ClearRestraints();}
    SDcorrection(const double& SDfac, const double& SDb, const double& SDadd,
		 const bool& fixSDb=false);

    void SetFixSDb(const bool& fixSDb) {fixsdb = fixSDb;
      SetRefineParameters();}
    void FixSDb() {fixsdb = true;
      SetRefineParameters();}
    void UnFixSDb() {fixsdb = false;
      SetRefineParameters();}

    // Reset minimum & maximum
    void ResetRange() {mincorr = +1.0e10; maxcorr = -mincorr;}

    // Multiply SDfac by update
    void UpdateFactor(const double& update) {sdfac *= update;}

    void Set(const double& SDfac, const double& SDb, const double& SDadd);

    double SDfac() const {return sdfac;} // get SDfac
    double SDb() const {return sdb;} // get SDb
    double SDadd() const {return sdadd;} // get SDadd

    // Correct sigma
    // sigma  uncorrected sigma(Ihl)
    // gscale inverse scale for Ihl
    // Iav    average <Ih>
    double SigmaPrime(const double& sigma, const double& gscale,
		      const double& Iav) const;

    //! in-place correction, using Iav as intensity, returns original uncorrected sd(I)
    float Correct(observation& Observation, const float& Iav) const;

    // return corrected sd as IsigI, using Iav as intensity (scaled by Gscale)
    IsigI Correct(const IsigI& Is, const double& gscale, const float& Iav) const;

    std::string format() const;
    // Range of values applied
    std::pair<double,double> MinMax() const
    {return std::pair<double,double>(mincorr, maxcorr);}

    // = = = for optimisation = = =
    // NB Optimisation parameters are not necessarily SdFac, SdB, SdAdd
    //  Use p1 = Sdfac^2; p2 = Sdfac^2 SdB; p3 = Sdfac^2 Sdadd^2
    //  Sdadd^2 may be negative
    //! return number of parameters
    int Nparams() const {return (fixsdb) ? 2 : 3;}

    // Get vector of parameters (2 or 3)
    std::vector<double> GetParameters() const;

    // Set all parameters from vector
    void SetParameters(const std::vector<double>& params);

    // Initial shifts for each parameter type, scaled by "scale", for Simplex
    std::vector<double> GetShifts(const double& scale) const;

    // Return d(sigma'f)/dp vector (2 or 3 parameters)
    //  sigma    uncorrected sigma
    //  Iav      average I <Ih>
    std::vector<double> GetDerivatives(const float& sigma,
				       const float& Iav) const;

    //! Get vector of lower bounds (0.0 means no bound)
    std::vector<double> LowerBounds() const;
    //! Get vector of upper bounds (0.0 means no bound)
    std::vector<double> UpperBounds() const;
    //! Get vector of "large shifts"
    std::vector<double> LargeShifts() const;

    // - - Restraints
    //! clear all restraints
    void ClearRestraints();
    //! set targets & SD (=0 for no restraint), for 3 parameters always
    //! weight = 1/SD^2
    void SetRestraints(const std::vector<double>& Targets,
		       const std::vector<double>& SDtarget);
    //! set default values for restraints (on SdB only)
    void SetDefaultRestraints();
    //! return restraint values
    void GetRestraints(std::vector<double>& Targets,
		        std::vector<double>& SDtarget) const;

    //! return contribution to restraint residual R2, summed over parameters
    double RestraintR() const;
    //! return derivative vector dR2/dp and Hessian, 2 or 3 parameters
    void RestraintDerivatives
    (std::vector<double>& dp, clipper::Array2d<double>& H) const;

    // = = = 
    
  private:
    double sdfac;  // NB number of parameters in Npar() above
    double sdb;
    double sdadd;
    double sdadd2;  // sdadd^2 negated if necessary
    bool fixsdb;   // true if SdB is fixed at 0.0 and p2 omitted
    mutable std::vector<double> p;  // parameter vector p1, [p2,] p3, p2 may be omitted
    mutable double mincorr, maxcorr;  // minimum & maximum values 

    // restraint target values for sdfac, sdb, sdadd^2, always dimension 3
    // sdadd is restrained as sdadd^2 to allow it to be negative
    std::vector<double> targets; // modified for sdadd^2
    std::vector<double> weights; // target weights, = 0 for no restraint
    std::vector<double> rawtargets; // original targets, unmodified
    std::vector<double> sdtargets;  // target SDs, = 0 for no restraint

    static const double MINVARINFRAC;  // minimum fraction of input variance for correction
    static const double MINSDFAC;      // minimum SDfac
    static const int NPARALL = 3;      // number of parameters = 3

    // Set internal vector p (2 or 3) from Sdfac etc
    //  vector p is (p1,p2,p3), p2 may be missing 
    void SetRefineParameters();


  };
}

#endif
