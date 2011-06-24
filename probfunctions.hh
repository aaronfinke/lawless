// probfunctions.hh

#ifndef PROBFUNC_HEADER
#define PROBFUNC_HEADER

#include <clipper/clipper.h>

#include "util.hh"


namespace scala {
  //--------------------------------------------------------------
  double GaussProb(const double& val,
		   const double& mean,
		   const double& sd);
  // Gaussian probability density function
  // returns p = (1/(sd*sqrt(twopi))) exp (-1/2 z^2)
  //   where z = (val-mean)/sd
  //--------------------------------------------------------------
  double TruncatedGaussProb(const double& val,
			    const double& mean,
			    const double& sd,
			    const double& minval,
			    const double& maxval
			    );
  // Truncated Gaussian probability density function
  // returns p = T * (1/(sd*sqrt(twopi))) exp (-1/2 z^2)
  //   where z = (val-mean)/sd
  //   and T allows for the fact that minval =< val < maxval
  //     zmin = (minval - mean)/sd
  //     zmax = (maxval - mean)/sd
  //     T = 1/(Phi(zmax) - Phi(zmin))
  //     Phi(z) = 1/2*[1 + erf(z/Sqrt2)]
  //                     ie the cumulative distribution function
  //--------------------------------------------------------------
  double TruncatedLorentzianProb(const double& val,
				 const double& mean,
				 const double& sd,
				 const double& minval,
				 const double& maxval
				 );
  // Truncated Lorentzian probability density function
  // returns p = T *(1/pi)*[1/(sd*(1+z^2))]
  //   where z = (val-mean)/sd
  //   and T is a normalisation factor
  //     zmin = (minval - mean)/sd
  //     zmax = (maxval - mean)/sd
  //     T = pi/(atan(zmax)-atan(zmin))
  //--------------------------------------------------------------
  class MeanModelBase
  // Base class for declining model probability functions
  {
  public:
    MeanModelBase() : High(0.0), Low(1.0) {}
    // Maps HighProb->LowProb on to 0->1
    MeanModelBase(const double& HighProb, const double& LowProb) 
      : High(HighProb), Low(LowProb) {}

    virtual ~MeanModelBase(){}

    virtual double dprob(const double& am) const {return -1.0;} // dummy

    double HighProb() const {return High;}
    double LowProb() const {return Low;}

  protected:
    // Returns v mapped on to range 0->1
    double nmap(const double& v) const
    {return Min(1.0,Max(0.0,(v-High)/(Low-High)));}

  private:
    double High;
    double Low;

  };
  //--------------------------------------------------------------
  class DM_power : public MeanModelBase
  //! (1 - x^k)^(1/k), default k = 2
  {
  public:
    DM_power(const double& HighProb, const double& LowProb)
      : MeanModelBase(HighProb, LowProb), k(2) {}

    double dprob(const double& am) const;

    void SetPower(const double& power) {k = power;}
  private:
    double k;  // power to use
  };
  //--------------------------------------------------------------
  class DM_2sqrt : public MeanModelBase
		   //!  sqrt(1 - m^2)
  {
  public:
    DM_2sqrt(const double& HighProb, const double& LowProb)
      : MeanModelBase(HighProb, LowProb) {}

    double dprob(const double& am) const;
  };
  //--------------------------------------------------------------
  class DM_lin : public MeanModelBase
  {
  public:
    DM_lin(const double& HighProb, const double& LowProb)
      : MeanModelBase(HighProb, LowProb) {}

    double dprob(const double& am) const;
  };

  class DM_cubic : public MeanModelBase
  {
  public:
    DM_cubic(const double& HighProb, const double& LowProb)
      : MeanModelBase(HighProb, LowProb) {}
  
    double dprob(const double& am) const;
  };

  class DM_cubicsu : public MeanModelBase
  {
  public:
    DM_cubicsu(const double& HighProb, const double& LowProb)
      : MeanModelBase(HighProb, LowProb) {}

    double dprob(const double& am) const;
  };

  class DM_1minusmSq : public MeanModelBase
  {
  public:
    DM_1minusmSq(const double& HighProb, const double& LowProb)
      : MeanModelBase(HighProb, LowProb) {}

    double dprob(const double& am) const;
  };
  class DM_1minusmCu : public MeanModelBase
  {
  public:
    DM_1minusmCu(const double& HighProb, const double& LowProb)
      : MeanModelBase(HighProb, LowProb) {}

    double dprob(const double& am) const;
  };
  //--------------------------------------------------------------
  class IntgrtProb
  // Probability of Gaussian integrated over prior for mean
  //   v = N(m,sigma)
  //   model for m is provided of function object derived
  //   from MeanModelBase
  //   m is integrated over the range m1->m2, where m1 corresponds
  //   to the higher probability, m2 to the lower (p(m2)=0)
  //   
  {
  public:
    IntgrtProb() : nd(10000) {}

    // Note that this stores a pointer to the external function
    // which therefore must continue to exist while this class
    // is used
    IntgrtProb(const MeanModelBase& DModel) : nd(10000)
    {DMeanModel = &DModel;}

    void init(const MeanModelBase& DModel) 
    {
      nd = 10000;
      DMeanModel = &DModel;
    }

    // Probability of Gaussian integrated over prior for mean
    //   v = N(m,sigma)
    //   model for m is provided of function object derived
    //   from MeanModelBase
    //   m is integrated over the range m1->m2, where m1 corresponds
    //   to the higher probability, m2 to the lower (p(m2)=0)
    //   Use Gaussian normalised to allow for truncation at 
    //   minval <= v <= maxval (unless both minval & maxval = 0.0)
    double Prob(const double& val, const double& sd,
			    const double& minval=0.0,
			    const double& maxval=0.0);

    // Probability of Lorentzian integrated over prior for mean
    //   v = L(m,sigma)
    //   model for m is provided of function object derived
    //   from MeanModelBase
    //   m is integrated over the range m1->m2, where m1 corresponds
    //   to the higher probability, m2 to the lower (p(m2)=0)
    //   Use Lorentzian normalised to allow for truncation at 
    //   minval <= v <= maxval (unless both minval & maxval = 0.0)
      double LorentzProb(const double& val, const double& sd,
				     const double& minval,
				     const double& maxval);

  private:
    const MeanModelBase* DMeanModel;
    int nd;  // number of integration intervals (divides range)
  };
//--------------------------------------------------------------
// Get probability of val(sd) being near posmean
// Assume Gaussian error around posmean
//  ie for case A
//  p(A|val,sd)  ~ N(posmean,sd)
//  p(!A|val,sd) = p(!A|val,sd,m) p(m)
// p(m) given by DModel
//
  double ProbBiassed(const double& val, const double& sd,
		   const double& posmean, const MeanModelBase& DModel);
//--------------------------------------------------------------
// Unbiassed normal probability
// if A, p(A|val,sd) ~ N(posmean,sd)
// if !A, p(!A|val,sd) ~ N(negmean,sd)
// Returns p(A) normalised
double ProbUnbiassed(const double& val, const double& sd,
		     const double& posmean,const double& negmean);
//--------------------------------------------------------------
}
#endif

