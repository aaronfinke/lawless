// probfunctions.cpp

#include "probfunctions.hh"
#include <scitbx/math/erf.h>


namespace scala{
//--------------------------------------------------------------
  double GaussProb(const double& val,
		   const double& mean,
		   const double& sd)
  // Gaussian probability density function
  // returns p = (1/(sd*sqrt(twopi))) exp (-1/2 z^2)
  //   where z = (val-mean)/sd
{
  if (sd <= 0.0)
    clipper::Message::message(clipper::Message_fatal
			      ("GaussProb: sd must be > 0"));
  double z = (val-mean)/sd;
  return (1./(sd*sqrt(clipper::Util::twopi()))) * exp (-0.5*z*z);
}
//--------------------------------------------------------------
  double PhiCDF(const double& z)
  // Phi(z) = 1/2*[1 + erf(z/Sqrt2)]
  //   ie the Gaussian cumulative distribution function
  {
    return 0.5*(1.0+scitbx::math::erf<double>(z/sqrt(2.0)));
  }
//--------------------------------------------------------------
  double TruncatedGaussProb(const double& val,
			    const double& mean,
			    const double& sd,
			    const double& minval,
			    const double& maxval
			    )
  // Truncated Gaussian probability density function
  // returns p = T * (1/(sd*sqrt(twopi))) exp (-1/2 z^2)
  //   where z = (val-mean)/sd
  //   and T allows for the fact that minval =< val < maxval
  //     zmin = (minval - mean)/sd
  //     zmax = (maxval - mean)/sd
  //     T = 1/(Phi(zmax) - Phi(zmin))
  //     Phi(z) = 1/2*[1 + erf(z/Sqrt2)]
  //                     ie the cumulative distribution function
{
  if (sd <= 0.0)
    clipper::Message::message(clipper::Message_fatal
			      ("GaussProb: sd must be > 0"));
  if (val < minval || val > maxval)
    clipper::Message::message(clipper::Message_fatal
			      ("GaussProb: value must lie between limits"));
  if (mean < minval || mean > maxval)
    clipper::Message::message(clipper::Message_fatal
			      ("GaussProb: mean must lie between limits"));
  double z = (val-mean)/sd;
  double z1 = (minval-mean)/sd;
  double z2 = (maxval-mean)/sd;
  double w = PhiCDF(z2) - PhiCDF(z1);
  w = 1.0/(w*(sd*sqrt(clipper::Util::twopi())));
  return w * exp (-0.5*z*z);
}
//--------------------------------------------------------------
double TruncatedLorentzianProb(const double& val,
			       const double& mean,
			       const double& sd,
			       const double& minval,
			       const double& maxval
			       )
// Truncated Lorentzian probability density function
// returns p = T *(1/pi)*[1/(sd*(1+z^2))]
//   where z = (val-mean)/sd
//   and T is a normalisation factor
//     zmin = (minval - mean)/sd
//     zmax = (maxval - mean)/sd
//     T = pi/(atan(zmax)-atan(zmin))    [note: pi cancels]
{
  if (sd <= 0.0)
    clipper::Message::message(clipper::Message_fatal
			      ("LorentzianProb: sd must be > 0"));
  if (val < minval || val > maxval)
    clipper::Message::message(clipper::Message_fatal
			      ("LorentzianProb: value must lie between limits"));
  if (mean < minval || mean > maxval)
    clipper::Message::message(clipper::Message_fatal
			      ("LorentzianProb: mean must lie between limits"));
  double sz = (val-mean);
  double z1 = (minval-mean)/sd;
  double z2 = (maxval-mean)/sd;
  double w = 1.0/(atan(z2) - atan(z1));
  double p = std::abs(w * (sd / (sd*sd + sz*sz)));
  //^
  //^  if (p < 0.0) {
  //^    std::cout << "Negative probability! " << p << " " << w << "\n";
  //^  }
  return p;
}
//--------------------------------------------------------------
double DM_power::dprob(const double& am) const
  //! (1 - x^k)^(1/k), default k = 2
{
  double m = nmap(am);
  return pow((1.0-pow(m,k)), 1.0/k);
}
//--------------------------------------------------------------
double DM_2sqrt::dprob(const double& am) const
{
  double m = nmap(am);
  return sqrt(1.0-m*m);
}
//--------------------------------------------------------------
double DM_lin::dprob(const double& am) const
{
  double m = nmap(am);
  return (1.0-m);
}
//--------------------------------------------------------------
double DM_cubic::dprob(const double& am) const
{
  double m = nmap(am);
  return (1.0-m*m*m);
}
//--------------------------------------------------------------
double DM_cubicsu::dprob(const double& am) const
{
  double m = nmap(am);
  return (1.0-m*m*m)*(1.0-m*m*m);
}
//--------------------------------------------------------------
double DM_1minusmSq::dprob(const double& am) const
{
  double m = nmap(am);
  return (1.0-m)*(1.0-m);
}
//--------------------------------------------------------------
double DM_1minusmCu::dprob(const double& am) const
{
  double m = nmap(am);
  return (1.0-m)*(1.0-m)*(1.0-m);
}
//--------------------------------------------------------------
// Probability of Gaussian integrated over prior for mean
//   v = N(m,sigma)
//   model for m is provided of function object derived
//   from MeanModelBase
//   m is integrated over the range m1->m2, where m1 corresponds
//   to the higher probability, m2 to the lower (p(m2)=0)
//   Use Gaussian normalised to allow for truncation at 
//   minval <= v <= maxval (unless both minval & maxval = 0.0)
  double IntgrtProb::Prob(const double& val, const double& sd,
			  const double& minval,
			  const double& maxval)
// Integrate probability over possible ideal values
{
  if (sd <= 0.0)
    clipper::Message::message(clipper::Message_fatal
			      ("IntgrtProb::Prob: sd must be > 0"));
  double m1 = DMeanModel->HighProb();
  double m2 = DMeanModel->LowProb();
  // Swap if wrong way round (make m2 > m1)
  // integration is independent of order
  if (m1 > m2)
    {
      double t = m1;
      m1 = m2;
      m2 = t;
    }
  double r = m2-m1;
  double d = r/double(nd);  // integration interval
  double sump = 0.0;
  double sumw = 0.0;
  double p;
  bool truncated = (minval != 0.0) || (maxval != 0.0);

  for (double m=m1;m<m2;m+=d)
    {
      if (truncated)
	{p = TruncatedGaussProb(val,m,sd,minval,maxval);}
      else
	{p = GaussProb(val,m,sd);}
      double w = DMeanModel->dprob(m);
      sump += p * w;
      sumw += w;
    }	  
  return sump/sumw;
}
//--------------------------------------------------------------
// Probability of Lorentzian integrated over prior for mean
//   v = L(m,sigma)
//   model for m is provided of function object derived
//   from MeanModelBase
//   m is integrated over the range m1->m2, where m1 corresponds
//   to the higher probability, m2 to the lower (p(m2)=0)
//   Use Lorentzian normalised to allow for truncation at 
//   minval <= v <= maxval (unless both minval & maxval = 0.0)
  double IntgrtProb::LorentzProb(const double& val, const double& sd,
				 const double& minval,
				 const double& maxval)
// Integrate probability over possible ideal values
{
  if (sd <= 0.0)
    clipper::Message::message(clipper::Message_fatal
			      ("IntgrtProb::LorentzProb: sd must be > 0"));
  double m1 = DMeanModel->HighProb();
  double m2 = DMeanModel->LowProb();
  // Swap if wrong way round (make m2 > m1)
  // integration is independent of order
  if (m1 > m2)
    {
      double t = m1;
      m1 = m2;
      m2 = t;
    }
  double r = m2-m1;
  double d = r/double(nd);  // integration interval
  double sump = 0.0;
  double sumw = 0.0;
  double p;
  bool truncated = (minval != 0.0) || (maxval != 0.0);

  for (double m=m1;m<m2;m+=d) {
    if (truncated)
	{p = TruncatedLorentzianProb(val,m,sd,minval,maxval);}
    else
      {p = GaussProb(val,m,sd);}
    double w = DMeanModel->dprob(m);
    sump += p * w;
    sumw += w;
  }	  
  //^
  if (sump/sumw < 0.0) {
    std::cout << "Negative probability! " << sump << " " << sumw << "\n";
  }
  return sump/sumw;
}
//--------------------------------------------------------------
  double ProbBiassed(const double& val, const double& sd,
		   const double& posmean, const MeanModelBase& DModel)
  {
    double pp = GaussProb(val, posmean, sd);
    double pm = IntgrtProb(DModel).Prob(val, sd);
    return pp/(pp+pm);
  }		   
//--------------------------------------------------------------
  double ProbUnbiassed(const double& val, const double& sd,
		    const double& posmean,const double& negmean)
  {
    double pp = GaussProb(val, posmean, sd);
    double pm = GaussProb(val, negmean, sd);
    return pp/(pp+pm);
  }
}
