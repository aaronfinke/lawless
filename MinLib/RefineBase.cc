//(c) 2000-2007 Cambridge University Technical Services Ltd
//All rights reserved
#include "RefineBase.h"

namespace phaser {

//------------------
// VIRTUAL FUNCTIONS
//------------------
void        RefineBase::cleanUp(outStream s,Output& o) {}
void        RefineBase::rejectOutliers(outStream s,Output& o) {}
repars1D    RefineBase::getRepar() { return repars1D(0); }
double      RefineBase::getMaxDistSpecial(TNT::Vector<double>&, TNT::Vector<double>&, TNT::Vector<double>&) { return std::numeric_limits<double>::max()/2.; }
bounds1D    RefineBase::getLowerBounds() { return bounds1D(numRefinePars()); }
bounds1D    RefineBase::getUpperBounds() { return bounds1D(numRefinePars()); }
void        RefineBase::logProtocolPars(outStream s,Output& o) { }

void RefineBase::setProtocol(protocolPtr protocol,Output& output)
{
  refinePar = getRefineMask(protocol);
  npars_ref = 0;
  for (int i = 0; i < refinePar.size(); i++) if (refinePar[i]) npars_ref++;
  npars_all = refinePar.size();
}

//----------
// FUNCTIONS
//----------
int RefineBase::numRefinePars() { return npars_ref; }

// Return maximum multiple of gradient that can be shifted
// before hitting bounds.  Used to define limits for line search.
// Also return vector of allowed shifts for each parameter, which
// can be used to decide which parameters can't be moved along the
// search direction
double      RefineBase::getMaxDist(TNT::Vector<double>& x, TNT::Vector<double>& g, TNT::Vector<double>& dist)
{
  floatType maxDist(std::numeric_limits<floatType>::max()/2.),ZERO(0);
  bounds1D Lower = getLowerBounds();
  reparBounds(Lower);
  bounds1D Upper = getUpperBounds();
  reparBounds(Upper);
  for (int i=0; i < numRefinePars(); i++)
  {
    // Test upper bound for negative gradient (positive shift),
    // lower bound for positive gradient
    if (g[i] < 0 && Upper[i].bounded)
    {
      dist[i] = std::max(ZERO,(x[i]-Upper[i].limit)/g[i]);
      maxDist = std::min(maxDist,dist[i]);
      // x[i] = std::min(x[i],Upper[i].limit); // Paranoia: keep in bounds
    }
    else if (g[i] > 0 && Lower[i].bounded)
    {
      dist[i] = std::max(ZERO,(x[i]-Lower[i].limit)/g[i]);
      maxDist = std::min(maxDist,dist[i]);
      // x[i] = std::max(x[i],Lower[i].limit); // Paranoia: keep in bounds
    }
    else
      dist[i] = -1.0; // Flag for parameter that is unbounded or not moving
  }
  // Normal Upper and Lower limits should not be set for correlated parameters,
  // which are checked here instead.
  maxDist = std::min(maxDist,getMaxDistSpecial(x,g,dist));
  return maxDist;
}

// Return maximum multiple of gradient that can be shifted before exceeding
// large shift limit. Used to limit size of first step in line search.
double RefineBase::getMaxStep(TNT::Vector<floatType>& x, TNT::Vector<floatType>& g)
{
  floatType maxStep(std::numeric_limits<floatType>::max());
  TNT::Vector<floatType> largeShifts = getLargeShifts();
  reparLargeShifts(largeShifts);
  for (int i=0; i < numRefinePars(); i++)
    if (g[i])
    {
      floatType step = std::fabs(largeShifts[i]/g[i]);
      maxStep = std::min(maxStep,step);
    }
  PHASER_ASSERT(maxStep > 0);
  return maxStep;
}

// Used when start of line search is at boundary, heading out.  Reset gradient
// for parameter(s) at boundary to zero. Return number of parameters that are still moving.
int RefineBase::filterGradient(TNT::Vector<floatType>& g, TNT::Vector<floatType>& dist, floatType dfilt)
{
  int nleft(numRefinePars());
  for (int i=0; i < numRefinePars(); i++)
  {
    if (dist[i] >= 0. && dist[i] < dfilt)
    {
      g[i] = 0.;
    }
    if (g[i] == 0.)  nleft--;
  }
  return nleft;
}

void RefineBase::logVector(outStream where,std::string what, Output& output,TNT::Vector<floatType>& vec)
{
  output.logTab(1,where,what);
  for (int i = 0; i < numRefinePars(); i++)
    output.logTabPrintf(1,where,"%3d [% 5.3e] : %s\n", i+1,vec[i],whatAmI(i).c_str());
  output.logBlank(where);
}

void RefineBase::logHessian(outStream where,Output& output,TNT::Fortran_Matrix<floatType>& Hessian)
{
  int max_dim(11);
  output.logTab(1,where,"Hessian for refined parameters - Matrix");
  if (Hessian.num_rows() >= max_dim)
    output.logTab(1,where,"Hessian is too large to write to log: size = "+itos(Hessian.num_rows()));
  for (int i=1; i<=Hessian.num_rows() && i < max_dim; i++)
  {
    output.logTabPrintf(1,where,"[");
    for (int j=1; j<=Hessian.num_cols() && j < max_dim; j++)
      if (Hessian(i,j)) output.logTabPrintf(0,where," % 1.0e", Hessian(i,j)) ;
      else output.logTabPrintf(0,where," --0-- ") ;
    (Hessian.num_rows() >= max_dim) ?  output.logTabPrintf(0,where," etc...\n") : output.logTabPrintf(0,where," ]\n");
  }
  if (Hessian.num_rows() >= max_dim) output.logTabPrintf(0,where," etc...\n");
  output.logBlank(where);
}

bool RefineBase::fixHessian(TNT::Fortran_Matrix<floatType>& Hessian)
// Eliminate non-positive curvatures.
// Return true if any curvatures were non-positive.
//
// Determine mean factor by which curvatures differ from 1/largeShift^2
// from parameters with positive curvatures.  This gives the average
// curvature that would be obtained if all the parameters were scaled
// so that their largeShift values were one.
// Expected curvatures (consistent with relative largeShift values) are
// computed from this and used to replace non-positive curvatures.
// Corresponding off-diagonal terms are set to zero.
// RJR Note 31/5/06: tried setting a minimum fraction of the expected
// curvature, but this degraded convergence and depended too much on
// accurate estimation of largeShift.
{
  floatType meanfac(0);
  bool needNewHessian(false);
  int npos(0);
  TNT::Vector<floatType> largeShifts = getLargeShifts();
  reparLargeShifts(largeShifts);
  for (int i = 0; i < numRefinePars(); i++)
  {
    if (Hessian(i+1,i+1) > 0.)
    {
      npos++;
      meanfac += Hessian(i+1,i+1)*fn::pow2(largeShifts[i]);
    }
  }
  if (npos)
  {
    meanfac /= npos;
    for (int i = 0; i < numRefinePars(); i++)
    {
      if (Hessian(i+1,i+1) <= 0)
      {
        Hessian(i+1,i+1) = meanfac/fn::pow2(largeShifts[i]);
        for (int j = 0; j < numRefinePars(); j++)
          if (i != j) Hessian(i+1,j+1) = Hessian(j+1,i+1) = 0.;
      }
    }
  }
  else // Fall back on diagonal matrix based on largeShifts shift values
  {
    for (int i = 0; i < numRefinePars(); i++)
    {
      Hessian(i+1,i+1) = 100./fn::pow2(largeShifts[i]);
      for (int j = 0; j < numRefinePars(); j++)
        if (i != j) Hessian(i+1,j+1) = Hessian(j+1,i+1) = 0.;
    }
  }
  if (numRefinePars()-npos) needNewHessian = true;
  return needNewHessian;
}

void RefineBase::studyParams(Output& output)
// Vary refined parameters and print out function, gradient and curvature
{
  std::string filename("studyParams");
  TNT::Vector<floatType> x(numRefinePars()),oldx(numRefinePars()),g(numRefinePars());
  TNT::Vector<floatType> unrepar_x(numRefinePars()),unrepar_oldx(numRefinePars()),unrepar_g(numRefinePars());
  TNT::Fortran_Matrix<floatType> Hessian(numRefinePars(),numRefinePars());
  bounds1D Lower = getLowerBounds();
  reparBounds(Lower);
  bounds1D Upper = getUpperBounds();
  reparBounds(Upper);
  TNT::Vector<floatType> largeShifts = getLargeShifts();
  reparLargeShifts(largeShifts);

  output.logBlank(LOGFILE);
  output.logTab(0,LOGFILE,"Study behaviour of parameters near current values\n");
  unrepar_x = getRefinePars();
  oldx = reparRefinePars(unrepar_x);
  x = oldx;
  applyShift(unrepar_x);
  floatType fmin(targetFn());

  // First check that function values from targetFn, gradientFn and hessianFn agree
  floatType gradLogLike = gradientFn(unrepar_g);
  g = reparGradient(unrepar_x,unrepar_g);
  if (std::abs(fmin-gradLogLike) >= 0.1)
  {
    output.logTab(1,LOGFILE,"Likelihoods from targetFn and gradientFn disagree");
    output.logTab(2,LOGFILE,"Likelihood from targetFn: " + dtos(-fmin));
    output.logTab(2,LOGFILE,"Likelihood from   gradientFn: " + dtos(-gradLogLike));
  }
  bool is_diagonal;
  floatType hessLogLike = hessianFn(Hessian,is_diagonal);
  reparHessian(unrepar_x,unrepar_g,Hessian);
  if (std::abs(fmin-hessLogLike) >= 0.1)
  {
    output.logTab(1,LOGFILE,"Likelihoods from targetFn and hessianFn disagree");
    output.logTab(2,LOGFILE,"Likelihood from  targetFn: " + dtos(-fmin));
    output.logTab(2,LOGFILE,"Likelihood from hessianFn: " + dtos(-hessLogLike));
  }
  PHASER_ASSERT(std::abs(fmin-gradLogLike) < 0.1);
  PHASER_ASSERT(std::abs(fmin-hessLogLike) < 0.1);
  fmin = gradLogLike; // Use function value from gradient for consistency below

  int ndiv(10);
  for (int i = 0; i < numRefinePars(); i++)
  {
    std::ofstream outstream;
    outstream.open(std::string(filename+"."+itos(i+1)).c_str());
    output.logBlank(LOGFILE);
    output.logTabPrintf(1,LOGFILE,"Refined Parameter #%d (%s)\n",i+1,whatAmI(i).c_str());
    output.logTabPrintf(1,LOGFILE,"Centered on %g\n",oldx[i]);
    floatType xmin = oldx[i] - 2.*largeShifts[i];
    if (Lower[i].bounded)
    {
      xmin = std::max(xmin,Lower[i].limit);
      output.logTabPrintf(1,LOGFILE,"Lower limit: %g\n",Lower[i].limit);
    }
    else // "Unbounded" parameters may have limits depending on correlations
    {
      TNT::Vector<floatType> g_fake(numRefinePars(),0.);
      TNT::Vector<floatType> d_fake(numRefinePars());
      g_fake[i] = largeShifts[i]/10.; // Put in range of plausible shift
      floatType specialLimit(oldx[i]-getMaxDist(oldx,g_fake,d_fake)*g_fake[i]);
      if (specialLimit > xmin)
      {
        xmin = specialLimit;
        output.logTabPrintf(1,LOGFILE,"Special lower limit: %g\n",specialLimit);
      }
    }
    floatType xmax = oldx[i] + 2.*largeShifts[i];
    if (Upper[i].bounded)
    {
      xmax = std::min(xmax,Upper[i].limit);
      output.logTabPrintf(1,LOGFILE,"Upper limit: %g\n",Upper[i].limit);
    }
    else
    {
      TNT::Vector<floatType> g_fake(numRefinePars(),0.);
      TNT::Vector<floatType> d_fake(numRefinePars());
      g_fake[i] = -largeShifts[i]/10.;
      floatType specialLimit(oldx[i]-getMaxDist(oldx,g_fake,d_fake)*g_fake[i]);
      if (specialLimit < xmax)
      {
        xmax = specialLimit;
        output.logTabPrintf(1,LOGFILE,"Special upper limit: %g\n",specialLimit);
      }
    }
    output.logTabPrintf(1,LOGFILE,"Large shift: %g\n",largeShifts[i]);
    output.logTab(1,LOGFILE," parameter   func-min    gradient   curvature  crv*lrg^2   FD grad     FD curv\n\n");
    float1D fval(ndiv+1);
    floatType dx((xmax-xmin)/ndiv),df,lastdf,lastg,lasth;
    for (int j = 0; j <= ndiv; j++)
    {
      x[i] = xmin + j*dx;
      TNT::Vector<double> unrepar_x = reparRefineParsInv(x);
      applyShift(unrepar_x);
      floatType gradLogLike = gradientFn(unrepar_g);
      g = reparGradient(unrepar_x,unrepar_g);
      df = gradLogLike - fmin;
      fval[j] = df;
      bool is_diagonal;
      hessianFn(Hessian,is_diagonal);
      reparHessian(unrepar_x,unrepar_g,Hessian);
      if (j > 1)
      {
        floatType fdgrad((fval[j]-fval[j-2])/(2.*dx));
        floatType fdhess((fval[j]+fval[j-2]-2.*fval[j-1])/fn::pow2(dx));
        output.logTabPrintf(1,LOGFILE,"%+.4e %+.3e %+.4e %+.4e %+.2e %+.4e %+.4e\n",x[i]-dx,lastdf,lastg,lasth,lasth*fn::pow2(largeShifts[i]),fdgrad,fdhess);
      }
      if (j == 0 || j == ndiv)
        output.logTabPrintf(1,LOGFILE,"%+.4e %+.3e %+.4e %+.4e %+.2e\n",x[i],df,g[i],Hessian(i+1,i+1),Hessian(i+1,i+1)*fn::pow2(largeShifts[i]));
      outstream << x[i] << " " << df << " " << g[i] << " " << Hessian(i+1,i+1) << " " << std::endl;
      lastdf = df;
      lastg = g[i];
      lasth = Hessian(i+1,i+1);
    }
    x[i] = oldx[i];
    outstream.close();
  }
  output.logBlank(LOGFILE);
  unrepar_oldx = reparRefineParsInv(oldx);
  applyShift(unrepar_oldx);
}

TNT::Vector<double> RefineBase::reparRefinePars(TNT::Vector<floatType> pars)
{
  std::vector<reparams> repar = getRepar();
  for (int i = 0; i < numRefinePars() && repar.size(); i++)
    if (repar[i].reparamed)
      pars[i] = std::log(pars[i] + repar[i].offset);
  return pars;
}

TNT::Vector<double> RefineBase::reparRefineParsInv(TNT::Vector<floatType> pars)
{
  std::vector<reparams> repar = getRepar();
  for (int i = 0; i < numRefinePars() && repar.size(); i++)
    if (repar[i].reparamed)
      pars[i] = std::exp(pars[i]) - repar[i].offset;
  return pars;
}

void RefineBase::reparBounds(bounds1D& bounds)
{
  TNT::Vector<floatType> x(npars_ref);
  for (int i = 0; i < numRefinePars(); i++) x[i] = bounds[i].limit;
  x = reparRefinePars(x);
  for (int i = 0; i < numRefinePars(); i++)  bounds[i].limit = x[i];
}

TNT::Vector<floatType> RefineBase::reparGradient(TNT::Vector<double>& pars,TNT::Vector<floatType> gradient)
{
  //This is only valid for y=ln(x+c) transformation.
  //f(y(x))
  //df    df   dx   df
  //-- =  -- * -- = -- (x+c)
  //dy    dx   dy   dx

  std::vector<reparams> repar = getRepar();
  for (int i = 0; i < numRefinePars() && repar.size(); i++)
    if (repar[i].reparamed)
      gradient[i] *= (pars[i]+repar[i].offset);
  return gradient;
}

void RefineBase::reparHessian(TNT::Vector<floatType>& pars,TNT::Vector<floatType>& gradient,TNT::Fortran_Matrix<floatType>& Hessian)
{
  //This is only valid for y=ln(x+c) transformation.
  std::vector<reparams> repar = getRepar();
  for (int i = 0; i < numRefinePars() && repar.size(); i++)
    for (int j = 0; j < numRefinePars(); j++)
      if (repar[i].reparamed || repar[j].reparamed)
      {
        floatType dxidyi(1),dxjdyj(1);
        if (repar[i].reparamed) dxidyi = pars[i]+repar[i].offset;
        if (repar[j].reparamed) dxjdyj = pars[j]+repar[j].offset;
        Hessian(i+1,j+1) *= dxidyi*dxjdyj;
        if (i==j) Hessian(i+1,i+1) += gradient[i]*dxidyi; //dxidyi=d2xidyi2
      }
}

void RefineBase::reparLargeShifts(TNT::Vector<floatType>& largeShifts)
{
  bounds1D oriLower = getLowerBounds();
  bounds1D Lower = getLowerBounds();
  reparBounds(Lower);
  std::vector<reparams> repar = getRepar();
  for (int i = 0; i < numRefinePars() && repar.size(); i++)
    if (repar[i].reparamed)
    {
      PHASER_ASSERT(Lower[i].bounded);
      largeShifts[i] = std::log(oriLower[i].limit + largeShifts[i] + repar[i].offset) - Lower[i].limit;
    }
}

} //phaser
