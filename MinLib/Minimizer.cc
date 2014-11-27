//(c) 2000-2007 Cambridge University Technical Services Ltd
//All rights reserved
#include "Minimizer.h"
#include "symm_eigen.h"
#include "Errors.h"
#include "jiffy.h"
#include "maths.h"
#include <cmath>

namespace phaser {

void Minimizer::run(RefineBase& target,protocolPtr p,Output& output)
{
  af::shared<protocolPtr> protocol;
  //^
  ///  std::cout << " Minimizer::run::Number of cycles " << p->getNCYC() << "\n";

  protocol.push_back(p);
  run(target,protocol,output);
}

void Minimizer::run(RefineBase& target,af::shared<protocolPtr> protocol,Output& output,bool clean_up,bool study_params)
{
  //this check is also done on Input
  if (!protocol.size()) PhaserError(FATAL,"No protocol for refinement");

  output.logTab(1,VERBOSE,"Values of parameters will be reported at end of each cycle");
  output.logBlank(VERBOSE);

  // Main loop
  bool calcNewHessian(false),needNewHessian(false);
  int ncyc_since(0);

  // Allocate memory for variables needed by the local minimizers
  floatType f(0),firstLogLike(0),oldLogLike(0),requiredGain(0);

  bool some_refinement(false);
  for (unsigned big_cyc = 0 ; big_cyc < protocol.size(); big_cyc++)
    if (!protocol[big_cyc]->is_off()) some_refinement = true;
  if (!some_refinement)
  {
    output.logTab(1,LOGFILE,"No refinement of parameters");
    output.logBlank(LOGFILE);
  }

  //for bfgs
  TNT::Vector<floatType> x_old;
  TNT::Vector<floatType> g_old;
  TNT::Fortran_Matrix<floatType> h_old;

  if (some_refinement)
  for (unsigned big_cyc = 0 ; big_cyc < protocol.size(); big_cyc++)
  {
    if (protocol.size() > 1) {
      output.logUnderLine(LOGFILE,"MACROCYCLE #" + itos(big_cyc+1) + " OF " + itos(protocol.size()));
    }
    target.setProtocol(protocol[big_cyc],output);
    if (!protocol[big_cyc]->getNCYC())
    {
      output.logTab(1,LOGFILE,"No cycles of refinement this macrocycle");
      output.logBlank(LOGFILE);
      continue; //skip rest of code in loop
    }
    else if (!target.numRefinePars())
    {
      output.logTab(1,LOGFILE,"No parameters to refine this macrocycle");
      output.logBlank(LOGFILE);
      continue; //skip rest of code in loop
    }
    else
    {
      output.logTabPrintf(0,LOGFILE,"%s",protocol[big_cyc]->logfile().c_str());
      output.logBlank(LOGFILE);
      target.logProtocolPars(LOGFILE,output);
    }

    target.rejectOutliers(LOGFILE,output);

    for (unsigned small_cyc = 0; small_cyc < protocol[big_cyc]->getNCYC(); small_cyc++)
    {
      output.logUnderLine(VERBOSE,"MACROCYCLE #" + itos(big_cyc+1) + ": CYCLE #" + itos(small_cyc+1));

      //store old value
      oldLogLike = f;
      //calculate new value
      f = target.targetFn();
      //Store value of refined parameters from start of this loop
      TNT::Vector<floatType> xbefore(target.numRefinePars());
      xbefore = target.getRefinePars();
      xbefore = target.reparRefinePars(xbefore);

      //if this is first call of the function, initialize reference function values
      if (!small_cyc) firstLogLike = oldLogLike = f;
      if (!small_cyc)
      {
        output.logTab(1,LOGFILE,"Optimization statistics macrocycle #" + itos(big_cyc+1));
        output.logTabPrintf(1,LOGFILE,"Cycle  %-15s %-18s %-18s\n","end-this-cycle","change-from-start","change-from-last");
        if (std::fabs(firstLogLike) > 10e-3)
        output.logTabPrintf(1,LOGFILE,"start  %13.3f %13.3s %18.3s\n",-firstLogLike,"","");
        else
        output.logTabPrintf(1,LOGFILE,"start  %13.3e %13.3s %18.3s\n",-firstLogLike,"","");
      }
      //determine gain required not to converge in this cycle
      const floatType EPS(1.0e-10);
      const floatType FTOL(1.0e-6);
      const floatType TWO(2.0);
      // Convergence test based on numerical precision
      requiredGain = FTOL*(fabs(f)+fabs(oldLogLike)+EPS)/TWO;

      // Call the minimizer, return the new function value, parameter shifts implicit

      bool hitBound(false),tooSmallShift(false);
      if (protocol[big_cyc]->getMINIMIZER() == "NEWTON")
      {
        f = newton(target,output,requiredGain,hitBound,tooSmallShift);
      }
      else if (protocol[big_cyc]->getMINIMIZER() == "BFGS")
      {
        if (!small_cyc) calcNewHessian = true;
        f = bfgs(target,output,requiredGain,calcNewHessian,needNewHessian,hitBound,tooSmallShift,x_old,g_old,h_old);
        if (calcNewHessian)
        {
          ncyc_since = 0;
          calcNewHessian = false;
        }
      }
      else if (protocol[big_cyc]->getMINIMIZER() == "DESCENT")
      {
        f = descent(target,output,requiredGain,hitBound,tooSmallShift);
      }

      target.logCurrent(VERBOSE,output); //current parameters

  //  Print out information detailing the minimization progress
      if (std::fabs(firstLogLike) > 10e-3)
      output.logTabPrintf(1,LOGFILE,"#%-5d %13.3f %13.3f %18.3f\n",small_cyc+1,-f,firstLogLike-f,oldLogLike-f);
      else
      output.logTabPrintf(1,LOGFILE,"#%-5d %13.3e %13.3e %18.3e\n",small_cyc+1,-f,firstLogLike-f,oldLogLike-f);
      TNT::Vector<floatType> largeShifts = target.getLargeShifts();
      // See whether parameters have moved
      TNT::Vector<floatType> xafter(target.numRefinePars());
      xafter = target.getRefinePars();
      xafter = target.reparRefinePars(xafter);
      bool zero_shift(true);
      const floatType ZERO(0.0);
      for (int i=0; i < target.numRefinePars(); i++)
        if ((xafter[i]-xbefore[i]) != ZERO) zero_shift = false;

      if ( zero_shift || tooSmallShift || (oldLogLike-f < requiredGain) )
      {
        if (hitBound)
        {
          output.logTab(1,DEBUG,"Hit bound in last search: carry on");        }
        else if (needNewHessian && ncyc_since)
        {
          output.logTab(1,DEBUG,"Recalculate Hessian and carry on");
          calcNewHessian = true;
          needNewHessian = false;
        }
        else
        {
          if (protocol.size() > 1) {
            output.logBlank(LOGFILE);
            output.logTab(1,LOGFILE,"---CONVERGENCE OF MACROCYCLE---");
            output.logBlank(LOGFILE);
            target.logCurrent(LOGFILE,output); //current parameters
          }
          break;
        }
      }
      if (small_cyc == protocol[big_cyc]->getNCYC()-1)
      {
        output.logBlank(LOGFILE);
        output.logTab(1,LOGFILE,"---ITERATION LIMIT OF MACROCYCLE---");
        output.logBlank(LOGFILE);
        target.logCurrent(LOGFILE,output); //current parameters
        break;
      }
      ncyc_since++;
      output.logBlank(VERBOSE);
    } //end cycle loop

    // study behaviour of function, gradient and curvature as parameters varied around minimum
    // using mathematica file studyParams_short.nb
    if (study_params || (getenv("PHASER_STUDY_PARAMS") != NULL)) target.studyParams(output);

    if (clean_up) target.cleanUp(LOGFILE,output); //clean up ready for next macrocycle

  } //end macrocycle loop
}


floatType Minimizer::descent(RefineBase& target,Output& output,floatType requiredGain,bool& hitBound, bool& tooSmallShift)
{
  floatType f(0.);
  const floatType ZERO(0.), MinDXoverESD(0.1);
  output.logTab(1,DEBUG,"== STEEPEST DESCENT ==");

  //get initial parameter values  x
  TNT::Vector<floatType> unrepar_x = target.getRefinePars();
  TNT::Vector<floatType> x = target.reparRefinePars(unrepar_x);

  // Compute gradient.
  // If it is zero (pathological case), return current likelihood.
  output.logTab(1,DEBUG,"Descent: Gradient");
  output.logEllipsisStart(DEBUG,"Descent: Calculating the gradient");
  TNT::Vector<floatType> unrepar_g(target.numRefinePars());
  floatType gradLogLike = target.gradientFn(unrepar_g);
  bool zero_gradient(true);
  for (int i = 0; i<target.numRefinePars(); i++)
    if (unrepar_g[i] != ZERO) zero_gradient = false;
  if (zero_gradient)
  {
    output.logTab(1,VERBOSE,"Descent: Zero gradient: exiting");
    hitBound = false;
    f = target.targetFn();
    return f;
  }

  TNT::Vector<floatType> g = target.reparGradient(unrepar_x,unrepar_g);
  output.logEllipsisEnd(DEBUG);
  output.logBlank(DEBUG);
  target.logVector(DEBUG,"Descent: Gradient",output,g);
  output.logBlank(DEBUG);

  // Scale the gradient by factor proportional to large shifts squared.
  // (Hessian is inversely proportional to scale of parameter squared.)
  // Additional factor of 1/100 approximates gradient over curvature,
  // with the largeShift values chosen by behaviour of function around
  // minimum.
  TNT::Vector<floatType> largeShifts = target.getLargeShifts();
  target.reparLargeShifts(largeShifts);
  TNT::Vector<floatType> gs(target.numRefinePars());
  for (int i =0; i < g.size(); i++) gs[i] = g[i]*0.01*fn::pow2(largeShifts[i]);

  output.logBlank(DEBUG);
  target.logVector(DEBUG,"Descent: Scaled Gradient",output,g);
  output.logBlank(DEBUG);

  // Figure out how far linesearch has to go to shift at least one parameter
  // by minimum shift/esd, using largeShift/10 as rough estimate of esd
  floatType requiredShift(0.);
  for (int i = 0; i < g.size(); i++)
    requiredShift = std::max(requiredShift,fabs(gs[i])/(largeShifts[i]/10.));
  requiredShift = MinDXoverESD/requiredShift;

  f = LineSearch(target,output,x,gradLogLike,g,gs,1.,requiredGain,requiredShift,hitBound,tooSmallShift);

  return f;
}

floatType Minimizer::bfgs(RefineBase& target,Output& output,floatType requiredGain,
           bool& calcNewHessian,bool& needNewHessian,bool& hitBound,bool& tooSmallShift, TNT::Vector<floatType>& x_old,
           TNT::Vector<floatType>& g_old,TNT::Fortran_Matrix<floatType>& h_old)
{
  floatType f;
  const floatType ZERO(0.), MinDXoverESD(0.1);
  output.logTab(1,DEBUG,"== BFGS ==");

  if (calcNewHessian && (x_old.size() != target.numRefinePars()))
  {
    // Resize if the number of parameters has changed
    x_old.newsize(target.numRefinePars());
    g_old.newsize(target.numRefinePars());
    h_old.newsize(target.numRefinePars(),target.numRefinePars());
  }

  //get initial parameter values x
  TNT::Vector<floatType> unrepar_x = target.getRefinePars();
  TNT::Vector<floatType> x = target.reparRefinePars(unrepar_x);

  // Compute gradient.
  // If gradient is zero (pathological case), return current likelihood.
  output.logTab(1,DEBUG,"BFGS: Gradient");
  output.logEllipsisStart(DEBUG,"BFGS: Calculating the gradient");
  TNT::Vector<floatType> unrepar_g(target.numRefinePars());
  floatType gradLogLike = target.gradientFn(unrepar_g);
  bool zero_gradient(true);
  for (int i = 0; i<target.numRefinePars(); i++)
    if (unrepar_g[i] != ZERO) zero_gradient = false;
  if (zero_gradient)
  {
    output.logTab(1,VERBOSE,"BFGS: Zero gradient: exiting");
    calcNewHessian = false;
    needNewHessian = false;
    hitBound = false;
    f = target.targetFn();
    return f;
  }

  TNT::Vector<floatType> g = target.reparGradient(unrepar_x,unrepar_g);
  output.logEllipsisEnd(DEBUG);
  output.logBlank(DEBUG);
  target.logVector(DEBUG,"BFGS: Gradient",output,g);
  output.logBlank(DEBUG);

  if (calcNewHessian)
  {
    output.logTab(1,DEBUG,"== Newton Step ==");
    output.logBlank(DEBUG);

    bool is_diagonal;
    target.hessianFn(h_old,is_diagonal);
    target.reparHessian(unrepar_x,unrepar_g,h_old);
    target.logHessian(DEBUG,output,h_old);
    needNewHessian = target.fixHessian(h_old);
    if (needNewHessian)
    {
      output.logTab(1,DEBUG,"BFGS::Newton: Fixed curvature(s) below expected minimum in Hessian");
      target.logHessian(DEBUG,output,h_old);
    }

    // Set up to repeat perturbed Newton step if insufficient progress
    int ntry(0);
    floatType hiiRMS(ZERO);
    TNT::Vector<floatType> hii(target.numRefinePars());
    for (int i = 0; i < g.size(); i++)
    {
      hii[i] = h_old(i+1,i+1);
      hiiRMS += fn::pow2(hii[i]);
    }
    if (g.size()) hiiRMS = std::sqrt(hiiRMS/g.size());
    int filtered_last(0);
    f = gradLogLike;
    hitBound = false;
    x_old = x;
    while (ntry < 3 && !hitBound && gradLogLike-f < requiredGain && filtered_last < target.numRefinePars())
    {
      ntry++;
      if (is_diagonal)
      {
        for (int i = 0; i < g.size(); i++)
          h_old(i+1,i+1) = 1/(hii[i]+(ntry-1)*hiiRMS); // downweight low curvature later
      }
      else
      {
        //  Scale Hessian by pre- and post-multiplication by diagonal matrix to
        //  have unit diagonal prior to computing pseudoinverse. Scale resulting
        //  pseudoinverse again to put back on original scale.
        TNT::Fortran_Matrix<floatType> h(target.numRefinePars(),target.numRefinePars());
        TNT::Vector<floatType> hscale(target.numRefinePars());
        for (int i = 0; i < g.size(); i++)
          hscale[i] = 1/sqrt(h_old(i+1,i+1));
        for (int i = 0; i < g.size(); i++)
          for (int j = 0; j < g.size(); j++)
            h(i+1,j+1) = hscale[i]*h_old(i+1,j+1)*hscale[j];
        int min_to_filter(0);
        if (ntry > 1) // First attempt did not minimize.  Perturb Hessian and try again.
        {
          min_to_filter = filtered_last + std::max(1,target.numRefinePars()/10);
          if (min_to_filter >= target.numRefinePars()) {
            min_to_filter = filtered_last + 1;
          }
          if (min_to_filter >= target.numRefinePars()) {
            min_to_filter = target.numRefinePars() - 1; // force less than number of parameters
          }
        }
        int filtered(0);
        h_old = SymmetricPseudoinverse<floatType>(h,filtered,false,min_to_filter).getInv();
        filtered_last = filtered;
        if (filtered) output.logTabPrintf(1,DEBUG,"Filtered %i small eigenvectors\n",filtered);
        for (int i = 0; i < g.size(); i++)
          for (int j = 0; j < g.size(); j++)
            h_old(i+1,j+1) = hscale[i]*h_old(i+1,j+1)*hscale[j];
        output.logTab(1,DEBUG,"Pseudoinverse Hessian");
        target.logHessian(DEBUG,output,h_old);
      }

      TNT::Vector<floatType> gs = h_old*g;

      target.logVector(DEBUG,"BFGS::Newton: Scaled Gradient",output,gs);
      output.logBlank(DEBUG);

      // Figure out how far linesearch has to go to shift at least one parameter
      // by minimum shift/esd, using inverse Hessian to estimate covariance matrix
      floatType requiredShift(0.);
      for (int i = 0; i < g.size(); i++)
      {
        if (h_old(i+1,i+1) > 0)
          requiredShift = std::max(requiredShift,fabs(gs[i])/std::sqrt(h_old(i+1,i+1)));
      }
      if (requiredShift > 0)
        requiredShift = MinDXoverESD/requiredShift;

      f = LineSearch(target,output,x,gradLogLike,g,gs,1.,requiredGain,requiredShift,hitBound,tooSmallShift);

      //LineSearch sets scaled gradient to zero on boundaries.
      //Do same to gradient before BFGS update.
      for (int i = 0; i < g.size(); i++) if (gs[i] == 0) g[i] = 0;
      g_old = g;
    }
  }
  else
  {
    output.logTab(1,DEBUG,"== BFGS Step ==");
    TNT::Vector<floatType> dx(target.numRefinePars());
    TNT::Vector<floatType> dg(target.numRefinePars());
    TNT::Vector<floatType> Hdg(target.numRefinePars());
    TNT::Vector<floatType> u(target.numRefinePars());
    dg  = g-g_old;
    dx  = x-x_old;
    Hdg = h_old*dg;
    target.logVector(DEBUG,"BFGS: dg",output,dg);
    target.logVector(DEBUG,"BFGS: dx",output,dx);
    target.logVector(DEBUG,"BFGS: Hdg",output,Hdg);
    floatType dx_dot_dg(0);
    floatType dg_dot_Hdg(0);
    for (int i = 0; i < g.size(); i++) dx_dot_dg += dx[i]*dg[i];
    for (int i = 0; i < g.size(); i++) dg_dot_Hdg += dg[i]*Hdg[i];
    output.logTab(1,DEBUG,"BFGS: dx_dot_dg, dg_dot_Hdg: " + dtos(dx_dot_dg) + "  " + dtos(dg_dot_Hdg));
    if (dx_dot_dg <= 0)  // Shouldn't happen, but we're cutting corners on line search convergence
    {                    // i.e. Wolfe condition considering new gradient
      output.logTab(1,VERBOSE,"BFGS: dx_dot_dg <= 0: restart with new Hessian");
      needNewHessian = true;
      f = target.targetFn();
      return f;
    }
    PHASER_ASSERT(dg_dot_Hdg);
    for (int i = 0; i < g.size(); i++)
      u[i] = dx[i]/dx_dot_dg - Hdg[i]/dg_dot_Hdg;
    for (int i = 0; i < g.size(); i++)
      for (int j = 0; j < g.size(); j++)
        h_old(i+1,j+1) += dx[i]*dx[j]/dx_dot_dg - Hdg[i]*Hdg[j]/dg_dot_Hdg + dg_dot_Hdg*u[i]*u[j];

    output.logBlank(DEBUG);
    target.logHessian(DEBUG,output,h_old);
    output.logBlank(DEBUG);

    g_old = g;
    x_old = x;
    TNT::Vector<floatType> gs = h_old*g;

    target.logVector(DEBUG,"BFGS: Scaled Gradient",output,gs);
    output.logBlank(DEBUG);

    // Figure out how far linesearch has to go to shift at least one parameter
    // by minimum shift/esd, using inverse Hessian to estimate covariance matrix
    floatType requiredShift(0.);
    for (int i = 0; i < g.size(); i++)
    {
      if (h_old(i+1,i+1) > 0)
        requiredShift = std::max(requiredShift,fabs(gs[i])/std::sqrt(h_old(i+1,i+1)));
    }
    if (requiredShift > 0)
      requiredShift = MinDXoverESD/requiredShift;

    f = LineSearch(target,output,x,gradLogLike,g,gs,1,requiredGain,requiredShift,hitBound,tooSmallShift);

    // Check for positive grad (possibly after modifying at bounds)
    floatType grad(-dot_prod(g,gs));
    if (grad >= 0.) needNewHessian = true; // Try again with new Hessian

    // Reset gradient before update, if any scaled gradient values were zeroed.
    for (int i = 0; i < g.size(); i++) if (gs[i] == ZERO) g[i] = ZERO;
  }
  return f;
}

floatType Minimizer::newton(RefineBase& target,Output& output,floatType requiredGain,bool& hitBound,bool& tooSmallShift)
{
  floatType f;
  const floatType ZERO(0.), MinDXoverESD(0.1);
  output.logTab(1,DEBUG,"==NEWTON==");

  //get initial parameter values x
  TNT::Vector<floatType> unrepar_x = target.getRefinePars();
  TNT::Vector<floatType> x = target.reparRefinePars(unrepar_x);

  // Compute gradient.
  // If gradient is zero (pathological case), return current likelihood.
  output.logTab(1,DEBUG,"Newton: Gradient");
  output.logEllipsisStart(DEBUG,"Newton: Calculating the gradient");
  TNT::Vector<floatType> unrepar_g(target.numRefinePars());
  floatType gradLogLike = target.gradientFn(unrepar_g);
  bool zero_gradient(true);
  for (int i = 0; i<target.numRefinePars(); i++)
    if (unrepar_g[i] != ZERO) zero_gradient = false;
  if (zero_gradient)
  {
    output.logTab(1,VERBOSE,"Newton: Zero gradient: exiting");
    hitBound = false;
    f = target.targetFn();
    return f;
  }

  TNT::Vector<floatType> g = target.reparGradient(unrepar_x,unrepar_g);
  output.logEllipsisEnd(DEBUG);
  output.logBlank(DEBUG);
  target.logVector(DEBUG,"Newton: Gradient",output,g);
  output.logBlank(DEBUG);

  // compute hessian
  bool is_diagonal;
  TNT::Fortran_Matrix<floatType> h;
  target.hessianFn(h,is_diagonal);
  target.reparHessian(unrepar_x,unrepar_g,h);
  target.fixHessian(h);
  output.logTab(1,DEBUG,"Newton: BEFORE INVERSE");
  target.logHessian(DEBUG,output,h);

  // invert hessian
  if (is_diagonal)
  {
    for (int i = 1; i <= g.size(); i++)
      h(i,i) = 1/h(i,i);
  }
  else
  {
// Scale Hessian by pre- and post-multiplication by diagonal matrix to
// have unit diagonal prior to computing pseudoinverse.  Scale resulting
// pseudoinverse again to put back on original scale.
    TNT::Vector<floatType> hscale(target.numRefinePars());
    for (int i = 0; i < g.size(); i++)
      hscale[i] = 1/sqrt(h(i+1,i+1));
    for (int i = 0; i < g.size(); i++)
      for (int j = 0; j < g.size(); j++)
        h(i+1,j+1) = hscale[i]*h(i+1,j+1)*hscale[j];
    int filtered;
    h = SymmetricPseudoinverse<floatType>(h,filtered).getInv();
    if (filtered) output.logTabPrintf(1,DEBUG,"Newton: Filtered %i small eigenvectors\n",filtered);
    for (int i = 0; i < g.size(); i++)
      for (int j = 0; j < g.size(); j++)
        h(i+1,j+1) = hscale[i]*h(i+1,j+1)*hscale[j];
    output.logTab(1,DEBUG,"Newton: AFTER COMPUTING PSEUDOINVERSE");
    target.logHessian(DEBUG,output,h);
  }

  // this is the scaled gradient (this SHOULD be exact near the minimum)
  TNT::Vector<floatType> gs = h*g;

  target.logVector(DEBUG,"Newton: Parameters",output,x);
  target.logVector(DEBUG,"Newton: Scaled Gradient",output,gs);
  output.logBlank(DEBUG);

  // Figure out how far linesearch has to go to shift at least one parameter
  // by minimum shift/esd, using inverse Hessian to estimate covariance matrix
  floatType requiredShift(0.);
  for (int i = 0; i < g.size(); i++)
  {
    if (h(i+1,i+1) > 0)
      requiredShift = std::max(requiredShift,fabs(gs[i])/std::sqrt(h(i+1,i+1)));
  }
  if (requiredShift > 0)
    requiredShift = MinDXoverESD/requiredShift;

  f = LineSearch(target,output,x,gradLogLike,g,gs,1.,requiredGain,requiredShift,hitBound,tooSmallShift);
  return f;
}


void Minimizer::ShiftX(floatType a, RefineBase& target, TNT::Vector<floatType> &x,
                                TNT::Vector<floatType> &oldx, TNT::Vector<floatType> &gs)
{
  for (int i = 0; i < gs.size(); i++)
    x[i] = oldx[i]- a*gs[i];

  TNT::Vector<floatType> unrepar_x = target.reparRefineParsInv(x);
  target.applyShift(unrepar_x);
}

floatType Minimizer::ShiftScore(floatType a, RefineBase& target, TNT::Vector<floatType> &x,
  TNT::Vector<floatType> &oldx, TNT::Vector<floatType> &gs,
  Output* poutput, bool bcount)
{
  ShiftX(a,target,x,oldx,gs);
  floatType f = target.targetFn();
  if (bcount)
    mincount++;

  if (poutput)
    poutput->logTabPrintf(1,VERBOSE,"f(%9.6f) = %10.6f\n", a, f);

  return f;
}


floatType Minimizer::LineSearch(RefineBase& target,Output& output,TNT::Vector<floatType> oldx,
   floatType f, TNT::Vector<floatType> g, TNT::Vector<floatType>& gs, floatType starting_distance,
   floatType requiredGain, floatType requiredShift, bool& hitBound, bool& tooSmallShift)
{
  hitBound = false;
  tooSmallShift = false;
  mincount = 0;

  TNT::Vector<floatType> x(target.numRefinePars());
  TNT::Vector<floatType> unrepar_x(target.numRefinePars());

  const floatType ZERO(0.0),HALF(0.5),ONE(1.0),TWO(2.0),FIVE(5.0);
  const floatType GOLDEN((std::sqrt(FIVE)+ONE)/TWO);
  const floatType GOLDFRAC((GOLDEN-ONE)/GOLDEN);
  const floatType WOLFEC1(1.E-4); // 1.E-4 suggested in Nocedal & Wright
  const floatType DTOL(1.e-5);
  const floatType DFILT(DTOL*starting_distance);

  output.logTab(1,DEBUG,"Determining Stepsize");
  output.logTabPrintf(1,DEBUG,"|");

  // By default, first test point in bracketing is current guess of step
  // Make sure that it will not go past bounds
  TNT::Vector<floatType> dist(target.numRefinePars());
  floatType maxDist = target.getMaxDist(oldx,gs,dist);
  if (maxDist < DFILT) // If on boundary and moving out, set component of gradient to zero
  {
    int nleft = target.filterGradient(gs,dist,DFILT);
    if (!nleft)
    {
      tooSmallShift = true;
      return f;
    }
    maxDist = target.getMaxDist(oldx,gs,dist);
    PHASER_ASSERT(maxDist >= DFILT);
  }
  if (starting_distance > maxDist)
  {
    output.logTab(1,VERBOSE,"To avoid exceeding bounds, starting distance reduced from "
         + dtos(starting_distance) + " to " + dtos(maxDist));
    starting_distance = maxDist;
  }

  // Check grad of function in direction of (possibly modified) line search
  floatType grad(-dot_prod(g,gs));
  if (grad >= ZERO)
  {
    output.logTab(1,VERBOSE,"Grad of target in search direction >= 0: grad = " +dtos(grad));
    tooSmallShift = true;
    return f;
  }

  // Make sure that first test point does not exceed any largeShift
  floatType maxStep = target.getMaxStep(oldx,gs);
  if (starting_distance > maxStep)
  {
    output.logTab(1,VERBOSE,"To avoid exceeding largeShift, starting distance reduced from "
         + dtos(starting_distance) + " to " + dtos(maxStep));
    starting_distance = maxStep;
  }

  floatType fk,flo,fhi,dk,dlo,dhi;
  PHASER_ASSERT(gs.size() == oldx.size());

  //Sample first test point
  output.logTab(1,DEBUG,"Unit distance = " + dtos(starting_distance));

  dk  = starting_distance;
  fk = ShiftScore(dk, target, x, oldx, gs, &output);

  const floatType WOLFEFRAC(HALF); // 0.5 slightly better in tests of 0,0.25,0.5,1
  if ((starting_distance >= std::min(maxDist,WOLFEFRAC)) // Significant fraction of (quasi-)Newton shift
      && (fk < f+WOLFEC1*dk*grad)) // Wolfe condition for first step
  {
    output.logTab(1,VERBOSE,"Satisfied Wolfe condition for first step");
    output.logTabPrintf(1,VERBOSE,"Final LLG: %12.6f distance: %12.6f \n", -fk, dk);
    output.logTab(1,VERBOSE,"Bracketing took " + itos(mincount) + " function evaluations");

    if (dk >= maxDist)
      hitBound = true;

    return fk;
  }

  dlo = ZERO;
  flo = f;

  if (f <= fk) // First step is too big
  {
    while (f <= fk)
    { // Shrink shift until we find a value lower than starting point
      // Fit quadratic to grad and f at 0 and dk, choose its minimum (but not too small)
      if (dk < starting_distance/1.E10) // Give up if shift too small
      {
        unrepar_x = target.reparRefineParsInv(oldx);
        target.applyShift(unrepar_x);  // Reset x before returning
        tooSmallShift = true;
        return f;
      }
      fhi = fk;
      dhi = dk;
      floatType denom(TWO*(grad*dk+(f-fk)));
      if (denom < 0)
        dk *= std::min(0.99,std::max(0.1,grad*dk/denom));
      else // unlikely but possible for vanishingly small gradient
        dk *= 0.1;
      fk = ShiftScore(dk, target, x, oldx, gs, &output);
    }
    if ((fk < f+WOLFEC1*dk*grad))
    {
      output.logTab(1,VERBOSE,"Satisfied Wolfe condition after backtracking");
      output.logTabPrintf(1,VERBOSE,"Final LLG: %12.6f distance: %12.6f \n", -fk, dk);
      output.logTab(1,VERBOSE,"Bracketing took " + itos(mincount) + " function evaluations");
      if (dk < requiredShift) tooSmallShift = true;
      return fk;
    }
  }
  else // First step goes down.
  {
    if (dk >= maxDist) // First step already hit bounds
    {
      dhi = maxDist;
      fhi = fk;
    }
    else
    {
      dhi = std::min(maxDist,(ONE+GOLDEN)*dk);
      fhi = ShiftScore(dhi, target, x, oldx, gs, &output);
    }
    if (flo > fk && fk >= fhi) // Still not coming up at second step
    {
      while (fk >= fhi && dhi < maxDist)
      {
        flo = fk;
        dlo = dk;

        fk = fhi;
        dk = dhi;

        dhi = std::min(maxDist,dk + GOLDEN*(dk-dlo));
        fhi = ShiftScore(dhi, target, x, oldx, gs, &output);
      }
      if (dhi >= maxDist && fk >= fhi)
      {
        // Reached boundary without defining bracket.
        // Use finite differences to test if still going down.
        // If so, stop this line search.  Otherwise, we've verified bracket.
        dk = (1.-DTOL)*dhi;
        fk = ShiftScore(dk, target, x, oldx, gs, &output);
        if (fk >= fhi)
        {
          for (int i = 0; i < x.size(); i++) x[i] = oldx[i]-maxDist*gs[i];
          unrepar_x = target.reparRefineParsInv(x);
          hitBound = true;
          return fhi;
        }
      }
    }
    if ((dk > WOLFEFRAC && fk < f+WOLFEC1*dk*grad))
    {
      output.logTab(1,VERBOSE,"Satisfied Wolfe condition for bigger than initial step");
      output.logTabPrintf(1,VERBOSE,"Final LLG: %12.6f distance: %12.6f \n", -fk, dk);
      output.logTab(1,VERBOSE,"Bracketing took " + itos(mincount) + " function evaluations");

      ShiftX(dk,target,x,oldx,gs);  // Set to current best before returning
      if (dk < requiredShift) tooSmallShift = true;
      return fk;
    }
  }

  //===================================  End of bracketing

  output.logTab(1,VERBOSE,"Bracketing took " + itos(mincount) + " function evaluations");
  output.logTab(1,VERBOSE,"xlow= " + dtos(dlo) + ", xhigh= " + dtos(dhi));
  floatType dmid(dk);
  floatType fmid(fk);

  // BISECTING/INTERPOLATION STARTS HERE (we now know that the value is between dlo and dhi)

  // let tolerance be small but not too close to zero as to ensure our algorithm tries a new value
  // distinctly different from existing one
  const floatType XTOL(0.01);

  floatType a, b, c;
  floatType stepsize = std::min(dhi-dmid,dmid-dlo);
  floatType lastlaststepsize, laststepsize;
  laststepsize = stepsize*TWO; // permit first step to actually happen
  floatType d1,f1,d2,f2;

  output.logTabPrintf(1,VERBOSE,"Line Search |");
  for (int i=0; i<20; i++) // fall-back upper limit on evaluations in interpolation stage
  {
    if (GetQuadraticCoefficients(flo, fmid, fhi, dlo, dmid, dhi, a, b, c)) // Make sure there will be quadratic with a minimum
    {
      output.logTab(1, DEBUG, "Quadratic interpolation used" );
      dk = -b/(TWO*c); // where the derivative of parabolic curve is zero
      dk = std::max(dk,dlo); // remain within high and low interval
      dk = std::min(dk,dhi); // (Paranoia: fmid < an end should be enough)
    }
    else
    {
      output.logTab(1, DEBUG, "No quadratic interpolation possible, taking interval midpoint..." );
      c = ZERO;
      dk = dmid;
    }

    lastlaststepsize = laststepsize;
    laststepsize = stepsize;
    stepsize = fabs(dk - dmid);

    output.logTabPrintf(0,VERBOSE," stepsize= %12.6f, laststepsize= %12.6f,  ",stepsize, laststepsize);
    if (c > ZERO  // i.e. the interval has a local minimum, not a local maximum so proceed.
        && dk-dlo > XTOL*(dhi-dlo) // New minimum is sufficiently far from previous points
        && dhi-dk > XTOL*(dhi-dlo) // to make quadratic convergence probable.
        && fabs(dk-dmid) > XTOL*(dhi-dlo)
        && stepsize <= HALF*lastlaststepsize // Steps are decreasing sufficiently in size
        )
    {
      fk = ShiftScore(dk, target, x, oldx, gs, &output);
      output.logTab(1,DEBUG,"Use quadratic step");
    }
    else // Can't find a quadratic solution. Maybe points are on a line or local maximum
    {    // Golden search instead
      dk = (dhi-dmid >= dmid-dlo) ? dmid + GOLDFRAC*(dhi-dmid) : dmid - GOLDFRAC*(dmid-dlo);
      fk = ShiftScore(dk, target, x, oldx, gs, &output);
      output.logTab(1,DEBUG,"Use golden search step");
    }
    if (dk < dmid) // Get data sorted by size of shift
    {
      d1 = dk;
      f1 = fk;
      d2 = dmid;
      f2 = fmid;
    }
    else
    {
      d1 = dmid;
      f1 = fmid;
      d2 = dk;
      f2 = fk;
    }
    if (fk < fmid) // Keep current best at dmid
    {
      dmid = dk;
      fmid = fk;
    }

    output.logTabPrintf(0,VERBOSE,"=");

    output.logTabPrintf(0,DEBUG,"[%12.6f,%12.6f,%12.6f,%12.6f]\n",dlo,d1,d2,dhi);
    output.logTabPrintf(0,DEBUG,"[    ----    ,%12.6f,%12.6f,    ----    ]\n",f1,f2);
    if (f1 < f2)
      output.logTabPrintf(0,DEBUG,"          <        ^^^^        >                     \n");
    else
      output.logTabPrintf(0,DEBUG,"                       <       ^^^^       >          \n");

    // convergence tests
    if ((fmid < f - std::max(requiredGain,-WOLFEC1*dmid*grad)) && (dmid >= requiredShift))
    {
      output.logTab(1,DEBUG,"Stop linesearch: good enough for next cycle");
      break;
    }
    else
    { // There may not be another cycle, so try harder
      PHASER_ASSERT(starting_distance);
      floatType fracStart((d2-d1)/starting_distance);
      if (std::max(fabs(f2-f),fabs(f1-f))<requiredGain && fracStart<0.1)
      {
        output.logTab(1,DEBUG,"Stop linesearch: no sign of improvement");
        break;
      }
      if (fabs(f2-f1)<requiredGain && fracStart<0.1)
      {
        output.logTab(1,DEBUG,"Stop linesearch: |f1-f2| < tolerance");
        break;
      }
      if (fracStart<0.01)
      {
        output.logTab(1,DEBUG,"Stop linesearch: |d2-d1|/startdist < 0.01");
        break;
      }
    }
    // Not converged, so prepare for next loop
    if (f1 < f2)
    {
      dmid = d1;
      fmid = f1;
      dhi = d2;
      fhi = f2;
    }
    else
    {
      dlo = d1;
      flo = f1;
      dmid = d2;
      fmid = f2;
    }
  }
  // converged or limit in steps

  if (f < fmid) // Shouldn't happen, but catch possibility that function didn't improve
  {
    unrepar_x = target.reparRefineParsInv(oldx);
    target.applyShift(unrepar_x);
  }
  else
  {
    ShiftX(dmid,target,x,oldx,gs);  // Make sure best shift has been applied.
    f = fmid;
  }

  output.logBlank(VERBOSE);
  output.logTab(1,VERBOSE,"This line search took " + itos(mincount) + " function evaluations");
  output.logBlank(VERBOSE);
  output.logTabPrintf(1,VERBOSE,"Final LLG: %12.6f distance: %12.6f \n", -f, dmid);
  output.logTab(1,DEBUG,"End distance = " +dtos(dmid));
  output.logTab(1,DEBUG,"Prediction ratio = " +dtos(dmid/starting_distance));
  output.logBlank(DEBUG);
  if (dmid < requiredShift) tooSmallShift = true;
  return f;
}


// Find minimum from parabolic interpolation. Return true as long as minimum is between xmin and xmax
bool Minimizer::GetQuadraticMinimum(floatType &xmin, floatType x1, floatType y1, floatType x2,
                             floatType y2, floatType x3, floatType y3, floatType xlow, floatType xhigh)
{
  floatType a,b,c;
  GetQuadraticCoefficients(y1, y2, y3, x1, x2, x3, a, b, c);
  xmin = -b/(2.0*c); // where the derivative of parabolic curve is zero

  if (xmin <xlow)
  {
    xmin = xlow;
    return false;
  }

  if (xmin >xhigh)
  {
    xmin = xhigh;
    return false;
  }

  return true;
}

// Given three points in the x-y plane get the coefficients corrsponding to
// the parabolic equation y(x) = a + b*x + c*x^2
bool Minimizer::GetQuadraticCoefficients(floatType y1, floatType y2, floatType y3, floatType x1,
                                    floatType x2, floatType x3, floatType& a, floatType& b, floatType& c)
{
  a = b = c = 0.0;

  floatType det = x2*x3*x3 - x3*x2*x2 - x1*(x3*x3 - x2*x2) + x1*x1*(x3 - x2);
  if (fabs(det) < std::numeric_limits<floatType>::epsilon()*3) // too close to zero within machine precision
    return false;

  floatType invdet = 1.0/det;

  //   a = (y1*(x2*x3*x3 - x2*x2*x3) + y2*(x1*x3*x3 - x3*x1*x1) + y3*(x1*x2*x2 - x2*x1*x1))*invdet; // correct but unstable
  b = (y1*(x2*x2 - x3*x3) + y2*(x3*x3 - x1*x1) + y3*(x1*x1 - x2*x2))*invdet;
  c = (y1*(x3 -x2) + y2*(x1 - x3) + y3*(x2 -x1))*invdet;
  a = y1 - (b*x1 +c*x1*x1);

  return true;
}


}//phaser
