//  refinereferencescale.cpp
//
//
//  Refine scale of reference intensities to test set, using BFGS minimiser
//
//  Minimise f(d) = f(w(Iobs - k Iref)) where k is a ReferenceScaleModel
//    d = w(Iobs - k Iref)
//    weight w is usually "unit" weights, actually = 2/Meanintensity
//           as a crude normalisation, but could be 1/sigma(Iobs)
//           sigma weighting seems to give biased results with <kIref> < <Iobs>
//
//  Function to minimise is either:
//    1) least squares, ie R = Sum [(wd)^2]  if targettype == RefineTargets::QUADRATIC,
//       quadratic == true, or
//    2) R = Sum(ln(cosh(d)) as recommended by Garib Murshudov as more robust
//       targettype = RefineTargets::LNCOSH, quadratic = false
//
//  First (gradient) and second (Hessian) derivatives:
//  1) Quadratic
//     dR/dk    = Sum[ -2 w^2 Iref d ]
//     d2R/dk2  = Sum[ 2 w^2 Iref^2]
//  2) ln cosh
//     dR/dk    = Sum[ -w Iref tanh(wd) ]
//     d2R/dk2  = Sum[ (w Iref)^2/cosh^2(wd) ] (strictly) ...
//            ... but Garib says this is better replaced by
//     d2R/dk2  = Sum[ (w Iref)^2 (1/wd) tanh(wd) ]
//     Note that if |wd| > ~3, ln(cosh(wd)) ~= |wd|, tanh(wd) ~= +-1
//               if |wd| very small (or 0), (1/wd) tanh(wd) = +1
//

#include <assert.h>
#define ASSERT assert

#include "refinereferencescale.hh"
#include "referencescalemodel.hh"
#include "timer.hh"

namespace scala {
  // ---------------------------------------------------------
  RefineReferenceScale::RefineReferenceScale
  (clipper::HKL_data<clipper::data32::I_sigI>& Isigi,
   const hkl_merge& Hklmergelist,
   ReferenceScaleModel& Referencescalemodel,
   const double& Meanintensity,
   const int& Nprocs)
  {
    // reflection data
    isigi = &Isigi;
    hklmergelist = &Hklmergelist;
    referencescalemodel = &Referencescalemodel; // pointer to model
    npar = referencescalemodel->Nparameters();
    params = referencescalemodel->GetParameters();  // initial parameters
    gradient.newsize(npar);
    gradientOK = false;
    meanintensity = Meanintensity;
    nprocs=Nprocs;
    targettype = RefineTargets::LNCOSH;  // default ln(cosh()) target
  }
  // ---------------------------------------------------------
  // set true to use quadratic residual, else ln(cosh(d))
  void RefineReferenceScale::setQuadratic(const bool& quadratic)
  {
    if (quadratic) {
      targettype = RefineTargets::QUADRATIC;
    } else {
      targettype = RefineTargets::LNCOSH;
    }
  }
  // ---------------------------------------------------------
  void RefineReferenceScale::setTargetType
  (const RefineTargets::REFINETARGETTYPES& targetType)
  {
    targettype = targetType;
  }
  // ---------------------------------------------------------
  // ---------------------------------------------------------
  floatType RefineReferenceScale::targetFn()
  {
    if (!gradientOK)  {
      TNT::Fortran_Matrix<floatType> H;
      TargetGradientHessian(false, false, H);
    }
    return target;
  }
  // ---------------------------------------------------------
  floatType    RefineReferenceScale::gradientFn(TNT::Vector<floatType>& grad)
  {
    if (!gradientOK) {
      TNT::Fortran_Matrix<floatType> H;
      TargetGradientHessian(true, false, H);
    }
    grad = gradient;
    return target;
  }
  // ---------------------------------------------------------
  floatType RefineReferenceScale::hessianFn(TNT::Fortran_Matrix<floatType>& H,
                                   bool& is_diagonal)
  {
    is_diagonal = false;
    TargetGradientHessian(true, true, H);
    return target;
  }
  // ---------------------------------------------------------
  void RefineReferenceScale::TargetGradientHessian(bool DoGradient,
                                           bool DoHessian,
                                           TNT::Fortran_Matrix<floatType>& H)
  // Merged data version
  // Private function to calculate target function, gradient & Hessian
  //  gradient is stored locally, Hessian is returned
  // DoHessian implies DoGradient
  // floatType == double
  //
  {
    const bool DEBUG = false;
    //    const bool DEBUG = true;
    // maximum value of argument to cosh(x) before using approximation
    //   set in refinetargets.cpp
    //    const double MAXCOSHARG = 3.0;
    // if wd < MINCOSHARG, 2nd derivative = 1
    //    const double MINCOSHARG = 0.01;
    bool quadratictarget = false;
    if (targettype == RefineTargets::QUADRATIC) {quadratictarget = true;}

    // it is slightly faster to accumulate Hessian in one-dimensional array Hv and
    // then copy it into the 2D array H
    // Accumulating the Hessian is the rate-limiting step
    std::vector<floatType> Hv(npar*npar*nprocs,0.0);

    //  Timer timer;
    target = 0.0;
    if (DoHessian) {
      // Set up & clear Hessian
      DoGradient = true;
      H.newsize(npar,npar);
      for (int i=0;i<npar;i++) { // loop parameters
        for (int j=0;j<npar;j++) // loop parameters
          {H(i+1,j+1) = 0.0;}
      }
    }
    if (DoGradient) {
      // Clear gradient
      for (int i=0;i<npar;i++) // loop parameters
        {gradient[i] = 0.0;}
    }

    // loop all test reflections
    nref = 0;
    int nbig = 0;
    // Set weight type
    std::vector<double> dkdp;

    typedef clipper::HKL_data_base::HKL_reference_index HRI;

    double unitweight = 2.0/meanintensity;

    MeanValue meanIref, meanIobs, meanDi, meanDi2; //^

    for ( HRI ih = isigi->first(); !ih.last(); ih.next() ) {
      clipper::HKL hkl = ih.hkl();
      // find equivalent in reference list, if present
      IsigI Isref = hklmergelist->Isig(hkl);
      if (Isref.sigI() > 0.0) {
        ftype I = (*isigi)[ih].I();     // Iobs
        ftype sd = (*isigi)[ih].sigI();
        if (sd > 0.0) {
          double scale =
            referencescalemodel->fderiv(true, hkl, dkdp);
          //double w = 1.0/sd;
          double w = unitweight;  // equal weights
          // target function to minimise is Sum(h){w *(Iav - k Iref)^2}
          double di = I - scale * Isref.I();
          double wdi = std::abs(w*di);
          if (quadratictarget) {
            meanDi.Add(wdi);
            w *= w;  // square
            target += w * di * di;
          } else {
            // ln cosh
            meanDi.Add(wdi);
            meanDi2.Add(wdi*wdi);
            if (wdi > RefineTargets::MAXCOSHARG) {
              target += wdi;
              nbig++;
            } else {
              target += log(cosh(wdi));
            }
          }
          nref++;

          meanIref.Add(scale * Isref.I());
          meanIobs.Add(I);   //^

          if (DoGradient) {
            // dR/dp =  for parameter p
            double dlncdx = 1.0;
            if (!quadratictarget) {
              if (wdi < RefineTargets::MAXCOSHARG) {
                dlncdx = tanh(w*di);  // signed
              } else if (di < 0.0) {
                dlncdx = -1.0;
              }
            }
            for (int ip=0;ip<npar;++ip) {
              if (quadratictarget) {
                gradient[ip] -= 2.0 * w * di * Isref.I() * dkdp[ip];
              } else {
                gradient[ip] -= dlncdx * Isref.I() * dkdp[ip];
              }
            }
            if (DoHessian) {
              double Iref2 = Isref.I() * Isref.I();
              if (quadratictarget) {
                for (int ip=0;ip<npar;++ip) {
                  int ip1 = ip*npar;
                  for (int jp=0;jp<=ip;jp++) { // loop parameters again (to ip)
                    Hv[(ip1+jp)] += 2.0 * w * Iref2 *
                      dkdp[ip]*dkdp[jp];
                  }
                }
              } else { // ln cosh target
                double d2dx = 1.0;
                // NB test against MAXCOSHARG already done for dlncdx
                if (std::abs(wdi) > RefineTargets::MINCOSHARG) {
                  d2dx = dlncdx / (w*di);  // (1/d) tanh(d)
                  // exact 2nd derivative
                  //d2dx = 1.0/(cosh(wdi)*cosh(wdi));
                }
                for (int ip=0;ip<npar;++ip) {
                  int ip1 = ip*npar;
                  for (int jp=0;jp<=ip;jp++) { // loop parameters again (to ip)
                    Hv[(ip1+jp)] += w * w * Iref2 *
                      dkdp[ip]*dkdp[jp] * d2dx;
                  }
                }
              }
            } // end loop parameters
          }  // DoGradient
        }  // sdI > 0
      }
    } // end loop reflections

    if(DoHessian) {
      // Symmetrise Hessian
      for (int i=0;i<npar;i++) {
        for (int j=0;j<=i;j++) {
          H(j+1,i+1) += Hv[(i*npar+j)];
          H(i+1,j+1) = H(j+1,i+1); // other half
        }
      }
      if (DEBUG) {
        printf("Target: %10.1f\n", target);
        printf("Gradient:\n");
        for (int i=0;i<npar;i++){
          printf(" %8.1f", gradient[i]);
        }
        printf("\n");

        printf("Full Hessian:\n");
        for (int j=0;j<npar;j++){
          for (int i=0;i<npar;i++){
            printf(" %.3lf", H(i+1,j+1));
          }
          printf("\n");
        }
        if (npar == 2) {
          double det = H(1,1) * H(2,2) - H(1,2) * H(2,1);
          printf("Determinant %.4f\n", det);
        }
      }

    } else { // Hessian
      if (DEBUG) {
        printf("Target: %10.1f\n", target);
      }
    }
    if (DEBUG) {
      printf("<Iref> <Iobs> <|wDi|> rms(wDi), nbig %8.1f %8.1f %8.2f %8.2f %7d\n",
             meanIref.Mean(), meanIobs.Mean(),
             meanDi.Mean(), sqrt(meanDi2.Mean()), nbig);
    }

    if (DoGradient) gradientOK = true;
  }
  // ---------------------------------------------------------
  void RefineReferenceScale::applyShift(TNT::Vector<floatType>& newpar)
  {
    //    std::cout << "applyShift: pars"; //^
    for (int i=0;i<npar;i++)  {
      params[i] = newpar[i];
      //      std::cout << " " << params[i]; //^
    }
    //    std::cout <<"\n"; //^
    referencescalemodel->SetParameters(params);
    gradientOK = false;
  }
  // ---------------------------------------------------------
  void RefineReferenceScale::logCurrent(outStream where, Output& output)
  {
    /*
      output.logTab(1,where,"\nScales:");
      for (int i=0;i<npar;i++)
      {output.logTabPrintf(2,where," %7.3f", params[i]);}
      output.logTab(1,where,"\n");
    */
  }
  // ---------------------------------------------------------
  std::string  RefineReferenceScale::whatAmI(int& i)
  {
    return "Scale "+itos(i);
  }
  // ---------------------------------------------------------
  std::vector<bounds>     RefineReferenceScale::getLowerBounds()
  {
    std::vector<bounds> Lower(npar);
    double bound;
    for (int i=0;i<npar;i++) {
      if (referencescalemodel->GetLowerBound(i, bound)) {
        Lower[i].on(bound);
      } else {
        Lower[i].off();
      }
    }
    return Lower;
  }
  // ---------------------------------------------------------
  std::vector<bounds>     RefineReferenceScale::getUpperBounds()
  {
    std::vector<bounds > Upper(npar);
    double bound;
    for (int i=0;i<npar;i++) {
      if (referencescalemodel->GetUpperBound(i, bound)) {
        Upper[i].on(bound);
      } else {
        Upper[i].off();
      }
    }
    return Upper;
  }
  // ---------------------------------------------------------
  TNT::Vector<floatType> RefineReferenceScale::getRefinePars()
  {
    TNT::Vector<floatType> pars(npar);
    for (int i=0;i<npar;i++) {
      pars[i] = params[i];
    }
    return pars;
  }
  // ---------------------------------------------------------
  TNT::Vector<floatType> RefineReferenceScale::getLargeShifts()
  {
    TNT::Vector<floatType> large(npar);
    for (int i=0;i<npar;i++) {
      ///      large[i] = 1.0;
      large[i] = referencescalemodel->GetLargeShift(i);
    }
    return large;
  }
  // ---------------------------------------------------------
  bool1D RefineReferenceScale::getRefineMask(protocolPtr protocol)
  {
    return bool1D(npar, true);
  }
  // ---------------------------------------------------------
}
