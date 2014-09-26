// fitresolutiondata.cpp
//

#include "fitresolutiondata.hh"

namespace scala {
  // ---------------------------------------------------------
  FitResolutionData::FitResolutionData(const std::vector<ResolutionData>& Data,
                                       RadialBase& Radialfunction)
  {
    rdata = &Data;
    radialfunction = &Radialfunction;
    npar = radialfunction->Nparameters();
    params = radialfunction->parameters();  // initial parameters
    //targettype = RefineTargets::LNCOSH;  // default ln(cosh()) target
    targettype = RefineTargets::QUADRATIC;  // default ln(cosh()) target
  }
  // ---------------------------------------------------------
  void FitResolutionData::ApplyShifts(const std::vector<double> shifts)
  {
    for (int i=0;i<npar;++i) {
      params[i] += shifts[i];
    }
    radialfunction->SetParameters(params);
  }
  // ---------------------------------------------------------
  TGH FitResolutionData::TargetGradientHessian()
  // Function to calculate target function, gradient & Hessian
  // returns target, gradient, Hessian
  {
    //    const bool DEBUG = true;
    const bool DEBUG = false;

    target = 0.0;
    bool quadratictarget = false;
    if (targettype == RefineTargets::QUADRATIC) {quadratictarget = true;}
    std::vector<double> Hv;

    // Set up & clear Hessian
    clipper::Matrix<double> H(npar,npar);
    for (int i=0;i<npar;i++) { // loop parameters
      for (int j=0;j<npar;j++) // loop parameters
        {H(i,j) = 0.0;}
    }
    // it is slightly faster to accumulate Hessian in
    // one-dimensional array Hv and then copy it into the 2D array H
    // Accumulating the Hessian is the rate-limiting step
    //   for large number of parameters (though probably not relevant here)
    Hv.assign(npar*npar,0.0);
    std::vector<double> gradient(npar, 0.0);   // Clear gradient

    std::vector<double> dvcdp(npar);

    int nbig = 0;
    nobs = 0;

    for (size_t k=0; k<rdata->size(); k++) { // loop data
      if ((*rdata)[k].w > 0.0) {
        double s = (*rdata)[k].s;
        double v = (*rdata)[k].v;
        double w = (*rdata)[k].w;

        double vc = radialfunction->value(s); // calculated value
        double di = (v - vc);
        double wdi = w * di;  // weighted residual

        nobs++;

        if (quadratictarget) {
          w *= w;  // square
          target += w * di * di;
        } else {
          // ln cosh
          if (wdi > RefineTargets::MAXCOSHARG) {
            target += wdi;
            nbig++;
          } else {
            target += log(cosh(wdi));
          }
        }
        dvcdp = radialfunction->deriv(s);
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
            gradient[ip] -= 2.0 * wdi * dvcdp[ip];
          } else {
            gradient[ip] -= dlncdx * dvcdp[ip];
          }
        }

        if (quadratictarget) {
          for (int ip=0;ip<npar;++ip) {
            int ip1 = ip*npar;
            for (int jp=0;jp<=ip;jp++) { // loop parameters again (to ip)
              //   note w^2
              Hv[(ip1+jp)] += 2.0 * w * dvcdp[ip] * dvcdp[jp];
              //              if (DEBUG) {
              //                printf("ip %3d jp %3d dvcdp %8.4f %8.4f H %8.4f\n",
              //                       ip, jp, dvcdp[ip], dvcdp[jp], Hv[(ip1+jp)]);
              //            }
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
              Hv[(ip1+jp)] += w * w * dvcdp[ip] * dvcdp[jp] * d2dx;
              if (DEBUG) {
                printf("ip %3d jp %3d dvcdp %8.4f %8.4f H %8.4f\n",
                       ip, jp, dvcdp[ip], dvcdp[jp], Hv[(ip1+jp)]);
              }
            }
          }
        }
      } // null data
    } // end loop data

    // Symmetrise Hessian
    for (int i=0;i<npar;i++) {
      for (int j=0;j<=i;j++) {
        H(j,i) += Hv[(i*npar+j)];
        H(i,j) = H(j,i); // other half
      }
    }

    if (DEBUG) {
      printf("Nobs: %5d\n", nobs);
      printf("Target: %10.1f\n", target);
      printf("Gradient:\n");
      for (int i=0;i<npar;i++){
        printf(" %8.1f", gradient[i]);
      }
      printf("\n");

      printf("Full Hessian:\n");
      for (int j=0;j<npar;j++){
        for (int i=0;i<npar;i++){
          printf(" %.3lf", H(i,j));
        }
        printf("\n");
      }
      if (npar == 2) {
        double det = H(1,1) * H(2,2) - H(1,2) * H(2,1);
        printf("Determinant %.4f\n", det);
      }
      std::cout << "Parameters: " << radialfunction->format() << std::endl;
    }

    return TGH(target, gradient, H);

  } // TargetGradientHessian
// ---------------------------------------------------------
}
