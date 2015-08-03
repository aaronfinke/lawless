// simpleminimise.cpp

// A simple Gauss-Newton minimiser, using Clipper classes

#include "simpleminimise.hh"

#include <assert.h>
#define ASSERT assert

namespace SimpleMinimise {
  // ---------------------------------------------------------
  DampedGaussNewton::DampedGaussNewton(FitBase& fitstuff,
                                       const int& Ncycles,
                                       const double& tolerance,
                                       const double& damp)
  {
    valid_ = false;
    //    std::cout << "DampedGaussNewton " << fitstuff.NvalidData()<<
    //      " "<< fitstuff.Nparameters() << "\n"; //^
    if (fitstuff.NvalidData() >= fitstuff.Nparameters()+2) {
      run(fitstuff, Ncycles, tolerance, damp);
    }
  }
  // ---------------------------------------------------------
  void dumpHmatrix(const clipper::Matrix<double>& H)
  // dump Hessian (DEBUG)
  {
    int npar = H.rows();
    std::cout << "\n";
    for (int i=0;i<npar;++i) {
      for (int j=0;j<npar;++j) {
        std::cout << " " << H(i,j);
      }
      std::cout << "\n";
    }
  }
  // ---------------------------------------------------------
  void DampedGaussNewton::run(FitBase& fitstuff,
                              const int& Ncycles,
                              const double& tolerance,
                              const double& damp)
  // On entry, FitBase contains initial parameters, updated on return
  // as well as the data
  {
    int npar = fitstuff.Nparameters();
    ASSERT (npar < 20);  // Clipper can't cope with large matrices

    int ndata = fitstuff.Ndata();
    std::vector<double> params0 = fitstuff.parameters();

    double target;
    double previoustarget;
    std::vector<double> gradient;
    clipper::Matrix<double> H;    // Hessian

    bool trying = true;
    double damping = damp;
    int kcycles = 0;  // big cycles
    const int MAXCYCLES = 3;
    const double DAMPSTEP = 0.1;
    bool OK;

    while (trying) {
      for (int icyc=0;icyc<Ncycles;++icyc) { // loop cycles
        TGH tgh = fitstuff.TargetGradientHessian();
        target = tgh.target;
        gradient = tgh.gradient;
        H = tgh.H;

        // Scale Hessian to make diagonals = 1.0
        std::vector<double> U(npar);    // diagonal of scaling matrix
        std::vector<double> Uinv(npar); // diagonal of inverse scaling matrix
        //      dumpHmatrix(H);
        for (int i=0;i<npar;++i) {
          ASSERT (H(i,i) > 0.0);  // positive definite
          U[i] = sqrt(H(i,i));
          Uinv[i] = 1.0/U[i];
        }

        clipper::Matrix<double> A(npar,npar);  // scaled Hessian
        // Scaled gradient vector
        std::vector<double> ugradient(npar, 0.0);
        for (int j=0;j<npar;++j) {
          ugradient[j] = - gradient[j] * Uinv[j]; // scaled gradient
          for (int i=0;i<npar;++i) {
            A(i,j) = H(i,j) * Uinv[i] * Uinv[j];
          }
          if (damping > 0.0) {
            A(j,j) += damping;
          }
        }

        // Get A^-1 via eigenvectors
        std::vector<double> ev = A.eigen();
        // Diagonal matrix of inverse eigenvalues
        //   L^-1(ii) = 1/(l(i) + damping)
        //   or if eigenvalue l(i) = ev[i] <= 0, L^-1ii = 0.0
        std::vector<double> Linv(npar, 0.0);
        for (int i=0;i<npar;++i) {
          if (ev[i] > 0.0) {
            Linv[i] = 1.0/(ev[i] + damping);
          }
        }

        clipper::Matrix<double> Ainv(npar,npar);
        clipper::Matrix<double> Hinv(npar,npar);
        // Matrix product:
        // [AB]ij = Sum(k) [A]ik [B]kj
        // [[A] [A]T ]ij = Sum(k) [A]ik [A]jk
        // [A]^-1 = [E] [L]^-1 [E]T  where [L] is diagonal eigenvalues
        //       [E] is matrix of eigenvectors, transpose [E]T = [E]^-1
        // A now contains matrix of eigenvectors
        for (int j=0;j<npar;++j) {
          for (int i=0;i<npar;++i) {
            Ainv(i,j) = 0.0;
            for (int k=0;k<npar;++k) {
              // [E] = A
              // [[E][L^-1]]ik = [E]ik (1/Lk)
              // [[E][L^-1][E]T]ij = Sum(k) [[E][L^-1]]ik [E]jk
              Ainv(i,j) += Linv[k] * A(i,k) * A(j,k);
            }
            // [H]^-1 = U^-1 A^-1 UT^-1
            Hinv(i,j) = Ainv(i,j) * Uinv[i] * Uinv[j];
          }
        }
        // scaled shifts
        std::vector<double> shifts = Ainv * ugradient;
        for (int i=0;i<npar;++i) {
          shifts[i] *= Uinv[i];  // unscale
        }

        // crudely finds best shift, updates parameters and target
        OK = linesearch(fitstuff, shifts, target);
        if (!OK) {
          break;
        }

        // Variance/covariance matrix = (R/(m-n)) H^-1
        //   for m observations, n variables
        double rmn = target/(double(ndata-npar));  // scaling factor
        std::vector<double> sdpkpl(npar);
        double bigrelshift = -1.0;

        for (int kpl=0;kpl<npar;++kpl) {  // Loop parameters kpl
          sdpkpl[kpl] = sqrt(rmn * Hinv(kpl,kpl));
          double relshift = shifts[kpl]/sdpkpl[kpl];
          bigrelshift = std::max(bigrelshift, std::abs(relshift));
          //^
          //      std::cout <<"Parameter "<<kpl<<" Shift "<<shifts[kpl]
          //                <<" sdpk " <<sdpkpl[kpl]<<" relshift " <<relshift
          //                <<" rmn " << rmn<<"\n"; //^-
        }

        if (bigrelshift < tolerance) {
          break;
        }
      } // end cycle
      if (OK) {break;} // OK so exit
      // failed, so reset, increase damping, try again
      kcycles++;
      if (kcycles > MAXCYCLES) {
        // fail
        valid_ = false;
        return;
      }
      damping += DAMPSTEP;
      fitstuff.SetParameters(params0); // reset parameters
      //^
      //      std::cout << damping << " new damp, Not OK\n";
      //^-
    } // end trying
    valid_ = true;
    return;
  }
  // ---------------------------------------------------------
  bool DampedGaussNewton::linesearch(FitBase& fitstuff,
                                     const std::vector<double>& shifts,
                                     double& target) const
  // returns false if can't find anything better than start, updates target
  // a VERY crude linesearch
  {
    // save initial parameters
    std::vector<double> params0 = fitstuff.parameters();
    // initial target value
    double target0 = target;

    // search fraction step
    const double FRACTSTEP = 0.2;

    int icyc = 0;
    double fract = 0.0;
    double bestfract = fract;
    double besttarget = target0;  // at lowfract

    while (fract <= 1.0001) {
      std::vector<double> trialparams = shift(params0, fract, shifts);
      fitstuff.SetParameters(trialparams);
      double newtarget = fitstuff.Target();
      //^
      //      std::cout << "DampedGaussNewton::linesearch fract, target "<<
      //        fract<<" "<<newtarget<<"\n";

      // is this better than what we have?
      if (newtarget < besttarget) {
        besttarget = newtarget;
        bestfract = fract;
      }
      fract += FRACTSTEP;
    }
    // reset to "best" values
    std::vector<double> trialparams = shift(params0, bestfract, shifts);
    fitstuff.SetParameters(trialparams);
    target = besttarget;

    //^^
    //    std::cout << "DampedGaussNewton::linesearch fract "<<
    //      bestfract <<"\n";
    //^-
    bool OK = false;
    if (bestfract > 0.0) {
      OK = true;
    } else {
      fitstuff.SetParameters(params0);
      target = target0;
    }
    return OK;
  }
  // ---------------------------------------------------------
  std::vector<double>  DampedGaussNewton::shift
  (const std::vector<double>& params0, const double& fract,
   const std::vector<double>& shifts) const
  {
    std::vector<double> params(params0.size());
    for (size_t i=0;i<shifts.size();++i) {
      params[i] = params0[i] + fract * shifts[i];
    }
    return params;
  }


}
