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
    run(fitstuff, Ncycles, tolerance, damp);
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

    double target;
    std::vector<double> gradient;
    clipper::Matrix<double> H;    // Hessian
    int kcycles = 0;

    for (int icyc=0;icyc<Ncycles;++icyc) { // loop cycles
      TGH tgh = fitstuff.TargetGradientHessian();
      target = tgh.target;
      gradient = tgh.gradient;
      H = tgh.H;

      // Scale Hessian to make diagonals = 1.0
      std::vector<double> U(npar);    // diagonal of scaling matrix
      std::vector<double> Uinv(npar); // diagonal of inverse scaling matrix
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
        if (damp > 0.0) {
          A(j,j) += damp;
        }
      }

      // Get A^-1 via eigenvectors
      std::vector<double> ev = A.eigen();
      // Diagonal matrix of inverse eigenvalues
      //   L^-1(ii) = 1/(l(i) + damp)
      //   or if eigenvalue l(i) = ev[i] <= 0, L^-1ii = 0.0
      std::vector<double> Linv(npar, 0.0);
      for (int i=0;i<npar;++i) {
        if (ev[i] > 0.0) {
          Linv[i] = 1.0/(ev[i] + damp);
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

      //Apply shifts
      fitstuff.ApplyShifts(shifts);

      // Variance/covariance matrix = (R/(m-n)) H^-1
      //   for m observations, n variables
      double rmn = target/(double(ndata-npar));  // scaling factor
      std::vector<double> sdpkpl(npar);
      double bigrelshift = -1.0;

      //      std::cout <<"Cycle "<< icyc <<", target " << target <<std::endl;

      for (int kpl=0;kpl<npar;++kpl) {  // Loop parameters kpl
        sdpkpl[kpl] = sqrt(rmn * Hinv(kpl,kpl));
        double relshift = shifts[kpl]/sdpkpl[kpl];
        bigrelshift = std::max(bigrelshift, std::abs(relshift));
        //^
        //      std::cout <<"Parameter "<<kpl<<" Shift "<<shifts[kpl]
        //                <<" sdpk " <<sdpkpl[kpl]<<" relshift " <<relshift
        //                        <<" rmn " << rmn<<"\n"; //^-
      }

      kcycles = icyc;
      if (bigrelshift < tolerance) {
        break;
      }
    } // end cycle
    //    std::cout << "Exit after "<<kcycles <<" cycles\n";
  }
// ---------------------------------------------------------
}
