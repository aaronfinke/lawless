//  RefineSDCorrection.cpp
//
// residual Sum(j) [ wj^2 (1 - sigma(delta(j)))^2 ] [QUADRATIC]
// or       Sum(j) [ ln cosh (wj (1 - sigma(delta(j)))) ] [LNCOSH]
// delta(hl) = (Ihl - <Ih>!l)/[sqrt(nh/nh-1) * sigma'(hl)]
// sigma'(hl) = SdFac *sqrt[sigma(hl)^2 + SdB <Ih> + (SdAdd * <Ih>)^2]
//            = p * sigma(hl)^2  + q * <Ih> + r * <Ih>^2
// p = SdFac^2 ; q = SdFac^2 SdB; r = SdFac^2 SdAdd^2
//
//  Note that parameters for each "class" (run+full/partial) can be optimised
//  separately, but sums over all reflections are done once/cycle for all groups.
//  For this reason, the BFGS minimiser cannot easily be used, so it's done
//  explicitly here, using clipper::Matrix with eigenvalue decomposition.
//  The matrix dimensions (nparam) are 2 or 3, so the clipper routines should be OK
//  (they are not suitable for larger matrices)
//
//

#include <assert.h>
#define ASSERT assert

#include "refinesdcorrection.hh"

namespace scala {
  // ---------------------------------------------------------
  // Refine SD correction model using LSQ minimiser for each "bin class" separately
  SDanalysis RefineSDcorrectionFactors(SDmodel& SDM,
                                       const hkl_unmerge_list& hkl_list,
                                       const all_controls& controls,
                                       IntensityBin& irange,
                                       const double& tolerance, const double& rtolerance,
                                       const int&  max_cycles,
                                       phaser_io::Output& output)
  {
    //^  std::cout << "SDM start: " << SDM.format() <<"\n"; //^

    //  int min_cycles = Min(max_cycles, 4); // at least 4 cycles unless < max
    int min_cycles = Min(max_cycles, 5); // at least N cycles unless < max
    double lastR = -1.0;
    double bestR = 10000.0;
    int bestcycle = -1;
    const int MAXNINCREASING = 3;  // allow 3 cycles with increasing R
    std::vector<double> bestsdmparams = SDM.GetParameters();

    output.logTab(0,LOGFILE,"\n");
    SDanalysis sdanal;

    double damp = SDM.Damp(); // damping factor for minimisation
    if (damp < 0.0) { // unset
      damp = 0.05;
    }
    output.logTabPrintf(0,LOGFILE,"Damping factor: %5.3f\n", damp);
    // print information about parameter restraints
    output.logTab(0,LOGFILE, SDM.formatTie());

    // Update targets and ties if similarity targets
    SDM.SetTargetsFromAverageParameters();

    int nmacrocycles = 0;
    const int MAXMACROCYCLES = 2;
    double R;
    std::vector<double> lastsdmparams;

    bool saveSdBfix = SDM.NoSDb(); // SdB fix value, true if fixed
    bool exit = true;  // to exit from macrocycles
    while (true) { // Macrocycles
      nmacrocycles++;
      if (nmacrocycles > MAXMACROCYCLES) {break;}
      int nincreasing = 0;
      lastR = -1.0;
      bestR = 10000.0;
      bestcycle = -1;

      for (int cyc=0;cyc<Max(1,max_cycles);++cyc) { // loop cycles
        // Accumulate all sums from data
        bool anomalous = controls.anomalouscontrol.AnomalousSDcorr;
        sdanal = SumsforSDcorrection(SDM, hkl_list, anomalous, irange);
        //^
        //    PrintSDanalysis(sdanal, SDanalysis(), RejectFlags(), irange, hkl_list.RunList(),
        //                  SDM, -1, PxdName(), false, output);
        //^-
        lastsdmparams = SDM.GetParameters();
        bool update = (max_cycles > 0); // don't update parameters if zero cycles
        TargetResiduals target = UpdateParameters(SDM, sdanal, tolerance, damp, update);
        R = target.R();
        if (SDM.restrained()) {  // restraints
          if (target.quadratic) {
            output.logTabPrintf(0,LOGFILE,
                                "Cycle %3d residual %10.5f   (main residual %8.5f restraint residual %8.5f)\n",
                                cyc+1, R, target.R1, target.R2);
          } else {
            output.logTabPrintf(0,LOGFILE,
                                "Cycle %3d residual %10.5f   (ln(cosh) residual %8.5f LSQresidual %8.5f restraint residual %8.5f)\n",
                                cyc+1, R, target.R1, target.R1lsq, target.R2);
          }
        } else { // no restraints
          if (target.quadratic) {
            output.logTabPrintf(0,LOGFILE,
                                "Cycle %3d residual %10.5f\n",
                                cyc+1, R);
          } else {
            output.logTabPrintf(0,LOGFILE,
                                "Cycle %3d ln(cosh) residual %10.5f LSQresidual %8.5f\n",
                                cyc+1, R, target.R1lsq);
          }
        }
        //^
        //      output.logTab(0,LOGFILE,
        //                    "\nSD correction parameters after cycle\n"+SDM.format()); //^-

        if (max_cycles <= 0) break;
        if (target.converged) {
          output.logTab(0,LOGFILE,"Convergence reached");
          break;
        }
        // Record the best so far
        if (R < bestR) {
          bestR = R;
          bestsdmparams = lastsdmparams;  // before last cycle
          bestcycle = cyc+1;
        }
        // Are the residuals increasing?
        if (lastR > 0.0 && R > lastR) {
          nincreasing++;
        } else {
          nincreasing = 0;
        }
        if (cyc+1 > min_cycles && lastR > 0.0) {
          // beyond minimum cycles, should we stop anyway?
          double diffR = std::abs(R - lastR);
          if (nincreasing > MAXNINCREASING) {
            // residual gone up, reinstate best parameter set
            SDM.SetParameters(bestsdmparams); // reset parameters
            if (nmacrocycles == 1) { // first macrocycle
              // If we are refining SdB, try again fixing it at the "best" values
              if (!saveSdBfix) {
                SDM.SetNoSDb(true);
                bestsdmparams = SDM.GetParameters();

                exit = false;  // to recycle
                output.logTabPrintf(0,LOGFILE,
                                    "Residual increasing, revert to best cycle %4d and try with SdB fixed\n",
                                    bestcycle);
                output.logTab(0,LOGFILE,
                              "\nSD correction parameters\n"+SDM.format()+"\n\n");
              }
            } else {
              output.logTabPrintf(0,LOGFILE,
                                  "Residual increasing, revert to best cycle %4d and exit\n",
                                  bestcycle);
            }
            break;
          }
          if (diffR < rtolerance) {
            // change in residual less than tolerance
            break;
          }
        }
        if (cyc < max_cycles-1) { // not last cycle
          lastR = R;
        }
      }  // loop cycles
      // Pick best cycle
      if (exit) {
        break;
      }
    } // loop macrocycles

    if (R > bestR) { //
      SDM.SetParameters(bestsdmparams); // reset parameters
      output.logTabPrintf(0,LOGFILE,
                          "Residual increased at end, revert to best cycle %4d and exit\n",
                          bestcycle);
    }
    SDM.SetNoSDb(saveSdBfix); // restore SdB fix value, true if fixed

    //^
    //  PrintSDanalysis(sdanal, SDanalysis(), RejectFlags(), irange, hkl_list.RunList(),
    //                      SDM, -1, PxdName(), false, output);
    //^-

    return sdanal;
  }
  // ---------------------------------------------------------
  SDanalysis SumsforSDcorrection(const SDmodel& SDM,
                                 const hkl_unmerge_list& hkl_list,
                                 const bool& anomalous,
                                 IntensityBin& irange)
  // Accumulate sums for SD correction refinement into SDanalysis object returned
  // anomalous   true to separate anomalous I+ & I- (usually true)
  {
    int Ndatasets = hkl_list.num_datasets();

    SDanalysis sdanal(irange, SDM, SDM.AllRunsSame(), true);
    // Set weight type
    WeightType::AverageWeightType weighttype = SDM.Weight();

    reflection this_refl;
    int nref = 0;
    hkl_list.rewind();
    std::vector<float> delta, delta2;

    //  deviations within each run & full/partial
    while (hkl_list.next_reflection(this_refl) >= 0)  {
      bool Centric = hkl_list.symmetry().is_centric(this_refl.hkl());
      // Correct sds in this_refl, return uncorrected scaled values
      std::vector<float> sig0 = SDM.CorrectReflection(this_refl);
      // Average I <I> over all observations
      SelectedObservations selobs(this_refl, -1, ALL, weighttype);
      //^
      //      std::cout << "rsd " << this_refl.hkl().format()
      //                << " N " << selobs.Number() << "\n";
      //^-
      float Iav = selobs.Average().I(); // average intensity for SD correction
      int mint = irange.bin(Iav);
      nref++;
      int nacc = 0; //^
      for (int id=0;id<Ndatasets;id++) {    // loop datasets
        if (Centric || !anomalous) {
          // No anomalous, treat all observations together
          if (Ndatasets > 1) {selobs.init(this_refl, id, ALL, weighttype);} // already done if 1 dataset
          if (selobs.Number() > 1) {
            sdanal.AddSelobsDelta2(selobs, mint);
            // partial derivatives, sums for d(delta)/dp
            sdanal.AddDerivatives(selobs, mint, SDM.GetDerivatives(selobs, sig0));
            nacc++;
          }
        } else {
          // Anomalous, treat I+ & I- separately
          selobs.init(this_refl, id, IPLUS, weighttype);
          if (selobs.Number() > 1) {
            sdanal.AddSelobsDelta2(selobs, mint);
            // partial derivatives, sums for d(delta)/dp
            sdanal.AddDerivatives(selobs, mint, SDM.GetDerivatives(selobs, sig0));
            nacc++;
          }
          selobs.init(this_refl, id, IMINUS, weighttype);
          if (selobs.Number() > 1) {
            sdanal.AddSelobsDelta2(selobs, mint);
            // partial derivatives, sums for d(delta)/dp
            sdanal.AddDerivatives(selobs, mint, SDM.GetDerivatives(selobs, sig0));
            nacc++;
          }
        } // end acentric
      } // end loop datasets
    } // end loop reflections

    return sdanal;
  }
  // ---------------------------------------------------------
  TargetResiduals UpdateParameters(SDmodel& SDM, const SDanalysis& sdanal,
                                   const double& tolerance, const double& damp,
                                   const bool& Update)
  // Solve LSQ equations for each "parameter class" to update SDmodel SDM,
  //    using sums in sdanal
  // Parameter classes are run and full/partial, subject to the Allrunssame flag and
  // few/no fulls/partials selections (useflags in SDM)
  //
  // Each parameter class is treated independently: this is possible since the
  //  weighted deviations delta within each class are independent, as they are calculated
  //  relative to <Ih>av weighted by 1/sqrt(scale), which is independent of sd(I)
  //  (at least for 1/sqrtscale weighting)
  //
  // Parameters may be tied to target values (set in SDmodel)
  // If Update false, just calculate residuals, no parameter update
  {
    //    RefineTargets::REFINETARGETTYPES targettype = RefineTargets::LNCOSH;
    RefineTargets::REFINETARGETTYPES targettype = RefineTargets::QUADRATIC;
    bool quadratictarget = false;
    if (targettype == RefineTargets::QUADRATIC) {quadratictarget = true;}

    // Number of parameter groups in analysis
    int npargroups = sdanal.NumberOfParameterGroups();

    // SD(delta(jc)) for each bin class jc (over all parameter classes)
    std::vector<double> sddelta = sdanal.SDdelta();
    std::vector<int> ninclass = sdanal.NumberinClass();  // number in each class

    //  dsigDeldp[jc][inb][k] partial derivatives d(sigma(delta(jc)))/dp(k) for
    //       parameter class jc, intensity bin inb, parameter k, k is index local to class jc
    std::vector<std::vector<std::vector<double> > > dsigDeldp = sdanal.Derivatives();

    int nintbins = sdanal.NumberIntensityBins();

    std::vector<double> sdmparams = SDM.GetParameters();

    // Restraint R2 for each parameter group
    std::vector<double> sdmrestraintR = SDM.GetRestraintR();
    ASSERT (int(sdmrestraintR.size()) == npargroups);
    // gradients and Hessians for each parameter group for restraints
    std::vector <std::vector<double> > dr2dp;
    std::vector <clipper::Array2d<double> > H2;
    SDM.GetRestraintDerivatives(dr2dp, H2);

    TargetResiduals target;    // R1, R2, R1lsq
    target.quadratic = false;
    if (quadratictarget) {
      target.quadratic = true;
    }
    double bigrelshift = -1.0;
    // relative weight for main residual compared to restraint residual
    const double WTREL = 0.1;
    // "SD" of residual in each intensity bin, for weighting by 1/SD^2
    const double SDRESID = 0.04;

    // for each parameter, true if values updated
    std::vector<bool> parameterupdated(sdmparams.size(), false);

    for (int jpc=0;jpc<npargroups;++jpc) { // loop parameter classes
      // weight for each intensity bin, equal (unit) weights
      double w1 = sqrt(WTREL)/(SDRESID*nintbins);
      //double w1 = sqrt(WTREL)/SDRESID;
      std::vector<double> wib(nintbins, w1);
      // Test! double weight on top bin
      wib.back() *= sqrt(2);

      double R1 = 0.0; // main residual
      double R1lsq = 0.0; // main LSQ residual
      double sumw = 0.0;
      int npar = sdanal.Nparam(jpc); // number of parameters for this parameter class
      ASSERT (int(dr2dp[jpc].size()) == npar);

      // index to 1st parameter in global list for class jpc
      int idxpar = sdanal.IdxParam(jpc);
      int number = 0;

      for (int mint=0;mint<nintbins;++mint) { // loop intensity bins
        int jc = sdanal.BinClass(jpc, mint);
        if (sddelta[jc] != 0.0) {
          number += ninclass[jc];
          double r = 1.0 - sddelta[jc]; // deviation
          // wib is sqrt of LSQ weight, ie weight on Delta
          double abswr = std::abs(wib[mint] * r);
          sumw += wib[mint]*wib[mint];
          R1lsq += abswr * abswr;;       // target residual ( * 2)
          if (quadratictarget) { // LSQ
            R1 += R1lsq;;       // target residual ( * 2)
          } else { // ln cosh
            if (abswr > RefineTargets::MAXCOSHARG) {
              R1 += abswr;
            } else {
              R1 += log(cosh(abswr));
            }
          }
        }
      }
      if (sumw > 0.0) {
        R1 = 0.5 * R1/sumw;
        R1lsq = 0.5 * R1lsq/sumw;
      }
      else {
        R1 = 0.0;
        R1lsq = 0.0;
      }  // target residual

      //^
      //      std::cout <<"Target contribution from group " << jpc<< " " << R1
      //                << " number " << number <<"\n"; //^-

      if (number > 0) {
        std::vector<double> gradient(npar, 0.0); // gradient vector dR/dp
        clipper::Matrix<double> H(npar,npar,0.0);  // Hessian ~= d2R/dp2

        // kpl is parameter number local to parameter class
        // kpg is global parameter number
        for (int kpl=0;kpl<npar;++kpl) {  // Loop parameters kpl for this class
          // dR/dpk = Sum[jc] (w (1 - SD(delta)) d(sigma(delta(jc)))/dpk)
          for (int mint=0;mint<nintbins;++mint) { // loop intensity bins
            int jc = sdanal.BinClass(jpc, mint);  // bin class number from param class & intbin
            if (sddelta[jc] != 0.0) {
              double r = 1.0 - sddelta[jc]; // deviation
              if (quadratictarget) { // LSQ
                double w2 = wib[mint] * wib[mint];
                gradient[kpl] += w2 * r * dsigDeldp[jpc][mint][kpl]; //??
                //^
                //              std::cout << "Gradient k, jc " << kpl << " " << jc
                //                        << " r " << r << " w " << wib[mint]
                //                        << " sddelta " << sddelta[jc]
                //                        << " Number " << sdanal.NumberinClass()[jc]
                //                        << " d(sig(delta))/dp " << dsigDeldp[jpc][mint][kpl]
                //                        << " g(k)(jc) " << - wib[mint] * r * dsigDeldp[jpc][mint][kpl]
                //                        << "\n";
                //^-
                for (int lp=0;lp<=kpl;++lp) {  // Loop local parameters lp, for half matrix
                  H(kpl,lp) += w2 * dsigDeldp[jpc][mint][kpl] * dsigDeldp[jpc][mint][lp];
                }
              } else { // ln cosh
                double dlncdwr = 1.0;  // d(ln cosh (wr))/d(wr)
                double abswr = std::abs(wib[mint] * r);
                if (abswr < RefineTargets::MAXCOSHARG) {
                  dlncdwr = tanh(wib[mint] * r);  // signed wr
                } else if (r < 0.0) {
                  dlncdwr = -1.0;
                }
                gradient[kpl] += wib[mint] * dlncdwr * dsigDeldp[jpc][mint][kpl];

                double d2dx = 1.0;
                // NB test against MAXCOSHARG already done for dlndx
                if (abswr > RefineTargets::MINCOSHARG) {
                  d2dx = dlncdwr / (wib[mint] * r);  // (1/d) tanh(d)
                  // exact 2nd derivative
                  //d2dx = 1.0/(cosh(wr)*cosh(wr));
                }
                for (int lp=0;lp<=kpl;++lp) {  // Loop local parameters lp, for half matrix
                  H(kpl,lp) += wib[mint] * d2dx * dsigDeldp[jpc][mint][kpl] * dsigDeldp[jpc][mint][lp];
                }
              }
            }
          } // end loop intensity bins
          //^
          //      std::cout << "k, Grad(k) no R2 " << kpl <<" "<<gradient[kpl] <<"\n";

          // Restraints: local parameter kpg, parameter group jpc
          gradient[kpl] += dr2dp[jpc][kpl]; // gradient
          //^     std::cout << "k, Grad(k)  R2   " << kpl <<" "<<gradient[kpl] <<"\n";
          for (int jpl=0;jpl<=kpl;++jpl) {  // Loop parameters jpl for this class
            H(jpl,kpl) += H2[jpc](jpl,kpl); // Hessian
          }
        }  // end loop local parameters k
        double R2 = sdmrestraintR[jpc];  // restraint residual R2
        double R = R1 + R2; // total residual
        double Rlsq = R1lsq + R2; // total residual

        //^
        //      std::cout << "Main residual: " << R1 << " Restraint Residual " << R2
        //                << " Total " << R <<"\n"; //^-

        // Symmetrise Hessian
        for (int k=0;k<npar-1;k++) {
          for (int l=k+1;l<npar;l++)
            {H(k,l) = H(l,k);} // other half
        }
        //^
        //      std::cout <<"Total Hessian:\n";
        //      for (int k=0;k<npar;k++) {
        //        for (int l=0;l<npar;l++)
        //          {std::cout <<" "<<H(l,k);}
        //        std::cout <<"\n";
        //      }
        //      // for npar = 2
        //      if (npar == 2) {
        //        double det = H(0,0)*H(1,1) - H(0,1)*H(1,0);
        //        std::cout << "Det [2x2]" << det <<"\n";
        //      }
        //^-

        // Variance/covariance matrix = (R/(m-n)) H^-1
        //   for m observations, n variables
        double rmn = R/(double(nintbins-npar));  // scaling factor

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
          ugradient[j] = gradient[j] * Uinv[j]; // scaled gradient
          for (int i=0;i<npar;++i) {
            A(i,j) = H(i,j) * Uinv[i] * Uinv[j];
          }}

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
        //^
        //      std::cout <<"Scaled Hessian:\n";
        //      for (int k=0;k<npar;k++) {
        //        for (int l=0;l<npar;l++)
        //          {std::cout <<" "<<A(l,k);}
        //        std::cout <<"\n";
        //      }
        //      std::cout << "Eigenvalues: ";
        //      for (int k=0;k<npar;k++) {
        //        std::cout <<" " << ev[k];
        //      }
        //      std::cout <<"\nDamp: " << damp <<"\n";;
        //^-1
        // Apply shifts
        for (int kpl=0;kpl<npar;++kpl) {  // Loop parameters kpl
          int kpg = kpl + idxpar;
          //^
          //      std::cout << "Shift: " << kpl <<" "<<kpg<<" "<<shifts[kpl]
          //                << " old parameter " <<sdmparams[kpg]
          //                << " new parameter " <<sdmparams[kpg]+shifts[kpl]<<"\n"; //^

          sdmparams[kpg] += shifts[kpl];
          parameterupdated[kpg] = true;
        }

        std::vector<double> sdpkpl(npar);
        for (int kpl=0;kpl<npar;++kpl) {  // Loop parameters kpl
          sdpkpl[kpl] = rmn * sqrt(Hinv(kpl,kpl));
          double relshift = shifts[kpl]/sdpkpl[kpl];
          bigrelshift = Max(bigrelshift, std::abs(relshift));
          //^
          //      std::cout <<"Parameter "<<kpl<<" sdpk " <<sdpkpl[kpl]<<" relshift " <<relshift
          //                <<" rmn " << rmn<<"\n"; //^-
        }
        target.Add(R1, R2);
        target.AddLsq(R1lsq);
      }
    } // end loop parameter classes

    target.converged = false;
    if (Update) {
      SDM.SetParameters(sdmparams, parameterupdated); // update parameters
      //      std::cout << "Biggest relative shift = " <<bigrelshift<<"\n";
      //      std::cout << "SDM: " << SDM.format() <<"\n"; //^
      if (bigrelshift < tolerance) {
        target.converged = true;
      }
    }
    return target;
  }
  // ---------------------------------------------------------
}
