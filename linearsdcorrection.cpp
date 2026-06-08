//  linearsdcorrection.cpp
//
// DO NOT USE, doesn't work (yet)
//
//
// Observational equation
//   Eps = (Ihl - <Ih>)^2 - sigma'(hl)^2
//
// sigma'(hl)   = SdFac *sqrt[sigma(hl)^2 + SdB <Ih> + (SdAdd * <Ih>)^2]
// sigma'(hl)^2 = p * sigma(hl)^2  + q * <Ih> + r * <Ih>^2
// p = SdFac^2 ; q = SdFac^2 SdB; r = SdFac^2 SdAdd^2
//
// xT(hl) = (sigma(hl)^2  <Ih> <Ih>^2)
// (XT W X)ij = Sum(hl) (whl  xT(i) xT(j) )   (ie W matrix diagonal)
// yT(hl)  = (n/n-1)(Ihl - <Ih>)^2
//
//  p = (XT W X)^-1 XT W y
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

#include "linearsdcorrection.hh"
#include "string_util.hh"

namespace scala {
  // ---------------------------------------------------------
  // Refine SD correction model using linear LSQ, for each "bin class" separately
  LinearSDcorrection::LinearSDcorrection(SDmodel& SDM,
                                         const hkl_unmerge_list& hkl_list,
                                         const all_controls& controls,
                                         IntensityBin& irange,
                                         const double& tolerance, const double& rtolerance,
                                         const int&  max_cycles,
                                         phaser_io::Output& output)
  {
    init(SDM, hkl_list, controls, irange, tolerance, rtolerance, max_cycles, output);
  }
  // ---------------------------------------------------------
    void LinearSDcorrection::init(SDmodel& SDM,
                                  const hkl_unmerge_list& hkl_list,
                                  const all_controls& controls,
                                  IntensityBin& irange,
                                  const double& tolerance, const double& rtolerance,
                                  const int&  max_cycles,
                                  phaser_io::Output& output)
  {
    //^  std::cout << "SDM start: " << SDM.format() <<"\n"; //^

    std::vector<double> bestsdmparams = SDM.GetParameters();  // starting parameters
    output.logTab(0,LOGFILE,"\n");

    // No ties
    // Number of parameter groups in analysis
    npargroups =  SDM.Ngroups();
    linearlsq.resize(npargroups);

    int nmacrocycles = 0;
    // only one needed for linear method, but leave it here for testing stability
    const int MAXMACROCYCLES = 1;
    double R;
    std::vector<double> lastsdmparams;

    while (true) { // Macrocycles (maybe just one)
      nmacrocycles++;
      if (nmacrocycles > MAXMACROCYCLES) {break;}
      // Accumulate all sums from data
      bool anomalous = controls.anomalouscontrol.AnomalousSDcorr;
      for (size_t k=0; k<linearlsq.size(); k++) {
        linearlsq[k].clear(SDM.NparamsPerGroup());
      }
      SumsforSDcorrection(SDM, hkl_list, anomalous, irange);
      lastsdmparams = SDM.GetParameters();
      bool update = (max_cycles > 0); // don't update parameters if zero cycles
      UpdateParameters(SDM);
    } // loop macrocycles

  }
  // ---------------------------------------------------------
  void LinearSDcorrection::SumsforSDcorrection(const SDmodel& SDM,
                                               const hkl_unmerge_list& hkl_list,
                                               const bool& anomalous,
                                               IntensityBin& irange)
  // Accumulate sums for SD correction refinement
  // anomalous   true to separate anomalous I+ & I- (usually true)
  {
    // Set weight type
    //    WeightType::AverageWeightType weighttype = SDM.Weight();
    WeightType::AverageWeightType weighttype = WeightType::SCALE;
    int Ndatasets = hkl_list.num_datasets();

    reflection this_refl;
    int nref = 0;
    hkl_list.rewind();
    std::vector<float> delta, delta2;
    //^
    //    std::cout << "SDcorr\na b d^2/var I^2/var weight\n";

    //  deviations within each run & full/partial
    while (hkl_list.next_reflection(this_refl) >= 0)  {
      bool Centric = hkl_list.symmetry().is_centric(this_refl.hkl());
      // Correct sds in this_refl, return uncorrected scaled values
      //      std::vector<float> sig0 = SDM.CorrectReflection(this_refl);
      // Average I <I> over all observations
      SelectedObservations selobs(this_refl, -1, ALL, weighttype);
      //^
      //      std::cout << "rsd " << this_refl.hkl().format()
      //                << " N " << selobs.Number() << "\n";
      //^-
      nref++;
      int nacc = 0; //^
      for (int id=0;id<Ndatasets;id++) {    // loop datasets
        if (Centric || !anomalous) {
          // No anomalous, treat all observations together
          if (Ndatasets > 1) {selobs.init(this_refl, id, ALL, weighttype);} // already done if 1 dataset
          if (selobs.Number() > 1) { // at least 2 observations
            // average intensity for SD correction
            double Iav = selobs.Average().I();
            addContribution(selobs, Iav, SDM);
            nacc++;
          }
        } else {
          // Anomalous, treat I+ & I- separately
          selobs.init(this_refl, id, IPLUS, weighttype);
          if (selobs.Number() > 1) {
            // average intensity for SD correction
            double Iav = selobs.Average().I();
            addContribution(selobs, Iav, SDM);
            nacc++;
          }
          selobs.init(this_refl, id, IMINUS, weighttype);
          if (selobs.Number() > 1) {
            // average intensity for SD correction
            double Iav = selobs.Average().I();
            addContribution(selobs, Iav, SDM);
            nacc++;
          }
        } // end acentric
      } // end loop datasets
    } // end loop reflections
    //    std::cout << "end sdcorr\n";
  }
  //---------------------------------------------------------------
  void LinearSDcorrection::addContribution(SelectedObservations& selobs,
                                           const double& Iav,
                                           const SDmodel& SDM)
  // Add in contribution to LSQ objects for this observation set
  {
    std::vector<float> delI = selobs.DelI();
    std::vector<double> weights = selobs.weights();
    //std::vector<double> weights(selobs.Nobs(),1.0);
    // "coordinate" from SDmodel, ie (var(I), [<I>,] <I>2)
    std::vector<double> xv;
    int idxgroup;
    int nused = selobs.Number();
    if (nused < 2) {return;}
    double fac = double(nused)/double(nused-1);
    const double MAXVAR = 10.0;
    //^
    reflection refl = selobs.Reflection();
    std::vector<float> delta2= selobs.Delta2(false);
    std::vector<float> sigmai = selobs.sigmaI();
    std::vector<IsigI> isigilist = selobs.IsigIlist();
    //    std::cout <<"ref "<<refl.hkl().format() <<" sig1 "<<sigmai[0] <<"\n";
    //^-
    for (size_t i=0;i<delI.size();++i) {
      if (delI[i] != 0.0) { // valid delta
        xv = SDM.coordinate(selobs, i, Iav, idxgroup); // get coordinate and group index
        //^
        //      xv[0] = sigmai[i]*sigmai[i];
        double y = fac*delI[i]*delI[i];
        if (y/xv[0] < MAXVAR) {
          //^
          double sdprime = SDM.sdCorrected(selobs, i, Iav);
          std::string flag = "  ";
          if (y > 6.0e+9) {
            flag = " *";
          }
          double chisq = y/xv[0];
          double diff = y - xv[0];
          printf("@@ %9.1f %9.4f %9.4f %9.0f %9.4f %9.1f %4d %s %s\n",
                 Iav, y, xv[0], xv.back(), weights[i],
                 isigilist[i].I(), nused,
                 refl.hkl().format().c_str(),
                 flag.c_str());
          //          printf("&$ %9.0f %9.4f %9.4f %9.0f %9.4f %9.4f %9.4f %s %s\n",
          //                 Iav, y, xv[0], xv.back(), weights[i], chisq, diff,
          //             refl.hkl().format().c_str(),
          //                 flag.c_str());
          //             Iav, sqrt(y), sqrt(xv[0]), sdprime, delta2[i]);
          //          printf("d/v, I2/v %9.4f %9.4f %9.4f %12.4g\n",
          //             delta2[i]*delta2[i], y/xv[0], xv.back()/xv[0], weights[i]);
          //    std::cout <<"addContribution "<<y<<" "<<weights[i];
          //    for (size_t k=0; k<xv.size(); k++) {
          //      std::cout <<" "<<xv[k];
          //    }
          //    std::cout <<std::endl;
          //^-
          double w = weights[i] * 0.00001;
          linearlsq.at(idxgroup).add(y, xv, w);
          //          linearlsq.at(idxgroup).add(y, xv, weights[i]);

        }
      }
    }
  }
  // ---------------------------------------------------------
  void LinearSDcorrection::UpdateParameters(SDmodel& SDM)
  // Solve LSQ equations for each "parameter class" to update SDmodel SDM,
  // Parameter classes are run and full/partial, subject to the Allrunssame flag and
  // few/no fulls/partials selections (useflags in SDM)
  //
  // Each parameter class is treated independently: this is possible since the
  //  weighted deviations delta within each class are independent, as they are calculated
  //  relative to <Ih>av, which is independent of sd(I)
  //  (at least for 1/sqrtscale or 1/scale weighting)
  //
  {
    //    std::cout <<"npargroups "<<npargroups<<std::endl;
    std::vector<double> allparams;
    for (int jpc=0;jpc<npargroups;++jpc) { // loop parameter classes
      std::vector<double> groupparams = linearlsq[jpc].solve();
      //^
      std::cout <<"Parameter group "<<jpc;
      for (size_t k=0; k<groupparams.size(); k++) {
        std::cout <<" "<<groupparams[k];
      }
      std::cout <<std::endl;
      std::cout << "Nobs "<< linearlsq[jpc].Nobs() << "\n";
      //^-
      allparams.insert(allparams.end(), groupparams.begin(), groupparams.end());
    } // end loop parameter classes

    SDM.SetParameters(allparams);

  }
  // ---------------------------------------------------------
}
