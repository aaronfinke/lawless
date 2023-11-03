// selectscalingreflections.cpp
//
// Mark reflections to be used for scaling with accept/reject flags
//
// Possible criteria:
//
//  1) accept if all (accepted) observations have I/sd'(I) > IovSDmin
//     sd'(I) after correction by SDmodel
//
//  2) accept if |E^2| > E2min and < E2max after normalisation
//

#include "selectscalingreflections.hh"
#include "selectedobservations.hh"


namespace scala {

  // ------------------------------------------------------------
  std::vector<int> SelectScalingReflections(hkl_unmerge_list& hkl_list,
					    const SDmodel& SDM,
					    const ScaleModel& AllScales,
					    double& IovSDmin,
					    const Normalise& NormRes,
					    const double& E2min, const double& E2max)
  // On entry:
  //   hkl_list    reflection list, scales applied if needed
  //   SDM         sd correction model
  //   AllScales   scale model just used to set Phi bins for I/sd cutoff
  //   IovSDmin    minimum value for  <I>/sd'(<I>), == 0 no test
  //                < 0 negative value from default, to be reset here
  //               [default set to -3.0 in controls.cpp]
  //   NormRes     normalisation object
  //   E2min       minimum |E^2|, <= 0 no test
  //   E2max       maximum |E^2|, <= 0 no test
  //
  // On exit:
  //   hkl_list    reflection list, reflection accept flags updated
  //
  // returns number of reflections accepted, total and nskip
  //
  {
    hkl_list.ResetReflAccept();  // set to accept everything
    if (IovSDmin == 0.0 && E2min <= 0.0) {
      return std::vector<int> {0,0,0};
    }

    reflection this_refl, temp_refl;
    clipper::Array2d<int>      nI;
    int nskip = 1;

    if (IovSDmin < -0.001) {
      // If doing SD cutoff, do some histogramming on rotation ranges to get
      // a suitable sample
      std::vector<Run> runlist = hkl_list.RunList();   // runs
      int nruns = runlist.size();
      // Number of ranges
      std::vector<int> nranges_run(nruns,0);
      int nrotranges = 0; // total number of ranges
      std::vector<int> idxrun(nruns); // index to 1st rotation range for each run

      // Set up rotation ranges for each run
      for (int irun=0;irun<nruns;++irun) {
        // Store number of rotation ranges
        nranges_run[irun] = AllScales.primary_scale(irun).Nintervals();
        idxrun[irun] = nrotranges;
        nrotranges += nranges_run[irun];
        runlist[irun].PhiRange().SetNbin(Max(1,nranges_run[irun])); // set up binning on phi
      } // end loop runs

      const int NIOVSBINS = 13;  // bins on I/sig(I)
      //std::vector<double> IovSbins(NIOVSBINS);
      // lower bin limits
      double ibins[] = {0.0, 0.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 17.0,  20.0, 25.0};
      std::vector<double> IovSbins(ibins, ibins+NIOVSBINS);

      nrotranges = Max(1,nrotranges);
      nI.resize(nrotranges, NIOVSBINS);
      for (int i=0;i<nrotranges;++i) {
        for (int j=0;j<NIOVSBINS;++j) {nI(i,j) = 0;}}

      double IovSDmax = -100000.;
      observation this_obs;
      hkl_list.rewind();
      int index;
      while (hkl_list.next_reflection(this_refl) >= 0)  {  // loop reflections
        SDM.CorrectReflection(this_refl);
        while ((index = this_refl.next_observation(this_obs)) >= 0) {
          // loop observations
          int irun = this_obs.run();
          Rtype phi = this_obs.phi();
          int irot = runlist[irun].PhiRange().bin(phi) + idxrun[irun];
          double IovS = this_obs.I()/this_obs.sigI();
          IovSDmax = Max(IovSDmax, IovS);
          int isbin = -1;
          for (int i=NIOVSBINS-1;i>=0;--i) {
            if (IovS >= IovSbins[i]) {
              isbin = i;
              break;
            }
          }
          if (isbin < 0) isbin = 0;
          nI(irot, isbin)++;  // histogram on I/sigma for each rotation bin
        }
      }  // end loop reflections
      // Target number of observations in each rotation bin
      const int NOBSINBIN_SIGM = 600;   // for I/sigI cutoff
      const int NOBSINBIN_SKIP = 200;   // for skipping a fraction of reflections
      const double MINIOVSDMIN = 3.0;    // minimum value for I/sd for making nskip > 1
      const double STARTIOVSDMIN = 100000.0;
      IovSDmin = STARTIOVSDMIN;
      int nmin = 0;
      int nmin2 = 0;

      int nminall = 100000000;
      for (int ir=0;ir<nrotranges;++ir) { // loop rotation ranges
        int n = 0;
        int isbin = -1;
        for (int i=NIOVSBINS-1;i>=0;--i) { // loop bins backwards
          //std::cout << "ir, i, nI(ir,i) " <<ir<<" "<<i<<" "<<nI(ir,i)<<"\n";
          nmin2 = n;
          n += nI(ir,i);
          if (n >= NOBSINBIN_SIGM) {
            isbin = i;
            nmin = n;
            break;
          }
        }
        double cutoffr = 0.0;
        if (isbin > 0) {
          cutoffr = IovSbins[isbin];
          IovSDmin = Min(IovSDmin, cutoffr);
          nminall = Min(nminall, nmin);  // overall minimum number above cutoff
          //^
          //      std::cout << "Reset IovS " << IovSDmin << " ir " << ir
          //                << " nmin " <<nmin <<" isbin " << isbin
          //                << " nmin2 " << nmin2
          //                <<"\n";
          //^-
        }
      } // end loop rotation ranges

      nskip = 1;
      if (IovSDmin > STARTIOVSDMIN - 1.0) {  // minimum I/sd is too small
        // Set to default
        const double IOVSDMINDEFAULT = 2.0;
        IovSDmin = Min(IOVSDMINDEFAULT, 0.1*IovSDmax);
      } else {
        // if we still have many reflections above IovSDmin, then just use every
        // nskip'th reflection
        if (IovSDmin > MINIOVSDMIN) {
          nskip = Max(1,nminall/NOBSINBIN_SKIP);  // skip to every nskip'th reflection
        }
      }
      //^
      //      std::cout << "Revised cutoff I/sigI = " << IovSDmin
      //                << " nskip " << nskip << "\n";
    } // end doing SD cutoff


    // * * * * Loop all reflections
    int refRejected = 0;
    int nc = 0;  // counting for skip
    hkl_list.rewind();
    int naccepted = 0;
    while (hkl_list.next_reflection(this_refl) >= 0)  {
      //  Apply current SD correction to copy of reflection (all observations)
      temp_refl = this_refl;
      SDM.CorrectReflection(temp_refl);
      // Select all (I+ & I-) accepted observations for all datasets
      SelectedObservations allobs(temp_refl, -1, ALL);
      IsigI AvIsig = allobs.Average();  // average I, 1/variance weight
      int raccept = 0; // accept
      if (IovSDmin > 0.0) {
        // I/sd test
        if (AvIsig.sigI() <= 0.0) {
          raccept = +1; //reject
        } else if (AvIsig.I()/AvIsig.sigI() < IovSDmin) {
          raccept = +1; //reject
        } else if (nc%nskip != 0) { // skipping as requested
          raccept = +1; //reject
        }
      }
      if (E2min > 0.0) {
        double E2 = NormRes.apply(AvIsig.I(), temp_refl.invresolsq());
        if (E2 < E2min || E2 > E2max) {
          raccept = +1; //reject
        }
      }
      if (raccept != 0) {
        this_refl.SetStatus(raccept);
        hkl_list.replace_reflection(this_refl); // store updated reflection
        refRejected++;
      } else {
	naccepted++;
      }
      nc++;
    }  // end loop reflections
    std::vector<int> refnums = {naccepted, naccepted+refRejected, nskip};
    return refnums;
  }
}
