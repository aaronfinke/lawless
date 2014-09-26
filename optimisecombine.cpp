// optimisecombine.cpp
//
// For Mosflm data with two intensity estimates, summation & profile fitted,
// Determine the best one to use or the best cross-over point for the combine option
//
//  Combine integrated intensity Iint and profile-fitted intensity Ipr as
//         I = w * Ipr  +  (1-w) * Iint
//       where w = 1/(1+(Iint/Imid)**IPOWER)
//       Imid (= IcolFlag)  is the point of equal weight, and
//       should be about the mean intensity.
//  Ipower default = 3
//


#include "optimisecombine.hh"
#include "scala_util.hh"
#include "string_util.hh"

using phaser_io::LOGFILE;

namespace scala {
  //--------------------------------------------------------------
  OptimiseCombine::OptimiseCombine(const hkl_unmerge_list& hkl_list,
                                   phaser_io::Output& output)
  {
    init(hkl_list, output);
  }
  //--------------------------------------------------------------
  void OptimiseCombine::init(const hkl_unmerge_list& hkl_list,
                             phaser_io::Output& output)
  {
    /// If not both profile & summation, then return!
    bothIpresent = SelectI::IsIprPresent();
    if (!bothIpresent) {return;}

    output.logTab(0,LOGFILE,
      "\n\n========= Optimising selection of intensity estimate =========\n");

    output.logTab(0,LOGFILE,
                  std::string("\nThe input HKLIN file contains two estimates of intensity,\n")+
                  "a summation integration value Isum and a profile-fitted value Ipr\n"+
                  "The optimisation here chooses the value which gives the smallest overall Rmeas:\n"+
                  "either Ipr, Isum or a combination of the two based on the Iraw value, where\n"+
                  "Iraw is the Isum value back-corrected for Lorentz and Polarisation\n");

    //Mean raw I    (ie mean I/LP)
    meanI = RawIntensityDistribution(hkl_list);
    Itop = meanI * 1.6;  // top of search range

    ResRange = hkl_list.ResLimRange();
    nresbin =  ResRange.Nbins();


    // We always have 3 scores
    std::vector<double> scores(3);
    // ... are R by resolution
    std::vector<std::vector<Rfactor> > Rfacs(3);
    // ... and their type: = 0 summation, -1 profile, > 0 imid combine
    std::vector<double> imid(3);

    output.logTabPrintf(0, LOGFILE,
                        "\nMean Iraw for all data: %10.1f\n", meanI);

    RPair inner = ResRange.boundsA(0);          // resolution range for inner bin
    std::string s1 = StringUtil::Strip(StringUtil::ftos(inner.first,8,1)+"-"+
                                       StringUtil::ftos(inner.second,8,2));
    RPair outer = ResRange.boundsA(nresbin-1);  // resolution range for outer bin
    std::string s2 = StringUtil::Strip(StringUtil::ftos(outer.first,8,2)+"-"+
                                       StringUtil::ftos(outer.second,8,2));

    output.logTab(0, LOGFILE,
                  "\nRmeas also printed for inner and outer resolution ranges");

    output.logTab(0, LOGFILE,
          "\n                         Intensity type        Rmeas       Inner       Outer");
    output.logTab(0, LOGFILE,
                  "                     Resolution range (A)        All"+
                  StringUtil::RightString(s1,12)+StringUtil::RightString(s2,12));
    output.logTab(0, LOGFILE," ");

    // summation I
    imid[0] = 0.0;
    scores[0] = TestValue(imid[0], hkl_list, Rfacs[0]);
    PrintR(imid[0], scores[0], Rfacs[0], output);

    // profile I
    imid[1] = -1.0;
    scores[1] = TestValue(imid[1], hkl_list, Rfacs[1]);
    PrintR(imid[1], scores[1], Rfacs[1], output);

    // Now search for best Combine::Imid value
    // Start at Itop/2
    imid[2] = Itop/2.0;
    scores[2] = TestValue(imid[2], hkl_list, Rfacs[2]);
    PrintR(imid[2], scores[2], Rfacs[2], output);

    bool KeepGoing = true;
    double TolFrac = 0.04;       // fractional tolerance
    double TolDelScore = 0.0;    // tolerance on score difference, see below
    double TolLimit = 0.1*Itop;  // end limits to select extremes
    int maxcycles = 5;           // maximum number of cycles
    int ib = 0;
    int nc = 0;
    int lbest = -1; // previous best
    double bestscore;

    while (KeepGoing) {
      // find worst and best (low is good)
      double scw = scores[0];
      int iw = 0;
      double scb = scores[0];
      ib = 0;
      for (int i=1;i<3;++i) {
        if (scores[i] > scw) {
          scw = scores[i];
          iw = i;
        }
        if (scores[i] < scb) {
          scb = scores[i];
          ib = i;
        }
      }
      int im = 3 - (iw+ib); // middle one
      if (nc == 0) { // 1st cycle, work out convergence criterion
        bestscore = scores[ib];
        double dsc = scores[im] - scores[ib];  //  difference of best to middle
        TolDelScore = dsc * TolFrac;
        //^
        //      std::cout << "Best, med, worst " <<
        //        imid[ib]<<" "<<scores[ib] << " " <<
        //        imid[im]<<" "<<scores[im] << " " <<
        //        imid[iw]<<" "<<scores[iw] << "\n";
      } else {
        //^
        //      std::cout << "Best, med, worst " <<
        //        imid[ib]<<" "<<scores[ib] << " " <<
        //        imid[im]<<" "<<scores[im] << " " <<
        //        imid[iw]<<" "<<scores[iw] << "\n";
        //      std::cout << "lbest, ib, imid[im], scores[im]-bestscore "
        //                <<lbest<<" "<<ib<<" "<<imid[im]<<" "<<scores[im]-bestscore<<"\n";
        //      std::cout << "best, diff, tol "
        //                << bestscore <<" "<<bestscore - scores[ib]
        //                <<" "<<TolDelScore<<"\n"; //^
        // Criteria for convergence:
        //  if we have a new best score, then test improvement from last best
        //  if the best one hasn't change, test improvement from the middle one
        if (((lbest != ib) && (bestscore - scores[ib]) < TolDelScore) ||
            ((lbest == ib) && ((scores[im] - bestscore) < TolDelScore))) {
          //  quit
          KeepGoing = false;
        }
        if (scores[ib] < bestscore) {bestscore = scores[ib];}
      }
      if (KeepGoing) {
        if (imid[ib] > 0.0) { // best is Combine
          if (imid[im] == 0.0) { // middle is Isum
            // replace worst with half-way between 0 & best
            imid[iw] = 0.5*imid[ib];
          } else if (imid[im] < 0.0) { // middle is Ipr
            // replace worst with half-way between top & best
            imid[iw] = 0.5*(imid[ib] + Itop);
          } else {
            // replace worst with half-way between middle & best
            imid[iw] = 0.5*(imid[ib] + imid[im]);
          }
        } else if (imid[ib] < 0.0) { // best is Profile
          if (imid[im] == 0.0) { // middle is Isum
            // replace worst with half-way between top & worst
            imid[iw] = 0.5*(imid[iw] + Itop);
          } else { // middle is combine
            if ((Itop-imid[im]) < TolLimit) { // hit top limit
              KeepGoing = false;
            }
            // replace worst with half-way between middle & top
            imid[iw] = 0.5*(Itop + imid[im]);
          }
        } else if (imid[ib] == 0.0) { // best is Isum
          if (imid[im] < 0.0) { // middle is Ipr
            // replace worst with half-way between 0 & worst
            imid[iw] = 0.5*imid[iw];
          } else { // middle is combine
            if (imid[im] < TolLimit) { // hit bottom limit
              KeepGoing = false;
            }
            // replace worst with half-way between middle & 0
            imid[iw] = 0.5*imid[im];
          }
        }
        scores[iw] = TestValue(imid[iw], hkl_list, Rfacs[iw]); // rescore
        if (KeepGoing) {
          nc++;
          PrintR(imid[iw], scores[iw], Rfacs[iw], output);
          if (nc > maxcycles) {KeepGoing = false;}
        }
      }
      lbest = ib;
    } // end while
    output.logTab(0,LOGFILE,"\nBest value:");
    PrintR(imid[ib], scores[ib], Rfacs[ib], output);
    SelectI::SetIcolFlag(int(imid[ib]), imid[ib]);
    output.logTab(0, LOGFILE, "\n"+SelectI::format());
  }
  //--------------------------------------------------------------
  double OptimiseCombine::TestValue(const double& imid, const hkl_unmerge_list& hkl_list,
                                    std::vector<Rfactor>& Rmeas)
  // Set imid value (= 0 Isum, < 0 Ipr, > 0 combine), set Rmeas by resolution &
  // return overall Rmeas
  {
    SelectI::SetIcolFlag(int(imid), imid);  // combine
    Rmeas = GetScores(hkl_list);
    Rfactor R;
    for (size_t i=0;i<Rmeas.size();++i) {
      R += Rmeas[i];
    }
    return R.R();
  }
  //--------------------------------------------------------------
  void OptimiseCombine::PrintR(const double& flag,
                               const double& score, const std::vector<Rfactor>& Rmeas,
                               phaser_io::Output& output)
  // flag = 0 Isum, < 0 Ipr, > 0 combine
  {
    std::string s;
    if (flag == 0) {
      s = "Summation intensities";
    } else if (flag < 0) {
      s = "Profile intensities";
    } else {
      s = "Combined intensities Imid = "+StringUtil::Strip(StringUtil::ftos(flag,10,0));
    }
    output.logTabPrintf(0,LOGFILE,
                        "%40s %11.4f %11.4f %11.4f\n",
                        s.c_str(), score, Rmeas[0].R(), Rmeas.back().R());
  }
//--------------------------------------------------------------
  std::vector<Rfactor> OptimiseCombine::GetScores(const hkl_unmerge_list& hkl_list) const
  // Only comes in here if we have both Isummation and Ipr
  // Get various scores
  {
    reflection this_refl;
    observation this_obs;

    int ndatasets = hkl_list.num_datasets();
    std::vector<std::vector<double> > intensities(ndatasets);

    std::vector<Rfactor> rmeasRes(nresbin);  // Rmeas, resolution bins, all datasets

    hkl_list.rewind();

    while (hkl_list.next_reflection(this_refl) >= 0)  {   // * * * * Loop accepted reflections
      for (int idts=0;idts<ndatasets;++idts) { // loop datasets to clear
        intensities[idts].clear();
      }
      Rtype invresolsq = this_refl.invresolsq();
      int mres = ResRange.bin(invresolsq);

      //  Loop all observations
      for (int i=0;i<this_refl.num_observations();++i) {
        this_obs = this_refl.get_observation(i);
        if (this_obs.IsAccepted()) {
          this_obs.sum_partials();  // with current SelectI settings
          if (this_obs.sigI() > 0.0) {
            double I = this_obs.kI();   // scaled I
            int didx = this_obs.datasetIndex();     // dataset
            intensities[didx].push_back(I);
          }
        }
      }
      for (int idts=0;idts<ndatasets;++idts) { // loop datasets for statistics
        if (intensities[idts].size() > 1) {
          double AvI = MnSd(intensities[idts]).first; // mean I
          double an = intensities[idts].size();
          double w = sqrt(an/(an-1.0));
          for (size_t i=0;i<intensities[idts].size();++i) { // loop observations
            double delI = intensities[idts][i] - AvI;
            rmeasRes[mres].add(delI, AvI, w);  // Rmeas in resolution bins
          }
        }
      } // end loop datasets
    } // end loop reflections
    return rmeasRes;
  }
  //--------------------------------------------------------------
  double OptimiseCombine::RawIntensityDistribution(const hkl_unmerge_list& hkl_list) const
  // get mean of raw intensity distribution
  {
    reflection this_refl;
    observation this_obs;
    hkl_list.rewind();
    MeanSD mnI;
    bool hasLP = hkl_list.DataFlags().is_LP;

    while (hkl_list.next_reflection(this_refl) >= 0)  {   // * * * * Loop accepted reflections
      //  Loop all observations
      for (int i=0;i<this_refl.num_observations();++i) {
        this_obs = this_refl.get_observation(i);
        if (this_obs.IsAccepted()) {
          this_obs.sum_partials();  // with current SelectI settings
          double I = this_obs.kI();   // scaled I
          double lp = this_obs.LP();
          if (hasLP) {
            if (lp > 0.0) {
              double Iraw = I/lp;  // raw intensity
              mnI.Add(Iraw);
            }
          } else {
            mnI.Add(I);
          }
        }
      }  // end loop observations
    } // end loop reflections
    return mnI.Mean();
  }
  //--------------------------------------------------------------
}
