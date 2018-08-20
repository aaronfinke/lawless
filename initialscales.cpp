//
// initialscales.cpp
//
// Initial rough scaling by making average intensities equal
//

#include "initialscales.hh"
#include "foxholmes.hh"
#include "string_util.hh"

#include <assert.h>
#define ASSERT assert

using phaser_io::LOGFILE;
using phaser_io::LXML;

namespace scala {
  // ---------------------------------------------------------
  InitialData::InitialData(const clipper::Array2d<double>& AvI)
  {
    // Store address of data array
    // Data is a 2D array(nrotranges, nresbins)
    avi = &AvI;
    npall = avi->rows();  // number of parameters = number of rotation ranges
    np = 0;
    validranges.assign(npall, true);
    rangecount.assign(npall, 0);
    for (int i=0;i<avi->rows();++i) {  // loop parameters (rotation ranges)
      next = -1;
      bool empty = true;
      while (++next < avi->cols()) { // loop resolution ranges
        if ((*avi)(i,next) > 0.0) {
          empty = false;
          rangecount[i]++;
        }
      }
      if (empty) {
        // this rotation range is empty
        validranges[i] = false;
      } else {
        np++;    // count valid rotation ranges
      }
    }
    next = -1;
  }
  // ---------------------------------------------------------
  bool InitialData::ObsArray(std::vector<DPair>& obs) const
  // return one observation, I, sigma pair, length fnp
  // Return false if end of data
  {
    ASSERT (int(obs.size()) >= avi->rows());
    double sd1 = 1.0;  // dummy sigma
    double I, sd;
    while (++next < avi->cols()) {
      for (int i=0;i<avi->rows();++i) {
        I = (*avi)(i,next);
        sd = (I>0.0) ? sd1  : 0.0;  // sd = 0 if I = 0
        obs[i] = DPair(I, sd);
        //^
        //std::cout <<"ObsArray i, obs " <<i<<" "<<obs[i].first<<" "<<obs[i].second<<"\n";
      }
      return true;
    }
    next = -1;
    return false;
  }
  // ---------------------------------------------------------
  // ---------------------------------------------------------
  InitialScales::InitialScales(hkl_unmerge_list& hkl_list,
                               ScaleModel& AllScales,
                               const bool& determineScales,
                               const all_controls& controls,
                               phaser_io::Output& output)
  // Get initial estimates of primary scales, from making intensity
  // averages equal
  // if !determineScales, just check for sufficient overlaps etc
  {
    status = 0;
    minimumoverlap = +1.00001;
    if (determineScales) {
      output.logTab(0,LOGFILE,
                    "\n========= Initial scaling =========\n\n");
    } else {
      output.logTab(0,LOGFILE,
                    "\n========= Checking for scaling overlaps =========\n\n");
    }
    std::vector<Run> runlist = hkl_list.RunList();   // runs
    int nruns = runlist.size();
    // true if a run has batch scales
    std::vector<bool> batch_scale_run(nruns, false);
    // Number of ranges, for batch scale = Nbatches, else number of scales-1
    std::vector<int> nranges_run(nruns,0);

    nrotranges = 0; // total number of ranges
    std::vector<int> idxrun(nruns); // index to 1st rotation range for each run
    // Set up rotation ranges for each run
    for (int irun=0;irun<nruns;++irun) {
      // Store number of rotation ranges
      //   for batch mode = Nbatches (excluding rejected ones)
      nranges_run[irun] = AllScales.primary_scale(irun).Nintervals();
      if (AllScales.primary_scale(irun).IsBatchScale()) {
        // Batch scale
        batch_scale_run[irun] = true;
      }
      // idxrun is 1st index in list for this run
      idxrun[irun] = nrotranges;
      nrotranges += nranges_run[irun];
      runlist[irun].PhiRange().SetNbin(nranges_run[irun]); // set up binning on phi
      //^
      //      std::cout << "Phirange run " << irun << " " << runlist[irun].PhiRange().format();
      //^-
    }


    // Resolution ranges
    ResoRange resrange = hkl_list.ResRange();
    if (controls.analysis.NresoBins() > 0) {
      // set number of ranges
      resrange.SetNbins(controls.analysis.NresoBins());
    }
    reflection this_refl;
    observation this_obs;
    int index;

    int nrbins =  resrange.Nbins();
    // sum->mean I (rotation, resolution)
    // try weighted mean
    clipper::Array2d<MeanValue> meanI(nrotranges, nrbins);
    // numbers
    clipper::Array2d<int>      nI(nrotranges, nrbins);

    //  Noverlap    number of observations which overlap with at least
    //              one other rotation range
    std::vector<int> noverlaprot(nrotranges, 0);
    //  NotOverlap  number of observations which are singletons or only
    //              overlap within the rotation range
    std::vector<int> notoverlap(nrotranges, 0);
    // number of reflections with multiple observations only within same rotrange
    nummultipleobsrotrange.assign(nrotranges, 0);
    // list of rotation range indices for each reflection
    std::vector<int> irotlist; // faster to create once and clear each iteration

    for (int i=0;i<nrotranges;++i) {
      for (int j=0;j<nrbins;++j) {
        nI(i,j) = 0;}
    }

    hkl_list.rewind();
    int irot, irot1;
    while (hkl_list.next_reflection(this_refl) >= 0)  {  // loop reflections
      int ires = resrange.bin(this_refl.invresolsq());
      int nobs = 0;
      irotlist.clear();
      irot1 = -1;
      while ((index = this_refl.next_observation(this_obs)) >= 0) {
        // loop observations
        int irun = this_obs.run();
        Rtype phi = this_obs.phi();
        if (AllScales.primary_scale(irun).IsBatchScale()) {
          // batch scale
          int batchnum = this_obs.Batch();
          // serial index in run
          int batchserial = AllScales.primary_scale(irun).batchSerialIndex(batchnum);
          irot = batchserial + idxrun[irun];
        } else {
          irot = runlist[irun].PhiRange().bin(phi) + idxrun[irun];
        }
        double weight = 1.0;
        double sigi = this_obs.sigI();
        if (sigi > 0.0) {
          weight = 1.0/(sigi*sigi);
        }
        meanI(irot, ires).Add(double(this_obs.I()), weight);
        nI(irot, ires)++;
        nobs++;
        irotlist.push_back(irot);
        if (irot1 == -1) {
          irot1 = irot;
        } else if (irot1 >= 0) {
          if (irot != irot1) {
            irot1 = -2; // flag for different irot values
          }
        }
      } // observation loop

      //^^
      //      std::cout << this_refl.hkl().format() <<" "<<irotlist.size()<<
      //        " "<<nobs<<" "<<irot1<<":";
      //      for (size_t k=0; k<irotlist.size(); k++) {
      //        std::cout <<" "<< irotlist[k];
      //      }
      //      std::cout <<"\n";
      //^-

      if (irot1 == -2) {
        // there are observations from different rotation ranges
        // count them
        for (size_t k=0; k<irotlist.size(); k++) {
          noverlaprot[irotlist[k]]++;
        }
      } else {
        // there are only observations from same rotation range
        // count them
        for (size_t k=0; k<irotlist.size(); k++) {
          notoverlap[irotlist[k]]++;
        }
        if (nobs > 1) {
          // number of observations with multiple observations only within same rotrange
          nummultipleobsrotrange[irotlist[0]] += nobs;
        }
      }
    } // reflection

    // Compute average intensities
    numobsrotrange.assign(nrotranges, 0); // number of observations for rotrange
    clipper::Array2d<double> sumI(nrotranges, nrbins);
    for (int i=0;i<nrotranges;++i) {
      for (int j=0;j<nrbins;++j) {
        if (nI(i,j) > 0) {
          sumI(i,j) = meanI(i,j).Mean();
          numobsrotrange[i] += nI(i,j);
        }
      }
    }

    // average overlap fraction by rotation range
    std::vector<bool> validranges(nrotranges, false);
    int nvalidranges = 0;
    fractionaloverlapbyrotrange.assign(nrotranges, 0.0);
    int noverlap = 0;
    int nnot = 0;
    for (int i=0;i<nrotranges;++i) {
      int ntot = noverlaprot[i]+notoverlap[i];
      if (ntot > 0) {
        fractionaloverlapbyrotrange[i] =
          double(noverlaprot[i])/double(noverlaprot[i]+notoverlap[i]);
        validranges[i] = true;
        nvalidranges++;
        minimumoverlap = std::min(minimumoverlap, fractionaloverlapbyrotrange[i]);
      }

      noverlap += noverlaprot[i];
      nnot += notoverlap[i];
    }
    averageoverlap = 0.0;
    if (nnot > 0) {
      averageoverlap = double(noverlap)/double(noverlap+nnot);
    }

    if (nrotranges == 1) {
      // special for one range, no scaling
      fractionaloverlapbyrotrange[0] =
         double(nummultipleobsrotrange[0])/double(nummultipleobsrotrange[0]+notoverlap[0]);
      averageoverlap = fractionaloverlapbyrotrange[0];
      output.logTab(0, LOGFILE,
                    std::string("\nOnly one rotation range, initial scales set to 1.0\n")+
                    "  Fractional overlap = Noverlapped/Ntotal = "+
                    StringUtil::ftos(averageoverlap,5,2)+
                    "\n   where Noverlapped is the number of observations"+
                    " with equivalent observations "+
                    "("+StringUtil::itos(nummultipleobsrotrange[0])+"),\n"+
                    "   and Ntotal is the total number of observations\n");

      // Store initial scales
      std::vector<double> gscales(nrotranges, 1.0);  // scales
      AllScales.SetInitialScales(gscales, numobsrotrange);
      status = nrotranges;
      return;
    }

    output.logTab(0, LOGFILE,
                  std::string("\nThe average fractional overlap = Noverlapped/Ntotal, ")+
                  "where Noverlapped is the number of observations\n"+
                  "with equivalent observations in a different rotation range, "+
                  "and Ntotal is the total number of observations\n");

    const int NPERLINE = 10;
    for (int irun=0;irun<nruns;++irun) {
      bool empty = false;
      int i1 = idxrun[irun];
      int i2 = nrotranges;
      if (irun < nruns-1) i2 = idxrun[irun+1];  // not last run
      output.logTabPrintf(0, LOGFILE,
         "\nAverage fractional overlap between rotation ranges for run %5d\n",
                          runlist[irun].RunNumber());
      for (int i=i1;i<i2;++i) {
        if ((i-i1) > 0 && (i-i1)%NPERLINE == 0) output.logTabPrintf(0,LOGFILE,"\n");
        output.logTabPrintf(0,LOGFILE," %9.2f", fractionaloverlapbyrotrange[i]);
      }
      output.logTabPrintf(0,LOGFILE,"\n");
    }
    output.logTabPrintf(0, LOGFILE,
       "\nOverall fractional overlap between rotation ranges %5.2f, minimum %5.2f\n",
                        averageoverlap, minimumoverlap);

    if (!determineScales) {return;}

    InitialData data(sumI);  // make data accessible to scale refinement

    FoxHolmes fh(data);
    int Ncyc =  5;  // number of cycles
    phaser::protocolPtr cPtr(new phaser::ProtocolScale(Ncyc));   // default protocols
    //^   std::cout << " InitScales::Minimizer::run::Number of cycles " << cPtr->getNCYC() << "\n";
    phaser::Minimizer Min;

    Min.run(fh, cPtr, output);

    // Print initial scales

    gscales = fh.getGscales();  // inverse scales (g)
    ASSERT (gscales.size() == size_t(nrotranges));
    std::vector<double> scales(nrotranges, 0.0);  // scales

    validranges = data.validRanges();

    for (int irun=0;irun<nruns;++irun) {
      output.logTabPrintf(0, LOGFILE,
                          "\nInitial scales for run %5d\n", runlist[irun].RunNumber());
      bool empty = false;
      int i1 = idxrun[irun];
      int i2 = nrotranges;
      if (irun < nruns-1) i2 = idxrun[irun+1];  // not last run
      std::vector<bool> validinrun(i2-i1, true);
      // scales for irun go from i1 to i2-1
      for (int i=i1;i<i2;++i) {
        if (gscales[i] == 0.0) {
          empty = true;  // empty slot
          validinrun[i] = false;
        }
      }

      if (empty) {
        // this run is missing at least one scale
        int k = i1;
        while (k < i2) {
          if (!validinrun[k]) { // k'th scale is missing
            // try to find valid ones to fill in
            // hunt backwards to k1
            int k1 = -1;
            if (k > i1) {
              for (int j=k-1;j>=i1;j--) {
                if (validinrun[j]) {
                  k1 = j;  // k1 is previous valid scale
                  break;
                }
              }
            }
            int k2 = -1;  // hunt forwards
            if (k < i2-1) {
              for (int j=k+1;j<i2;j++) {
                if (validinrun[j]) {
                  k2 = j;  // k2 is next valid scale
                  break;
                }
              }
            }
            // we need to fill in from k to k2-1
            floatType g = 0.0;
            int ng = 0;
            if (k1 >= 0) {g += gscales[k1]; ng++;}
            if (k2 >= 0) {g += gscales[k2]; ng++;}
            else {k2 = i2;}
            if (ng > 0) {
              g /= floatType(ng);
              for (int j=k;j<k2;j++) {  // to k2-1 or i2-1
                gscales[j] = g;
                //std::cout <<"fill in " << j <<" with "<<1.0/g<<"\n"; //^
              }
            }
            k = k2-1;
          }
          k++;
        }
      }
      for (int i=i1;i<i2;++i) {
        if ((i-i1) > 0 && (i-i1)%NPERLINE == 0) output.logTabPrintf(0,LOGFILE,"\n");
        if (gscales[i] != 0.0) {
          scales[i] = 1./gscales[i];
        }
        output.logTabPrintf(0,LOGFILE," %9.3f", scales[i]);
      }
      output.logTabPrintf(0,LOGFILE,"\n");

    }  // end loop runs
    // Store initial scales
    AllScales.SetInitialScales(gscales, numobsrotrange);
    status = nrotranges;
  }  // InitialScales
  // ---------------------------------------------------------
  // return true if there seems to be enough data to refine scales
  bool InitialScales::enoughData(const double& minimum_overlap,
                                 const int& maximum_gap)
  // Count for each rotation range:
  //  Noverlap    number of observations which overlap with at least
  //              one other rotation range
  //  NotOverlap  number of observations which are singletons or only
  //              overlap within the rotation range
  // then fractionaloverlap = Noverlap / (Noverlap + NotOverlap)
  // This should be > abs(minimum_overlap)
  // sets overlapthreshold & allowedgap (mutable)
  // If minimum_overlap <= 0.0, no check is made
  {
    enoughdata = true;
    overlapthreshold = std::abs(minimum_overlap);
    allowedgap = maximum_gap;  // allow [default two] below threshold
    //    if (overlapthreshold <= 0.0) {
    //      return enoughdata;
    //    }
    ASSERT (nrotranges == int(fractionaloverlapbyrotrange.size()));
    if (nrotranges == 1) {
      // special for 1 range, always allow here
      if (fractionaloverlapbyrotrange[0] < overlapthreshold) {
        enoughdata = false;
      }
    } else if (nrotranges > 1) {
      int ngap = 0;
      bool ingap = false;
      int nbad = 0;
      // count ranges with low multiplicity
      for (size_t ir=0; ir<nrotranges; ir++) {
        if (fractionaloverlapbyrotrange[ir] < overlapthreshold) {
          // too few
          nbad++; // count contiguous bad ranges
          if (!ingap) {
            ingap = true;
          }
        } else {
          if (ingap && (nbad > allowedgap)) {
            // found a range longer than allowedgap, count them
            ngap++;
          }
          ingap = false;
          nbad = 0;
        }
      }
      if (ingap && (nbad > allowedgap)) { // possible gap at end
        // found a range longer than allowedgap, count them
        ngap++;
      }
      if (ngap > 0) {
        enoughdata = false;
      }
    }
    return enoughdata;
  }
  // ---------------------------------------------------------
  void InitialScales::reportOverlapXML(phaser_io::Output& output,
                                       const bool& allowgap) const
  // Report overlap status information to XML
  //  allowgap true if overlap gap check suppressed
  {
    if (overlapthreshold < 0.0) {
      return;  // no information
    }
    output.logTab(0, LXML,"<RotationalOverlap>");

    if (enoughdata) {
      // Yes sufficient data overlap
      output.logTab(1, LXML,
                    StringUtil::MakeXMLtag("EnoughOverlap","True"));
    } else {
      // No insufficient data overlap
      output.logTab(1, LXML,
                    StringUtil::MakeXMLtag("EnoughOverlap","False"));
      if (allowgap) {
        output.logTab(1, LXML,
                      StringUtil::MakeXMLtag("GapAllowed","True"));
      }
    }
    output.logTab(1, LXML,
                  StringUtil::MakeXMLtag("AverageOverlap",
                                         averageoverlap,5,2));
    output.logTab(1, LXML,
                  StringUtil::MakeXMLtag("Overlapthreshold",
                                         overlapthreshold,5,2));
    output.logTab(1, LXML,
                  StringUtil::MakeXMLtag("AllowedGap",
                                         allowedgap,5,2));
    output.logTab(1, LXML,
                  StringUtil::MakeXMLtag("MinimumOverlap",
                                         minimumoverlap,5,2));
    output.logTab(0, LXML,"</RotationalOverlap>");

  }
} // namespace scala
