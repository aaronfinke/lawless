//
// batchgroup.cpp
//

#include "batchgroup.hh"
#include "util.hh"
#include "string_util.hh"

namespace scala {
  // ------------------------------------------------------------
  // allscales             scale model, for maximum primary scale smoothing range
  // analysiscontrols  for batchgroupwidth, <= 0.0 for individual
  // hkl_list          address store for batch number -> batch serial lookup
  // datasetIndex      < 0 for all datasets, see statistics
  Batchgroup::Batchgroup(const ScaleModel& allscales,
                         const AnalysisControls& analysiscontrols,
                         const hkl_unmerge_list& hkl_list,
                         const int& datasetIndex)
  {
    init(allscales, analysiscontrols, hkl_list, datasetIndex);
  }
  // ------------------------------------------------------------
  void Batchgroup::init(const ScaleModel& allscales,
                        const AnalysisControls& analysiscontrols,
                        const hkl_unmerge_list& hkl_list,
                        const int& datasetIndex)
  // allscales             scale model, for maximum primary scale smoothing range
  // analysiscontrols      for user specifications
  // hkl_list              for header information: a pointer to this is
  //                       stored for the batch number/serial lookup
  // datasetIndex          < 0 for all datasets
  {
    std::vector<Batch> batches = hkl_list.Batches();
    nbatches = batches.size();
    batchgroupwidth = analysiscontrols.Batchgroupwidth();
    datasetindex = datasetIndex;
    hkl_list_pointer = &hkl_list;

    groupindex.assign(nbatches, -2); // group number for each batch

    validPhi = true;  // false if any batch is invalid
    float phirangemin = +10000.;
    float phirangemax = -phirangemin;

    individualbatches = false;
    if (batchgroupwidth <= 0.0) {individualbatches = true;}
    //    std::cout << "batchgroupwidth "<<batchgroupwidth <<std::endl;

    // Get minimum and maximum phi range
    for (int ib=0;ib<nbatches;++ib) {
      if ((datasetIndex < 0) || (batches[ib].datasetindex() == datasetIndex)) {
        if (batches[ib].Accepted()) {
          if (batches[ib].ValidPhi()) {
            phirangemin = std::min(phirangemin, batches[ib].PhiRange());
            phirangemax = std::max(phirangemax, batches[ib].PhiRange());
            groupindex[ib] = -1;  // in dataset and valid
          } else {
            validPhi = false;
            groupindex[ib] = -1;  // in dataset and valid
          }
        }
      }
    }
    if (phirangemin < 0.00001) {
      validPhi = false;
    }

    differentPhiRanges = false;
    batchindex.clear();
    nbatchingroup.clear();
    std::vector<Run> runList = runlist();
    ngroupsinrun.assign(runList.size(), 0);

    //    std::cout << "valid " <<validPhi<<" "<< individualbatches<<std::endl;
    if (individualbatches) {
      // explicit individual batches
      individualbatches = true;
      numberinGroup = 1;
      // batches will be individually analysed
      int igroup = 0;
      for (int ib=0;ib<nbatches;++ib) {
        if (groupindex[ib] == -1) { // valid and accepted batches
          groupindex[ib] = igroup;      // individual batches
          batchindex.push_back(ib); //
          nbatchingroup.push_back(1); // number of batches in group
          int irun = batches[ib].RunIndex();
          grouprun.push_back(irun);
          igroup++;
        }
      }
      for (size_t irun=0; irun<runList.size(); irun++) {
        if (runList[irun].DatasetIndex() == datasetIndex) {
          // count of accepted batches
          ngroupsinrun[irun] = runList[irun].BatchList(true).size();
        }
      }
    } else {
      if (validPhi) {
        // valid Phi values
        if (!Close<float, float>(phirangemin, phirangemax, 0.00001)) {
          // different phi ranges, should we do something about this?
          differentPhiRanges = true;
        }
        // determine a suitable group width
        numberinGroup = std::max(1, Nint(batchgroupwidth/phirangemax));
      } else {
        // not valid Phi ranges, but still group batches
        const int DEFAULT_NINGRP = 5;
        numberinGroup = std::max(1, std::min(DEFAULT_NINGRP, nbatches/10));
      }
      //      std::cout <<"batchgroup, width, maxphi, numberingroup "<<
      //        batchgroupwidth<<" "<<phirangemax<<" "<<numberinGroup<<std::endl; //^

      for (size_t irun=0; irun<runList.size(); irun++) {
        // Runs in dataset
        if (runList[irun].DatasetIndex() == datasetIndex) {
          // accepted batches
          std::vector<int> batchnums = runList[irun].BatchList(true);
          //^^
          //      std::cout <<"\n!!! Run "<<irun<<"\n";
          //      for (size_t i=0; i<batchnums.size(); i++) {
          //        std::cout << "  Batch: "<<i<<" "<<batchnums[i]<<"\n";
          //      } //^-
          // make groups, store count of groups in run
          ngroupsinrun[irun]  = makeGroups(batchnums, irun);
        }
      }
      // groupindex and batchindex arrays are now set up
      //^^
      //      std::cout << "Group number for each batch serial, NB "
      //                << groupindex.size()<<"\n";
      //      for (size_t ib=0; ib<groupindex.size(); ib++) {
      //        int ibatch = hkl_list_pointer->batch(ib).num();
      //        std::cout << ib<<" "<<ibatch <<" "<<groupindex[ib]<<"\n";
      //      }
      //      std::cout << "1st batch serial in each group\n";
      //      for (size_t j=0; j<batchindex.size(); j++) {
      //        int ibatch = hkl_list_pointer->batch(batchindex[j]).num();
      //        int ibser = hkl_list_pointer->batch_serial(ibatch);
      //        std::cout <<j<<" "<< batchindex[j]<<" "<<ibatch<<" "<<ibser
      //                  <<" run "<<grouprun[j]<<" nbg "<<nbatchingroup[j]<<"\n";
      //      } //^-

    } // not individual

    // Smoothing
    double spacing = 0.0;
    for (int k=0; k<hkl_list.num_runs(); k++) {
      spacing = std::max(spacing, allscales.primary_scale(k).Spacing());
    }

    numgroupsmooth = 5;
    double smoothwidth = analysiscontrols.SmoothStatisticsRange();
    if (smoothwidth <= 0.0) {
      smoothwidth = spacing;
    }
    if (phirangemax > 0.0 && smoothwidth > 0.0) {
      int bgw = 1.0;
      if (batchgroupwidth > 0.0) {
        bgw = batchgroupwidth;
      }
      numgroupsmooth = (Nint(smoothwidth/bgw)/2)*2 + 1; // force to be odd
    }
    //    std::cout << "numgroupsmooth "<< numgroupsmooth<<std::endl;
  }
  // ------------------------------------------------------------
  // Dummy for testing, just set ngroups
  //  Assume batch numbers are 1 to Ngroups
  Batchgroup::Batchgroup(const int& Ngroups)
  {
    ngroup = Ngroups;
    nbatches = -1;
  }
  // ------------------------------------------------------------
  // return group number for batch number, < 0 if not set
  int Batchgroup::batchgroup(const int& batchnumber) const
  //   within dataset!!!
  {
    // batch serial number for Batch
    int jgroup = std::max(0, batchnumber-1);
    if (hkl_list_pointer != 0) {
      int ib = hkl_list_pointer->batch_serial(batchnumber);
      jgroup = groupindex[ib];
    }
    return jgroup;
  }
  // ------------------------------------------------------------
  // return central batch serial for group
  int Batchgroup::batchserial(const int& jgroup) const
  {
    int serial = jgroup;
    if (hkl_list_pointer != 0) {
      serial = std::min(batchindex.at(jgroup) +
                        (nbatchingroup.at(jgroup)-1)/2,
                        nbatches-1);
    }
    //    std::cout <<"*&* batchserial "<<serial<<"\n";
    return serial;
  }
  // ------------------------------------------------------------
  // Return batch list
  std::vector<Batch> Batchgroup::batches() const {
    if (hkl_list_pointer != 0) {
      return hkl_list_pointer->Batches();
    }
    return std::vector<Batch>(ngroup, Batch());
  }
  // ------------------------------------------------------------
  // format for printing
  std::string Batchgroup::format() const
  {
    std::string s;
    if (individualbatches) {
      s = "Analysis by Batch uses separate batches";
      return s;
    }
    s = "Analysis by Batch is in groups of batches, group width "+
      StringUtil::ftos(batchgroupwidth,7,2)+" degrees, "+
      StringUtil::itos(numberinGroup,3)+" batches/group";
    return s;
  }
  // ------------------------------------------------------------
  // format for XML
  std::string Batchgroup::formatXML() const
  {
    std::string s = StringUtil::MakeXMLtag("BatchGroupMessage", format());
    s += "\n"+StringUtil::MakeXMLtag("BatchGroupWidth",batchgroupwidth,5,2);
    s += "\n"+StringUtil::MakeXMLtag("BatchGroupNumber",numberinGroup);
    return s;
  }
  //--------------------------------------------------------------
  int Batchgroup::makeGroups(const std::vector<int>& batchnumberlist,
                             const int& irun)
  // Given a list of contiguous batch numbers for run irun,
  //   assign them to groups
  //   returns number of groups
  {
    int nb = batchnumberlist.size();
    int ngroups = 1;
    int lastgroupsize = 1;
    if (numberinGroup > 0) {
      ngroups = nb/numberinGroup;
      lastgroupsize = nb%numberinGroup;  // number in last group if > 0
    }

    if (lastgroupsize > 0) {ngroups++;}
    // If last group is "too small", put into previous group
    const float MINSIZE = 0.25;  // allowed fraction of group size
    // don't bother to check if numberinGroup is small
    const int MINNINGROUP = 3;
    if (lastgroupsize > 0 && numberinGroup > MINNINGROUP) {
      if (lastgroupsize < Nint(MINSIZE*numberinGroup)) {
        ngroups = std::max(1, ngroups-1);
        lastgroupsize = nb - (ngroups - 1)*numberinGroup;
      }}
    //    std::cout << "Batchgroup::makeGroups NiG, nb, ngroups, lastgroupsize, irun " <<
    //      numberinGroup<<" "<<
    //      nb<<" "<<ngroups <<" "<< lastgroupsize <<" "<<irun<<"\n";

    // Loop batches in list
    int lastgroup = -1;
    int k = 0;
    int igroup = -1;
    // number of groups so far
    int jgroup0 = batchindex.size();
    int ib0 = -1;  // first batch serial in group
    int ningroup = 0;  // number of batches in each group
    int ib;
    int ng = 0;
    for (int i=0; i<nb; ++i) {  // loop batches in group
      // batch serial
      ib = hkl_list_pointer->batch_serial(batchnumberlist[k]);
      if (ib0 < 0) {ib0 = ib;}
      if (groupindex[ib] == -1) { // accepted batch
        // group assignment within this list
        igroup = k/numberinGroup;
        if (igroup > ngroups-1) {
          igroup = ngroups-1;  // put excess into last group
        }
        // global group number
        int jgroup = igroup + jgroup0;
        groupindex[ib] = jgroup;
        //      std::cout << "batch ib, jgroup "<<ib<<" "<<jgroup<<"\n";
        if (lastgroup >= 0 && igroup != lastgroup) {
          batchindex.push_back(ib0);
          nbatchingroup.push_back(ib-ib0); // number of batches in group
          //      std::cout <<" new group, ib0, ib, irun, bisize "<<
          //        ib0<<" "<<ib<<" "<<irun<<" "<<batchindex.size()<<"\n";
          ib0 = ib;
          grouprun.push_back(irun);
          ng++;
        }
        lastgroup = igroup;
        k++;
      }
    }
    ib = hkl_list_pointer->batch_serial(batchnumberlist.back());
    if (ib >= ib0) {
      batchindex.push_back(ib0);
      // number of batches in group
      nbatchingroup.push_back(std::max(1,(ib-ib0)));
      grouprun.push_back(irun);
      ng++;
      //      std::cout <<" last group, ib, ib0, irun, bisize "<<
      //        ib<<" "<<ib0<<" "<<irun<<" "<<batchindex.size()<<"\n";
    }
    ASSERT (ng == ngroups);

    return ngroups;
  }
  //--------------------------------------------------------------
  // List of groups in run
  std::vector<int> Batchgroup::groupsinrun(const int& irun) const
  {
    std::vector<int> grouplist;
    for (size_t jgroup=0; jgroup<grouprun.size(); jgroup++) {
      if (grouprun[jgroup] == irun) {
        grouplist.push_back(jgroup);
      }
    }
    return grouplist;
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  Xbreaks::Xbreaks(const std::vector<Batch>& batches,
                   const int& datasetIndex)
  // Find all breaks in batch number list, for X axis in plots
  // datasetIndex < 0 for all datasets
  {
    init(batches, datasetIndex);
  }
  // ------------------------------------------------------------
  Xbreaks::Xbreaks(const Batchgroup& batchgroup,
                   const int& datasetIndex)
  // datasetIndex = -1 for all datasets
  {
    std::vector<Batch> batches = batchgroup.batches();
    init(batches, datasetIndex);
  }
  //--------------------------------------------------------------
  void Xbreaks::init(const std::vector<Batch>& batches,
                     const int& datasetIndex)
  // Find all breaks in batch number list, for X axis in plots
  // datasetIndex < 0 for all datasets
  {
    validbatchnumbers.clear();
    int lastbatchnum = -1;
    int mingap = 2; // don't break with fewer than mingap missing
    for (size_t i=0;i<batches.size();++i) {  // even batches that have no reflections
      // ... but not rejected batches
      if (((datasetIndex < 0) || (batches[i].datasetindex() == datasetIndex))
          && batches[i].Accepted()) {
        if (lastbatchnum >= 0) {
          int gap = batches[i].num() - lastbatchnum;
          if (gap > mingap) { // we have a break
            breaks.push_back(Range(lastbatchnum, batches[i].num()));
          }
        }
        lastbatchnum = batches[i].num();
        validbatchnumbers.update(batches[i].num());
      }
    }
  }
  //--------------------------------------------------------------
  int GetSymbolSizeforNpoints(const int& npoints)
  //  symbol size  = 0 if "too many" points
  {
    // no symbol if too many points
    int symbolsize = -1; // default symbol size
    const int MAXGROUPSFORSYMBOL = 90;
    if (npoints > MAXGROUPSFORSYMBOL) {symbolsize = 0;}
    return symbolsize;
  }
}
