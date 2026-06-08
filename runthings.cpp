// runthings.cpp

#include <algorithm>

#include "runthings.hh"
#include "string_util.hh"
#include "dataset.hh"

namespace scala {
  //--------------------------------------------------------------
  Run::Run() {clear();}
  //--------------------------------------------------------------
  Run::Run(const int& DatasetIndex, const int& DatasetID,
           const PxdName& Pxdname)
  {
    clear();
    dataset_index= DatasetIndex;
    datasetID = DatasetID;
    pxdname = Pxdname;
  }
  //--------------------------------------------------------------
  void Run::clear()
  {
    batch_number_list.clear();
    runnumber = 0;
    dataset_index = -1;
    datasetID = -1;
    pxdname = PxdName();
    batch_number_offset = 0;
    batchserial0 = 0;
    file_number = 0;
    phirange = Range();
    timerange = Range();
    validtime = false;
    validorientation = false;
    nfulls = 0;
    npartials = 0;
    fullsandpartials = FULLSANDPARTIALS;
    resrangeset = false;
    resrange = ResoRange();
    use_ = true;
  }
  //--------------------------------------------------------------
  void Run::clearCounts()
  {
    nfulls = 0;
    npartials = 0;
    fullsandpartials = FULLSANDPARTIALS;
  }
  //--------------------------------------------------------------
  void Run::AddBatch(const int& Batch, const bool& accept)
  {
    batch_number_list.push_back(std::pair<int,bool>(Batch,accept));
  }
  //--------------------------------------------------------------
  int Run::indexInList(const int& Batchnum) const
  // returns index into list if found, else -1
  {
    int jb = -1;
    for (size_t ib=0;ib<batch_number_list.size();++ib) {
      if (batch_number_list[ib].first == Batchnum) {
        jb = ib;
        break;
      }
    }
    return jb;
  }
  //--------------------------------------------------------------
  bool Run::IsInList(const int& Batchnum) const
  // true if Batch is in list
  {
    return (indexInList(Batchnum) >= 0);
  }
  //--------------------------------------------------------------
  std::vector<int> Run::BatchList(const bool& Accepted) const
  // Return batch number list:
  //  if Accepted == true, only return accepted batches
  {
    std::vector<int> bn;
    for (size_t ib=0;ib<batch_number_list.size();++ib) {
      if (Accepted) {
        if (batch_number_list.at(ib).second) bn.push_back(batch_number_list[ib].first);
      } else {
        bn.push_back(batch_number_list[ib].first);
      }
    }
    return bn;
  }
  //--------------------------------------------------------------
  void Run::SetBatchAccept(const int& Batchnum, const bool& Accepted)
  // Set batch accept flag
  {
    int jb = indexInList(Batchnum);
    if (jb < 0) {return;} // not found
    batch_number_list[jb].second = Accepted;
  }
  //--------------------------------------------------------------
  bool Run::IsBatchAccepted(const int& Batchnum) const
  // return batch accepted flag
  {
    int jb = indexInList(Batchnum);
    if (jb < 0) {return false;} // not found
    return batch_number_list[jb].second;
  }
  //--------------------------------------------------------------
  void Run::OffsetBatchNumbers(const int& offset)
  {
    for (size_t i=0;i<batch_number_list.size();i++) {
      batch_number_list[i].first += offset;
    }
    batch_number_offset += offset;
  }
  //--------------------------------------------------------------
  std::string Run::formatPrint(const std::vector<Dataset>& datasets) const
  {
    std::string s =
      FormatOutput::logTabPrintf(0,
              "\nRun number: %3d  Dataset: %3d %s consists of batches:-",
                         runnumber, dataset_index+1,
                                 datasets[dataset_index].pxdname().format().c_str());
    const int nperline = 15;
    for (size_t i=0;i<batch_number_list.size();i++)  {
      if (batch_number_list[i].second) {
        if (i%nperline == 0) s += FormatOutput::logTabPrintf(0,"\n");
        s += FormatOutput::logTabPrintf(0," %6d", batch_number_list[i].first);
      }
    }
    if (latnum > 0) {
      s += " Lattice number " + StringUtil::itos(latnum, 2);
    }
    if (!use_) {
      s += " [Not Used]";
    }

    s += FormatOutput::logTabPrintf(0,"\n");
    if (resrangeset) {
      s +=  FormatOutput::logTabPrintf(0,"\n   Resolution range for run: %8.2f    %8.2f\n",
                          resrange.ResLow(), resrange.ResHigh());
    }
    s += FormatOutput::logTabPrintf(3,
                        "Phi range: %8.2f to %8.2f   Time range: %8.2f to %8.2f\n",
                                    phirange.first(),phirange.last(),
                                    timerange.first(),timerange.last());
    if (validorientation) {
      s += FormatOutput::logTab(1,
                    "Closest reciprocal axis to spindle: "+
                    spindletoprincipleaxis);
    }
    return s;
  }
  //--------------------------------------------------------------
  //! return list of batch number ranges
  std::vector<scala::IntRange> Run::BatchNumberRanges() const
  {
    std::vector<scala::IntRange> ranges;
    int i1 = -1;
    int i2;
    for (size_t ib=0;ib<batch_number_list.size();++ib) { // loop batches
      if (batch_number_list[ib].second) {
        if (i1 < 0) { // first
          i1 = batch_number_list[ib].first;
          i2 = i1;
        } else { // not first
          if ((batch_number_list[ib].first - i2) == +1) {
            // contiguous batch numbers
            i2 = batch_number_list[ib].first;
          } else {
            // break in batch number series, batches i1 to i2
            ranges.push_back(IntRange(i1,i2));
            i1 = batch_number_list[ib].first;
            i2 = i1;
          }
        }
      }
    } // end loop batches
    ranges.push_back(IntRange(i1,i2));
    return ranges;
  }
  //--------------------------------------------------------------
  std::string Run::formatPrintBrief() const
  {
    // Only list accepted batches
    std::vector<scala::IntRange> batchnumberranges = BatchNumberRanges();

    std::string s;

    if (batchnumberranges.size() > 0) {
      for (size_t i=0;i<batchnumberranges.size();++i) {
        if (s != "") s += ", ";
        s += StringUtil::Strip(clipper::String(batchnumberranges[i].min()))+" - "+
          StringUtil::Strip(clipper::String(batchnumberranges[i].max()));
      }


      /*
        int firstbatch = -1;;
        int lastbatch = -1;
        for (size_t i=0;i<batch_number_list.size();++i) {
        if (batch_number_list[i].second) {
        firstbatch = batch_number_list[i].first;
        break;
        }
        }
        if (firstbatch >= 0) {
        for (int i=int(batch_number_list.size())-1;i>=0;i--) {
        if (batch_number_list[i].second) {
        lastbatch = batch_number_list[i].first;
        break;
        }
        }
      */

      s = FormatOutput::logTab(2,
                               "Run number: "+StringUtil::itos(runnumber, 3)+
                               " consists of batches "+s, false);

      if (latnum > 0) {
        s += " Lattice number " + StringUtil::itos(latnum, 2);
      }
      if (!use_) {
        s += " [Not Used]";
      }
      s += "\n";
      if (resrangeset) {
        s +=  FormatOutput::logTabPrintf(3,"Resolution range for run: %8.2f    %8.2f\n",
                                         resrange.ResLow(), resrange.ResHigh());
      }
      s += FormatOutput::logTabPrintf(3,
                                      "Phi range: %8.2f to %8.2f   Time range: %8.2f to %8.2f\n",
                                      phirange.first(),phirange.last(),
                                      timerange.first(),timerange.last());
      if (validorientation) {
        s += FormatOutput::logTab(3,
                                  "Closest reciprocal axis to spindle: "+
                                  spindletoprincipleaxis);
      }
    } else {
      s = "No accepted batches in run";
    }
    return s;
  }
  //--------------------------------------------------------------
  const int Run::MaxBatchNumber = 9999999; // should be OK in 32bit float
  //--------------------------------------------------------------
  void Run::SortList()  // call this after last AddBatch
  {
    // sort list, sorted on batch number (.first), carries along accepted flag
    std::sort(batch_number_list.begin(), batch_number_list.end());
  }
  //--------------------------------------------------------------
  std::pair<int,int> Run::BatchRange() const
  // return minimum & maximum batch number
  {
    return std::pair<int,int>(batch_number_list[0].first, batch_number_list.back().first);
  }
  //--------------------------------------------------------------
  // return first batch number
  int Run::Batch0() const {
    return batch_number_list[0].first;
  }
  //--------------------------------------------------------------
  //! Set FullsAndPartials flag based on numbers
  void Run::SetFullsAndPartials()
  {
    fullsandpartials = FULLSANDPARTIALS;
    if (nfulls == 0 && npartials == 0) return;
    if (nfulls == 0) fullsandpartials = ONLYPARTIALS;
    else if (npartials == 0) fullsandpartials = ONLYFULLS;
    else {
      // Check for very few fulls or partials: if so combine them
      // Criteria for only a "few" fulls or partials
      //  "Few" if either
      //     1. fraction (f = full/partial or 1/f)  .lt. FMINFEW
      //     2. number .lt. MINFEW   and f .lt. 0.5
      int combine = 0;
      const float FMINFEW = 0.01; // 3%
      const int MINFEW  = 100;
      float f = float(nfulls)/float(npartials);
      if (f < FMINFEW) combine = -1;  // few fulls
      else if (1.0/f < FMINFEW) combine = +1;  // few partials
      if (nfulls < MINFEW && f < 0.5) combine = -1;  // few fulls
      if (npartials < MINFEW && 1.0/f < 0.5) combine = +1;  // few partials
      if (combine < 0) fullsandpartials = FEWFULLS;
      else if (combine > 0) fullsandpartials = FEWPARTIALS;
    }
  }
  //--------------------------------------------------------------
  // format for saving essentials
  std::string Run::FormatSave() const
  {
    std::string dump = "Run V1 {\n";
    int nb = batch_number_list.size();
    std::vector<int> numbers(nb);
    for (size_t i=0;i<batch_number_list.size();++i) {
      numbers[i] = batch_number_list[i].first;
    }
    dump += " Batch_number_list "+
      clipper::String(nb)+"\n"+
      StringUtil::FormatSaveVector(numbers);
    return dump+"}\n";
  }
  //--------------------------------------------------------------
  //! set resolution range
  void Run::StoreResoRange(const ResoRange& Resrange)
  {
    resrange = Resrange;
    resrangeset = true;
  }
  //--------------------------------------------------------------
  //! return true if observation is in range
  // note that the batch number is not used at present but may be in future
  bool Run::InResoRange(const Rtype& invresolsq, const int& batch)
  {
    if (!resrangeset) return true;
    if (resrange.tbin(invresolsq) < 0) {return false;}
    return true;
  }
  //--------------------------------------------------------------
  RunRange::RunRange(const int& Low, const int& High)
    : LowBatchNumber(Low), HighBatchNumber(High), Offset(0)
  {
    if (LowBatchNumber > HighBatchNumber) {
      std::swap(LowBatchNumber, HighBatchNumber);
    }
  }
  //--------------------------------------------------------------
  RunRange::RunRange(const std::pair<int,int>& LowHigh)
    :   LowBatchNumber(LowHigh.first), HighBatchNumber(LowHigh.second), Offset(0)
  {
    if (LowBatchNumber > HighBatchNumber) {
      std::swap(LowBatchNumber, HighBatchNumber);
    }
  }
  //--------------------------------------------------------------
  bool RunRange::Adjacent(const RunRange& test) const
  // returns true if test run and this one are adjacent
  // test both ways round
  {
    return ((LowBatchNumber == test.HighBatchNumber+1) ||
            (test.LowBatchNumber == HighBatchNumber+1));
  }
  //--------------------------------------------------------------
  bool RunRange::Encloses(const RunRange& test) const
  // returns true if test run overlaps with this one, or are adjacent
  // test both ways round
  {
    return (Encloses(test.LowBatchNumber) || Encloses(test.HighBatchNumber) ||
            test.Encloses(LowBatchNumber) || test.Encloses(HighBatchNumber));
  }
  //--------------------------------------------------------------
  bool RunRange::Encloses(const int& testN) const
  // returns true if testN is inside this one, or are adjacent
  {
    return (testN >= LowBatchNumber && testN <= HighBatchNumber+1);
  }
  //--------------------------------------------------------------
  void RunRange::IncrementOffset(const int& offset)
  {
    LowBatchNumber += offset;
    HighBatchNumber += offset;
    Offset += offset;
  }
  //--------------------------------------------------------------
  std::vector<int> CompareRunRanges(const std::vector<Run>& refruns,
                                    const std::vector<Run>& testruns)
  // Determine suitable batch number offsets for each run in testruns
  // such that they do not clash with any runs in refruns nor with other
  // runs in testruns
  // Returns vector of offsets, or empty vector if fails
  {
    int OffsetIncrement = 1000;
    int MaxOffset = scala::Run::MaximumBatchNumber() - 2*OffsetIncrement;

    std::vector<RunRange> ref_rr(refruns.size());
    for (size_t i=0;i<refruns.size();i++) {
      ref_rr[i] = RunRange(refruns[i].BatchRange());
    }
    std::vector<RunRange> test_rr(testruns.size());
    for (size_t i=0;i<testruns.size();i++) {
      test_rr[i] = RunRange(testruns[i].BatchRange());
    }

    std::vector<int> offsets(test_rr.size(), 0);
    bool overflow = false; // true if we run off end of MaxBatch

    for (size_t it=0;it<test_rr.size();it++) { // Loop each test run
      bool found = false;
      while (test_rr[it].Offset <= MaxOffset) {
        bool OK = true;
        for (size_t ir=0;ir<ref_rr.size();ir++) {
          // Loop each reference run
          if (ref_rr[ir].Encloses(test_rr[it])) {
            // but it might be OK if they are adjacent before any offset
            if (!((test_rr[it].Offset == 0) && ref_rr[ir].Adjacent(test_rr[it]))) {
              OK = false;
              break;
            }
          }
        } // end loop reference runs
        if (OK) {
          // this test run does not overlap with any reference run
          // Test it against the other test runs
          for (size_t it2=0;it2<test_rr.size();it2++) {
            if (it2 != it) {
              if (test_rr[it2].Encloses(test_rr[it])) {
                OK = false;
                break;
              }
            }
          }
        }
        if (OK) {
          // No overlap, offset OK
          found = true;
          break;
        }
        test_rr[it].IncrementOffset(OffsetIncrement);
      }  // test next offset
      if (found) {
        // We have found a suitable offset for this testrun
        offsets[it] = test_rr[it].Offset;
      }
      else
        {overflow = true;}
    } // end loop test_rr
    if (overflow) offsets.clear();  // clear vector to indicate failure
    return offsets;
  }
  //--------------------------------------------------------------
  // Get run index number for a given run number, from a list of runs
  // returns -1 if not found
  int FindRunIndex(const int& runNumber, const std::vector<Run>& runList)
  {
    for (size_t irun=0;irun<runList.size();irun++) {
      if (runNumber == runList[irun].RunNumber()) {
        return irun;
      }
    }
    return -1;
  }
  //--------------------------------------------------------------
}
