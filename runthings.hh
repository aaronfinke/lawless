// runthings.hh
//
// Things concerning runs & their unique batch numbers

#ifndef RUNTHINGS_HEADER
#define RUNTHINGS_HEADER

#include <vector>
#include "hkl_datatypes.hh"
#include "range.hh"

namespace scala {

  //--------------------------------------------------------------
  class Run
  // A run is a group of contiguous batches belonging to the same dataset
  {
  public:
    enum FullsAndPartials {FULLSANDPARTIALS, ONLYFULLS, ONLYPARTIALS, FEWFULLS, FEWPARTIALS};

    Run();
    Run(const int& DatasetIndex);
    void clear();
    void clearCounts();  // clear reflection/observation counts

    void AddBatch(const int& Batch, const bool& accept=true);

    void SortList();  // call this after last AddBatch
    
    // true if Batch is in list
    bool IsInList(const int& Batch) const;
    // Return batch number list:
    //  if Accepted == true, only return accepted batches
    std::vector<int> BatchList(const bool& Accepted=false) const;
    int DatasetIndex() const {return dataset_index;}
    // return minimum & maximum batch number
    std::pair<int,int> BatchRange() const;
    // return first batch number
    int Batch0() const;
    //  First batch serial number
    void SetBatchSerial0(const int& BatchSerial)
      {batchserial0 = BatchSerial;} // set
    int  BatchSerial0() const {return batchserial0;} // get

    // Number of batches
    int Nbatches() const {return batch_number_list.size();}

    // offset all batch numbers & increment (store) total offset
    void OffsetBatchNumbers(const int& offset);
    int& BatchNumberOffset() {return batch_number_offset;} // set
    int  BatchNumberOffset() const {return batch_number_offset;} // get
    // Run number
    int& RunNumber() {return runnumber;} // set
    int  RunNumber() const {return runnumber;} // get
    // Filenumber
    int& FileNumber() {return file_number;} // set
    int  FileNumber() const {return file_number;} // get
    // Phi start & finish, Phi1, Phi2
    Range& PhiRange() {return phirange;}  // set
    Range  PhiRange() const {return phirange;}  // get
    // Time start & finish, Time1, Time2
    Range& TimeRange() {return timerange;}  // set
    Range  TimeRange() const {return timerange;}  // get
    // Time is invalid if it is just a copy of Phi
    bool& ValidTime() {return validtime;}  // set
    bool ValidTime() const {return validtime;}  // get
    void SetValidOrientation(const bool& validOrientation)
      {validorientation = validOrientation;}
    bool ValidOrientation() const {return validorientation;}
    void StoreSpindleToPrincipleAxis(const std::string& Spindletoprincipleaxis)
    {spindletoprincipleaxis = Spindletoprincipleaxis;}
    std::string SpindleToPrincipleAxis() const
    {return spindletoprincipleaxis;}

    static int MaximumBatchNumber() {return MaxBatchNumber;}
    
    std::string formatPrint(const std::vector<Xdataset>& datasets) const;
    std::string formatPrintBrief(const std::vector<Xdataset>& datasets) const;

    int& Nfulls() {return nfulls;}   //!< set number of fully recorded observations
    int& Npartials() {return npartials;}  //!< set number of partial observations

    int Nfulls() const {return nfulls;}   //!< get number of fully recorded observations
    int Npartials() const {return npartials;}  //!< get number of partial observations

    //! Set FullsAndPartials flag based on numbers
    void SetFullsAndPartials();
    //! Return FullsAndPartials flag
    FullsAndPartials fullsAndPartials()const {return fullsandpartials;}

    // format for saving essentials
    std::string FormatSave() const;

    //! set resolution range
    void StoreResoRange(const ResoRange& Resrange);
    //! return resolution range
    ResoRange GetResoRange() const {return resrange;}
    //! return true if there is a resolution cutoff
    bool IsResoRange() const {return resrangeset;}
    //! return true if observation is in range
    // note that the batch number is not used at present but may be in future
    bool InResoRange(const Rtype& invresolsq, const int& batch);

  private:
    int runnumber;   // run number, an arbitrary integer (typically 1,2,3...)
    static const int MaxBatchNumber;
    int dataset_index; // index into dataset list
    mutable std::vector<int> batch_number_list;
    std::vector<bool> batch_accepted;
    int batch_number_offset;
    int batchserial0;
    int file_number;
    Range phirange;
    Range timerange;
    bool validtime;
    bool validorientation;
    std::string spindletoprincipleaxis; // name of reciprocal axis closest to spindle
    int nfulls;     // number of fully recorded observations
    int npartials;  // number of partial observations
    FullsAndPartials fullsandpartials;
    ResoRange resrange; // optional resolution limit for this run
    bool resrangeset;   // true if there is a run resolution limit
  };
  //--------------------------------------------------------------
  class RunRange 
  {
  public:
    RunRange() : Offset(0) {}
    RunRange(const int& Low, const int& High)
      : LowBatchNumber(Low), HighBatchNumber(High), Offset(0)  {}
    RunRange(const std::pair<int,int>& LowHigh)
      :	LowBatchNumber(LowHigh.first), HighBatchNumber(LowHigh.second), Offset(0) {}
    
    bool Encloses(const RunRange& test) const;
    bool Encloses(const int& testN) const;
    
    void IncrementOffset(const int& offset);
    
    // Range of batch numbers
    int LowBatchNumber;
    int HighBatchNumber;
    int Offset;   // offset, initially 0
  };
  //--------------------------------------------------------------
  // Determine suitable batch number offsets for each run in testruns
  // such that they do not clash with any runs in refruns nor with other
  // runs in testruns
  // Returns vector of offsets, or empty vector if fails
  // Note the offsets which are tried are all multiples of 1000
  std::vector<int> CompareRunRanges(const std::vector<Run>& refruns,
				    const std::vector<Run>& testruns);
  //--------------------------------------------------------------
  // Get run index number for a given run number, from a list of runs
  // returns -1 if not found
  int FindRunIndex(const int& runNumber, const std::vector<Run>& runList);
} // namespace scala
#endif
