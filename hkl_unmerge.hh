// hkl_unmerge.hh
// Phil Evans August 2003-2004 etc
//
// Classes for unmerged hkl lists
//  at present just for the items used by Mosflm & Scala

// These data structures are primarily directed at Scala so are
// in its namespace

#ifndef HKL_UNMERGE_HEADER
#define HKL_UNMERGE_HEADER

#include <map>
#include "hkl_datatypes.hh"
#include "hkl_controls.hh"
#include "controls.hh"
#include "hkl_symmetry.hh"
#include "icering.hh"
#include "observationflags.hh"
#include "range.hh"
#include "hash.hh"
#include "runthings.hh"

namespace scala {
  //==============================================================
  class data_flags
  //! Data column flags to indicate which columns are present in the file
  // True if column is present
  {
  public:
    data_flags();  //!< Compulsory flags set true, others false
    void print() const; // for debugging
    
    bool is_h, is_k, is_l, is_misym, is_batch,
      is_I, is_sigI, is_Ipr, is_sigIpr, is_fractioncalc,
      is_Xdet, is_Ydet, is_Rot, is_Width, is_LP, is_Mpart,
      is_ObsFlag, is_BgPkRatio, is_scale, is_sigscale, is_time;
  };  // data_flags
  //===================================================================
  // *******************   observation, reflections etc
  class SelectI 
  {
    //!< A static class to select profile-fitted, summeation integrated, or combined intensity
    //

  public:
    SelectI(){}
    //  Algorithm for selection as in Scala
    //   INTEGRATED   integrated intensity I                  (SelectIcolFlag=0)
    //   PROFILE      profile-fitted intensity IPR            (SelectIcolFlag=-1)
    //   COMBINE      weighted mean after combining partials  (SelectIcolFlag=+1)
    //       combine integrated intensity Iint and
    //       profile-fitted intensity Ipr as
    //         I = w * Ipr  +  (1-w) * Iint
    //       where w = 1/(1+(Iint/Imid)**IPOWER)
    //       Imid (= IcolFlag)  is the point of equal weight, and
    //       should be about the mean intensity.
    //  Ipower default = 3

    // return true if COMBINE
    static bool Combine() {return (selecticolflag > 0);}
    // If IcolFlag > 0, set selecticolflag and set imid = Imid, ie COMBINE
    // If IcolFlag = 0, select Iint, < 0 select Ipr
    static void SetIcolFlag(const int& IcolFlag, const double& Imid, const int Ipower=3);
    // Set SelectIcolFlag
    static int& SelectIcolFlag() {return selecticolflag;}
    // store imid = overall <I> if needed and not set
    static void SetAverageIntensity(const double& meanI);
    // store imid = overall <I> if needed
    static void ResetAverageIntensity(const double& meanI);

    static IsigI GetCombinedI(const Rtype& Iraw, const Rtype& Ic, const Rtype& varIc,
			      const Rtype& Ipr, const Rtype& varIpr);

    static IsigI GetCombinedI(const Rtype& Iraw, const IsigI& Isc, const IsigI& Ispr);

    //! set true if we have a second intensity Ipr stored
    static void SetIprPresent(const bool& Isiprpresent) {iprpresent = Isiprpresent;}
    //! return true if we have a second intensity Ipr stored
    static bool IsIprPresent() {return iprpresent;}

    static std::string format();

  private:
    static int selecticolflag;
    static int ipowercomb;
    static double imid;
    static bool iprpresent; // true if we have a second intensity Ipr stored
  //===================================================================
  };  // SelectI

// PartFlagSwitch
//   EMPTY  nothing stored
//   FULL   fully recorded reflection (1 part)
//   COMPLETE_CHECKED  complete, all Mpart flags consistent
//   COMPLETE passed total fraction checks
//   SCALE passed checks to be scaled
//   INCOMPLETE
//   EXCLUDED
 enum PartFlagSwitch
   {EMPTY, FULL, COMPLETE_CHECKED, COMPLETE, SCALE, INCOMPLETE, EXCLUDED};
  //===================================================================
  //! An observation of an intensity: may be a partial
  class observation_part
  // hkl is the reduced indices in the "current" spacegroup
  // run is set by call from hkl_unmerge_list::organise
  {
  public:

    // constructors
    observation_part();  // don't use 
    observation_part(const Hkl& hkl_in,
                     const int& isym_in, const int& batch_in,
                     const Rtype& I_in, const Rtype& sigI_in,
                     const Rtype& I_pr, const Rtype& sigIpr_in,
                     const Rtype& Xdet_in, const Rtype& Ydet_in,
                     const Rtype& phi_in, const Rtype& time_in,
                     const Rtype& fraction_calc_in, const Rtype& width_in,
                     const Rtype& LP_in,
                     const int& Npart_in, const int& Ipart_in,
                     const ObservationFlag& ObsFlag_in);
    // Npart is number of parts, from input (MPART)
    //       = 1 for full, = -1 unknown but partial
    // Ipart serial number in parts, from input (MPART)
    //       = 1 for full or unknown

    void set_run(const int& run) {run_ = run;}  //!< set by call from hkl_unmerge_list::organise


    inline  IsigI I_sigI() const {return IsigI(I_, sigI_);}
    inline  Rtype Ic() const {return I_;}       // actual I column
    inline  Rtype sigIc() const {return sigI_;} 
    inline  IsigI I_sigIpr() const {return IsigI(Ipr_, sigIpr_);}
    inline  Rtype Ipr() const {return Ipr_;}       //  Ipr column
    inline  Rtype sigIpr() const {return sigIpr_;} 

    void StoreIsigI(const IsigI& Is) {I_ = Is.I(); sigI_ = Is.sigI();}
    void StoreIsigIpr(const IsigI& Ipr) {Ipr_ = Ipr.I(); sigIpr_ = Ipr.sigI();}
    void StoreLP(const Rtype& LP) {LP_ = LP;}

    // selected I column
    // selIcolFlag = 0 Ic, = -1 Ipr
    inline  Rtype Ic(const int& selIcolFlag) const
    {return (selIcolFlag == 0) ? I_ : Ipr_;}
    inline  Rtype sigIc(const int& selIcolFlag) const
    {return (selIcolFlag == 0) ? sigI_ : sigIpr_;}


    inline  Rtype Xdet() const {return Xdet_;}
    inline  Rtype Ydet() const {return Ydet_;}
    inline  Rtype phi()  const {return phi_;}
    inline  Rtype time() const {return time_;}
    inline  Rtype fraction_calc() const {return fraction_calc_;}
    inline  Rtype width() const {return width_;}
    inline  Rtype LP() const {return LP_;}   // divide by this to get raw intensity

    inline  int isym() const {return isym_;}
    inline  int batch() const {return batch_;}
    inline  int Npart() const {return Npart_;}
    inline  int Ipart() const {return Ipart_;}
    inline  ObservationFlag ObsFlag() const {return ObsFlag_;}
    inline  int run() const {return run_;}
    inline  Hkl hkl()  const {return hkl_;};

    void set_isym(const int& isym) {isym_ = isym;}  // set isym
    void set_hkl(const Hkl& hkl) {hkl_ = hkl;} // set hkl
    void set_batch(const int& Batch) {batch_ = Batch;} // set batch number
    void set_phi(const Rtype& Phi) {phi_ = Phi;} // set Phi
    void offset_phi(const Rtype& Offset) {phi_ += Offset;} // offset Phi
    void offset_time(const Rtype& Offset) {time_ += Offset;} // offset Time
    void negate_time() {time_ = -time_;} // negate Time

  private:
    Hkl hkl_;
    int isym_, batch_;
    Rtype I_, sigI_;           // from column I
    Rtype Ipr_, sigIpr_;           // from column Ipr (if present)
    Rtype Xdet_, Ydet_, phi_, time_;
    Rtype fraction_calc_, width_, LP_;
    int  Npart_, Ipart_;
    ObservationFlag ObsFlag_;
    int run_;
  };  // class observation_part

  //===================================================================
  //! An observation of a reflection, which may consist of one or more parts
  class observation
  {
  public:
    observation();
    observation(const Hkl hkl_in,
                const int& isym_in,
                const int& run_in,
		const int& datasetIndex_in,
		const int& Npart_in,
		observation_part ** const part1_in, 
		const Rtype& TotFrac,
		const PartFlagSwitch& partialstatus_in,
		const ObservationFlag& obsflag_in);

    // Accessors
    Rtype I() const {return I_;}    //!< return I
    Rtype sigI() const {return sigI_;} //!< return sigI
    IsigI I_sigI() const {return IsigI(I_,sigI_);}    //!< return I, sigI

    Rtype kI() const {return I_/gscale;}    //!< return scaled I
    Rtype ksigI() const {return sigI_/gscale;} //!< return scaled sigI
    IsigI kI_sigI() const {return IsigI(I_/gscale,sigI_/gscale);}    //!< return scaled I, sigI

    //! Return "summation" integration IsigI, summed over partials if necessary
    // This is also the sole intensity if there is only one 
    // Also sets mean phi, time, LP
    IsigI IsigIsummation();

    // Return "profile" integration I sigI, summed over partials if necessary
    IsigI IsigIpr() const;

    // Sum (or scale) all partials for this observation
    // Assumes that SelectI has been set up correctly to choose
    // either summation, profile or combined intensity measurements
    // Sets I_, sigI_, phi_, time_, LP_, batch_
    void sum_partials();

    Rtype phi() const {return phi_;}  //!< return rotation angle "phi"
    Rtype time() const {return time_;} //!< return "time"
    int Isym() const {return isym_;}  //!< return symmetry number ISYM
    Rtype LP() const {return LP_;}

    //! return run number
    int run() const {return run_;};

    int datasetIndex() const {return datasetIndex_;}    //!< return dataset index
    bool IsAccepted() const {return obs_status.IsAccepted();} //!< return "accepted" flag

    int num_parts() const;  //!< return number of parts
    observation_part get_part(const int& kpart) const; //!< return kpart'th part
    void replace_part(const int& kpart, const observation_part& obs_part); //!< replace kpart'th part

    Hkl hkl_original() const {return hkl_original_;} //!< return original indices hkl
    bool IsFull() const {return part_flag == FULL;} //!< return true if fully recorded
    Rtype Gscale() const {return gscale;}  //!< return stored inverse (dividing) scale 
    std::pair<float,float> XYdet() const; //!< return average detector coordinates 
    float TotalFraction() const {return totalfraction;} //!< return total fraction
    int Batch() const;  //!< return central batch number

    // Partial status flag
    PartFlagSwitch PartFlag() const {return part_flag;} //!< return partial status

    //! Store I, sigI  incomplete partials are scaled if necessary, return true if scaled
    bool set_IsigI_phi_time(const Rtype& I, const Rtype& sigI,
 			    const Rtype& phi, const Rtype& time,
			    const Rtype& LP);
    //! Set batch number for "central" batch
    void set_batch(const int& batch) {batch_ = batch;}
    //! Store updated sigI
    void set_sigI(Rtype& sigI) {sigI_ = sigI;}
    //! Store inverse scale g
    void SetGscale(const Rtype& g) {gscale = g;}
    //! Store secondary beam directions
    void StoreS2(const float& Thetap, const float& Phip)
    {thetap=Thetap; phip=Phip;}
    //! Return secondary beam directions
    void GetS2(double& Thetap, double& Phip) const 
    {Thetap=thetap; Phip=phip;}
    //! Store diffraction vector
    void StoreS(const FVect3& dStarvec) {s_dif = dStarvec;}
    //! Return diffraction vector at this phi, diffractometer frame
    FVect3 GetS() const {return s_dif;}
    //! Reset status flag from ObservationFlag according to control settings
    /*!  ObservationFlag obs_flag   is the set of bit flags read from the input file
             eg overload: flagged observations may be conditionally accepted by setting
             obs_status.
             ObservationStatus obs_status  is a volatile flag indicating the
             current status of the observation
     */
    void ResetObsAccept(ObservationFlagControl& ObsFlagControl);	
    //! return observation flag
    ObservationFlag Observationflag() const {return obs_flag;}

    //! return current value of volatile status flag
    ObservationStatus ObsStatus() const {return obs_status;}
    //! set status flag
    void UpdateStatus(const ObservationStatus& Status) {obs_status = Status;}

  private:
    Hkl hkl_original_;
    int isym_;
    int run_;
    int datasetIndex_;
    int Npart_;
    observation_part ** part1;
    int batch_;  // batch number of "central" batch
    Rtype totalfraction;
    PartFlagSwitch part_flag;
    ObservationFlag obs_flag;
    Rtype gscale;         // inverse scale g
    Rtype I_, sigI_;
    Rtype phi_;
    Rtype time_;
    Rtype LP_;
    ObservationStatus obs_status;
    Rtype thetap, phip;   // secondary beam direction polar angles
    FVect3 s_dif;             // diffraction vector at Phi setting, 1/A units
  }; // class observation

  //===================================================================
  //! A unique reflection, containing a list of observations
  class reflection
  {
  public:
    reflection();  //!< dummy, don't use

    //! constructor with index into observation_part list from hkl_unmerge_list object
    /*! \param hkl         hkl indices
        \param index_obs1  index of first part in observation_part list 
        \param s           1/d^2  (invresolsq) 
     */
    reflection(const Hkl& hkl, const int& index_obs1, const Dtype& s);
    reflection(const reflection& refl); //!< copy constructor, resets NextObs

    reflection& operator= (const reflection& refl); //!< copy, resets NextObs

    // Storage
    //! store index of last part in observation_part list
    void store_last_index(const int& index_obs2);
    inline int first_index() const {return index_first_obs_;} //!< return 1st index
    inline int last_index() const {return index_last_obs_;} //!< return last index
    inline Hkl hkl() const {return hkl_reduced_;} //!< return reduced hkl

    //! add in an observation
    void add_observation_list(const std::vector<observation>& obs_in);

    //! add all partials, no scales, returns smallest sigma found, return = -1 if no valid observations
    /*! On exit:
      Nfull, number of fulls; Npart, number of partials; Nscaled, number scaled */
    Rtype sum_partials(int& Nfull, int& Npart, int& Nscaled);

    // Retrieval
    int num_observations() const; //!< return number of observations
    int NvalidObservations() const; //!< return number of valid observations
    observation get_observation(const int& lobs) const; //!< return lobs'th observation

    //! Return next valid observation in obs, returns index number, = -1 if end
    int next_observation(observation& obs) const; 
    //! Reset NextObs count 
    void reset() const {NextObs = -1;} // (NextObs is mutable)
    //! replace current observation with updated version
    void replace_observation(const observation& obs);
    //! replace lobs'th observation with updated version
    void replace_observation(const observation& obs, const int& lobs);
    // Count valid observations
    void CountNValid();

    Rtype invresolsq() const {return invresolsq_;} //!< return 1/d^2

    //! Reset observation accepted flags to allow for acceptance of observations flagged as possible errors
    // Counts observations reclassified in ObsFlagControl
    void ResetObsAccept(ObservationFlagControl& ObsFlagControl);

    void SetStatus(const int& istat) {statusflag = istat;} //!< set status (internal use)
    int Status() const {return statusflag;} //!< return status (internal use)

  private:
    Hkl hkl_reduced_; // reduced indices
    std::vector<observation> observations;
    int index_first_obs_, index_last_obs_;
    Rtype invresolsq_;
    mutable int NextObs;
    int NvalidObs;
    int statusflag;  // 0 OK accept, !=0  reject
  }; //class reflection

  //===================================================================
  class hkl_unmerge_list
  //! This is the top-level unmerged reflection object
  /*!
   It contains :- 
\li     1) obs_part_list  a list of raw observations
       (almost a mirror of the file) created by a series of
       calls to store_part followed by close_part_list
\li     2) refl_list  a list of unique reflections, created by
       a call to "organise"
\li     3) within each reflection, a list of observations which
       may have multiple parts (spots), assembled by a call to "partials"
   ie the hierarchy is 
\n   part < observation < reflection
  
    To fill this object:
\li     1) Construct as dummy
\li     2) call "init" to set up
         init may be called with or without dataset & batch
         information: if without, then this must be added later
         with a call to either (a) StoreDatasetBatch, or
         (b) a series of calls to AddDatasetBatch
\li     3) call "store_part" to add each observation part
\li     4) call "close_part_list" to finish list
\li     5) call "StoreDataFlags" to marked which data items are actually present
\li     6) call "prepare" to classify into observations &
        reflections etc

   Note that objects of this class may not be copied, unless they are empty
   (ie created with the default constructor), and also a copy constructor will fail.
   However, you may "append" to an empty object, as a way of copying.
   These objects may be very large, so copying should be avoided where possible.
  
   No reading of external files (eg MTZ file) is done within this class
  */
  {
  public:
    //! construct empty object, must be followed by init
    hkl_unmerge_list();
    //! Initialise for internal writing, with dataset and batch information
    void init(const std::string& Title,
	      const int& NreflReserve, 
	      const hkl_symmetry& symmetry,
	      const all_controls& controls,
	      const std::vector<Xdataset>& DataSets,
	      const std::vector<Batch>& Batches);
    //! Initialise for internal writing,dataset and batch information added later
    void init(const std::string& Title,
	      const int& NreflReserve, 
	      const hkl_symmetry& symmetry,
	      const all_controls& controls);
    //! Clear out list ready for new init
    void clear();

    //! true if list is empty
    bool IsEmpty() const {return (status == EMPTY);}

    //!true if list processed ready for use
    bool IsReady() const {return (status == SUMMED);}

    //! Store datasets & batch info following previous call to init
    void StoreDatasetBatch(const std::vector<Xdataset>& DataSets,
			   const std::vector<Batch>& Batches);

    //! append datasets & batch info following previous call to init
    /*! This may be one of several calls*/
    void AddDatasetBatch(const std::vector<Xdataset>& DataSets,
					   const std::vector<Batch>& Batches);

    //! Add in another hkl_unmerge_list to this one
    /*! Returns status  =  0  OK
                        = +1  different symmetry (point group) */
    int append(const hkl_unmerge_list& OtherList);

    //----
    //! Copy constructor throws exception unless object is EMPTY
    hkl_unmerge_list(const hkl_unmerge_list& List);
    //! Copy operator throws exception unless object is EMPTY
    hkl_unmerge_list& operator= (const hkl_unmerge_list& List);

    // Set controls, limits, etc ----------
    void SetResoLimits(const float& LowReso, const float& HighReso); //!< set resolution

    //! reset to file limits
    void ResetResoLimits();
    //! store ice rings, only store ones flagged as "reject"
    void SetIceRings(const Rings& rings);

    //! If there are any resolution limits set by run, go through the observation
    // list and flag observations which are outside these limits
    void ImposeResoByRunLimits();

    // Retrieve information  -------------------------------
    //! total number of reflections
    int num_reflections() const;
    //! number of valid reflections
    int num_reflections_valid() const {return Nref_valid;}
    //! number of reflections in ice rings
    int num_reflections_icering() const {return Nref_icering;}

    //! total number of observations
    int num_observations() const;
    //! number of fully-recorded observations
    int num_observations_full() const {return Nobs_full;}
    //! number of partially-recorded observations
    int num_observations_part() const {return Nobs_partial;}
    //! number of scaled partial observations
    int num_observations_scaled() const {return Nobs_scaled;}
    //! number of partial observations rejected for too small fraction
    int num_observations_rejected_FracTooSmall() const
    {return partial_flags.NrejFractionTooSmall();}
    //! number of partial observations rejected for too large fraction
    int num_observations_rejected_FracTooLarge() const
    {return partial_flags.NrejFractionTooLarge();}
    //! number of partial observations rejected for gap
    int num_observations_rejected_Gap() const {return partial_flags.NrejGap();}
    //! Minimum sigma(I) of summed partials
    Rtype MinSigma() const {return sigmamin;}
    //! return unit cell for named dataset: if name is blank return average cell for all datasets
    Scell cell(const PxdName& PXDsetName = PxdName()) const;
    Scell Cell() const {return cell(PxdName());}  //!< returns average cell

    //! returns symmetry
    hkl_symmetry symmetry() const {return refl_symm;}

    //! Store MTZ symmetry (just for checking if next file has same symmetry)
    void SetMtzSym(const CMtz::SYMGRP& Mtzsym) {mtzsym = Mtzsym;}
    //! Return MTZ symmetry (for file output)
    CMtz::SYMGRP MtzSym() const {return mtzsym;}

    //! Total reindexing so far (cumulative)
    ReindexOp TotalReindex() const {return totalreindex;}
    // Return resolution
    ResoRange ResRange() const {return ResolutionRange;}      //!< Resolution
    Range Srange() const {return ResolutionRange.SRange();}   //!< smin, smax (1/d^2)
    Range RRange() const {return ResolutionRange.RRange();}   //!< low, high, A
    ResoRange ResLimRange() const {return ResoLimRange;}      //!< Resolution limits
    Rtype DstarMax() const;                                   //!< maximum d* = lambda/d
    // Return dataset stuff
    int num_datasets() const {return ndatasets;} //!< number of datasets
    Xdataset xdataset(const int& jset) const {return datasets.at(jset);}  //!< jset'th dataset
    std::vector<Xdataset> AllXdatasets() const {return datasets;} //!< all datasets
    // Return batch stuff
    int num_batches() const {return nbatches;} //!< number of batches
    int num_accepted_batches() const;          //!< number of accepted batches
    std::vector<Batch> Batches() const {return batches;}  //!< all batches

    //! Actual batch number for batch serial
    Batch batch(const int& jbat) const {return batches.at(jbat);}
    //! Batch serial for given batch
    int batch_serial(const int& batch) const {return batch_lookup.lookup(batch);}

    // Return batch serial number for batch batchnum, or if this one is
    // not present, search upwards until one is found or maxbatchnum is reached.
    // Return -1 if nothing found
    int NextBatchSerial(const int& batchnum, const int& maxbatchnum) const;
    // Return batch serial number for batch batchnum, or if this one is
    // not present, search backwards until one is found or 0 is reached.
    // Return -1 if nothing found
    int LastBatchSerial(const int& batchnum) const;

    // Run stuff
    int num_runs() const {return runlist.size();} //!< number of runs
    std::vector<Run> RunList() const {return runlist;}  //!< all runs

    //! Apply offset to batch numbers, one offset for each run
    void OffsetBatchNumbers(const std::vector<int>& runOffsets);
    //! Mark batch number ibatch as not accepted
    /*! Data records are not changed */
    void RejectBatch(const int& batch);
    //! Mark batch with serial number jbat as not accepted
    /*! Data records are not changed */
    void RejectBatchSerial(const int& jbat);

    //! Remove all observation parts belonging to rejected batches: returns number of parts rejected
    int PurgeRejectedBatches();
    //! Remove all observation parts belonging to specified rejected batches and remove them entirely from the list
    /*! Returns number of parts rejected */
    int EliminateBatches(const std::vector<int> RejectedBatches);

    //! Title from file
    std::string Title() const {return FileTitle;}

    //! Append filename to list of names
    void AppendFileName(const std::string& Name);
    std::string Filename() const {return filename;} //!< returns filename list

    //! return one reflection for index jref, unconditional, no checks
    reflection get_reflection(const int& lref) const;
    //! Return reflection for given (reduced) hkl in refl. Returns index number or -1 if missing
    /*! Records current reflection */
    int get_reflection(reflection& refl, const Hkl& hkl) const;
    //! Return next accepted reflection in argument refl, returns index number  = -1 if end
    int next_reflection(reflection& refl) const;
    //! replace current reflection with updated version
    void replace_reflection(const reflection& refl);

    //!  Get reflection jref, if accepted, else next acceptable one
    //!  Return index of returned reflection, = -1 if end of list
    //! Thread safe, no change in mutable data
    int get_accepted_reflection(const int& jref, reflection& this_refl) const;
      
    //! replace lobs'th observation in current reflection with updated version
    void replace_observation(const observation& obs, const int& lobs);

    //! Reset reflection counter to beginning for "next_reflection"
    /*! not always needed since next_reflection resets it at end */
    void rewind() const;

    //! Retrieve i'th obs_part using pointer list
    observation_part& find_part(const int& i) const;
    //! Total number of parts. Required for MTZ dump, otherwise for internal use
    int num_parts() const {return int(N_part_list);} // length of part list

    // Do things             -------------------------------

    //! Change symmetry or reindex
    /*! If AllowFractIndex true, allow discarding of fractional index
      observations after reindexing, otherwise this is a fatal error
      Returns number of fractional index reflections discarded */
    int change_symmetry(const hkl_symmetry& new_symm,
			 const ReindexOp& reindex_op,
			 const bool& AllowFractIndex = false);

    //! Prepare list for reflection processing
    /*!  if required, sort, organise, partials
      returns number of unique reflections */
    int prepare();
    //! Sum partials
    /*!  return number of partials */
    /*! If forcesum is true, sum them even if already summed */
    int sum_partials(const bool& forcesum=false);

    //! Calculate all secondary beam directions, in chosen frame
    /*! On entry:
      \param pole  =  0 SECONDARY  camera frame
                  !=  0 ABSORPTION, crystal frame = 1,2,3 for h,k,l, 
 		   = -1 unspecified, use closest reciprocal axis for each run */
    void CalcSecondaryBeams(const int& pole);

    //! Set up poles for Absorption
    void SetPoles(const int& pole);

    //! Calculate secondary beam direction for one observation
    /*! On entry
      \param  batchNum     [in]  batch number
      \param  hkl_original [in]  original hkl
      \param  phi          [in]  incident beam rotation, degrees
      \param  sPhi         [out] diffraction vector at actual phi position

     Returns: pair(thetap, phip)   secondary beam directions, radians 
    */
    std::pair<float, float> CalcSecondaryBeamPolar 
    (const int& batchNum, const Hkl& hkl_original, const float& phi,
     DVect3& sPhi) const;

    //! Calculate secondary beam directions
    /*! On entry
      \param  batchNum     [in]  batch number
      \param  hkl_original [in]  original hkl
      \param  phi          [in]  incident beam rotation, degrees
      \param  sPhi         [out] diffraction vector at actual phi position

     Returns  DVect3    secondary beam directions, direction cosines
    */
    DVect3 CalcSecondaryBeam
    (const int& batchNum, const Hkl& hkl_original, const float& phi,
     DVect3& sPhi) const;
      
    //! Reset all observation accepted flags
    /*! Reset all observation accepted flags to allow for acceptance of observations
      flagged as possible errors
      Counts observations reclassified them in ObsFlagControl */
    void ResetObsAccept(ObservationFlagControl& ObsFlagControl);
    //! Reset reflection accepted flags to accept all (subject to resolution checks etc)
    void ResetReflAccept();

  //! Update polarisation corrections for all parts, returns range of corrections
  Range UpdatePolarisationCorrections(const bool& Total,
				     const double& polarisationfactor);
    // If Total == true, then apply complete correction
    //   else assume the unpolarised correction is already applied, apply
    //   additional correction for polarised incident beam
    // polarisationfactor if fraction polarised, = 0 for unpolarised, ~ 0.9 for synchrotrons

    //! Adding observation parts (spots)
    void store_part(const Hkl& hkl_in,
		    const int& isym_in, const int& batch_in,
		    const Rtype& I_in, const Rtype& sigI_in,
		    const Rtype& Ipr, const Rtype& sigIpr,
		    const Rtype& Xdet_in, const Rtype& Ydet_in,
		    const Rtype& phi_in, const Rtype& time_in,
		    const Rtype& fraction_calc_in, const Rtype& width_in,
		    const Rtype& LP_in, 
		    const int& Npart_in, const int& Ipart_in,
		    const ObservationFlag& ObsFlag_in);
    //! finish adding parts
    int close_part_list(const ResoRange& RRange,
			const bool& Sorted);  // returns number of observations

    //! Put data flags to indicate which data items are actually present
    void StoreDataFlags(const data_flags& DataFlags) {dataflags = DataFlags;}
    //! Get data flags to indicate which data items are actually present
    data_flags DataFlags() const {return dataflags;}

    // dump given reflection, for debugging
    void dump_reflection(const Hkl& hkl) const;

    // Set runs
    void AutoSetRun();

  private:
    void initialise(const int NreflReserve, 
		    const hkl_symmetry& symmetry);

    // organise observations into reflections
    int organise();   // returns number of reflections
    // assemble partials into observations, returns number of observations
    int partials();
    void MakeHklLookup() const;

    void SetBatchList();
    void AverageBatchData();

    // Status:-
    //   EMPTY          initial state
    //   RAWLIST        raw observation list read in
    //   SORTED         sorted
    //   ORGANISED      organised into reflections
    //   PREPARED       reflections split into observations
    //                   (ie partials checked and assembled,
    //                   but not summed)
    //   SUMMED         partials summed, reflection.observations valid
    enum rfl_status
      {EMPTY, RAWLIST, SORTED, ORGANISED, PREPARED, SUMMED};
    rfl_status status;
    bool ChangeIndex;
    ReindexOp totalreindex;   // cumulative from original HKLIN setting

    std::vector <observation_part> obs_part_list;
    std::vector <observation_part *> obs_part_pointer; 
    std::vector <reflection> refl_list;
    Rtype sigmamin;
    size_t N_part_list;
    int Nref;
    int Nref_valid;

    int Nobservations;
    int Nobs_full;
    int Nobs_partial;
    int Nobs_scaled;
    // hkl index list: only generated if needed
    mutable clipper::HKL_lookup hkl_lookup;
    mutable bool is_hkl_lookup;   // initially false, true if hkl_lookup has been generated

    mutable int Nref_icering;

    mutable int NextRefNum;  // for next_reflection

    hkl_symmetry refl_symm;
    CMtz::SYMGRP mtzsym;     // Mtz-style symmetry, may be null
			     //  null is mtzsym.spcgrp = -1

    // Average cell (over all datasets)
    Scell averagecell;
    ResoRange ResolutionRange; // low, high, in stored array
    ResoRange ResoLimRange;    // resolution range limits
			       // reflections outside this range will
			       // not be "accept"ed

    // Controls
    int run_set;   // has run be set in part list?
		   // = -1 no controls, = 0 not set, = +1 set
    run_controls run_flags;

    bool partial_set;  // Has partial_flags been set?
    partial_controls partial_flags;

    // Ice rings
    Rings Icerings;

    // Datasets & batches
    std::vector<Xdataset> datasets;
    int ndatasets;
    std::vector<Batch> batches;      // list of batches
    int nbatches;
    hash_table batch_lookup;

    // IsPhiOffset initially false, reset in AutoSetRun
    // true if batch list set up & offset needed
    // Used in store_part to apply offset to data as it is stored,
    // or if StoreDatasetBatch (which calls AutoSetRun) is called
    // after storing the data, then offsets are applied in AutoSetRun
    // Reset to false in close_part_list
    bool IsPhiOffset;

    // Runs
    std::vector<Run> runlist;

    std::string FileTitle;
    std::string filename;

    // Flags for which data items are actually present
    data_flags dataflags;

    // ****  Private member functions
    // sort part list
    void sort();
    // for the organise method
    void add_refl(const Hkl& hkl_red, const int& index, const Dtype& s);
    void end_refl(const int& index);


    void set_run();   // set all runs in part list

    bool accept() const;    // true if reflection(NextRefNum) is acceptable
                            // also counts ice ring
    // private method, no counting or other internal storage
    bool acceptableReflection(const int& RefNum) const;

    //! Set runs
    void SetUpRuns();
    void CheckAllRuns();
    void SetRunInput();

    // Check that all batches in each run have the same or similar orientations
    // Return false if different
    // Store orientation in all runs
    bool StoreOrientationRun();

    // Are new datasets the same as any old ones?
    // For each dataset from other list, store equivalent dataset
    // index in present list, if it is the same dataset
    void MergeDatasetLists(const std::vector<Xdataset>& otherDatasets,
			   const std::vector<Batch>& otherBatches);

    struct ComparePartOrder
    // This construct seems to be a way of getting pointer-to-function
    // into the argument for sort. Copied from the Web.
    {
	bool operator()(const observation_part * part1,
			const observation_part * part2);
    };
  };  // class hkl_unmerge_list
} // end namespace scala

#endif
// end hkl_unmerge
