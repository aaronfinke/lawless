// hkl_unmerge.cpp
//
// Phil Evans August 2003
//
// Classes for unmerged hkl lists
//  at present just for the items used by Mosflm & Scala

// These data structures are primarily directed at Scala so are
// in its namespace

#include <algorithm>

#define ASSERT assert
#include <assert.h>

#include "hkl_unmerge.hh"
#include "scala_util.hh"
#include "string_util.hh"

// Clipper
#include <clipper/clipper.h>
#include "clipper/clipper-ccp4.h"
using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;


namespace scala {
  // -------------------------------------------------------------
  // -------------------------------------------------------------
  data_flags::data_flags()
  {
    is_h = true;
    is_k = true;
    is_l = true;
    is_misym = true;
    is_batch = true;
    is_I = true;
    is_sigI = true;
    // Optional
    is_Ipr = false;
    is_sigIpr = false;
    is_fractioncalc = false;
    is_Xdet = false;
    is_Ydet = false;
    is_Rot = false;
    is_Width = false;
    is_LP = false;
    is_Mpart = false;
    is_ObsFlag = false;
    is_BgPkRatio = false;
    is_scale = false;
    is_sigscale = false;
    is_time = false;
  }
  // ------------------------------------------------------------------
  void data_flags::print() const // for debugging
  {
    std:: cout << "\nDataFlags:\n"
	       << "\nis_h  " << is_h
	       << "\nis_k  " << is_k
	       << "\nis_l  " << is_l
	       << "\nis_misym  " << is_misym
	       << "\nis_batch  " << is_batch
	       << "\nis_I  " << is_I
	       << "\nis_sigI  " << is_sigI
	       << "\nis_Ipr  " << is_Ipr
	       << "\nis_sigIpr  " << is_sigIpr
	       << "\nis_fractioncalc  " << is_fractioncalc
	       << "\nis_Xdet  " << is_Xdet
	       << "\nis_Ydet  " << is_Ydet
	       << "\nis_Rot  " << is_Rot
	       << "\nis_Width  " << is_Width
	       << "\nis_LP  " << is_LP
	       << "\nis_Mpart  " << is_Mpart
	       << "\nis_ObsFlag  " << is_ObsFlag
	       << "\nis_BgPkRatio  " << is_BgPkRatio
	       << "\nis_scale  " << is_scale
	       << "\nis_sigscale  " << is_sigscale
	       << "\nis_time  " << is_time  << "\n\n";
  }
  // ------------------------------------------------------------------
  // ******************************************************************
  //--------------------------------------------------------------
  // Initialise static members
  int SelectI::selecticolflag = 0;  // No Ipr
  int SelectI::ipowercomb = 3;      // Ipower =3
  double SelectI::imid = -1.0;         // unset
  bool SelectI::iprpresent = false; // true if we have a second intensity Ipr stored
  //--------------------------------------------------------------
  void SelectI::SetIcolFlag(const int& IcolFlag, const double& Imid, const int Ipower)
  {
    selecticolflag = IcolFlag;
    ipowercomb = Ipower;
    if (IcolFlag > 0) {
      selecticolflag = 1;
      imid = Imid;
    }
  }
  //--------------------------------------------------------------
  void SelectI::SetAverageIntensity(const double& meanI)
  // store imid = overall <I> if needed and not set
  {
    if (selecticolflag > 0 && imid < 0.1) {
      imid = meanI;
    }
  }
  //--------------------------------------------------------------
  void SelectI::ResetAverageIntensity(const double& meanI)
  // store imid = overall <I> if needed
  {
    if (selecticolflag > 0) {
      imid = meanI;
    }
  }
  //--------------------------------------------------------------
  IsigI SelectI::GetCombinedI(const Rtype& Iraw, const Rtype& Ic, const Rtype& varIc,
			      const Rtype& Ipr, const Rtype& varIpr)
  // Return combined I, sigI 
  {
    // COMBINE option, weighted mean of I & Ipr
    double w = 1.0/(1.0 + pow((std::abs(Iraw)/imid), ipowercomb));
    return IsigI((w*Ipr+(1.0-w)*Ic),
		 sqrt(w*varIpr + (1.0-w)*varIc));
  }
  //--------------------------------------------------------------
  IsigI SelectI::GetCombinedI(const Rtype& Iraw, const IsigI& Isc, const IsigI& Ispr)
  {
    Rtype varIc  = Isc.sigI()*Isc.sigI();
    Rtype varIpr = Ispr.sigI()*Ispr.sigI();
    return GetCombinedI(Iraw, Isc.I(), varIc, Ispr.I(), varIpr);
  }
  //--------------------------------------------------------------
  std::string SelectI::format() 
  {
    std::string s;
    if (selecticolflag < 0) {
      s = "Profile-fitted intensities will be used";
    } else if (selecticolflag == 0) {
      s = "Summation-integration (or sole) intensities will be used";
    } else {
      s = "Combined intensities will be used:\n";
      s += "  weighted mean of profile-fitted (Ipr) & summation (Isum) intensities\n";
      s += "    I = w * Ipr + (1-w) * Isum\n";

      s += "    w = 1/(1+(Iraw/"+
	StringUtil::Strip(StringUtil::ftos(imid,8,1))+
	")^"+StringUtil::Strip(clipper::String(ipowercomb))+")";
    }
    return s;
  }
  // ******************  observation_part  *******************

  observation_part::observation_part() {}  // dummy default constructor
  //              **************
  observation_part::observation_part(const Hkl& hkl_in,
                                     const int& isym_in, const int& batch_in,
                                     const Rtype& I_in, const Rtype& sigI_in,
				     const Rtype& Ipr_in, const Rtype& sigIpr_in,
                                     const Rtype& Xdet_in, const Rtype& Ydet_in,
                                     const Rtype& phi_in, const Rtype& time_in,
                                     const Rtype& fraction_calc_in, const Rtype& width_in,
                                     const Rtype& LP_in, 
                                     const int& Npart_in, const int& Ipart_in,
				     const ObservationFlag& ObsFlag_in)
    :     hkl_(hkl_in),
          isym_(isym_in), batch_(batch_in),
          I_(I_in), sigI_(sigI_in),
          Ipr_(Ipr_in), sigIpr_(sigIpr_in),
          Xdet_(Xdet_in), Ydet_(Ydet_in), phi_(phi_in), time_(time_in),
          fraction_calc_(fraction_calc_in), width_(width_in),
          LP_(LP_in),
          Npart_(Npart_in), Ipart_(Ipart_in), ObsFlag_(ObsFlag_in),
	  run_(1)
  {}
  // ******************  observation  *******************

  observation::observation()
    //        *********
    : part_flag(EMPTY)   // initialise as unfilled
  {}
  //--------------------------------------------------------------
  observation::observation(const Hkl hkl_in,
     //       **********
                           const int& isym_in, const int& run_in,
			   const int& datasetIndex_in,
                           const int& Npart_in,
			   observation_part ** const part1_in,
			   const Rtype& TotFrac,
                           const PartFlagSwitch& partialstatus_in,
			   const ObservationFlag& obsflag_in)
    : hkl_original_(hkl_in),
      isym_(isym_in), run_(run_in), datasetIndex_(datasetIndex_in),
      Npart_(Npart_in), part1(part1_in), batch_(0),
      totalfraction(TotFrac), part_flag(partialstatus_in), obs_flag(obsflag_in),
      gscale(1.0)
  {
    // By default here reject if any flag set
    // This observation may be accepted later if the flags pass a conditional test
    // (see ResetObsAccept)
    if (! obs_flag.OK()) obs_status.SetObsFlag();    
  }

  //--------------------------------------------------------------
  int observation::num_parts() const {return Npart_;}
  //               ^^^^^^^^
  //--------------------------------------------------------------
  observation_part observation::get_part(const int& kpart) const
    //                         ^^^^^^^
  {
    return **(part1+kpart);
  }
  // -------------------------------------------------------------------
  void observation::replace_part(const int& kpart, const observation_part& obs_part)
  {
    **(part1+kpart) = obs_part;
  }
  //--------------------------------------------------------------
  // Central batch number
  int observation::Batch() const
  {return batch_;}
  //--------------------------------------------------------------
  void observation::ResetObsAccept(ObservationFlagControl& ObsFlagControl)
  // Reset observation accepted flags to allow for acceptance of
  // observations flagged as possible errors
  // Counts observations reclassified them in ObsFlagControl
  {
    // Update acceptance flag
    if (ObsFlagControl.IsAccepted(obs_flag)) // true if allowed
      {obs_status.UnsetObsFlag();}
    else
      {obs_status.SetObsFlag();}
  }
  //--------------------------------------------------------------
  bool observation::set_IsigI_phi_time(const Rtype& I, const Rtype& sigI,
 				       const Rtype& phi, const Rtype& time,
				       const Rtype& LP)
  // Store I, sigI  incomplete partials are scaled if necessary
  // Assumes that totalfraction has already been checked
  // return true if scaled
  {
    I_ = I; sigI_ = sigI; phi_ = phi; time_ = time; LP_ = LP;
    if (part_flag == SCALE) {
      I_ /= totalfraction;
      sigI_ /= totalfraction;
      return true;
    }
    return false;
  }
  //--------------------------------------------------------------
  std::pair<float,float> observation::XYdet() const
  // Return average detector coordinates
  {
    float Xd = 0.0;
    float Yd = 0.0;
    if (Npart_ <= 0) return std::pair<float,float>(0.,0.);
    for (int i=0;i<Npart_;++i) {
      Xd += get_part(i).Xdet();
      Yd += get_part(i).Ydet();
    }
    return std::pair<float,float>(Xd/float(Npart_),Yd/float(Npart_));
  }
  //--------------------------------------------------------------
  IsigI observation::IsigIsummation()
  {
    // Return "summation" integration I sigI, summed over partials if necessary
    // This is also the sole intensity if there is only one 
    // Also sets mean phi, time, LP
    Rtype Itot = 0.0;
    Rtype varItot = 0.0;
    IsigI Is;
    observation_part this_part;
    // Stored values
    phi_ = 0.0;
    time_ = 0.0;
    LP_ = 0.0;
    
    if (Npart_ == 1) {
      // Full
      Is = get_part(0).I_sigI();
      phi_ = get_part(0).phi();
      time_ = get_part(0).time();
      LP_ = get_part(0).LP();
      batch_ = get_part(0).batch();
    } else {  // partial
      Rtype max_bit = -1.0;
      batch_ = 0;
      for (int kpart = 0; kpart < Npart_; kpart++) { // loop parts
	this_part = get_part(kpart);
	Itot += this_part.Ic();
	varItot += this_part.sigIc()*this_part.sigIc();
	phi_  += this_part.phi();
	time_  += this_part.time();
	LP_ += this_part.LP();
	// find biggest bit to mark as central batch
	if (this_part.fraction_calc() > max_bit) {
	  max_bit = this_part.fraction_calc();
	  batch_ = this_part.batch();
	}
      } // end loop parts 
      phi_ = phi_/Npart_;  // average phi over all parts
      time_ = time_/Npart_;  // average time over all parts
      LP_ = LP_/Npart_;
      Rtype sigItot = sqrt(varItot);
      if (part_flag == SCALE) {
	Itot /= totalfraction;
	sigItot /= totalfraction;
      }
      // Central batch
      if (batch_ == 0) {
	// Not set, use one in the middle
	batch_ = get_part(Npart_/2).batch();
      }
      Is = IsigI(Itot, sigItot);
    } // end if partial
    return Is;
  }
  //--------------------------------------------------------------
  IsigI observation::IsigIpr() const
  {
    // Return "profile" integration I sigI, summed over partials if necessary
    Rtype Itot = 0.0;
    Rtype varItot = 0.0;
    IsigI Is;
    observation_part this_part;
      
    if (Npart_ == 1) {
      // Full
      Is = get_part(0).I_sigIpr();
    } else {  // partial
      for (int kpart = 0; kpart < Npart_; kpart++) { // loop parts
	this_part = get_part(kpart);
	Itot += this_part.Ipr();
	varItot += this_part.sigIpr()*this_part.sigIpr();
      } // end loop parts 
      Rtype sigItot = sqrt(varItot);
      if (part_flag == SCALE) {
	Itot /= totalfraction;
	sigItot /= totalfraction;
      }
      Is = IsigI(Itot, sigItot);
    } // end if partial
    return Is;
  }
  //--------------------------------------------------------------
  void observation::sum_partials()
  {
    // Sum (or scale) all partials for this observation
    // Assumes that SelectI has been set up correctly to choose
    // either summation, profile or combined intensity measurements
    //
    // Sets I_, sigI_, phi_, time_, LP_, batch_
    // - - - -

    // Get summation integration or sole intensity, sum parts, set phi, time, LP
    IsigI Ic = IsigIsummation();
    IsigI Ipr(0.0, 0.0);
    if (SelectI::IsIprPresent()) {
      // ... and for Ipr if present
      Ipr = IsigIpr();
    }
    IsigI Isum;
    if (SelectI::Combine()) {
      Rtype Iraw = Ic.I();
      if (LP_ > 0.0) Iraw /= LP_;  // raw intensity back-corrected for LP
      Isum = SelectI::GetCombinedI(Iraw, Ic, Ipr);
    } else if (SelectI::SelectIcolFlag() < 0) { // profile
      Isum = Ipr;
    } else {
      Isum = Ic;
    }
    I_ = Isum.I();
    sigI_ = Isum.sigI();
  }
  // ****************** reflection   *******************
  reflection::reflection()  {}   // dummy
  // Normal constructor
  reflection::reflection(const Hkl& hkl, const int& index_obs1,
                         const Dtype& s) 
    : hkl_reduced_(hkl), index_first_obs_(index_obs1),
      invresolsq_(s), statusflag(0)
  {
    NextObs = -1;
  }
  //--------------------------------------------------------------
  reflection::reflection(const reflection& refl)
  // copy constructor, resets NextObs
  {
    hkl_reduced_ = refl.hkl_reduced_;
    observations = refl.observations;
    index_first_obs_ = refl.index_first_obs_;
    index_last_obs_ = refl.index_last_obs_;
    invresolsq_ = refl.invresolsq_;
    NextObs = -1;
    NvalidObs = refl.NvalidObs;
    statusflag = refl.statusflag;
  }
  //--------------------------------------------------------------
  reflection& reflection::operator= (const reflection& refl)
  // Copy resets NextObs
  {
    hkl_reduced_ = refl.hkl_reduced_;
    observations = refl.observations;
    index_first_obs_ = refl.index_first_obs_;
    index_last_obs_ = refl.index_last_obs_;
    invresolsq_ = refl.invresolsq_;
    NextObs = -1;
    NvalidObs = refl.NvalidObs;
    statusflag = refl.statusflag;
    return *this;
  }
  //--------------------------------------------------------------
  void reflection::store_last_index(const int& index_obs2)
    //             ^^^^^^^^^^^^^^
  {index_last_obs_ = index_obs2;}
  
  //--------------------------------------------------------------
  void reflection::add_observation_list(const std::vector<observation>& obs_list)
    //             ^^^^^^^^^^^^^^^^^^
  {
    observations = obs_list;
  }
  //--------------------------------------------------------------
  Rtype reflection::sum_partials(int& Nfull, int& Npart, int& Nscaled)
  // add all partials, no scales, returns smallest sigma found
  //  return = -1 if no valid observations
  //
  // On exit:
  //  Npart   = number of partials
  //  Nscaled = number scaled
  {
    Rtype sdmin0 = +10000000.;
    Rtype sdmin = sdmin0;
    Nfull = 0;
    Npart = 0;
    Nscaled = 0;
    NvalidObs = 0; // counts valid and accepted

    for (int lobs = 0; lobs < num_observations(); lobs++)  {
      observations[lobs].sum_partials();
      if (observations[lobs].IsFull()) {
	if (observations[lobs].IsAccepted()) {NvalidObs++;}
	Nfull++;
      } else {  // partial
	if (observations[lobs].PartFlag() == SCALE) Nscaled++;
	if (observations[lobs].IsAccepted()) {NvalidObs++;}
	Npart++;
      } // end if partial
      if (observations[lobs].sigI() > 0.0) {
	sdmin = Min(sdmin, observations[lobs].sigI());
      }
    } // end loop observations
    if (sdmin > sdmin0*0.9) sdmin = -1.0;  
    return sdmin;
  }
  //--------------------------------------------------------------
  // Methods to return information
  int  reflection::num_observations() const
    //             ^^^^^^^^^^^^^^
    // Return number of observation in reflection
  {
    return observations.size();
  }
  //--------------------------------------------------------------
  int reflection::NvalidObservations() const {return NvalidObs;}
  //--------------------------------------------------------------
  observation reflection::get_observation(const int& lobs) const
    //                    ^^^^^^^^^^^^^^
    // Return lobs'th observation
  {
    NextObs = lobs;  // reset current observation
    return observations[lobs];
  }
  //--------------------------------------------------------------
  int reflection::next_observation(observation& obs) const
  // Return next valid observation, returns -1 if end
  //  On exit: NextObs is index of current observation 
  {
    while (++NextObs < int(observations.size()))  {
      if (observations[NextObs].IsAccepted()) {
	obs = observations[NextObs];
	return NextObs;
      }
    }
    NextObs = -1;
    return NextObs;
  }
  //--------------------------------------------------------------
  void reflection::replace_observation(const observation& obs)
  // replace current observation with updated version
  {
    observations.at(NextObs) = obs;
  }
  //--------------------------------------------------------------
  void reflection::replace_observation(const observation& obs, const int& lobs)
  // replace lobs'th observation with updated version
  {
    observations.at(lobs) = obs;
  }
//--------------------------------------------------------------
  void reflection::ResetObsAccept (ObservationFlagControl& ObsFlagControl)
  // Reset observation accepted flags to allow for acceptance of
  // observations flagged as possible errors
  // Counts observations reclassified in ObsFlagControl
  {
    // count valid observations
    NvalidObs = 0;
    NextObs = -1;
    while (++NextObs < int(observations.size())) {
      observations[NextObs].ResetObsAccept(ObsFlagControl);
      if (observations[NextObs].IsAccepted()) {NvalidObs++;}
    }
    NextObs = -1;
  }
//--------------------------------------------------------------
  void reflection::CountNValid()
  // Count valid observations
  {
    NvalidObs = 0;
    NextObs = -1;
    while (++NextObs < int(observations.size())) {
      if (observations[NextObs].IsAccepted()) {NvalidObs++;}
    }
    NextObs = -1;
  }
//--------------------------------------------------------------

  // ****************** hkl_unmerge_list   *******************
  //--------------------------------------------------------------
  //! construct empty object, must be followed by init
  hkl_unmerge_list::hkl_unmerge_list(): status(EMPTY)
  {
    clear();
  }
  //--------------------------------------------------------------
  // construct list of all reflections and observations from
  // internal calls
  void hkl_unmerge_list::init(const std::string& Title,
			      const int& NreflReserve, 
			      const hkl_symmetry& symmetry,
			      const all_controls& controls,
			      const std::vector<Xdataset>& DataSets,
			      const std::vector<Batch>& Batches)
  {
    // Initialise reflection list with number of reflections and spacegroup
    int Nobspart = NreflReserve;
    initialise(Nobspart, symmetry);

    filename = "";
    FileTitle = Title;

    // store controls
    run_flags = controls.runs;
    run_set = 0;
    partial_flags = controls.partials;
    partial_set = true;

    StoreDatasetBatch(DataSets, Batches);
  }
  //--------------------------------------------------------------
  // construct list of all reflections and observations from
  // internal calls, datasets & batches to be added later 
  void hkl_unmerge_list::init(const std::string& Title,
			      const int& NreflReserve, 
			      const hkl_symmetry& symmetry,
			      const all_controls& controls)
  {
    // Initialise reflection list with number of reflections and spacegroup
    int Nobspart = NreflReserve;
    initialise(Nobspart, symmetry);

    filename = "";
    FileTitle = Title;

    // store controls
    run_flags = controls.runs;
    run_set = 0;
    partial_flags = controls.partials;
    partial_set = true;
    datasets.clear();
    batches.clear();
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::clear()
  // Clear out list ready for new init
  {
    init("Empty list",0,hkl_symmetry(),all_controls());
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::SetBatchList()  // private
  {
    int idataset;
    for (size_t i=0;i<batches.size();i++) {
      if (in_datasets(batches[i].DatasetID(), datasets, idataset)) {
	datasets[idataset].add_batch(batches[i].num());
      }	
    }
    nbatches = batches.size();
    // Sort batch list
    std::sort(batches.begin(), batches.end());
    // Make lookup table (hash table)
    // Setup up hash lookup table a bit larger than required
    batch_lookup.set_size( int(1.2 * nbatches));
    for (int i = 0; i < nbatches; i++)  {
      batch_lookup.add(batches[i].num(), i);
    }
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::AverageBatchData()
  {
    // Average batch cells for each dataset
    // also mosaicity & wavelength, and store in dataset object
    std::vector<float> averageMosaicity, averageWavelength;
    std::vector<Scell> avbcell = AverageBatchCell(batches, ndatasets,
				    averageMosaicity, averageWavelength);
    for (int j=0;j<ndatasets;j++) {
      if (avbcell[j].null()) {
	// Average batch cell is invalid, use dataset cell instead
	avbcell[j] = datasets[j].cell();
      }
      datasets[j].cell() = avbcell[j];
      datasets[j].Mosaicity() = averageMosaicity[j];
      float wvl = averageWavelength[j];
      if (wvl > 0.001) {	datasets[j].wavelength() = wvl;}
    }
    averagecell = AverageDsetCell(datasets);
  }
  //--------------------------------------------------------------
  // Store datasets & batch info following previous call to init
  void hkl_unmerge_list::StoreDatasetBatch(const std::vector<Xdataset>& DataSets,
					   const std::vector<Batch>& Batches)
  {
    // datasets
    datasets.clear();
    for (size_t i=0;i<DataSets.size();i++)
      {datasets.push_back(DataSets[i]);}
    
    // Average unit cells over all datasets & store average
    ndatasets = datasets.size();

    // Average cell, mosaicity & wavelength over all batches for each dataset
    // & store in dataset
    AverageBatchData();

    // Batches
    batches = Batches;
    SetBatchList();

    // Set up run definitions
    SetUpRuns();
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::AddDatasetBatch(const std::vector<Xdataset>& Datasets,
					 const std::vector<Batch>& Batches)
  // append datasets & batch info following previous call to init
  // This may be one of several, terminated by a call to CloseDatasetBatch
  {
    // append datasets & batches if not the same as existing ones
    MergeDatasetLists(Datasets, Batches);
    ndatasets = datasets.size();
    nbatches = batches.size();
    SetBatchList();

    // Average batch cells for each dataset, so far, so that
    // averagecell is available
    AverageBatchData();
    SetUpRuns();
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::OffsetBatchNumbers(const std::vector<int>& runOffsets)
  // Apply offset to batch numbers, one offset for each run
  {
    if (!run_flags.Set()) {
      Message::message(Message_fatal
		       ("hkl_unmerge_list::OffsetBatchNumbers - no runs set"));
    }
    if (run_set == 0) set_run();  // setup runs in list if not already done
    ASSERT (runOffsets.size() == runlist.size());

    bool allzero = true;
    for (size_t irun=0;irun<runlist.size();irun++) {   // Loop runs
      if (runOffsets[irun] != 0) {allzero = false;}
    }
    if (allzero) {return;}  // if all offsets are zero, don't do anything

    for (size_t irun=0;irun<runlist.size();irun++) {   // Loop runs
      std::vector<int> batchlist = runlist[irun].BatchList();  // all batches
      for (size_t ib=0;ib<batchlist.size();ib++) {  // apply offsets to batch list
	batches[batch_lookup.lookup(batchlist[ib])].OffsetNum(runOffsets[irun]);
      }
      // apply offset to batch list in runlist
      runlist[irun].OffsetBatchNumbers(runOffsets[irun]);
    }
    // Remake batch lookup & batch list for each dataset
    batch_lookup.Clear();
    int idataset;
    for (size_t i=0;i<datasets.size();i++) {
      datasets[i].ClearBatchList();
    }
    for (int i = 0; i < nbatches; i++)   {
      batch_lookup.add(batches[i].num(), i);
      if (in_datasets(batches[i].DatasetID(), datasets, idataset)) {
	datasets[idataset].add_batch(batches[i].num());
      }	
    }
    // Reset all batch number in observation part list
    observation_part part;
    for (size_t i = 0; i < N_part_list; i++) {  // loop all raw observations
      part = find_part(i);
      find_part(i).set_batch(part.batch() + runOffsets[part.run()]);
    }
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::RejectBatch(const int& ibatch)
  // Mark batch number ibatch as not accepted
  // Data records are not changed
  {
    batches[batch_lookup.lookup(ibatch)].SetAccept(false);
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::RejectBatchSerial(const int& jbat)
  // Mark batch with serial number jbat as not accepted
  // Data records are not changed
  {
    batches.at(jbat).SetAccept(false);
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::AppendFileName(const std::string& Name)
  {
    if (filename != "") filename += " + ";
    filename += Name;
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::MergeDatasetLists(const std::vector<Xdataset>& otherDatasets,
					   const std::vector<Batch>& otherBatches)
  // Are new datasets (in otherDatasets) the same as any old ones? Append new ones to list
  // For each dataset from otherDataset list, store equivalent dataset
  // index in present list, if it is the same dataset. Return index list
  // Append batches with dataset references
  {
    //    const double Tolerance = 1.0;
    std::vector<int> DtsIndex(otherDatasets.size());
    int ndts = datasets.size(); // number of current datasets
    int kd = ndts - 1;  // index for new datasets appended to old ones
    int id = -1;        // dataset ID for new datasets, largest current id
    for (int j=0;j<ndts;j++) {id = Max(id, datasets[j].setid());}

    for (size_t i=0;i<otherDatasets.size();i++) { // loop new datasets
      DtsIndex[i] = -1;
      for (int j=0;j<ndts;j++) { // loop current datasets
	if (otherDatasets[i] == datasets[j]) {
	  DtsIndex[i] = j;  // i'th "Other" dataset has same name as j'th
	  // add new cell and wavelength into list
	  datasets[j].AddCellWavelength(otherDatasets[i].cell(),
					otherDatasets[i].wavelength());
	  // Check for similar unit cell & wavelength
	  //	  if (!datasets[j].cell().equalsTol(otherDatasets[i].cell(), Tolerance)) {
	  //	    Message::message(Message_warn
	  //	      ("\nWARNING: Datasets with same name have different unit cells, 1st one used"));
	  //	  }
	  break;
	}
      }
      if (DtsIndex[i] < 0) {
	// New dataset
	Xdataset OtherDataset = otherDatasets[i];
	DtsIndex[i] = ++kd;  // new index for i'th dataset
	OtherDataset.setid() = ++id;  // new setid
	datasets.push_back(OtherDataset);
      }
    }
    // Largest file number so far
    int filenum = -1;
    for (int i=0;i<nbatches;i++) {
      filenum = Max(filenum, batches[i].FileNumber());
    }
    filenum++;  // new file number for new batches

    // Append batch list
    for (size_t i=0;i<otherBatches.size();i++) {
      Batch OtherBatch = otherBatches[i];
      // Check that we don't already have this batch number
      if (batch_lookup.lookup(OtherBatch.num()) >= 0) {
	Message::message
	  (Message_fatal("Non-unique batch number"+
			 clipper::String(otherBatches[i].num())));
      }
      // Update dataset index
      int otherIndex = OtherBatch.index(); // old dataset index
      OtherBatch.index() = DtsIndex[otherIndex];  // new index
      OtherBatch.DatasetID() = datasets[DtsIndex[otherIndex]].setid();
      OtherBatch.FileNumber() = filenum;  // store filenumber
      batches.push_back(OtherBatch);
    }
  }
  //--------------------------------------------------------------
  // Add in another hkl_unmerge_list to this one
  int hkl_unmerge_list::append(const hkl_unmerge_list& OtherList)
  // Returns status  =  0  OK
  //                 = +1  different symmetry
  // OtherList is assumed to be compatible:
  // MakeHKL_listscompatible should be run first!
  // Compatibiity means that the two lists have:
  //   1) same point group
  //   2) equivalent indexing schemes if there is ambiguity
  //   3) similar unit cells (within tolerance)
  //   4) unique batch numbers
  {
    // Fail on self-append
    if (&OtherList == this) {
      Message::message(Message_fatal
		       ("hkl_unmerge_list::Append - cannot append to self"));
    }
    // otherwise construct it

    // special for appending to an EMPTY list, initialise it first
    if (status == EMPTY) {
      all_controls othercontrols;
      othercontrols.runs = OtherList.run_flags;
      othercontrols.partials = OtherList.partial_flags;
      init(OtherList.Title(), OtherList.num_parts(), OtherList.symmetry(), othercontrols);
    }

    // First some sanity checks to see if it is allowed
    // These are not exhaustive
    int istat = 0;
    // Same symmetry (ignoring translations)
    if (!refl_symm.equals_r(OtherList.refl_symm)) {
      Message::message(Message_warn
		       ("\nWARNING: Cannot combine reflection lists with different symmetry"));
      istat = 1;
      return istat;
    }

    // Copy all parts
    for (size_t i=0;i<OtherList.N_part_list;i++) {
      obs_part_list.push_back(OtherList.find_part(i));
      N_part_list++;
    }
    if (obs_part_list.size() != N_part_list)
      Message::message
	(Message_fatal("hkl_unmerge_list::close_part - Wrong length list") );
    obs_part_list.resize(N_part_list);
    obs_part_pointer.resize(N_part_list);
    // Set pointer list
    for (size_t i=0;i<N_part_list;i++) {
      obs_part_pointer[i] = &(obs_part_list[i]);
    }
    status = RAWLIST;
    ResolutionRange = ResolutionRange.MaxRange(OtherList.ResRange());
    ResoLimRange = ResolutionRange;  // FIXME?


    // Merge dataset & batch lists
    MergeDatasetLists(OtherList.datasets, OtherList.batches);
    // For each dataset from other list, store equivalent dataset
    // index in present list, if it is the same dataset
    ndatasets = datasets.size();

    // Average cell, mosaicity & wavelength over all batches for each dataset
    // & store in dataset
    AverageBatchData();

    SetBatchList();

    is_hkl_lookup = false;
    // Set up run definitions
    run_set = 0;
    SetUpRuns();
    return istat;
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::initialise(const int NreflReserve,
				    const hkl_symmetry& symmetry)
  //                         ^^^^^^^^^
  {
    refl_symm = symmetry;
    obs_part_list.reserve(NreflReserve);   // reserve space for all observations
    Nref = 0;
    obs_part_list.clear();
    obs_part_pointer.clear();
    refl_list.clear();
    status = EMPTY;
    ndatasets = 0;
    nbatches = 0;
    NextRefNum = -1;
    sigmamin = 0.0;
    IsPhiOffset = false;
    N_part_list = 0;
    Nref = 0;
    Nref_valid = 0;
    Nobservations = 0;
    Nobs_full = 0;
    Nobs_partial = 0;
    Nobs_scaled = 0;
    runlist.clear();
    ResolutionRange = ResoRange();
    ResoLimRange = ResolutionRange;
    mtzsym.spcgrp = -1;
    ChangeIndex = false;
    totalreindex = ReindexOp();
    is_hkl_lookup = false;
    hkl_lookup = clipper::HKL_lookup();
    averagecell = Scell();
    run_set = 0;
    partial_set = false;
    Icerings = Rings();
    datasets.clear();
    batches.clear();
    batch_lookup.Clear();
    filename = "";
    FileTitle = "";
    dataflags = data_flags();
  } // initialise
  //--------------------------------------------------------------
  // Store a raw observation part in list
  void hkl_unmerge_list::store_part(const Hkl& hkl,
  //                     ^^^^^^^^^
                                    const int& isym, const int& batch,
                                    const Rtype& I, const Rtype& sigI,
                                    const Rtype& Ipr, const Rtype& sigIpr,
                                    const Rtype& Xdet, const Rtype& Ydet,
                                    const Rtype& phi, const Rtype& time,
                                    const Rtype& fraction_calc, const Rtype& width,
                                    const Rtype& LP,
                                    const int& Npart, const int& Ipart,
				    const ObservationFlag& ObsFlag)
  {
    Rtype Phi = phi;
    Rtype Time = time;
    if (IsPhiOffset) {   // set if batch list has been set up & phi offsets are needed
      // Apply Phi offset for batch
      float offset = batches.at(batch_lookup.lookup(batch)).PhiOffset();
      Phi += offset;
      if (batches.at(batch_lookup.lookup(batch)).IsTimePhi()) {
	Time += offset;  // also offset time if it is a copy of phi
      }
    }
    obs_part_list.push_back(observation_part(hkl, isym, batch,
                                             I, sigI, Ipr, sigIpr,
                                             Xdet, Ydet, Phi, Time,
                                             fraction_calc, width, LP, 
                                             Npart, Ipart, ObsFlag));
    // Don't set pointer list obs_part_pointer until end (in close_part)
    // in case vector gets extended
    N_part_list++;
  } // store_part

  //--------------------------------------------------------------
  // Close raw observation part list, return number of parts
  int hkl_unmerge_list::close_part_list(const ResoRange& RRange,
					const bool& Sorted)
  {
    if (obs_part_list.size() != N_part_list) {
      Message::message(Message_fatal("hkl_unmerge_list::close_part - Wrong length list") );}
    obs_part_list.resize(N_part_list);
    obs_part_pointer.resize(N_part_list);
    // Set pointer list
    for (size_t i=0;i<N_part_list;i++)  {
      obs_part_pointer[i] = &(obs_part_list[i]);
    }

    // Phi offset already applied, if needed: don't do it again
    IsPhiOffset = false;
    if (N_part_list > 0)  ResolutionRange = RRange.MaxRange(ResolutionRange);
    ResoLimRange = ResolutionRange;  // FIXME?
    status = RAWLIST;
    // List is already sorted if sorted in input file & no change of asu
    if (Sorted) status = SORTED;
    ChangeIndex = false;

    return  N_part_list;
  } // close_part
  //--------------------------------------------------------------
  observation_part& hkl_unmerge_list::find_part(const int& i) const
  //                                  ^^^^^^^^
  // Retrieve i'th obs_part using pointer list
  {
    return  *obs_part_pointer[i];
  }
  //--------------------------------------------------------------
  int hkl_unmerge_list::PurgeRejectedBatches()
  // Remove all observation parts belonging to rejected batches
  // Returns number of parts rejected
  {
    if (status == EMPTY) return 0;
    bool reject = false;
    // Are there any rejected batches?
    for (int ib=0;ib<nbatches;ib++) {
      if (!batches[ib].Accepted()) {
	reject = true;
	break;
      }
    }
    if (!reject) return 0;   // Nothing to do
    unsigned int nrej = 0;
    unsigned int k = 0;
    for (size_t i=0;i<N_part_list;i++) {
      if (batches[batch_lookup.lookup(obs_part_list[i].batch())].Accepted()) {
	// This part belongs to an accepted batch, copy it
	if (k != i) {
	  obs_part_list[k] = obs_part_list[i];
	}
	k++;
      } else {
	// Skip reject observation part
	nrej++;
      }
    }
    ASSERT (N_part_list == k + nrej);
    if (nrej > 0) {
      N_part_list = k;
      obs_part_list.resize(N_part_list);
      obs_part_pointer.resize(N_part_list);
      // Reset pointer list
      for (size_t i=0;i<N_part_list;i++) {
	obs_part_pointer[i] = &(obs_part_list[i]);
      }
      status = RAWLIST;
    }
    return nrej;
  }
  //--------------------------------------------------------------
  int hkl_unmerge_list::EliminateBatches
  (const std::vector<int> RejectedBatches)
  // Remove all observation parts belonging to specified rejected batches
  // and remove them entirely from the list
  // Returns number of parts rejected
  {
    if (status == EMPTY) return 0;
    if (RejectedBatches.size() == 0) return 0;   // Nothing to do

    // Make temporary hash table
    hash_table rejLookUp(int(RejectedBatches.size()*1.5));
    for (size_t i=0;i<RejectedBatches.size();++i) {
      rejLookUp.add(RejectedBatches[i], i);
    }

    unsigned int nrej = 0;
    unsigned int k = 0;
    for (size_t i=0;i<N_part_list;i++) {
      if (rejLookUp.lookup(obs_part_list[i].batch()) < 0) {
	// This part belongs to an accepted batch, copy it
	if (k != i) {
	  obs_part_list[k] = obs_part_list[i];
	}
	k++;
      } else {
	// Skip reject observation part
	nrej++;
      }
    }
    ASSERT (N_part_list == k + nrej);
    if (nrej > 0) {
      N_part_list = k;
      obs_part_list.resize(N_part_list);
      obs_part_pointer.resize(N_part_list);
      // Reset pointer list
      for (size_t i=0;i<N_part_list;i++) {
	obs_part_pointer[i] = &(obs_part_list[i]);
      }
      status = RAWLIST;
      run_set = 0;
    }
    // Remove batches
    k=0;
    for (int ib=0;ib<nbatches;++ib) {
      if (rejLookUp.lookup(batches[ib].num()) < 0) {
	// batch not rejected, copy it
	batches[k] = batches[ib];
	k++;
      }
    }
    batches.resize(k);
    // Clear batch lists for all datasets
    for (size_t id=0;id<datasets.size();id++) {
      datasets[id].ClearBatchList();
    }
    SetBatchList();  // reset batch lists in datasets & batch lookup tables
    return nrej;
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::SetUpRuns()
  // Setup up runs either using automatic criteria or explicit input definitions
  // in run_flags
  {
    // true if runs were specified on input
    bool explicitruns = run_flags.Explicit();
    // Always do the autoset first, since this also does phi offsets etc
    AutoSetRun();
    if (explicitruns) SetRunInput(); // reset runs from input specifications
    CheckAllRuns(); // finish run specification
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::AutoSetRun()
  // Divide batches up into runs
  {
    Run ThisRun;
    int PreviousBatNum = -1;
    float PreviousPhi = 0.0;
    float phioffset;
    double tolerance = 0.01;
    int offset = 0;
    int filenum = -1;
    int batNgap;
    float delPhi = 0.0; 
    float gap = 0.0;
    int nbatAccepted = 0; // number of accepted batches in run

    runlist.clear();
    // Clear dataset run indices
    for (size_t id=0;id<datasets.size();id++) {
      datasets[id].ClearRunList();
    }
    // Loop batches
    for (size_t ib=0;ib<batches.size();ib++) {
      if (ib == 0) {
	// First batch, start run, store dataset index in run
	ThisRun = Run(batch(ib).index());
	// Store run index in dataset
	datasets[batch(ib).index()].AddRunIndex(runlist.size());
	offset  = batch(ib).BatchNumberOffset();
	filenum = batch(ib).FileNumber();
	phioffset = 0.0;
	batches[ib].OffsetPhi(phioffset);
	if (batches[ib].IsTimePhi()) {
	  // if time == phi offset
	  batches[ib].OffsetTime(phioffset);
	}
       } else {
	// Compare this batch with last one
	bool newgroup = false;
	// Conditions for being in the same group (run):
	// same dataset
	if (batch(ib).index() != ThisRun.DatasetIndex()) {newgroup = true;}
	// contiguous batch numbers
	batNgap = batch(ib).num() - (PreviousBatNum+1);  // gap in batch numbers eg 0
	if (batNgap != 0) {newgroup = true;}
	// No phi gap from previous, but allow mod(360)
	gap = batch(ib).Phi1() - PreviousPhi;
	// phioffset must always be a multiple of 360.0
	phioffset = double(Nint(-gap/360.))*360.0;
	if (std::abs(gap) > tolerance) {
	  if (std::abs(gap + phioffset) < tolerance) {
	    // gap = 0(modulo 360)
	    gap = gap + phioffset;
	  } else {
	    newgroup = true;
	    // Allow for gap of one batch
	    if (batNgap == 1) {
	      float dgap = gap - delPhi;
	      if (std::abs(dgap - float(Nint(dgap/360.))*360.) < tolerance) {
		// dgap = 0(modulo 360)
		newgroup = false;
		phioffset = double(Nint(-dgap/360.))*360.0;
		//^
		//		std::cout << "AutoSetRun: not newgroup, batch, gap, phioffset "
		//			  << batch(ib).num() << " " << gap <<" " << phioffset << "\n";
		//^-
	      }
	    }
	  }
	}
	if (newgroup) {
	  //^
	  //	  std::cout << "AutoSetRun: newgroup, batch, batNgap, gap "
	  //		    << batch(ib).num() << " " << batNgap << " " << gap << "\n";
	  //^-
	  // Store completed run as long as it has some accepted batches
	  if (nbatAccepted > 0) {
	    ThisRun.BatchNumberOffset() = offset;
	    ThisRun.FileNumber() = filenum;
	    ThisRun.SortList();
	    ThisRun.RunNumber() = runlist.size()+1;
	    runlist.push_back(ThisRun);
	  }
	  // Start new group, store dataset index
	  ThisRun = Run(batch(ib).index());
	  // Store run index
	  datasets[batch(ib).index()].AddRunIndex(runlist.size());
	  offset  = batch(ib).BatchNumberOffset();
	  filenum = batch(ib).FileNumber();
	  phioffset = 0.0;
	  batches[ib].OffsetPhi(phioffset);
	  if (batches[ib].IsTimePhi()) {
	    // if time == phi offset
	    batches[ib].OffsetTime(phioffset);
	  }
	  nbatAccepted = 0;
	} else {
	  // still in same group, apply offset to stored phi values
	  batches[ib].OffsetPhi(phioffset);
	  //^
	  //	  std::cout << "ib, offset " << ib <<" "<< phioffset <<"\n";
	  //^-
	  if (batches[ib].IsTimePhi()) {
	    // if time == phi offset
	    batches[ib].OffsetTime(phioffset);
	  }
	  if (std::abs(phioffset) > 0.1) {IsPhiOffset = true;}
	  // Sanity checks
	  ASSERT (offset == batch(ib).BatchNumberOffset());
	  //^		ASSERT (filenum == batch(ib).FileNumber());
	}
      }
      // Add batch to run even if not accepted
      ThisRun.AddBatch(batch(ib).num(), batch(ib).Accepted());
      if (batch(ib).Accepted()) {
	nbatAccepted++; // count accepted batches
	// Store run index in batch: this = current size of runlist
	batches[ib].SetRunIndex(runlist.size());
      } else {
	// Store null run index in batch
	batches[ib].SetRunIndex(-1);
      }
      PreviousBatNum = batch(ib).num();
      PreviousPhi = batch(ib).Phi2();
      delPhi =  batch(ib).Phi2() - batch(ib).Phi1();
      // Store run index in batch: this = current size of runlist
      batches[ib].SetRunIndex(runlist.size());
    }  // end loop batches
    // Store last run as long as it has some accepted batches
    if (nbatAccepted > 0) {
      ThisRun.BatchNumberOffset() = offset;
      ThisRun.FileNumber() = filenum;
      ThisRun.SortList();
      ThisRun.RunNumber() = runlist.size()+1;
      runlist.push_back(ThisRun);
    }
    run_flags.SetStatus(0);  // status set to indicate auto run assignment

    // All runs now defined
  }
  //--------------------------------------------------------------
  int hkl_unmerge_list::num_accepted_batches() const          //!< number of accepted batches
  {
    int n = 0;
    for (size_t i=0;i<batches.size();++i) {
      if (batches[i].Accepted()) {n++;}
    }
    return n;
  }
  //--------------------------------------------------------------
  int hkl_unmerge_list::NextBatchSerial(const int& batchnum, const int& maxbatchnum) const
  // Return batch serial number for batch batchnum, or if this one is
  // not present, search upwards until one is found or maxbatchnum is reached.
  // Return -1 if nothing found
  {
    int ibatch = batchnum;
    while (ibatch <= maxbatchnum) {
      // return serial number if found
      if (batch_serial(ibatch) >= 0) return batch_serial(ibatch);
      ibatch++;
    }
    return -1; // not found
  }
  //--------------------------------------------------------------
  int hkl_unmerge_list::LastBatchSerial(const int& batchnum) const
  // Return batch serial number for batch batchnum, or if this one is
  // not present, search backwards until one is found or 0 is reached.
  // Return -1 if nothing found
  {
    int ibatch = batchnum;
    while (ibatch >= 0) {
      // return serial number if found
      if (batch_serial(ibatch) >= 0) return batch_serial(ibatch);
      ibatch--;
    }
    return -1; // not found
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::SetRunInput()
  // Set up runs from input specification
  {
    // List of run numbers specified
    std::vector<int> runnumberlist = run_flags.RunNumberList();
    int nruns = runnumberlist.size();
    if (nruns <= 0) return; // Null list, leave as autoset

    // Clear run specs
    runlist.clear();
    // Clear dataset run indices
    for (size_t id=0;id<datasets.size();id++) {
      datasets[id].ClearRunList();
    }
    //  mark all batches as omitted
    int maxbatchnum = -1; // maximum batch number
    for (size_t ib=0;ib<batches.size();ib++) {
      // Store null run index in batch
      batches[ib].SetRunIndex(-1);
      maxbatchnum = Max(maxbatchnum, batches[ib].num());
    }

    for (int irun=0;irun<nruns;++irun) {    // Loop runs
      int runnum = runnumberlist[irun];  // run number
      // list of batch ranges for this run
      std::vector<IntRange> batchranges = run_flags.BatchRanges(runnum);
      if (batchranges.size() > 0) {
	int nbatAccepted = 0;
	// First batch number in first range
	int ib0 = NextBatchSerial(batchranges[0].min(), maxbatchnum);
	if (ib0 < 0) {
	  // Batch range not found
	  clipper::String br = clipper::String(batchranges[0].min(),6)+
	    " to "+clipper::String(batchranges[0].max(),6);
	  Message::message(Message_fatal
		("hkl_unmerge_list:: batch range not found "+br));
	}
	//	start run, store dataset index in run	
	Run ThisRun = Run(batch(ib0).index());      
	for (size_t i=0;i<batchranges.size();++i) { // loop batch ranges
	  int ibs1 = NextBatchSerial(batchranges[i].min(), maxbatchnum);
	  if (ibs1 < 0) {
	    // Batch range not found
	    clipper::String br = clipper::String(batchranges[i].min(),6)+
	      " to "+clipper::String(batchranges[i].max(),6);
	    Message::message(Message_fatal
			     ("hkl_unmerge_list:: batch range not found "+br));
	  }
	  // batch serial for end of range
	  int ibs2 = LastBatchSerial(batchranges[i].max());
	  if (ibs2 < 0) {
	    // reset to last batch serial
	    ibs2 = num_batches()-1;
	  }
	  for (int ib=ibs1;ib<=ibs2;++ib) {
	    // Add batch to run even if not accepted
	    ThisRun.AddBatch(batch(ib).num(), batch(ib).Accepted());
	    if (batch(ib).Accepted()) {
	      nbatAccepted++; // count accepted batches
	      // Store run index in batch: this = current size of runlist
	      batches[ib].SetRunIndex(runlist.size());
	    } else {
	      // Store null run index in batch
	      batches[ib].SetRunIndex(-1);
	    }
	  } // end loop batches in range
	} // end loop ranges
	if (nbatAccepted > 0) {
	  // Store run index in dataset
	  datasets[batch(ib0).index()].AddRunIndex(runlist.size());
	  ThisRun.BatchNumberOffset() = batch(ib0).BatchNumberOffset();
	  ThisRun.FileNumber() = batch(ib0).FileNumber();
	  ThisRun.SortList();
	  ThisRun.RunNumber() = runlist.size()+1;
	  runlist.push_back(ThisRun);
	}
      }
    } // end loop runs
    //  mark batches not in a run as omitted
    int nbrej = 0;
    for (size_t ib=0;ib<batches.size();ib++) {
      if (batches[ib].RunIndex() < 0) {
	batches[ib].SetAccept(false);
	nbrej++;
      }
    }
    if (nbrej > 0) {
      PurgeRejectedBatches(); // flag observations for rejected batches
    }
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::CheckAllRuns() {
    // Check each run for status of time information
    // Also add in any resolution range cutoffs from run_flags (run_controls)
    std::vector<bool> negateTimeInBatch(nbatches, false); // true if time negated in batch
    bool negateTime = false; // true if any to be negated
    for (size_t irun=0;irun<runlist.size();++irun) {
      // Just look at the first batch in run for time status
      //  batch serial number for first batch
      int ib0 = batch_lookup.lookup(runlist[irun].Batch0());
      runlist[irun].SetBatchSerial0(ib0);  // store
      runlist[irun].ValidTime() = batches[ib0].IsValidTime();
      if (batches[ib0].IsValidTime()) {
	// batch has time or phi information
	if (batches[ib0].IsTimePhi()) {
	  //  time is actually phi
	  if (batches[ib0].PhiRange() < 0.0) {
	    //  phi is decreasing?	    
	    // We need to loop all batches in run to fix up things 
	    //  batch number list
	    std::vector<int> runbatchlist = runlist[irun].BatchList();
	    for (size_t batch=0;batch<runbatchlist.size();++batch) {
	      int ib = batch_lookup.lookup(runbatchlist[batch]);
	      negateTimeInBatch[ib] = true;  // negate time in this batch
	      batches[ib].StoreTimeRange(-batches[ib].Time1(), -batches[ib].Time2());
	      batches[ib].SetTimeReversedFromPhi();
	      negateTime = true;
	    }	    
	  }
	}}
    }

    // If phioffsets are needed and we have already read in the data, then we need
    // to offset phi & possibly time = phi
    // Time from phi may also be negated so that it is increasing
    if ((IsPhiOffset || negateTime) && N_part_list > 0) {
      // If there are phi offsets, apply offsets to all observation parts
      observation_part part;
      for (size_t i = 0; i < N_part_list; i++) { // loop all raw observations
	int batch = find_part(i).batch();
	int ib = batch_lookup.lookup(batch);
	find_part(i).offset_phi(batches[ib].PhiOffset());
	if (batches[ib].IsTimePhi()) {
	  // If time is just a copy of phi, then offset this too
	  find_part(i).offset_time(batches[ib].PhiOffset());
	  if (negateTimeInBatch[ib]) find_part(i).negate_time();
	}
      }
    }
    StoreOrientationRun();

    // resolution ranges, if any
    if (run_flags.IsResoByRun()) {
      // run number, resorange pairs from input
      std::vector<std::pair<int,ResoRange> > resobyrun = run_flags.GetResoByRun();
      for (size_t i=0;i<resobyrun.size();++i) { // loop runs with specified limits
	int irun = FindRunIndex(resobyrun[i].first, runlist); // run index
	if (irun < 0) { // specified run not found
	  clipper::String s = "\n**** Run "+clipper::String(resobyrun[i].first)+
	    " specified on RESOLUTION RUN command does not exist ****\n";
	  Message::message(Message_fatal(s));
	}
	// Store resolution range limit for this run, forced to be within file range
	runlist[irun].StoreResoRange(resobyrun[i].second.MinRange(ResolutionRange));
      }  // end loop specified limits
    }
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::set_run()
  //                     ^^^^^^^
  // Set run numbers into obs_part list using run index from batch
  // 
  {
    if (status != RAWLIST && status != SORTED) 
      Message::message(Message_fatal("hkl_unmerge_list::set_run - not RAWLIST or SORTED") );

    for (size_t i = 0; i < N_part_list; i++) {  // loop all raw observations
      // Assign run
      int irun = batches[batch_lookup.lookup(find_part(i).batch())].RunIndex();
      find_part(i).set_run(irun);
    }
    run_set = +1;
  } // set_run
  //--------------------------------------------------------------
  bool hkl_unmerge_list::StoreOrientationRun()
  // Check that all batches in each run have the same or similar orientations
  // Return false if different
  // Store orientation in all runs
  {
    bool sameOrientation = true;
    DMat33 DUB0;
    DVect3 s0;
    clipper::Rotation U0inv, U;
    // Limit for monitoring large orientation change, 3 degrees
    double rotlim = clipper::Util::d2rad(3.0);

    // Loop runs
    for (size_t irun=0;irun<runlist.size();++irun) {
      bool allvalid = true;
      bool firstBatch = true;
      int lastBatchIdx = 0; // index to last accepted batch

      // Loop batches in run
      std::vector<int> batchlist = runlist[irun].BatchList();
      for (size_t ib=0;ib<batchlist.size();ib++) {  // apply offsets to batch list etc
	Batch batch = batches[batch_lookup.lookup(batchlist[ib])];
	if (batch.Accepted()) {
	  lastBatchIdx = batch_lookup.lookup(batchlist[ib]);
	  // Note that any phi offsets have already been applied in AutoSetRun
	  if (firstBatch) {
	    // 1st batch
	    runlist[irun].PhiRange().first() = batch.Phi1();
	    if (batch.IsValidTime()) {runlist[irun].TimeRange().first() = batch.Time1();}
	    else {runlist[irun].TimeRange().first() = batch.Phi1();}   // substitute Phi for time if missing
	  }
	  if (batch.ValidOrientation()) {
	    if (firstBatch) {
	      // 1st batch, store things
	      U0inv = clipper::Rotation(batch.Umat().inverse());
	      runlist[irun].StoreSpindleToPrincipleAxis(batch.SpindleToPrincipleAxis());
	      runlist[irun].SetValidOrientation(true);
	    } else {
	      // test for change of orientation
	      double rot = (clipper::Rotation(batch.Umat()) * U0inv).abs_angle();
	      if (rot > rotlim) {
		Message::message(Message_warn
				 ("WARNING: problem in run "+clipper::String(int(irun+1))+
				  " batch "+clipper::String(batchlist[ib])+
				  "  has a different orientation that of initial batch "+
				  clipper::String(batchlist[0])+
				  "\n Orientation difference = "+
				  clipper::String(clipper::Util::rad2d(rot))+"\n"));
		sameOrientation = false;
	      }
	    }  // end valid orientation
	  } else {
	  // No valid orientation
	    allvalid = false;
	  }
	  firstBatch = false;
	} // end batch accepted
      } // end loop batches

      // Update information from last accepted batch
      // last batch
      Batch batch = batches[lastBatchIdx];
      runlist[irun].PhiRange().last() = batch.Phi2();
      if (batch.IsValidTime()) {runlist[irun].TimeRange().last() = batch.Time2();}
      else {runlist[irun].TimeRange().last() = batch.Phi2();} // substitute Phi for time if missing

      if (!allvalid) {
	runlist[irun].SetValidOrientation(false);
      }
      runlist[irun].PhiRange().AllowDescending(); // allow negative phi range 
    } // end loop runs
    return sameOrientation;
  }
  //--------------------------------------------------------------
  int hkl_unmerge_list::organise()
    //                  ^^^^^^^
    // Set up reflection list refl_list from raw observation list
    // return number of reflections
  {
    if (status == RAWLIST)
      sort();
    if (status != SORTED) 
      Message::message(Message_fatal(
                                     "hkl_unmerge_list::organise - not SORTED") );

    if (run_set < 0) 
      Message::message(Message_fatal("hkl_unmerge_list::organise - run not set") );
    if (run_set == 0) set_run();  // setup runs in list if not already done
    
    const int MaxIndex = 99999999;
    Hkl lhkl = Hkl(MaxIndex,MaxIndex,MaxIndex);
    Hkl hkl;

    if (Nref !=0) {
      refl_list.resize(0);      // List has been used before, clear & allocate memory
      refl_list.reserve(Nref);  // for same size as before
    }
    Nref = 0;
    
    for (size_t i = 0; i < N_part_list; i++) {  // loop all raw observations
      hkl = find_part(i).hkl();  // reduced indices from list
      if (hkl != lhkl) {
	// New (or first) reflection
	if (i != 0) end_refl(i-1);
	add_refl(hkl, i, hkl.invresolsq(averagecell));
      }
      lhkl = hkl;
    }
    end_refl(N_part_list-1);
    refl_list.resize(Nref);
    status = ORGANISED;
    return Nref;
  } // organise
  //--------------------------------------------------------------
  // start a new reflection in list: called by organise
  void hkl_unmerge_list::add_refl(const Hkl& hkl_red,
  //                     ^^^^^^^
                                  const int& index, const Dtype& s)
  {
    // Create new reflection, store hkl, & index to first observation
    // add it to list    
    refl_list.push_back(reflection(hkl_red, index, s));
    Nref++;
  } // add_refl

  //--------------------------------------------------------------
  void hkl_unmerge_list::end_refl(const int& index)
    //                   ^^^^^^^
    // Mark end of this reflection, ie change of hkl
    //   called by organise
  {
    refl_list[Nref-1].store_last_index(index);
  }


  //--------------------------------------------------------------
  int hkl_unmerge_list::partials()
    //                   ^^^^^^^
    // Allocate observations within reflections to partials
    // ie create observation list for each reflection
  // Returns number of observations
  {
    if (status != ORGANISED) 
      Message::message(Message_fatal("hkl_unmerge_list::partials - not ORGANISED") );
    if (!partial_set)
      Message::message(Message_fatal("hkl_unmerge_list::partials - no partial selection information") );
    Nobservations = 0;
    int batchgap;
    partial_flags.Clear();
		
    // Temporary store for observation list for each reflection
    //   this can expand beyond allocated length if necessary
    std::vector <observation> obs_list(100);

    // Clear run reflection counts
    for (size_t irun=0;irun<runlist.size();++irun) {
      runlist[irun].clearCounts();
    }
    // bool combine = SelectI::Combine(); // true if we want average I for combination
    Rtype avI; // for each observation
    MeanSD meanI;

    for (size_t j = 0; j < refl_list.size(); j++) {  // loop all reflections
      obs_list.clear();   // clear temporary list
      int i = refl_list[j].first_index();
      bool obsOK = false; // true if at least one accepted observations
      while (i <= refl_list[j].last_index())  {
	// start possible observation
	//  Set values for first part or full
	int i1 = i;
	int Nfound = 1;
	int isym1 = find_part(i).isym();
	int batch1 = find_part(i).batch();
	int run1 = find_part(i).run();
	int datasetIndex = runlist[run1].DatasetIndex();
	Rtype total_fraction = find_part(i).fraction_calc();
	bool check_ok = partial_flags.check(); // false if no check on Mpart
	bool scale_frac = false; // don't scale incomplete partial
	PartFlagSwitch partial_status = FULL;
	ObservationFlag obsflag(find_part(i).ObsFlag());
	avI = find_part(i).Ic();
	  
	int kpart = 1;
	int Npart = find_part(i).Npart();
	
	if (Npart != 1) {
	  // Started a Partial (not classified as full)
	  if (partial_flags.check()) 
	    // Check 1st part for consistency of Mpart flags
	    {if (find_part(i).Ipart() != 1) check_ok = false;}
	  // gap flag, store maximum gap between batches, should = 1
	  
	  batchgap = 0;  
	  while (++i <= refl_list[j].last_index()) {
	    // Loop through parts 2->EndObs
	    // Tests for still same observation
	    // same symmetry
	    if (isym1 != find_part(i).isym()) break;
	    //same run
	    if (run1 != find_part(i).run()) break;
	    // check contiguous batches: count gaps, should == 0
	    int gap = (find_part(i).batch() - (batch1+kpart));
	    if (gap > 2) break;
	    batchgap += gap;
	    // Found another part belonging to this observation
	    // Check for rejection
	    obsflag.AddFlag(find_part(i).ObsFlag());
	    kpart++; // kpart counts accepted part
	    if (partial_flags.check()) 
	      // Check for consistency of Mpart flags
	      {if (find_part(i).Ipart() != kpart) check_ok = false;}
	    
	    total_fraction += find_part(i).fraction_calc();
	    avI += find_part(i).Ic();	    
	  }
	  Nfound = kpart;
	  avI /= Rtype(Nfound);

	  if (partial_flags.check()) {
	    // Check for consistency of Mpart flags
	    if (Npart == Nfound) 
	      {partial_status = COMPLETE_CHECKED;}
	    else
	      {check_ok = false;}
	  }
	  if (! check_ok) {
	    // Check for acceptability & completeness unless
	    // Mpart check is all OK
	    if (total_fraction < partial_flags.accept_fract_min()) {
	      if (partial_flags.correct_fract_min() > 0.0001
		  && total_fraction >= partial_flags.correct_fract_min()) {
		check_ok = true;   // accept &
		scale_frac = true; // scale incomplete partial
		partial_status = SCALE; 
	      } else
		{partial_flags.IncrementNrejFractionTooSmall();}
	    } else if (total_fraction > partial_flags.accept_fract_max())
	      {partial_flags.IncrementNrejFractionTooLarge();}
	    else {
	      check_ok = true; // total fraction in range
	      partial_status = COMPLETE; 
	    }
	  }
	  // Gap check
	  if (batchgap > partial_flags.maxgap()) {
	    check_ok = false;
	    partial_flags.IncrementNrejGap();
	  }
	  // end of observation, all parts
	} // partial
	else {
	  // full
	  check_ok = true;
	  i++;  //i indexes next part after this observation
	}
	if (check_ok) {
	  // Store observation
	  // If any part of obsflag is set, mark observation as REJECTED for now
	  Nobservations += 1;
	  obs_list.push_back(observation(
			 refl_symm.get_from_asu(refl_list[j].hkl(), isym1),
			 isym1, run1, datasetIndex, Nfound,
			 &obs_part_pointer[i1],
			 total_fraction, partial_status, obsflag));
	  obsOK = true;
	  if (partial_status == FULL) {
	    runlist[run1].Nfulls()++;
	  } else {
	    runlist[run1].Npartials()++;
	  }
	  meanI.Add(avI);
	}
      } // observation loop
      refl_list[j].add_observation_list(obs_list);
    } // reflection loop

    // Set flags into runs for ony||few fulls||partials
    for (size_t irun=0;irun<runlist.size();++irun) {
      runlist[irun].SetFullsAndPartials();
    }

    if (SelectI::Combine()) {
      // if we want average I for combination
      if (meanI.Count() > 0) SelectI::SetAverageIntensity(meanI.Mean());
    }
    status = PREPARED;
    ImposeResoByRunLimits();  // mark observations if outside run limits

    return Nobservations;
  } // end partials
  //--------------------------------------------------------------
  int hkl_unmerge_list::sum_partials(const bool& forcesum)
  //                   ^^^^^^^^^^^
  // Sum partials within each observation for all reflections
  // If forcesum is true, sum them even if already summed
  // Returns number of valid partials
  {
    if (!forcesum && (status == SUMMED)) return Nobs_partial;
    if (!(status == PREPARED) && (status != SUMMED)) 
      Message::message(Message_fatal("hkl_unmerge_list::sum_partials - not PREPARED") );

    sigmamin = +1000000.;
    double sm;
    Nref_valid = 0;
    Nobs_full = 0;
    Nobs_partial = 0;
    Nobs_scaled = 0;
    int Nfull, Npart, Nscaled;

    for (size_t j = 0; j < refl_list.size(); j++) {  // loop all reflections
      // Sum partials
      //    reflection.sum_partials returns min sigma found
      //    (excluding zeroes)
      sm = refl_list[j].sum_partials(Nfull, Npart, Nscaled);
      if (sm > 0.0) {
	sigmamin = Min(sigmamin, sm);
	Nref_valid++;
	Nobs_full += Nfull;
	Nobs_partial += Npart;
	Nobs_scaled += Nscaled;
      }
    }
    status = SUMMED;
    NextRefNum = 0;        // point to first reflection in list
    return Nobs_partial;
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::ImposeResoByRunLimits()
  // If there are any resolution limits set by run, go through the observation
  // list and flag observations which are outside these limits
  {
    if (!((status == SUMMED) || (status == PREPARED))) {
      Message::message(Message_fatal
		       ("hkl_unmerge_list::ImposeResoByRunLimits - not PREPARED or SUMMED") );
    }
    if (!run_flags.IsResoByRun()) {return;}

    observation this_obs;
    ObservationStatus obs_status;

    // * * * * Loop reflections
    for (size_t j=0;j<refl_list.size();++j) {
      for (int iobs=0;iobs<refl_list[j].num_observations();++iobs) { // loop observations
	this_obs = refl_list[j].get_observation(iobs);
	int irun = this_obs.run();  // run index
	int ibatch = this_obs.Batch(); // batch number (central slot)
	obs_status =  this_obs.ObsStatus();
	// resolution limit to test
	if (runlist[irun].IsResoRange() && 
	    !runlist[irun].InResoRange(refl_list[j].invresolsq(), ibatch)) {
	  // outside limits
	  obs_status.SetResolution(); // set resolution reject flag
	} else {
	  obs_status.UnSetResolution(); // unset resolution reject flag
	}
	this_obs.UpdateStatus(obs_status);
	refl_list[j].replace_observation(this_obs);
      } // end loop observations
    } // end loop reflections
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::SetResoLimits(const float& LowReso,
                                       const float& HighReso)
  {
    // resolution range limits & bins
    ResoLimRange.SetRange(LowReso, HighReso, Nobservations); 
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::ResetResoLimits()
  {
    // resolution range limits rest to file range
    ResoLimRange = ResolutionRange;
  }
  //--------------------------------------------------------------
  Rtype hkl_unmerge_list::DstarMax() const
  // maximum d* = lambda/d
  {
    Rtype dsmax = 0.0;
    for (int id=0;id<ndatasets;++id) {
      Rtype ds = datasets[id].wavelength()/ResolutionRange.ResHigh();
      dsmax = Max(dsmax, ds);
    }
    return dsmax;
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::SetIceRings(const Rings& rings)
  {
    // Only store rings flagged as "reject"
    Icerings.CopyRejRings(rings);
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::MakeHklLookup() const
  // Make hkl lookup table to lookup reflections
  //  (not really const, but [is_]hkl_lookup mutable
  {
    std::vector<clipper::HKL> hkls(Nref); // list of all hkl (as clipper class)
    for (int i=0;i<Nref;i++)
      {hkls[i] = get_reflection(i).hkl().HKL();}
    hkl_lookup.init(hkls);
    is_hkl_lookup = true;
  }
  //--------------------------------------------------------------
  // Methods to return information
  int hkl_unmerge_list::num_reflections() const
    //                  ^^^^^^^^^^^^^
    // Return number of reflections in list
  {
    return Nref;
  }
  //--------------------------------------------------------------
  int hkl_unmerge_list::num_observations() const
    //                  ^^^^^^^^^^^^^
    // Return number of observations in list
  {
    return Nobservations;
  }
  //--------------------------------------------------------------
  Scell hkl_unmerge_list::cell(const PxdName& PXDsetName) const
    {
      // retrieve cell for a dataset
      // If dataset name blank, get average over all datasets
      if (ndatasets <= 0) 
        {
          Message::message(Message_fatal(
                                     "hkl_unmerge_list::cell - no datasets") );
        }

      if (PXDsetName.is_blank())
        {
            return averagecell;
        }
      // Find dataset
      for (int k=0; k<ndatasets; k++)
        if (PXDsetName == datasets[k].pxdname())
            return datasets[k].cell();
      Message::message(Message_fatal(
                      "hkl_unmerge_list::cell - dataset not found "+PXDsetName.format()) );
      return Scell(); // dummy
    }
  //--------------------------------------------------------------
  void hkl_unmerge_list::rewind() const
  {
    NextRefNum = -1;
  }
  //--------------------------------------------------------------
  bool hkl_unmerge_list::accept() const
    // private method, called from next_reflection
  {
    // Accepted?
    if (refl_list[NextRefNum].Status() != 0) {
      return false;
    }
    // Reject if no observations
    if (refl_list[NextRefNum].NvalidObservations() <= 0)
      return false;
    // Test resolution limits
    Rtype s2 = refl_list[NextRefNum].invresolsq();
    if (ResoLimRange.tbin(s2) < 0)
      return false;
    // Check ice rings: only "reject" rings are stored
    if (Icerings.Nrings() > 0) {
      int ir = Icerings.InRing(s2);
      if (ir >= 0) {
	Nref_icering++; // mutable
	return false; 
      }
    }
    return true;
  }
  //--------------------------------------------------------------
  bool hkl_unmerge_list::acceptableReflection(const int& RefNum) const
  // private method, no counting or other internal storage
  {
    // Accepted?
    if (refl_list[RefNum].Status() != 0) {
      return false;
    }
    // Reject if no observations
    if (refl_list[RefNum].NvalidObservations() <= 0)
      return false;
    // Test resolution limits
    Rtype s2 = refl_list[RefNum].invresolsq();
    if (ResoLimRange.tbin(s2) < 0)
      return false;
    // Check ice rings: only "reject" rings are stored
    if (Icerings.Nrings() > 0) {
      int ir = Icerings.InRing(s2);
      if (ir >= 0) {
	return false; 
      }
    }
    return true;
  }
  //--------------------------------------------------------------
  int hkl_unmerge_list::next_reflection(reflection& refl) const
  {
    if (status != SUMMED) 
      Message::message(Message_fatal(
	        "hkl_unmerge_list::next_reflection - not SUMMED") );

    if (NextRefNum < 0) {
      // First time, some initialisations
      Nref_icering = 0;
    }

    while (++NextRefNum < Nref) {
      if (accept()) {   // test acceptance criteria (eg resolution)
	refl = get_reflection(NextRefNum);
	return NextRefNum;
      }
    }
    NextRefNum = -1;
    return NextRefNum;
  }
  //--------------------------------------------------------------
  int hkl_unmerge_list::get_accepted_reflection(const int& jref, reflection& this_refl) const
  //  Get reflection jref, if accepted, else next acceptable one
  //  Return index of returned reflection, = -1 if end of list
  // Thread safe, no change in mutable data
  {
    int jr = jref;
    while (jr < Nref) { 
      if (acceptableReflection(jr)) {
	this_refl = refl_list[jr];
	return jr;
      }
      jr++;
    }
    return -1;
  }
  //--------------------------------------------------------------
  reflection hkl_unmerge_list::get_reflection(const int& jref) const
  //  unconditional                  ^^^^^^^^^^^^
  {
    if (status != PREPARED && status != SUMMED) 
      Message::message(Message_fatal(
                  "hkl_unmerge_list::get_reflection - not PREPARED or SUMMED") );
    NextRefNum = jref;      // record current reflection
    return refl_list[jref];
  }
  //--------------------------------------------------------------
  // Return reflection for given (reduced) hkl
  // returns index number or -1 if missing
  int hkl_unmerge_list::get_reflection(reflection& refl, const Hkl& hkl) const
  {
    if (status != SUMMED) 
      Message::message(Message_fatal(
		"hkl_unmerge_list::reflection - not SUMMED") );
    if (!is_hkl_lookup) {MakeHklLookup();}  // make lookup table if needed
    int idx = hkl_lookup.index_of(hkl.HKL());
    if (idx >= 0) {refl = refl_list[idx];}
    NextRefNum = idx;
    return idx;  // -1 if missing
  }     
  //--------------------------------------------------------------
  void hkl_unmerge_list::replace_reflection(const reflection& refl)
  // replace current reflection with updated version
  {
    refl_list.at(NextRefNum) = refl;
    refl_list[NextRefNum].CountNValid(); // update number of valid observations
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::replace_observation(const observation& obs,
					     const int& lobs) 
  // replace lobs'th observation in current reflection with updated version
  {
    refl_list.at(NextRefNum).replace_observation(obs, lobs);
    refl_list[NextRefNum].CountNValid(); // update number of valid observations
  }
  //--------------------------------------------------------------
  bool hkl_unmerge_list::ComparePartOrder::operator()(const observation_part * part1,
                                                      const observation_part * part2)
    // Return true if part1 is before part2 in sort order
    // Comparison function for sorting on H,K,L,M/ISYM,BATCH
  {
    // Compare most significant keys first
    if (part1->hkl().h() < part2->hkl().h())
      return true;
    else if (part1->hkl().h() > part2->hkl().h())
      return false;

    if (part1->hkl().k() < part2->hkl().k())
      return true;
    else if (part1->hkl().k() > part2->hkl().k())
      return false;

    if (part1->hkl().l() < part2->hkl().l())
      return true;
    else if (part1->hkl().l() > part2->hkl().l())
      return false;

    int M1 = Min(part1->Npart(), 2);
    int M2 = Min(part2->Npart(), 2);
    if (M1 < M2)
      return true;
    else if (M1 > M2)
      return false;

    if (part1->isym() < part2->isym())
      return true;
    else if (part1->isym() > part2->isym())
      return false;

    if (part1->batch() < part2->batch())
      return true;
    else if (part1->batch() > part2->batch())
      return false;

    return false;
  }
  

  //--------------------------------------------------------------
  void hkl_unmerge_list::sort()
    //                   ^^^^
    // Sort Observation part list
  {
    std::sort(obs_part_pointer.begin(), obs_part_pointer.end(), ComparePartOrder());
    status = SORTED;
  }
  //--------------------------------------------------------------
  int hkl_unmerge_list::change_symmetry(const hkl_symmetry& new_symm,
                                         const ReindexOp& reindex_op,
					 const bool& AllowFractIndex)
    // Change symmetry in all internal lists to new_spgp,
    // ie
    // 1. for each obs_part, get original indices 
    // 2. reindex (change basis)
    // 3. rereduce to new asymmetric unit
    // 4. flag as RAWLIST (unsorted)
  //
  // If AllowFractIndex true, allow discarding of fractional index
  // observations after reindexing, otherwise this is a fatal error
  {
    Hkl hkl, hkl_reindex, hkl_new;
    int new_isym;
    int NfractIdx = 0;

    // all Batches, fix up cell constraint flags
    std::vector<int> lbcell = new_symm.CellConstraint();
    ASSERT (lbcell.size() == 6);
    for (int ib=0; ib<nbatches; ib++)
      {batches[ib].SetCellConstraint(lbcell);}

    bool reindex = true;
    if (reindex_op.IsIdentity()) reindex = false;

    // change_basis unit cell
    if (reindex) {
      // Overall cell
      Scell tmpcell = averagecell.change_basis(reindex_op);
      averagecell = tmpcell;
      // all Batch info
      for (int ib=0; ib<nbatches; ib++)
	batches[ib].change_basis(reindex_op);
      // All datasets
      for (int id=0;id<ndatasets;id++)  {                             
	datasets[id].change_basis(reindex_op);
      }
      // Accumulate total reindexing
      totalreindex = totalreindex * reindex_op;
    }
    //^
    //^    std::cout << "**** Change symmetry from  " << refl_symm.symbol_xHM() << " to "
    //^	      << new_symm.symbol_xHM() << "  reindex " << reindex_op.a  s_hkl()
    //^	      << "\n";

    for (size_t i = 0; i < N_part_list; i++) {  // loop all raw observations 
      // Original indices
      hkl = refl_symm.get_from_asu(find_part(i).hkl(), find_part(i).isym());
      
      if (reindex)  {
	// returns false if non-integral indices
	bool HklOK = hkl.change_basis(hkl_reindex, reindex_op);
	if (!HklOK) {
	  NfractIdx++;
	  obs_part_pointer[i] = NULL;  // Clear pointer
	  continue;
	}
	hkl_new = new_symm.put_in_asu(hkl_reindex, new_isym);
      } else {
	hkl_new = new_symm.put_in_asu(hkl, new_isym);
      }
      //^
      //	std::cout << find_part(i).hkl().format() << " " << find_part(i).isym()
      //		  << " " << hkl.format()
      //		  << " " << hkl_new.format() << " " << new_isym
      //		  << " " << find_part(i).batch() << "\n";
      //^-
      find_part(i).set_hkl(hkl_new);
      find_part(i).set_isym(new_isym);
    }
    // Reset symmetry
    refl_symm = new_symm;
    // Store reindex operator to get back
    refl_symm.set_reindex(reindex_op.inverse());

    if (NfractIdx > 0) {
      // Some fractional indices have been found & discarded
      // Is this allowed?
      if (!AllowFractIndex) {
	Message::message(Message_fatal
	    ("hkl_unmerge_list::change_symmetry: illegal fractional indices generated by reindex operator") );
      }
	// Pack down pointer list
      int j = 0;
      for (size_t i=0;i<N_part_list;i++)  {
	if (obs_part_pointer[i])
	  {obs_part_pointer[j++] = obs_part_pointer[i];}
      }
      N_part_list = j;
      obs_part_pointer.resize(N_part_list);
    }
    status = RAWLIST;
    return NfractIdx;
  }
  //--------------------------------------------------------------
  int hkl_unmerge_list::prepare()
    // sort organise partials as required
    // returns number of reflections
  {
    if (status == EMPTY) 
      Message::message(Message_fatal(
                                     "hkl_unmerge_list::prepare - EMPTY") );
    if (N_part_list == 0) {
      Message::message(Message_fatal(
                     "hkl_unmerge_list::prepare  No observations in list") );
    }
    int n = Nref;
    if (status == RAWLIST) {
      // sort
      sort();
    }
    if (status == SORTED) {
      // organise
      n = organise();
    }
    if (status == ORGANISED) {
      // set partials
      partials();
    }
    // Reset resolution bin defaults using number of reflections
    ResolutionRange.SetRange(Nref);
    ResoLimRange.SetRange(Nref);
    return n;
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::SetPoles(const int& pole)
  // Set up poles for Absorption
  // maybe should be done by run?
  {
    for (int ib=0; ib<nbatches; ib++) {
      batches[ib].SetPole(pole);
    }
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::CalcSecondaryBeams(const int& pole)
  // Calculate all secondary beam directions, in chosen frame, also
  // diffraction vectors d*vec
  // On entry:
  //  pole  =  0 SECONDARY  camera frame
  //       !=  0 ABSORPTION, crystal frame = 1,2,3 for h,k,l, 
  //        = -1 unspecified, use closest reciprocal axis for each run
  {
    SetPoles(pole);
    reflection this_refl;
    observation this_obs;
    DVect3 sPhi;
    FVect3 fsPhi;

    // * * * * Loop reflections
    for (size_t j=0;j<refl_list.size();++j) {
      for (int iobs=0;iobs<refl_list[j].num_observations();++iobs) { // loop observations
	this_obs = refl_list[j].get_observation(iobs);
	int batchNum = this_obs.Batch();
	std::pair<float,float> thphi =
	  CalcSecondaryBeamPolar(batchNum, this_obs.hkl_original(),
				 this_obs.phi(), sPhi);
	this_obs.StoreS2(thphi.first, thphi.second);
	for (int i=0;i<3;++i) {fsPhi[i] = sPhi[i];}  // float not double for storage
	this_obs.StoreS(fsPhi);
	refl_list[j].replace_observation(this_obs);
      }
    }
  }
  //--------------------------------------------------------------
  std::pair<float, float> hkl_unmerge_list::CalcSecondaryBeamPolar
  (const int& batchNum, const Hkl& hkl_original, const float& phi,
   DVect3& sPhi) const
  // Calculate secondary beam directions
  //
  // On entry
  //  batchNum       batch number
  //  hkl_original   original hkl
  //  phi            incident beam rotation, degrees
  //
  // Returns
  //  sPhi      diffraction vector at actual phi position 
  //  pair(thetap, phip)   secondary beam directions, radians
  {
    DVect3 s2 = CalcSecondaryBeam(batchNum, hkl_original, phi, sPhi);
    // z is polar axis
    double r = sqrt(s2[0]*s2[0] + s2[1]*s2[1]);
    // thetap in range 0 -> pi
    float thetap = atan2(r, s2[2]);
    // phip 
    float phip = 0.0;
    if (r != 0.0) {
      phip = atan2(s2[1], s2[0]);
    }
    return std::pair<float, float>(thetap,phip);
  }
  //--------------------------------------------------------------
  DVect3 hkl_unmerge_list::CalcSecondaryBeam
  (const int& batchNum, const Hkl& hkl_original, const float& phi,
   DVect3& sPhi) const
  //
  // Calculate secondary beam directions
  //
  //  geometry calculations mostly done in Batch class
  //
  //     s = [E1] [E2] [E3] [UB] h
  //  
  //  1) outer rotation, 1-axis or omega scan
  //     s = [Phi/Omega] [D] [UB] h
  //     [D] = [E1] [E2] [E3] = [Omega0][Chi/Kappa][Phi]
  //
  //  2) 3-axis, phi scan
  //     s = [Omega] [Chi/Kappa] [Phi] [D] [UB] h
  //     [D] = [Phi0]      [E1E2] = [Omega] [Chi/Kappa]
  //
  //  s(phi0 frame) = [DUB] h
  //      
  // On entry
  //  irun           run serial number
  //  hkl_original   original hkl
  //  phi            incident beam rotation, degrees
  //
  // Returns
  //  sPhi      diffraction vector at actual phi position (camera frame)
  //  DVect3    secondary beam directions, direction cosines
  {
    int ib = batch_lookup.lookup(batchNum);  // batch serial
    // diffraction vector at phi = 0  zero rotation angle frame
    //    s(r0) = DUB h (dimensionless reciprocal lattice units)
    DVect3 sr0 = batches[ib].HtoSr0(hkl_original);
    // Source vector in zero rotation angle frame (unit vector)
    DVect3 s0r0 = batches[ib].SrtoSr0(batches[ib].Source(), phi);
    // secondary beam s2(r0) = s(r0) - s0(r0) in zero rotation angle frame
    DVect3 s2 = sr0 - s0r0;
    if (batches[ib].Pole() > 0) {
      // ABSORPTION, crystal frame = 1,2,3 for h,k,l
      // back rotate s2 into permuted crystal frame
      //   s2' = [P][DU]^-1 s2  
      s2 = batches[ib].Sr0toP(s2);
    }
    // s(r) = [R][D][U][B]h  camera frame
    sPhi =  batches[ib].HtoSr(hkl_original, phi);
    //^ sanity check
    //    float wvl = batches[ib].Wavelength();
    //    std::cout << "\ns(r)   = " << sPhi.format() << " d " << wvl/sqrt(sPhi*sPhi)<<"\n";
    //    std::cout << "s0(r0) = " << s0r0.format() <<"\n";
    //    DVect3 s2r = sPhi - batches[ib].Source();
    //    std::cout << "s2(r)  = " << s2r.format() << " " << sqrt(s2r*s2r)<<"\n";
    //    std::cout << "s2(r0) = " << s2.format() <<" "<< sqrt(s2*s2) << "\n";
    //    // s2(r0) = 
    //    DVect3 s2r0 = batches[ib].SrtoSr0(s2r, phi);
    //    std::cout << "s2(r0) = " << s2r0.format() << " " << sqrt(s2r0*s2r0) <<"\n";
    //^-
    return s2.unit();  // unit vector along s2
  
    /*

    // Phi rotation matrix
    DMat33 Phi = clipper::Rotation
      (clipper::Polar_ccp4(0.0,0.0,+clipper::Util::d2rad(phi))).matrix();
    // Inverse Phi rotation matrix
    DMat33 PhiInv = Phi.transpose();
    // diffraction vector at phi = 0   s = DUB h
    //    dimensionless reciprocal lattice units
    DVect3 sPhi0 = batches[ib].DUBmat() * hkl_original.real();
    // s0 (Phi=0) = [E1E2]^-1  [Phi]^-1 s0
    DVect3 s0 = batches[ib].Source();
    DVect3 s00 = s0;
    if (batches[ib].PhiScan()) {
      // 3-axis Phi scan
      s00 = batches[ib].E1E2().transpose() * s00;
      // total goniostat rotation
      Phi = batches[ib].E1E2() * Phi;
      //^
      std::cout << "Phi\n"<<Phi.format()
		<<"\nDet = " << Phi.det() <<"\n";
      //^-
    }
    s00 = PhiInv * s00;
    // secondary beam s2 = s - s0(0)  (phi = 0 camera frame)
    DVect3 s2 = sPhi0 - s00;
    if (batches[ib].Pole() > 0) {
      // ABSORPTION, crystal frame = 1,2,3 for h,k,l
      // back rotate s2 into permuted crystal frame
      //   s2' = [DUP]^-1 s2  
      s2 = batches[ib].DUPinv() * s2;
    }
    // Rotate s(phi=0) to s: s = [Phi] s(phi=0)
    sPhi = Phi * sPhi0;    // returned
    */
  }
//--------------------------------------------------------------
  void hkl_unmerge_list::ResetObsAccept (ObservationFlagControl& ObsFlagControl)
  // Reset observation accepted flags to allow for acceptance of
  // observations flagged as possible errors
  // Counts observations reclassified them in ObsFlagControl
  {
    ObsFlagControl.Clear();  // clear counts

    int Next = 0;
    while (Next < Nref) {
      // Process all reflections unconditionally
      reflection this_refl = get_reflection(Next++);
      this_refl.ResetObsAccept(ObsFlagControl);
      replace_reflection(this_refl);
    }
  }
//--------------------------------------------------------------
  void hkl_unmerge_list::ResetReflAccept()
  // Reset reflection accepted flags to accept all
  //  (subject to resolution checks etc)
  {
    int Next = 0;
    while (Next < Nref) {
      // Process all reflections unconditionally
      reflection this_refl = get_reflection(Next++);
      if (this_refl.Status() != 0) {
	this_refl.SetStatus(0);  // accept
	replace_reflection(this_refl);
      }
    }
  }
  //--------------------------------------------------------------
  Range hkl_unmerge_list::UpdatePolarisationCorrections(const bool& Total,
					  const double& polarisationfactor)
  // Update polarisation corrections for all parts
  // If Total == true, then apply complete correction
  //   else assume the unpolarised correction is already applied, apply
  //   additional correction for polarised incident beam
  // polarisationfactor if fraction polarised, = 0 for unpolarised, ~ 0.9 for synchrotrons
  //
  // Returns range of correction factors
  //
  // see Kahn, Fourme, Gadet, Janin, Dumas, & Andre,
  // J. Appl. Cryst. (1982). 15, 330-337
  {
    Range PFrange;
    Hkl hkl_original;
    observation_part part;
    // loop parts
    for (size_t i = 0; i < N_part_list; i++) {  // loop all raw observations
      part = find_part(i);
      // Original indices
      hkl_original = refl_symm.get_from_asu(part.hkl(), part.isym());
      int ib = batch_lookup.lookup(part.batch());  // batch serial
      // s(r) = [R][D][U][B]h  camera frame, reciprocal lattice units
      DVect3 sPhi =  batches[ib].HtoSr(hkl_original, part.phi());
      double z = sPhi[2]; // z coordinate
      double sinSqtheta = 0.25 * (sPhi * sPhi);    // |s| = 2 sin theta; sin^2 theta = 0.25*|s|^2
      double cos2theta = 1. - 2.0 * sinSqtheta;    // cos 2theta = cos^2 theta - sin^2 theta
						   //  = 1 - 2 sin^2 theta
      double cosSq2theta = cos2theta * cos2theta;  // cos^2 2theta
      double sinSq2theta = 1.0 - cosSq2theta;      // sin^2 2theta
      double cosrho = z/sqrt(sinSq2theta);         // cos rho = z/sine 2theta
      // P0 = 1/2 [ 1 + cos^2 2 theta]
      double P0 = 0.5*(1.0 + cosSq2theta);  // unpolarised part
      // P' = 1/2 (-Xsi') cos 2rho sin^2 2theta
      double PP = 0.5 * polarisationfactor * (2.0*cosrho*cosrho - 1.0) * sinSq2theta;
      double PolFac = 1.0/P0;
      if (!Total) {
	// P0 already applied, so just correct it
	//  complete PolFac = P0 - P', negative sign because polarisationfactor should be negative
	PolFac = P0/(P0 - PP); // inverse
      }
      // Apply it, dividing intensities
      find_part(i).StoreIsigI(part.I_sigI().scaleIs(PolFac));
      find_part(i).StoreIsigIpr(part.I_sigIpr().scaleIs(PolFac));

      // stored LP factor is 1/LP, update it
      find_part(i).StoreLP(part.LP()*PolFac);
      PFrange.update(PolFac);
      //^
      //      std::cout << hkl_original.format() <<" "<< P0 <<" "<< PP
      //		<<" "<< PolFac <<" " << clipper::Util::rad2d(acos(cosrho))<<"\n";
      //^-
    } // end loop parts
    //^
    //    std::cout << "Range of polarisation corrction factors "
    //	      << PFrange.min() <<" " << PFrange.max() <<"\n";
    return PFrange;
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::dump_reflection(const Hkl& hkl) const
  // dump given reflection, for debugging
  {
    reflection this_refl;
    int idx = get_reflection(this_refl, hkl);
    if (idx < 0) {
      std::cout << "Reflection missing " << hkl.format() << "\n";
    } else {
      int nobs = this_refl.num_observations();
      std::string acc = " accepted";
      std::string rej = " rejected";
      std::string refacc = (this_refl.Status()==0) ? acc : rej;
      std::cout << "** Reflection " << hkl.format()
		<< " resolution " << sqrt(1./this_refl.invresolsq())
		<< " number observations " << nobs
		<< " valid " << this_refl.NvalidObservations()
		<< refacc << "\n";
      int Next = 0;
      while (Next < nobs) {
      // all observations unconditionally
	observation obs = this_refl.get_observation(Next++);
	std::string obsacc = obs.IsAccepted()  ? acc : rej;
	std::string fp =  obs.IsFull() ? " full" : " partial";
	std::cout << " * Observation " << obs.hkl_original().format()
		  <<  obsacc
		  << " run " << obs.run() << " dataset " << obs.datasetIndex()
		  << " batch " << obs.Batch() << " isym " << obs.Isym()
		  << fp << "\n";
      }
    }
  }
//--------------------------------------------------------------
// Copy constructor throws exception unless object is EMPTY
  hkl_unmerge_list::hkl_unmerge_list(const hkl_unmerge_list& List)
  {
    if (List.status != EMPTY) {
      Message::message(Message_fatal
		       ("hkl_unmerge_list: illegal copy constructor"));
    }
    clear();
  }
//--------------------------------------------------------------
  // Copy operator throws exception unless object is EMPTY
  hkl_unmerge_list& hkl_unmerge_list::operator= (const hkl_unmerge_list& List)
  {
    if (List.status != EMPTY) {
      Message::message(Message_fatal
		       ("hkl_unmerge_list: illegal copy operation"));
    }
    clear();
    return *this; 
  }

} // namespace 
