// hkl_unmerge.cpp
//
// Phil Evans August 2003
//
// Classes for unmerged hkl lists
//  at present just for the items used by Mosflm & Scala

// These data structures are primarily directed at Scala so are
// in its namespace

#include <algorithm>

#include <assert.h>
#define ASSERT assert

#include "hkl_unmerge.hh"
#include "scala_util.hh"
#include "string_util.hh"
#include "report_errors.hh"

// Clipper
#include <clipper/clipper.h>
#include "clipper/clipper-ccp4.h"


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
    is_latnum = false;
    is_lathkl = false;
    is_latinfo = false;
    is_latscale = false;
    n_latinfo = 0;
  }
  // ------------------------------------------------------------------
  void data_flags::print() const // for debugging
  {
    std::cout << "\nDataFlags:\n"
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
              << "\nis_time  " << is_time
              << "\nis_latnum " << is_latnum
              << "\nis_lathkl " << is_lathkl
              << "\nis_latinfo " << is_latinfo
              <<"\nn_latinfo " << n_latinfo
              << "\n\n";
  }
  // ------------------------------------------------------------------
  // set up flags for multilattice output, for maxNoverlap sets of columns
  // if scalecolumn true, include a column for the scale
  void data_flags::SetMultilatticeFlags(const int& maxNoverlap,
                                        const bool& scalecolumn)
  {
    is_latnum = true;
    is_latinfo = true;
    is_lathkl = false;
    n_latinfo = maxNoverlap;
    is_latscale = scalecolumn;
  }
  // -------------------------------------------------------------
  void data_flags::combineFlags(const data_flags& other)
  // set each flag to true if either are true
  {
    is_h |= other.is_h;
    is_k |= other.is_k;
    is_l |= other.is_l;
    is_misym |= other.is_misym;
    is_batch |= other.is_batch;
    is_I |= other.is_I;
    is_sigI |= other.is_sigI;
    // Optional
    is_Ipr |= other.is_Ipr;
    is_sigIpr |= other.is_sigIpr;
    is_fractioncalc |= other.is_fractioncalc;
    is_Xdet |= other.is_Xdet;
    is_Ydet |= other.is_Ydet;
    is_Rot |= other.is_Rot;
    is_Width |= other.is_Width;
    is_LP |= other.is_LP;
    is_Mpart |= other.is_Mpart;
    is_ObsFlag |= other.is_ObsFlag;
    is_BgPkRatio |= other.is_BgPkRatio;
    is_scale |= other.is_scale;
    is_sigscale |= other.is_sigscale;
    is_time |= other.is_time;
    is_latnum |= other.is_latnum;
    is_lathkl |= other.is_lathkl;
    is_latinfo |= other.is_latinfo;
    is_latscale |= other.is_latscale;
    n_latinfo = Max(n_latinfo, other.n_latinfo);
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
  // return true if COMBINE
  bool SelectI::Combine() {
    return ((selecticolflag > 0) && (imid > 0.001));
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
  //--------------------------------------------------------------
  LatticeIndexInfo::LatticeIndexInfo(const int& Latnum, const Hkl& jhkl,
                           const Rtype& Gscale)
  {
    init(Latnum, jhkl, Gscale);
  }
  //--------------------------------------------------------------
  void LatticeIndexInfo::init(const int& Latnum, const Hkl& jhkl,
                         const Rtype& Gscale)
  {
    latnum = Latnum;
    hkl = jhkl;
    gscale = Gscale;
  }
  //--------------------------------------------------------------
  // ******************  observation_part  *******************

  observation_part::observation_part(const Hkl& hkl_in,
                                     const int& isym_in, const int& batch_in,
                                     const Rtype& I_in, const Rtype& sigI_in,
                                     const Rtype& Ipr_in, const Rtype& sigIpr_in,
                                     const Rtype& Xdet_in, const Rtype& Ydet_in,
                                     const Rtype& phi_in, const Rtype& time_in,
                                     const Rtype& fraction_calc_in, const Rtype& width_in,
                                     const Rtype& LP_in,
                                     const int& Npart_in, const int& Ipart_in,
                                     const ObservationFlag& ObsFlag_in,
                                     const int& latnum_in,
                                     const std::vector<LatticeIndexInfo>& lathkl_in)
    :     hkl_(hkl_in),
          isym_(isym_in), batch_(batch_in),
          I_(I_in), sigI_(sigI_in),
          Ipr_(Ipr_in), sigIpr_(sigIpr_in),
          Xdet_(Xdet_in), Ydet_(Ydet_in), phi_(phi_in), time_(time_in),
          fraction_calc_(fraction_calc_in), width_(width_in),
          LP_(LP_in),
          Npart_(Npart_in), Ipart_(Ipart_in), ObsFlag_(ObsFlag_in),
          run_(1), latnum_(latnum_in), lathkl_(lathkl_in)
  {}
  //--------------------------------------------------------------
  std::string observation_part::format() const // for debugging, incomplete
  {
    std::string s = "observation_part: " + hkl_.format();
    s += " Isym "+ clipper::String(isym_) + " batch " + clipper::String(batch_);
    s += " Npart " + clipper::String(Npart_) + " Ipart " + clipper::String(Ipart_);
    if (latnum_ != 0) {
      s += "\n  Lattice " + clipper::String(latnum_) +"\n";
      for (size_t i=0; i<lathkl_.size(); i++) {
        s += "     Overlap " + clipper::String(lathkl_[i].latnum) +
          " " + lathkl_[i].hkl.format() +
          " gscale " + clipper::String(lathkl_[i].gscale) + "\n";
      }
    }
    s += "\n";
    return s;
  }
  //--------------------------------------------------------------
  // apply scale to both IsigIs
  void observation_part::ScaleIsigI(const Rtype& scale)
  {
    I_ *= scale;
    sigI_ *= scale;
    Ipr_ *= scale;
    sigIpr_ *= scale;
  }
  //--------------------------------------------------------------
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
                           const ObservationFlag& obsflag_in,
                           const int& latnum_in,
                           const std::vector<LatticeIndexInfo>& lathkl_in)
    : hkl_original_(hkl_in),
      isym_(isym_in), run_(run_in), datasetIndex_(datasetIndex_in),
      Npart_(Npart_in), part1(part1_in), batch_(0),
      totalfraction(TotFrac), part_flag(partialstatus_in), obs_flag(obsflag_in),
      gscale(1.0), vargscale(-1.0), latnum(latnum_in), lathkl_(lathkl_in)
  {
    // By default here reject if any flag set
    // This observation may be accepted later if the flags pass a conditional test
    // (see ResetObsAccept)
    if (! obs_flag.OK()) obs_status.SetObsFlag();
  }

  //--------------------------------------------------------------
  int observation::num_parts() const
  //               ^^^^^^^^
  {
    return Npart_;
  }
  //--------------------------------------------------------------
  int observation::num_parts_mpart() const
  //               ^^^^^^^^
  // Return number of parts as recorded in the MPART column of the first part
  // This will get the number of parts from observations where the parts have already
  // been summed
  {
    int npart = Npart_;
    if (npart == 1) {
      npart = get_part(0).Npart();
    }
    return npart;
  }
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
  // return scaled sigI
  Rtype observation::ksigI() const
  //
  {
    if (vargscale <= 0.0) {
      return sigI_/gscale;
    }
    Rtype scI = kI();
    return sqrt(sigI_*sigI_ + scI*scI*vargscale)/gscale;
  }
  //--------------------------------------------------------------
  IsigI observation::kI_sigI() const
  // return scaled I, sigI
  {
    if (vargscale <= 0.0) {
      {return IsigI(I_/gscale,sigI_/gscale);}
    }
    return IsigI(kI(), ksigI());
  }
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
  bool observation::hasIpr() const
  {
    // Return true if there is an IPR value
    // only need to test 1st or only part
    return (get_part(0).I_sigIpr().sigI() > 0.0);
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
      if (Ipr.sigI() <= 0.0) {
        Isum = Ic;  // no valid Ipr for this observation
      } else {
        Rtype Iraw = Ic.I();
        if (LP_ > 0.0) Iraw /= LP_;  // raw intensity back-corrected for LP
        Isum = SelectI::GetCombinedI(Iraw, Ic, Ipr);
      }
    } else if (SelectI::SelectIcolFlag() < 0) { // profile
      if (Ipr.sigI() > 0.0) {
        Isum = Ipr;
      } else {
        Isum = Ic;
      }
    } else {
      Isum = Ic;
    }
    I_ = Isum.I();
    sigI_ = Isum.sigI();
  }
  //--------------------------------------------------------------
  Rtype observation::width() const  //! return reflection width (degrees, from input)
  {
    // Average over parts
    if (Npart_ == 1) {
      return get_part(0).width();
    }
    // partial
    Rtype width = 0.0;
    for (int kpart = 0; kpart < Npart_; kpart++) { // loop parts
      width += get_part(kpart).width();
    }
    return width/Rtype(Npart_);
  }
  //--------------------------------------------------------------
  // return range of batches for this observation
  IntRange observation::BatchRange() const
  {
    int npart = num_parts();
    IntRange batchrange;
    for (int jp=0;jp<npart;++jp) {
      batchrange.update(get_part(jp).batch());
    }
    return batchrange;
  }
  //--------------------------------------------------------------
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
                              const std::vector<Dataset>& DataSets,
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
    nlattices = 0;
    nlatticesall = 0;

    filename = "";
    FileTitle = Title;
    historylines.clear();

    // store controls
    run_flags = controls.runs;
    run_set = 0;
    partial_flags = controls.partials;
    partial_set = true;
    datasets.clear();
    batches.clear();
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::setNoPartials(const bool& nopartials)
  {
    // set NoPartial flag
    bool currentNoPartial = partial_flags.noPartials();
    partial_flags.setNoPartials(nopartials);  // reset
    rfl_status currentstatus = status;
    if ((status == PREPARED) || (status == SUMMED)) {
      if (!nopartials && currentNoPartial) {
        // we are switching from noPartials to Partials, so re-organise
        status = ORGANISED;
        prepare();
        if (currentstatus == SUMMED) {
          sum_partials();
        }
      }
    }
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::clear()
  // Clear out list ready for new init
  {
    status = EMPTY;
    init("Empty list",0,hkl_symmetry(),all_controls());
  }
  //--------------------------------------------------------------
  //! append to MTZ history
  void hkl_unmerge_list::addHistory(const std::vector<std::string>& history)
  {
    historylines.insert(historylines.end(), history.begin(), history.end());
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::SetBatchList()  // private
  {
    int idataset;
    for (size_t i=0;i<batches.size();i++) {
      if (in_datasets(batches[i].PXDname(), datasets, idataset)) {
        datasets[idataset].add_batch(batches[i].PXDname(), batches[i].num());
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
      ///      datasets[j].SetCell(avbcell[j]);
      datasets[j].SetMosaicity(averageMosaicity[j]);
      float wvl = averageWavelength[j];
      datasets[j].SetCellWavelength(avbcell[j], wvl);
    }
    averagecell = AverageDsetCell(datasets);
  }
  //--------------------------------------------------------------
  Scell hkl_unmerge_list::AverageOtherBatchData(const std::vector<Batch>& batches,
                                                const int& ndatasets,
                                                std::vector<float>& averageMosaicity,
                                                std::vector<float>& averageWavelength,
                                                std::vector<Scell>& avbcell) const
  {
    // Average batch cells for each dataset
    //  (like AverageBatchData only for other data)
    // also mosaicity & wavelength
    // Returns overall average cell
    // Sets averageMosaicity, averageWavelength and avbcell

    avbcell = AverageBatchCell(batches, ndatasets,
                               averageMosaicity, averageWavelength);
    return AverageDsetCell(datasets);
  }
  //--------------------------------------------------------------
  // Store datasets & batch info following previous call to init
  void hkl_unmerge_list::StoreDatasetBatch(const std::vector<Dataset>& DataSets,
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
  void hkl_unmerge_list::AddDatasetBatch(const std::vector<Dataset>& Datasets,
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
      ReportErrors::printFatalError
        ("hkl_unmerge_list::OffsetBatchNumbers - no runs set");
    }
    if (run_set == 0) set_run();  // setup runs in part list if not already done
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
      if (in_datasets(batches[i].PXDname(), datasets, idataset)) {
        datasets[idataset].add_batch(batches[i].PXDname(), batches[i].num());
      }
    }
    // Reset all batch numbers in observation part list
    observation_part part;
    for (size_t i = 0; i < N_part_list; i++) {  // loop all raw observations
      part = find_part(i);
      find_part(i).set_batch(part.batch() + runOffsets[part.run()]);
    }
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::OffsetLatticeNumbers(const int& latticeoffset)
  // Apply offset to lattice numbers
  {
    if (latticeoffset == 0) {return;} // don't bother if 0

    //std::cout <<"OffsetLatticeNumbers " << latticeoffset <<std::endl; //^
    if (!run_flags.Set()) {
      ReportErrors::printFatalError
        ("hkl_unmerge_list::OffsetLatticeNumbers - no runs set");
    }
    if (run_set == 0) set_run();  // setup runs in part list if not already done

    for (size_t irun=0;irun<runlist.size();irun++) {   // Loop runs
      runlist[irun].SetLatticeNumber(runlist[irun].LatticeNumber() + latticeoffset);
    }
    // Offset lattice numbers in batches
    for (size_t ib=0; ib<batches.size(); ib++) {
      batches[ib].SetLatticeNumber(batches[ib].LatticeNumber() + latticeoffset);
    }

    // Reset all lattice numbers in observation part list
    observation_part part;
    for (size_t i = 0; i < N_part_list; i++) {  // loop all raw observations
      part = find_part(i);
      find_part(i).set_latnum(part.latnum() + latticeoffset);
      std::vector<LatticeIndexInfo> lathkl = part.lathkl();
      if (lathkl.size() > 0) {
        for (size_t j=0; j<lathkl.size(); j++) {
          lathkl[j].latnum += latticeoffset;
        }
        find_part(i).set_lathkl(lathkl);
      }
    }
  }
  //--------------------------------------------------------------
  //! set range of lattice numbers either as main lattice or secondary
  void hkl_unmerge_list::SetLatticeNumberRange(const IntRange& latticenumberRange)
  {
    latticenumberrange = latticenumberRange;
    nlatticesall = latticenumberrange.AbsRange()+1;
  }
  //--------------------------------------------------------------
  //! range of lattice numbers as main lattice
  void hkl_unmerge_list::SetMainLatticeNumberRange(const IntRange& mainlatticenumberRange)
  {
    mainlatticenumberrange = mainlatticenumberRange;
    nlattices = mainlatticenumberrange.AbsRange()+1;
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::RejectBatch(const int& ibatch, const bool& fromrun)
  // Mark batch number ibatch as not accepted
  // Data records are not changed, runlist is updated if fromrun true
  {
    RejectBatchSerial(batch_lookup.lookup(ibatch), fromrun);
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::RejectBatchSerial(const int& jbat, const bool& fromrun)
  // Mark batch with serial number jbat as not accepted
  // Data records are not changed, runlist is updated if fromrun true
  {
    batches.at(jbat).SetAccept(false);
    if (fromrun) {
      int irun = batch(jbat).RunIndex();
      int batchnum = batch(jbat).num();
      runlist.at(jbat).SetBatchAccept(batchnum, false);
    }
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::AppendFileName(const std::string& Name)
  {
    if (filename != "") filename += " + ";
    filename += Name;
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::MergeDatasetLists(const std::vector<Dataset>& otherDatasets,
                                           const std::vector<Batch>& otherBatches)
  // Are new datasets (in otherDatasets) the same as any old ones? Append new ones to list
  // For each dataset from otherDataset list, store equivalent dataset
  // index in present list, if it is the same dataset. Return index list
  // Append batches with dataset references
  {
    //    const double Tolerance = 1.0;

    // Average cell, mosaicity & wavelength over all batches for each dataset
    // & store in dataset
    std::vector<float> averageMosaicity;  // size otherDatasets.size()
    std::vector<float> averageWavelength;
    std::vector<Scell> avbcell;
    Scell otheraveragecell =
      AverageOtherBatchData(otherBatches, otherDatasets.size(),
                            averageMosaicity, averageWavelength, avbcell);

    std::vector<int> DtsIndex(otherDatasets.size());
    int ndts = datasets.size(); // number of current datasets
    int kd = ndts - 1;  // index for new datasets appended to old ones
    int id = -1;        // dataset ID for new datasets, largest current id
    for (int j=0;j<ndts;j++) {id = datasets[j].MaxID(id);}

    for (size_t i=0;i<otherDatasets.size();i++) { // loop new (other) datasets
      DtsIndex[i] = -1;
      for (int j=0;j<ndts;j++) { // loop current datasets
        if (otherDatasets[i] == datasets[j]) {
          DtsIndex[i] = j;  // i'th "Other" dataset has same names as j'th
          // add new cell and wavelength into list (put into 1st Xdataset in Dataset)
          datasets[j].AddCellWavelength(avbcell[i],
                                        averageWavelength[i]);
          // Check for similar unit cell & wavelength
          //      if (!datasets[j].cell().equalsTol(otherDatasets[i].cell(), Tolerance)) {
          //        Message::message(Message_warn
          //          ("\nWARNING: Datasets with same name have different unit cells, 1st one used"));
          //      }
          break;
        } else {
          // Add any matching Xdatasets into current one
          bool added = datasets[j].AddDataset(otherDatasets[i]);
          if (added) {
            DtsIndex[i] = j;  // i'th "Other" dataset has same names as j'th
            break;
          }
        }
      }
      // DtsIndex[i] >= 0 (= j) if the i'th new dataset has been added into the j'th old one
      if (DtsIndex[i] < 0) {
        // New dataset
        Dataset OtherDataset = otherDatasets[i];
        DtsIndex[i] = ++kd;  // new index for i'th dataset
        id++;
        // id returned incremented if OtherDataset contains > 1 Xdataset
        id = OtherDataset.StoreSetID(id);
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
        ReportErrors::printFatalError
          ("Non-unique batch number"+clipper::String(otherBatches[i].num()));
      }
      // Update dataset index
      int otherIndex = OtherBatch.datasetindex(); // old dataset index
      int idts = DtsIndex[otherIndex];  // new index
      ASSERT (idts >= 0);
      OtherBatch.datasetindex() = idts;  // store new index
      PxdName pxdname = OtherBatch.PXDname();  // name
      int ID = datasets[idts].GetID(pxdname);  // datasetID for this Xdataset
      OtherBatch.DatasetID() = ID;
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
      ReportErrors::printFatalError
        ("hkl_unmerge_list::Append - cannot append to self");
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
      ReportErrors::printWarning
        ("\nWARNING: Cannot combine reflection lists with different symmetry", "WarningMessage",false);
      istat = 1;
      return istat;
    }

    // Copy all parts
    for (size_t i=0;i<OtherList.N_part_list;i++) {
      obs_part_list.push_back(OtherList.find_part(i));
      N_part_list++;
    }
    if (obs_part_list.size() != N_part_list)
      ReportErrors::printFatalError("hkl_unmerge_list::close_part - Wrong length list");
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

    SetBatchList();
    for (int j=0;j<ndatasets;j++) {
      datasets[j].setAverageCell(AverageBatchCellforDataset(batches, j));
    }

    is_hkl_lookup = false;
    // Set up run definitions
    run_set = 0;
    SetUpRuns();

    // Flag column present if in either file
    dataflags.combineFlags(OtherList.dataflags);

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
    maxlatnum = 0;
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
                                    const ObservationFlag& ObsFlag,
                                    const int& latnum,
                                    const std::vector<LatticeIndexInfo>& lathkl)
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
                                             Npart, Ipart, ObsFlag,
                                             latnum, lathkl));
    // Don't set pointer list obs_part_pointer until end (in close_part)
    // in case vector gets extended
    N_part_list++;
    if (latnum > 0) {
      // get maximum lattice number
      maxlatnum = Max(latnum, maxlatnum);
    }
  } // store_part
  //--------------------------------------------------------------
  // Store a raw observation part in list
  void hkl_unmerge_list::store_part(const observation_part& part)
  {
    obs_part_list.push_back(part);
    // Don't set pointer list obs_part_pointer until end (in close_part)
    // in case vector gets extended
    N_part_list++;
    if (part.latnum() > 0) {
      // get maximum lattice number
      maxlatnum = Max(part.latnum(), maxlatnum);
    }
  } // store_part
  //--------------------------------------------------------------
  // Close raw observation part list, return number of parts
  int hkl_unmerge_list::close_part_list(const ResoRange& RRange,
                                        const bool& Sorted)
  {
    if (obs_part_list.size() != N_part_list) {
      ReportErrors::printFatalError("hkl_unmerge_list::close_part - Wrong length list");}
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
    if (explicitruns) {
      SetRunInput(); // reset runs from input specifications
    }
    CheckAllRuns(); // finish run specification
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::ResetRuns(const run_controls& Runcontrols)
  {
    run_flags = Runcontrols;
    SetRunInput(); // reset runs from input specifications
    CheckAllRuns(); // finish run specification
    run_set = 0;
    status = SORTED;
    prepare();
    sum_partials();
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::AutoSetRun()
  // Divide batches up into runs
  {
    Run ThisRun;
    int PreviousBatNum = -1;
    float PreviousPhi = 0.0;
    clipper::Rotation PreviousUinv;
    // Limit for monitoring large orientation change, 2 degrees
    double rotlim = clipper::Util::d2rad(3.0);
    float phioffset;
    double tolerance = 0.01;
    int offset = 0;
    int filenum = -1;
    int batNgap;
    float delPhi = 0.0;
    int latnum = 0;
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
        ThisRun = Run(batch(ib).datasetindex(), batch(ib).DatasetID(),
                      batch(ib).PXDname());
        // Store run index in dataset
        if (batch(ib).Accepted()) {
          datasets[batch(ib).datasetindex()].AddRunIndex(batch(ib).PXDname(), runlist.size());
        }
        offset  = batch(ib).BatchNumberOffset();
        filenum = batch(ib).FileNumber();
        phioffset = 0.0;
        batches[ib].OffsetPhi(phioffset);
        if (batches[ib].IsTimePhi()) {
          // if time == phi offset
          batches[ib].OffsetTime(phioffset);
        }
      } else { // not 1st batch
        // Compare this batch with last one
        bool newgroup = false;
        // Conditions for being in the same group (run):
        // same dataset ID (dataset & crystal)
        //      if (batch(ib).datasetindex() != ThisRun.DatasetIndex()) {newgroup = true;}
        if (batch(ib).DatasetID() != ThisRun.DatasetID()) {newgroup = true;}
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
                //              std::cout << "AutoSetRun: not newgroup, batch, gap, phioffset "
                //                        << batch(ib).num() << " " << gap <<" " << phioffset << "\n";
                //^-
              }
            }
          }
        }
        // Check for large change in orientation
        double rot = (clipper::Rotation(batch(ib).Umat()) * PreviousUinv).abs_angle();
        if (rot > rotlim) {
          newgroup = true;
        }
        if (newgroup) {
          //^
          //      std::cout << "AutoSetRun: newgroup, batch, batNgap, gap "
          //                << batch(ib).num() << " " << batNgap << " " << gap << "\n";
          //^-
          // Store completed run as long as it has some accepted batches
          if (nbatAccepted > 0) {
            ThisRun.BatchNumberOffset() = offset;
            ThisRun.FileNumber() = filenum;
            ThisRun.SortList();
            ThisRun.RunNumber() = runlist.size()+1;
            Batch bat = batch(ib-1);
            latnum = batch(ib-1).LatticeNumber();
            ThisRun.SetLatticeNumber(latnum);
            runlist.push_back(ThisRun);
          }
          // Start new group, store dataset index
          ThisRun = Run(batch(ib).datasetindex(), batch(ib).DatasetID(),
                        batch(ib).PXDname());
          // Store run index
          if (batch(ib).Accepted()) {
            datasets[batch(ib).datasetindex()].AddRunIndex(batch(ib).PXDname(), runlist.size());
          }
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
          //      std::cout << "ib, offset " << ib <<" "<< phioffset <<"\n";
          //^-
          if (batches[ib].IsTimePhi()) {
            // if time == phi offset
            batches[ib].OffsetTime(phioffset);
          }
          if (std::abs(phioffset) > 0.1) {IsPhiOffset = true;}
          // Sanity checks
          ASSERT (offset == batch(ib).BatchNumberOffset());
          //^           ASSERT (filenum == batch(ib).FileNumber());
        }
      } // not 1st batch
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
      PreviousUinv = clipper::Rotation(batch(ib).Umat().inverse());
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
      latnum = batches[batch_lookup.lookup(ThisRun.BatchList()[0])].LatticeNumber();
      ThisRun.SetLatticeNumber(latnum);
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
  void hkl_unmerge_list::SetRunByFile()
  // Set up runs from input specification
  {
    // Clear run specs
    runlist.clear();
    // Clear dataset run indices
    for (size_t id=0;id<datasets.size();id++) {
      datasets[id].ClearRunList();
    }
    int currentfilenumber = -1; // initial
    Run ThisRun;
    int runindex = -1;

    for (size_t i=0;i<batches.size();++i) { //  loop batches
      int filenumber = batches[i].FileNumber();
      if (filenumber != currentfilenumber) {
        if (currentfilenumber >= 0) {
          // not first, so close this one
          ThisRun.SortList();
          runlist.push_back(ThisRun);
        }
        // start new run
        ThisRun = Run(batches[i].datasetindex(), batches[i].DatasetID(),
                          batches[i].PXDname());
        runindex = runlist.size();
        ThisRun.BatchNumberOffset() = batches[i].BatchNumberOffset();
        ThisRun.FileNumber() = filenumber;
        ThisRun.RunNumber() = runindex + 1;
        int latnum = batches[i].LatticeNumber();
        ThisRun.SetLatticeNumber(latnum);
        currentfilenumber = filenumber;
      }
      ThisRun.AddBatch(batches[i].num(), batches[i].Accepted());
      if (batches[i].Accepted()) {
        // Store run index in batch: this = current size of runlist
        batches[i].SetRunIndex(runindex);
      } else {
        // Store null run index in batch
        batches[i].SetRunIndex(-1);
      }

    }  // end loop batches

    if (ThisRun.Nbatches() > 0) {
      ThisRun.SortList();
      runlist.push_back(ThisRun);  // close last one
    }
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  void hkl_unmerge_list::SetRunInput()
  // Set up runs from input specification
  {
    if (run_flags.Auto()) {return;}  // default AUTO
    if (run_flags.Byfile()) {SetRunByFile(); return;}  // BYFILE

    // set run_controls status
    run_flags.SetStatus(-2);

    // batch ranges
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
      //      std::cout << "Accepted batches " << ib <<  " "<<batches[ib].num()
      //                <<" "<<batches[ib].Accepted()<<"\n";
      // Store null run index in batch
      batches[ib].SetRunIndex(-1);
      maxbatchnum = Max(maxbatchnum, batches[ib].num());
    }
    int latnum = 0;
    Run this_run;

    for (int irun=0;irun<nruns;++irun) {    // Loop runs
      int runnum = runnumberlist[irun];  // run number
      bool firstbatchinrun = true;
      int nbatAccepted = 0;
      int ib0 = -1;
      //^^
      //      std::cout << "\n\nRun " <<runnum<<"\n"; //^-
      for (size_t ib=0;ib<batches.size();ib++) {
        // ib is batch serial number
        int batchnumber = batch(ib).num(); // batch number
        int filenum = batch(ib).FileNumber()+1; // file number (if relevant)
        // reconstruct original batch number if offset (else unchanged)
        int originalbatchnumber = batchnumber - batch(ib).BatchNumberOffset();
        //^^
        //      std::cout << "ib, runnum, bnum, origbnum, filenum "
        //                <<ib<<" "<<runnum<<" "<< batchnumber
        //                <<" "<< originalbatchnumber<<" "<< filenum<<"\n";
        if (run_flags.InSelection
            (runnum, batchnumber, originalbatchnumber, filenum)) {
          //      std::cout << "Accepted\n"; //^^
          // this batch is selected to be in this run
          if (firstbatchinrun) { // 1st batch, start a run
            ib0 = ib;
            this_run = Run(batch(ib).datasetindex(), batch(ib).DatasetID(),
                           batch(ib).PXDname());
          }
          firstbatchinrun = false;
          // Add batch to run even if not accepted
          this_run.AddBatch(batch(ib).num(), batch(ib).Accepted());
          if (batch(ib).Accepted()) {
            nbatAccepted++; // count accepted batches
            // Store run index in batch: this = current size of runlist
            batches[ib].SetRunIndex(runlist.size());
          } else {
            // Store null run index in batch
            batches[ib].SetRunIndex(-1);
          }
        }
      } // end loop batches in range
      //      std::cout <<"End Run " << runnum <<" Naccepted " << nbatAccepted <<"\n";
      if (nbatAccepted > 0) {
        // Store run index in dataset
        datasets[batch(ib0).datasetindex()].AddRunIndex(batch(ib0).PXDname(), runlist.size());
        this_run.BatchNumberOffset() = batch(ib0).BatchNumberOffset();
        this_run.FileNumber() = batch(ib0).FileNumber();
        this_run.SortList();
        this_run.RunNumber() = runnum;
        latnum = batches[batch_lookup.lookup(this_run.BatchList()[0])].LatticeNumber();
        this_run.SetLatticeNumber(latnum);
        runlist.push_back(this_run);
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
    // Remove any empty datasets
    for (size_t id=0;id<datasets.size();id++) {
      if (datasets[id].RunIndexList().size() == 0) {
        datasets.erase(datasets.begin()+id);
      }
    }
    ndatasets = datasets.size();
    // set run_controls status
    run_flags.SetStatus(+1);
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
      // Accumulate maximum resolution range for all runs, within file range
      ResoRange maxrange;
      for (size_t i=0;i<resobyrun.size();++i) { // loop runs with specified limits
        if (i == 0) {
          maxrange = resobyrun[i].second.MinRange(ResolutionRange);  // 1st range
        } else {
          maxrange = maxrange.MaxRange(resobyrun[i].second.MinRange(ResolutionRange));
        }
        int irun = FindRunIndex(resobyrun[i].first, runlist); // run index
        if (irun < 0) { // specified run not found
          clipper::String s = "\n**** Run "+clipper::String(resobyrun[i].first)+
            " specified on RESOLUTION RUN command does not exist ****\n";
          ReportErrors::printFatalError(s);
        }
        // Store resolution range limit for this run, forced to be within file range
        runlist[irun].StoreResoRange(resobyrun[i].second.MinRange(ResolutionRange));
      }  // end loop specified limits
      // If resolution ranges are defined for all runs, then reset overall limit to maximum range
      if (int(resobyrun.size()) == num_runs()) {
        ResoLimRange = maxrange;
      }
    }
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::set_run()
  //                     ^^^^^^^
  // Set run numbers into obs_part list using run index from batch
  //
  {
    if (status != RAWLIST && status != SORTED)
      ReportErrors::printFatalError("hkl_unmerge_list::set_run - not RAWLIST or SORTED");

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
              // Note that this shouldn't happen after automatic run generation, as
              // change of orientation should change run
              double rot = (clipper::Rotation(batch.Umat()) * U0inv).abs_angle();
              if (rot > rotlim) {
                ReportErrors::printWarning
                  ("WARNING: problem in run "+clipper::String(int(irun+1))+
                   " batch "+clipper::String(batchlist[ib])+
                   "  has a different orientation that of initial batch "+
                   clipper::String(batchlist[0])+
                   "\n Orientation difference = "+
                   clipper::String(clipper::Util::rad2d(rot))+"\n",
                   "WarningMessage",false);
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
      ReportErrors::printFatalError("hkl_unmerge_list::organise - not SORTED");

    if (run_set < 0)
      ReportErrors::printFatalError("hkl_unmerge_list::organise - run not set");
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
      {ReportErrors::printFatalError("hkl_unmerge_list::partials - not ORGANISED");}
    if (partial_flags.noPartials()) {
      return nopartials();  // simpler processing if no partials
    }
    if (!partial_set)
      ReportErrors::printFatalError
        ("hkl_unmerge_list::partials - no partial selection information");
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
    int latnum;
    std::vector<LatticeIndexInfo> lathkl;
    maxhkloverlappart = 0;  // maximum number of overlaps for any one part
    // count lattices
    nlattices = 0;
    nlatticesall = 0;
    maxhkloverlap = 0; // maximum number of overlapped hkl on any one observations

    // MAXNLATTICES is maximum number of lattices allowed
    // count of "main" lattice entries
    std::vector<int> numberinlattice(MAXNLATTICES+1,0); // +1 as lattices are numbered from 1
    // count of overlapped lattice entries
    std::vector<int> numberinlatticeall(MAXNLATTICES+1,0); // +1 as lattices are numbered from 1

    std::vector<Range> invresrangebydataset(ndatasets);
    std::vector<Range> invresrangebyrun(runlist.size());

    // Update lattice number ranges
    latticenumberrange.clear();
    mainlatticenumberrange.clear();

    for (size_t j = 0; j < refl_list.size(); j++) {  // loop all reflections
      obs_list.clear();   // clear temporary list
      int i = refl_list[j].first_index();
      bool obsOK = false; // true if at least one accepted observations
      int lastbatch = -1;

      while (i <= refl_list[j].last_index())  {
        // start possible observation
        //  Set values for first part or full
        int i1 = i;
        int Nfound = 1;
        int isym1 = find_part(i).isym();
        int batch1 = find_part(i).batch();
        lastbatch = batch1;
        int run1 = find_part(i).run();
        int datasetIndex = runlist[run1].DatasetIndex();
        Rtype total_fraction = find_part(i).fraction_calc();
        bool check_ok = partial_flags.check(); // false if no check on Mpart
        bool scale_frac = false; // don't scale incomplete partial
        PartFlagSwitch partial_status = FULL;
        ObservationFlag obsflag(find_part(i).ObsFlag());
        avI = find_part(i).Ic();
        if (dataflags.is_latnum) {
          latnum = find_part(i).latnum();
          mainlatticenumberrange.update(latnum);
          lathkl = find_part(i).lathkl();
          UpdateLatticeNumberRanges(lathkl);
          numberinlattice.at(latnum)++;    // count entries for each lattice
          numberinlatticeall.at(latnum)++;    // count entries for each lattice
        }

        int kpart = 1;
        //  Npart for 1st part: = 1 for a full,
        //     > 1 if extracted from MPART column,
        //     = -1 for a partial with no MPART column
        int Npart = find_part(i).Npart();

        // Test both fulls and partials, in case fulls run over > 1 part
        if (partial_flags.check()) {
          // Check 1st part for consistency of Mpart flags
          if (find_part(i).Ipart() != 1) check_ok = false;
        }
        // gap flag, store maximum gap between batches, should = 1

        batchgap = 0;
        while (++i <= refl_list[j].last_index()) { // loop parts
          // Loop through parts 2->EndObs
          // Tests for still same observation
          // same symmetry
          if (isym1 != find_part(i).isym()) break;
          //same run
          if (run1 != find_part(i).run()) break;

          if ((Npart == 1) && (find_part(i).Npart() == 1)) {
            // both parts marked as full, so not partial
            // could be partial if one part is marked as full
            break;
          }

          bool addingoverlaps = false;
          if (dataflags.is_latnum) {
            if (latnum <= 0 ){
              latnum = find_part(i).latnum();
            } else {
              // all parts should belong to the same basic lattice
              if (latnum != find_part(i).latnum()) {break;}
            }
            mainlatticenumberrange.update(latnum);
            // set flag to add them in, conditional on passing later tests
            addingoverlaps = true;
            if (lathkl.size() == 0) {
              lathkl = find_part(i).lathkl();
              UpdateLatticeNumberRanges(lathkl);
            }
          }

          // check contiguous batches: count gaps, should == 0
          int gap = 0;
          if (lastbatch >= 0) { // not 1st part
            // check for contiguous batches
            int gap = find_part(i).batch();  // this batch
            if ((gap - lastbatch) != +1) {
              //              std::cout << "Non-contiguous "<<refl_list[j].hkl().format()
              //                        <<" "<<lastbatch<<" "<<gap<<"\n"; //^
              break;   // not contiguous
            }
            gap = gap - (batch1+kpart);
          }
          lastbatch = find_part(i).batch();
          if (std::abs(gap) > 2) break;
          batchgap += gap;

          // add in any additional lathkl components
          if (addingoverlaps) {
            CombineLathkl(lathkl, find_part(i).lathkl());
            UpdateLatticeNumberRanges(find_part(i).lathkl());
          }

          // Found another part belonging to this observation
          // Check for rejection
          obsflag.AddFlag(find_part(i).ObsFlag());

          kpart++; // kpart counts accepted part
          if (partial_flags.check()) {
            // Check for consistency of Mpart flags
            if (find_part(i).Ipart() != kpart) check_ok = false;
            if (find_part(i).Npart() != Npart) check_ok = false;
          }

          total_fraction += find_part(i).fraction_calc();
          avI += find_part(i).Ic();
        } // end loop parts
        Nfound = kpart;

        // At this point we have identified Nfound parts from part i1
        //    as an observation
        // check_ok = true if MPART flags are being checked & are consistent
        //    else  = false
        // Npart is number of parts recorded for the 1st obs_part

        if (Nfound == 1) { // potential FULL
          if (Npart == 1) { // yes it is
            partial_status = FULL;
            if (dataflags.is_Mpart) {  // only check for Mosflm output
              if (total_fraction > 0.99) { // should be > 1 for full
                check_ok = true;  // OK
              }
              else {
                check_ok = false; // check total fraction
              }
            } else {
                check_ok = true;  // OK
            }
          } else {
            check_ok = false;  // check fraction
          }
        } else { // more than one part found
          if (partial_flags.check()) {
            // Check for consistency of Mpart flags
            if (Npart == Nfound) {
              partial_status = COMPLETE_CHECKED;  // unless !check_ok
            } else {
              check_ok = false;
            }
          }
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
            } else {
              partial_flags.IncrementNrejFractionTooSmall();
              //^
              //              std::cout << "hkl_unmerge_list::partials, rejected small "
              //                        << refl_symm.get_from_asu(refl_list[j].hkl(), isym1).format()
              //                        <<" fract "<< total_fraction <<"\n";
              //^-
            }
          } else if (Nfound > 1 && total_fraction > partial_flags.accept_fract_max()) {
            partial_flags.IncrementNrejFractionTooLarge();
            //^
            //      std::cout << "hkl_unmerge_list::partials, rejected large "
            //                << refl_symm.get_from_asu(refl_list[j].hkl(), isym1).format()
            //                <<" fract "<< total_fraction <<"\n";
            //^-
          } else {
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
        if (check_ok) {
          // Store observation
          // If any part of obsflag is set, mark observation as REJECTED for now
          if (Npart == 1) {
            partial_status = FULL;
          }
          Nobservations += 1;
          maxhkloverlap = Max(maxhkloverlap, int(lathkl.size()));
          obs_list.push_back(observation(
                         refl_symm.get_from_asu(refl_list[j].hkl(), isym1),
                         isym1, run1, datasetIndex, Nfound,
                         &obs_part_pointer[i1],
                         total_fraction, partial_status, obsflag, latnum, lathkl));
          //^
          //      bool DEBUG = true;
          //      if (DEBUG) {
          //        std::string s = "singleton";
          //        if (lathkl.size() > 0) {
          //          s = "multiple "+StringUtil::itos(int(lathkl.size()),1);
          //        }
          //        // print all
          //        std::cout << "\nhkl_unmerge_list::partials, " << s <<" "
          //                  << obs_list.back().hkl_original().format()<<" lattice "
          //                  << latnum << " I1 = " << find_part(i1).Ic() <<"\n   ";
          //        if (lathkl.size() > 0) {
          //          for (size_t jj=0; jj<lathkl.size(); jj++) {
          //            std::cout << " lat " <<lathkl[jj].latnum
          //                      << " " <<lathkl[jj].hkl.format();
          //          }
          //          std::cout <<"\n";
          //        }
          //        std::cout << "i, j, refl_list[j].last_index() "<< i
          //                  <<" " << j<<" "<< refl_list[j].last_index() <<"\n";
          //      }
          //^-

          invresrangebydataset[datasetIndex].update(refl_list[j].invresolsq());
          invresrangebyrun[run1].update(refl_list[j].invresolsq());
          obsOK = true;
          if (partial_status == FULL) {
            runlist[run1].Nfulls()++;
          } else {
            runlist[run1].Npartials()++;
          }
          meanI.Add(avI);

          if (dataflags.is_latnum) {
            // update counts for each lattice mentioned in lathkl list
            UpdateNumberInLattice(numberinlatticeall, lathkl);
          }

        } else {
          //^
          ///     std::cout << "Not OK " << total_fraction <<"\n"; //^-
        }
      } // observation loop
      if (obs_list.size() > 0) {
        refl_list[j].add_observation_list(obs_list);
      }
    } // reflection loop

    nlattices = 0;
    excludeoverlaps = true; // default for single lattice
    if (dataflags.is_latnum) {
      // count lattices with non-zero entries
      ASSERT (numberinlattice.size() == numberinlatticeall.size());
      for (size_t j=1; j<numberinlattice.size(); j++) { // loop from 1
        if (numberinlattice[j] > 0) {
          nlattices++;
        }
        if (numberinlatticeall[j] > 0) {
          nlatticesall++;
        }
      }
      excludeoverlaps = false;
    }

    // Set flags into runs for only||few fulls||partials
    for (size_t irun=0;irun<runlist.size();++irun) {
      runlist[irun].SetFullsAndPartials();
    }

    if (SelectI::Combine()) {
      // if we want average I for combination
      if (meanI.Count() > 0) SelectI::SetAverageIntensity(meanI.Mean());
    }
    status = PREPARED;
    updateResolutionranges(invresrangebydataset, invresrangebyrun);

    ImposeResoByRunLimits();  // mark observations if outside run limits
    return Nobservations;
  } // end ::partials
  //--------------------------------------------------------------
  int hkl_unmerge_list::nopartials()
  //                    ^^^^^^^
  //  create observation list for each reflection for nopartials
  // Returns number of observations
  {
    Nobservations = 0;
    partial_flags.Clear();

    // Temporary store for observation list for each reflection
    //   this can expand beyond allocated length if necessary
    std::vector <observation> obs_list(100);

    // Clear run reflection counts
    for (size_t irun=0;irun<runlist.size();++irun) {
      runlist[irun].clearCounts();
    }
    Rtype avI; // for each observation
    MeanSD meanI;
    int latnum;
    std::vector<LatticeIndexInfo> lathkl;
    maxhkloverlappart = 0;  // maximum number of overlaps for any one part
    // count lattices
    nlattices = 0;
    nlatticesall = 0;
    maxhkloverlap = 0; // maximum number of overlapped hkl on any one observations

    // MAXNLATTICES is maximum number of lattices allowed
    // count of "main" lattice entries
    std::vector<int> numberinlattice(MAXNLATTICES+1,0); // +1 as lattices are numbered from 1
    // count of overlapped lattice entries
    std::vector<int> numberinlatticeall(MAXNLATTICES+1,0); // +1 as lattices are numbered from 1

    std::vector<Range> invresrangebydataset(ndatasets);
    std::vector<Range> invresrangebyrun(runlist.size());

    // Update lattice number ranges
    latticenumberrange.clear();
    mainlatticenumberrange.clear();
    PartFlagSwitch partial_status = FULL; // always
    Rtype total_fraction = 0.0;

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
        ObservationFlag obsflag(find_part(i).ObsFlag());
        avI = find_part(i).Ic();
        if (dataflags.is_latnum) {
          latnum = find_part(i).latnum();
          mainlatticenumberrange.update(latnum);
          lathkl = find_part(i).lathkl();
          UpdateLatticeNumberRanges(lathkl);
          numberinlattice.at(latnum)++;    // count entries for each lattice
          numberinlatticeall.at(latnum)++;    // count entries for each lattice
        }

        //  Npart for 1st part: = 1 for a full, > 1 if extracted from MPART column,
        //     = -1 for a partial with no MPART column
        ASSERT (find_part(i).Npart() == 1);

        // Store observation
        // If any part of obsflag is set, mark observation as REJECTED for now
        Nobservations += 1;
        maxhkloverlap = Max(maxhkloverlap, int(lathkl.size()));
        obs_list.push_back(observation(
                         refl_symm.get_from_asu(refl_list[j].hkl(), isym1),
                         isym1, run1, datasetIndex, Nfound,
                         &obs_part_pointer[i1],
                         total_fraction, partial_status, obsflag, latnum, lathkl));
        invresrangebydataset[datasetIndex].update(refl_list[j].invresolsq());
        invresrangebyrun[run1].update(refl_list[j].invresolsq());
        obsOK = true;
        runlist[run1].Nfulls()++;
        meanI.Add(avI);
        if (dataflags.is_latnum) {
          // update counts for each lattice mentioned in lathkl list
          UpdateNumberInLattice(numberinlatticeall, lathkl);
        }
        i++;
      } // observation loop
      if (obs_list.size() > 0) {
        refl_list[j].add_observation_list(obs_list);
      }
    } // reflection loop

    nlattices = 0;
    excludeoverlaps = true; // default for single lattice
    if (dataflags.is_latnum) {
      // count lattices with non-zero entries
      ASSERT (numberinlattice.size() == numberinlatticeall.size());
      for (size_t j=1; j<numberinlattice.size(); j++) { // loop from 1
        if (numberinlattice[j] > 0) {
          nlattices++;
        }
        if (numberinlatticeall[j] > 0) {
          nlatticesall++;
        }
      }
      excludeoverlaps = false;
    }

    // Set flags into runs for only||few fulls||partials
    for (size_t irun=0;irun<runlist.size();++irun) {
      runlist[irun].SetFullsAndPartials();
    }

    status = PREPARED;
    updateResolutionranges(invresrangebydataset, invresrangebyrun);

    ImposeResoByRunLimits();  // mark observations if outside run limits
    return Nobservations;
  } // end ::nopartials
  //--------------------------------------------------------------
  void hkl_unmerge_list::updateResolutionranges
  (const std::vector<Range>& invresrangebydataset,
   const std::vector<Range>& invresrangebyrun)
  // Update dataset and run resolution ranges
  {
    ResoRange overallrange = ResoLimRange;
    for (int id=0;id<ndatasets;++id) { // Resolution range for each dataset
      ResoRange dtsresrange(invresrangebydataset[id]);
      // check all runs for this dataset
      std::vector<int> runindexlist = datasets[id].RunIndexList();
      ResoRange maxrunresorange = runlist.at(0).GetResoRange();
      bool runlimits = false;
      for (size_t jrun=0; jrun<runindexlist.size(); jrun++) {
        int irun = runindexlist[jrun];
        if (runlist[irun].IsResoRange()) {
          // defined resolution cutoff for this run
          maxrunresorange =
            maxrunresorange.MaxRange(runlist[irun].GetResoRange());
          //^^
          //      std::cout <<"\nrun " << irun << " res " << runlist[irun].GetResoRange().format() <<"\n";
          runlimits = true;
        }
      }
      //^^
      //      std::cout << "\nhkl_unmerge_list::updateResolutionranges "<<id <<" "
      //                <<runlimits<<" maxrange "<<maxrunresorange.format()<<
      //        " dtsresrange " <<dtsresrange.format() <<"\n";
      //^-
      if (runlimits) {
        dtsresrange = maxrunresorange; // reset dataset range
      }
      datasets[id].SetResRange(dtsresrange);
      // Overall
      overallrange = overallrange.MaxRange(datasets[id].ResRange());
    }
    for (size_t irun=0;irun<runlist.size();++irun) {
      if (!runlist[irun].IsResoRange()) {
        runlist[irun].StoreResoRange(invresrangebyrun[irun]);
      }
    }
    overallrange.ExtendRange();  // add a little tolerance
    ResoLimRange = overallrange;
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::UpdateLatticeNumberRanges(const std::vector<LatticeIndexInfo>& lathkl)
  // update maxhkloverlappart and latticenumberrange
  {
    if (lathkl.size() > 0) {
      maxhkloverlappart = Max(maxhkloverlappart, lathkl.size());
      for (size_t j=0;j<lathkl.size();++j) { // loop new lathkl
        latticenumberrange.update(lathkl[j].latnum);
      }
    }
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::UpdateNumberInLattice
  (std::vector<int>& numberinlatticeall,
   const std::vector<LatticeIndexInfo> lathkl) const
  // update counts for each lattice mentioned in lathkl list
  {
    for (size_t i=0; i<lathkl.size(); i++) {
      if (lathkl[i].latnum > 0) {
        numberinlatticeall.at(lathkl[i].latnum)++;
      }
    }
  }
  //--------------------------------------------------------------
  // add in any additional lathkl components from newlathkl into lathkl
  void hkl_unmerge_list::CombineLathkl
  (std::vector<LatticeIndexInfo>& lathkl,
   const std::vector<LatticeIndexInfo>& newlathkl) const
  {
    for (size_t j=0;j<newlathkl.size();++j) { // loop new lathkl
      // Do we have this one already?
      bool found = false;
      for (size_t i=0;i<lathkl.size();++i) { // loop lathkl
        if (newlathkl[j] == lathkl[i]) {
          found = true;
        }}
      if (!found) { // a new one, so add it
        //^
        //      std::cout << "CombineLathkl " <<j<<" "
        //                << newlathkl[j].latnum <<" "
        //                << newlathkl[j].hkl.format() <<"\n"; //^
        //      for (size_t jj=0;jj<newlathkl.size();++jj) { // loop new lathkl
        //        std::cout << "NewLathkl " << newlathkl[jj].latnum
        //                  <<" "<< newlathkl[jj].hkl.format() <<"\n"; //^
        //      }
        //      for (size_t i=0;i<lathkl.size();++i) { // loop lathkl
        //        std::cout << "Lathkl " << lathkl[i].latnum
        //                  <<" "<< lathkl[i].hkl.format() <<"\n"; //^
        //      }
        //^
        lathkl.push_back(newlathkl[j]);
      }
    }
  }
  //--------------------------------------------------------------
  int hkl_unmerge_list::sum_partials(const bool& forcesum)
  //                   ^^^^^^^^^^^
  // Sum partials within each observation for all reflections
  // If forcesum is true, sum them even if already summed
  // Returns number of valid partials
  {
    if (!forcesum && (status == SUMMED)) return Nobs_partial;
    if (!(status == PREPARED) && (status != SUMMED))
      ReportErrors::printFatalError("hkl_unmerge_list::sum_partials - not PREPARED");

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
  // Also generate resolution ranges for each dataset
  {
    if (!((status == SUMMED) || (status == PREPARED))) {
      ReportErrors::printFatalError
        ("hkl_unmerge_list::ImposeResoByRunLimits - not PREPARED or SUMMED");
    }
    if (!run_flags.IsResoByRun()) {return;}

    observation this_obs;
    ObservationStatus obs_status;

    std::vector<Range> invresrangebydataset(ndatasets);

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
          // inv resolution range by dataset
          invresrangebydataset[this_obs.datasetIndex()].update(refl_list[j].invresolsq());
        }
        this_obs.UpdateStatus(obs_status);
        refl_list[j].replace_observation(this_obs);
      } // end loop observations
    } // end loop reflections

    // Store resolution range for each dataset
    for (int id=0;id<ndatasets;++id) {
      datasets[id].SetResRange(ResoRange(invresrangebydataset[id]));
    }
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::SetResoLimits(const float& LowReso,
                                       const float& HighReso)
  {
    // resolution range limits & bins
    ResoLimRange.SetRange(LowReso, HighReso, Nobservations);
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::SetResoLimits(const ResoRange& resRange)
  {
    // resolution range limits & bins
    ResoLimRange = resRange;
    ResoLimRange.SetRange(Nobservations);
  }
  //--------------------------------------------------------------
  void hkl_unmerge_list::ResetResoLimits()
  {
    // resolution range limits reset to file range
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
      if (ndatasets <= 0) {
        ReportErrors::printFatalError("hkl_unmerge_list::cell - no datasets");
      }

      if (PXDsetName.is_blank()) {
        return averagecell;
      }
      // Find dataset
      for (int k=0; k<ndatasets; k++) {
        if (PXDsetName == datasets[k].pxdname()) {
            return datasets[k].cell();
        }}
      ReportErrors::printFatalError
        ("hkl_unmerge_list::cell - dataset not found "+PXDsetName.format());
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
    if (ResoLimRange.tbin(s2) < 0) {
      return false;
    }
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
      ReportErrors::printFatalError("hkl_unmerge_list::next_reflection - not SUMMED");

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
      ReportErrors::printFatalError(
                  "hkl_unmerge_list::get_reflection - not PREPARED or SUMMED");
    NextRefNum = jref;      // record current reflection
    return refl_list[jref];
  }
  //--------------------------------------------------------------
  // Return reflection for given (reduced) hkl
  // returns index number or -1 if missing
  int hkl_unmerge_list::get_reflection(reflection& refl, const Hkl& hkl) const
  {
    if (status != SUMMED)
      ReportErrors::printFatalError("hkl_unmerge_list::reflection - not SUMMED");
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
    // Comparison function for sorting on H,K,L,M/ISYM,[LATTNUM],BATCH
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

    if (part1->latnum() < part2->latnum())
      return true;
    else if (part1->latnum() > part2->latnum())
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
                                        const bool& reindexSecondaryLattices,
                                        const bool& AllowFractIndex)
  // Change symmetry in all internal lists to new_spgp,
  // ie
  // 1. for each obs_part, get original indices
  // 2. reindex (change basis)
  // 3. rereduce to new asymmetric unit
  // 4. flag as RAWLIST (unsorted)
  //
  // reindexSecondaryLattices true if secondary lattices from multilattice overlaps
  //  should be reindexed as well
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
    //^       << new_symm.symbol_xHM() << "  reindex " << reindex_op.a  s_hkl()
    //^       << "\n";

    if (dataflags.is_latnum) {
      ASSERT (nlattices > 0);    // should have already called partials() to set nlattices
    }

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
        // Multilattice
        if (dataflags.is_latnum && reindexSecondaryLattices) {
          observation_part& part = find_part(i);
          std::vector<LatticeIndexInfo> lathkl = part.lathkl();
          for (int j=0;j<nlattices;++j) {
            HklOK = HklOK || lathkl[j].hkl.change_basis(hkl_reindex, reindex_op);
          }
          if (!HklOK) {
            NfractIdx++;
            obs_part_pointer[i] = NULL;  // Clear pointer
            continue;
          }
          find_part(i).set_lathkl(lathkl);
        }
      } else {  // no reindex
        hkl_new = new_symm.put_in_asu(hkl, new_isym);
      }
      //^
      //        std::cout << find_part(i).hkl().format() << " " << find_part(i).isym()
      //                  << " " << hkl.format()
      //                  << " " << hkl_new.format() << " " << new_isym
      //                  << " " << find_part(i).batch() << "\n";
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
        ReportErrors::printFatalError
          ("hkl_unmerge_list::change_symmetry: illegal fractional indices generated by reindex operator");
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
      ReportErrors::printFatalError("hkl_unmerge_list::prepare - EMPTY");
    if (N_part_list == 0) {
      ReportErrors::printFatalError("hkl_unmerge_list::prepare  No observations in list");
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
    int pole0 = pole;
    if (pole < 0) {
      // pole unspecified, use value from 1st batch for consistency
      batches[0].SetPole(pole);
      pole0 = batches[0].Pole();
    }
    for (int ib=0; ib<nbatches; ib++) {
      batches[ib].SetPole(pole0);
    }
  }
  //--------------------------------------------------------------
  bool hkl_unmerge_list::validOrientation() const
  // Returns true if OK, false if not OK eg some batch does not have valid Umat
  {
    // Check that all batches have valid orientation
    bool OK = true;
    for (int i = 0; i < nbatches; i++)  {
      if (!batches[i].ValidOrientation()) {
        OK = false;
        break;
      }
    }
    return OK;
  }
  //--------------------------------------------------------------
  bool hkl_unmerge_list::CalcSecondaryBeams(const int& pole)
  // Calculate all secondary beam directions, in chosen frame, also
  // diffraction vectors d*vec
  // On entry:
  //  pole  =  0 SECONDARY  camera frame
  //       !=  0 ABSORPTION, crystal frame = 1,2,3 for h,k,l,
  //        = -1 unspecified, use closest reciprocal axis from 1st batch, use for all
  // not yet!        = -1 unspecified, use closest reciprocal axis for each run
  //
  // Returns true if OK, false if not OK eg some batch does not have valid Umat
  {
    SetPoles(pole);
    if(!validOrientation()) {return false;}

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
        this_obs.StoreS(fsPhi);   // rlu
        refl_list[j].replace_observation(this_obs);
      }
    }
    return true;
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
  //  sPhi      diffraction vector at actual phi position (rlu)
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
  //  sPhi      diffraction vector at actual phi position (camera frame) (rlu)
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
    // s(r) = [R][D][U][B]h  camera frame  (rlu)
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
  Range hkl_unmerge_list::UpdatePolarizationCorrections(const bool& Total,
                                        const PolarizationControl& polarizationcontrol)
  // Update polarization corrections for all parts
  // If Total == true, then apply complete correction
  //   else assume the unpolarised correction is already applied, apply
  //   additional correction for polarised incident beam
  // PolarizationControl contains:
  //  polarizationfactor if fraction polarised, = 0 for unpolarised, ~ +0.9 for synchrotrons
  //  direction of polarization in current coordinate frame
  // Returns range of correction factors
  //
  // see Kahn, Fourme, Gadet, Janin, Dumas, & Andre,
  // J. Appl. Cryst. (1982). 15, 330-337
  {
    double polarizationfactor = polarizationcontrol.Factor();
    // component of electric vector in synchrotron plane, = Pn x s0, unit vector
    clipper::Vec3<double> EprimePi = polarizationcontrol.EprimePi();


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
      // zp = projection of diffraction vector on to E'pi
      double zp = clipper::Vec3<double>::dot(sPhi, EprimePi);
      //^^
      //      double zd = clipper::Vec3<double>::dot(sPhi, clipper::Vec3<double>(0.0,0.0,1.0));
      //      std::cout << "z " << sPhi[2] << " zp " << zp <<" zd "<<zd<< "\n"; //^-
      //      //zp = sPhi[2]; // testing
      //^-
      double sinSqtheta = 0.25 * (sPhi * sPhi);    // |s| = 2 sin theta; sin^2 theta = 0.25*|s|^2
      double cos2theta = 1. - 2.0 * sinSqtheta;    // cos 2theta = cos^2 theta - sin^2 theta
                                                   //  = 1 - 2 sin^2 theta
      double cosSq2theta = cos2theta * cos2theta;  // cos^2 2theta
      double sinSq2theta = 1.0 - cosSq2theta;      // sin^2 2theta
      double cosrho = zp/sqrt(sinSq2theta);         // cos rho = zp/sin 2theta
      // P0 = 1/2 [ 1 + cos^2 2 theta]
      double P0 = 0.5*(1.0 + cosSq2theta);  // unpolarised part
      // P' = 1/2 (-Xsi') cos 2rho sin^2 2theta
      double PP = 0.5 * polarizationfactor * (2.0*cosrho*cosrho - 1.0) * sinSq2theta;
      double PolFac = 1.0/P0;
      if (!Total) {
        // P0 already applied, so just correct it
        //  complete PolFac = P0 - P', negative sign because polarizationfactor should be negative
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
      //                <<" "<< PolFac <<" " << clipper::Util::rad2d(acos(cosrho))<<"\n";
      //^-
    } // end loop parts
    //^
    //    std::cout << "Range of polarization corrction factors "
    //        << PFrange.min() <<" " << PFrange.max() <<"\n";
    return PFrange;
  }
  //--------------------------------------------------------------
  //! Store Space group status (CMtz::SYMGRP.spg_confidence) into mtzsym
  void hkl_unmerge_list::SetSpaceGroupStatus(const char& spg_status)
  {
    mtzsym.spg_confidence = spg_status;
  }
  //--------------------------------------------------------------
  //! Store use run flags
  void hkl_unmerge_list::StoreUseRun(const std::vector<bool>& userun)
  {
    ASSERT (userun.size() == runlist.size());
    run_flags.StoreUseRun(userun);
    for (size_t i=0; i<runlist.size(); i++) {
      runlist[i].StoreUse(userun[i]);  // store use flag in runs
    }
  }
  //--------------------------------------------------------------
  int hkl_unmerge_list::num_accepted_datasets() const
  //!< number of accepted datasets
  {
    int ndts = 0;
    for (size_t id=0; id<datasets.size(); id++) {
      if (datasets[id].accepted()) {
        ndts++;
      }
    }
    return ndts;
  }
  //--------------------------------------------------------------
  std::vector<Dataset> hkl_unmerge_list::AllAcceptedDatasets() const
  //!< all accepted datasets
  {
    std::vector<Dataset> valid_datasets;
    for (size_t id=0; id<datasets.size(); id++) {
      if (datasets[id].accepted()) {
        valid_datasets.push_back(datasets[id]);
      }
    }
    return valid_datasets;
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
      ReportErrors::printFatalError("hkl_unmerge_list: illegal copy constructor");
    }
    clear();
  }
//--------------------------------------------------------------
  // Copy operator throws exception unless object is EMPTY
  hkl_unmerge_list& hkl_unmerge_list::operator= (const hkl_unmerge_list& List)
  {
    if (List.status != EMPTY) {
      ReportErrors::printFatalError("hkl_unmerge_list: illegal copy operation");
    }
    clear();
    return *this;
  }

} // namespace
