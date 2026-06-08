//
// reject.cpp
//
// set or reset outlier flags for all reflections

#include "reject.hh"
#include "selectedobservations.hh"
#include "file_util.hh"
#include "observationstatuscontrol.hh"
#include "string_util.hh"
#include "weighttype.hh"

#include <assert.h>
#define ASSERT assert

namespace scala {
  // ------------------------------------------------------------
  void ClearOutlierFlags(hkl_unmerge_list& hkl_list)
  //  Clear outlier status flags for all observations
  //  Other flags are left unaltered
  {
    reflection this_refl;
    observation this_obs;
    hkl_list.rewind();

    // loop all reflections unconditionally
    for (int jref=0;jref<hkl_list.num_reflections();++jref) {
      this_refl = hkl_list.get_reflection(jref);
      // loop all observations, ignoring accept flag
      for (int lobs=0;lobs<this_refl.num_observations();++lobs) {
        this_obs = this_refl.get_observation(lobs);
        ObservationStatus status = this_obs.ObsStatus();
        status.UnsetOutlier();
        this_obs.UpdateStatus(status);
        this_refl.replace_observation(this_obs);
      }
      hkl_list.replace_reflection(this_refl);
    }
  }
  // ------------------------------------------------------------
  class RejectList
  {
  public:
    // Construct & store addresses
    RejectList(reflection& This_Refl,
               const OutlierControl& Outliercontrol)
      : this_refl(&This_Refl), outliercontrol(&Outliercontrol), discrepant(false)
    {
      // clear deviation list
      deviations.assign(This_Refl.num_observations(), 0.0);
      // Test for "weak" reflections in Emax test is based on E2/sigE2
      // Make the acceptance criterion stricter than the outlier test
      const float SDREJ_RATIO = 0.6;
      sdrej = SDREJ_RATIO*
	Outliercontrol.Reject(IPLUS, AnomalousClass()).sdrej;
    }

    // Check for outliers in specified class, flags with status, accumulate flags
    void Check(const AnomalousClass& selclass, const int& Dts_Index,
               const bool& anomoutlier);

    // Check all observations against Emax test
    // Accumulate indices to outlier rejects from vector "selected" into vectors
    // "rejected " and "statusflags". Status for rejects is in "status"
    void CheckEmax(const Normalise& NormRes, const EProb& eprobtest,
                   const bool& Centric);

    const std::vector<int>& Rejected() const  {return rejected;}
    const std::vector<ObservationStatus>& Statusflags() const {return statusflags;}
    std::vector<float> Deviations() const {return deviations;}
    //! Return true if outliers found even if not rejected
    bool Discrepant() const {return discrepant;}

    // Update reflection with new observation status flags
    void UpdateReflection();

    static void setWeightType(const WeightType::AverageWeightType& Weighttype)
    {weighttype = Weighttype;}

  private:
    reflection* this_refl;
    const OutlierControl* outliercontrol;
    float sdrej;  // for Emax test

    std::vector<int> rejected;
    std::vector<ObservationStatus> statusflags;
    std::vector<float> deviations;
    bool discrepant;
    // Weighting scheme for averaging observations in outlier testing
    static WeightType::AverageWeightType weighttype;   // type of weighting for average

    // Return list of index numbers for each Emax outlier observation, if any
    std::vector<int> EmaxRejectIndexList
    (SelectedObservations& selobs,
     const Normalise& NormRes, const EProb& eprobtest,
     const bool& Centric,
     const float& sdrej) const;

    // negate all indices unless already negated, -1 to distinguish +0 from -0
    std::vector<int> negateIndices(const std::vector<int>& idxv) const;
    int negateIndex(const int& idxin) const;

  };
  // Default weight type, mostly set via keywords and their default
  WeightType::AverageWeightType RejectList::weighttype = WeightType::SQRTSCALE;
  // ------------------------------------------------------------
  void RejectList::Check(const AnomalousClass& selclass, const int& dts_index,
                         const bool& anomoutlier)
  // Select observations according to dataset index dts_index and anomalous class
  // Accumulate indices to outlier rejects from vector "selected" into vectors
  // "rejected " and "statusflags".
  // anomoutlier true is comparing I+ with I-
  //
  {
    SelectedObservations sel(*this_refl, dts_index, selclass, weighttype);
    RejectFlags rejflags = outliercontrol->Reject(selclass, dts_index);
    // outlierindexlist, list of (indices+1) to rejects, or if negated, deviants
    std::vector<int> outlierindexlist =
      sel.OutlierIndexList(rejflags);
    discrepant = discrepant || sel.Discrepant();

    ObservationStatus status;
    // Append to rejected and statusflags arrays
    for (size_t k=0; k<outlierindexlist.size(); k++) {
      ASSERT (outlierindexlist[k] != 0);
      size_t kk = std::abs(outlierindexlist[k])-1;
      if (outlierindexlist[k] < 0) {
        status.SetDeviant(); // ObservationStatus::OBSSTAT_DEVIANT
      } else if (anomoutlier) {
        status.SetOutlierAnom(); // ObservationStatus::OBSSTAT_OUTLIERANOM
      } else if (outlierindexlist[k] > 0) {
        status.SetOutlier(); //ObservationStatus::OBSSTAT_OUTLIER
      }
      rejected.push_back(kk);
      statusflags.push_back(status);
    }
    std::vector<float> ddi = sel.DeltaAll();
    ASSERT (ddi.size() == deviations.size());
    for (size_t i=0;i<ddi.size();++i) {
      if (ddi[i] != 0.0) {deviations[i] = ddi[i];}  // copy deviations
    }
  }
  // ------------------------------------------------------------
  void RejectList::CheckEmax(const Normalise& NormRes, const EProb& eprobtest,
                             const bool& Centric)
  // Check all observations against Emax test
  // Accumulate indices to outlier rejects from vector "selected" into vectors
  // "rejected " and "statusflags".
  // Status for rejects are OBSSTAT_EMAX, or OBSSTAT_EMAX_OK if kept
  //
  {
    SelectedObservations sel(*this_refl, -1, ALL);
    // get index list of rejected observations, negated to keep
    std::vector<int> outlierindexlist =
      EmaxRejectIndexList(sel, NormRes, eprobtest, Centric, sdrej);
    if (outlierindexlist.size() > 0) {
      //      std::cout <<"CheckEmax number = " << outlierindexlist.size()
      //                <<" " << rejected.size() <<" "<<sel.Number()
      //                <<" "<<sel.hkl().format()
      //                <<std::endl; //^
      for (size_t k=0; k<outlierindexlist.size(); k++) {
        ObservationStatus status;
        int i = outlierindexlist[k];
        if (i >= 0) {
          status.SetEmax();
        } else {
          // keep
          i = -i-1;
          status.SetEmaxOK();
        }
        rejected.push_back(i);
        statusflags.push_back(status);
      }
      discrepant = true;
    }
  }
  // ------------------------------------------------------------
  void RejectList::UpdateReflection()
  // Update reflection with new observation status flags
  {
    if (rejected.size() > 0) {
      // Set outlier flags back into observations within reflections
      for (size_t i=0;i<rejected.size();++i) {
        observation this_obs = this_refl->get_observation(rejected[i]);
        ObservationStatus status = this_obs.ObsStatus();
        this_obs.UpdateStatus(status.mergestatus(statusflags[i]));
        this_refl->replace_observation(this_obs);
      }
    }
  }
  // ------------------------------------------------------------
  // ------------------------------------------------------------
  void RejectOutlier(hkl_unmerge_list& hkl_list,
                     const SDmodel& SDM,
                     const Normalise& NormRes,
                     const bool& anomOn,
                     const OutlierControl& outliercontrol,
                     WriteRogues& RoguesList)
  // Check for outliers in all reflections, using parameters in outliercontrol,
  // and set status flags as required on each observation
  //
  //  SDM              current sd correction model
  //  NormRes          normalisation for Emax test
  //  anomOn           true to do main outlier check only within the
  //                   I+ and I- sets
  //  outliercontrol   parameters for rejection
  //  RoguesList       optional rogues list output
  {
    reflection this_refl, temp_refl;
    observation this_obs;
    int ndsets = hkl_list.num_datasets();
    int ndsetscheck = ndsets;

    // First clear all outlier & other status flags (except ObsFlags)
    ClearObsStatus(hkl_list);

    // if OUTLIER COMBINE option, check outliers across all datasets
    if (outliercontrol.Combine()) ndsetscheck = 1;

    EProb eprobtest = outliercontrol.EMaxTest();  // Emax test

    // Weighting type for outliers
    RejectList::setWeightType(outliercontrol.weightType());

    hkl_list.rewind();

    while (hkl_list.next_reflection(this_refl) >= 0)  {  // loop reflections
      bool Centric = hkl_list.symmetry().is_centric(this_refl.hkl());
      temp_refl = this_refl;  // copy
      // start reject list, store pointers to temp_refl and outliercontrol
      RejectList rejlist(temp_refl, outliercontrol);
      //  Apply current SD correction to reflection (all observations)
      SDM.CorrectReflection(temp_refl);
      //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      //      if (this_refl.hkl() == Hkl(-32,-9,1)) {
      //	std::cout << this_refl.hkl().format() << " <<<<\n";
      //      }
      // Emax test, before outlier test
      if (!(eprobtest.Null())) {
        rejlist.CheckEmax(NormRes, eprobtest, Centric);
        rejlist.UpdateReflection();  // update flags
      }

      //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // Outlier check
      // Loop datasets
      for (int idts=0;idts<ndsetscheck;++idts) {
        int dts_index = idts;
        // if OUTLIER COMBINE option, check outliers across all datasets
        if (outliercontrol.Combine()) dts_index = -1;

        if (anomOn && !Centric) {
          // Anomalous, check I+ & I- seperately
          rejlist.Check(IPLUS, dts_index, false);
          rejlist.Check(IMINUS,dts_index, false);
          if (outliercontrol.Anom()) {
            // check for outliers between I+ & I-
            rejlist.UpdateReflection();  // update flags for rejections within I+/-
            rejlist.Check(BOTH, dts_index, true);
          }
        } else {
          // no anomalous, just check all observations
          rejlist.Check(ALL, dts_index, false);
        }
        rejlist.UpdateReflection();  // update flags for rejections within I+/-
      }  // end loop datasets
      //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      //  put rejected list in "rejected" & "status"
      //   list of rejected or deviant (negated) observations ...
      std::vector<int> rejected(rejlist.Rejected().begin(), rejlist.Rejected().end());
      //   ... and their status
      std::vector<ObservationStatus> statusflags(rejlist.Statusflags().begin(),
                                                 rejlist.Statusflags().end());

      ASSERT (rejected.size() == statusflags.size());
      if (rejlist.Discrepant() > 0) {
        if (rejected.size() > 0) {
          // Set outlier flags back into reflection
          for (size_t i=0;i<rejected.size();++i) {
            this_obs = this_refl.get_observation(rejected[i]);
            this_obs.UpdateStatus(statusflags[i]);
            this_refl.replace_observation(this_obs);
          }
          hkl_list.replace_reflection(this_refl);
        }
        // Optional output to ROGUES file & ROGUEPLOT
        if (RoguesList.Open()) {
          SDM.CorrectAllReflection(this_refl);
          RoguesList.RogueReflection(this_refl, rejlist.Deviations(), NormRes);
        }
      }
    } // end loop reflections
  }
  // ------------------------------------------------------------
  std::vector<int> CountOutliers(const hkl_unmerge_list& hkl_list)
  {
    std::vector<int> rejectedbatch;
    std::vector<int> rejecteddataset;
    return CountOutliers(hkl_list, rejectedbatch, rejecteddataset);
  }
  // ------------------------------------------------------------
  std::vector<int> CountOutliers(const hkl_unmerge_list& hkl_list,
                                 std::vector<int>& rejectedbatch,
                                 std::vector<int>& rejecteddataset)
  // Returns counts of flagged outliers within I+/- & between +/-
  // return[0] number of rejects [1] number on I+- [2] number on Emax
  //
  // On exit:
  //  rejectedbatch   count of rejected reflections for each batch
  //  rejecteddataset count of rejected reflections for each dataset
  {
    reflection this_refl;
    int n = 0;
    int na = 0;
    int nemax = 0;
    int nbatches = hkl_list.num_batches();
    observation this_obs;
    int ndatasets = hkl_list.num_datasets();
    rejectedbatch.assign(nbatches, 0);
    rejecteddataset.assign(ndatasets, 0);

    hkl_list.rewind();

    // loop all reflections unconditionally
    for (int jref=0;jref<hkl_list.num_reflections();++jref) {
      this_refl = hkl_list.get_reflection(jref);
      bool rejref = false;
      // loop all observations, ignoring accept flag
      for (int lobs=0;lobs<this_refl.num_observations();++lobs) {
        this_obs = this_refl.get_observation(lobs);
        int batchnum = this_obs.Batch();
        int jbatch = hkl_list.batch_serial(this_obs.Batch()); // batch serial
        int jdataset = this_obs.datasetIndex();
        ObservationStatus status = this_obs.ObsStatus();
        bool rej = false;
        if (status.TestOutlier()) {
          n++;
          rej = true;
        }
        if (status.TestOutlierAnom()) {
          na++;
          rej = true;
        }
        if (status.TestEmax()) {
          nemax++;
          rej = true;
        }
        if (rej) {
          rejectedbatch[jbatch]++;
          rejecteddataset[jdataset]++;
          rejref = true;
        }
      }
    }
    std::vector<int> rejv(3);
    rejv[0]=n; rejv[1]=na; rejv[2]=nemax;
    return rejv;
  }
  // ------------------------------------------------------------
  std::string CountOutliersXML(const std::vector<int>& nrejs)
  {
    std::string s = "<Outliers>\n";
    // within I+, I- sets
    s += StringUtil::MakeXMLtag("RejectNumberUnique", nrejs[0])+"\n";;
    // between I+ & I-
    s += StringUtil::MakeXMLtag("RejectNumberFriedel", nrejs[1])+"\n";
    // on Emax
    s += StringUtil::MakeXMLtag("RejectNumberEmax", nrejs[2])+"\n";
    s += "</Outliers>\n";
    return s;
  }
  // ------------------------------------------------------------
  // negate indices unless already negated, -1 to distinguish +0 from -0
  //  > 0 reject, < 0 keep
  std::vector<int> RejectList::negateIndices(const std::vector<int>& idxlist) const
  {
    std::vector<int> idxv = idxlist;
    for (size_t k=0; k<idxv.size(); k++) {
      idxv[k] = negateIndex(idxv[k]);
    }
    return idxv;
  }
  // ------------------------------------------------------------
  // negate index unless already negated, -1 to distinguish +0 from -0
  //  > 0 reject, < 0 keep
  int RejectList::negateIndex(const int& idxin) const
  {
    int idx = idxin;
    if (idx >= 0) {
        idx = -idx-1;
    }
    return idx;
  }
  // ------------------------------------------------------------
  std::vector<int> RejectList::EmaxRejectIndexList
  (SelectedObservations& selobs,
   const Normalise& NormRes, const EProb& eprobtest,
   const bool& Centric,
   const float& sdrej) const
  // Return list of index numbers for each Emax outlier observation, if any
  {
    std::vector<int> idxlist;
    reflection this_ref = selobs.Reflection();
    Rtype invresolsq = this_ref.invresolsq();
    DVect3 rhkl = this_ref.hkl().real();

    // This is a real fudge!
    float bigtestfactor = NormRes.mmratio(); // a measure of mean/median

    observation this_obs;
    int i;
    int n = 0;
    bool ilist = true;
    float corr = NormRes.Corr(invresolsq);
    Hkl hkl;
    bool allok = true;
    int nfail = 0;

    while ((i = selobs.next_observation(this_obs)) >= 0) {
      this_obs =this_ref.get_observation(i);
      float E2 = NormRes.apply(this_obs.kI(), invresolsq, rhkl);
      n++;
      if (eprobtest.TooBig(E2, Centric)) {
	
	//	if (nfail == 0) {
	//	  std::cout <<"\n"<< this_ref.hkl().format() << " N: " <<selobs.Number()
	//		    << " sdrej "<< sdrej <<  "\n";
	//	}

	allok = false;  // at least one marked as too big (for now)
	nfail++;  // count failures
	// is it weak?
	float sigE2 = NormRes.apply(this_obs.ksigI(), invresolsq, rhkl);
	// E^2 too big, positive or negative
	//  check E^2/sig(E^2), and keep weak observations, as the
	//  normalisation is unreliable in very weak shells
	//    (maybe pending improved normalisation)
	if (std::abs(E2/sigE2) > sdrej) {
	  idxlist.push_back(i);
	  //	  std::cout << "**Reject "<<" "<<E2<<" "<<sigE2<<
	  //	    " "<< E2/sigE2<<" "<<" "<<std::endl;
	} else { 
	  //	  std::cout << "**Keep "<<" "<<E2<<" "<<sigE2<<
	  //	  " "<< E2/sigE2<<" "<<" "<<std::endl;
	  idxlist.push_back(-i-1);  // flag for accepted Emax fail
	}
      }
    }
    return idxlist;
  }
  /*
    if (allok) {
      // all observations pass the Emax test, quit returning empty list
      return idxlist;
    }

    std::cout <<"\n"<< this_ref.hkl().format() << " N: " <<selobs.Number()
	      << " sdrej "<< sdrej <<  "\n";

    // At least one observation failed initial Emax test

    // Count the number of weak observations
    bool unflaggedweak = true;  // if all the unflagged observations are weak
    int nweak = 0;
    while ((i = selobs.next_observation(this_obs)) >= 0) {
      this_obs =this_ref.get_observation(i);
      float E2 = NormRes.apply(this_obs.kI(), invresolsq, rhkl);
      float sigE2 = NormRes.apply(this_obs.ksigI(), invresolsq, rhkl);
      // E^2 too big, positive or negative
      //  check E^2/sig(E^2), and keep weak observations, as the
      //  normalisation is unreliable in very weak shells
      //    (maybe pending improved normalisation)
      if (std::abs(E2/sigE2) <= sdrej) {
	nweak++;
	std::cout << "**Keep "<<" "<<E2<<" "<<sigE2<<
	  " "<< E2/sigE2<<" "<<" "<<std::endl;
      } else {
	// not weak, is it flagged?
	if (flags[i]) {
	  // not flagged but not weak
	  unflaggedweak = false;  // not all unflagged obs are weak
	} else {
	  // flagged and not weak
	  std::cout << "**Reject "<<" "<<E2<<" "<<sigE2<<
	  " "<< E2/sigE2<<" "<<" "<<std::endl;
	}
      }
    }  // end obs loop
    // If all observations are weak, keep all
    if (nweak == nobs) {
      idxlist = negateIndices(idxlist);
    }
    return idxlist;
  }
    */

  /*
  //  older program versions had more exceptions to the reject/keep options
  //   0.7.12 skip this
    // If we have at least one reject, check that we really want to reject them
    // if all or most observations are rejected, then we probably want to keep them
    // if (false) {  //
    if (idxlist.size() > 0) {  // TESTING
      // fraction rejected
      double fracrejected = double(idxlist.size())/double(n);
      const double REJFRAC = 0.8;
      if (fracrejected > REJFRAC && n > 2) {
        // More than REJFRAC rejected, and more than 2 observations
        // Keep all, negate indices (-1 to distinguish 0), if more than 2
        idxlist = negateIndices(idxlist);
	std::cout << "*^* Emax test, keeping all observations for "
		  << selobs.hkl().format()<<std::endl;
      } else {
        // Not all rejected, but are all the observations close to the limit?
        // First test the average I
        IsigI E2av = NormRes.apply(selobs.Average(), invresolsq, rhkl);
        const float AVTESTFACTOR = 0.8;
        if ((eprobtest.TooBig(E2av.I(), Centric, AVTESTFACTOR)) &&
            (std::abs(E2av.I()/E2av.sigI()) > sdrej)) {
          // ... but not really too big
          if (!eprobtest.TooBig(E2av.I(), Centric, bigtestfactor)) {
            // Keep all, negate list indices
            idxlist = negateIndices(idxlist);
	    std::cout << "*^* Emax test on average, E2av = "<<E2av.I()<<" "
		      <<E2av.sigI() <<", Av(I) " <<selobs.Average().I()<<" "
		      <<selobs.Average().sigI()
		      <<", keeping all observations for "
		      << selobs.hkl().format()<<std::endl;
	  } else {
	    std::cout << "*^* Emax test on average, E2av = "<<E2av.I()<<" "
		      <<E2av.sigI() <<", Av(I) " <<selobs.Average().I()<<" "
		      <<selobs.Average().sigI()
		      <<", REJECTING all observations for "
		      << selobs.hkl().format()<<std::endl;
          }
        } else {
          // Test individual observations against smaller test value
          //  count how many are above that
          const float TESTFACTOR = 0.64;  // 0.8^2
          int nabovetest = 0;
          n = 0;
          while ((i = selobs.next_observation(this_obs)) >= 0) {
            this_obs =this_ref.get_observation(i);
            IsigI E2sigE2 =
	      NormRes.apply(this_obs.kI_sigI(), invresolsq, rhkl);
            n++;
	    std::cout << "*** Emax test, E2 "<<E2sigE2.I()
		      <<" "<<E2sigE2.sigI()<<" "<<sdrej<<std::endl; //^^
            if (eprobtest.TooBig(E2sigE2.I(), Centric, TESTFACTOR)) {
              // E^2 too big, positive or negative
              //  check E^2/sig(E^2), and keep weak observations, as the
              //  normalisation is unreliable in very weak shells (pending
              //  improved normalisation
              if (std::abs(E2sigE2.I()/E2sigE2.sigI()) > sdrej) {
                nabovetest++;
              }
            }
          }
          fracrejected = double(nabovetest)/double(n);
          if (fracrejected > REJFRAC) {
            // most above smaller limit
            // Keep all, negate indices for previously rejected observations
            idxlist = negateIndices(idxlist);
	    std::cout << "*^* Emax test on individuals, nabove="<<
	      nabovetest<<" of "<<n<<
	      ", keeping all observations for "
		      << selobs.hkl().format()<<std::endl;
	  } else {
	    std::cout << "*!* Emax test on individuals, nabove="<<
	      nabovetest<<" of "<<n<<
	      ", rejecting some observations for "
		      << selobs.hkl().format()<<std::endl;
          }
        }
      }
    }
    return idxlist;
    }
    */
  // ------------------------------------------------------------
  // ------------------------------------------------------------
  WriteRogues::WriteRogues(const std::string& filename,
                           const bool& Start, const bool& Plot,
                           const bool& multilattice,
                           const std::string& title, const float& dstarMax,
                           const float& wavelength,
			   const Rings& icerings,
			   const OutlierControl& outliercontrol,
                           const bool& xmgraceoutput)
  // Open ROGUES file & write header if Start true
  // Open ROGUESPLOT file & write header if Plot true
  // If filename = "", no xmgr output
  // multilattice = true is there are multiple lattices
  // title & maximum resolution d* = lambda/d
  //  outliercontrol   parameters for rejection
  // xmgraceoutput if true, write xmgr file ROGUEPLOT
  {
    if (Start) {
      rogues = OpenFile(filename, true);  // open ROGUES file
      std::string s = "The ROGUES file contains all monitored outliers";
      if (outliercontrol.Reject(ALL).formatReject2Policy() ==
          "KEEP") {
        s += ", including unrejected outliers measured twice";
      }

      fprintf(rogues, "%s", s.c_str());
      std::string rs =
        std::string("\nRej = '*', '@' for I+- rejects, '#' for Emax rejects, ")+
        "'x' for accepted flagged observation,\n"+
        "   'd' for deviant but kept, '$' for >Emax but kept because the observation is small";
      if (multilattice) {
        rs += ",\n      'M' for multiple lattice overlaps";
      }
      rs += "\nFlag is taken from the Mosflm FLAG column or XDS MISFIT";
      rs += "\n";
      fprintf(rogues, "%s", rs.c_str());
      fprintf(rogues,
              "TotFrc = total fraction, fulls (f) or partials (p),");
      fprintf(rogues,
              " Bijv I+ or I- for Bijvoet classes\n");
      fprintf(rogues,
              "DelI/sd = (Ihl - Mn(I)others)/sqrt[sd(Ihl)**2 + sd(Mn(I))**2]\n\n");
      if (multilattice) {
        rs = std::string("Note that multilattice overlapped observations are not used in outlier calculation nor in means,\n")+
          "  and are listed here only under one of their hkl indices\n\n";
        fprintf(rogues, "%s", rs.c_str());
      }
      fprintf(rogues,
 "Flagged observations kept are labelled as: B BGratio; P PKratio; N TooNeg; G BGgradient; O Overload; E Edge; X XDS MISFIT\n");
      fprintf(rogues,
 "Deviant reflections with two measurements are always listed. Policy for deviant reflections measured twice: %s\n\n",
              outliercontrol.Reject(ALL).formatReject2Policy().c_str());


      fprintf(rogues,
              "   h   k   l     h   k   l  Batch      I  sigI    E  TotFrc ");
      fprintf(rogues,
              "Bijv  Scale DelI/sd d(A)   Xdet   Ydet    Phi   LP   Rej Flag\n");
      fprintf(rogues,
              "   (measured)     (unique)\n");
    } else {
      rogues = NULL;
    }
    if (Plot) {
      std::string roguefilename = "ROGUEPLOT";
      if (!xmgraceoutput) {roguefilename = "";}
      rogueplot = RoguePlot(roguefilename, title, dstarMax,
			    wavelength, icerings);
      rogueplot.Start();
    }
  }
  // ------------------------------------------------------------
  void WriteRogues::RogueReflection(const reflection& this_refl,
                                    const std::vector<float>& deviations,
                                    const Normalise& NormRes)
  // Write rogues entry
  {
    bool outlier;
    char partial;
    std::string reject;
    std::string PlusMinus;
    std::string flagtype = "";
    fprintf(rogues,"\n");
    float invresolsq = this_refl.invresolsq();
    float d = 1./sqrt(invresolsq);
    DVect3 rhkl = this_refl.hkl().real();

    for (int lobs=0;lobs<this_refl.num_observations();++lobs) {
      outlier = false;
      observation obs = this_refl.get_observation(lobs);
      ObservationStatus status = obs.ObsStatus();
      // List observations which are accepted or outliers, not
      // non-accepted non-outliers (ie observations rejected on
      // observationflags
      if (!status.TestObsFlag()) {
        partial = 'p';
        if (obs.IsFull()) {partial = 'f';}

        reject = ' ';
        if (!obs.Observationflag().OK()) {
          reject += 'x';
        }
        if (status.TestOutlier()) {
          reject += '*';
          outlier = true;
        } else if (status.TestOutlierAnom()) {
          reject += '@';
          outlier = true;
        } else if (status.TestEmax()) {
          reject += '#';
          outlier = true;
        } else if (status.TestEmaxOK()) {
          reject += '$';
        } else if (status.TestDeviant()) {
          reject += 'd';
        }

        if (obs.Isym()%2 == 0) {
          PlusMinus = "I-";
        } else {
          PlusMinus = "I+";
        }
        if (!obs.IsSingleton()) {
          // multiple lattice overlap
          reject += 'M';
        }
        float scale = obs.Gscale();
        if (scale != 0.0) scale = 1./scale;
        std::pair<float,float> XY = obs.XYdet();

        float E = NormRes.apply(obs.kI(), invresolsq, rhkl);
        if (E > 0.0) {E = sqrt(E);}
        else {E = -sqrt(-E);}
        std::string flagtype = obs.Observationflag().format(); // 6 characters

        fprintf(rogues,
        "%4d%4d%4d  %4d%4d%4d%6d%8d%6d%6.2f%6.1f%c   %2s%7.3f%7.1f%6.2f%7.1f%7.1f%7.1f %5.4f %s %s\n",
                obs.hkl_original()[0],obs.hkl_original()[1],obs.hkl_original()[2],
                this_refl.hkl()[0],this_refl.hkl()[1],this_refl.hkl()[2],
                obs.Batch(), Nint(obs.kI()), Nint(obs.ksigI()),
                E, obs.TotalFraction(),
                partial,  PlusMinus.c_str(), scale, deviations[lobs],
                d, XY.first, XY.second, obs.phi(), obs.LP(), reject.c_str(), flagtype.c_str());
        // ROGUEPLOT? Classes Outlier, OutlierAnom, Emax
        // Not kept ones (Deviant, EmaxOK)
        if (rogueplot.IsPlot() && outlier) {
          ObservationStatus status = obs.ObsStatus();
          int pclass = 0;
          if (status.TestOutlier())     {pclass = 1;}
          if (status.TestOutlierAnom()) {pclass = 2;}
          if (status.TestEmax())        {pclass = 3;}
          rogueplot.PlotOutlier(obs.GetS(), pclass);  // diffraction vector (rlu)
        }
      }
    }  // end loop observations
    SelectedObservations allobs(this_refl, -1, ALL);
    SelectedObservations plusobs(this_refl, -1, IPLUS);
    SelectedObservations minusobs(this_refl, -1, IMINUS);
    fprintf(rogues,"              Weighted mean, sd%9d%6d  I+ %9d%6d  I- %9d%6d\n",
            Nint(allobs.Average().I()), Nint(allobs.Average().sigI()),
            Nint(plusobs.Average().I()), Nint(plusobs.Average().sigI()),
            Nint(minusobs.Average().I()), Nint(minusobs.Average().sigI()));
  }
}  // namespace scala
