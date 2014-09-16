//
// reject.cpp
//
// set or reset outlier flags for all reflections

#include "reject.hh"
#include "selectedobservations.hh"
#include "file_util.hh"
#include "observationstatuscontrol.hh"
#include "string_util.hh"

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
    {deviations.assign(This_Refl.num_observations(), 0.0);}  // clear deviation list

    // Check for outliers in specified class, flags with status, accumulate flags
    void Check(const AnomalousClass& selclass, const int& Dts_Index,
               const ObservationStatus& status);

    // Check all observations against Emax test
    // Accumulate indices to outlier rejects from vector "selected" into vectors
    // "rejected " and "statusflags". Status for rejects is in "status"
    void CheckEmax(const Normalise& NormRes, const EProb& eprobtest,
                   const bool& Centric, const ObservationStatus& status);

    const std::vector<int>& Rejected() const  {return rejected;}
    const std::vector<ObservationStatus>& Statusflags() const {return statusflags;}
    std::vector<float> Deviations() const {return deviations;}
    //! Return true if outliers found even if not rejected
    bool Discrepant() const {return discrepant;}

    // Update reflection with new observation status flags
    void UpdateReflection();

  private:
    reflection* this_refl;
    const OutlierControl* outliercontrol;

    std::vector<int> rejected;
    std::vector<ObservationStatus> statusflags;
    std::vector<float> deviations;
    bool discrepant;
  };
  // ------------------------------------------------------------
  void RejectList::Check(const AnomalousClass& selclass, const int& dts_index,
                         const ObservationStatus& status)
  // Select observations according to dataset index dts_index and anomalous class
  // Accumulate indices to outlier rejects from vector "selected" into vectors
  // "rejected " and "statusflags". Status for rejects is in "status"
  //
  {
    SelectedObservations sel(*this_refl, dts_index, selclass);
    RejectFlags rejflags = outliercontrol->Reject(selclass, dts_index);
    std::vector<int> outlierindexlist =
      sel.OutlierIndexList(rejflags);
    discrepant = discrepant || sel.Discrepant();

    if (rejected.size() == 0) {
      rejected.assign(outlierindexlist.begin(),
                      outlierindexlist.end());
      statusflags.assign(outlierindexlist.size(), status);
    } else {
      rejected.insert(rejected.end(),
                      outlierindexlist.begin(),
                      outlierindexlist.end());
      statusflags.insert(statusflags.end(), outlierindexlist.size(),
                         status);
    }
    std::vector<float> ddi = sel.DeltaAll();
    ASSERT (ddi.size() == deviations.size());
    for (size_t i=0;i<ddi.size();++i) {
      if (ddi[i] != 0.0) {deviations[i] = ddi[i];}  // copy deviations
    }
  }
  // ------------------------------------------------------------
  void RejectList::CheckEmax(const Normalise& NormRes, const EProb& eprobtest,
                             const bool& Centric, const ObservationStatus& status)
  // Check all observations against Emax test
  // Accumulate indices to outlier rejects from vector "selected" into vectors
  // "rejected " and "statusflags". Status for rejects is in "status"
  //
  {
    SelectedObservations sel(*this_refl, -1, ALL);
    std::vector<int> outlierindexlist =
      EmaxRejectIndexList(sel, NormRes, eprobtest, Centric);
    if (outlierindexlist.size() > 0) {
      if (rejected.size() == 0) {
        rejected.assign(outlierindexlist.begin(),
                        outlierindexlist.end());
        statusflags.assign(outlierindexlist.size(), status);
      } else {
        rejected.insert(rejected.end(),
                        outlierindexlist.begin(),
                        outlierindexlist.end());
        statusflags.insert(statusflags.end(), outlierindexlist.size(),
                           status);
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
        this_obs.UpdateStatus(statusflags[i]);
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

    hkl_list.rewind();

    while (hkl_list.next_reflection(this_refl) >= 0)  {  // loop reflections
      bool Centric = hkl_list.symmetry().is_centric(this_refl.hkl());
      temp_refl = this_refl;  // copy
      // start reject list, store pointers to temp_refl and outliercontrol
      RejectList rejlist(temp_refl, outliercontrol);
      //  Apply current SD correction to reflection (all observations)
      SDM.CorrectReflection(temp_refl);

      //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // Outlier check
      // Loop datasets
      for (int idts=0;idts<ndsetscheck;++idts) {
        int dts_index = idts;
        // if OUTLIER COMBINE option, check outliers across all datasets
        if (outliercontrol.Combine()) dts_index = -1;

        if (anomOn && !Centric) {
          // Anomalous, check I+ & I- seperately
          rejlist.Check(IPLUS, dts_index, ObservationStatus::OBSSTAT_OUTLIER);
          rejlist.Check(IMINUS,dts_index, ObservationStatus::OBSSTAT_OUTLIER);
          if (outliercontrol.Anom()) {
            // check for outliers between I+ & I-
            rejlist.UpdateReflection();  // update flags for rejections within I+/-
            rejlist.Check(BOTH, dts_index, ObservationStatus::OBSSTAT_OUTLIERANOM);
          }
        } else {
          // no anomalous, just check all observations
          rejlist.Check(ALL, dts_index, ObservationStatus::OBSSTAT_OUTLIER);
        }
      }  // end loop datasets
      //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // Emax test
      if (!(eprobtest.Null())) {
        rejlist.CheckEmax(NormRes, eprobtest, Centric, ObservationStatus::OBSSTAT_EMAX);
      }
      //  put rejected list in "rejected" & "status"
      //   list of rejected observations ...
      std::vector<int> rejected(rejlist.Rejected().begin(), rejlist.Rejected().end());
      //   ... and their status
      std::vector<ObservationStatus> statusflags(rejlist.Statusflags().begin(),
                                                 rejlist.Statusflags().end());

      ASSERT (rejected.size() == statusflags.size());
      if (rejlist.Discrepant() > 0) {
        if (rejected.size() > 0) {
          // Set outlier flags back into observations within reflections
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
  std::vector<int> EmaxRejectIndexList
  (const SelectedObservations& selobs,
   const Normalise& NormRes, const EProb& eprobtest,
   const bool& Centric)
  // Return list of index numbers for each Emax outlier observation, if any
  {

    reflection this_ref = selobs.Reflection();
    Rtype invresolsq = this_ref.invresolsq();
    observation this_obs;

    std::vector<int> idxlist;
    int i;
    while ((i = selobs.next_observation(this_obs)) >= 0) {
      this_obs =this_ref.get_observation(i);
      float E2 = NormRes.applyAvg(this_obs.kI(), invresolsq);
      if (eprobtest.TooBig(E2, Centric)) {
        // reject
        idxlist.push_back(i);
      }
    }
    return idxlist;
  }
  // ------------------------------------------------------------
  // ------------------------------------------------------------
  WriteRogues::WriteRogues(const bool& Start, const bool& Plot,
                           const bool& multilattice,
                           const std::string& title, const float& dstarMax,
                           const float& wavelength,
                           const OutlierControl& outliercontrol)
  // Open ROGUES file & write header if Start true
  // Open ROGUESPLOT file & write header if Plot true
  // multilattice = true is there are multiple lattices
  // title & maximum resolution d* = lambda/d
  //  outliercontrol   parameters for rejection
  {
    if (Start) {
      rogues = OpenFile("ROGUES", true);  // open ROGUES file
      fprintf(rogues,
              "The ROGUES file contains all rejected reflections ");
      std::string rs =
        std::string("\nRej = '*', '@' for I+- rejects, '#' for Emax rejects, ")+
        "'x' for accepted flagged observation";
      if (multilattice) {
        rs += ",\n      'M' for multiple lattice overlaps";
      }
      rs += "\n";
      fprintf(rogues, rs.c_str());
      fprintf(rogues,
              "TotFrc = total fraction, fulls (f) or partials (p),");
      fprintf(rogues,
              " Bijv I+ or I- for Bijvoet classes\n");
      fprintf(rogues,
              "DelI/sd = (Ihl - Mn(I)others)/sqrt[sd(Ihl)**2 + sd(Mn(I))**2]\n\n");
      if (multilattice) {
        rs = std::string("Note that multilattice overlapped observations are not used in outlier calculation nor in means,\n")+
          "  and are listed here only under one of their hkl indices\n\n";
        fprintf(rogues, rs.c_str());
      }
      fprintf(rogues,
 "Flagged observations kept are labelled as: B BGratio; P PKratio; N TooNeg; G BGgradient; O Overload; E Edge\n");
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
      rogueplot = RoguePlot("ROGUEPLOT", title, dstarMax, wavelength);
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
    char reject;
    std::string PlusMinus;
    std::string flagtype = "";
    fprintf(rogues,"\n");
    float invresolsq = this_refl.invresolsq();
    float d = 1./sqrt(invresolsq);

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
          reject = 'x';
        }
        if (status.TestOutlier()) {
          reject = '*';
          outlier = true;
        } else if (status.TestOutlierAnom()) {
          reject = '@';
          outlier = true;
        } else if (status.TestEmax()) {
          reject = '#';
          outlier = true;
        }
        if (obs.Isym()%2 == 0) {
          PlusMinus = "I-";
        } else {
          PlusMinus = "I+";
        }
        if (!obs.IsSingleton()) {
          // multiple lattice overlap
          reject = 'M';
        }
        float scale = obs.Gscale();
        if (scale != 0.0) scale = 1./scale;
        std::pair<float,float> XY = obs.XYdet();

        float E = NormRes.applyAvg(obs.kI(), invresolsq);
        if (E > 0.0) {E = sqrt(E);}
        else {E = -sqrt(-E);}
        std::string flagtype = obs.Observationflag().format();

        fprintf(rogues,
        "%4d%4d%4d  %4d%4d%4d%6d%8d%6d%6.2f%6.1f%c   %2s%7.3f%7.1f%6.2f%7.1f%7.1f%7.1f %5.4f %c  %s\n",
                obs.hkl_original()[0],obs.hkl_original()[1],obs.hkl_original()[2],
                this_refl.hkl()[0],this_refl.hkl()[1],this_refl.hkl()[2],
                obs.Batch(), Nint(obs.kI()), Nint(obs.ksigI()),
                E, obs.TotalFraction(),
                partial,  PlusMinus.c_str(), scale, deviations[lobs],
                d, XY.first, XY.second, obs.phi(), obs.LP(), reject, flagtype.c_str());
        // ROGUEPLOT?
        if (rogueplot.IsPlot() && outlier) {
          rogueplot.PlotOutlier(obs.GetS());  // diffraction vector (rlu)
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
