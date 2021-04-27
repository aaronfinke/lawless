//  sdmodel.cpp

#define ASSERT assert
#include <assert.h>

#include "sdmodel.hh"
#include "selectedobservations.hh"
#include "string_util.hh"
#include "restore.hh"

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;

namespace scala
{
//--------------------------------------------------------------
SDmodel CreateSDmodel(const phaser_io::InputAll& input,
                      const std::vector<Run>& runlist,
                      const bool& setnull)
// Create SDmodel for each run from input or defaults
{
  SDcorrection sdcdeffull(1.0,0.0,0.02);
  if (setnull) {
    sdcdeffull = SDcorrection(1.0,0.0,0.0);
  }
  sdcdeffull.SetDefaultRestraints();

  SDcorrection sdcdefpartial = sdcdeffull;  // default values
  SDcorrection sdcfull = sdcdeffull;
  SDcorrection sdcpartial = sdcdefpartial;  // actual values

  SDmodel SDM;
  int Nruns = runlist.size();
  int nrunsused = Nruns;

  // Refine flag
  SDM.SetRefine(input.SDC_Refine());
  SDM.SetDamp(input.SDCdamp());

  const int MINIMUMSAMPLE = 6; // minimum number for sample SD
  SDM.SetSampleSD(input.SampleSD(), MINIMUMSAMPLE);
  // If SampleSD, switch off refinement unless explicit  DON'T
  //  if(SDM.SampleSD()) {
  //    if (!input.SDC_RefineSet()) {
  //      SDM.SetRefine(false);
  //    }}

  std::vector<std::pair<SDcorrection,SDcorrection> >
    sdcval = input.SDC_SDcorrections();
  std::vector<int> runnumbers = input.SDC_RunNumbers();
  ASSERT (sdcval.size() == runnumbers.size());

  if (input.SDC_NumberInput() < 0) {
    // Overall values for all runs given
    sdcfull = sdcval[0].first;
    sdcpartial = sdcval[0].second;
  }

  // All runs same flag
  bool allrunssame = input.SDC_AllRunsSame();
  if (Nruns == 1) allrunssame = false;
  if (allrunssame) nrunsused = 1;

  std::vector<Run::FullsAndPartials> FandP(Nruns);

  for (int irun=0;irun<Nruns;irun++) { // loop all runs
    if (input.SDC_NumberInput() > 0) {
      bool found = false;
      for (size_t i=0;i<runnumbers.size();++i) {
        if (runnumbers[i] == runlist[irun].RunNumber()) {
          sdcfull = sdcval[i].first;
          sdcpartial = sdcval[i].second;
          found = true;
          break;
        }
      }
      if (!found) {
        sdcfull = sdcdeffull;
        sdcpartial = sdcdefpartial;
      }
    }
    SDM.AddRun(runlist[irun].RunNumber(), sdcfull, sdcpartial);
    if (!allrunssame) { // not all same, set flags for all runs
      SetSdmFullPartialFlags(runlist[irun].fullsAndPartials(), SDM, irun);
    } else { // all the same, store flags
      FandP[irun] = runlist[irun].fullsAndPartials();
    }
  } // end loop runs
  SDM.SetAllRunsSame(allrunssame);
  // Fix SdB flag
  bool fixSDb = input.SDC_NoSDb();
  SDM.SetNoSDb(fixSDb);

  //^^^^
  //      SDM.FullAsPartial(0, true);   // full as partial for test FIXME
  //      std::cout << "SDM full as partial!!\n";
  //  SDM.PartialAsFull(0, true);   //  partial as full for test FIXME
  //  std::cout << "SDM partial as full!!\n";
  //^^^^-
  if (allrunssame) { // all runs same, find overall setting for ONLY|FEW Fulls|PARTIALS
    // only here if >1 run
    Run::FullsAndPartials FP1;
    bool sameflag = true; // all have the same flag
    bool anyboth = false;
    bool anyFewFull = false;
    bool anyFewPartial = false;

    for (int irun=0;irun<Nruns;irun++) { // loop runs
      if (irun == 0) {
        FP1 = FandP[irun];
      }
      if (FandP[irun] != FP1) sameflag = false;
      if (FandP[irun] == Run::FULLSANDPARTIALS) anyboth = true;
      if (FandP[irun] == Run::FEWFULLS) anyFewFull = true;
      if (FandP[irun] == Run::FEWPARTIALS) anyFewPartial = true;
    } // end loop runs
    // Set flags for all runs, not just the first, even though this is allrunssame
    for (int irun=0;irun<Nruns;irun++) { // loop runs
      if (sameflag) {
        // all same flags, use that one
        SetSdmFullPartialFlags(FP1, SDM, irun);
      } else if (anyboth) {
        // any have both fulls & partials (not few), leave as that default
      } else if (anyFewFull && !anyFewPartial) {
        // use FEWFULLS unless FEWPARTIALS is also set
        SetSdmFullPartialFlags(Run::FEWFULLS, SDM, irun);
      } else if (!anyFewFull && anyFewPartial) {
        // use FEWPARTIALS unless FEWFULLS is also set
        SetSdmFullPartialFlags(Run::FEWPARTIALS, SDM, irun);
      }
    }  // end loop runs
  } // allrunssame

  std::vector<double> targets(3);    // 3 targets
  std::vector<double> sdtargets(3);  // ... and their SDs (= 0 no target)
  // = 0 no tie, = -1 defaults, = +1 set from input, = +2 similarity tie
  int tietype = input.SDCties(targets, sdtargets);
  SDM.SetTies(tietype, targets, sdtargets);  // for all SD corrections

  return SDM;
}
//-------------------------------------------------------------
  void SetSdmFullPartialFlags(const Run::FullsAndPartials& FandP,
                              SDmodel& SDM, const int& irun)
  // Set appropriate flags into SDM...[irun] according to FandP
  {
    if (FandP == Run::ONLYFULLS) {
      // no partials
      SDM.PartialAsFull(irun, false); // use values from fulls for partials
    } else if (FandP == Run::FEWFULLS) {
      // few fulls
      SDM.FullAsPartial(irun, true); // use values from partials for fulls
    } else if (FandP == Run::ONLYPARTIALS) {
      // no fulls
      SDM.FullAsPartial(irun, false); // use values from partials for fulls
    } else if (FandP == Run::FEWPARTIALS) {
      // few partials
      SDM.PartialAsFull(irun, true); // use values from fulls for partials
    }
  }
//-------------------------------------------------------------
  void SDmodel::init()
  {
    nosdb = false;
    allrunssame = false;
    nsets = 0;
    ties.tietype = 0; // no ties
    ties.targets.assign(3,0.0);
    ties.sdtargets.assign(3,0.0);
    SetVarianceWeights();    // default weighting scheme
    sampleSD = false;
    minimumsample = 10;
}
//-------------------------------------------------------------
  void SDmodel::ResetRange()
  // Reset all minimum & maximum values
  {
    for (size_t i=0;i<sdc_full_run.size();++i) {
      sdc_full_run[i].ResetRange();
    }
    for (size_t i=0;i<sdc_partial_run.size();++i) {
      sdc_partial_run[i].ResetRange();
    }
  }
//-------------------------------------------------------------
  void SDmodel::AddRun(const int& runNum,
                       const SDcorrection& SDCfull, const SDcorrection& SDCpartial)
 {
   sdc_full_run.push_back(SDCfull);
   sdc_full_run.back().SetFixSDb(nosdb);
   sdc_partial_run.push_back(SDCpartial);
   sdc_partial_run.back().SetFixSDb(nosdb);
   usetype.push_back(0); // fulls and partials
   runnumbers.push_back(runNum);
   ASSERT (sdc_full_run.size() == sdc_partial_run.size());
   nsets = sdc_full_run.size();
   SetIdxParam();  // set up index list
 }
  //-------------------------------------------------------------
  //! Store ties for all SD corrections
  // tietype = 0 no tie, = -1 defaults, = +1 set from parameters, = +2 similarity
  void SDmodel::SetTies(const int& Tietype,
                        const std::vector<double>& Targets,
                        const std::vector<double>& SDtargets)
  {
    ties.tietype = Tietype;
    ties.targets = Targets;
    ties.sdtargets = SDtargets;
    SetTies();
  }
  //-------------------------------------------------------------
  //! Store one tie for all SD corrections, for parameter ipar (0-2)
  void SDmodel::ResetTie(const int& ipar,
                        const double& Target,
                        const double& SDtarget)
  {
    ties.tietype = +1;
    ASSERT (ipar < int(ties.targets.size()));
    ties.targets.at(ipar) = Target;
    ties.sdtargets.at(ipar) = SDtarget;
    SetTies();
  }
//-------------------------------------------------------------
  //! Store ties for all SD corrections
  // tietype = 0 no tie, = -1 defaults, = +1 set from parameters
  void SDmodel::SetTies()
  {
    //^
    //    std::cout << "SDmodel::SetTies " <<
    //      " " << ties.targets[0]<< " " << ties.targets[1]
    //        << " " << ties.targets[2] <<"\n"
    //        << " " << ties.sdtargets[0]<< " " << ties.sdtargets[1]
    //        << " " << ties.sdtargets[2] <<"\n";
    //^-
    for (int irun=0;irun<Nruns();++irun) {
      if (ties.tietype == 0) { // no ties
        sdc_full_run[irun].ClearRestraints();
        sdc_partial_run[irun].ClearRestraints();
      } else if (ties.tietype < 0) { // use defaults
        sdc_full_run[irun].SetDefaultRestraints();
        sdc_partial_run[irun].SetDefaultRestraints();
        // Store default parameters for printing
        sdc_full_run[irun].GetRestraints(ties.targets, ties.sdtargets);
      } else { // use input values or similarity
        sdc_full_run[irun].SetRestraints(ties.targets, ties.sdtargets);
        sdc_partial_run[irun].SetRestraints(ties.targets, ties.sdtargets);
      }
    }
  }
  //-------------------------------------------------------------
  //! true if there are active ties
  bool SDmodel::restrained() const
  {
    bool restraints = false;  // true if we have some restraints
    if (ties.tietype != 0) {  // ties.tietype = 0 for no restraints
      // If we are not refining SdB and the only restraint is on SdB, then no restraints
      if (nosdb) {
        // No SdB refinement
        if ((ties.sdtargets[0] != 0.0) || (ties.sdtargets[2] != 0.0)) {
          restraints = true;
        }
      } else { // SdB refinement
        if ((ties.sdtargets[0] != 0.0) || (ties.sdtargets[1] != 0.0) ||
            (ties.sdtargets[2] != 0.0)) {
          restraints = true;
        }
      }
    }
    return restraints;
  }
  //-------------------------------------------------------------
  //! format tie information
  std::string SDmodel::formatTie() const
  {
    // NB all the same
    std::string s;
    bool restraints = restrained();  // true if we have some restraints

    if (ties.tietype == +2) {
      s = "SD parameters tied to average across all runs";
    } else if (restraints) {
      s = "Restraints on SD correction parameters (target (+-SD)):";
      if (ties.sdtargets[0] != 0.0) { // SdAdd
        s += " SdAdd "+StringUtil::Strip(StringUtil::ftos(ties.targets[0],5,1))+
          " (+-"+StringUtil::Strip(StringUtil::ftos(ties.sdtargets[0],5,1))+")";
      }
      if (!nosdb && ties.sdtargets[1] != 0.0) { // SdB
        s += " SdB "+StringUtil::Strip(StringUtil::ftos(ties.targets[1],8,1))+
          " (+-"+StringUtil::Strip(StringUtil::ftos(ties.sdtargets[1],8,1))+")";
      }
      if (ties.sdtargets[2] != 0.0) { // SdAdd
        s += " SdAdd "+StringUtil::Strip(StringUtil::ftos(ties.targets[2],8,3))+
          " (+-"+StringUtil::Strip(StringUtil::ftos(ties.sdtargets[2],8,3))+")";
      }
    } else {
      s = "No restraints on SD correction parameters";
    }
    return s;
  }
  //-------------------------------------------------------------
  //! Set flag = true to fix SDb = 0.0
  void SDmodel::SetNoSDb(const bool& NoSDb)
  {
    nosdb = NoSDb;
    for (int irun=0;irun<Nruns();++irun) {
      if (nosdb) {
        sdc_full_run[irun].FixSDb();
        sdc_partial_run[irun].FixSDb();
      } else {
        sdc_full_run[irun].UnFixSDb();
        sdc_partial_run[irun].UnFixSDb();
      }
    }
    SetIdxParam();  // set up index list
  }
  //-------------------------------------------------------------
  //! Set flag = true to use same parameters for all runs
  void SDmodel::SetAllRunsSame(const bool& flag)
  {
    allrunssame = flag;
    if (allrunssame) {
      nsets = 1;
    } else {
      nsets = Nruns();
    }
  }
  //-------------------------------------------------------------
  //! Number of parameters, depends on allrunsame flag and full/partial usage
  int SDmodel::Nparams() const
  {
    int nparams = 0;
    for (int irun=0;irun<nsets;++irun) {
      if (usetype[irun] >= 0) { // use parameters for fulls
        nparams += sdc_full_run[irun].Nparams();
      }
      if (usetype[irun] <= 0) { // use parameters for partials
        nparams += sdc_partial_run[irun].Nparams();
      }
    }
    return nparams;
  }
  //-------------------------------------------------------------
  //! Number of parameter groups, depends on allrunsame flag and full/partial usage
  int SDmodel::Ngroups() const
  {
    int ngroups = 0;
    for (int irun=0;irun<nsets;++irun) {
      if (usetype[irun] >= 0) { // use parameters for fulls
        ngroups++;
      }
      if (usetype[irun] <= 0) { // use parameters for partials
        ngroups++;
      }
    }
    return ngroups;
  }
  //-------------------------------------------------------------
  //! Set index list for each run
  void SDmodel::SetIdxParam()
  {
    int nparams = 0;
    idxfullparam.assign(Nruns(), 0);
    idxpartialparam.assign(Nruns(), 0);
    idxparamgroups.assign(Nruns(), std::pair<int,int>(-1,-1));
    std::vector<int> nparamsinallgroups;
    int jpargroup = 0;

    for (int irun=0;irun<nsets;++irun) {
      if (usetype[irun] >= 0) { // use parameters for fulls
        idxfullparam[irun] = nparams; // index to first parameter in run
        nparams += sdc_full_run[irun].Nparams();
        nparamsinallgroups.push_back(sdc_full_run[irun].Nparams());
        idxparamgroups[irun].first = jpargroup++;
      }
      if (usetype[irun] <= 0) { // use parameters for partials
        idxpartialparam[irun] = nparams; // index to first parameter in run
        nparams += sdc_partial_run[irun].Nparams();
        nparamsinallgroups.push_back(sdc_partial_run[irun].Nparams());
        idxparamgroups[irun].second = jpargroup++;
      }
    }
    if (nsets < Nruns()) {
      // if allsame, copy parameter indices from 1st run
      for (int irun=nsets;irun<Nruns();++irun) {
        if (usetype[irun] >= 0) { // use parameters for fulls
          idxfullparam[irun] = idxfullparam[0]; // index to first parameter in run
          idxparamgroups[irun].first = idxparamgroups[0].first;
        }
        if (usetype[irun] <= 0) { // use parameters for partials
          idxpartialparam[irun] = idxpartialparam[0]; // index to first parameter in run
          idxparamgroups[irun].second = idxparamgroups[0].second;
        }
      }
    }
    // Check all groups have the same number of parameters
    ASSERT (nparamsinallgroups.size() > 0);
    for (size_t k=1; k<nparamsinallgroups.size(); k++) {
      ASSERT (nparamsinallgroups[k] == nparamsinallgroups[0]);
    }
    nparamspergroup = nparamsinallgroups[0];
  }
//-------------------------------------------------------------
  //! use partial correction for fulls
  // if few == true, there are few fulls, if false none
  void SDmodel::FullAsPartial(const int& RunIndex, const bool& few)
  {
    int flag = few ? -2 : -1;
    usetype[RunIndex] = flag;
    SetIdxParam();
  }
//-------------------------------------------------------------
  //! use full correction for partials
  // if few == true, there are few partials, if false none
  void SDmodel::PartialAsFull(const int& RunIndex, const bool& few)
  {
    int flag = few ? +2 : +1;
    usetype[RunIndex] = flag;
    SetIdxParam();
  }
  //-------------------------------------------------------------
  //! return true if at least one run has few fulls||partials
  bool SDmodel::SomeFewRuns() const
  {
    for (int irun=0;irun<Nruns();++irun) {
      if (std::abs(usetype[irun]) == 2) return true;
    }
    return false;
  }
  //-------------------------------------------------------------
  void SDmodel::UpdateFactor(const int& RunIndex, const bool& copy,
                             const float& UpdateFull, const float& UpdatePartial)
  // Update SDC for run RunIndex by multiplying SDfac values, only if > 0
  //    UpdateFull     update factor for fulls
  //    UpdatePartial  update factor for partials
  // If copy is true, then copy parameters  full <-> partial according to usetype
  // If allrunssame, apply to all runs
  {
    ASSERT (RunIndex < int(sdc_full_run.size()));
    int ir1 = RunIndex;
    int ir2 = RunIndex+1;
    if (allrunssame) {
      ir1 = 0; ir2 = Nruns();
    }
    for (int ir=ir1;ir<ir2;++ir) { // just RunIndex or loop all if allrunssame
      if (UpdateFull > 0.0) {sdc_full_run[ir].UpdateFactor(UpdateFull);}
      if (UpdatePartial > 0.0) {sdc_partial_run[ir].UpdateFactor(UpdatePartial);}
      if (copy) {
        // copy values if needed for usetype flags
        if (usetype[ir] < 0) {
          sdc_full_run[ir] = sdc_partial_run[ir]; // full as partial
        } else if (usetype[ir] > 0) {
          sdc_partial_run[ir] = sdc_full_run[ir]; // partial as full
        }
      }
    }
  }
  //-------------------------------------------------------------
  void SDmodel::SetSDadd(const double& SDadd)
  // Reset SDadd for all runs
  {
    for (int ir=0;ir<Nruns();++ir) { // loop all runs
      sdc_full_run[ir].SetSDadd(SDadd);
      sdc_partial_run[ir].SetSDadd(SDadd);
    }
  }
//-------------------------------------------------------------
  std::vector<float> SDmodel::CorrectReflection(reflection& Ref) const
  // Apply appropriate SD correction to all valid observations in reflection
  // Return vector of uncorrected sigmas, for all observations, even unselected ones
  // Use full/partial values irrespective of usetype, assuming that values
  // have been duplicated if necessary
  {
    // number of observations in this including unselected ones
    int nobs = Ref.num_observations();
    std::vector<float> sig0(nobs,0.0); // uncorrected sigma(I)
    SelectedObservations selobs(Ref, -1, ALL, weighttype);
    float Iav = selobs.Average().I();  // average intensity for SD correction

    observation this_obs;
    Ref.reset(); // reset next_observation count
    int iobs; // observation index in reflection

    while ((iobs=Ref.next_observation(this_obs)) >= 0) {
      // correct sigI for observation, return uncorrected value
      if (this_obs.IsFull())
        {sig0[iobs] = sdc_full_run[this_obs.run()].Correct(this_obs, Iav);}
      else
        {sig0[iobs] = sdc_partial_run[this_obs.run()].Correct(this_obs, Iav);}
      // store updated observation
      Ref.replace_observation(this_obs);
    }
    return sig0;
  }
//-------------------------------------------------------------
  void SDmodel::CorrectAllReflection(reflection& Ref) const
  // Apply appropriate SD correction to all observations in reflection, including outliers
  // Use full/partial values irrespective of usetype, assuming that values
  // have been duplicated if necessary
  {
    // number of observations in this including unselected ones
    int nobs = Ref.num_observations();
    SelectedObservations selobs(Ref, -1, ALL, weighttype);
    float Iav = selobs.Average().I();  // average intensity for SD correction (omitting rejects

    observation this_obs;

    for (int iobs=0;iobs<nobs;++iobs) {
      this_obs = Ref.get_observation(iobs);
      ObservationStatus obsstatus = this_obs.ObsStatus() ;
      // OK is Accepted or Outlier
      if (obsstatus.IsOKforRogues()) {
        // correct sigI for observation, return uncorrected value
        if (this_obs.IsFull())
          {sdc_full_run[this_obs.run()].Correct(this_obs, Iav);}
        else
          {sdc_partial_run[this_obs.run()].Correct(this_obs, Iav);}
        // store updated observation
        Ref.replace_observation(this_obs);
      } // accepted
    }
  }
//-------------------------------------------------------------
  // Return overall minimum & maximum values
  void SDmodel::GetSDcorrectionRanges(float& minSDcorrFulls, float& maxSDcorrFulls,
                                      float& minSDcorrPartials, float& maxSDcorrPartials) const
  {
    std::pair<float,float> minmax;
    for (int i=0;i<Nruns();++i) {  // loop runs for fulls
      if (i == 0) {
        minmax = sdc_full_run[i].MinMax();
      } else {
        minmax.first = Min(minmax.first, sdc_full_run[i].MinMax().first);
        minmax.second = Max(minmax.second, sdc_full_run[i].MinMax().second);
      }
    }
    minSDcorrFulls = 0.0;
    maxSDcorrFulls = 0.0;
    if (minmax.first < 1.0e+9) {
      minSDcorrFulls = minmax.first;
    }
    if (minmax.second > -1.0e+9) {
      maxSDcorrFulls = minmax.second;
    }

    for (size_t i=0;i<sdc_partial_run.size();++i) { // loop runs for partials
      if (i == 0) {
        minmax = sdc_partial_run[i].MinMax();
      } else {
        minmax.first = Min(minmax.first, sdc_partial_run[i].MinMax().first);
        minmax.second = Max(minmax.second, sdc_partial_run[i].MinMax().second);
      }
    }
    minSDcorrPartials = 0.0;
    maxSDcorrPartials = 0.0;
    if (minmax.first < 1.0e+9) {
      minSDcorrPartials = minmax.first;
    }
    if (minmax.second > -1.0e+9) {
      maxSDcorrPartials = minmax.second;
    }
  }
  //-------------------------------------------------------------
  // Get vector of parameters
  std::vector<double> SDmodel::GetParameters() const
  {
    std::vector<double> params;
    std::vector<double> p;

    // nsets = Nruns, or 1 if allrunssame
    for (int irun=0;irun<nsets;++irun) {
      if (usetype[irun] >= 0) { // values for fulls
        p = sdc_full_run[irun].GetParameters();
        params.insert(params.end(), p.begin(), p.end());
      }
      if (usetype[irun] <= 0) { // values for partials
        p = sdc_partial_run[irun].GetParameters();
        params.insert(params.end(), p.begin(), p.end());
      }
    }
    return params;
  }
  //-------------------------------------------------------------
  // Get vector of parameter "shifts" ie starting values for varying
  // the parameters in simplex optimisation
  std::vector<double> SDmodel::GetShifts(const double& scale) const
  {
    std::vector<double> shifts;
    std::vector<double> p;
    for (int irun=0;irun<nsets;++irun) {
      if (usetype[irun] >= 0) { // values for fulls
        p = sdc_full_run[irun].GetShifts(scale);
        shifts.insert(shifts.end(), p.begin(), p.end());
      }
      if (usetype[irun] <= 0) { // values for partials
        p = sdc_partial_run[irun].GetShifts(scale);
        shifts.insert(shifts.end(), p.begin(), p.end());
      }
    }
    return shifts;
  }
//-------------------------------------------------------------
  void SDmodel::AddToAverages(std::vector<MeanValue>& averagerealparameters,
                              const std::vector<double>& realparameters) const
  // add 2 or 3 parameters into averages
  {
    ASSERT (averagerealparameters.size() == realparameters.size());
    for (size_t i=0; i<averagerealparameters.size(); i++) {
      averagerealparameters[i].Add(realparameters[i]);
    }
  }
  //-------------------------------------------------------------
  void SDmodel::SetTargetsFromAverageParameters()
  {
    if (ties.tietype != +2) {return;} // not similarity

    int npc = Max(sdc_full_run[0].Nparams(), sdc_partial_run[0].Nparams());
    // NB refined parameters are not SdFac, SDb, SDadd but
    //  averages and restraints are
    std::vector<MeanValue> averagerealparameters(npc); // SdFac, [SDb,] SDadd

    for (int irun=0;irun<Nruns();++irun) { // loop all runs
      AddToAverages(averagerealparameters, sdc_full_run[irun].GetRealParameters());
      AddToAverages(averagerealparameters, sdc_partial_run[irun].GetRealParameters());
    } // end loop runs

    // Average parameter values
    std::vector<double> newtargets(3);  // always 3
    for (int j=0;j<npc;++j) {
      newtargets[j] = averagerealparameters[j].Mean();
    }
    if (npc == 2) { // fixSDb
      newtargets[2] = newtargets[1];
      newtargets[1] = 0.0;
    }
    ties.targets = newtargets;
    //^
    //    std::cout << "Updating SDcorrection targets "
    //        << " " << newtargets[0]
    //        << " " << newtargets[1]
    //        << " " << newtargets[2] <<"\n";
    //    std::cout << "SDmodel::SetTies " <<
    //      " " << ties.targets[0]<< " " << ties.targets[1]
    //        << " " << ties.targets[2] <<"\n"
    //        << " " << ties.sdtargets[0]<< " " << ties.sdtargets[1]
    //        << " " << ties.sdtargets[2] <<"\n";
    //^-
  }
//-------------------------------------------------------------
  // Set all parameters from vector
  void SDmodel::SetParameters(const std::vector<double>& params)
  {
    std::vector<bool> parameterupdated(params.size(), true);
    SetParameters(params, parameterupdated);
  }
//-------------------------------------------------------------
  // Set all parameters from vector
  // parameterupdated true if this parameter has been updated
  void SDmodel::SetParameters(const std::vector<double>& params,
                              const std::vector<bool>& parameterupdated)
  {
    int npc = Max(sdc_full_run[0].Nparams(), sdc_partial_run[0].Nparams());
    if (int(params.size()) != Nparams()) {
      Message::message(Message_fatal("SDmodel::SetParameters: wrong number of parameters"));
    }
    std::vector<double> pars(npc);

    int k=0;
    for (int irun=0;irun<Nruns();++irun) { // loop all runs
      if (allrunssame) k=0; // if allrunssame, one set of parameters for all runs
      int k1=k; // save
      // Always put in values for fulls even if not refined
      bool updatedfull = false;
      for (int j=0;j<sdc_full_run[irun].Nparams();++j) {
        pars[j] = params[k];
        if (parameterupdated[k]) {updatedfull = true;}
        k++;
      }
      sdc_full_run[irun].SetParameters(pars); // set values for fulls anyway
      if (usetype[irun] != 0) { // no actual values for fulls or partials
        // pick up same parameters for partials
        k = k1;
      }
      // Always put in values for partials even if not refined
      bool updatedpartial = false;
      for (int j=0;j<sdc_partial_run[irun].Nparams();++j) {
        pars[j] = params[k];
        if (parameterupdated[k]) {updatedpartial = true;}
        k++;
      }
      sdc_partial_run[irun].SetParameters(pars); // set values for partials
      // Copy full <-> partial if needed
      if (updatedfull && !updatedpartial) {
        //  partial from full
        sdc_partial_run[irun].SetParameters(sdc_full_run[irun].GetParameters());
      }
      if (!updatedfull && updatedpartial) {
        // full from partial
        sdc_full_run[irun].SetParameters(sdc_partial_run[irun].GetParameters());
      }
    } // end loop runs

    SetTargetsFromAverageParameters();
  }
  //-------------------------------------------------------------
  std::vector <std::vector<double> >
  SDmodel::GetDerivatives(SelectedObservations& selobs,
                          const std::vector<float>& sigmaI) const
  // uncorrected scaled sigma(I) for each observation (including unselected ones)
  // return vector elements for each observation in selobs
  // each element is vector of elements for each parameter
  //  elements for each parameter are d(delta(iobs))/dp(k)
  {
    ASSERT (selobs.Nobs() == int(sigmaI.size()));
    // corrected sigma' for each observation (scaled)
    std::vector<float> sigmaprime = selobs.sigmaI();
    // delI for each observation (scaled by 1/g) = Ihl - <Ih>
    std::vector<float> delI = selobs.DelI();
    // number of observations in selobs including unused slots
    int nobs = delI.size();
    ASSERT (nobs == int(sigmaprime.size()));
    // d(delta(i))/d(p(j)) for all observations i
    //  empty observation slots will contain empty vectors
    std::vector <std::vector<double> >  ddeltadp(nobs);

    std::vector<double> ddeltaidp; // d(delta(i))/dp for i'th observation
    float Iav = selobs.Average().I();  // average intensity for SD correction
    std::vector<double> dp;  // for each param set
    double an = selobs.Number();  // number used
    double fac = -sqrt(an/(an-1.0));

    for (int iobs=0;iobs<nobs;++iobs) { // loop observations in selobs
      if (delI[iobs] != 0.0) {  // a valid delta
        ddeltaidp.assign(Nparams(),0.0);
        //--- calculate d(sigma')/dp vector for all parameters
        int irun = selobs.Run(iobs);
        int idx; // first parameter index
        if (selobs.Full(iobs)) {
          // Full
          if (usetype[irun] >= 0) { // values for fulls
            dp = sdc_full_run[irun].GetDerivatives(sigmaI[iobs], Iav);
            idx = idxfullparam[irun];
          } else if (usetype[irun] < 0) { // full as partial
            dp = sdc_partial_run[irun].GetDerivatives(sigmaI[iobs], Iav);
            idx = idxpartialparam[irun];
          }
        } else {
          // Partial
          if (usetype[irun] > 0) { // partial as full
            dp = sdc_full_run[irun].GetDerivatives(sigmaI[iobs], Iav);
            idx = idxfullparam[irun];
          } else if (usetype[irun] <= 0) { // partial
            dp = sdc_partial_run[irun].GetDerivatives(sigmaI[iobs], Iav);
            idx = idxpartialparam[irun];
          }
        }
        for (size_t i=0;i<dp.size();++i) {
          ddeltaidp[idx++] = dp[i];
        }
        // all parameters d(sigma')/dp done
        // uncorrected sigma(I) for observation
        //      float sigma = sigmaI[iobs];
        // d(delta)/d(sigma') = -sqrt(n/n-1) delI / (sigma')^2

        double dddsp = fac * delI[iobs] / (sigmaprime[iobs]*sigmaprime[iobs]);

        for (size_t i=0;i<ddeltaidp.size();++i) {
          ddeltaidp[i] *= dddsp; // d(delta)/dp = d(delta)/d(sigma') d(sigma')/dp
        }
        ddeltadp[iobs] = ddeltaidp;  // store d(delta(i))/dp vector
      } // end valid delta
    } // end loop observations in selobs
    return ddeltadp;
  }
  //-------------------------------------------------------------
  std::vector<double> SDmodel::LowerBounds() const
  //! lower bounds for each parameter (0.0 means no bound)
  {
    std::vector<double> bounds;
    std::vector<double> p;
    // nsets = Nruns, or 1 if allrunssame
    for (int irun=0;irun<nsets;++irun) {
      if (usetype[irun] >= 0) { // values for fulls
        p = sdc_full_run[irun].LowerBounds();
        bounds.insert(bounds.end(), p.begin(), p.end());
      }
      if (usetype[irun] <= 0) { // values for partials
        p = sdc_partial_run[irun].LowerBounds();
        bounds.insert(bounds.end(), p.begin(), p.end());
      }
    }
    return bounds;
  }
  //-------------------------------------------------------------
  std::vector<double> SDmodel::UpperBounds() const
  //! upper bounds for each parameter (0.0 means no bound)
  {
    std::vector<double> bounds;
    std::vector<double> p;
    // nsets = Nruns, or 1 if allrunssame
    for (int irun=0;irun<nsets;++irun) {
      if (usetype[irun] >= 0) { // values for fulls
        p = sdc_full_run[irun].UpperBounds();
        bounds.insert(bounds.end(), p.begin(), p.end());
      }
      if (usetype[irun] <= 0) { // values for partials
        p = sdc_partial_run[irun].UpperBounds();
        bounds.insert(bounds.end(), p.begin(), p.end());
      }
    }
    return bounds;
  }
  //-------------------------------------------------------------
  std::vector<double> SDmodel::LargeShifts()
  //! large shifts  for each parameter
  {
    std::vector<double> large;
    std::vector<double> p;
    // nsets = Nruns, or 1 if allrunssame
    for (int irun=0;irun<nsets;++irun) {
      if (usetype[irun] >= 0) { // values for fulls
        p = sdc_full_run[irun].LargeShifts();
        large.insert(large.end(), p.begin(), p.end());
      }
      if (usetype[irun] <= 0) { // values for partials
        p = sdc_partial_run[irun].LargeShifts();
        large.insert(large.end(), p.begin(), p.end());
      }
    }
    return large;
  }
  //-------------------------------------------------------------
  std::vector<double> SDmodel::GetRestraintR() const
  // Return restraint R2 for each parameter group
  {
    std::vector<double> R2;

    // nsets = Nruns, or 1 if allrunssame
    for (int irun=0;irun<nsets;++irun) {
      if (usetype[irun] >= 0) { // values for fulls
        R2.push_back(sdc_full_run[irun].RestraintR());
      }
      if (usetype[irun] <= 0) { // values for partials
        R2.push_back(sdc_partial_run[irun].RestraintR());
      }
    }
    ASSERT (int(R2.size()) == Ngroups());
    return R2;
  }
  //-------------------------------------------------------------
    //! return derivatives dR2/dp and Hessian for each parameter group
  void SDmodel::GetRestraintDerivatives
  (std::vector <std::vector<double> >& dr2dp,
   std::vector <clipper::Array2d<double> >& H) const
  {
    dr2dp.clear();
    H.clear();
    std::vector<double> dp;
    clipper::Array2d<double> h;

    for (int irun=0;irun<nsets;++irun) {
      if (usetype[irun] >= 0) { // values for fulls
        sdc_full_run[irun].RestraintDerivatives(dp,h);
        dr2dp.push_back(dp); // gradients
        H.push_back(h);      // Hessian
      }
      if (usetype[irun] <= 0) { // values for partials
        sdc_partial_run[irun].RestraintDerivatives(dp,h);
        dr2dp.push_back(dp); // gradients
        H.push_back(h);      // Hessian
      }
    }
    ASSERT (int(dr2dp.size()) == Ngroups());
  }
  //-------------------------------------------------------------
  //! return the "coordinate" of the iobs'th entry in selobs
  // Also return parameter group index idxgroup
  //   Iav = average I
  std::vector<double> SDmodel::coordinate(const SelectedObservations& selobs,
                                          const int& iobs, const double& Iav,
                                          int& idxgroup) const
  {
    double variance = selobs.sigmaI()[iobs]*selobs.sigmaI()[iobs];
    // Get relevant SDCmodel
    return SDCmodel(selobs, iobs, idxgroup).coordinate(variance, Iav);
  }
  //-------------------------------------------------------------
  //! return the corrected SD (from parameters) for the iobs'th entry in selobs
  //   Iav = average I
  double SDmodel::sdCorrected(const SelectedObservations& selobs,
                              const int& iobs, const double& Iav) const
  {
    double sigma = selobs.sigmaI()[iobs];
    // Get relevant SDCmodel
    int idxgroup;
    return SDCmodel(selobs, iobs, idxgroup).SigmaPrime(sigma, 1.0, Iav);
  }
  //-------------------------------------------------------------
  //! return a reference to the relevant SDC model (run, full/partial
  // Also return parameter group index idxgroup
  const SDcorrection& SDmodel::SDCmodel(const SelectedObservations& selobs,
                                        const int& iobs,
                                        int& idxgroup) const
  {
    int irun = selobs.Run(iobs);
    if (selobs.Full(iobs)) {
      // Full
      if (usetype[irun] >= 0) { // values for fulls
        idxgroup = idxparamgroups[irun].first;
        return sdc_full_run[irun];
      } else if (usetype[irun] < 0) { // full as partial
        idxgroup = idxparamgroups[irun].second;
        return sdc_partial_run[irun];
      }
    } else {
      // Partial
      if (usetype[irun] > 0) { // partial as full
        idxgroup = idxparamgroups[irun].first;
        return sdc_full_run[irun];
      } else if (usetype[irun] <= 0) { // partial
        idxgroup = idxparamgroups[irun].second;
        return sdc_partial_run[irun];
      }
    }
    ASSERT (false);  // crash
  }
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  void SDmodel::dump() const
  {
    std::cout << "Nruns: " << sdc_full_run.size() <<"\n";
    for (size_t ir=0;ir<sdc_full_run.size();++ir) {
      std::cout << "Run " << ir+1 << " useflag " << usetype[ir] << "\n";
      std::cout << "Fulls:    " << sdc_full_run[ir].format() << "\n";
      std::cout << "Partials: " << sdc_partial_run[ir].format() << "\n";
    }
  }
  //-------------------------------------------------------------
  //! return formatted only||few + fulls||partials information for all runs
  std::string SDmodel::formatFullPartialInfo() const
  {
    std::string ss;
    int nr = Nruns();
    if (allrunssame) {nr = 1;}

    for (int irun=0;irun<nr;++irun) {
      std::string label = "fulls & partials";
      if (usetype[irun] > 0) {
        label = (usetype[irun] == +1) ? "only fulls" : "relatively few partials";
      } else if (usetype[irun] < 0){
        label = (usetype[irun] == -1) ? "only partials" : "relatively few fulls";
      }
      if (allrunssame) {
      ss += FormatOutput::logTabPrintf(0,"All runs have %s\n", label.c_str());
      } else {
        ss += FormatOutput::logTabPrintf(0,"Run %4d has %s\n",
                                         irun+1, label.c_str());
      }
    }
    return ss;
  }
  //-------------------------------------------------------------
  //! return formatted use flag (only||few + fulls||partials)
  std::string SDmodel::formatUseFlag(const int& RunIndex) const
  {
    std::string label = "Fulls & partials";
    int irun = RunIndex;
    if (usetype[irun] > 0) {
      label = (usetype[irun] == +1) ? "  OnlyFulls   " : " FewPartials  ";
    } else if (usetype[irun] < 0){
      label = (usetype[irun] == -1) ? " OnlyPartials " : "  FewFulls    ";
    }
    return label;
  }
  //-------------------------------------------------------------
  std::string SDmodel::format() const
  {
    std::string ss;
    bool both = true;
    int fullpart = 0;  // >1 fulls only, < 0 partials only

    if (usetype[0] == +1) { // fulls only
      fullpart = +1;
    } if (usetype[0] == -1) { // partials only
      fullpart = -1;
    }


    both = false;
    bool somefulls = false;
    bool somepartials = false;
    for (int irun=0;irun<nsets;++irun) { // loop sets
      if (usetype[irun] > 0) { // some fulls
        somefulls = true;
      } else if (usetype[irun] < 0) { // partials
        somepartials = true;
      } else if (usetype[irun] == 0) { // both
        both = true;
      }
    }  // end loop run sets
    if (somefulls && somepartials) {both = true;}

    std::string label;
    if (both) {
      ss += FormatOutput::logTab(0,
         std::string(
           "                                  Fulls                        Partials\n")+
           "    Run                    SdFac    SdB   SdAdd  ISa    SdFac    SdB   SdAdd  ISa\n");
    } else {
      label = (fullpart > 0) ? "Fulls" : "Partials";
      ss += FormatOutput::logTab(0,
         std::string("                                  ")+label+"\n"+
                     "    Run                  SdFac    SdB   SdAdd  ISa\n");
    }

    for (int irun=0;irun<nsets;++irun) { // loop runs
      std::string runnum = clipper::String(runnumbers[irun], 7);
      if (allrunssame) {
        runnum = "AllRuns";
      }

      if (both) {
        double sdaf = 0.0;
        double sdbf = 0.0;
        double sdcf = 0.0;
        double sdap = 0.0;
        double sdbp = 0.0;
        double sdcp = 0.0;
        SDcorrection sdc1, sdc2;
        if (usetype[irun] >= 0) {
          sdaf = sdc_full_run[irun].SDfac();
          sdbf = sdc_full_run[irun].SDb();
          sdcf = sdc_full_run[irun].SDadd();
          sdc1  = sdc_full_run[irun];
        }
        if (usetype[irun] <= 0) {
          sdap = sdc_partial_run[irun].SDfac();
          sdbp = sdc_partial_run[irun].SDb();
          sdcp = sdc_partial_run[irun].SDadd();
          sdc2  = sdc_partial_run[irun];
        }
        label = formatUseFlag(irun);
        ss += FormatOutput::logTabPrintf(0,
                    "%s %s %7.2f %6.2f %7.4f %4.1f  %7.2f %6.2f %7.4f %4.1f\n",
                                         runnum.c_str(), label.c_str(),
                                         sdaf, sdbf, sdcf, ISa(sdc1),
                                         sdap, sdbp, sdcp, ISa(sdc2));
      } else if (fullpart > 0) {
        // Fulls
        label = formatUseFlag(irun);
        ss += FormatOutput::logTabPrintf(0,"%s %s %7.2f %6.2f %7.4f %4.1f\n",
                                         runnum.c_str(), label.c_str(),
                                         sdc_full_run[irun].SDfac(),
                                         sdc_full_run[irun].SDb(),
                                         sdc_full_run[irun].SDadd(),
                                         ISa(sdc_full_run[irun]));
      } else {
        label = formatUseFlag(irun);
        ss += FormatOutput::logTabPrintf(0,"%s %s %7.2f %6.2f %7.4f %4.1f\n",
                                         runnum.c_str(),label.c_str(),
                                         sdc_partial_run[irun].SDfac(),
                                         sdc_partial_run[irun].SDb(),
                                         sdc_partial_run[irun].SDadd(),
                                         ISa(sdc_partial_run[irun]));
      }
    }  // end loop runs

    // sample SD stuff
    if (sampleSD) {
      ss += std::string("\nFinal sigma(I) estimates will be calculated from the sample")+
      " variance of each reflection,\n"+
      "instead of from the individual SD(I), for those reflections with more than "+
        StringUtil::itos(minimumsample) + " observations\n";
    }

    return ss;
  }
  //-------------------------------------------------------------
  double SDmodel::ISa(const SDcorrection& sdc) const
  // ISa = 1/(Sdfac*SDadd)   =~ (I/sig(I))asymtotic for large I
  // see K.Diederichs, Acta Cryst. D66,733
  {
    double ISa = 0.0;
    double sdfacsdadd = sdc.SDfac()*sdc.SDadd();
    if (sdfacsdadd > 0.0) {
      ISa = 1.0/sdfacsdadd;
    }
    return ISa;
  }
  //-------------------------------------------------------------
  std::string SDmodel::asXML() const
  {
    std::string s = "<SDcorrection>\n";

    for (int irun=0;irun<nsets;++irun) { // loop runs
      if (allrunssame) {
        s += "  <AllRuns>\n";
      } else {
        s += "<Run> <number>"+ StringUtil::itos(runnumbers[irun],4)+" </number>\n";
      }
      if (usetype[irun] >= 0) {
        s += "  <Fulls>\n";
        s += "    "+StringUtil::MakeXMLtag("SDfac", sdc_full_run[irun].SDfac(),6,2);
        s += StringUtil::MakeXMLtag("SDb",   sdc_full_run[irun].SDb(),6,2);
        s += StringUtil::MakeXMLtag("SDadd", sdc_full_run[irun].SDadd(),8,4);
        s += StringUtil::MakeXMLtag("ISa", ISa(sdc_full_run[irun]),5,1);
        s += "\n  </Fulls>\n";
      }
      if (usetype[irun] <= 0) {
        s += "  <Partials>\n";
        s += "    "+StringUtil::MakeXMLtag("SDfac", sdc_partial_run[irun].SDfac(),6,2);
        s += StringUtil::MakeXMLtag("SDb",   sdc_partial_run[irun].SDb(),6,2);
        s += StringUtil::MakeXMLtag("SDadd", sdc_partial_run[irun].SDadd(),8,4);
        s += StringUtil::MakeXMLtag("ISa", ISa(sdc_partial_run[irun]),5,1);
        s += "\n  </Partials>\n";
      }
      if (allrunssame) {
        s += "  </AllRuns>\n";
      } else {
        s += "</Run>\n";
      }
    }  // end loop runs
    // Sample SD
    if (sampleSD) {
      s += "<SampleSD> ";
      s += StringUtil::MakeXMLtag("MinimumSample", minimumsample);
      s += " </SampleSD>\n";
    }
    s += "</SDcorrection>\n";

    return s;
  }
//-------------------------------------------------------------
  //! return formatted weight information
  std::string SDmodel::formatWeightType() const
  {
    return WeightType::formatWeightType(weighttype);
  }
//-------------------------------------------------------------
  std::string SDmodel::FormatSave() const
  {
    const std::string SDMODELVERSION = "V1.1";
    std::string ds = "SDModel "+SDMODELVERSION+" {\n";

    int nruns = Nruns();
    ASSERT (int(sdc_full_run.size()) == nruns);
    ASSERT (int(sdc_partial_run.size()) == nruns);
    ASSERT (int(runnumbers.size()) == nruns);
    //    ASSERT (idxfullparam.size() == nruns);
    //    ASSERT (idxpartialparam.size() == nruns);
    ASSERT (int(usetype.size()) == nruns);

    ds += "Nruns "+clipper::String(nruns)+"\n";
    for (int i=0;i<nruns;++i) { // loop runs
      ds += "SDC {\n";
      ds += "RunNumber "+clipper::String(runnumbers[i])+"\n";
      ds += "SDCfulls    " + sdc_full_run[i].FormatSave() + "\n";
      ds += "SDCpartials " + sdc_partial_run[i].FormatSave() + "\n";
      ds += "Usetype " + clipper::String(usetype[i]) + "\n";
      ds += "}\n";
    }
    ds += "NoSDB " + clipper::String(nosdb)  +"\n";
    ds += "Allrunssame " + clipper::String(allrunssame)  +"\n";
    ds += "Refine " + clipper::String(refine)  +"\n";
    ds += "Nsets " + clipper::String(nsets)  +"\n";
    ds += "Damp " + clipper::String(damp)  +"\n";

    ds += "Tietype " + clipper::String(ties.tietype)  +"\n";
    ds += "Ntargets " + clipper::String(int(ties.targets.size()))  +"\n";
    ds += "Targets " + StringUtil::FormatSaveVector(ties.targets);
    ds += "SDtargets " + StringUtil::FormatSaveVector(ties.sdtargets);

    ds += "}\n";
    return ds;  // null for now
  }
//-------------------------------------------------------------
  void SDmodel::Restore(const std::string& restorefilename,
                        const std::vector<Run>& runlist)
  {
    std::ifstream scalesin(restorefilename.c_str());
    Fileread FR(scalesin, restorefilename, "RESTORE");

    // Get run definitions from Scalemodel part of file
    FR.ReadTag("ScaleModel"); // Note ReadTag fails if tag is wrong
    if (FR.GetTag() != "V1.1") {  // version check
      clipper::Message::message(Message_fatal
                ("RESTORE incompatible version in "+restorefilename));
    }
    FR.Skip(); // skip "{"

    RunsFromSavefile savefileruns(FR, runlist);
    if (savefileruns.NumberRunsFound() < int(runlist.size())) {
      clipper::Message::message(Message_fatal
                        ("RESTORE not all runs found in save file"));
    }

    int svnruns = savefileruns.NumberRunsInSaveFile();
    std::vector<int> runsfromsavefile = savefileruns.Runsfromsavefile();
    // runsfromsavefile now contains the index in runlist
    // for each run in save file, or -1 if not wanted

    FR.SkipToTag("SDModel"); // Skip to tag
    if (FR.GetTag() != "V1.1") {  // version check
      clipper::Message::message(Message_fatal
        ("RESTORE SDmodel incompatible version in "+restorefilename));
    }
    FR.Skip(); // skip "{"

    FR.ReadTag("Nruns");
    int nr = FR.Int();
    if (nr != svnruns) {
      clipper::Message::message(Message_fatal
                        ("RESTORE SDmodel inconsisent Nruns"));
    }

    int jpr = 0; // index to accepted runs
    for (int ipr = 0;ipr<svnruns;++ipr) { // loop runs in save file
      FR.ReadTag("SDC");
      // Do we want this one?
      if (runsfromsavefile[ipr] < 0) {
        FR.SkipSection(0);  // no, skip it
      } else {
        FR.Skip(); // skip "{"
        FR.ReadTag("RunNumber");
        int runnum = FR.Int();  // this should match run number in runlist
        if (runnum != runlist[jpr].RunNumber()) {
          clipper::Message::message(Message_fatal
                    ("RESTORE SDC: mismatch run number "+
                     clipper::String(runnum)+" "+
                     clipper::String(runlist[ipr].RunNumber())));
        }
        FR.ReadTag("SDCfulls");
        sdc_full_run[jpr].Restore(FR);
        FR.ReadTag("SDCpartials");
        sdc_partial_run[jpr].Restore(FR);
        FR.ReadTag("Usetype");
        usetype[jpr] = FR.Int();
        jpr++;
        if (!FR.CheckEnd()) {
          clipper::Message::message(Message_warn
                                    ("SDmodel Restore unexpected tag "+FR.Tag()));
        }
      }  // end run
    }  // end loop runs

    FR.ReadTag("NoSDB");
    nosdb = FR.Int();
    FR.ReadTag("Allrunssame");
    allrunssame = FR.Int();
    FR.ReadTag("Refine");
    refine = FR.Int();
    FR.ReadTag("Nsets");
    nsets = FR.Int();
    FR.ReadTag("Damp");
    damp = FR.Double();
    FR.ReadTag("Tietype");
    ties.tietype = FR.Int();
    FR.ReadTag("Ntargets");
    int ntargets = FR.Int();
    FR.ReadTag("Targets");
    ties.targets = FR.DoubleVec(ntargets);
    FR.ReadTag("SDtargets");
    ties.sdtargets = FR.DoubleVec(ntargets);

    scalesin.close();

    nsets = Min(nsets, jpr); // number of accepted sets

    SetIdxParam(); // set index list
    SetTies();
  }
//-------------------------------------------------------------
}
