// statistics.cpp
//
// Print all merging statistics etc
//

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;

#include "statistics.hh"
#include "score_datatypes.hh"
#include "printing.hh"
#include "selectedobservations.hh"
#include "reject.hh"
#include "intensitybin.hh"
#include "halfdataset.hh"
#include "Output.hh"
#include "summarystatistics.hh"
#include "cumulativecompleteness.hh"
#include "sdanalysis.hh"
#include "tile.hh"
#include "anisotropy.hh"
#include "timer.hh"
#include "radiationdamageanalysis.hh"
#include "comparesds.hh"
#include "string_util.hh"

using phaser_io::LOGFILE;
using phaser_io::LXML;

namespace scala {
  // ------------------------------------------------------------
  void BatchScales0(const std::vector<Batch>& batches,
                    const Batchgroup& batchgroup,
                    const int& datasetIndex, const ScaleModel& AllScales,
                    std::vector<float>& scale0batch,
                    std::vector<float>& bfacbatch)
  // Set arrays of Primary scales at theta=0 for centre of each batch group,
  // & Bfactor, for selected dataset
  //
  // On entry:
  //  batches        list of all batches (including those not in this dataset)
  //  batchgroup     batch grouping
  //  datasetIndex   dataset index number to select dataset
  //  AllScales      scales
  //
  // On exit:
  //  scale0batch    primary scales at theta=0 for centre of each batch group
  //  bfacbatch      Bfactor for each batch
  {
    int nbatches = batches.size();
    size_t ngroups = batchgroup.numberofgroups();
    if (ngroups <= 0) {ngroups = nbatches;}
    scale0batch.assign(ngroups, 0.0);
    bfacbatch.assign(ngroups, 0.0);
    float ps;

    for (size_t i=0;i<ngroups;++i) {  // print even batches that have no reflections
      int ib = batchgroup.batchserial(i);
      // Is it this dataset?
      if ((datasetIndex < 0) || (batches[ib].datasetindex() == datasetIndex)) {
        if (batches[ib].Accepted()) {
          int irun = batches[ib].RunIndex();
          if (irun >= 0) {
            PrimaryScale pscale = AllScales.primary_scale(irun);
            if (pscale.IsBatchScale()) {
              if (pscale.ValidScale(batches[ib].num())) {
                ps = pscale.Scale(batches[ib].num());
              } else {
                ps = 0.0;
              }
            } else {
              ps = pscale.Scale(batches[ib].MidPhi());
            }
            if (ps != 0.0) {
              scale0batch[i] = 1./ps;
            } else {
              scale0batch[i] = 0.0;
            }

            RelativeBfactor bfac = AllScales.Bfactor(irun);

            if (bfac.IsBatchBfactor()) {
              int batchN = batches[ib].num();
              ps = bfac.BfactorValueB(batchN);
            } else {
              ps = bfac.BfactorValue(batches[ib].MidTime());
              //^^
              //              std::cout << "Bfacs: irun, ib, t, B " <<irun<<" "
              //                        <<ib<<" "<<batches[ib].MidTime()
              //                        <<" "<<ps<<"\n";
              //^-
            }
            bfacbatch[i] = ps;
          }
        }
      }
    } // end loop batches
  }
  // ------------------------------------------------------------
  void AddDelStats(const float& delI, const float& AvI, const int& nmult,
                   const int& jbatchgroup,  std::vector<Rfactor>& rmergebatch,
                   const int& mres, const bool& isfull,
                   std::vector<Rfactor>& rmergeRes,
                   std::vector<Rfactor>& rmergeResFull,
                   std::vector<Rfactor>& rmeasRes,
                   std::vector<Rfactor>& rpimRes,
                   const int& mint,
                   std::vector<Rfactor>& rmergeInt,
                   std::vector<Rfactor>& rmeasInt,
                   std::vector<Rfactor>& rpimInt)
  // Add in deviations to Rmerge statistics etc
  // by batch, resolution, intensity
  //
  // On entry:
  //  obs       observation
  //  delI      Ihl - <I>
  //  AvI       <I>
  //  nmult     multiplicity for this observation
  //  jbatchgroup    batch group number, < 0 don't use
  //  mres      resolution bin, < 0 don't use
  //  mint      intensity bin, < 0 don't use
  //
  // On exit:
  //  rmergebatch  updated Rmerge by batch
  //  rmergeRes    updated Rmerge by resolution
  //  rmergeResFull updated Rmerge by resolution for fulls only
  //  rmeasRes     updated Rmeas by resolution
  //  rpimRes      updated Rpim by resolution
  //  rmergeInt    updated Rmerge by intensity
  //  rmeasInt     updated Rmeas by intensity
  //  rpimInt      updated Rpim by intensity
  //
  {
    double unitw = 1.0;
    if (jbatchgroup >= 0) rmergebatch[jbatchgroup].add(delI, AvI, unitw);
    double an = nmult;
    if (mres >= 0) {
      rmergeRes[mres].add(delI, AvI, unitw); // Rmerge
      if (isfull) {
        rmergeResFull[mres].add(delI, AvI, unitw); // Rmerge for fulls
      }
      double w = sqrt(an/(an-1.0));
      rmeasRes[mres].add(delI, AvI, w);  // Rmeas
      w = sqrt(1.0/(an-1.0));
      rpimRes[mres].add(delI, AvI, w);  // Rpim
    }
    if (mint >= 0) {
      rmergeInt[mint].add(delI, AvI, unitw); // Rmerge
      double w = sqrt(an/(an-1.0));
      rmeasInt[mint].add(delI, AvI, w);  // Rmeas
      w = sqrt(1.0/(an-1.0));
      rpimInt[mint].add(delI, AvI, w);  // Rpim
    }
  }
  // ------------------------------------------------------------
  void AddDelStatsOv(const float& delI, const float& AvI, const int& nmult,
                     const int& jbatchgroup,  std::vector<Rfactor>& rmergebatchOv,
                     const int& mres, const bool& isfull, const int& irun,
                     std::vector<Rfactor>& rmergeResOv,
                     std::vector<Rfactor>& rmergeResFullOv,
                     std::vector<Rfactor>& rmeasResOv,
                     std::vector<Rfactor>& rpimResOv,
                     std::vector<std::vector<Rfactor> >& rmeasRun,
                     const int& mint,
                     std::vector<Rfactor>& rmergeIntOv,
                     std::vector<Rfactor>& rmeasIntOv,
                     std::vector<Rfactor>& rpimIntOv)
  // Add in deviations to Rmerge statistics etc
  // by resolution, over all I+ & I- observations together
  //
  // On entry:
  //  obs       observation
  //  delI      Ihl - <I>
  //  AvI       <I>
  //  nmult     multiplicity for this observation
  //  jbatchgroup    batch group number, < 0 don't use
  //  mres      resolution bin, < 0 don't use
  //  isfull    true if fully recorded
  //  mint      intensity bin, < 0 don't use
  //
  // On exit:
  //  rmergebatchOv    updated Rmerge by batch
  //  rmergeResOv      updated Rmerge by resolution
  //  rmergeResFullOv  updated Rmerge by resolution
  //  rmeasResOv       updated Rmeas by resolution
  //  rpimResOv        updated Rpim by resolution
  //  rmergeIntOv      updated Rmerge by intensity
  //  rmeasIntOv       updated Rmeas by intensity
  //  rpimIntOv        updated Rpim by intensity
  //
  {
    double unitw = 1.0;
    if (jbatchgroup >= 0) rmergebatchOv[jbatchgroup].add(delI, AvI, unitw);
    double an = nmult;
    if (mres >= 0) {
      rmergeResOv[mres].add(delI, AvI, unitw); // Rmerge
      if (isfull) {
        rmergeResFullOv[mres].add(delI, AvI, unitw); // Rmerge fulls
      }
      double w = sqrt(an/(an-1.0));
      rmeasResOv[mres].add(delI, AvI, w);  // Rmeas
      rmeasRun[irun][mres].add(delI, AvI, w);  // Rmeas
      w = sqrt(1.0/(an-1.0));
      rpimResOv[mres].add(delI, AvI, w);  // Rpim
    }
    if (mint >= 0) {
      rmergeIntOv[mint].add(delI, AvI, unitw); // Rmerge
      double w = sqrt(an/(an-1.0));
      rmeasIntOv[mint].add(delI, AvI, w);  // Rmeas
      w = sqrt(1.0/(an-1.0));
      rpimIntOv[mint].add(delI, AvI, w);  // Rpim
    }
  }
  // ------------------------------------------------------------
  void AddDelStats(const float& delI, const float& AvI,
                   const int& jbatchgroup,  std::vector<Rfactor>& rmergebatch)
  {
    double unitw = 1.0;
    if (jbatchgroup >= 0) rmergebatch[jbatchgroup].add(delI, AvI, unitw);
  }
  // ------------------------------------------------------------
  void AddChiSqBatch(SelectedObservations& selobs,
                     const Batchgroup& batchgroup,
                     std::vector<MeanSD>& meanChiSqBatch)
  // Mean ChiSq by batch
  {
    int idx;
    observation this_obs;
    std::vector<float> del2 = selobs.Delta2();
    while ((idx=selobs.next_observation(this_obs)) >= 0) {  // loop all valid observations
      int batchn = this_obs.Batch();  // batch number
      int jbatchgroup = batchgroup.batchgroup(batchn);
      if (jbatchgroup >= 0) {
        if (del2[idx] != 0) {
            meanChiSqBatch[jbatchgroup].Add(del2[idx]*del2[idx]);
        }
      }
    }
  }
  // ------------------------------------------------------------
  void AddChiSq(SelectedObservations& selobs,
                MeanSD& meanChiSq)
  // Mean ChiSq by resolution, if != 0
  {
    if (selobs.chiSq(false) != 0.0) {
      meanChiSq.Add(selobs.chiSq(false));
      //^
      //      std::cout <<"AddChiSq " <<  selobs.Reflection().hkl().format()
      //                <<"  "<<selobs.chiSq(false)<<"\n";
    }
  }
  // ------------------------------------------------------------
  void BiasSums(const SelectedObservations& allobs, const IsigI& AvI,
                MeanSD& biasRes,
                MeanSD& biasIRes)
  //
  //  Bias calculation
  //  This compairs each "partial" observation  Ihl with the mean of
  //  the "fulls" <Ifull>. For this purpose, "fulls" are considered as
  //  all observations with the minimum number of parts for this reflection
  //  (=1 for true fulls), and "partials" are all observations with more parts
  //
  // On entry:
  //  allobs    observations for this reflection (I+ & I-)
  //  AvI       average IsigI
  //
  // On exit:
  //  biasRes   update Sum(<Ifull> - Ipartial)
  //  biasIRes  update Sum(<I>)
  {
    observation this_obs;
    int idx;
    MeanSD AvIsmall;

    int minparts = 1000000;
    // Find smallest width (fulls if present)
    while ((idx=allobs.next_observation(this_obs)) >= 0) {
      minparts = Min(minparts, this_obs.num_parts_mpart());
    }
    while ((idx=allobs.next_observation(this_obs)) >= 0) {
      if (this_obs.num_parts_mpart() == minparts) {
        AvIsmall.Add(this_obs.kI());
      }
    }
    float Ifull = AvIsmall.Mean();
    //^
    //    std::cout << "minparts, Ifull " <<minparts<<" "<< Ifull <<"\n"; //^-
    while ((idx=allobs.next_observation(this_obs)) >= 0) {
      if (this_obs.num_parts_mpart() > minparts) {
        biasRes.Add(Ifull - this_obs.kI());
        //      std::cout <<"bias " << Ifull - this_obs.kI() <<"\n"; //^
        biasIRes.Add(AvI.I());
      }
    }
  }
  // ------------------------------------------------------------
  void SmoothStatisticsByBatch(const std::vector<std::vector<MeanSD> >& mnIsdResBatch,
                               const Batchgroup& batchgroup,
                               std::vector<double>& maxresbatchsmoothed,
                               std::vector<Rfactor>& rmergebatch,
                               std::vector<Rfactor>& rreferencebatchsmoothed,
                               std::vector<MeanValue>& ccreferencebatchsmoothed,
                               const std::vector<Batch>& batches,
                               const std::vector<Run>& runlist,
                               const ResoRange& ResRange,
                               const double& MinimumIoverSigmaBatch,
                               const int& NbatchSmooth)
  // Arguments:
  // mnIsdResBatch by batch group for each resolution bin = Mean(<I>/sd(<I>))
  // maxresbatchsmoothed (returned) "maximum resolution" by batch smoothed over NbatchSmooth batches
  // rmergebatch by batch group, replaced by smooth version
  // rreferencebatchsmoothed  vs. reference, replaced by smooth version
  // ccreferencebatchsmoothed vs. reference, replaced by smooth version
  // MinimumIoverSigmaBatch     threshold for resolution
  // NbatchSmooth          number of groups to smooth over
  //                       should be odd, if not will be forced to be odd here
  {
    int nbatchgroups = mnIsdResBatch.size();
    maxresbatchsmoothed.assign(nbatchgroups, 0.0);
    std::vector<Rfactor> Rsmooth(nbatchgroups);
    bool hklref = false;
    std::vector<Rfactor> Rrefsmooth;
    std::vector<MeanValue> ccrefsmooth;
    if (rreferencebatchsmoothed.size() > 0) {
      hklref = true;
      Rrefsmooth.resize(nbatchgroups);
      ccrefsmooth.resize(nbatchgroups);
    }

    int nbs = (NbatchSmooth/2)*2 + 1; // force odd

    // Loop runs
    for (size_t irun=0; irun<runlist.size(); irun++) {
      // list of groups in this run
      std::vector<int> grouplist = batchgroup.groupsinrun(irun);
      // all in same or all datasets, in batchgroup
      int nbgrouprun = batchgroup.nGroupsinRun(irun); // number of batch groups in run
      int nbsr = Min(nbgrouprun, nbs); // smoothing range for this run (nb may be even)
      if (nbsr >= 3) { // don't bother smoothing if very few batches
        nbsr = (nbsr/2)*2 + 1;        // Make it odd
        if (nbsr > nbgrouprun) nbsr -= 2; // // ... but not larger than number in run
        int half = nbsr/2;

        int jgroup0 = grouplist.at(0);  // first group in this run

        for (int k=0; k<grouplist.size(); k++) {
          int jgroup = k + jgroup0;
          // smooth from batch group i1 to i2
          int i1 = jgroup-half;
          int i2 = jgroup+half;
          if (i1 < jgroup0) {
            i1 = jgroup0;
            i2 = i1 + nbsr - 1;
          }
          if (i2 >= jgroup0+nbgrouprun) {
            i2 = jgroup0+nbgrouprun-1;
            i1 = std::max(i2-nbsr+1, jgroup0);
          }

          //      std::cout <<"Smoothing "<<i1<<" "<<i2<<" "<<nbsr<<
          //        " "<<nbatchgroups<<" "<<jgroup<<" "<<nbgrouprun<<
          //        " "<<jgroup0<<" "<<irun<<"\n";

          ASSERT (i2 < nbatchgroups);
          // Number of resolution bins
          int nrbins = mnIsdResBatch[jgroup].size();
          std::vector<MeanSD> msd(nrbins); // for each resolution bin
          for (int j=i1;j<=i2;++j) { // loop nbsr batchgroups
            Rsmooth[jgroup] += rmergebatch[j];  // Rmerge
            for (int i=0;i<nrbins;++i) { // loop resolution bins
              msd[i] += mnIsdResBatch[j][i];
            }
            // R and CC against reference, if present
            if (hklref) {
              Rrefsmooth[jgroup] += rreferencebatchsmoothed[j];
              ccrefsmooth[jgroup] += ccreferencebatchsmoothed[j];
            }
          }
          //      std::cout << i1 <<" "<<i2
          //                <<"  "<<irun
          //                <<" " << rmergebatch[jgroup].R()
          //                <<" " << Rsmooth[jgroup].R()
          //                << " i1,i2\n"; //^

          // assign resolution limit for this group to batch ib

          ResolutionLimit batchreslimit(msd, ResRange,
                                        MinimumIoverSigmaBatch,
                                        ResolutionLimit::NONE);
          maxresbatchsmoothed[jgroup] = batchreslimit.HighResolution();
        } // smoothing or not
      } // end loop batch groups
    } // end loop runs
    rmergebatch = Rsmooth; // return overwriting input
    rreferencebatchsmoothed = Rrefsmooth;
    ccreferencebatchsmoothed = ccrefsmooth;
    return;
  }
  // ------------------------------------------------------------
  std::vector<MeanValue> AverageCCoverresolution
  (const std::vector<std::vector<correl_coeff> >& ccreferencebatch,
   std::vector<int>& numberinCC)
  // average of CC over resolution bins, weighted by number of contributions,
  // also returns numberinCC containing total count in average
  {
    std::vector<MeanValue> averageccoverresolution(ccreferencebatch.size());
    if (ccreferencebatch.size() == 0) return averageccoverresolution;
    int nresbins = ccreferencebatch[0].size();
    numberinCC.assign(ccreferencebatch.size(), 0);
    for (size_t ib=0; ib<ccreferencebatch.size(); ib++) {
      for (int mres=0;mres<nresbins;++mres) {
        double weight = ccreferencebatch[ib][mres].Number();
        averageccoverresolution[ib].Add(ccreferencebatch[ib][mres].CC(), weight);
        numberinCC[ib] += ccreferencebatch[ib][mres].Number();
      }
    }
    return averageccoverresolution;
  }
  // ------------------------------------------------------------
  SummaryStatistics Statistics(const ScaleModel& AllScales,
                               const hkl_unmerge_list& hkl_list,
                               SDmodel& SDM,
                               const all_controls& controls,
                               const int& datasetIndex,
                               const ResoRange& ResRange,
                               const Normalise& NormRes,
                               const AnomDistribution& anomDistribution,
                               const float& anomProbSlope,
                               const ReferenceList& hklreflist,
                               phaser_io::Output& output)
  //
  // Statistics for within a dataset datasetIndex
  //
  //  On entry:
  //   AllScales    scale model
  //   hkl_list     list with scales applied & outliers rejected
  //   SDM          Sd correction model
  //   controls     all controls
  //   datasetIndex dataset index
  //   ResRange     resolution range with bins
  //   NormRes      normalisation object
  //   anomProbSlope slope of anomalous normal probability plot
  //   hklreflist   reference data for analysis, if present
  //   output
  //
  {
    // Check valid datasetIndex
    if (datasetIndex < 0 || datasetIndex >= hkl_list.num_datasets()) {
      Message::message(Message_fatal("Statistics: datasetIndex "+
                                     clipper::String(datasetIndex)+" out of range"));
    }

    // Summary data
    SummaryStatistics summaryStatistics;

    // Project/Crystal/Dataset for this dataset
    PxdName dataset_pxd = hkl_list.dataset(datasetIndex).pxdname();
    summaryStatistics.StorePXDname(dataset_pxd);
    std::vector<Run> runlist = hkl_list.RunList();
    Dataset this_dataset =  hkl_list.dataset(datasetIndex);
    // Statistics by batch
    //   batches in whole file, including other datasets
    int nbatches = hkl_list.num_batches();

    std::string tt = "Merging statistics for dataset "+dataset_pxd.format();
    std::string marks(tt.size()+4,'*');

    output.logTab(0,LOGFILE,"\n\n"+marks);
    output.logTab(0,LOGFILE,"* "+tt+" *");
    output.logTab(0,LOGFILE,marks+"\n");

    // Resolution ranges, overall, inner, outer
    int nresbin =  ResRange.Nbins();
    // Reset low resolution limit to real one
    //  overall
    ResoRange resrangedataset = hkl_list.dataset(datasetIndex).ResRange();
    // inner
    ResoRange resrange0 = ResRange;
    resrange0.SetRange(resrangedataset.ResLow(), ResRange.BinRange(0).ResHigh());
    summaryStatistics.StoreResRanges(resrangedataset, resrange0,
                                     ResRange.BinRange(nresbin-1));

    // Smoothing and grouping
    Batchgroup batchgroup(AllScales, controls.analysis, hkl_list, datasetIndex);
    output.logTab(0,LOGFILE,"\n"+batchgroup.format()+"\n");
    output.logTab(0,LXML,batchgroup.formatXML());
    // number of batch groups for batch analysis
    int nbatchgroups = batchgroup.numberofgroups();

    // Intensity bins etc
    int NintBin = controls.analysis.NiBins();
    //   Number of bins, number of "reference" bin,
    //   intensity at "reference" bin, maximum intensity
    float Iav = NormRes.Imean();
    float Jmax = NormRes.Imax();
    IntensityBin Irange(NintBin, NintBin/2, Iav, Jmax);

    //  these determined for each batch
    std::vector<float> scale0batch(nbatchgroups);  // mean scale at theta=0
    std::vector<float> bfacbatch(nbatchgroups);    // B-factor
    std::vector<Batch> batches = hkl_list.Batches();  // all batches
    BatchScales0(batches, batchgroup, datasetIndex, AllScales, scale0batch, bfacbatch);
    std::vector<int> rejectedbatch(nbatchgroups);   // count of outliers
    std::vector<int> rejecteddataset(hkl_list.num_datasets());   // count of outliers
    // NOT DONE //    std::vector<int> overloadsbatch(nbatchgroups);  // count of overloads
    std::vector<std::vector<MeanSD> >  mnIsdResBatch(nbatchgroups);  // Mean(<I>/sd(<I>))
    for (int i=0;i<nbatchgroups;++i) {  // ... by resolution for each batch
      mnIsdResBatch[i].assign(nresbin,MeanSD());
    }

    // Accumulated over all data
    std::vector<MeanSD>  scalebatch(nbatchgroups);     // mean scale overall
    std::vector<Rfactor> rmergebatch(nbatchgroups);    // Rmerge within I+/I-
    std::vector<Rfactor> rmergebatchOv(nbatchgroups);  // Rmerge (all I+, I-)
    std::vector<MeanSD>  imeanbatch(nbatchgroups);     // Imean (all I+, I-)
    std::vector<MeanSD>  rmsDbatch(nbatchgroups);      // RMS scatter from mean (all I+,I-)
    std::vector<int>     NumObsBatch(nbatchgroups,0);  // Number of observations
    std::vector<MeanSD>  meanChiSqBatch(nbatchgroups); // mean Chi^2
    std::vector<MeanSD>  meanChiSqBatch2(nbatchgroups);// mean Chi^2 excluding outliers

    // statistics relative to reference dataset
    bool hklref = !(hklreflist.IsEmpty());
    std::vector<Rfactor> rreferencebatch;    // Rfactor to reference data
    // CC to reference data by resolution & batch
    std::vector<std::vector<correl_coeff> > ccreferencebatch;
    std::vector<MeanValue> meanIrefbatch;
    std::vector<MeanValue> meanIobsbatch;
    if (hklref) {
      rreferencebatch.resize(nbatchgroups);
      ccreferencebatch.resize(nbatchgroups);
      for (int i=0;i<nbatchgroups;++i) {  // ... by resolution for each batch
        ccreferencebatch[i].assign(nresbin,correl_coeff());
      }
      meanIrefbatch.resize(nbatchgroups);
      meanIobsbatch.resize(nbatchgroups);
    }

    // by resolution
    // within I+/I- sets
    std::vector<Rfactor> rmergeRes(nresbin); // Rmerge
    std::vector<Rfactor> rmergeResFull(nresbin); // Rmerge for fulls
    std::vector<Rfactor> rmeasRes(nresbin);  // Rmeas
    // Rmeas by run & resolution (not all runs may be in this dataset)
    int nruns = hkl_list.num_runs();
    std::vector<std::vector<Rfactor> > rmeasRun(nruns);
    std::vector<int> nbfacrun(nruns);  // number of Bfactors in run
    for (int i=0;i<nruns;++i) {
      rmeasRun[i].assign(nresbin, Rfactor());
      nbfacrun[i] = AllScales.Bfactor(i).Number();
    }
    std::vector<Rfactor> rpimRes(nresbin);   // Rpim
    // over all I+ & I- sets
    std::vector<Rfactor> rmergeResOv(nresbin); // Rmerge
    std::vector<Rfactor> rmeasResOv(nresbin);  // Rmeas
    std::vector<Rfactor> rpimResOv(nresbin);   // Rpim
    std::vector<Rfactor> rmergeResFullOv(nresbin); // Rmerge for fulls

    std::vector<MeanSD>  imeanRes(nresbin);  // <I>
    std::vector<MeanSD>  rmsDRes(nresbin);   // RMS scatter from mean (all I+,I-)

    std::vector<MeanSD>  avSdRes(nresbin);   // Average corrected SD
    std::vector<MeanSD>  mnIsdRes(nresbin);  // Mean(<I>/sd(<I>))
    std::vector<MeanSD>  biasRes(nresbin);   // bias Mean (<I"full"> - Ihl(partial))
    std::vector<MeanSD>  biasIRes(nresbin);  // Mean <I> for fractional bias
    std::vector<int>     NumRef(nresbin,0);  // Number of unique reflections
    std::vector<int>     NumObs(nresbin,0);  // Number of observations
    std::vector<MeanSD>  meanChiSqRes(nresbin);  // mean Chi^2
    std::vector<MeanSD>  meanChiSqRes2(nresbin); // mean Chi^2 excluding outliers
    int NumRefAll = 0;
    int NumObsAll = 0;
    int NumObsFull = 0;
    int NumObsPart = 0;
    int NumObsScaled = 0;
    std::vector<int> NumRefSphere(nresbin,0);  // Number unique in sphere
    std::vector<int> NumCentric(nresbin,0);    // Number unique centric
    std::vector<int> NumACentric(nresbin,0);   // Number unique acentric
    std::vector<int> NumAnom(nresbin,0);       // number unique anomalous
    std::vector<int> NumAnomSphere(nresbin,0); // number unique in sphere
    std::vector<double> SNumAnomPairs(nresbin,0.0); // anomalous pairs
    // by intensity
    //  within I+/I- sets
    std::vector<Rfactor> rmergeInt(NintBin); // Rmerge
    std::vector<Rfactor> rmeasInt(NintBin);  // Rmeas
    std::vector<Rfactor> rpimInt(NintBin);   // Rpim
    // over all I+ & I- sets
    std::vector<Rfactor> rmergeIntOv(NintBin); // Rmerge
    std::vector<Rfactor> rmeasIntOv(NintBin);  // Rmeas
    std::vector<Rfactor> rpimIntOv(NintBin);   // Rpim
    std::vector<MeanSD>  imeanInt(NintBin);  // <I>
    std::vector<MeanSD>  rmsDInt(NintBin);   // RMS scatter from mean (all I+,I-)
    std::vector<MeanSD>  avSdInt(NintBin);   // Average corrected SD
    std::vector<MeanSD>  mnIsdInt(NintBin);  // Mean(<I>/sd(<I>))
    std::vector<MeanSD>  biasInt(NintBin);   // bias Mean (<I"full"> - Ihl(partial))
    std::vector<MeanSD>  biasIInt(NintBin);  // Mean <I> for fractional bias

    // ---- Set up anisotropy directions (axes or planes)
    // I/sd analysis in up to 3 anisotropic directions, by resolution
    std::vector<std::vector<MeanSD> > mnIsdResAniso(3);  // Mean(<I>/sd(<I>))
    for (int i=0;i<3;++i) {mnIsdResAniso[i].resize(nresbin);}

    Timer anisotime;
    // Get principal axes of anisotropy depending on symmetry and data
    AnisotropicAnalysis anisoanal(hkl_list, datasetIndex, SDM);
    anisoanal.SetConeAngle(controls.analysis.ConeAngle());  // store cone angle
    summaryStatistics.StoreAnisoDeltaB(anisoanal.BfactorDifference());
    output.logTab(0, LOGFILE,
         "\nTime for determination of anisotropic axes: "+anisotime.format(true));
    // ----

    // Sample SD option
    CompareSDs comparesds;
    if (SDM.SampleSD()) {
      comparesds.init(ResRange, Irange, SDM.MinimumSample(), batchgroup);
    }

    // Half dataset correlations etc, by resolution
    HalfDataset halfDatasetScores(nresbin, dataset_pxd);
    // Store relevant anomalous statistics
    halfDatasetScores.StoreAnomStats(anomDistribution);

    // Number of resolution ranges for cumulative analysis
    int nCumulativeResoBins = -1;  // automatic setting
    CumulativeCChalf cumulativecchalf(ResRange, nCumulativeResoBins, batchgroup);

    // SD analysis by intensity, runs, full/partial
    //  Separate run analysis even if SDM.AllRunsSame()
    SDanalysis sdanalysis(Irange, SDM, false, false);
    //^
    //    sdanalysis.SetDump("sddump.dat");  //^ dumping analysis data
    // "core" data only, ie within smaller limits on delta
    SDanalysis sdanalysiscore(Irange, SDM, false, false);
    //^    sdanalysiscore.SetDump("sdanal.dat"); //^^-
    // Analysis on detector
    DetectorAnalysis detectoranalysis;
    if (controls.analysis.DetectorAnalysis()) {
      detectoranalysis.init(hkl_list);
    }

    reflection this_refl;
    observation this_obs;

    CumulativeCompleteness cumulativecompleteness(nbatchgroups);

    bool Anom = controls.anomalouscontrol.Anomalous;
    summaryStatistics.SetAnom(Anom);

    SelectedObservations allobs;     // all I+ and I-
    SelectedObservations obsplus;    // just I+
    SelectedObservations obsminus;   // just I-
    std::vector<float> delI;
    std::vector<float> delIplus;
    std::vector<float> delIminus;
    std::vector<IvarI> AvIothers;
    IsigI AvIsig, AvIsigplus, AvIsigminus;
    SDM.ResetRange();  // range of sd correction values

    // For SD analysis
    float sdrej = 5.0;     // for now, FIXME
    float sdrej2 = sdrej;
    scala::RejectFlags::Reject2Policy Rej2policy = scala::RejectFlags::KEEP;
    double multiplicity = double(hkl_list.num_observations())/
      double(hkl_list.num_reflections_valid());
    const double MINMULTFORKEEP2 = 1.5;
    if (multiplicity < MINMULTFORKEEP2) {
        output.logTab(0,LOGFILE,
                      std::string("WARNING: multiplicity is low, ")+
                      StringUtil::Strip(StringUtil::ftos(multiplicity,8,1))+
                      " (below threshold "+
                      StringUtil::Strip(StringUtil::ftos(MINMULTFORKEEP2,8,1))+
                      ") so deviant reflections measured twice are KEPT in SD and Chi^2 analysis");
    } else {
        Rej2policy = scala::RejectFlags::REJECT;
        output.logTab(0,LOGFILE,
                      std::string("Multiplicity ")+
                      StringUtil::Strip(StringUtil::ftos(multiplicity,8,1))+
                      " is above threshold "+
                      StringUtil::Strip(StringUtil::ftos(MINMULTFORKEEP2,8,1))+
                      " so deviant reflections measured twice are REJECTED in SD and Chi^2 analysis");
    }
    RejectFlags rejflags(sdrej, sdrej2, Rej2policy);

    // Count outliers/batch
    std::vector<int> outliercount = CountOutliers(hkl_list, rejectedbatch, rejecteddataset);
    double maxinvresolsq = 0.0; // actual maximum resolution
    // number of symmetry operators including lattice centering (since epsilon allows
    // for lattice centering)
    float NumSymm = hkl_list.symmetry().Nsym();
    hkl_list.rewind();

    // * * * * Loop reflections
    while (hkl_list.next_reflection(this_refl) >= 0)  {
      Rtype invresolsq = this_refl.invresolsq();
      bool Centric = hkl_list.symmetry().is_centric(this_refl.hkl());

      // Multiplicity is number of times this reflection will occur in a
      //     complete sphere of data. This is Nsym/Epsilon, multiplied by 2
      //     for acentric reflections (since Nsym symmetry operations generate
      //     a hemisphere of acentric data, but all centric reflections)
      // NOTE this is not the same as the multiplicity used for weighting
      //     <I> in Wilson plot calculations (see Iwasaki & Ito, Acta Cryst,
      //     A33,227-229(1977))
      float epsiln = hkl_list.symmetry().epsilon(this_refl.hkl());
      int multcy;
      if (epsiln == 0.0) {
        // systematic absence
        multcy = 0;
      } else {
        multcy = Nint(float(NumSymm)/epsiln);
        if (!Centric) multcy = multcy*2;
      }

      // Resolution bin
      int mres = ResRange.bin(invresolsq);
      maxinvresolsq = Max(maxinvresolsq, invresolsq);

      //  Apply current SD correction to reflection (all observations)
      SDM.CorrectReflection(this_refl);

      // Select all (I+ & I-) accepted observations for this dataset
      allobs.init(this_refl, datasetIndex, ALL);
      AvIsig = allobs.Average();  // average I, 1/variance weight
      delI = allobs.DelI();
      AvIothers = allobs.MeanIothers();

      // Intensity bins
      int mint = Irange.bin(AvIsig.I());

      // Anisotropic analysis, get index for direction (2 or 3 directions), = -1 if none
      std::pair<int,double> axisw = anisoanal.Axis(this_refl.hkl(), invresolsq);
      int jconeaxis = axisw.first;        // axis
      double wtaniso = axisw.second;  // weight for CCs

      // Counts
      if (allobs.Number() > 0) {
        float IovsigI = AvIsig.I()/AvIsig.sigI();
        mnIsdRes[mres].Add(IovsigI);
        mnIsdInt[mint].Add(IovsigI);
        // by cone
        if (mres == 0) {
          for (int j=0;j<3;++j) {  // lowest res bin, add into all directions
            mnIsdResAniso[j][mres].Add(IovsigI);
          }
        } else if (jconeaxis >= 0) {
          mnIsdResAniso[jconeaxis][mres].Add(IovsigI);
        }

        NumRef[mres]++;                    // Number unique
        NumObs[mres] += allobs.Number();   // Number observed
        NumRefAll++;
        NumObsAll += allobs.Number();

        // Counts for completeness & multiplicity
        //  Total in sphere allowing for symmetry multiplicity
        NumRefSphere[mres] += multcy;      // Total unique in sphere
        if (Centric) {
          NumCentric[mres]++;
        } else {
          NumACentric[mres]++;
        }

      }
      if (allobs.Number() > 1) {
        BiasSums(allobs, AvIsig, biasRes[mres], biasIRes[mres]); // Bias
        BiasSums(allobs, AvIsig, biasInt[mint], biasIInt[mint]); // Bias
      }
      // Store multiplicity of reflection for cumulative completeness
      cumulativecompleteness.StartReflection(multcy);
      int idx;

      while ((idx=allobs.next_observation(this_obs)) >= 0) {  // loop all valid observations
        int batchn = this_obs.Batch();  // batch number
        int jbatchgroup = batchgroup.batchgroup(batchn);
        bool isfull = (this_obs.PartFlag() == FULL);  // true if fully recorded, false for partial
        NumObsBatch[jbatchgroup] += multcy;              // Number observed (in sphere)

        if (isfull) {
          NumObsFull++;
        } else {
          NumObsPart++;
          if (this_obs.PartFlag() == SCALE) {NumObsScaled++;}  // scaled partial
        }

        // record an observation for cumulative completeness
        cumulativecompleteness.AddObservationBatch(batchn, jbatchgroup);
        scalebatch[jbatchgroup].Add(1./this_obs.Gscale());    // actual scale = 1/g

        imeanbatch[jbatchgroup].Add(this_obs.kI());        // Imean (all I+-)
        imeanRes[mres].Add(this_obs.kI());
        avSdRes[mres].Add(this_obs.ksigI());
        imeanInt[mint].Add(this_obs.kI());
        avSdInt[mint].Add(this_obs.ksigI());
        // by resolution for each batch
        mnIsdResBatch[jbatchgroup][mres].Add(this_obs.kI()/this_obs.ksigI());
        // Analysis against reference data
        if (hklref) {
          // find matching hkl, if present

          IsigI Isref = hklreflist.Isig(this_refl.hkl());
          if (Isref.sigI() > 0.0) {
            rreferencebatch[jbatchgroup].add(this_obs.kI()-Isref.I(), this_obs.kI(), 1.0);
            ccreferencebatch[jbatchgroup][mres].add(this_obs.kI(), Isref.I(), 1.0);
            meanIrefbatch[jbatchgroup].Add(Isref.I());
            meanIobsbatch[jbatchgroup].Add(this_obs.kI());
          }
        }
        if (allobs.Number() > 1) {
          rmsDbatch[jbatchgroup].Add(delI[idx]*delI[idx]);  // Sum(DelI^2) (all I+-)
          // Rmerge(batch) all I+ and I-
          AddDelStats(delI[idx], AvIsig.I(), jbatchgroup, rmergebatch);
          rmsDRes[mres].Add(delI[idx]*delI[idx]);  // Sum(DelI^2) (all I+-)
          rmsDInt[mint].Add(delI[idx]*delI[idx]);  // Sum(DelI^2) (all I+-)
          // Detector analysis
          if (controls.analysis.DetectorAnalysis()) {
            // AvIothers   <I> of other observations
            float xd = this_obs.XYdet().first;
            float yd = this_obs.XYdet().second;
            detectoranalysis.AddStats(this_obs.kI(), AvIothers[idx].I(),
                                      this_obs.run(), xd, yd);
          }
          if (Centric) {
            // No anomalous
            // Rmerge etc  (added in later for acentrics)
            AddDelStats(delI[idx], AvIsig.I(), allobs.Number(),
                        jbatchgroup, rmergebatch,
                        mres, isfull,
                        rmergeRes, rmergeResFull, rmeasRes, rpimRes,
                        mint, rmergeInt, rmeasInt, rpimInt);
          }
          // over all I+ & I- sets
          AddDelStatsOv(delI[idx], AvIsig.I(), allobs.Number(),
                        jbatchgroup, rmergebatchOv,
                        mres, isfull, this_obs.run(),
                        rmergeResOv, rmergeResFullOv, rmeasResOv, rpimResOv, rmeasRun,
                        mint, rmergeIntOv, rmeasIntOv, rpimIntOv);
        }
      }  // end loop observations


      // Correlations on <I>
      halfDatasetScores.AddMean(mres, allobs);
      // Anisotropic analysis
      // Halfdataset correlations by cone etc: for inner resolution bin, use all data
      if (jconeaxis >= 0 || mres == 0) {
        halfDatasetScores.AddAniso(mres, jconeaxis, wtaniso, allobs);
      }
      cumulativecchalf.addreflection(mres, allobs);

      // ---- For anisotropic analysis on projections, expand symmetry
      int nsymp = hkl_list.symmetry().NsymP(); // number of primitive operations
      //float normscale = NormRes.Corr(invresolsq); // Normalisation factor (multiplying)
      float normscale = 1.0;

      for (int isym=1;isym<=nsymp*2;isym+=2) { // loop odd ISYM, for I+
        DVect3 projection =
          anisoanal.Projection(hkl_list.symmetry().get_from_asu(this_refl.hkl(), isym), true);
        IVect3 anisores;
        for (int i=0;i<3;++i) {
          projection[i] *= projection[i];  // square component
          anisores[i] = ResRange.bin(projection[i]);
        }
        halfDatasetScores.AddAnisoProjection(anisores, allobs, normscale);
      }
      // ----

      if (Centric) {
        // Dummy anomalous as control
        halfDatasetScores.AddAnomCentric(mres, allobs);
      } else {
        // Always get I+ & I- sets (unless centric) for anomalous analysis
        obsplus.init(this_refl, datasetIndex, IPLUS);
        obsminus.init(this_refl, datasetIndex, IMINUS);
        halfDatasetScores.AddAnom(mres, obsplus, obsminus);

        // Statistics within  the I+/I-1 sets even if Anomalous Off
        AvIsigplus = obsplus.Average();  // average I+, 1/variance weight
        delIplus = obsplus.DelI();
        AvIsigminus = obsminus.Average();  // average I+, 1/variance weight
        delIminus = obsminus.DelI();

        // Counts for anomalous completeness & multiplicity
        bool both = false;
        if (obsplus.Number() > 0 && obsminus.Number() > 0) {
          both = true;
          NumAnom[mres]++;           // number unique
          NumAnomSphere[mres] += multcy;   // number unique in sphere
          // Multiplicity = Min(n+, n-) + Dn/(Dn+1) where Dn = ||n+ - n-||
          float Dn = std::abs(obsplus.Number() - obsminus.Number());
          SNumAnomPairs[mres] += Min(obsplus.Number(), obsminus.Number()) + Dn/(Dn+1.0f);
        }

        while ((idx=obsplus.next_observation(this_obs)) >= 0) {  // loop I+ observations
          int batchn = this_obs.Batch();  // batch number
          int jbatchgroup = batchgroup.batchgroup(batchn);
          // record an observation for cumulative completeness
          if (both) cumulativecompleteness.AddObservationBatch(batchn, jbatchgroup, IPLUS);
          if (obsplus.Number() > 1) {
            // Rmerge etc
            bool isfull = (this_obs.PartFlag() == FULL);
            ///     jbatchgroup = -1;  // rmergebatch already done for all I+ and I-, switch off here
            AddDelStats(delIplus[idx], AvIsigplus.I(), obsplus.Number(),
                        jbatchgroup, rmergebatch,
                        mres, isfull,
                        rmergeRes, rmergeResFull, rmeasRes, rpimRes,
                        mint, rmergeInt, rmeasInt, rpimInt);
          }
        }

        while ((idx=obsminus.next_observation(this_obs)) >= 0) {  // loop I- observations
          int batchn = this_obs.Batch();  // batch number
          int jbatchgroup = batchgroup.batchgroup(batchn);
          // record an observation for cumulative completeness
          if (both) cumulativecompleteness.AddObservationBatch(batchn, jbatchgroup, IMINUS);
          if (obsminus.Number() > 1) {
            // Rmerge etc
            bool isfull = (this_obs.PartFlag() == FULL);
            ///     jbatchgroup = -1;  // rmergebatch already done for all I+ and I-, switch off here
            AddDelStats(delIminus[idx], AvIsigminus.I(), obsminus.Number(),
                        jbatchgroup, rmergebatch,
                        mres, isfull,
                        rmergeRes, rmergeResFull, rmeasRes, rpimRes,
                        mint, rmergeInt, rmeasInt, rpimInt);
          }
        }
      } // end acentric

      // Analysis of delta = deviation/sigma in intensity bins
      // intensity bins mint / NintBin
      // runs                  nruns
      // full/partial
      if (!Anom || Centric) {
        // No anomalous, selectedobservations are in allobs
        sdanalysis.AddSelobsDelta2(allobs, mint);
        AddChiSq(allobs, meanChiSqRes[mres]);
        AddChiSqBatch(allobs, batchgroup, meanChiSqBatch);
        allobs.Outliers(rejflags);
        sdanalysiscore.AddSelobsDelta2(allobs, mint);
        AddChiSq(allobs, meanChiSqRes2[mres]);
        AddChiSqBatch(allobs, batchgroup, meanChiSqBatch2);
        if (SDM.SampleSD()) {
          comparesds.add(allobs, mres, mint);
        }
      } else {
        // Anomalous
        sdanalysis.AddSelobsDelta2(obsplus, mint);
        sdanalysis.AddSelobsDelta2(obsminus, mint);
        AddChiSq(obsplus, meanChiSqRes[mres]); // if not 0.0
        AddChiSq(obsminus, meanChiSqRes[mres]);
        AddChiSqBatch(obsplus, batchgroup, meanChiSqBatch);
        AddChiSqBatch(obsminus, batchgroup, meanChiSqBatch);
        obsplus.Outliers(rejflags);
        obsminus.Outliers(rejflags);
        sdanalysiscore.AddSelobsDelta2(obsplus, mint);
        sdanalysiscore.AddSelobsDelta2(obsminus, mint);
        AddChiSq(obsplus, meanChiSqRes2[mres]); // if not 0.0
        AddChiSq(obsminus, meanChiSqRes2[mres]);
        AddChiSqBatch(obsplus, batchgroup, meanChiSqBatch2);
        AddChiSqBatch(obsminus, batchgroup, meanChiSqBatch2);
        if (SDM.SampleSD()) {
          comparesds.add(obsplus, mres, mint);
          comparesds.add(obsminus, mres, mint);
        }
      }
      cumulativecompleteness.EndReflection();
    }  // end loop reflections

    // ================================================================

    std::vector<float> batchcompleteness =
      cumulativecompleteness.BatchCompleteness
      (ResRange, hkl_list.symmetry(), hkl_list.Cell());
    std::vector<float> batchanomcompleteness =
      cumulativecompleteness.BatchAnomCompleteness
      (ResRange, hkl_list.symmetry(), hkl_list.Cell());
    std::vector<float> batchmultiplicity =
      cumulativecompleteness.BatchMultiplicity
      (NumObsBatch, ResRange, hkl_list.symmetry(), hkl_list.Cell());

    // Actual maximum invresolsq
    summaryStatistics.StoreMaxinvresolsq(maxinvresolsq);

    // Estimates of "maximum resolution" for each batch, based on MinimumIoverSigma
    std::vector<double> maxresbatch(nbatchgroups);
    double MinimumIoverSigmaBatch = controls.analysis.MinimumBatchIoverSigma();
    for (int i=0;i<nbatchgroups;++i) {  // ... by resolution for each batch
      int ib = batchgroup.batchserial(i);
      if (batches[ib].datasetindex() == datasetIndex) {
        //std::cout <<"Entering Batchlimit " << ib<<"\n";
        ResolutionLimit batchreslimit(mnIsdResBatch[i], ResRange,
                                      MinimumIoverSigmaBatch,
                                      ResolutionLimit::NONE);
        maxresbatch[i] = batchreslimit.HighResolution();
        //std::cout <<"Batchlimit " << ib<<" "<<maxresbatch[ib]<<"\n";
      } else {
        maxresbatch[i] = 0.0;
      }
    }
    // and generate a smoothed version of this, as well as Rmerge
    std::vector<double> maxresbatchsmoothed = maxresbatch;
    std::vector<Rfactor> rmergebatchsmoothed = rmergebatch;

    std::vector<Rfactor> rreferencebatchsmoothed;    // Rfactor to reference data

    std::vector<int> numberinCC; // number in CC for each batch
    std::vector<MeanValue> averageccbatch = AverageCCoverresolution(ccreferencebatch,
                                                                    numberinCC);
    std::vector<MeanValue> averageccbatchsmoothed;    // CC to reference data
    if (hklref) {
      // copy, to be overwritten
      rreferencebatchsmoothed = rreferencebatch;    // Rfactor to reference data
      averageccbatchsmoothed  = averageccbatch;     // CC to reference data
    }

    if (batchgroup.numGroupSmooth() > 1) {
      SmoothStatisticsByBatch(mnIsdResBatch, batchgroup, maxresbatchsmoothed,
                              rmergebatchsmoothed,
                              rreferencebatchsmoothed,
                              averageccbatchsmoothed,
                              batches,
                              runlist,
                              ResRange,
                              MinimumIoverSigmaBatch,
                              batchgroup.numGroupSmooth());
    }

    // Print stuff
    output.logTabPrintf(0,LOGFILE,
                "\nAccepted data:\nNumber of unique reflections                  %9d\n",
                        NumRefAll);
    output.logTabPrintf(0,LOGFILE,"Number of observations                        %9d\n",
                        NumObsAll);
    if (NumObsPart > 0 && NumObsFull > 0) {
      output.logTabPrintf(0,LOGFILE,"Number of fully-recorded observations         %9d\n",
                          NumObsFull);
      output.logTabPrintf(0,LOGFILE,"Number of partially-recorded observations     %9d\n",
                          NumObsPart);
      output.logTabPrintf(0,LOGFILE,"Number of scaled partial observations         %9d\n",
                          NumObsScaled);
    }
    if (hkl_list.num_datasets() > 1) {
      output.logTabPrintf(0,LOGFILE,"\nNumber of rejected outliers (this dataset)    %9d\n",
                          rejecteddataset.at(datasetIndex));
    } else {
      output.logTabPrintf(0,LOGFILE,"\nNumber of rejected outliers                   %9d\n",
                          outliercount.at(0)+outliercount.at(1));
    }

    output.logTabPrintf(0,LOGFILE,"Number of observations rejected on Emax limit %9d\n\n",
                        outliercount.at(2));

    PrintScalesByBatch(dataset_pxd, batches, batchgroup, runlist, datasetIndex,
                       scale0batch, bfacbatch, nbfacrun, scalebatch,
                       output);
    PrintDeviationsByBatch(dataset_pxd, batches, batchgroup, datasetIndex,
                           imeanbatch, rmsDbatch, rmergebatch, rmergebatchsmoothed,
                           rejectedbatch,
                           batchcompleteness, batchanomcompleteness, batchmultiplicity,
                           maxresbatch, maxresbatchsmoothed,
                           meanChiSqBatch, meanChiSqBatch2,
                           MinimumIoverSigmaBatch, batchgroup.numGroupSmooth(),
                           ResRange,rejflags, output);
    // Analysis by batch against reference
    if (hklref) {
      PrintComparisonToReferenceByBatch(dataset_pxd, batches, batchgroup, datasetIndex, batchgroup.numGroupSmooth(),
                                        rreferencebatch, averageccbatch, numberinCC,
                                        rreferencebatchsmoothed, averageccbatchsmoothed,
                                        meanIrefbatch, meanIobsbatch,
                                        output);
    }

    // process halfdataset scores, work out resolution "limits"
    halfDatasetScores.Analyse(ResRange,
                              controls.analysis.MinimumHalfdatasetCC(),
                              controls.analysis.MinimumHalfdatasetAnomCC());

    PrintHalfDatasetCorrelations(dataset_pxd,
                                 ResRange,
                                 halfDatasetScores,
                                 summaryStatistics,
                                 output);

    PrintAnisotropyAnalysis(dataset_pxd,
                            ResRange, halfDatasetScores, mnIsdResAniso,
                            anisoanal, controls.analysis.MinimumIoverSigma(),
                            summaryStatistics, output);

    cumulativecchalf.printcumulativeCC(datasetIndex, runlist, output);

    // Statistics within I+/I- sets
    //  rmergeRes, rmeasRes, rpimRes, rmergeResFull
    // or overall I+-
    //  rmergeResOv, rmeasResOv, rpimResOv, rmergeResFullOv
    if (Anom) { // anomalous on, print statistics with I+/I- sets
      PrintDeviationsByResolution(dataset_pxd, ResRange, Anom,
                                  rmergeRes, rmergeResFull, rmeasRes,
                                  rpimRes, imeanRes, rmsDRes, avSdRes, mnIsdRes,
                                  biasRes, biasIRes, meanChiSqRes, meanChiSqRes2,
                                  controls.analysis.MinimumIoverSigma(),
                                  summaryStatistics, rejflags, output);
    } else { // no anomalous, use overall statistics
      PrintDeviationsByResolution(dataset_pxd, ResRange, Anom,
                                  rmergeResOv, rmergeResFullOv, rmeasResOv,
                                  rpimResOv, imeanRes, rmsDRes, avSdRes, mnIsdRes,
                                  biasRes, biasIRes, meanChiSqRes, meanChiSqRes2,
                                  controls.analysis.MinimumIoverSigma(),
                                  summaryStatistics, rejflags, output);
    }

    PrintDeviationsByResolutionOv(dataset_pxd, ResRange,
                                  rmergeRes, rmeasRes, rpimRes,
                                  rmergeResOv, rmeasResOv, rpimResOv,
                                  summaryStatistics, output);


    // Within I+/I- sets      rmergeInt, rmeasInt, rpimInt
    // overall I+/I-          rmergeIntOv, rmeasIntOv, rpimIntOv
    if (Anom) { // anomalous on, print statistics with I+/I- sets
      PrintDeviationsByIntensity(dataset_pxd, Irange, Anom, rmergeInt, rmeasInt,
                               rpimInt, imeanInt, rmsDInt, avSdInt, mnIsdInt,
                               biasInt, biasIInt, output);
      summaryStatistics.StoreRtopI(rmergeInt[NintBin-1]);
    } else { // no anomalous
      PrintDeviationsByIntensity(dataset_pxd, Irange, Anom, rmergeIntOv, rmeasIntOv,
                               rpimIntOv, imeanInt, rmsDInt, avSdInt, mnIsdInt,
                               biasInt, biasIInt, output);
      summaryStatistics.StoreRtopI(rmergeIntOv[NintBin-1]);
    }

    PrintDeviationsByRun(dataset_pxd, ResRange, hkl_list.RunList(), rmeasRun, output);

    PrintCompletenessMultiplicity(dataset_pxd, ResRange,
                                  hkl_list.symmetry(), hkl_list.Cell(),
                                  NumRef, NumObs, NumRefSphere, NumCentric, NumACentric,
                                  NumAnom, NumAnomSphere, SNumAnomPairs,
                                  summaryStatistics, output);


    PrintSDanalysis(sdanalysis, sdanalysiscore,
                    rejflags, Irange, runlist,
                    SDM, datasetIndex, dataset_pxd, true, output);


    if (SDM.SampleSD()) {
      comparesds.printByResolution(output);
      comparesds.printByIntensity(Irange, output);
      comparesds.printByBatch(datasetIndex, output);
    }

    // Correlplot
    std::string s = halfDatasetScores.PlotCorrel();
    if (s != "") {
      output.logTab(0,LXML,s);
    }

    // Analysis on the detector
    if (controls.analysis.DetectorAnalysis()) {
      detectoranalysis.WriteImages("DETECTORIMAGE");  // write out analyses as images
    }

    // Radiation damage analysis
    //  only if one run, and not Batch scaling
    if (runlist.size() == 1 && !AllScales.isAllBatch()
        && nbatchgroups > 1) {
      RadiationDamageAnalysis radiationdamageanalysis(hkl_list, 0,
                                                      batchgroup);
      radiationdamageanalysis.plot(batchcompleteness, output);
    }

    // Other things for summary
    Scell avcell = hkl_list.cell(dataset_pxd);
    summaryStatistics.StoreAverageCell(avcell);
    //    summaryStatistics.StoreAverageCell(hkl_list.cell(dataset_pxd));
    summaryStatistics.StoreSpaceGroupName(hkl_list.symmetry().symbol_xHM());
    float minsdcorrfulls, maxsdcorrfulls, minsdcorrpartials, maxsdcorrpartials;
    SDM.GetSDcorrectionRanges(minsdcorrfulls, maxsdcorrfulls,
                              minsdcorrpartials, maxsdcorrpartials);
    summaryStatistics.StoreSDcorrectioRange(minsdcorrfulls, maxsdcorrfulls,
                                            minsdcorrpartials, maxsdcorrpartials);
    summaryStatistics.StoreAnomNPslope(anomProbSlope);
    summaryStatistics.StoreAverageMosaicity(hkl_list.dataset(datasetIndex).Mosaicity());
    summaryStatistics.StoreAnisoAxisLabels(anisoanal.Axesformat());
    summaryStatistics.StoreNlattices(hkl_list.NumberofLattices());
    return summaryStatistics;
  }  // Statistics
  // ------------------------------------------------------------
}
