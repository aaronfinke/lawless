// halfdataset.hh

#include "halfdataset.hh"
#include "file_util.hh"
#include "plotfiles.hh"
#include "anomdistribution.hh"
#include "tablegraph.hh"
#include "string_util.hh"

#define ASSERT assert

namespace scala {
  // ------------------------------------------------------------
  HalfDataset::HalfDataset(const int& NresBins, const PxdName& Dataset_pxd)
    : iscorrelplot(false)
  {
    init(NresBins, Dataset_pxd);
  }
  // ------------------------------------------------------------
  void HalfDataset::init(const int& NresBins, const PxdName& Dataset_pxd)
  {
    dataset_pxd = Dataset_pxd;
    nresbin = NresBins;
    ccanomreso.resize(nresbin);
    ccanomresoCen.resize(nresbin);
    ccIreso.resize(nresbin);
    meanIreso.resize(nresbin);
    meanVarianceImeanreso.resize(nresbin);
    rsplitreso.resize(nresbin);
    rmsCorrel.resize(nresbin);
    rmsError.resize(nresbin);
    rmsCorrelCen.resize(nresbin);
    rmsErrorCen.resize(nresbin);
    maxDelAnom = 5.0;
    ccaniso.resize(3);
    for (int i=0;i<3;++i) {ccaniso[i].resize(nresbin);}
    // Storage for projections, divided into the same bins as resolution,
    // but along the principle axes
    ccanisoprj.resize(3);
    for (int i=0;i<3;++i) {ccanisoprj[i].resize(nresbin);}
    anisoresolimit.resize(3);
  }
  // ------------------------------------------------------------
  void HalfDataset::StoreAnomStats(const AnomDistribution& anomDistribution)
  // Store relevant anomalous statistics
  // Store RMS DelAnom for this dataset, for each resolution bin
  //  & and number of reflections with 2 delAnom
  {
    std::vector<MeanSD> RMSdelanom = anomDistribution.RMSdelAnom();
    std::vector<int> nAnomPairsRes = anomDistribution.NdelAnomHalf();
    ASSERT (nresbin <= int(RMSdelanom.size()) && nresbin <= int(nAnomPairsRes.size()));
    rmsdelanom.resize(nresbin);
    nAnomPairs = 0;

    MeanSD rmsOverall;
    for (int i=0;i<nresbin;++i) {
      rmsdelanom[i] = RMSdelanom[i].SD();
      rmsOverall += RMSdelanom[i];
      nAnomPairs += nAnomPairsRes[i];
    }
    rmsdelanomOverall = rmsOverall.SD();
    correlplot.init("DelAnom/RMS scatter plot",
                    dataset_pxd.format(), nresbin,
                    rmsdelanomOverall, nAnomPairs);
    iscorrelplot = true;
  }
  // ------------------------------------------------------------
  void HalfDataset::StoreRMS(std::vector<MeanSD>& RMSdelanom)
  {
    ASSERT (int(RMSdelanom.size()) == nresbin);
    rmsdelanom.resize(nresbin);
    MeanSD rmsOverall;
    for (int i=0;i<nresbin;++i) {
      rmsdelanom[i] = RMSdelanom[i].SD();
      rmsOverall += RMSdelanom[i];
    }
  }
  // ------------------------------------------------------------
  void HalfDataset::AddMean(const int& mres, SelectedObservations& allobs)
  // Add into sums, for Imean
  {
    float I1, I2;
    allobs.SetNpart(2);
    if (allobs.HalfAverages(I1, I2)) {
      ccIreso[mres].add(I1, I2);
      rsplitreso[mres].add((I1-I2), 0.5*(I1+I2), 1.0);
    }

    // Calculate CC(1/2) from variances as in
    //  Assam, Brehm & Diederichs, (2016), J.Appl.Cryst. 49, 1021-1028
    // For each resolution range we want:
    //  1) Var(<I>) = sigma(y)^2 (related to Var(Jtrue) - error)
    //  2) Mean Var(<Iobs>) = sigma(eps_all)^2 = 0.5 * sigma(eps_half)^2
    double Iav = allobs.Average().I(); // also get sample variance of mean
    meanIreso[mres].Add(Iav);
    double varI = allobs.Variance(); // variance of mean
    if (varI > 0.0) {
      meanVarianceImeanreso[mres].Add(varI);
    }
  }
  // ------------------------------------------------------------
  void HalfDataset::AddAnomCentric(const int& mres,
                                   SelectedObservations& allobs)
  // Add into sums, centric without anomalous
  {
    const float cossin45 = 0.707106781;   // cos 45 = sin 45
    std::vector<float> Is;
    allobs.SetNpart(4);
    if (allobs.PartAverages(Is)) {
      float del1 = Is[0] - Is[2];  // two estimates of DelAnom
      float del2 = Is[1] - Is[3];
      // Exclude unfeasibly large differences
      if (del1 < maxDelAnom*rmsdelanom[mres] && del2 < maxDelAnom*rmsdelanom[mres]) {
        ccanomresoCen[mres].add(del1, del2);
        // Sums for RMS correlation ratio
        rmsCorrelCen[mres].Add(cossin45*(del1+del2));
        rmsErrorCen[mres].Add(cossin45*(del1-del2));      }
    }
  }
  // ------------------------------------------------------------
  void HalfDataset::AddAnom(const int& mres,
                            SelectedObservations& obsplus,
                            SelectedObservations& obsminus)
  // Add into sums, acentric with anomalous
  {
    const float cossin45 = 0.707106781;   // cos 45 = sin 45
    float I1p, I2p, I1m, I2m;
    if ((obsplus.HalfAverages(I1p, I2p)) && (obsminus.HalfAverages(I1m, I2m))) {
      float del1 = I1p - I1m;  // two estimates of DelAnom
      float del2 = I2p - I2m;
      // Exclude unfeasibly large differences
      if (std::abs(del1) < maxDelAnom*rmsdelanom[mres] &&
          std::abs(del2) < maxDelAnom*rmsdelanom[mres]) {
        ccanomreso[mres].add(del1, del2);
        // Sums for RMS correlation ratio
        rmsCorrel[mres].Add(cossin45*(del1+del2));
        rmsError[mres].Add(cossin45*(del1-del2));
        // Add point for correlation plot, with sampling if necessary
        if (iscorrelplot) {
          correlplot.AddPoint(mres, del1, del2);
        }
      }
    }
  }
  // ------------------------------------------------------------
  correl_coeff HalfDataset::CCanom() const
  // Overall
  {
    correl_coeff CC;
    for (int i=0;i<nresbin;++i) {
      CC += ccanomreso[i];
    }
    return CC;
  }
  // ------------------------------------------------------------
  correl_coeff HalfDataset::CCanomCen() const
  // Overall
  {
    correl_coeff CC;
    for (int i=0;i<nresbin;++i) {
      CC += ccanomresoCen[i];
    }
    return CC;
  }
  // ------------------------------------------------------------
  correl_coeff HalfDataset::CC_Imean() const
  // Overall
  {
    correl_coeff CC;
    for (int i=0;i<nresbin;++i) {
      CC += ccIreso[i];
    }
    return CC;
  }
  // ------------------------------------------------------------
  // CC(1/2) from variances
  double HalfDataset::CC_half(const int& mres) const
  {
    double cc = 0.0;
    double varI = meanIreso[mres].SampleVariance();
    double meanVarI = meanVarianceImeanreso[mres].Mean();
    if ((varI > 0.0) && (meanVarI > 0.0)) {
      cc = (varI - meanVarI)/(varI + meanVarI);
    }
    //std::cout << "&& "<<mres<<" "<<varI<<" "<<meanVarI<<" "<<cc<<"\n";
    return cc;
  }
  // ------------------------------------------------------------
  // CC(1/2) from variances, overall
  double HalfDataset::CC_half() const
  {
    MeanValue CC;
    for (int mres=0;mres<nresbin;++mres) {
      double varI = meanIreso[mres].SampleVariance();
      double meanVarI = meanVarianceImeanreso[mres].Mean();
      double cc = 0.0;
      if ((varI > 0.0) && (meanVarI > 0.0)) {
        cc = (varI - meanVarI)/(varI + meanVarI);
      }
      CC.Add(cc);
    }
    return CC.Mean();
  }
  // ------------------------------------------------------------
  // CC(1/2) from variances, all
  std::vector<double> HalfDataset::allCC_half() const
  {
    std::vector<double> CC;
    for (int mres=0;mres<nresbin;++mres) {
      double varI = meanIreso[mres].SampleVariance();
      double meanVarI = meanVarianceImeanreso[mres].Mean();
      double cc = 0.0;
      if ((varI > 0.0) && (meanVarI > 0.0)) {
        cc = (varI - meanVarI)/(varI + meanVarI);
      }
      CC.push_back(cc);
    }
    return CC;
  }
  // ------------------------------------------------------------
  Rfactor HalfDataset::rsplit(const int& mres) const
  {
    Rfactor R = rsplitreso.at(mres);
    R.scale(1.0/sqrt(2.0));
    return R;
  }
  // ------------------------------------------------------------
  Rfactor HalfDataset::rsplit() const
  // Overall
  {
    Rfactor rsplit;
    for (int i=0;i<nresbin;++i) {
      rsplit += rsplitreso[i];
    }
    rsplit.scale(1.0/sqrt(2.0));
    return rsplit;
  }
  // ------------------------------------------------------------
  std::vector<Rfactor> HalfDataset::allrsplit() const
  // shells
  {
    std::vector<Rfactor> rsplitv;
    for (int i=0;i<nresbin;++i) {
      rsplitv[i] = rsplitreso[i];
      rsplitv[i].scale(1.0/sqrt(2.0));
    }
    return rsplitv;
  }
  // ------------------------------------------------------------
  double HalfDataset::RMScorrelRatio(const int& mres) const
  {
    if (rmsError.at(mres).Count() <= 2) return 0.0;
    return rmsCorrel.at(mres).SD()/rmsError.at(mres).SD();
  }
  // ------------------------------------------------------------
  double HalfDataset::RMScorrelRatio() const
  // Overall
  {
    MeanSD C, E;
    for (int i=0;i<nresbin;++i) {
      C += rmsCorrel.at(i);
      E += rmsError.at(i);
    }
    if (E.Count() <= 2)  return 0.0;
    return C.SD()/E.SD();
  }
  // ------------------------------------------------------------
  double HalfDataset::RMScorrelRatioCen(const int& mres) const
  {
    if (rmsErrorCen.at(mres).Count() <= 2) return 0.0;
    return rmsCorrelCen.at(mres).SD()/rmsErrorCen.at(mres).SD();
  }
  // ------------------------------------------------------------
  double HalfDataset::RMScorrelRatioCen() const
  // Overall
  {
    MeanSD C, E;
    for (int i=0;i<nresbin;++i) {
      C += rmsCorrelCen.at(i);
      E += rmsErrorCen.at(i);
    }
    if (E.Count() <= 2)  return 0.0;
    return C.SD()/E.SD();
  }
  // ------------------------------------------------------------
  std::string HalfDataset::PlotCorrel(const bool& plotxmgr) const
  // plot stuff, returns XML plot
  {
    std::string s;
    if (iscorrelplot) {
      FILE* correlplotfile = NULL;
      if (plotxmgr) {
        correlplotfile = OpenFile("CORRELPLOT", true);
      }
      s = correlplot.Plot(correlplotfile);
      if (plotxmgr) {fclose (correlplotfile);}
    }
    return s;
  }
  // ------------------------------------------------------------
  void HalfDataset::AddAniso(const int& mres, const int& jaxis,
                             const double& wt,
                             SelectedObservations& allobs)
  // Add into sums, for anisotropy analysis along three directions
  {
    float I1, I2;
    if (mres == 0) {
      // inner resolution bin, use all data for all directions
      if (allobs.HalfAverages(I1, I2)) {
        for (int j=0;j<3;++j) {
          ccaniso[j][mres].add(I1,I2,wt);  // analyis by axis and resolution
        }
      }
    } else if (jaxis >= 0) { // near axis
      if (allobs.HalfAverages(I1, I2)) {
        ccaniso[jaxis][mres].add(I1,I2,wt);  // analysis by axis and resolution
      }
    }
  }
  // ------------------------------------------------------------
  void HalfDataset::AddAnisoProjection(const IVect3& anisores,
                                       SelectedObservations& allobs,
                                       const float& normscale)
  // Add into sums, for anisotropy analysis by projection
  // along three directions
  // anisores are 3 projected resolution bins along the principle axes
  // normscale is scale to multiply I to E^2
  {
    float I1, I2;
    if (allobs.HalfAverages(I1, I2)) {
      I1 *= normscale;
      I2 *= normscale;
      //^      std::cout << "I12 " << I1 <<"  " <<I2<<"\n"; //^
      for (int i=0;i<3;++i) {
        ccanisoprj[i][anisores[i]].add(I1,I2);  // analyis by axis and resolution
      }
    }
  }
  // ------------------------------------------------------------
  // for axis and resolution
  correl_coeff HalfDataset::CCaniso(const int& jaxis, const int& mres) const
  {
    return ccaniso.at(jaxis).at(mres);
  }
  // ------------------------------------------------------------
  correl_coeff HalfDataset::CCaniso(const int& jaxis) const  // overall
  {
    correl_coeff CC;
    for (int i=0;i<nresbin;++i) {
      CC += ccaniso.at(jaxis)[i];
    }
    return CC;
  }
  // ------------------------------------------------------------
  // for axis and resolution, projection
  correl_coeff HalfDataset::CCanisoProjection(const int& jaxis, const int& mres) const
  {
    return ccanisoprj.at(jaxis).at(mres);
  }
  // ------------------------------------------------------------
  // overall
  correl_coeff HalfDataset::CCanisoProjection(const int& jaxis) const
  {
    correl_coeff CC;
    for (int i=0;i<nresbin;++i) {
      CC += ccanisoprj.at(jaxis)[i];
    }
    return CC;
  }
  // ------------------------------------------------------------
  void HalfDataset::Analyse(const ResoRange& ResRange,
                            const double& MinimumHalfdatasetCC,
                            const double& MinimumHalfdatasetCCanom)
  // Determine resolution "limits" from half-dataset CCs
  // Also clear directions if there is only a value in the 1st bin
  {
    // Anisotropic analysis along 3 axes
    for (int jax=0;jax<3;++jax) {  // Loop directions
      ASSERT (int(ccaniso[jax].size()) == ResRange.Nbins());
      std::vector<double> cc(ccaniso[jax].size(),0.0);
      bool OK = false;
      for (size_t i=0;i<ccaniso[jax].size();++i) {
        if (i>0 && ccaniso[jax][i].result().count > 0) {
          OK = true;
        }
        cc[i] = ccaniso[jax][i].result().val;
      }
      if (!OK) {
        // clear directions if there is only a value in the 1st bin
        for (size_t i=0;i<ccaniso[jax].size();++i) {
          ccaniso[jax][i].zero();
        }
      } else {
        anisoresolimit[jax].init(cc, ResRange, MinimumHalfdatasetCC,
                                 ResolutionLimit::TANH);
      }
    }

    // Overall values
    std::vector<double> cc(ccIreso.size(),0.0);
    for (size_t i=0;i<ccIreso.size();++i) {
      cc[i] = ccIreso[i].result().val;
    }
    overallresolimit.init(cc, ResRange, MinimumHalfdatasetCC,
                          ResolutionLimit::TANH);

    // for anomalous
    cc.assign(ccanomreso.size(),0.0);
    for (size_t i=0;i<ccanomreso.size();++i) {
      cc[i] = ccanomreso[i].result().val;
    }
    anomresolimit.init(cc, ResRange, MinimumHalfdatasetCCanom,
                          ResolutionLimit::TANH);
  }
  // ------------------------------------------------------------
  // ------------------------------------------------------------
  CumulativeCChalf::CumulativeCChalf(const ResoRange& ResRange,
                                     const int& nCumulativeResoBins,
                                     const Batchgroup& batchgroup)
  {
    resrange = ResRange;
    int nresbin =  resrange.Nbins();  // number of bins for initial analysis
    pbatchgroup = &batchgroup;

    // Set number of resolution ranges: the full number nresbin
    // are used to accumulate data, but then averaged down to a smaller number
    // at the end, since initially CCs need to be in reasonably narrow resolution ranges
    if (nCumulativeResoBins > 0) {
      ncumulativeresobins = nCumulativeResoBins;
    } else {
      // typically 5 bin
      const int NRESCUMBIN = 5;
      ncumulativeresobins = std::min(nresbin, NRESCUMBIN);
    }

    // for each batch group
    meanI.resize(nresbin);
    meanVarianceImean.resize(nresbin);
    meanIdelta.resize(nresbin);
    meanVarianceImeandelta.resize(nresbin);
    int ngroups = std::max(1, pbatchgroup->numberofgroups());
    for (int k=0;k<nresbin;++k) {
      meanI[k].resize(ngroups);
      meanVarianceImean[k].resize(ngroups);
      meanIdelta[k].resize(ngroups);
      meanVarianceImeandelta[k].resize(ngroups);
    }
  }
  // ------------------------------------------------------------
  bool compareIwgrp(const Iwgrp& p1, const Iwgrp& p2)
  {return (p1.batch < p2.batch);}
  // ------------------------------------------------------------
  void CumulativeCChalf::addreflection(const int& mres,
                                       const SelectedObservations& selobs)
  // add in variance contributions for one reflection
  {
    if (selobs.Number() < 2) {return;}

    if ((mres < 0) || (mres >= meanI.size())) {
      //      std::cout << "CCcum mres "<<mres <<"\n";
      ASSERT (!((mres < 0) || (mres >= meanI.size())));
    }
    // List of (I, weight, Batch)
    std::vector<Iwgrp> iwgp = selobs.iwgroup();
    // sort into order of batch (same order as batchgroup)
    std::sort(iwgp.begin(), iwgp.end(), compareIwgrp);

    // Get successive variances for items 1 to N,
    // store them according to last batch group index
    MeanVariance mvp;
    for (size_t k=0; k<iwgp.size(); k++) {
      mvp.Add(iwgp[k].I, iwgp[k].w);
      int batchgroup = pbatchgroup->batchgroup(iwgp[k].batch);
      meanI[mres][batchgroup].Add(mvp.Mean());
      //^^
      //      if (mres == meanI.size()-1 &&  k == iwgp.size()-1) {
      //        // print last bit
      //        std::cout <<"@@ "<<mres<<" "<<k<<" "<<mvp.Mean()<<"\n";
      //      } //^-
      if (k > 0) {
        if ((batchgroup > 0) && (batchgroup < meanI[mres].size())) {
          //^
          //      if (mres == meanI.size()-3) {
          //        std::cout <<"CCcum bgp " <<k<<" "<<batchgroup<<
          //          " "<<mvp.Mean()<<" "<<mvp.VarianceofMean()<<"\n";
          //      }
          //^
          meanVarianceImean[mres][batchgroup].Add(mvp.VarianceofMean()); // variance so far
        }
      }
    }
    // Sums for Delta(CC(1/2))
    //  mvp is MeanVariance for all observations of this reflection
    MeanVariance delmvp;
    for (size_t k=0; k<iwgp.size(); k++) {
      int batchgroup = pbatchgroup->batchgroup(iwgp[k].batch);
      // remove this observation
      delmvp = mvp;
      delmvp.Subtract(iwgp[k].I, iwgp[k].w);
      meanIdelta[mres][batchgroup].Add(delmvp.Mean());
      meanVarianceImeandelta[mres][batchgroup].Add(delmvp.VarianceofMean());
    }
    //    std::cout << "CumVar " << mvp.VarianceofMean()<<"\n";
  }
  // ------------------------------------------------------------
  std::vector<std::vector<MeanValue> > CumulativeCChalf::getcumulativeCC()
  // cc(1/2) = (varI - meanVarI)/(varI + meanVarI);

  {
    // Map accumulated resolution bins on to smaller number for output ncumulativeresobins)
    int nbatchgroups = pbatchgroup->numberofgroups();
    resrangecoarse = resrange;
    resrangecoarse.SetNbins(ncumulativeresobins);

    // for output
    cumulativeCC.resize(ncumulativeresobins);
    for (int k=0;k<ncumulativeresobins;++k) {
      cumulativeCC[k].resize(nbatchgroups);
    }

    for (size_t j=0; j<meanVarianceImean.size(); j++) { // loop "fine" reso bins
      double smid = resrange.middle(j);  // middle of resolution range
      int mrescoarse = resrangecoarse.bin(smid);
      for (int i=0;i<nbatchgroups;++i) { // loop batch groups
        double varI = meanI[j][i].SampleVariance(); // variance of mean <I>
        double meanVarI = meanVarianceImean[j][i].Mean(); // mean sample variance
        //      std::cout << "j,i, vI, n, mVI, n "<<j<<" "<<i<<" "
        //                << varI<<" "<<meanI[j][i].Count()<<" "
        //                <<meanVarI<<" "<<meanVarianceImean[j][i].Count() <<std::endl;
        double cc = 0.0;
        if ((varI > 0.0) && (meanVarI > 0.0)) {
          cc = (varI - meanVarI)/(varI + meanVarI);
        }
        // CC for this batch group i, coarse resolution
        cumulativeCC[mrescoarse][i].Add(cc);
        //      std::cout << " BG "<<i<< "  "<<cc;
      }
      //      std::cout<<"\n";

    }

    //^
    //    for (int i=0;i<nbatchgroups;++i) { // loop batch groups
    //      std::cout << "CC cum " << i;
    //      for (size_t mrescoarse=0; mrescoarse<ncumulativeresobins; mrescoarse++) {
    //        std::cout << "  "<<cumulativeCC[mrescoarse][i].Mean();
    //      }
    //      std::cout <<"\n";
    //    }

    return cumulativeCC;
  }
  // ------------------------------------------------------------
  void CumulativeCChalf::printcumulativeCC(const int& datasetIndex,
                                           const std::vector<Run>& RunList,
                                           phaser_io::Output& output)
  {
    if (cumulativeCC.size() == 0) {
      getcumulativeCC();
    }
    output.logTab(0,LOGFILE,
                  std::string
                  ("\nCumulative CC(1/2) analysed by Batch for each dataset\n")+
                     "=====================================================\n\n");
    output.logTab(0,LOGFILE, pbatchgroup->format());

    int nbatchgroups = pbatchgroup->numberofgroups();
    std::vector<double> cc(ncumulativeresobins);

    // no symbol if too many points
    int symbolsize = GetSymbolSizeforNpoints(nbatchgroups);

    // $TABLE  start
    TableGraph table(" Cumulative CC(1/2) in resolution ranges vs. Batch");
    table.StoreID("Graph-CumulativeCC(1/2)");

    TableGraphPlot graph("Cumulative CC(1/2)");
    std::string description =
      "Cumulative CC(1/2) in coarse resolution ranges";
    graph.SetDescription(description);

    std::vector<std::string> collabels;
    std::string format = "%5d%5d%8d";
    collabels.push_back("N");         // 1
    collabels.push_back("Run");       // 2
    collabels.push_back("Batch");     // 3
    int jcol = 4;

    // Resolution ranges
    int nresobins = resrangecoarse.Nbins();
    int m0 = std::min(nresobins-1, 1); // usually skip first bin
    for (int mres=m0;mres<nresobins;++mres) {
      RPair reslowhigh = resrangecoarse.boundsA(mres);
      std::string sr = StringUtil::ftos(reslowhigh.first,7,2)+"-"+
        StringUtil::ftos(reslowhigh.second,7,2);
      const int FWIDTH=10;
      sr = StringUtil::CentreString(StringUtil::Strip(sr), FWIDTH);
      collabels.push_back(sr);
      format += "%"+StringUtil::itos(FWIDTH,2)+".2f";
      graph.AddLine(TableGraphPlotline(1,jcol,"","",symbolsize));  //
      jcol++;
    }

    graph.SetYaxis("", true);  // Y from zero
    // Breaks in X axis
    int xcolbr = 4;  // column for real batch number
    Xbreaks xbreaks(*pbatchgroup, datasetIndex);
    std::vector<Range> xbreaklist = xbreaks.get_breaks();
    graph.SetXbreak(xcolbr, xbreaklist);

    Range xrange(xbreaks.get_batchnumberrange());  // overall batch number range
    graph.SetXaxis("", false, xrange, true);
    table.AddGraph(graph);

    int nc = collabels.size();
    std::vector<bool> Zero(nc, true);
    table.StoreColumnFields(collabels, Zero,format+"\n");
    int ncr = 3;

    std::vector<Batch> batches = pbatchgroup->batches();

    for (int i=0;i<nbatchgroups;++i) { // loop batch groups
      int ib = pbatchgroup->batchserial(i);
      std::vector<double> ccr;
      for (int mres=m0;mres<nresobins;++mres) {
        ccr.push_back(cumulativeCC[mres][i].Mean());
      }
      table.Line(ccr, ncr, ib+1,
                 RunList[batches[ib].RunIndex()].RunNumber(),
                 batches[ib].num());
    }
    table.CloseTable();
    output.logTab(0,LOGFILE, "\n"+table.format());
    output.logTab(0,LOGFILE,table.RawLabels());
    output.logTab(0,LXML,table.XMLformat());
  }
  // ------------------------------------------------------------
}
