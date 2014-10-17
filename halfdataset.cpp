// halfdataset.hh


#include "halfdataset.hh"
#include "file_util.hh"
#include "plotfiles.hh"
#include "anomdistribution.hh"

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
      // Sample if necessary
      ccIreso[mres].add(I1, I2);
      rsplitreso[mres].add((I1-I2), 0.5*(I1+I2), 1.0);
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
  std::string HalfDataset::PlotCorrel() const
  // plot stuff
  {
    std::string s;
    if (iscorrelplot) {
      FILE* correlplotfile = OpenFile("CORRELPLOT", true);
      s = correlplot.Plot(correlplotfile);
      fclose (correlplotfile);
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
        //^^    std::cout << "HDS aniso " << jax <<std::endl;
        anisoresolimit[jax].init(cc, ResRange, MinimumHalfdatasetCC,
                                 ResolutionLimit::TANH);
      }
    }

    // Overall values
    std::vector<double> cc(ccIreso.size(),0.0);
    for (size_t i=0;i<ccIreso.size();++i) {
      cc[i] = ccIreso[i].result().val;
    }
    //^^   std::cout << "HDS overall\n";
    overallresolimit.init(cc, ResRange, MinimumHalfdatasetCC,
                          ResolutionLimit::TANH);

    // for anomalous
    cc.assign(ccanomreso.size(),0.0);
    for (size_t i=0;i<ccanomreso.size();++i) {
      cc[i] = ccanomreso[i].result().val;
    }
    //^^    std::cout << "HDS anom\n";
    anomresolimit.init(cc, ResRange, MinimumHalfdatasetCCanom,
                          ResolutionLimit::TANH);
  }
  // ------------------------------------------------------------
}
