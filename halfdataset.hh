// halfdataset.hh
//
// Analyses in half datasets

#ifndef HALFDATASET_HEADER
#define HALFDATASET_HEADER

#include "selectedobservations.hh"
#include "score_datatypes.hh"
#include "plotfiles.hh"
///#include "anomdistribution.hh"
#include "resolutionlimit.hh"
#include "batchgroup.hh"
#include "Output.hh"

namespace scala {

  class AnomDistribution;

  // ------------------------------------------------------------
  class HalfDataset {
    // Half-datasets correlation scores etc
  public:
    HalfDataset() : iscorrelplot(false) {} 
    HalfDataset(const int& NresBins, const PxdName& Dataset_pxd);

    void init(const int& NresBins, const PxdName& Dataset_pxd);

    // Store relevant anomalous statistics
    // Store RMS DelAnom for each resolution bin
    //  & and number of reflections with 2 delAnom
    void StoreAnomStats(const AnomDistribution& anomDistribution);

    void StoreRMS(std::vector<MeanSD>& RMSdelanom);

    // Add into sums, for Imean
    void AddMean(const int& mres, SelectedObservations& allobs);

    // Add into sums,  without anomalous
    void AddAnomCentric(const int& mres, SelectedObservations& allobs);

    // Add into sums, acentric with anomalous
    void AddAnom(const int& mres,
		 SelectedObservations& obsplus, SelectedObservations& obsminus);

    int NresBin() const {return nresbin;}

    correl_coeff CCanom(const int& mres) const {return ccanomreso.at(mres);}
    correl_coeff CCanom() const;  // overall
    std::vector<correl_coeff> allCCanom() const {return ccanomreso;}

    correl_coeff CCanomCen(const int& mres) const {return ccanomresoCen.at(mres);}
    correl_coeff CCanomCen() const;  // overall

    correl_coeff CC_Imean(const int& mres) const {return ccIreso.at(mres);}
    correl_coeff CC_Imean() const;  // overall
    std::vector<correl_coeff> allCC_Imean() const {return ccIreso;}

    // CC(1/2) from variances
    double CC_half(const int& mres) const;
    double CC_half() const;  // overall average
    std::vector<double> allCC_half() const;

    Rfactor rsplit(const int& mres) const;
    Rfactor rsplit() const; // overall
    std::vector<Rfactor> allrsplit() const;


    double RMScorrelRatio(const int& mres) const;
    double RMScorrelRatio() const;  // overall

    double RMScorrelRatioCen(const int& mres) const;
    double RMScorrelRatioCen() const;  // overall

    // plot stuff, returns XML, writes xmgr file if plotxmgr true
    std::string PlotCorrel(const bool& plotxmgr) const;

    // Add into sums, for anisotropy analysis along three directions
    void AddAniso(const int& mres, const int& jaxis,
    			     const double& wt, SelectedObservations& allobs);
    // Add into sums, for anisotropy analysis by projection along three directions
    // anisores are 3 projected resolution bins along the principle axes
    // normscale is scale to multiply I to E^2
    void AddAnisoProjection(const IVect3& anisores,
			    SelectedObservations& allobs,
			    const float& normscale);
    
    // for axis and resolution
    correl_coeff CCaniso(const int& jaxis, const int& mres) const;
    correl_coeff CCaniso(const int& jaxis) const;  // overall
    correl_coeff CCanisoProjection(const int& jaxis, const int& mres) const;
    correl_coeff CCanisoProjection(const int& jaxis) const;  // overall

    //! Determine resolution "limits" from half-dataset CCs
    void Analyse(const ResoRange& ResRange,
		 const double& MinimumHalfdatasetCC,
		 const double& MinimumHalfdatasetAnomCC);

    // return overall resolution limit calculated by Analyse
    ResolutionLimit OverallResoLimit() const {return overallresolimit;}

    // return anisotropic resolution limit calculated by Analyse
    std::vector<ResolutionLimit> AnisoResoLimits() const
    {return anisoresolimit;}

    // return anomalous resolution limit calculated by Analyse
    ResolutionLimit AnomalousResoLimit() const {return anomresolimit;}

  private:
    int nresbin;   // number of resolution bins
    PxdName dataset_pxd;
    // DelAnom correlations by resolution
    std::vector<correl_coeff> ccanomreso;
    std::vector<correl_coeff> ccanomresoCen;
    // MeanI correlations by resolution
    std::vector<correl_coeff> ccIreso;
    std::vector<MeanVariance> meanIreso; // for Var(<I>) by resolution
    std::vector<MeanValue> meanVarianceImeanreso; // <var(Isample)>

    // Rsplit by resolution
    std::vector<Rfactor> rsplitreso;
    // RMS DelAnom for each resolution bin
    std::vector<float> rmsdelanom;
    float rmsdelanomOverall;  // overall value
    // number of delAnom pairs in both half sets
    int nAnomPairs;
    // Limit on delAnom as multiple of RMS in resolution range
    float maxDelAnom;
    // Sums for RMS correlation ratio
    std::vector<MeanSD> rmsCorrel;  // cos45 * (Del1 + Del2)
    std::vector<MeanSD> rmsError;   // cos45 * (Del1 - Del2)
    // Sums for RMS correlation ratio
    std::vector<MeanSD> rmsCorrelCen;  // cos45 * (Del1 + Del2)
    std::vector<MeanSD> rmsErrorCen;   // cos45 * (Del1 - Del2)

    CorrelPlot correlplot;
    bool iscorrelplot;

    // For anisotropy analysis
    // outer: cone axis 0-2 ; inner, resolution bins
    std::vector<std::vector<correl_coeff> > ccaniso;
    // outer: projection axis 0-2 ; inner, projected resolution bins
    std::vector<std::vector<correl_coeff> > ccanisoprj;

    ResolutionLimit overallresolimit;
    std::vector<ResolutionLimit> anisoresolimit; // for axes 0-2
    // for anomalous
    ResolutionLimit anomresolimit;

  };
  // ------------------------------------------------------------
  class CumulativeCChalf {
  public:
    CumulativeCChalf (){}
    CumulativeCChalf(const ResoRange& ResRange,
		     const int& nCumulativeResoBins,
		     const Batchgroup& batchgroup);

    void addreflection(const int& mres,
		       const SelectedObservations& selobs);

    // calculate results
    std::vector<std::vector<MeanValue> > getcumulativeCC();

    // print them
    void printcumulativeCC(const int& datasetIndex,
			   const std::vector<Run>& RunList,
			   phaser_io::Output& output);

  private:
    ResoRange resrange;  // local copy
    ResoRange resrangecoarse;  // coarse ranges
    int ncumulativeresobins;
    const Batchgroup* pbatchgroup;
    // <var(Iav)> for each resolution range for each batch
    // for Var(<I>) by resolution
    std::vector<std::vector<MeanVariance> > meanI;
    // <var(Isample)> for each resolution range for each batch
    std::vector<std::vector<MeanValue> > meanVarianceImean;
    // for Delta(CC(1/2))
    // <var(Iav)> for each resolution range for each batch
    // for Var(<I>) by resolution
    std::vector<std::vector<MeanVariance> > meanIdelta;
    // <var(Isample)> for each resolution range for each batch
    std::vector<std::vector<MeanValue> > meanVarianceImeandelta;

    std::vector<std::vector<MeanValue> > cumulativeCC;

  };
  // ------------------------------------------------------------
}

#endif
