// halfdataset.hh
//
// Analyses in half datasets

#ifndef HALFDATASET_HEADER
#define HALFDATASET_HEADER

#include "selectedobservations.hh"
#include "score_datatypes.hh"
#include "plotfiles.hh"
#include "anomdistribution.hh"
#include "resolutionlimit.hh"

namespace scala {

  // ------------------------------------------------------------
  class HalfDataset {
    // Half-datasets correlation scores etc
  public:
    HalfDataset(){} 
    HalfDataset(const int& NresBins, const PxdName& Dataset_pxd);

    // Store relevant anomalous statistics
    // Store RMS DelAnom for each resolution bin
    //  & and number of reflections with 2 delAnom
    void StoreAnomStats(const AnomDistribution& anomDistribution);

    // Add into sums, for Imean
    void AddMean(const int& mres, SelectedObservations& allobs);

    // Add into sums,  without anomalous
    void AddAnomCentric(const int& mres, SelectedObservations& allobs);

    // Add into sums, acentric with anomalous
    void AddAnom(const int& mres,
		 SelectedObservations& obsplus, SelectedObservations& obsminus);

    // Divide into Npart parts
    void SetNpart(const int& Npart);

    int NresBin() const {return nresbin;}

    correl_coeff CCanom(const int& mres) const {return ccanomreso.at(mres);}
    correl_coeff CCanom() const;  // overall

    correl_coeff CCanomCen(const int& mres) const {return ccanomresoCen.at(mres);}
    correl_coeff CCanomCen() const;  // overall

    correl_coeff CC_Imean(const int& mres) const {return ccIreso.at(mres);}
    correl_coeff CC_Imean() const;  // overall

    double RMScorrelRatio(const int& mres) const;
    double RMScorrelRatio() const;  // overall

    double RMScorrelRatioCen(const int& mres) const;
    double RMScorrelRatioCen() const;  // overall

    void PlotCorrel();  // plot stuff

    // Add into sums, for anisotropy analysis along three directions
    void AddAniso(const int& mres, const int& jaxis, SelectedObservations& allobs);

    // for axis and resolution
    correl_coeff CCaniso(const int& jaxis, const int& mres) const;
    correl_coeff CCaniso(const int& jaxis) const;  // overall

    //! Determine resolution "limits" from half-dataset CCs
    void Analyse(const ResoRange& ResRange,
		 const double& MinimumHalfdatasetCC);

    // return overall resolution limit calculated by Analyse
    ResolutionLimit OverallResoLimit() const {return overallresolimit;}

    // return anisotropic resolution limit calculated by Analyse
    std::vector<ResolutionLimit> AnisoResoLimits() const
    {return anisoresolimit;}

  private:
    int nresbin;   // number of resolution bins
    PxdName dataset_pxd;
    // DelAnom correlations by resolution
    std::vector<correl_coeff> ccanomreso;
    std::vector<correl_coeff> ccanomresoCen;
    // MeanI correlations by resolution
    std::vector<correl_coeff> ccIreso;
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

    // For anisotropy analysis
    // outer: axis 0-2; inner, resolution bins
    std::vector<std::vector<correl_coeff> > ccaniso;

    ResolutionLimit overallresolimit;
    std::vector<ResolutionLimit> anisoresolimit; // for axes 0-2

  };

}

#endif
