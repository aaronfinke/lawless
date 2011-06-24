//
// anomdistribution.hh
//

#ifndef ANOMDISTRIBUTION_HEADER
#define ANOMDISTRIBUTION_HEADER

#include "hkl_unmerge.hh"
#include "controls.hh"
#include "sdmodel.hh"
#include "range.hh"
#include "normalise.hh"
#include "scala_util.hh"

namespace scala {

  class AnomDistribution {
    // Distribution of anomalous differences to get estimate
    // of maximum likely values etc
    // for each dataset
  public:
    AnomDistribution(){}

    void SetNresbin(const int& Nresbin);

    // Store delAnom, & count reflections used for
    // half-dataset correlations (ie with n+ & n- > 1, correlAnom true)
    //  for resolution range mres
    void Add(const int& mres,
	   const float& delAnom, const bool& correlAnom);

    // RMS values for datasets Idts
    std::vector<MeanSD> RMSdelAnom() const {return rmsDelAnom;}

    // Number of reflections with half-dataset anomalous correlation info
    std::vector<int> NdelAnomHalf() const {return nDelAnom;}

  private:
    int nresbin;
    // by resolution
    std::vector<MeanSD> rmsDelAnom;
    std::vector<int> nDelAnom;
    
  };
  // ------------------------------------------------------------
class AllAnomDistributions
// Anomalous distributions for all datasets
{
public:
  AllAnomDistributions(){}
  AllAnomDistributions(const hkl_unmerge_list& hkl_list,
		       const SDmodel& SDM,
		       const all_controls& controls,
		       const ResoRange& ResRange,
		       const Normalise& NormRes);

  // Store delAnom, & count reflections used for
  // half-dataset correlations (ie with n+ & n- > 1, correlAnom true)
  //  for dataset idts & resolution range mres
  void Add(const int& idts, const int& mres,
	   const float& delAnom, const bool& correlAnom)
  {anomdistributions.at(idts).Add(mres, delAnom, correlAnom);}

  //! Add in to correlation sums, resolution bin mres
  void AddCorrelations(const std::vector<float>& danomdts,
		       const std::vector<float>& Imeandts,
		       const int& mres);

  AnomDistribution Anomdistribution(const int& idts) const
  {return anomdistributions.at(idts);}

  // Print correlation tables
  void Print(phaser_io::Output& output) const;

private:
  int ndatasets;
  std::vector<PxdName> pxdnames;  // datasets
  std::vector<float> wavelengths; // for each dataset
  std::vector<std::string> dnames; // dataset names
  std::vector<AnomDistribution> anomdistributions;
  int basedataset;  // base dataset index
  ResoRange resrange;
  std::vector<std::vector<correl_coeff> > cca;  // anom differences
  std::vector<std::vector<correl_coeff> > ccd;  // dispersive differences
  // dataset index pairs for each CC in cca
  std::vector<std::pair<int,int> > ccadtsindex;
  // dataset index pairs for each CC in ccd
  std::vector<std::pair<int,int> > ccddtsindex;

  // Private methods
  // diff = true for dispersive differences
  std::string FormatTable(const std::string& title,
			  const std::string& graphtitle,
			  const std::string& ccl1,
			  const std::string& ccl2,
			  const std::vector<std::vector<correl_coeff> >& cc,
			  const std::vector<std::pair<int,int> >& ccidx,
			  const bool& diff,
			  std::vector<correl_coeff>& allcc) const;

  // format CC of anomalous differences between datasets as table
  std::string CrossCorrelation(const std::string& title,
			       const std::vector<correl_coeff>& allcc,
			       const bool& diff) const;


};
}
#endif
