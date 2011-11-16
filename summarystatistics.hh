//
// summarystatistics.hh
//
// Overall summary data class (== "Table 1")
//

#ifndef SUMMARYSTATISTICS_HEADER
#define SUMMARYSTATISTICS_HEADER

#include <vector>
#include "range.hh"
#include "score_datatypes.hh"
#include "Output.hh"
#include "resolutionlimit.hh"

namespace scala {
  class SummaryStatistics
  // A class to accumulate all the things needed for a final summary
  // even if they are also stored elsewhere
  {
  public:
    SummaryStatistics():Anom(false){}
    void SetAnom(const bool& anom) {Anom = anom;}

    // store PXD name
    void StorePXDname(const PxdName& PXDname) {pxdname = PXDname;}
    // store resolution ranges, overall, inner shell, outer shell
    void StoreResRanges(const ResoRange& overall, const ResoRange& inner,
			const ResoRange& outer);
    // store Rmerge, overall, inner shell, outer shell
    void StoreRmergeReso(const Rfactor& overall, const Rfactor& inner,
			 const Rfactor& outer);
    // store Rmeas, overall, inner shell, outer shell
    void StoreRmeasReso(const Rfactor& overall, const Rfactor& inner,
			 const Rfactor& outer);
    // store Rpim, overall, inner shell, outer shell
    void StoreRpimReso(const Rfactor& overall, const Rfactor& inner,
			 const Rfactor& outer);
    // store Rmerge, overall, inner shell, outer shell, v. overall mean I+-
    void StoreRmergeResoOv(const Rfactor& overall, const Rfactor& inner,
			   const Rfactor& outer);
    // store Rmeas, overall, inner shell, outer shell, v. overall mean I+-
    void StoreRmeasResoOv(const Rfactor& overall, const Rfactor& inner,
			  const Rfactor& outer);
    // store Rpim, overall, inner shell, outer shell, v. overall mean I+-
    void StoreRpimResoOv(const Rfactor& overall, const Rfactor& inner,
			 const Rfactor& outer);
    // Rmerge in top intensity bin
    void StoreRtopI(const Rfactor& R) {RmergeTopI = R;}

    // Numbers
    void StoreNumbers(const int& overallNobs,const int& innerNobs,const int& outerNobs,
		      const int& overallNuniq,const int& innerNuniq,const int& outerNuniq);
    // Mean(I/sd)
    void StoreMnIsd(const float& overallMnIsd, const float& innerMnIsd,
		    const float& outerMnIsd);
    // <I> half-dataset CC
    void StoreImeanCorrel(const float& overall,
			  const float& inner,
			  const float& outer);
    // Completeness & multiplicity
    void StoreCmplMult(const float& overallComplete, const float& innerComplete,
		       const float& outerComplete,
		       const float& overallMult, const float& innerMult,
		       const float& outerMult);
    // Anomalous Completeness & multiplicity
    void StoreAnomCmplMult(const float& overallComplete, const float& innerComplete,
			   const float& outerComplete,
			   const float& overallMult, const float& innerMult,
			   const float& outerMult);
    // Anomalous correlation of half-datasets
    void StoreAnomCorrel(const float& overall,
			 const float& inner,
			 const float& outer);
    // Anomalous RMS correlation ratio of half-datasets
    void StoreAnomRCR(const float& overall,
		      const float& inner,
		      const float& outer);
    // slope of anomalous normal probability
    void StoreAnomNPslope(const float& anomnpslope);
    // Average unit cell
    void StoreAverageCell(const Scell& cell) {averageCell = cell;}
    // Space group name
    void StoreSpaceGroupName(const std::string& sgname)
    {spacegroupname = sgname;}
    // Average mosaicity
    void StoreAverageMosaicity(const float& mosaic)
    {averageMosaicity = mosaic;}
    // Range of SD corrections
    void StoreSDcorrectioRange(const float& minsdcorrfulls, const float& maxsdcorrfulls,
			       const float& minsdcorrpartials, const float& maxsdcorrpartials);

    // resolution limit estimates
    //  Store
    // overall limit from half-dataset CCs
    void StoreHalfdatsetCCresolimit(const ResolutionLimit& OverallResoLimitCC);
    // overall limit from Mn(I/sd)
    void StoreMnIsigresolimit(const ResolutionLimit& OverallResoLimitIsig);
    // anisotropic limits from half-dataset CCs
    void StoreHalfdatsetCCAnisoresolimit
      (const std::vector<ResolutionLimit>& AnisoresolimitCC);
    // anisotropic limits from Mn(I/sd)
    void StoreMnIsigAnisoresolimit
      (const std::vector<ResolutionLimit>& AnisoresolimitIsig);
    // store anisotropic axis labels
    void StoreAnisoAxisLabels(const std::vector<std::string>& Anisoaxislabels);

    // print the final summary table as RESULT if Result true
    void PrintSummaryTable(const bool& Result, phaser_io::Output& output);

    PxdName pxdname;                // project, crystal, dataset
    // Various statistics for overall, inner, outer shells(3-vectors)
    std::vector<ResoRange> resRange;  // overall, inner, outer
    std::vector<Rfactor> rmerge;    // Rmerge
    std::vector<Rfactor> rmeas;     // Rmeas
    std::vector<Rfactor> rpim;      // Rpim
    std::vector<Rfactor> rmergeOv;  // Rmerge overall I+ & I-
    std::vector<Rfactor> rmeasOv;   // Rmeas
    std::vector<Rfactor> rpimOv;    // Rpim
    Rfactor RmergeTopI;             // Rmerge in top intensity bin
    std::vector<int> Nobs;          // Total number of observations
    std::vector<int> Nuniq;         // Number of unique reflections
    std::vector<float> MnIsd;       // Mn(I/sd)
    std::vector<float> Icorrelation;// <I> correlation in halfdatasets
    std::vector<float> complete;        // completeness
    std::vector<float> multiplicity;    // multiplicity
    bool Anom;                          // true if anomalous information
    std::vector<float> anomcomplete;    // Anomalous completeness
    std::vector<float> anommultiplicity; // Anomalous multiplicity
    std::vector<float> anomcorrelation; // Anomalous correlation in halfdatasets
    std::vector<float> anomRCR;         // Anomalous RMS correlation ratio in halfdatasets
    float anomNPslope;                  // slope of anomalous normal probability
    Scell averageCell;                  // average unit cell
    std::string spacegroupname;
    float averageMosaicity;
    float minSDcorrFulls, maxSDcorrFulls, minSDcorrPartials, maxSDcorrPartials;

    // Resolution limit estimates
    ResolutionLimit overallresolimitCC;   // overall, from half-dataset CCs
    ResolutionLimit overallresolimitIsig; // overall, from Mn(I/sd)
    // anisotropic, from half-dataset CCs
    std::vector<ResolutionLimit> anisoresolimitCC;
    // anisotropic, from Mn(I/sd)
    std::vector<ResolutionLimit> anisoresolimitIsig;
    std::vector<std::string> anisoaxislabels;
  };  //   class SummaryStatistics
  // ------------------------------------------------------------
  template<class T> std::vector<T> Store3val(const T& overall, const T& inner, const T& outer)
  // Store 3 values in vector 
  {
    std::vector<T> v(3);
    v[0] = overall;
    v[1] = inner;
    v[2] = outer;
    return v;
  }
  // ------------------------------------------------------------
  class AllSummaryStatistics {
    //! Summary statistics for all datasets (only relevant if > 1)
  public:
    AllSummaryStatistics(){}

    //! store statistics for one dataset
    void AddSummaryStatistics(const SummaryStatistics& summarystatistics);

    //! print the final summary table for all datasets as RESULT if Result true
    void PrintSummaryTable(const bool& Result,
			   const bool& Anom, phaser_io::Output& output);

    //! print the final summary table for one dataset (idts), as RESULT if Result true
    void PrintOneSummaryTable(const int& idts, const bool& Result, phaser_io::Output& output);    

  private:
    std::vector<SummaryStatistics> allsummarystatistics;
  };
}  // namespace scala
#endif
