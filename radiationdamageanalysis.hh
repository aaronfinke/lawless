//
// radiationdamageanalysis.hh
//

#ifndef RADIATIONDAMAGEANALYSIS_HEADER
#define RADIATIONDAMAGEANALYSIS_HEADER

#include "hkl_unmerge.hh"
#include "score_datatypes.hh"
#include "batchgroup.hh"

namespace scala {
//--------------------------------------------------------------
  class RadiationDamageAnalysis {

    // At present this analysis is done only if there is a single run
    // Note that this analysis will not be useful if the multiplicity is low

    // Rcp is the cumulative pairwise residual devised by Graeme Winter for the program CHEF,
    // inspired by Diederichs Rd statistic (Acta Cryst D62, 96-101 (2005))

    //         Rcp(k) = Sum(||Ii - Ij||)/Sum(0.5*(Ii + Ij))
    // where i & j are the batch group numbers (proxy for radiation dose) and k = Max(i, j)
    // ie a pairwise R-factor up to batch group  k

    //CmPoss is cumulative completeness

  public:
    RadiationDamageAnalysis(){}

    RadiationDamageAnalysis(const hkl_unmerge_list& hkl_list,
			    const int& jrun,
			    const Batchgroup& batchgroup);

    void init(const hkl_unmerge_list& hkl_list,
	      const int& jrun,
	      const Batchgroup& batchgroup);

    void plot(const std::vector<float>& batchcompleteness,
	      phaser_io::Output& output) const;

  private:
    int irun;    // run serial number
    int runnum;  // run number
    const Batchgroup* pbatchgroup;

    double phibinsize;  // bin size in degrees
    int batch0;     // 1st batch number
    int nresbin;
    int ntimebin;
    std::vector<std::vector<Rfactor> > rfactor;  // resolution, time
    std::vector<std::vector<correl_coeff> > cc;

    // get index of this batch in the run, allowing for rejected batches
    int batchIndex(const int& batchnum) const;
  };
//--------------------------------------------------------------
} // namespace scala
#endif
