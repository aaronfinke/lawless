//
// radiationdamageanalysis.hh
//

#ifndef RADIATIONDAMAGEANALYSIS_HEADER
#define RADIATIONDAMAGEANALYSIS_HEADER

#include "hkl_unmerge.hh"
#include "score_datatypes.hh"

namespace scala {
//--------------------------------------------------------------
  class RadiationDamageAnalysis {
  public:
    RadiationDamageAnalysis(){}

    RadiationDamageAnalysis(const hkl_unmerge_list& hkl_list,
			    const int& jrun, const int& nbatchgroup=0);

    void init(const hkl_unmerge_list& hkl_list,
	      const int& jrun, const int& nbatchgroup=1);

    void plot(const std::vector<float>& batchcompleteness,
	      phaser_io::Output& output) const;

  private:
    int irun;    // run serial number
    int runnum;  // run number
    int batchgroup; // number of batches in group (usually 1)

    hash_table batch_lookup;  // batch lookup for index into run

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
