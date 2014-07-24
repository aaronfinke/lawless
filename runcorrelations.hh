// runcorrelations.hh
//
// Matrix of correlations between runs


#ifndef RUNCORRELATIONS_HEADER
#define RUNCORRELATIONS_HEADER

#include "hkl_unmerge.hh"
#include "sdmodel.hh"
#include "normalise.hh"
#include "score_datatypes.hh"

namespace scala {

  class RunCorrelations {
  public:
    RunCorrelations(){}
    RunCorrelations(const hkl_unmerge_list& hkl_list,
		    const SDmodel& SDM,
		    const Normalise& NormRes,
		    const int& nresbins);

    void init(const hkl_unmerge_list& hkl_list,
	      const SDmodel& SDM,
	      const Normalise& NormRes,
	      const int& nresbins);

    void add(const reflection& refl);

    void formatTable(phaser_io::Output& output) const;


  private:
    int nruns;
    const hkl_unmerge_list* hkl_list_p;
    const Normalise* normres;  // pointer to Normalisation object
    std::vector<double> mnCC_E2; // CC(E^2) for each pair
    std::vector<int> nmeanCC;

    std::vector<std::vector<correl_coeff> > CC_E2_res;  // on E^2

    // Diifferent weightings for average over resolution bins (hard-wired)
    int resweighttype;  // 0 = unit, +1 Var(CC), +2 Number

    ResoRange resrange;

    void resolutiongraph(const int& nrungraph,
			 phaser_io::Output& output) const;

    // output all CC by resolution to XML
    void writeCCtoXML(phaser_io::Output& output) const;

  };

} // namespace scala

#endif
