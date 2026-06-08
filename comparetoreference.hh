#ifndef COMPARETOREFERENCE
#define COMPARETOREFERENCE
// comparetoreference.hh

#include <vector>
#include <string>

#include "score_datatypes.hh"
#include "range.hh"
#include "hkl_datatypes.hh"
#include "hkl_unmerge.hh"
#include "referencelist.hh"
 
namespace scala {
  class CompareToReference {
    // Analyse agreementof multiple datasets with reference set
  public:
    CompareToReference() : ndatasets(-1) {}

    CompareToReference(const hkl_unmerge_list& hkl_list,
		       const ReferenceList& hklreflist,
		       const ResoRange& ResRange);

    void printTable(phaser_io::Output& output);

  private:
    int ndatasets;
    std::vector<PxdName> pxdnames;  // datasets
    std::vector<std::string> dnames; // dataset names
    ResoRange resrange;

    std::vector<std::vector<Rfactor> > rfref;
    std::vector<std::vector<correl_coeff> > ccref;
    std::vector<Rfactor>  totalrfref;
    std::vector<correl_coeff>  totalccref;

  };
}
#endif
