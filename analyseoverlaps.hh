// analyseoverlaps.hh
//

#include <vector>

#include "hkl_unmerge.hh"
#include "Output.hh"

namespace scala {

  class Analyseoverlaps {
  public:
    Analyseoverlaps() : numbermerged(-1) {}

    //! init does the work
    void init(hkl_unmerge_list& hkl_list, const int& datasetIndex,
	      const bool& verbose, phaser_io::Output& output);
    //<! Reclassify self-overlaps as singletons
    //<! Accumulate statistics on overlaps for this dataset (-1 for all)

    //! return true if initialised
    bool Valid() {return (numbermerged >= 0);}

    void PrintOverlapTable(phaser_io::Output& output) const;

    //! return number of observations merged
    int NumberMerged() const {return numbermerged;}

    //! return number of observations merged by resolution
    std::vector<int> NumberMergedbyResolution() const {return numbermergedres;}
    //! number of unique reflections not represented by a singleton
    std::vector<int> NumberNosingleton() const {return numbernosingletonsres;}
    //! number of acentric reflections not represented by a singleton for both I+ and I-
    std::vector<int> NumberNoanomSingletons() const {return numbernoanomsingletonsres;}

  private:
    int datasetindex;
    ResoRange resrange;
    //  number of observations merged into singletons, constructed to -1
    int numbermerged;
    std::vector<int> numbermergedres;  // by resolution
    // number of unique reflections not represented by a singleton
    std::vector<int> numbernosingletonsres;  // by resolution
    // number of unique acentric reflections not represented by a
    //  singleton for both I+ and I-
    std::vector<int> numbernoanomsingletonsres;  // by resolution
    // total number of unique reflections
    std::vector<int> numberres;
    // total number of acentricreflections
    std::vector<int> numberacentricres;
  };



} // namespace scala
