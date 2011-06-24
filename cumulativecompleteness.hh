// cumulativecompleteness.hh
//
#ifndef CUMULATIVECOMPLETENESS_HEADER
#define CUMULATIVECOMPLETENESS_HEADER

#include <vector>

#include "range.hh"
#include "hkl_symmetry.hh"
#include "hkl_datatypes.hh"

namespace scala {
  class CumulativeCompleteness {
  public:
    CumulativeCompleteness(){}
    // nbatches  number of batches in all datasets
    // ResRange  resolution range
    CumulativeCompleteness(const int& Nbatches);

    // Store multiplicity of reflection, ie the number of times this reflection
    // will occur in a complete sphere of data
    void StartReflection(const int& Multcy);
    // Record an observation in the current reflection
    //  batchn  batch number
    //  jbatch  batch serial number
    //  Anomclass ALL, IPLUS, IMINUS
    void AddObservationBatch(const int& batchn, const int& jbatch,
			     const AnomalousClass& Anomclass=ALL);
    // End of a series of AddObservationBatch calls for the current reflection
    void EndReflection();

    // return cumulative completeness for each batch serial
    std::vector<float> BatchCompleteness(const ResoRange& ResRange,
					 const hkl_symmetry& symmetry,
					 const Scell& cell);
    // anomalous completeness
    std::vector<float> BatchAnomCompleteness(const ResoRange& ResRange,
					     const hkl_symmetry& symmetry,
					     const Scell& cell);

  private:
    int nbatches;  // in all datasets
    int multcy;
    int maxBatSer;
    std::vector<int> batchserials;
    std::vector<int> batchserialsPlus;
    std::vector<int> batchserialsMinus;
    // count of reflections first appearing in this batch
    std::vector<int> batchserialcount;
    std::vector<int> batchserialcountanom;
    float nref_sphere;     // all
    float nrefAcen_sphere; // acentrics

    int LowestSerial(const std::vector<int>& bs);
    // Totals in  sphere
    void CalcSphere(const ResoRange& ResRange,
		    const hkl_symmetry& symmetry,
		    const Scell& cell);


  };
}

#endif
