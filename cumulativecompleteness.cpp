// cumulativecompleteness.cpp
//

#include "cumulativecompleteness.hh"
#include "numbercomplete.hh"

namespace scala {
  // ------------------------------------------------------------
  CumulativeCompleteness::CumulativeCompleteness(const int& Nbatches)		  
    : nbatches(Nbatches)
  {
    batchserialcount.assign(nbatches, 0);
    batchserialcountanom.assign(nbatches, 0);
    nref_sphere = -1;
    nrefAcen_sphere = -1;
    maxBatSer = -1;
  }
  // ------------------------------------------------------------
  void  CumulativeCompleteness::StartReflection(const int& Multcy)
  // Store multiplicity of reflection, ie the number of times this reflection
  // will occur in a complete sphere of data
  {
    multcy = Multcy;
    // Clear arrays
    batchserials.clear();
    batchserialsPlus.clear();
    batchserialsMinus.clear();
  }
  // ------------------------------------------------------------
  void CumulativeCompleteness::AddObservationBatch
  (const int& batchn, const int& jbatch, const AnomalousClass& Anomclass)
  // Record an observation in the current reflection
  //  batchn  batch number
  //  jbatch  batch serial number
  //  Anomclass ALL, IPLUS, IMINUS
  {
    ASSERT (jbatch <= nbatches);
    maxBatSer = Max(maxBatSer, jbatch);
    if (Anomclass == ALL) {
      batchserials.push_back(jbatch);
    } else if (Anomclass == IPLUS) {
      batchserialsPlus.push_back(jbatch);
    } else if (Anomclass == IMINUS) {
      batchserialsMinus.push_back(jbatch);
    }
  }
  // ------------------------------------------------------------
  int CumulativeCompleteness::LowestSerial(const std::vector<int>& bs)
  {
    // Find lowest batch serial, returns index 
    if (bs.size() <= 0) return -1;
    int lowb = 1000000000;
    int k = -1;
    for (size_t i=0;i<bs.size();++i) {
      if (bs[i] < lowb) {
	lowb = bs[i];
	k = bs[i];
      }
    }
    return k;
  }
  // ------------------------------------------------------------
  void CumulativeCompleteness::EndReflection()
  // End of a series of AddObservationBatch calls for the current reflection
  {
    if (batchserials.size() <= 0) return;
    int k = LowestSerial(batchserials); // ALL
    if (k >= 0) {
      batchserialcount.at(k) += multcy; // increment count for first appearance
    }
    int kp = LowestSerial(batchserialsPlus); // IPLUS
    int km = LowestSerial(batchserialsMinus); // IMINUS
    if (kp>=0 && km>=0) {
      k = Max(kp,km); // wait for the second of the Bijvoet pair
      batchserialcountanom.at(k) += multcy; // increment count for first appearance
    }
  }
  // ------------------------------------------------------------
  void CumulativeCompleteness::CalcSphere(const ResoRange& ResRange,
					 const hkl_symmetry& symmetry,
					 const Scell& cell)
  // Totals in  sphere
  {
    ResoRange resrange = ResRange;
    resrange.SetNbins(1);  // one bin
    std::vector<int> Nrefres(1);
    std::vector<int> Nrefacen(1);
    NumberComplete(resrange, symmetry, cell, Nrefres, Nrefacen);
    nref_sphere = Nrefres[0];  // number of reflections in sphere
    nrefAcen_sphere = Nrefacen[0];  // number of reflections in sphere
  }
  // ------------------------------------------------------------
  std::vector<float> CumulativeCompleteness::BatchCompleteness(const ResoRange& ResRange,
							       const hkl_symmetry& symmetry,
							       const Scell& cell)
  // return cumulative completeness for each batch serial
  {
    if (nref_sphere <= 0) CalcSphere(ResRange, symmetry, cell);
    std::vector<float> complete(nbatches, 0.0);
    int sofar = 0;
    for (int i=0;i<nbatches;++i) {
      sofar += batchserialcount[i];
      if (nref_sphere == 0) {
	complete.at(i) = 0.0;
      } else {
	complete.at(i) = float(sofar)/nref_sphere;
      }
    }
    return complete;
  }
  // ------------------------------------------------------------
  std::vector<float> CumulativeCompleteness::BatchAnomCompleteness
  (const ResoRange& ResRange,
   const hkl_symmetry& symmetry,
   const Scell& cell)
  // return cumulative anomalous completeness for each batch serial
  {
    if (nrefAcen_sphere <= 0) CalcSphere(ResRange, symmetry, cell);
    std::vector<float> complete(nbatches, 0.0);
    int sofar = 0;
    for (int i=0;i<nbatches;++i) {
      sofar += batchserialcountanom[i];
      if (nrefAcen_sphere == 0) {
	complete.at(i) = 0.0;
      } else {
	complete.at(i) = float(sofar)/nrefAcen_sphere;
      }
    }
    return complete;
  }

}
