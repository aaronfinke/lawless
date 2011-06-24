//
// optimisesdcorr.hh
//

#ifndef OPTIMISESDCORR_HEADER
#define OPTIMISESDCORR_HEADER

#include <vector>

#include "hkl_unmerge.hh"
#include "sdmodel.hh"
#include "intensitybin.hh"
#include "simplex-lib.h"
#include "selectedobservations.hh"
#include "analysesd.hh"
#include "sdanalysis.hh"

namespace scala {
  //---------------------------------------------------------------
  class SDresiduals
  //! a group of residuals
  {
  public:
    SDresiduals(const double& R, const int& NR,
		const std::vector<double>& cR, const std::vector<int>& cNR)
      : overallR(R), nR(NR), classRs(cR), cnR(cNR) {}

    double overallR;
    int nR;
    std::vector<double> classRs; // for each parameter group
    std::vector<int> cnR;
  };
  //---------------------------------------------------------------
class SDcorrResidual {
  // Determine residual for optimisation of SD corrections
public:
  SDcorrResidual(){}
  SDcorrResidual (SDmodel& sdm,
		  const hkl_unmerge_list& Hkl_list,
		  IntensityBin& Irange,
		  const bool& Anomalous);
  void init(SDmodel& sdm,
	    const hkl_unmerge_list& Hkl_list,
	    IntensityBin& Irange,
	    const bool& Anomalous);

  SDcorrResidual(const SDanalysis& Sdanal) {init(Sdanal);}
  void init(const SDanalysis& Sdanal);

  // Return function value given the parameter vector args
  double operator() (const std::vector<double>& args) const;

  SDresiduals Residual() const;

  void PrintTable(phaser_io::Output& output);
    
private:
  SDmodel* SDM;                      // SD model (may change)
  const hkl_unmerge_list* hkl_list;  // reflection data (const, doesn't change)
  mutable IntensityBin irange;
  bool anomalous;
  std::vector<bool> combine;  // combine fulls & partials, for each run
  int nrunsused;   // = number of runs or +1 if SDM.AllRunsSame
  mutable SDanalysis sdanal;
};  // SDcorrResidual
  //---------------------------------------------------------------
  class SDcorrRefine : public Target_fn_order_zero
  {
  public:
    SDcorrRefine(){}
    SDcorrRefine(SDmodel& sdm,
		 const hkl_unmerge_list& Hkl_list,
		 IntensityBin& Irange,
		 const bool& anomalous);

    int num_params() const {return nparams;}
    double operator() (const std::vector<double>& args) const;

    void PrintTable(phaser_io::Output& output) {SDCresid.PrintTable(output);}

  private:
    SDcorrResidual SDCresid;
    int nparams;

  }; // SDcorrRefine
  //---------------------------------------------------------------
  // Get starting values for simplex optimisation
  std::vector<std::vector<double> > StartValues(const SDmodel& SDM,
						const double& scale);
  //---------------------------------------------------------------
  // Optimise SD correction model using Simplex minimiser
  void OptimiseSDcorr(SDmodel& SDM,
		      const hkl_unmerge_list& hkl_list,
		      const all_controls& controls,
		      IntensityBin& Irange,
		      const double& tolerance,
		      const int&  max_cycles,
		      phaser_io::Output& output);

}
#endif
