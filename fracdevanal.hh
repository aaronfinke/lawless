// fracdevanal.hh


//  Analysis of fractional deviations in intensity bins


#ifndef FRACDEVANAL_HEADER
#define FRACDEVANAL_HEADER

#include <vector>
#include "Output.hh"
#include "range.hh"
#include "hkl_datatypes.hh"
#include "score_datatypes.hh"
#include "intensitybin.hh"

namespace scala
{
class FracDeviationAnalysis
{
public:
  FracDeviationAnalysis() {}
  FracDeviationAnalysis(const IntensityBin& Irange);

  // Add in one deviation for bin corresponding to intensity I
  //  D = I - <I>
  //  S = sd(I)
  //  w = weight
  void AddDeviation(const float& I, const float& D, const float& S, const float& w);
  void AddDeviation(const float& I, const RPair& DS, const float& w)
  {AddDeviation(I,DS.first,DS.second,w);}
  void AddDeviation(const float& I, const float& delta);

  void clear();

  MeanSD TotalDeviation() const {return total_deviation;}
  std::vector<MeanSD> Deviations() const {return deviations;}
  int Number() const {return n_resid;}

  float MeanHalfDelta2() const
  {return (n_resid > 0) ? 0.5*float(sum_delta2/double(n_resid)) : 0.0;}

  float MeanLnS() const
  {return (n_resid > 0) ? float(sum_lnS/double(n_resid)) : 0.0;}

  // Slope of SD(delta) v. intensity
  float Slope() const;

  // Print results
  void Print(const std::string& label, phaser_io::Output& output) const;

private:
  IntensityBin irange;
  int Nbin;
  // Totals
  MeanSD total_deviation;
  double sum_delta2;
  double sum_lnS;
  double sum_w;
  int n_resid;

  // for each intensity bin
  std::vector<MeanSD> deviations;

};
}
#endif
