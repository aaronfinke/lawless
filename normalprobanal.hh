// normalprobanal.hh

// Normal probability analysis

// Normal probability analysis of list of normalized deltas (delta/sd(delta))

// David Smith and Lynne Howell ( J.Appl. Cryst (1992) 25, 81-86) 
// David Smith, CCP4 Study Weekend (1993), page 99-106

#ifndef NORMALPROBANAL_HEADER
#define NORMALPROBANAL_HEADER

#include <vector>
#include "hkl_datatypes.hh"
#include "score_datatypes.hh"
#include "plotfiles.hh"

class NormalProbAnal
{
public:
  NormalProbAnal()
    : dltlim(0.9), maxptd(25), Fitted(false) {}

  void Clear() {deltalist.clear(); Fitted=false;}  
  void AddDelta(const float& delta);
  void AddDelta(const std::vector<float>& deltas);
  void StoreDelta(const std::vector<float>& deltas);

  // Returns slope for central part within dltlim,
  // sets dltlim if > 0
  //   or for all data if All = true
  //  = 0,0 if no data
  float Slope(const float& DltLim=0.9, const bool& All=false);
  // Returns intercept for central part within dltlim,
  // sets dltlim if > 0
  //   or for all data if All = true
  //  = 0,0 if no data
  float Intercept(const float& DltLim=0.9, const bool& All=false);
  // Numbers
  int Number(const float& DltLim=0.0, const bool& All=false);

  // Write delta(expected), delta(obs) pairs to plotfile  
  void Plot(NormalProbPlot& NormPlot,
	    const std::string& legend);

  float DltLim() const { return dltlim;}   // current value of limit

private:
  void Fit();  // fit the lines

  // probability limit for central part of distribution,
  //                    used for slope
  //                    default 0.9, c. 63% of data
  float dltlim; 
  //  maximum point density for plot: above this number,
  //  only a sample of points will be written to the plot file
  //  default = 25
  int maxptd;
  std::vector<float> deltalist;  // delta = normalised deviations

  bool Fitted;

  float slope_all;
  float slope_sel;
  float intercept_all;
  float intercept_sel;
  int num_all;
  int num_sel;
};


#endif

