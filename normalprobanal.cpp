// normalprobanal.cpp

#include "normalprobanal.hh"
#include "normprobfunc.hh"

//-------------------------------------------------------------
void NormalProbAnal::AddDelta(const float& delta)
{
  if (delta != 0.0) deltalist.push_back(delta);
}
//-------------------------------------------------------------
void NormalProbAnal::AddDelta(const std::vector<float>& deltas)
{
  for (std::vector<float>::const_iterator i=deltas.begin();i!=deltas.end();i++)
    {if (*i != 0.0) deltalist.push_back(*i);}
}
//-------------------------------------------------------------
void NormalProbAnal::StoreDelta(const std::vector<float>& deltas)
{
  deltalist.clear();
  AddDelta(deltas);
}
//-------------------------------------------------------------
void NormalProbAnal::Fit()
// Fit points to line
{
  int Nobs = deltalist.size();
  if (Nobs == 0)
    {
      slope_all = 0.0;
      intercept_all = 0.0;
      num_all = 0;
      slope_sel = 0.0;
      intercept_sel = 0.0;
      num_sel = 0;
      return;
    }
  // Sort list
  std::sort(deltalist.begin(), deltalist.end());
  scala::LinearFit LineAll;
  scala::LinearFit LineSel;
  NormalProbability NormProb;
  float w = 1.0; // unit weights

  for (int i=0;i<Nobs;i++)
    {
      // Expected delta
      // rank is 1->Nobs = i+1
      float DeltaExp = NormProb.ExpectedDelta(i+1, Nobs);
      LineAll.add(DeltaExp, deltalist[i], w);
      if (std::abs(DeltaExp) < dltlim)
        LineSel.add(DeltaExp, deltalist[i], w);
    }

  RPair FitAll = LineAll.result();
  slope_all = FitAll.first;
  intercept_all = FitAll.second;
  num_all = LineAll.Number();
  RPair FitSel = LineSel.result();
  slope_sel = FitSel.first;
  intercept_sel = FitSel.second;
  num_sel = LineSel.Number();

  Fitted = true;
}
//-------------------------------------------------------------
  //  = 0,0 if no data
float NormalProbAnal::Slope(const float& DltLim, const bool& All)
// Returns slopes for central part within dltlim
//   or for all data if All = true
//  = 0,0 if no data
{
  if (DltLim > 0.0) dltlim = DltLim;
  if (!Fitted) {Fit();}  // do the line fitting if not done
  if (All) {return slope_all;}
  else {return slope_sel;}
}
//-------------------------------------------------------------
float NormalProbAnal::Intercept(const float& DltLim, const bool& All)
// Returns intercept for central part within dltlim
//   or for all data if All = true
{
  if (DltLim > 0.0) dltlim = DltLim;
  if (!Fitted) {Fit();}  // do the line fitting if not done
  if (All) {return intercept_all;} else {return intercept_sel;}
}
//-------------------------------------------------------------
int NormalProbAnal::Number(const float& DltLim, const bool& All)
// Returns numbers for all data or central part within dltlim
//  = 0,0 if no data
{
  if (DltLim > 0.0) dltlim = DltLim;
  if (!Fitted) {Fit();}  // do the line fitting if not done
  if (All) {return num_all;} else {return num_sel;}
}
//-------------------------------------------------------------
void NormalProbAnal::Plot(NormalProbPlot& NormPlot,
                          const std::string& legend)
{
  int Nobs = deltalist.size();
  if (Nobs > 0)  {
    if (!Fitted) Fit();
    NormPlot.NewLine(legend);
    NormPlot.SetSample(deltalist.size());  // set sampling
    NormalProbability NormProb;
    for (size_t i=0;i<deltalist.size();i++) {
      float DeltaExp = NormProb.ExpectedDelta(i+1, Nobs);
      NormPlot.OutputPoint(DeltaExp, deltalist[i]);
    }
    NormPlot.EndLine();
  }
}
//-------------------------------------------------------------
//-------------------------------------------------------------
