
// fracdevanal.cpp

#include "fracdevanal.hh"
using phaser_io::LOGFILE;

namespace scala {
//--------------------------------------------------------------
FracDeviationAnalysis::FracDeviationAnalysis(const IntensityBin& Irange)
{
  irange = Irange;  // copy intensity bin object
  Nbin = irange.NumberBins();  // number of bins
  // resize all arrays
  deviations = std::vector<MeanSD>(Nbin); // zeroes sums too
  sum_delta2 = 0.0;
  sum_lnS = 0.0;
  sum_w = 0.0;
  n_resid = 0;
}
//--------------------------------------------------------------
  void FracDeviationAnalysis::clear()
  {
    irange.Clear();
    total_deviation.clear();
    for (int i=0;i<Nbin;i++)
      {
	deviations[i].clear();
      }
    sum_delta2 = 0.0;
    sum_lnS = 0.0;
    sum_w = 0.0;
    n_resid = 0;
  }
//--------------------------------------------------------------
  void FracDeviationAnalysis::AddDeviation(const float& I, const float& D,
  const float& S, const float& w)
{
  if (S > 0.0001)
    {
      int Ibin = irange.bin(I);
      double delta = D/double(S);
      deviations[Ibin].Add(delta);
      total_deviation.Add(delta);
      sum_delta2 += w*(delta)*(delta);
      sum_lnS += w*log(S);
      sum_w +=w;
      n_resid++;
    }
}
//--------------------------------------------------------------
  void FracDeviationAnalysis::AddDeviation(const float& I, const float& delta)
  {
    int Ibin = irange.bin(I);
    deviations[Ibin].Add(delta);
    total_deviation.Add(delta);
    n_resid++;
}
//--------------------------------------------------------------
  float FracDeviationAnalysis::Slope() const
  // Slope of SD(deviation) v. Intensity, from bins
  {
    LinearFit line;
    float w = 1.0;
    for (int i=0;i<Nbin;i++)
      {
	if (irange.Count(i) > 0)
	  {
	    float x = irange.middle(i);
	    float y = deviations[i].SD();
	    line.add(x,y,w);
	  }
      }
    return line.result().first;
  }
//--------------------------------------------------------------
// Print results
void FracDeviationAnalysis::Print(const std::string& label, phaser_io::Output& output) const
{
  output.logTabPrintf(0,LOGFILE,
      "\n\nTotal SD of fractional deviation = %12.4g for %d observations\n\n",
		      total_deviation.SD(), total_deviation.Count());

  output.logTab(0,LOGFILE,"\n$TABLE: "+label+" standard deviation v. Intensity:\n");
  output.logTab(0,LOGFILE,"$GRAPHS: Sigma(scatter/SD) :N:4,7: $$\n");
  output.logTab(0,LOGFILE,
	" Range    Imin     Imax     Imean      Number    Mean   Sigma  $$ $$\n");
  
  
  for (int i=0;i<Nbin;i++)
    {
      RPair bounds = irange.bounds(i);

      output.logTabPrintf(0,LOGFILE,
			  "%5d%9.0f%9.0f  %9.0f  %9d%8.2f%8.2f\n",
			  i+1, bounds.first, bounds.second,
			  irange.middle(i), irange.Count(i),
			  deviations[i].Mean(), deviations[i].SD());
    }

  output.logTab(0,LOGFILE,"$$\n\n");

  output.logTabPrintf(0,LOGFILE,"Slope x 10000: %8.3f\n", 10000.*Slope());
}
//--------------------------------------------------------------
//--------------------------------------------------------------
}
