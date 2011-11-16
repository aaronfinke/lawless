// plotfiles.hh

// Classes for writing plotfiles
//   currently these are all for Grace (alias xmgr)

#ifndef PLOTFILES_HEADER
#define PLOTFILES_HEADER

#include <string>
#include <vector>
#include "CCP4base.hh"

#include "range.hh"

//--------------------------------------------------------------
class PlotSample
// Class to optionally sample points in a scatter or normal probability plot,
// to avoid the need to plot every possible point
// Based on an expected Gaussian distribution of point density
//   cf scala s/r lsample
{
public:
  PlotSample(){}

  // Npoints            total number of points to sample
  // maxPointDensity    approximate maximum number of points to plot in each bin
  // numberBins         number of sampling bins
  // binWidth           width of bins
  PlotSample(const int& Npoints, const int& maxPointDensity,
	     const int& numberBins, const float& binWidth);

  // Return true if we want this point
  bool Keep(const float& value) const;

private:
  bool lsample;
  float maxsample;
  float binwidth;
  std::vector<float> sampleBins;
};
//--------------------------------------------------------------
class XMGRACE
{
public:
  XMGRACE(){}

  // Return true if plotting is turned on
  bool IsPlot() const {return (file != NULL);}


  // File            output file
  // title1, title2  titles
  // xmin, xmax      range of coordinates
  // ymin, ymax      range of coordinates
  // xv, yv          size of viewport
  // xtick, ytick    tick intervals
  // ticklabel       true to label ticks
  // xlabel, ylabel  axis labels
  // xlegend, ylegend position of legends
  // WriteLegend     true to write legend
  void Header(FILE* file, const std::string& title1, const std::string& title2,
	      const float& xmin, const float& xmax,
	      const float& ymin, const float& ymax,
	      const float& xv, const float& yv,
	      const float& xtick, const float& ytick, const bool& ticklabel,
	      const std::string& xlabel, const std::string& ylabel,
	      const float& xlegend, const float& ylegend, 
	      const bool& WriteLegend);

  // start new line
  // Symbol    !=0 draw symbols, if < 0 fill symbol
  // Join = true to draw line as well as points
  void Line(const std::string& legend,
	    const int& lcolor, const int& Symbol, const bool& Join);

  
  void Point(const std::string& fmt,
	     const float& x, const float& y);

  void EndLine();

  // plot diagonal line if Nlines > 0
  // File not closed
  void ClosePlot(const float& range, const bool& diagonal);

private:
  FILE* file;
  float symbolSize;
  int objectNumber;  // track graph object number, from 0
};
//--------------------------------------------------------------
class NormalProbPlot
{
public:
  NormalProbPlot() :  file(0), Nlines(0) {}
  // Create and write header
  NormalProbPlot(const std::string& FileName,
		 const bool& WriteLegend,
		 const std::string& title1,
		 const std::string& title2);

  void init(const std::string& FileName,
	    const bool& WriteLegend,
	    const std::string& title1,
	    const std::string& title2);

  // Set up sampling for plotting central region
  // Because the plot can contain a ridiculously large number of points,
  // if there are a lot, only a sample are output. maxPointDensity is the maximum
  // point density, above this a random number is used to select points.
  // Setup up sampling array at intervals of dx from -maxSample to +maxSample
  // Outside this range, all points will be used
  void SetSample(const int& Npoints);

  // Start new line
  void NewLine(const std::string& legend);
  // Output one point (x,y)
  void OutputPoint(const float& x, const float& y);
  // Write end plot marker & increment line number
  void EndLine();
  // Write diagonal line & close
  void ClosePlot(const float& range=4.0);


private:
  FILE* file;
  int Nlines;
  XMGRACE xmgrplot;
  PlotSample sample;  // smapling of points for plotting
  float limit;
};
// ------------------------------------------------------------
class CorrelPlot {
  // correlation plot
public:
  CorrelPlot(){}
  
// title         graph title
// Pxd_title     dataset title
// NresBins      number of resolution bins
// UnitValue     RMS value, ie value to plot as 1.0
// Npoints       total number of points to plot (including those
//               omitted by sampling)
  CorrelPlot(const std::string& Title, const std::string& Pxd_title,
	     const int& NresBins,
	     const float& UnitValue,
	     const int& Npoints);
  void init(const std::string& Title, const std::string& Pxd_title,
	    const int& NresBins,
	    const float& UnitValue,
	    const int& Npoints);
  
  void SetEqualLimit(const bool& EqualLimit) {equallimit = EqualLimit;}

  void SetNbins(const int& NresBins);
  
  void AddPoint(const int& mres, const float& I1, const float& I2);
  
  void Plot(FILE* plotfile) const;
  
private:
  int nresbin;   // number of resolution bins
  std::string title;
  std::string pxd_title;
  std::vector<float> x;
  std::vector<float> y;
  std::vector<int> resbin;
  scala::Range valrange;
  float scale;  // scale factor for plot x,y values
  float limit;  // ranges for plotting
  PlotSample sample;  // sampling of points for plotting
  bool equallimit;  // true for equal +/- limits
};
//--------------------------------------------------------------
class RoguePlot {
  // Rogue plot
public:
  RoguePlot(){}

  // Smax  maximum 4(sin theta/lambda)**2
  //                      = (d*max)**2 = 1/dmin**2
  RoguePlot(const std::string& FileName,
	    const std::string& Title, const float& Smax);

  // Return true if plotting is turned on
  bool IsPlot() const {return xmgrplot.IsPlot();}

  // Start points, symbols, no line
  void Start();
  // End plot
  void End();

  // Plot outlier point
  //  d    d spacing
  //  s    diffraction vector in diffratometer frame, 1/A units
  void PlotOutlier(const float& d, const FVect3& s);

private:
  FILE* file;
  XMGRACE xmgrplot;
};
//--------------------------------------------------------------
#endif
