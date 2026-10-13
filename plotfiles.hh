// plotfiles.hh

// Classes for writing plotfiles
//   for Grace (alias xmgr), and as self-contained gnuplot scripts (.plt)

#ifndef PLOTFILES_HEADER
#define PLOTFILES_HEADER

#include <string>
#include <vector>
#include "CCP4base.hh"

#include "range.hh"
#include "icering.hh"

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
class XMLplot
// XML plot for Pimple/Qloggraph
{
public:
  XMLplot(){}

  // Return true if plotting is turned on
  bool IsPlot() const {return (headerstring != "");}

  // title           title
  // xmin, xmax      range of coordinates
  // ymin, ymax      range of coordinates
  // xtick, ytick    tick intervals  (ignored)
  // ticklabel       true to label ticks  (ignored)
  // xlabel, ylabel  axis labels
  // xlegend, ylegend position of legends
  // WriteLegend     true to write legend
  void Header(const std::string& title,
	      const float& xmin, const float& xmax,
	      const float& ymin, const float& ymax,
	      const float& xtick, const float& ytick, const bool& ticklabel,
	      const std::string& xlabel, const std::string& ylabel,
	      const float& xlegend, const float& ylegend, 
	      const bool& WriteLegend);

  // start new line
  // legend    also used as ID for associated dataset
  // linestyle = 0 no line, else 1,2,3,4 = 'Solid','Dashed','Dash-dot','Dotted'
  // linesize  line width, default = 1
  // lcolour    <= 0 use default, else colour number
  // Symbol    !=0 draw symbols, if < 0 use default symbol,
  //            >0 use specified symbol
  void StartLine(const std::string& legend,
		 const int& linestyle,
		 const int& linesize=1,
		 const int& lcolour=-1,
		 const int& Symbol=-1 );

  // write point x,y
  // fw = field width, fd = #decimal
  void Point(const float& x, const float& y, const int& fw, const int& fd);
  void EndLine();

  // Draw a line
  // linestyle = 0 no line, else 1,2,3,4 = 'Solid','Dashed','Dash-dot','Dotted'
  // lcolour    <= 0 use default, else colour number
  void DrawLine(const float& xmin, const float& xmax,
		const float& ymin, const float& ymax,
		const int& fw, const int& fd,
		const int& linestyle=1,
		const int& linesize=1,
		const int& lcolour=-1);
  
  // Draw a circle
  void DrawCircle(const float& xcen, const float& ycen,
		  const float& radius,
		  const int& fw, const int& fd,
		  const int& linestyle=1,
		  const int& linesize=1,
		  const int& lcolour=-1,
		  const int& fillcolour=-1);

  void ClosePlot();

  // Return formatted XML
  std::string format() const;

private:
  std::string titlestring;   // Table title
  std::string headerstring;  // for header stuff (<plot> etc)
  std::string datastring;    // for the data tables (<data>)
  std::vector<std::string> dataIDs; // ID string for each set of data/plotline
  int symbolSize;
  
  // return colour string
  std::string Colour(const int& lcolour) const;
  // return linestyle string
  std::string LineStyle(const int& linestyle) const;


};
//--------------------------------------------------------------
class XMGRACE
{
public:
  XMGRACE() : file(NULL) {}

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
	    const int& lcolour, const int& Symbol, const bool& Join);

  
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
class GNUPLOT
// Self-contained gnuplot script, written as a parallel to XMGRACE.
// Points are buffered, since a gnuplot script needs its data (as inline
// datablocks) before the plot command; the file is written in ClosePlot.
// Colour numbers are the xmgr ones, so both files look alike.
{
public:
  GNUPLOT() : diagonal(false), square(false), ticklabel(true) {}

  // Return true if plotting is turned on
  bool IsPlot() const {return (filename != "");}

  // FileName        output file, written in ClosePlot
  // title1, title2  titles
  // xmin, xmax      range of coordinates
  // ymin, ymax      range of coordinates
  // xtick, ytick    tick intervals
  // Ticklabel       true to label ticks
  // xlabel, ylabel  axis labels
  // WriteLegend     true to write legend
  // Square          true for equal x and y scales
  void Header(const std::string& FileName,
	      const std::string& title1, const std::string& title2,
	      const float& xmin, const float& xmax,
	      const float& ymin, const float& ymax,
	      const float& xtick, const float& ytick, const bool& Ticklabel,
	      const std::string& xlabel, const std::string& ylabel,
	      const bool& WriteLegend, const bool& Square);

  // start new line, arguments as XMGRACE::Line
  void Line(const std::string& legend,
	    const int& lcolour, const int& Symbol, const bool& Join);

  void Point(const std::string& fmt,
	     const float& x, const float& y);

  void EndLine() {}

  // Draw a circle centred on 0,0
  void Circle(const int& lcolour, const std::string& legend,
	      const float& Radius, const int& Npoints);

  // plot diagonal line and axes if diagonal true, then write the file
  void ClosePlot(const bool& Diagonal);

private:
  struct Series {
    std::string legend;
    int colour;
    int symbol;
    bool join;
    std::string data;
  };
  std::string filename;
  std::string title1, title2;
  float xmin, xmax, ymin, ymax, xtick, ytick;
  std::string xlabel, ylabel;
  bool writelegend, diagonal, square, ticklabel;
  std::vector<Series> series;
};
//--------------------------------------------------------------
class NormalProbPlot
{
public:
  NormalProbPlot() :  file(0), Nlines(0), xmgr(true), gnu(false) {}
  // Create and write header
  // FileName     xmgr file, "" for none
  // GnuFileName  gnuplot file, "" for none
  NormalProbPlot(const std::string& FileName,
		 const bool& WriteLegend,
		 const std::string& title1,
		 const std::string& title2,
		 const std::string& GnuFileName="");

  void init(const std::string& FileName,
	    const bool& WriteLegend,
	    const std::string& title1,
	    const std::string& title2,
	    const std::string& GnuFileName="");

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
  void ClosePlot();

  std::string formatXML() const {return xmlplot.format();}

private:
  FILE* file;
  int Nlines;
  float xmax, ymax;
  bool xmgr;  // true to write xmgr file
  bool gnu;   // true to write gnuplot file
  XMGRACE xmgrplot;
  GNUPLOT gnuplot;
  XMLplot xmlplot;
  PlotSample sample;  // sampling of points for plotting
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
  
  // returns XML; plotfile = NULL for no xmgr file, gnufile = "" for no gnuplot
  std::string Plot(FILE* plotfile, const std::string& gnufile="") const;
  
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
  RoguePlot() : file(NULL), xmgr(false), gnu(false) {}

  // Smax  maximum 4(sin theta/lambda)**2
  //                      = (d*max)**2 = 1/dmin**2
  // wavelength  
  // FileName     xmgr file, "" for none
  // GnuFileName  gnuplot file, "" for none
  RoguePlot(const std::string& FileName,
	    const std::string& Title, const float& Smax,
	    const float& wavelength, const scala::Rings& icerings,
	    const std::string& GnuFileName="");

  // Return true if plotting is turned on
  bool IsPlot() const {return xmgrplot.IsPlot() || gnuplot.IsPlot();}

  // Start points, symbols, no line
  void Start();
  // End plot
  void End();

  // Plot outlier point
  //  s       diffraction vector in diffratometer frame, 1/A units
  //  pclass   = 1 outlier, = 2 outlierAnom, = 3 Emax
  void PlotOutlier(const FVect3& s, const int& pclass=0);

  std::string formatXML() const {return xmlplot.format();}

private:
  FILE* file;
  bool xmgr;  // true to write xmgr file
  bool gnu;   // true to write gnuplot file
  XMGRACE xmgrplot;
  GNUPLOT gnuplot;
  XMLplot xmlplot;
  std::vector<float> x;
  std::vector<float> y;
  std::vector<int> sclass;  // = 1 outlier, = 2 outlierAnom, = 3 Emax
};
//--------------------------------------------------------------
#endif
