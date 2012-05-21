// plotfiles.cpp

#include "plotfiles.hh"
#include "scala_util.hh"
#include "string_util.hh"
#include "file_util.hh"
#include "jiffy.hh"
#include "icering.hh"

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;

//--------------------------------------------------------------
PlotSample::PlotSample(const int& Npoints, const int& maxPointDensity,
		       const int& numberBins, const float& binWidth)
// Npoints            total number of points to sample
// maxPointDensity    approximate maximum number of points to plot in each bin
// numberBins         number of sampling bins
// binWidth           width of bins
//   cf scala s/r lsample
{
  // 1/sqrt(2 pi)
  float rrt2pi = 1./sqrt(clipper::Util::twopi());
  if (Npoints < maxPointDensity) {
    //  Few points, no sampling needed
    lsample = false;
    return;
  }
  lsample = true;
  sampleBins.resize(numberBins+1);
  maxsample = numberBins * binWidth;
  binwidth = binWidth;
  //^
  //^  std::cout << "Npoints " << Npoints << "\n"; //^!

  for (int i=0;i<=numberBins;++i) {
    float x = float(i) * binWidth;
    // Expected number for this x = dP(x)/dx
    int nexpect = Nint(binWidth * float(Npoints) * rrt2pi *exp(-0.5*x*x));
    if (nexpect > maxPointDensity) {
      sampleBins[i] = float(maxPointDensity)/float(nexpect);
    } else {
      sampleBins[i] = 1.0;
    }
    //^    std::cout << "SampleBin " << i << " " << sampleBins[i]
    //^	      << " Nexpect " << nexpect << "\n"; //^!
  }
}
//--------------------------------------------------------------
bool PlotSample::Keep(const float& value) const
// Return true if we want this point
{
  if (lsample) {
    if (value < maxsample) {
      int k = int(value/binwidth);
      if (scala::FRandom(1.0) > sampleBins.at(k)) {
	//^	std::cout << "NotKept " << value << "\n";
	return false;
      }
    }
  }
  return true;
}
//--------------------------------------------------------------
void DrawCircle(XMGRACE& xmgrplot, const int& lcolour, const std::string& legend,
		const float& Radius, const int& Npoints)
// Draw circle centred on 0,0, radius Radius, sampled into Npoints
// xmgrace format
{
  const int maxpnt = 200;
  const int minpnt = 16;
  int npoints = Min(Max(minpnt, Npoints), maxpnt);
  double delang = 8.*atan(1.0)/double(npoints);
  xmgrplot.Line(legend, lcolour, 0, true);

  for (int i=0;i<=npoints;++i) {
    double a = delang * double(i);
    xmgrplot.Point("%8.4f %8.4f\n",
		   Radius*cos(a), Radius*sin(a));
  }
  xmgrplot.EndLine();
}
//--------------------------------------------------------------
//--------------------------------------------------------------
NormalProbPlot::NormalProbPlot(const std::string& FileName,
			       const bool& WriteLegend,
			       const std::string& title1,
			       const std::string& title2)
  : Nlines(0)
{
  init(FileName, WriteLegend, title1, title2);
}
//--------------------------------------------------------------
void NormalProbPlot::init(const std::string& FileName,
			  const bool& WriteLegend,
			  const std::string& title1,
			  const std::string& title2)
{
  file = OpenFile(FileName, true);  // open write file
  xmgrplot.Header(file, title1, title2, -4, +4, -5, +5,
		  0.6, 0.7,
		  1.0, 1.0, true,
		  "Delta(expected)", "Delta(observed)", 0.8, 0.0, true);
}
//--------------------------------------------------------------
void NormalProbPlot::SetSample(const int& Npoints)
// Set up sampling for plotting central region
// Because the plot can contain a ridiculously large number of points,
// if there are a lot, only a sample are output. maxPointDensity is the maximum
// point density, above this a random number is used to select points.
// Setup up sampling array at intervals of dx from -maxSample to +maxSample
// Outside this range, all points will be used
{
  limit = 5.0;  // range = -limit -> +limit on x & y
  // Set up sampling
  int maxPointDensity = 200;
  float maxSample = 3.0;   // sample values less than this
  int numberBins = 10;
  float binWidth = maxSample/float(numberBins);
  sample = PlotSample(Npoints, maxPointDensity, numberBins, binWidth);
}
//--------------------------------------------------------------
void NormalProbPlot::NewLine(const std::string& legend)
// Start new line, set style parameters
// Colour is incremented through line number (from 0)
{
  //  Colours: (maximum 15 in xmgr) 6 & 7 are not very good
  int lcol[] = {1,2,3,4,5,8,9,10,11,12,13,14,15,6,7};
  std::vector<int> lcolor(lcol, lcol+15);
  int maxlin = lcolor.size();
  int icol = Nlines%maxlin;  // wrap round after maxlin
  xmgrplot.Line(legend, lcolor[icol], +2, false);
}

//--------------------------------------------------------------
void NormalProbPlot::OutputPoint(const float& x, const float& y)
{
  if (std::abs(x) < limit && std::abs(y) < limit) {
    if (sample.Keep(Max(std::abs(x), std::abs(y)))) {
      xmgrplot.Point("%10.4f %9.4f\n", x, y);
    }
  }
}
//--------------------------------------------------------------
void NormalProbPlot::EndLine()
{
  xmgrplot.EndLine();
  Nlines++;  // increment line number
}
//--------------------------------------------------------------
void NormalProbPlot::ClosePlot(const float& range)
// Draw diagonal line
{
  xmgrplot.ClosePlot(range, true);
}
  // ------------------------------------------------------------
  // ------------------------------------------------------------
CorrelPlot::CorrelPlot(const std::string& Title, const std::string& Pxd_title,
		       const int& NresBins,
		       const float& UnitValue,
		       const int& Npoints)
// title         graph title
// Pxd_title     dataset title
// NresBins      number of resolution bins
// UnitValue     RMS value, ie value to plot as 1.0
// Npoints       total number of points to plot (including those
//               omitted by sampling)
{
  init(Title, Pxd_title, NresBins, UnitValue, Npoints);
}
// ------------------------------------------------------------
void CorrelPlot::init(const std::string& Title, const std::string& Pxd_title,
		      const int& NresBins,
		      const float& UnitValue,
		      const int& Npoints)
// title         graph title
// Pxd_title     dataset title
// NresBins      number of resolution bins
// UnitValue     RMS value, ie value to plot as 1.0
// Npoints       total number of points to plot (including those
//               omitted by sampling)
{
  title = Title;
  pxd_title = Pxd_title;
  SetNbins(NresBins);
  scale = 1.0;
  if (UnitValue != 0.0) {scale = 1.0/UnitValue;}
  limit = 10.0;  // range = -limit -> +limit on x & y
  // Set up sampling
  int maxPointDensity = 200;
  float maxSample = 3.0;   // sample values less than this
  int numberBins = 10;
  float binWidth = maxSample/float(numberBins);
  sample = PlotSample(Npoints, maxPointDensity, numberBins, binWidth);
  equallimit = true;  // equal +/- by default
}
// ------------------------------------------------------------
void CorrelPlot::SetNbins(const int& NresBins)
{
  nresbin = NresBins;
}
// ------------------------------------------------------------
void CorrelPlot::AddPoint(const int& mres, const float& I1, const float& I2)
{
  const float cossin45 = 0.707106781;   // cos 45 = sin 45
  if (std::abs(I1) < limit/scale && std::abs(I2) < limit/scale) {
    float xx = I1*scale;
    float yy = I2*scale;
    float samplevalue;
    samplevalue = Max(std::abs(xx), std::abs(yy));
    if (sample.Keep(samplevalue)) {
      x.push_back(xx);
      y.push_back(yy);
      resbin.push_back(mres);
      valrange.update(xx);
      valrange.update(yy);   // accumulate range of scaled values
    }
  }
}
// ------------------------------------------------------------
void CorrelPlot::Plot(FILE* plotfile) const
{
  if (x.size() <= 0) return;
  XMGRACE xmgr;
  // Limits on axes
  float lowlimit = valrange.min();
  if (equallimit) {
    lowlimit = -limit;
  }
  xmgr.Header(plotfile, title, pxd_title,
	      lowlimit, limit, lowlimit, limit,
	      0.6, 0.7,
	      2.0, 2.0, true,
	      "", "", 0.0, 0.0, false);
  // Split into nres parts by resolution bin
  int nres = 2;
  //  Colours: (maximum 15 in xmgr) 6 & 7 are not very good
  //  Not many used here, only nres
  int lcol[] = {1,2,3,4,5,8,9,10,11,12,13,14,15,6,7};
  std::vector<int> lcolour(lcol, lcol+15);
  int iline = 0;
  for (int iset=0;iset<nres;++iset) {  // loop nres sets of ranges
    std::string legend = "ResRange "+phaser_io::itos(iset+1);
    if (nres == 2) {
      legend = "Low resolution";
      if (iset > 0) {
	legend = "High resolution";
      }
    }
    xmgr.Line(legend, lcolour[iset], +3, false);
    int iresmin = iset*(nresbin/nres);
    int iresmax = (iset+1)*(nresbin/nres);
    if (iset == nres-1) {iresmax = nresbin;}  // last set
    int np = 0;
    for (size_t i=0;i<x.size();++i) {
      if (resbin[i] >= iresmin && resbin[i] < iresmax) {
	xmgr.Point("%8.3f %8.3f\n", x[i], y[i]);
	np++;
      }
    }
    if (np > 0) iline++;
    xmgr.EndLine();
  }
  xmgr.ClosePlot(limit, true);
}
// ------------------------------------------------------------
// ------------------------------------------------------------
double RingRadius(const double& dstar)
// dstar in rlu
// Radius of ring = tan(2theta) = d* cos(theta)/cos(2theta)
// d* = lambda/d
{
  double theta = asin(0.5*dstar);
  return dstar * cos(theta)/cos(2.*theta);
}
// ------------------------------------------------------------
RoguePlot::RoguePlot(const std::string& FileName,
		     const std::string& Title, const float& Smax,
		      const float& wavelength)
// Smax  maximum 4(sin theta/lambda)**2
//                      = (d*max)**2 = 1/dmin**2
{
  file = OpenFile(FileName, true);  // open write file
  
  double dstar = wavelength*sqrt(Smax);   // rlu
  float radius = RingRadius(dstar);
  float legx = radius*0.5;
  float legy = -radius*0.9;
  xmgrplot.Header(file, "Outliers on detector (horizontal rotation axis)",
		  Title, -radius, +radius, -radius, +radius,
		  0.6, 0.75,
		  radius, radius, false,
		  "", "", legx, legy, true);
  
  // Resolution ring
  int npoint = 96;  // sampling of circle
  int lcol = 1;     // black
  DrawCircle(xmgrplot, lcol, "", radius, npoint);
  
  // Draw ice rings
  scala::Rings icerings;
  icerings.DefaultIceRings();  // set default ice rings
  lcol = 2;     // probably red
  bool first = true;
  for (int ir=0;ir<icerings.Nrings();++ir) {
    float rad = RingRadius(wavelength * icerings.Dstar(ir));
    if (rad < radius) {
      std::string label = "";
      if (first) {
	label = "ice rings";
	first = false;
      }
      DrawCircle(xmgrplot, lcol, label, rad, npoint);
    }
  }
  // Axis lines
  lcol = 1;     // black
  xmgrplot.Line("", lcol, 0, true); // X axis
  xmgrplot.Point("%8.4f %8.4f\n", -radius, 0.0);
  xmgrplot.Point("%8.4f %8.4f\n", +radius, 0.0);
  xmgrplot.EndLine();
  xmgrplot.Line("", lcol, 0, true); // Y axis
  xmgrplot.Point("%8.4f %8.4f\n", 0.0, -radius);
  xmgrplot.Point("%8.4f %8.4f\n", 0.0, +radius);
  xmgrplot.EndLine();
}
// ------------------------------------------------------------
void RoguePlot::Start()
{
  // Start points, symbols, no line
  int lcol = 1;
  xmgrplot.Line("", lcol, -1, false);
}
// ------------------------------------------------------------
void RoguePlot::End()
{
  xmgrplot.EndLine();
  xmgrplot.ClosePlot(0.0, false);
}
// ------------------------------------------------------------
void RoguePlot::PlotOutlier(const FVect3& s)
// Plot outlier point
//  s    diffraction vector in diffractometer frame, rlu
{
  double dstar = sqrt(s*s);
  double theta = asin(0.5*dstar);
  float sc = 1.0/cos(2.*theta);

  float ydn = sc * s[1];  // y
  float zdn = sc * s[2];  // z
  xmgrplot.Point("%10.4f %10.4f\n", zdn, ydn);
}
//--------------------------------------------------------------
void XMGRACE::Header(FILE* File, const std::string& title1, const std::string& title2,
		     const float& xmin, const float& xmax,
		     const float& ymin, const float& ymax,
		     const float& xv, const float& yv,
		     const float& xtick, const float& ytick, const bool& ticklabel,
		     const std::string& xlabel, const std::string& ylabel,
		     const float& xlegend, const float& ylegend, 
		     const bool& WriteLegend)
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
{
  objectNumber = 0;  // count objects
  file = File;
  // Calculate symbol size
  float r = Max(xmax-xmin, ymax-ymin);
  symbolSize = 0.5/r;

  fprintf(file, "@version 40102\n");
  if (title1.size() > 0)
    {fprintf(file,"@ title \"%s\"\n", title1.c_str());}

  if (title2.size() > 0)
    {fprintf(file, "@ subtitle \"%s\"\n", title2.c_str());}

  fprintf(file, "@with G0\n");
  fprintf(file, "@  world xmin %9.4f\n",xmin);
  fprintf(file, "@  world xmax %9.4f\n",xmax);
  fprintf(file, "@  world ymin %9.4f\n",ymin);
  fprintf(file, "@  world ymax %9.4f\n",ymax);
  // Viewport
  float xvmin = 0.05;
  float yvmin = 0.15;
  fprintf(file, "@  view xmin %9.4f\n", xvmin);
  fprintf(file, "@  view xmax %9.4f\n", xvmin+xv);
  fprintf(file, "@  view ymin %9.4f\n", yvmin);
  fprintf(file, "@  view ymax %9.4f\n", yvmin+yv);
  fprintf(file, "@  xaxis  tick on\n");
  fprintf(file, "@  xaxis  tick major %9.4f\n", xtick);
  if (xlabel.size() > 0) {
    fprintf(file, "@  xaxis  label \"%s\"\n", xlabel.c_str());
  }
  fprintf(file, "@  xaxis  tick major grid on\n");
  fprintf(file, "@  xaxis  tick minor ticks 0\n");

  fprintf(file, "@  yaxis  tick on\n");
  fprintf(file, "@  yaxis  tick major %9.4f\n", ytick);
  if (ylabel.size() > 0) {
    fprintf(file, "@  yaxis  label \"%s\"\n", ylabel.c_str());
  }
  fprintf(file, "@  yaxis  tick major grid on\n");
  fprintf(file, "@  yaxis  tick minor ticks 0\n");
  if (!ticklabel) {
    fprintf(file, "@  xaxis  ticklabel off\n");
    fprintf(file, "@  yaxis  ticklabel off\n");
  }
  fprintf(file, "@  legend loctype world\n");
  fprintf(file, "@  legend x1 %9.3f\n", xlegend);
  fprintf(file, "@  legend y1 %9.3f\n", ylegend);
  fprintf(file, "@  legend color 1\n");
  fprintf(file, "@  legend vgap 2\n");
  fprintf(file, "@  legend hgap 1\n");
  fprintf(file, "@  legend length 1\n");
  fprintf(file, "@  legend font 4\n");
  fprintf(file, "@  legend char size 0.8\n");

  if (WriteLegend)
    {fprintf(file, "@  legend on\n");}
  else
    {fprintf(file, "@  legend off\n");}
  fprintf(file, "@G0 on\n");
}
//--------------------------------------------------------------
void XMGRACE::Line(const std::string& legend,
		   const int& lcolor,
		   const int& Symbol, const bool& Join)
// Symbol    !=0 draw symbols, if < 0 fill symbol, >1 scale symbol size
// Join = true to draw line as well as points
{
  std::string set = "s"+StringUtil::Strip(clipper::String(objectNumber));  // "Sn"
  int fill = 0;
  if (Symbol < 0) fill = 1;

  fprintf(file, "@  %s type xy\n", set.c_str());
  if (Symbol != 0) {
    fprintf(file, "@  %s  symbol size %9.4f\n", set.c_str(), Symbol*symbolSize);
    fprintf(file, "@  %s  symbol fill %2d\n", set.c_str(), fill);
    // symbol
    fprintf(file, "@  %s  symbol 2\n", set.c_str());
    // Colour
    fprintf(file, "@  %s  symbol color %3d\n", set.c_str(), lcolor);
  }
  if (Join) {
    fprintf(file, "@  %s  linestyle 1\n", set.c_str());
  } else {
    fprintf(file, "@  %s  linestyle 0\n", set.c_str());
  }
  // Colour
  fprintf(file, "@  %s  color %3d\n", set.c_str(), lcolor);
  // Legend
  if (legend.size() > 0) {
    fprintf(file, "@  legend string %3d \"%s\"\n", objectNumber,
	    legend.c_str());
  }
  objectNumber++;
}
//--------------------------------------------------------------
void XMGRACE::Point(const std::string& fmt,
		    const float& x, const float& y)
{
  fprintf(file,fmt.c_str(), x, y);
}
//--------------------------------------------------------------
void XMGRACE::EndLine()
{ 
  fprintf(file,"&\n");
}
//--------------------------------------------------------------
void XMGRACE::ClosePlot(const float& range,
			const bool& diagonal)
// plot diagonal line if diagonal true
{
  if (diagonal) {
    std::string set = "s"+StringUtil::Strip(clipper::String(objectNumber));
    fprintf(file, "@ %s color 1\n", set.c_str());
    fprintf(file, "@ %s linestyle 1\n", set.c_str());
    fprintf(file,"%10.4f %9.4f\n", -range, -range);
    fprintf(file,"%10.4f %9.4f\n", range, range);
    EndLine();
    objectNumber++;
  }
}
//--------------------------------------------------------------

