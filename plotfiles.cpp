// plotfiles.cpp

#include "plotfiles.hh"
#include "scala_util.hh"
#include "string_util.hh"
#include "file_util.hh"
#include "jiffy.hh"
#include "icering.hh"
#include "report_errors.hh"
#include "version.hh"

// Clipper
#include <clipper/clipper.h>
//using clipper::Message;
//using clipper::Message_fatal;


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
    //^       << " Nexpect " << nexpect << "\n"; //^!
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
        //^     std::cout << "NotKept " << value << "\n";
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
                               const std::string& title2,
                               const std::string& GnuFileName)
  : Nlines(0)
{
  init(FileName, WriteLegend, title1, title2, GnuFileName);
}
//--------------------------------------------------------------
void NormalProbPlot::init(const std::string& FileName,
                          const bool& WriteLegend,
                          const std::string& title1,
                          const std::string& title2,
                          const std::string& GnuFileName)
{
  xmgr = true;
  if (FileName == "") {
    // No XMGR output
    xmgr = false;
  }
  gnu = (GnuFileName != "");
  xmax = 4.0;
  ymax = 5.0;

  // Initialise file for XMgrace
  if (xmgr) {
    file = OpenFile(FileName, true);  // open write file for xmgrace
    xmgrplot.Header(file, title1, title2, -xmax, xmax, -ymax, ymax,
                    0.6, 0.7,
                    1.0, 1.0, true,
                    "Delta(expected)", "Delta(observed)", 0.8, 0.0, true);
  }
  if (gnu) {
    gnuplot.Header(GnuFileName, title1, title2, -xmax, xmax, -ymax, ymax,
                   1.0, 1.0, true,
                   "Delta(expected)", "Delta(observed)", true, false);
  }
  // Always initialise XML output
  xmlplot.Header(title1, -xmax, xmax, -ymax, ymax,
                 1.0, 1.0, true,
                 "Delta(expected)", "Delta(observed)", 0.6, 0.0, true);
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
  std::vector<int> lcolour(lcol, lcol+15);
  int maxlin = lcolour.size();
  int icol = Nlines%maxlin;  // wrap round after maxlin
  if (xmgr) {xmgrplot.Line(legend, lcolour[icol], +2, false);}
  if (gnu) {gnuplot.Line(legend, lcolour[icol], +2, false);}
  // start XML line, with a legend, not connected
  xmlplot.StartLine(legend, 0, 1, lcolour[icol]);
}

//--------------------------------------------------------------
void NormalProbPlot::OutputPoint(const float& x, const float& y)
{
  if (std::abs(x) < limit && std::abs(y) < limit) {
    if (sample.Keep(Max(std::abs(x), std::abs(y)))) {
      if (xmgr) {xmgrplot.Point("%10.4f %9.4f\n", x, y);}
      if (gnu) {gnuplot.Point("%10.4f %9.4f\n", x, y);}
      xmlplot.Point(x, y, 9, 4);
    }
  }
}
//--------------------------------------------------------------
void NormalProbPlot::EndLine()
{
  if (xmgr) {xmgrplot.EndLine();}
  if (gnu) {gnuplot.EndLine();}
  xmlplot.EndLine();
  Nlines++;  // increment line number
}
//--------------------------------------------------------------
void NormalProbPlot::ClosePlot()
// Draw diagonal line
{
  float range = xmax;
  if (xmgr) {xmgrplot.ClosePlot(range, true);}
  if (gnu) {gnuplot.ClosePlot(true);}
  // XML
  xmlplot.DrawLine(-range, +range, -range, +range, 9,4,
                   1,1,7); // black
  xmlplot.DrawLine(-range, +range, 0.0, 0.0, 9,4,
                   1,1,7); // black
  xmlplot.DrawLine(0.0, 0.0, -ymax, +ymax, 9,4,
                   1,1,7); // black

  xmlplot.ClosePlot();
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
std::string CorrelPlot::Plot(FILE* plotfile, const std::string& gnufile) const
// returns XML representation
// If plotfile = NULL, don't write xmgr file
// If gnufile = "", don't write gnuplot file
{
  if (x.size() <= 0) return "";
  XMLplot xmlplot;
  XMGRACE xmgr;
  GNUPLOT gnuplot;
  bool bxmgr = (plotfile != NULL);
  bool bgnu = (gnufile != "");

  // Limits on axes
  float lowlimit = valrange.min();
  if (equallimit) {
    lowlimit = -limit;
  }
  if (bxmgr) {
    xmgr.Header(plotfile, title, pxd_title,
                lowlimit, limit, lowlimit, limit,
                0.6, 0.7,
                2.0, 2.0, true,
                "", "", 0.6, 0.0, true);
  }
  if (bgnu) {
    gnuplot.Header(gnufile, title, pxd_title,
                   lowlimit, limit, lowlimit, limit,
                   2.0, 2.0, true,
                   "", "", true, equallimit);
  }
  xmlplot.Header(title,
              lowlimit, limit, lowlimit, limit,
              2.0, 2.0, true,
              "", "", 0.6, 0.0, true);
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
    if (bxmgr) {xmgr.Line(legend, lcolour[iset], +3, false);}
    if (bgnu) {gnuplot.Line(legend, lcolour[iset], +3, false);}
    xmlplot.StartLine(legend, 0, 1, lcolour[iset]);
    int iresmin = iset*(nresbin/nres);
    int iresmax = (iset+1)*(nresbin/nres);
    if (iset == nres-1) {iresmax = nresbin;}  // last set
    int np = 0;
    for (size_t i=0;i<x.size();++i) {
      if (resbin[i] >= iresmin && resbin[i] < iresmax) {
        if (bxmgr) {xmgr.Point("%8.3f %8.3f\n", x[i], y[i]);}
        if (bgnu) {gnuplot.Point("%8.3f %8.3f\n", x[i], y[i]);}
        xmlplot.Point( x[i], y[i], 8, 3);
        np++;
      }
    }
    if (np > 0) iline++;
    if (bxmgr) {xmgr.EndLine();}
    if (bgnu) {gnuplot.EndLine();}
    xmlplot.EndLine();
  }
  if (bxmgr) {xmgr.ClosePlot(limit, true);}
  if (bgnu) {gnuplot.ClosePlot(true);}

  // XML
  xmlplot.DrawLine(-limit, +limit, -limit, +limit, 9,4,
                   1,1,7); // black
  xmlplot.DrawLine(-limit, +limit, 0.0, 0.0, 9,4,
                   1,1,7); // black
  xmlplot.DrawLine(0.0, 0.0, -limit, +limit, 9,4,
                   1,1,7); // black
    xmlplot.ClosePlot();
  return xmlplot.format();
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
		     const float& wavelength, const scala::Rings& icerings,
                     const std::string& GnuFileName)
// Smax  maximum 4(sin theta/lambda)**2
//                      = (d*max)**2 = 1/dmin**2
// FileName     xmgr file, "" for none
// GnuFileName  gnuplot file, "" for none
{
  xmgr = true;
  if (FileName == "") {xmgr = false;}
  gnu = (GnuFileName != "");
  if (xmgr) {file = OpenFile(FileName, true);}  // open write file

  x.clear(); y.clear(); sclass.clear();  // clear arrays

  double dstar = wavelength*sqrt(Smax);   // rlu
  float radius = RingRadius(dstar);
  float legx = 0.0;
  float legy = 0.0;
  if (xmgr) {
    xmgrplot.Header(file, "Outliers on detector (horizontal rotation axis)",
                    Title, -radius, +radius, -radius, +radius,
                    0.6, 0.75,
                    radius, radius, false,
                    "", "", legx, legy, true);
  }
  if (gnu) {
    gnuplot.Header(GnuFileName, "Outliers on detector (horizontal rotation axis)",
                   Title, -radius, +radius, -radius, +radius,
                   radius, radius, false,
                   "", "", true, true);
  }
  xmlplot.Header("Outliers on detector (horizontal rotation axis)",
                 -radius, +radius, -radius, +radius,
                 radius, radius, false,
                 "", "", legx, legy, true);

  // Resolution ring
  int npoint = 96;  // sampling of circle
  int lcol = 1;     // black
  if (xmgr) {DrawCircle(xmgrplot, lcol, "", radius, npoint);}
  if (gnu) {gnuplot.Circle(lcol, "", radius, npoint);}
  xmlplot.DrawCircle(0.0, 0.0, radius, 8, 4, 1, 1, lcol);

  // Draw ice rings
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
      if (xmgr) {DrawCircle(xmgrplot, lcol, label, rad, npoint);}
      if (gnu) {gnuplot.Circle(lcol, label, rad, npoint);}
      xmlplot.DrawCircle(0.0, 0.0, rad, 8, 4, 1, 1, lcol);
    }
  }
  // Axis lines
  lcol = 1;     // black
  if (xmgr) {
    xmgrplot.Line("", lcol, 0, true); // X axis
    xmgrplot.Point("%8.4f %8.4f\n", -radius, 0.0);
    xmgrplot.Point("%8.4f %8.4f\n", +radius, 0.0);
    xmgrplot.EndLine();
    xmgrplot.Line("", lcol, 0, true); // Y axis
    xmgrplot.Point("%8.4f %8.4f\n", 0.0, -radius);
    xmgrplot.Point("%8.4f %8.4f\n", 0.0, +radius);
    xmgrplot.EndLine();
  }
  if (gnu) {
    gnuplot.Line("", lcol, 0, true); // X axis
    gnuplot.Point("%8.4f %8.4f\n", -radius, 0.0);
    gnuplot.Point("%8.4f %8.4f\n", +radius, 0.0);
    gnuplot.EndLine();
    gnuplot.Line("", lcol, 0, true); // Y axis
    gnuplot.Point("%8.4f %8.4f\n", 0.0, -radius);
    gnuplot.Point("%8.4f %8.4f\n", 0.0, +radius);
    gnuplot.EndLine();
  }
  xmlplot.DrawLine(-radius, +radius, 0.0, 0.0, 8,4);
  xmlplot.DrawLine(0.0, 0.0, -radius, +radius, 8,4);
}
// ------------------------------------------------------------
void RoguePlot::Start()
{
  // Nothing done here, all plotting in function "End"

}
// ------------------------------------------------------------
void RoguePlot::End()
{
  ASSERT ((x.size() == y.size()) &&  (y.size() == sclass.size()));
  // Do plots for each class
  // Do we have members of each class?
  std::vector<int> numinclass(3,0);
  for (size_t k=0; k<sclass.size(); k++) {
    if (sclass[k] > 0 && sclass[k] <= 3) {
      numinclass[sclass[k]-1]++;
    } else {
      // shouldn't happen
      ReportErrors::printFatalError("RoguePlot::End invalid class");
    }
  }

  // Loop classes 1,2,3
  std::string legend;
  for (int j=1;j<=3;++j){
    if (numinclass[j-1] > 0) {
      // Start points, symbols, no line
      int lcol = j;  // colour
      if (j == 1) {
        legend = "Outlier";
      } else if (j == 2) {
        legend = "AnomOutlier";
      } else if (j == 3) {
        legend = "Emax";
      }
      if (xmgr) {xmgrplot.Line(legend, lcol, -1, false);}
      if (gnu) {gnuplot.Line(legend, lcol, -1, false);}
      xmlplot.StartLine(legend, 0, 0, lcol);
      for (size_t k=0; k<sclass.size(); k++) {
        if (sclass[k] == j) {
          // a point in this class
          if (xmgr) {xmgrplot.Point("%10.4f %10.4f\n", x[k], y[k]);}
          if (gnu) {gnuplot.Point("%10.4f %10.4f\n", x[k], y[k]);}
          xmlplot.Point(x[k], y[k], 10, 4);
        }
      }
      if (xmgr) {xmgrplot.EndLine();}
      if (gnu) {gnuplot.EndLine();}
      xmlplot.EndLine();
    }
  }
  if (xmgr) {xmgrplot.ClosePlot(0.0, false);}
  if (gnu) {gnuplot.ClosePlot(false);}
  xmlplot.ClosePlot();
}
// ------------------------------------------------------------
void RoguePlot::PlotOutlier(const FVect3& s, const int& pclass)
// Plot outlier point
//  s    diffraction vector in diffractometer frame, rlu
//  pclass   = 1 outlier, = 2 outlierAnom, = 3 Emax
{
  double dstar = sqrt(s*s);
  double theta = asin(0.5*dstar);
  float sc = 1.0/cos(2.*theta);

  float ydn = sc * s[1];  // y
  float zdn = sc * s[2];  // z
  x.push_back(zdn);
  y.push_back(ydn);
  sclass.push_back(pclass);
}
//--------------------------------------------------------------
void XMLplot::Header(const std::string& title,
                     const float& xmin, const float& xmax,
                     const float& ymin, const float& ymax,
                     const float& xtick, const float& ytick, const bool& ticklabel,
                     const std::string& xlabel, const std::string& ylabel,
                     const float& xlegend, const float& ylegend,
                     const bool& WriteLegend)
// title           title
// xmin, xmax      range of coordinates
// ymin, ymax      range of coordinates
// xtick, ytick    tick intervals (not used)
// ticklabel       true to label ticks (not used)
// xlabel, ylabel  axis labels
// xlegend, ylegend position of legends
// WriteLegend     true to write legend
{
  datastring = "<headers separator=\" \">\n";
  datastring += " X Y\n";
  datastring += "</headers>\n";

  dataIDs.clear();

  // Title
  titlestring = "<CCP4Table title=\""+title+"\">\n";
  headerstring += "<plot>\n";  // start plot
  headerstring += "<title>"+title+"</title>\n";

  symbolSize = 2;

  // Axes
  headerstring += "<xrange min=\""+StringUtil::ftos(xmin,8,4)+
    "\" max=\""+StringUtil::ftos(xmax,8,4)+"\"/>\n";
  headerstring += "<yrange min=\""+StringUtil::ftos(ymin,8,4)+
    "\" max=\""+StringUtil::ftos(ymax,8,4)+"\"/>\n";
  headerstring += "<fixaspectratio>true</fixaspectratio>\n";

  if (xlabel != "") {headerstring += StringUtil::MakeXMLtag("xlabel", xlabel)+"\n";}
  if (ylabel != "") {headerstring += StringUtil::MakeXMLtag("ylabel", ylabel)+"\n";}

  // legends
  if (WriteLegend) {
    headerstring += StringUtil::MakeXMLtag("showlegend","true")+"\n";
    headerstring += "<legendposition x=\""+StringUtil::ftos(xlegend,8,3)+
      "\" y=\""+StringUtil::ftos(ylegend,8,3)+"\"/>\n";
  } else {
    headerstring += StringUtil::MakeXMLtag("showlegend","false")+"\n";
  }
}
//--------------------------------------------------------------
void XMLplot::StartLine(const std::string& legend,
                        const int& linestyle,
                        const int& linesize,
                        const int& lcolour,
                        const int& Symbol)
// legend    also used as ID for associated dataset
// linestyle <= 0 no line, else 1,2,3,4 = 'Solid','Dashed','Dash-dot','Dotted'
// linesize  line width, default = 1
// lcolour    <= 0 use default, else colour number
// Symbol    !=0 draw symbols, if < 0 use default symbol,
//            >0 use specified symbol
{
  // If legend is blank don't add a dataid to plotline
  dataIDs.push_back(legend);
  if (legend == "") {
    headerstring += "<plotline xcol=\"1\" ycol=\"2\">\n";
  } else {
    headerstring += "<plotline dataid=\""+legend+
    "\" xcol=\"1\" ycol=\"2\">\n";
    headerstring += StringUtil::MakeXMLtag("label", legend)+"\n";
  }
  int symbolsize = symbolSize;
  if (Symbol == 0) { // no symbol
    symbolsize = 0;
  }
  if (Symbol > 0) { // use chosen symbol
    std::string symbols[] = {"o","s","d","^",">","<","p","h","*"};
    const int nsymbols = 9;  // number in list
    int index = (Symbol-1)%nsymbols;  // wrap around
    headerstring += StringUtil::MakeXMLtag("symbol", symbols[index])+"\n";
  }
  headerstring += StringUtil::MakeXMLtag("symbolsize",
                                        StringUtil::itos(symbolsize,2))+"\n";
  if (lcolour > 0) { // colour
    headerstring += StringUtil::MakeXMLtag("colour", Colour(lcolour))+"\n";
  }
  headerstring += StringUtil::MakeXMLtag("markeredgewidth","0.0")+"\n"; // no symbol border
  // linestyle = no line, else 1,2,3,4 = 'Solid','Dashed','Dash-dot','Dotted'
  headerstring += StringUtil::MakeXMLtag("linestyle",LineStyle(linestyle))+"\n";
  if (linesize > 1) {
    headerstring += StringUtil::MakeXMLtag("linesize", StringUtil::itos(linesize,3));
  }
  if (legend == "") {
    datastring += "<data>\n";
  } else {
    datastring += "<data id=\""+legend+"\">\n";
  }
  headerstring += "</plotline>\n";
}
//--------------------------------------------------------------
void XMLplot::Point(const float& x, const float& y, const int& fw, const int& fd)
{
  datastring += StringUtil::ftos(x, fw, fd) + " " + StringUtil::ftos(y, fw, fd)+"\n";
}
//--------------------------------------------------------------
void XMLplot::EndLine()
{
  datastring += "</data>\n";
}
//--------------------------------------------------------------
void XMLplot::ClosePlot()
{
    headerstring += "</plot>\n";
}
//--------------------------------------------------------------
// Draw a line
// linestyle = 0 no line, else 1,2,3,4 = 'Solid','Dashed','Dash-dot','Dotted'
// lcolour    <= 0 use default, else colour number
void XMLplot::DrawLine(const float& xmin, const float& xmax,
                       const float& ymin, const float& ymax,
                       const int& fw, const int& fd,
                       const int& linestyle,
                       const int& linesize,
                       const int& lcolour)
{
  headerstring += "<line x1=\"" + StringUtil::ftos(xmin,fw,fd)+"\"" +
    " x2=\"" + StringUtil::ftos(xmax,fw,fd)+"\"" +
    " y1=\"" + StringUtil::ftos(ymin,fw,fd)+"\"" +
    " y2=\"" + StringUtil::ftos(ymax,fw,fd)+"\"";
  headerstring += " linestyle=\""+LineStyle(linestyle)+"\"";
  if (linesize > 1) {
    headerstring += " linesize=\""+StringUtil::itos(linesize,3)+"\"";
  }
  if (lcolour > 0) {
    headerstring += " linecolour=\""+Colour(lcolour)+"\"";
  }
  headerstring += "/>\n";
}
//--------------------------------------------------------------
// Draw a circle
void XMLplot::DrawCircle(const float& xcen, const float& ycen,
                         const float& radius,
                         const int& fw, const int& fd,
                         const int& linestyle,
                         const int& linesize,
                         const int& lcolour,
                         const int& fillcolour)
{
  headerstring += "<circle xpos=\"" + StringUtil::ftos(xcen,fw,fd)+"\"" +
    " ypos=\"" + StringUtil::ftos(ycen,fw,fd)+"\"" +
    " radius=\"" + StringUtil::ftos(radius,fw,fd)+"\"";
  if (linestyle > 1) {
    headerstring += " linestyle=\""+LineStyle(linestyle)+"\"";
  }
  if (linesize > 1) {
    headerstring += " linesize=\""+StringUtil::itos(linesize,3)+"\"";
  }
  if (lcolour > 0) {
    headerstring += " linecolour=\""+Colour(lcolour)+"\"";
  }
  if (fillcolour > 0) {
    headerstring += " fillcolour=\""+Colour(fillcolour)+"\"";
  }
  headerstring += "/>\n";
}
//--------------------------------------------------------------
std::string XMLplot::Colour(const int& lcolour) const
// return colour string
{
  std::string s = "";
  if (lcolour > 0) { // colour
    std::string colours[] = {"red", "green", "blue", "yellow",
                             "magenta", "cyan", "black"};
    const int ncol = 7;
    int index = (lcolour-1)%ncol;
    s = colours[index];
  }
  return s;
}
//--------------------------------------------------------------
std::string XMLplot::LineStyle(const int& linestyle) const
// return linestyle string
{
  // linestyle = no line, else 1,2,3,4 = 'Solid','Dashed','Dash-dot','Dotted'
    std::string s = "";
  if (linestyle <= 0) {
    s = ".";
  } else {
    std::string styles[] = {"-","--","-.",":"};
    const int nstyles = 4;
    int index = (linestyle-1)%nstyles;
    s = styles[index];
  }
  return s;
}
//--------------------------------------------------------------
std::string XMLplot::format() const
// return formatted XML
{
  return titlestring + datastring + headerstring + "\n</CCP4Table>\n";
}
//--------------------------------------------------------------
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
  //  float r = Max(xmax-xmin, ymax-ymin);
  // symbolSize = 0.5/r;
  // Fixed symbol size
  symbolSize = 0.25;

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
                   const int& lcolour,
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
    fprintf(file, "@  %s  symbol color %3d\n", set.c_str(), lcolour);
  }
  if (Join) {
    fprintf(file, "@  %s  linestyle 1\n", set.c_str());
  } else {
    fprintf(file, "@  %s  linestyle 0\n", set.c_str());
  }
  // Colour
  fprintf(file, "@  %s  color %3d\n", set.c_str(), lcolour);
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
//--------------------------------------------------------------
//--------------------------------------------------------------
namespace {
  // gnuplot double-quoted string
  std::string GnuQuote(const std::string& str)
  {
    std::string q = "\"";
    for (size_t i=0;i<str.size();++i) {
      if (str[i] == '\n') {q += "\\n"; continue;}
      if (str[i] == '"' || str[i] == '\\') q += '\\';
      q += str[i];
    }
    return q + "\"";
  }
  // xmgr colour number to rgb
  std::string GnuColour(const int& lcolour)
  {
    static const char* colours[] = {
      "#ffffff", "#000000", "#e41a1c", "#00a000", "#1f4ed8", "#d4a017",
      "#8b4513", "#808080", "#8a2be2", "#00b4c8", "#d020d0", "#ff8c00",
      "#4b0082", "#800000", "#40c8b0", "#008b00"};
    int ic = (lcolour >= 0 && lcolour <= 15) ? lcolour : 1;
    return GnuQuote(colours[ic]);
  }
}
//--------------------------------------------------------------
void GNUPLOT::Header(const std::string& FileName,
                     const std::string& Title1, const std::string& Title2,
                     const float& Xmin, const float& Xmax,
                     const float& Ymin, const float& Ymax,
                     const float& Xtick, const float& Ytick,
                     const bool& Ticklabel,
                     const std::string& Xlabel, const std::string& Ylabel,
                     const bool& WriteLegend, const bool& Square)
{
  filename = FileName;
  title1 = Title1;  title2 = Title2;
  xmin = Xmin;  xmax = Xmax;  ymin = Ymin;  ymax = Ymax;
  xtick = Xtick;  ytick = Ytick;  ticklabel = Ticklabel;
  xlabel = Xlabel;  ylabel = Ylabel;
  writelegend = WriteLegend;
  square = Square;
  diagonal = false;
  series.clear();
}
//--------------------------------------------------------------
void GNUPLOT::Line(const std::string& legend,
                   const int& lcolour, const int& Symbol, const bool& Join)
// Symbol    !=0 draw symbols, if < 0 fill symbol, >1 scale symbol size
// Join = true to draw line as well as points
{
  Series sr;
  sr.legend = legend;
  sr.colour = lcolour;
  sr.symbol = Symbol;
  sr.join = Join;
  series.push_back(sr);
}
//--------------------------------------------------------------
void GNUPLOT::Point(const std::string& fmt,
                    const float& x, const float& y)
{
  if (series.empty()) return;
  char buf[64];
  snprintf(buf, sizeof(buf), fmt.c_str(), x, y);
  series.back().data += buf;
}
//--------------------------------------------------------------
void GNUPLOT::Circle(const int& lcolour, const std::string& legend,
                     const float& Radius, const int& Npoints)
// Draw circle centred on 0,0, radius Radius, sampled into Npoints
{
  const int maxpnt = 200;
  const int minpnt = 16;
  int npoints = Min(Max(minpnt, Npoints), maxpnt);
  double delang = 8.*atan(1.0)/double(npoints);
  Line(legend, lcolour, 0, true);
  for (int i=0;i<=npoints;++i) {
    double a = delang * double(i);
    Point("%8.4f %8.4f\n", Radius*cos(a), Radius*sin(a));
  }
  EndLine();
}
//--------------------------------------------------------------
void GNUPLOT::ClosePlot(const bool& Diagonal)
// plot diagonal line and axes if Diagonal true, then write the file
{
  if (filename == "") return;
  diagonal = Diagonal;

  std::string s;
  s += "# " + filename + " - " + title1 + "\n";
  s += "# Generated by " + PROGRAM_NAME + " " + PROGRAM_VERSION + "\n";
  if (title2 != "") {
    std::string t2 = title2;  // keep the comment on one line
    for (size_t i=0;i<t2.size();++i) {if (t2[i] == '\n') t2[i] = ' ';}
    s += "# " + t2 + "\n";
  }
  s += "#\n";
  s += "# View this plot with gnuplot (5.0 or later):\n";
  s += "#     gnuplot -p " + filename + "\n";
  s += "# (the -p flag keeps the plot window open after drawing)\n";
  s += "#\n";
  s += "set termoption noenhanced\n";
  std::string title = title1;
  if (title2 != "") title += "\n" + title2;
  if (title != "") s += "set title " + GnuQuote(title) + "\n";
  if (xlabel != "") s += "set xlabel " + GnuQuote(xlabel) + "\n";
  if (ylabel != "") s += "set ylabel " + GnuQuote(ylabel) + "\n";
  char buf[200];
  snprintf(buf, sizeof(buf), "set xrange [%g:%g]\nset yrange [%g:%g]\n",
           xmin, xmax, ymin, ymax);
  s += buf;
  if (xtick > 0.0) {snprintf(buf, sizeof(buf), "set xtics %g\n", xtick); s += buf;}
  if (ytick > 0.0) {snprintf(buf, sizeof(buf), "set ytics %g\n", ytick); s += buf;}
  if (!ticklabel) s += "set format x \"\"\nset format y \"\"\n";
  s += "set grid\n";
  if (square) s += "set size ratio -1\n";
  if (diagonal) s += "set xzeroaxis lt -1\nset yzeroaxis lt -1\n";
  bool anylegend = false;
  for (size_t i=0;i<series.size();++i) {
    if (series[i].legend != "" && series[i].data != "") anylegend = true;
  }
  if (writelegend && anylegend) {
    s += "set key outside right top\n";
  } else {
    s += "unset key\n";
  }
  s += "\n";

  // Data, one datablock per series
  std::vector<std::string> clauses;
  for (size_t i=0;i<series.size();++i) {
    const Series& sr = series[i];
    if (sr.data == "") continue;  // gnuplot complains about empty data
    std::string name = "$S" + StringUtil::Strip(clipper::String(int(i)));
    s += name + " << EOD\n" + sr.data + "EOD\n";

    std::string c = name + " using 1:2 with ";
    if (sr.symbol != 0) {
      c += sr.join ? "linespoints" : "points";
      int pt = (sr.symbol < 0) ? 5 : 4;  // filled or open square
      float ps = Max(0.25f, 0.25f*float(std::abs(sr.symbol)));
      snprintf(buf, sizeof(buf), " pt %d ps %.2f", pt, ps);
      c += buf;
    } else {
      c += "lines";
    }
    c += " lc rgb " + GnuColour(sr.colour);
    c += (sr.legend != "") ? " title " + GnuQuote(sr.legend) : " notitle";
    clauses.push_back(c);
  }
  if (diagonal) {
    clauses.push_back("x with lines lc rgb \"#000000\" notitle");
  }
  s += "\n";
  if (!clauses.empty()) {
    s += "plot ";
    for (size_t i=0;i<clauses.size();++i) {
      if (i > 0) s += ", \\\n     ";
      s += clauses[i];
    }
    s += "\n";
  }

  FILE* file = OpenFile(filename, true);
  fputs(s.c_str(), file);
  fclose(file);
}
//--------------------------------------------------------------
