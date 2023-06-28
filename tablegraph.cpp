//
// tablegraph.cpp
//
//  TableGraph class for writing table for loggraph
//
// Graph syntax (ccp4 6.1)
// ============
//
//  $TABLE :table name:
//  $GRAPHS :graph1 name:graphtype:column_list:
//          :graph2 name:graphtype:column_list:
//          :graph 3 ...: ... $$
//  column1_name column2_name ... $$ any_characters $$
//   numbers $$
//
//  graphtype is
//
//  A[UTO]
//     for fully automatic scaling (e.g. ... :A:1,2,4,5:)
//  N[OUGHT]
//    for automatic y coordinate scaling, where y lowest limit is 0
//    (e.g. ... :N:1,2,4,5:)
//  XMIN|XMAXxYMIN|YMAX
//    for user defined scaling where XMIN ... are axis limits
//    (e.g. ... :0|100x-1|1:1,2,4,5:)
//
//

#include "tablegraph.hh"

#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;
#include <cstdarg>
#include <stdlib.h>
#include <stdio.h>

#include "string_util.hh"

//--------------------------------------------------------------
// construct as explicit X,Y ranges
GraphAxesType::GraphAxesType(const scala::Range& Xrange,
                             const scala::Range& Yrange,
                             const bool& ZeroY)
{
  init(Xrange, Yrange, ZeroY);
}
//--------------------------------------------------------------
// initialise as explicit X,Y ranges
void GraphAxesType::init(const scala::Range& Xrange,
                         const scala::Range& Yrange,
                         const bool& ZeroY)
{
  graphtype = XY_SPECIFIED;
  xrange = Xrange;
  xinvresolsq = false;
  yrange = Yrange;
  FixYrange(ZeroY);
  zeroy_RH = false;
}
//--------------------------------------------------------------
// X-axis range and type flag (true for 1/d^2)
void GraphAxesType::SetXaxis(const scala::Range& Xrange,
                             const bool& isinvresolsq)
{
  xrange = Xrange;
  xinvresolsq = isinvresolsq;
}
//--------------------------------------------------------------
// Y-axis range and ZeroY true to start Y at 0
void GraphAxesType::SetYaxis(const scala::Range& Yrange, const bool& ZeroY)
{
  yrange = Yrange;
  FixYrange(ZeroY);
}
//--------------------------------------------------------------
// Right Y-axis range and ZeroY true to start Y at 0
void GraphAxesType::SetRightYaxis(const scala::Range& Yrange,
                                  const bool& ZeroY)
{
  yrange_RH = Yrange;
  FixYrangeRH(ZeroY);
}
//--------------------------------------------------------------
void GraphAxesType::FixYrange(const bool& ZeroY)
// Fix Yrange to be "sensible"
// If ZeroY true, then y range should start at 0
{
  if (ZeroY) {
    yrange.first() = 0.0;
  }
  zeroy = ZeroY;
}
//--------------------------------------------------------------
void GraphAxesType::FixYrangeRH(const bool& ZeroY)
// Fix Yrange to be "sensible"
// If ZeroY true, then y range should start at 0
{
  if (ZeroY) {
    yrange_RH.first() = 0.0;
  }
  zeroy_RH = ZeroY;
}
//--------------------------------------------------------------
std::string GraphAxesType::FormatType() const
// return formatted for loggraph
{
  GraphType graphtype1 = graphtype;
  if (xrange.Valid() && yrange.Valid()) {
    graphtype1 = XY_SPECIFIED;
  }
  if (graphtype1 == AUTO_Y) {
    if (zeroy) {
      return "N";
    }
    return "A";
  } else if (graphtype1 == NOUGHT_Y) {
    return "N";
  } else if (graphtype1 ==  XY_SPECIFIED) {
    std::string s =
      clipper::String(xrange.min())+"|"+clipper::String(xrange.max())+
      "x"+
      clipper::String(yrange.min())+"|"+clipper::String(yrange.max());
    return StringUtil::Strip(s);
  }
  return "A";
}
//--------------------------------------------------------------
//--------------------------------------------------------------
TablegraphLineStyle::TablegraphLineStyle(const std::string& colr,
                             const std::string& linestyle,
                             const int& linewidth)
{
  init(colr, linestyle, linewidth);
}
//--------------------------------------------------------------
void TablegraphLineStyle::init(const std::string& colr,
                         const std::string& linestyle,
                         const int& linewidth)
{
  SetColour(colr);
  SetLine(linestyle, linewidth);
}
//--------------------------------------------------------------
//! Set colour
void TablegraphLineStyle::SetColour(const std::string& colr)
{
  colour_ = colr;
}
//--------------------------------------------------------------
void TablegraphLineStyle::SetLine(const std::string& linestyle,
                            const int& width)
{
  slinestyle = Style(linestyle);  // standard value
  linesize = width; // default = -1, unspecified
  // The style of the line, allowed values:
  // '-','--','-.',':','.',
  // corresponding to: 'Solid','Dashed','Dash-dot','Dotted','Blank'.
  linestylevalue_ = "";
  if (slinestyle == "Solid") {
    linestylevalue_ = "-";
  } else if (slinestyle == "Dashed") {
    linestylevalue_ = "--";
  } else if (slinestyle == "Dash-dot") {
    linestylevalue_ = "-.";
  } else if (slinestyle == "Dotted") {
    linestylevalue_ = ":";
  } else if (slinestyle == "Blank") {
    linestylevalue_ = ".";
  }
}
//--------------------------------------------------------------
// convert string to standard linestyle string (static)
std::string TablegraphLineStyle::Style(const std::string& style)
{
  std::string upperstyle = StringUtil::ToUpper(style);
  if (upperstyle == "" || upperstyle == "DEFAULT") {
    return "Solid";
  }
  if (upperstyle == "SOLID") {
    return "Solid";
  }
  if (upperstyle == "DASHED") {
    return "Dashed";
  }
  if (upperstyle == "DASH-DOT" || upperstyle == "DASH_DOT") {
    return "Dash-dot";
  }
  if (upperstyle == "DOTTED") {
    return "Dotted";
  }
  if (upperstyle == "BLANK") {
    return "Blank";
  }
  return "Solid";
}
//--------------------------------------------------------------
//--------------------------------------------------------------
TableGraphPlotline::TableGraphPlotline() {
  init();
}
//--------------------------------------------------------------
TableGraphPlotline::TableGraphPlotline
(const int& Xcol, const int& Ycol, const std::string& colr,
 const std::string& symb, const int& symbsize, const bool& symbedge,
 const std::string& linestyle, const int& linewidth)
{
  init(Xcol, Ycol, colr, symb, symbsize, symbedge,
       linestyle, linewidth);
}
//--------------------------------------------------------------
void TableGraphPlotline::init() {
  init(1,2);
}
//--------------------------------------------------------------
// set xcol, ycol from 1, etc
void TableGraphPlotline::init(const int& Xcol, const int& Ycol,
                              const std::string& colr,
                              const std::string& symb, const int& symbsize,
                              const bool& symbedge,
                              const std::string& linestyle,
                              const int& linewidth)
{
  xcol = Xcol;
  ycol = Ycol;
  SetColour(colr);
  SetSymbol(symb, symbsize, symbedge);
  tablegraphlinestyle.init(colr, linestyle, linewidth);
  rhaxis = false;
  // +1 log file only, 0 [default] both, -1 XML only
  logorXML = 0;
}
//--------------------------------------------------------------
void TableGraphPlotline::SetSymbol(const std::string& symb,
                                   const int& size,
                                   const bool& edge) {
  // size default = -1 ie unspecified
  symbol = symb;
  symbolsize = size;
  symboledge = edge;
}
//--------------------------------------------------------------
void TableGraphPlotline::SetColour(const std::string& colr)
{
  colour = colr;
}
//--------------------------------------------------------------
// return formatted XML block
std::string TableGraphPlotline::XMLformat(const int& xcolbreak) const
// if xcolbreak >= 0, then use this column for x axis instead of xcol
{
  if (logorXML > 0) {return "";}
  std::string s;
  // <plotline xcol="ix" ycol="iy">
  int xc = xcol;
  if (xcolbreak >= 0) {xc = xcolbreak;}
  std::string sxcol = StringUtil::itos(xc,3);
  std::string sycol = StringUtil::itos(ycol,3);
  s += "\n<plotline xcol=\""+sxcol+"\" ycol=\""+sycol+"\"";
  if (rhaxis) {
    s += " rightaxis=\"true\"";
  }
  s += " >\n";
  if (symbol != "") {
    s += StringUtil::MakeXMLtag("symbol", symbol)+"\n";
  }
  if (symbolsize >= 0) {
    s += StringUtil::MakeXMLtag("symbolsize",StringUtil::itos(symbolsize,3))+"\n";
  }
  if (!symboledge) {
    s += StringUtil::MakeXMLtag("markeredgewidth",StringUtil::ftos(0.0,4,1))+"\n";
  }
  if (tablegraphlinestyle.linestylevalue() != "") {
    s += StringUtil::MakeXMLtag("linestyle",
                tablegraphlinestyle.linestylevalue())+"\n";
  }
  if (tablegraphlinestyle.linewidth() > 0) {
    s += StringUtil::MakeXMLtag("linesize",
                StringUtil::itos(tablegraphlinestyle.linewidth(),3))+"\n";
  }
  if (colour != "") {
    s += StringUtil::MakeXMLtag("colour", colour)+"\n";
  }
  if (label != "") {
    s += StringUtil::MakeXMLtag("label", label)+"\n";
  }
  s += "</plotline>\n"; // end block
  return s;
}
//--------------------------------------------------------------
// convert string to standard linestyle string
std::string TableGraphPlotline::Style(const std::string& style)
{
  return TablegraphLineStyle::Style(style);
}
//--------------------------------------------------------------
//--------------------------------------------------------------
TableGraphPlot::TableGraphPlot() {
  init("");
}
//--------------------------------------------------------------
TableGraphPlot::TableGraphPlot(const std::string& ptitle) {
  init(ptitle);
}
//--------------------------------------------------------------
void TableGraphPlot::init(const std::string& ptitle) {
  plottype = "xy";
  title = ptitle;
  description = "";
  xlabel = "";
  ylabel = "";
  xscale = "";
  yscale = "";
  xrange.clear();
  yrange.clear();
  yrange_RH.clear();
  isRHyaxis = true;
  xbreaks.clear();
  xcolbreak = -1;
  xinvresolsq = false;
  zeroy = false;
  zeroy_RH = false;
  ybreaks.clear();
  xintegral = false;
  yintegral = false;
  yintegral_RH = false;
  axistypes.init(GraphAxesType::AUTO_Y); // default Y axis type
  plotlines.clear();
  lines.clear();
}
//--------------------------------------------------------------
//--------------------------------------------------------------
TablegraphLine::TablegraphLine(const std::pair<double, double> XY1,
                               const std::pair<double, double> XY2,
                               const int& fw, const int& fd,
                               const std::string& colr,
                               const std::string& linestyle,
                               const int& linewidth)
{
  init(XY1, XY2, fw, fd, colr, linestyle, linewidth);
}
//--------------------------------------------------------------
void TablegraphLine::init(const std::pair<double, double> XY1,
                          const std::pair<double, double> XY2,
                          const int& fw, const int& fd,
                          const std::string& colr,
                          const std::string& linestyle,
                          const int& linewidth)
{
  xy1 = XY1;
  xy2 = XY2;

  fw_ = fw;
  fd_ = fd;

  tablegraphlinestyle.init(colr, linestyle, linewidth);
}
//--------------------------------------------------------------
//! Return XML format for Pimple
std::string TablegraphLine::XMLformat() const
{
  std::string s = "<line x1=\"" +
    StringUtil::ftos(xy1.first,fw_,fd_)+"\"" +
    " x2=\"" + StringUtil::ftos(xy2.first,fw_,fd_)+"\"" +
    " y1=\"" + StringUtil::ftos(xy1.second,fw_,fd_)+"\"" +
    " y2=\"" + StringUtil::ftos(xy2.second,fw_,fd_)+"\"";
  s += " linestyle=\""+tablegraphlinestyle.linestylevalue()+"\"";
  if (tablegraphlinestyle.linewidth() > 1) {
    s += " linesize=\""+
      StringUtil::itos(tablegraphlinestyle.linewidth(),3)+"\"";
  }
  if (tablegraphlinestyle.colour() != "") {
    s += " linecolour=\""+tablegraphlinestyle.colour()+"\"";
  }
  s += "/>\n";
  return s;
}
//--------------------------------------------------------------
//--------------------------------------------------------------
void TableGraphPlot::SetXaxis(const std::string& label,
                              const bool& isinvresolsq,
                              const scala::Range& range,
                              const bool& integral)
//  Define X-axis:
//  label    for axis, "" to get from data table
//  isinvresolsq true if x axis is 1/d^2
//  range    axis range, null for auto determination
//  integral true if axis values are integral
{
  if (label != "") {
    xlabel = label;
  }
  xinvresolsq = isinvresolsq;
  if (range.Valid()) {
    xrange = range;
  }
  xintegral = integral;

  axistypes.SetXaxis(xrange, xinvresolsq);
}
//--------------------------------------------------------------
void TableGraphPlot::SetXbreak(const int& xcolbr,
                               const std::vector<scala::Range>& xbreak,
                               const scala::Range& rangexbreak)
{
  // all x-breaks must refer to the same column, check
  if (xcolbreak < 0) {
    xcolbreak = xcolbr;
  } else if (xcolbreak != xcolbr) {
    Message::message(Message_fatal("TableGraphPlot::SetXbreak: different xcolbreak"));
  }
  if (xbreak.size() > 0) {
    for (size_t i=0;i<xbreak.size();++i) {
      xbreaks.push_back(xbreak[i]);
    }}
  if (rangexbreak.Valid()) {
    xbreakrange = rangexbreak;
  }
}
//--------------------------------------------------------------
void TableGraphPlot::SetYaxis(const std::string& label,
                              const bool& ZeroY,
                              const scala::Range& range,
                              const bool& integral)
// Define Y-axis:
//  label    for axis, "" to get from data table
//  ZeroY    true to run y from zero
//  range    axis range, null for auto determination
//  integral true if axis values are integral
{
  if (label != "") {
    ylabel = label;
  }
  if (range.Valid()) {
    yrange = range;
    if (ZeroY) {yrange.first() = 0.0;}
    if (yrange.last() == yrange.first()) {
      yrange.last() = yrange.first() + 1.0;
    }
  }
  yintegral = integral;
  zeroy = ZeroY;
  axistypes.SetYaxis(yrange, zeroy);
}
//--------------------------------------------------------------
void TableGraphPlot::SetRightYaxis(const std::string& label,
                                   const bool& ZeroY,
                                   const scala::Range& range,
                                   const bool& integral)
// Define right-hand Y-axis:
//  label    for axis, "" to get from data table
//  ZeroY    true to run y from zero
//  range    axis range, null for auto determination
//  integral true if axis values are integral
{
  if (label != "") {
    ylabel_RH = label;
  }
  if (range.Valid()) {
    yrange_RH = range;
    if (ZeroY) {yrange_RH.first() = 0.0;}
  }
  yintegral_RH = integral;
  zeroy_RH = ZeroY;
  axistypes.SetRightYaxis(yrange_RH, zeroy_RH);
  isRHyaxis = true;
}
//--------------------------------------------------------------
// Add a line to the plot
void TableGraphPlot::AddLine(const TableGraphPlotline& pltline)
{
  // Check for consistent RH axis specification
  if (pltline.IsRHaxis() && !isRHyaxis) {
    Message::message(Message_fatal
                     ("TableGraphPlot::AddLine: must specify RH yaxis to add RH axis line"));
  }
  plotlines.push_back(pltline);
}
//--------------------------------------------------------------
// Add a simple line to the plot
void TableGraphPlot::AddPlainLine(const TablegraphLine& plainline)
{
  lines.push_back(plainline);
}
//--------------------------------------------------------------
std::string TableGraphPlot::formatXbreaks() const
{
  std::string s = "<xbreaks>\n";
  for (size_t i=0;i<xbreaks.size();++i) {
    std::string symin;
    std::string symax;
    if (xintegral) {
      symin = StringUtil::itos(Nint(xbreaks[i].min()));
      symax = StringUtil::itos(Nint(xbreaks[i].max()));
    } else {
      symin = StringUtil::ftos(xbreaks[i].min());
      symax = StringUtil::ftos(xbreaks[i].max());
    }
    s += "<break min=\""+symin+"\" max=\""+symax+"\"/>\n";
  }
  s += "</xbreaks>\n";
  return s;
}
//--------------------------------------------------------------
// XML format for Pimple
std::string TableGraphPlot::XMLformat() const
{
  std::string s = "<plot>\n";
  if (title != "") {
    s += StringUtil::MakeXMLtag("title", title)+"\n";
  }

  if (description != "") {
    s += StringUtil::MakeXMLtag("description", description)+"\n";
  }

  if (xlabel != "") {
    s += StringUtil::MakeXMLtag("xlabel", xlabel)+"\n";
  }
  if (ylabel != "") {
    s += StringUtil::MakeXMLtag("ylabel", ylabel)+"\n";
  }
  if (ylabel_RH != "") {
    s += StringUtil::MakeXMLtag("rylabel", ylabel_RH)+"\n";
  }
  if (xinvresolsq) {
    s += StringUtil::MakeXMLtag("xscale", "oneoversqrt")+"\n";
  }
  if (xrange.Valid()) {
    // <xrange min="xmin" max="xmax"\>
    std::string sxmin = StringUtil::ftos(xrange.min());
    std::string sxmax = StringUtil::ftos(xrange.max());
    s += "<xrange min=\""+sxmin+"\" max=\""+sxmax+"\"/>\n";
  } else if (xbreaks.size() > 0) {
    if (xbreakrange.Valid()) {
      // no range, but some xbreaks, use alternative xrange
      std::string sxmin = StringUtil::ftos(xbreakrange.min());
      std::string sxmax = StringUtil::ftos(xbreakrange.max());
      s += "<xrange min=\""+sxmin+"\" max=\""+sxmax+"\"/>\n";
    }
  }

  if (yrange.Valid()) {
    // <yrange min="ymin" max="ymax"\>
    std::string symin = StringUtil::ftos(yrange.min());
    std::string symax = StringUtil::ftos(yrange.max());
    s += "<yrange min=\""+symin+"\" max=\""+symax+"\"/>\n";
  } else if (zeroy) {
    // <yrange min="0" max="None"\>
    s += "<yrange min=\"0\" max=\"None\"/>\n";
  }
  if (yrange_RH.Valid()) {
    // <yrange min="ymin" max="ymax"\>
    std::string symin = StringUtil::ftos(yrange_RH.min());
    std::string symax = StringUtil::ftos(yrange_RH.max());
    s += "<yrange min=\""+symin+"\" max=\""+symax+"\" rightaxis=\"true\"/>\n";
  } else if (zeroy_RH) {
    // <yrange min="0" max="None"\>
    s += "<yrange min=\"0\" max=\"None\" rightaxis=\"true\"/>\n";
  }

  if (xbreaks.size() > 0) {
    s += formatXbreaks();
  }
  if (xintegral) {
    s += StringUtil::MakeXMLtag("xintegral", "true");
  }
  if (yintegral) {
    s += StringUtil::MakeXMLtag("yintegral", "true");
  }
  // Reverse order of plotlines so that the first is last and
  // therefore on top
  if (plotlines.size() > 0) {
    for (int i=int(plotlines.size())-1;i>=0;--i) {
      s += plotlines[i].XMLformat(xcolbreak);
    }
  }

  if (lines.size() > 0) {
    for (int i=int(lines.size())-1;i>=0;--i) {
      s += lines[i].XMLformat();
    }
  }

  s += "</plot>\n";
  return s;
}
//--------------------------------------------------------------
std::string TableGraphPlot::format(const bool& first) const
// format a graph for loggraph
// first  true for first graph
{
  std::string text;
  // Check type
  //  if (axistypes != "A" && axistypes != "N" &&
  //      (axistypes.find("|") == std::string::npos)) {
  //    Message::message(Message_fatal("TableGraph::Graph: invalid graph type:"+
  //                               Graphtype));
  //  }
  if (first) {text += "$GRAPHS";}
  text += ":"+title+":"+axistypes.FormatType()+":";
  if (plotlines.size() <= 0) {
    Message::message(Message_fatal("TableGraph::Graph: no lines"));
  }
  std::vector<int> columnNumbers;
  int xcol = plotlines[0].Xcol();
  for (size_t i=0;i<plotlines.size();++i) {
    if (plotlines[i].LogorXML() >= 0) {
      if (plotlines[i].Xcol() != xcol) {
        Message::message(Message_fatal
                         ("TableGraph::Graph: all lines must have same x column"));
      }
      if (i==0) {columnNumbers.push_back(xcol);}
      columnNumbers.push_back(plotlines[i].Ycol());
    }
  }
  std::string numbers;
  for (size_t i=0;i<columnNumbers.size();++i) {
    numbers += clipper::String(columnNumbers[i]);
    if (int(i)<int(columnNumbers.size())-1) {numbers+=",";}  // comma separated
  }
  text += StringUtil::Strip(numbers)+":\n";
  return text;
}
//--------------------------------------------------------------
TableGraph::TableGraph(const std::string& Title)
// construct & store title
  : title(Title), ngraphs(0), id("")  {}
//--------------------------------------------------------------
// Store title
void TableGraph::init(const std::string& Title)
{
  title = Title;
  extratitle = "";
  ngraphs = 0;
  id = "";
  sdatatable.clear();
}
//--------------------------------------------------------------
// Add a graph, return graph header
std::string TableGraph::Graph(const std::string& GraphTitle,
                              const GraphAxesType& Graphaxestype,
                              const std::vector<int>& columnNumbers)
{
  return Graph(GraphTitle, Graphaxestype.FormatType(), columnNumbers);
}
//--------------------------------------------------------------
std::string TableGraph::Graph(const std::string& GraphTitle,
                              const std::string& Graphtype,
                              const std::vector<int>& columnNumbers)
// Add a graph
{
  std::string text;
  // Check type
  if (Graphtype != "A" && Graphtype != "N" &&
      (Graphtype.find("|") == std::string::npos)) {
    Message::message(Message_fatal("TableGraph::Graph: invalid graph type:"+
                                   Graphtype));
  }
  if (ngraphs == 0) {text += "$GRAPHS";}
  ngraphs++;
  text += ":"+GraphTitle+":"+Graphtype+":";
  std::string numbers;
  for (size_t i=0;i<columnNumbers.size();++i) {
    numbers += clipper::String(columnNumbers[i]);
    if (int(i)<int(columnNumbers.size())-1) {numbers+=",";}  // comma separated
  }
  text += StringUtil::Strip(numbers)+":\n";
  return text;
}
//--------------------------------------------------------------
void TableGraph::AddGraph(const TableGraphPlot& tgplot)
// Add and store a graph
{
  graphs.push_back(tgplot);
  ngraphs = graphs.size();
}
//--------------------------------------------------------------
std::string TableGraph::AddInLabel(const std::string& label,
                                   const int& ifw_in, const int& ifd,
                                   int& overhang) const
// returns overhang = number of characters after field (maximum =1)
// Private function
{
  int icp;          // centre position in field
  overhang = 0;
  int len = label.size();
  int ifw = ifw_in;
  if (ifd == 0) {
    // Integer, try to right-justify, but allow overhang
    if (ifw-len < 2) {
      overhang = 1;
      ifw += overhang;
    }
    //^      std::cout << "\nAILd " << label << " " << ifw << "\n";  //^
    return StringUtil::RightString(label, ifw);  // integer, right justify
  } else {
    icp = ifw - ifd/2 -2;  // real, centre in field
    icp = clipper::Util::min(icp, ifw-(len+1)/2);
    overhang = clipper::Util::max(0,icp - ifw + (len+1)/2);
  }
  //^    std::cout << "\nAILf " << label << " " << ifw << " " << ifd << " " << icp << " "<< overhang<<"\n";//^
  //^    std::cout<< CentreString(label, ifw, icp) << "\n";; //^
  return StringUtil::CentreString(label, ifw, icp);
}
//--------------------------------------------------------------
std::string TableGraph::ColumnFields(const std::vector<std::string>& Labels,
                                     const std::vector<bool>& ZeroMark,
                                     const std::string& pformat,
                                     const bool& lastmark)
// Define format, & write out labels
// ZeroMark = true to replace zero value with "-"
// if lastmark true [default] add final "$$" after headers
//
{
  StoreColumnFields(Labels, ZeroMark, pformat);
  return GetLabels(lastmark);
}
//--------------------------------------------------------------
std::string TableGraph::GetLabels(const bool& lastmark) const
// format label string for loggraph
{
  const std::string LABELLEADER  = "$$\n";
  const std::string LABELTRAILER = "  $$";
  const std::string LABELFINAL   = " $$\n";

  std::string lglabels = LABELLEADER+labels+LABELTRAILER;
  if (lastmark) lglabels += LABELFINAL; // optional final "$$" mark
  else lglabels += "\n";
  return lglabels;
}
//--------------------------------------------------------------
void TableGraph::StoreColumnFields(const std::vector<std::string>& Labels,
                                   const std::vector<bool>& ZeroMark,
                                   const std::string& pformat)
// Define format for data table
// ZeroMark = true to replace zero value with "-"
//
{
  ASSERT (Labels.size() == ZeroMark.size());
  ncolumns = Labels.size();
  labels = "";
  labelarray = Labels;
  prtf_format = pformat;
  //  char bs = '\\';
  //^    std::cout << pformat << "\n"; //^

  // Parse format
  int i=-1;
  int ifield = -1;  // field count
  int i0=0;  // start of field
  int ifw = 0;
  int ifd = 0;      // number after decimal point
  int overhang=0;
  while (++i < int(pformat.size())) {
    if (pformat[i] == '%') {
      ifield++; // new field
      // number format field
      i++;
      std::string fw;
      while (std::isdigit(pformat[i])) {
        fw += pformat[i++];
      }
      ifw += atoi(fw.c_str());
      ifd = 0;
      int tp = 0;
      if (pformat[i] == '.') {    // should be either '.' or 'd'
        // Assume single digit decimal count
        ifd = atoi(std::string(1,pformat[++i]).c_str());
        i++;
        tp = +1;
      }
      int ifww = ifw+overhang;  // add previous overhang
      // Store label
      //^
      //      std::cout << ifield << " " << ifw << " " << ifd << " " << overhang << "\n"; //^
      labels += AddInLabel(Labels.at(ifield), ifw, ifd, overhang);
      // Store field information
      //   field width
      //   position for "-" character: right-justified for integer
      std::string sfmt = pformat.substr(i0, i-i0+1);
      int dp = ifww - ifd - 1;
      if (ifd == 0) {dp = ifww-1;}
      if (!ZeroMark[ifield]) dp = -1;  // don't replace zeroes
      fields.push_back(fieldInfo(ifww, dp, tp, sfmt));
      i0 = i+1;  // start of next format field
      // next char [i] is 'd' or 'f', will be ignored
      ifw = -overhang;
    } else {
      ifw++; // count leading characters before '%'
    }
  }
  if (ncolumns != ifield+1) {
    Message::message(Message_fatal
                     ("TableGraph: wrong number of fields in format"));
  }
  lablen = labels.size();
}
//--------------------------------------------------------------
// Column labels string, omitting the leader and trailer (if present)
std::string TableGraph::RawLabels() const
{
  return labels;
}
//--------------------------------------------------------------
std::string TableGraph::NumberLine(const int nc, ...) const
// Write nc numbers, using predefined format
// This will probably fail if the number of arguments doesn't match
// the format
{
  ASSERT (ncolumns == nc);
  // cf Output::logTabPrintf
  static const std::size_t temp_size = 8192;
  char temp[temp_size];
  temp[temp_size-1] = '\0';
  va_list arglist;
  va_start(arglist, nc);
  vsnprintf(temp,8192,prtf_format.c_str(),arglist);
  va_end(arglist);
  assert(temp[temp_size-1] == '\0');
  sdatatable.push_back(std::string(temp));
  return std::string(temp);
}
//--------------------------------------------------------------
std::string TableGraph::Line(const int nc, ...) const
// Write nc numbers, using predefined format, replacing zeroes by "-"
// This will probably fail if the number of arguments doesn't match
// the format
{
  ASSERT (ncolumns == nc);
  va_list arglist;
  va_start(arglist, nc);
  static const std::size_t buf_size = 8192;
  char buf[buf_size];
  buf[buf_size-1] = '\0';
  int iv;
  float fv;
  std::string sfld;
  line = "";

  for (int i=0;i<nc;++i) {
    sfld.assign(fields[i].fieldwidth,' ');
    // What type?
    if (fields[i].type == 0) {
      // Integer
      iv = va_arg(arglist, int);
      if (iv == 0 && fields[i].dashpos >= 0) {
        sfld[fields[i].dashpos] = '-';
      } else {
        snprintf(buf, buf_size, fields[i].fmt.c_str(), iv);
        sfld.assign(buf, fields[i].fieldwidth);
      }
    } else {
      // real
      fv = va_arg(arglist, double);
      if (fv == 0.0 && fields[i].dashpos >= 0) {
        sfld[fields[i].dashpos] = '-';
      } else {
        snprintf(buf, buf_size, fields[i].fmt.c_str(), fv);
        sfld.assign(buf, fields[i].fieldwidth);
      }
    }
    if (sfld[0] != ' ') {
      sfld = ' ' + sfld;
    }
    line += sfld;
  }
  va_end(arglist);
  sdatatable.push_back(line); // store without lf
  line += "\n";
  return line;
}
//--------------------------------------------------------------
std::string TableGraph::Line(const std::vector<double>& val, const int nc, ...) const
// Write nc numbers, then vector val
// using predefined format, optionally replacing zeroes by "-"
// This will probably fail if the number of arguments doesn't match
// the format
{
  va_list arglist;
  va_start(arglist, nc);
  static const std::size_t buf_size = 8192;
  char buf[buf_size];
  buf[buf_size-1] = '\0';
  int iv;
  float fv;
  std::string sfld;
  line = "";
  int nval = val.size();  // length of vector
  int ntot = nc + nval;
  ASSERT (ncolumns == ntot);
  int k = 0;

  for (int i=0;i<ntot;++i) {
    sfld.assign(fields[i].fieldwidth,' ');
    // What type?
    if (fields[i].type == 0) {
      // Integer
      if (i >= nc) { // take from vector
        iv = Nint(val[k++]);
      } else {
        iv = va_arg(arglist, int);
      }
      if (iv == 0 && fields[i].dashpos >= 0) {
        sfld[fields[i].dashpos] = '-';
      } else {
        snprintf(buf, buf_size, fields[i].fmt.c_str(), iv);
        sfld.assign(buf, fields[i].fieldwidth);
      }
    } else {
      // real
      if (i >= nc) { // take from vector
        fv = val[k++];
      } else {
        fv = va_arg(arglist, double);
      }
      if (fv == 0.0 && fields[i].dashpos >= 0) {
        sfld[fields[i].dashpos] = '-';
      } else {
        snprintf(buf, buf_size, fields[i].fmt.c_str(), fv);
        sfld.assign(buf, fields[i].fieldwidth);
      }
    }
    if (sfld[0] != ' ') {
      sfld = ' ' + sfld;
    }
    line += sfld;
  }
  va_end(arglist);
  sdatatable.push_back(line);
  line += "\n";
  return line;
}
//--------------------------------------------------------------
void TableGraph::StartLine()
// clear line
{line = ""; kfield=0;}
//--------------------------------------------------------------
void TableGraph::AddToLine(const int& iv)
// Add integer to next field in line
{
  if (fields[kfield].type != 0) {
    Message::message(Message_fatal
                     ("TableGraph: adding integer to non-integer field "));
  }
  std::string sfld(fields[kfield].fieldwidth,' ');
  if (iv == 0 && fields[kfield].dashpos >= 0) {
    sfld[fields[kfield].dashpos] = '-';
  } else {
    char buf[256];
    snprintf(buf, 256, fields[kfield].fmt.c_str(), iv);
    sfld.assign(buf, fields[kfield].fieldwidth);
  }
  if (sfld[0] != ' ') {
    sfld = ' ' + sfld;
  }
  line += sfld;
  kfield++;
}
//--------------------------------------------------------------
void TableGraph::AddToLine(const float& v)
// Add float to next field in line
{
  if (fields[kfield].type != +1) {
    Message::message(Message_fatal
                     ("TableGraph: adding float to non-float field "));
  }
  std::string sfld(fields[kfield].fieldwidth,' ');
  if (v == 0.0 && fields[kfield].dashpos >= 0) {
    sfld[fields[kfield].dashpos] = '-';
  } else {
    char buf[256];
    snprintf(buf, 256, fields[kfield].fmt.c_str(), v);
    sfld.assign(buf, fields[kfield].fieldwidth);
  }
  if (sfld[0] != ' ') {
    sfld = ' ' + sfld;
  }
  line += sfld;
  kfield++;
}

//--------------------------------------------------------------
void TableGraph::patchLine(const std::string& patch, const int& firstchar)
// replace part of "line" with patch, starting at character firstchar
{
  int plen = patch.size(); // length of patch string
  int llen = line.size();  // length of current line
  int len = plen;
  if (firstchar + plen > llen) {
    len = llen - firstchar;
  }
  if (len > 0) {
    line.replace(firstchar, len, patch, 0, len);
  }
}
//--------------------------------------------------------------
std::string TableGraph::GetLine()
// return line assembled in AddToLine calls
{
  sdatatable.push_back(line);
  line += "\n";
  return line;
}
//--------------------------------------------------------------
std::string TableGraph::XMLformat() const
{
  // Title
  std::string s = "\n<CCP4Table ";
  // Always label as a Graph
  s += "groupID=\"Graph\" ";
  if (id != "") {
    s += "id=\""+id+"\" ";
  }
  s += "title=\""+StringUtil::XMLstring(title)+"\">\n";
  // graph headers
  for (int igr=0;igr<ngraphs;++igr) {
    s += graphs[igr].XMLformat();
  }
  // Column labels
  s += "<headers separator=\" \">\n";
  s += StringUtil::XMLstring(labels)+"\n";
  s += "</headers>\n";
  // Data lines
  s += "<data>\n";
  for (size_t i=0;i<sdatatable.size();++i) {
    s += sdatatable[i]+"\n";
  }
  s += "</data>\n";
  // Terminate the table for XML
  return s+"</CCP4Table>\n";
}
//--------------------------------------------------------------
std::string TableGraph::formatTitle() const
{
  std::string s = "\n$TABLE: "+title+":\n";
  return s;
}
//--------------------------------------------------------------
std::string TableGraph::format() const
{
  std::string s = "\n$TABLE: "+title+":\n";
  // graph headers
  for (int igr=0;igr<ngraphs;++igr) {
    bool first = (igr == 0);
    s += graphs[igr].format(first);
  }
  if (extratitle != "") {
    s += extratitle+"\n";
  }
  s += " $$\n";
  // Column labels
  s += labels+"   $$ $$\n";
  // Data lines
  for (size_t i=0;i<sdatatable.size();++i) {
    s += sdatatable[i]+"\n";
  }
  return s+CloseTable();
}
//--------------------------------------------------------------
std::string TableGraph::CloseTable() const
{
  // Terminate the table
  return "$$\n";
}
//--------------------------------------------------------------
