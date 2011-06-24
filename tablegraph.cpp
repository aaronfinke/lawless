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
//! construct as explicit X,Y ranges
GraphAxesType::GraphAxesType(const scala::Range& Xrange, const scala::Range& Yrange,
			     const bool& ZeroY)
{
  init(Xrange, Yrange, ZeroY);
}
//--------------------------------------------------------------
//! initialise as explicit X,Y ranges
void GraphAxesType::init(const scala::Range& Xrange, const scala::Range& Yrange,
			 const bool& ZeroY)
{
  xrange = Xrange;
  yrange = Yrange;
  graphtype = XY_SPECIFIED;
  FixYrange(ZeroY);
}
//--------------------------------------------------------------
void GraphAxesType::FixYrange(const bool& ZeroY)
// Fix Yrange to be "sensible"
// If ZeroY true, then y range should start at 0
{
  if (ZeroY) {
    yrange.first() = 0.0;
  }
}
//--------------------------------------------------------------
std::string GraphAxesType::FormatType() const
  //! return formatted for graph
{
  if (graphtype == AUTO_Y) {
    return "A";
  } else if (graphtype == NOUGHT_Y) {
    return "N";
  } else if (graphtype ==  XY_SPECIFIED) {
    std::string s =
      clipper::String(xrange.min())+"|"+clipper::String(xrange.max())+
      "x"+
      clipper::String(yrange.min())+"|"+clipper::String(yrange.max());
    return StringUtil::Strip(s);
  }
  return "A";
}
//--------------------------------------------------------------
const std::string TableGraph::LABELLEADER  = "$$\n";
const std::string TableGraph::LABELTRAILER = "  $$";
const std::string TableGraph::LABELFINAL   = " $$\n";
//--------------------------------------------------------------
TableGraph::TableGraph(const std::string& Title)
// construct & store title
  : title(Title), ngraphs(0)  {}
//--------------------------------------------------------------
// Store title
void TableGraph::init(const std::string& Title)
{
  title = Title;
  ngraphs = 0;
}
//--------------------------------------------------------------
std::string TableGraph::formatTitle() const
{ 
  return "\n$TABLE: "+title+":\n";
}
//--------------------------------------------------------------
//! Add a graph, return graph header
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
  ASSERT (Labels.size() == ZeroMark.size());
  ncolumns = Labels.size();
  labels = "";
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
	ifd = atoi(&pformat[++i]);
	i++;
	tp = +1;
      }
      int ifww = ifw+overhang;  // add previous overhang
      // Store label
      labels += AddInLabel(Labels[ifield], ifw, ifd, overhang);
      //^
      //^	std::cout << ifield << " " << ifw << " " << ifd << " " << overhang << "\n"; //^
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
  labels = LABELLEADER+labels+LABELTRAILER;
  if (lastmark) labels += LABELFINAL; // optional final "$$" mark
  else labels += "\n";
  return labels;
}
//--------------------------------------------------------------
//! Column labels string, omitting the leader and trailer (if present)
std::string TableGraph::RawLabels() const
{
  return labels.substr(LABELLEADER.size(), lablen);
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
  vsprintf(temp,prtf_format.c_str(),arglist);
  va_end(arglist);
  assert(temp[temp_size-1] == '\0');
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
	sprintf(buf, fields[i].fmt.c_str(), iv);
	sfld.assign(buf, fields[i].fieldwidth);
      }
    } else {
      // real
      fv = va_arg(arglist, double);
      if (fv == 0.0 && fields[i].dashpos >= 0) {
	sfld[fields[i].dashpos] = '-';
      } else {
	sprintf(buf, fields[i].fmt.c_str(), fv);
	sfld.assign(buf, fields[i].fieldwidth);
      }
    }
    line += sfld;
  }
  va_end(arglist);
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
    sprintf(buf, fields[kfield].fmt.c_str(), iv);
    sfld.assign(buf, fields[kfield].fieldwidth);
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
    sprintf(buf, fields[kfield].fmt.c_str(), v);
    sfld.assign(buf, fields[kfield].fieldwidth);
  }
  line += sfld;
  kfield++;
}
//--------------------------------------------------------------
std::string TableGraph::GetLine()
// return line assembled in AddToLine calls
{
  line += "\n";
  return line;
}
//--------------------------------------------------------------
std::string TableGraph::CloseTable() const
// Terminate the table
{
  return "$$\n";
}
//--------------------------------------------------------------

