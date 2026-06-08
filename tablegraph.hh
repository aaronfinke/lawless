//
//   tablegraph.hh
//
//!  TableGraph class for writing table for loggraph or qloggraph (aka Pimple)
//
// LogGraph syntax (ccp4 6.1)
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

#ifndef TABLEGRAPH_HEADER
#define TABLEGRAPH_HEADER

#include <assert.h>
#define ASSERT assert

#include <string>
#include <vector>
#include "range.hh"

class GraphAxesType {
  //! graph axis type for both x and y axes
  //
  //!  [Q]Loggraph options are:
  //!  A[UTO]
  //!     for fully automatic scaling of y axis (e.g. ... :A:1,2,4,5:)
  //!  N[OUGHT]
  //!    for automatic y coordinate scaling, where y lowest limit is 0
  //!    (e.g. ... :N:1,2,4,5:)
  //!  XMIN|XMAXxYMIN|YMAX
  //!    for user defined scaling where XMIN ... are axis limits
  //!    (e.g. ... :0|100x-1|1:1,2,4,5:)
  //!  ONEOVERSQRT for resolution x-axis
  //!
  //!  these may be extended in future
public:
  enum GraphType {AUTO_Y, NOUGHT_Y, XY_SPECIFIED};
  //! construct as default
  GraphAxesType() : graphtype(AUTO_Y), xinvresolsq(false),
		    zeroy(false) {} // default
  //! construct as specified
  GraphAxesType(const GraphType& gt) : graphtype(gt){}
  //! construct as explicit X,Y ranges, ZeroY true to start Y at 0
  GraphAxesType(const scala::Range& Xrange, const scala::Range& Yrange,
		const bool& ZeroY);
  //! initialise as specified
  void init(const GraphType& gt)  {graphtype = gt;}
  //! initialise as explicit X,Y ranges, ZeroY true to start Y at 0
  void init(const scala::Range& Xrange, const scala::Range& Yrange,
	    const bool& ZeroY);
  //! X-axis range and type flag (true for 1/d^2)
  void SetXaxis(const scala::Range& Xrange, const bool& isinvresolsq);
  //! Y-axis range and ZeroY true to start Y at 0
  void SetYaxis(const scala::Range& Yrange, const bool& ZeroY);
  //! Right Y-axis range and ZeroY true to start Y at 0
  void SetRightYaxis(const scala::Range& Yrange, const bool& ZeroY);

  //! return type
  GraphType Graphtype() const {return graphtype;}

  //! return formatted for loggraph
  std::string FormatType() const;

private:
  GraphType graphtype;
  scala::Range xrange;
  bool xinvresolsq;  // true if x axis is 1/d^2
  scala::Range yrange;
  bool zeroy;  // true if Y axis should start at zero
  scala::Range yrange_RH;  // RH axis
  bool zeroy_RH;  // true if Y axis should start at zero

  // Fix Yrange to be "sensible"
  // If ZeroY true, then y range should start at 0
  void FixYrange(const bool& ZeroY);
  void FixYrangeRH(const bool& ZeroY);
};


class TablegraphLineStyle {
  //! Linestyles are specified as a string, one of 
  //! 'Solid','Dashed','Dash-dot','Dotted','Blank' (case-insensitive)
  //! corresponding to: '-','--','-.',':','.' as used by Pimple
  //!  also line width (integer) and colour 
  //! Colour: allowed values red, green, blue, yellow, magenta, cyan, black,
  //! or r, g, b, y, m, c, k,
  //! or any other colour specification understood by matplotlib,
  //! e.g. orange, #ff7700.
public:
  TablegraphLineStyle() : linesize(-1) {}
  
  TablegraphLineStyle(const std::string& colr, const std::string& linestyle="",
		const int& linewidth=-1);
  
  void init(const std::string& colr="", const std::string& linestyle="",
	    const int& linewidth=-1);

  //! Set line type
  void SetLine(const std::string& linestyle, const int& width=-1);
  //! Set colour
  void SetColour(const std::string& colr);

  std::string colour() const {return colour_;}
  std::string linestylevalue() const {return linestylevalue_;}
  int linewidth() const {return linesize;}
  
  static std::string Style(const std::string& style);

private:
  std::string colour_;
  // 'Solid','Dashed','Dash-dot','Dotted','Blank' as input
  std::string slinestyle;
  // The style of the line, allowed values for Pimple: '-','--','-.',':','.',
  std::string linestylevalue_;
  int linesize; // line width
};

// Definitions in TableGraphPlotline & TableGraphPlot follow Pimple,
// though not all options are currently supported
// Many defaults can be left blank

//! a plotline corresponds to a set of data points drawn with a style
class TableGraphPlotline {
  //! Only Xcol, Ycol are relevant for loggraph
public:
  TableGraphPlotline();
  //! Construct with properties

  //! \param Xcol, Ycol column numbers for x, y (from 1)
  //! \param colour string (see below)
  //! \param symbol a character to define the symbol drawn (see below)
  //! \param symbolsize = -1 for default
  //! \param linestyle string, eg "SOLID"
  //! \param linewidth = -1 for default
  //!
  //! Colour: allowed values red, green, blue, yellow, magenta, cyan, black,
  //! or r, g, b, y, m, c, k,
  //! or any other colour specification understood by matplotlib,
  //! e.g. orange, #ff7700.
  //!
  //! Symbol strings are 
  //! 'o','s','d','^','>','<','p','h','*','$\\lambda$',
  //! '$\\bowtie$', '$\\circlearrowleft$', '$\\clubsuit$', '$\\checkmark$'.
  //! These values correspond to:
  //! 'Circle','Square','Diamond','Up arrow','Right arrow','Left arrow',
  //! 'Pentagon','Hexagon','Star', 'Lambda',
  //! 'Bow tie','Circle arrow left','Clubs (as in playing cards)','Tick'. 
  //!
  //! Linestyles are
  //! 'Solid','Dashed','Dash-dot','Dotted','Blank' (case-insensitive)
  //! corresponding to: '-','--','-.',':','.',

  TableGraphPlotline(const int& Xcol, const int& Ycol,
		     const std::string& colour="",
		     const std::string& symbol="",
		     const int& symbolsize=-1,
		     const bool& symboledge=true,
		     const std::string& linestyle="",
		     const int& linewidth=-1);

  void init();
  //! Initialise with all properties (or xcol,ycol), cf constructor
  void init(const int& Xcol, const int& Ycol,
	    const std::string& colr="",
	    const std::string& symbol="",
	    const int& symbolsize=-1,
	    const bool& symboledge=true,
	    const std::string& linestyle="",
	    const int& linewidth=-1);
  //! Set line symbol type
  void SetSymbol(const std::string& symb, const int& size=-1,
		 const bool& edge=true);
  //! Set line type
  void SetLine(const std::string& linestyle, const int& width=-1);
  //! Set colour
  void SetColour(const std::string& col);
  //! Set right-hand axis
  void SetRHaxis() {rhaxis = true;}
  //! True if RH axis specified
  bool IsRHaxis() const {return rhaxis;}

  //! Set logfile or XML
  void SetLogorXML(const int& logorxml) {logorXML = logorxml;}
  // +1 log file only, 0 [default] both, -1 XML only

  int LogorXML() const {return logorXML;}

  int Xcol() const {return xcol;} //!< return x-column
  int Ycol() const {return ycol;} //!< return y-column

  //! return formatted XML block
  std::string XMLformat(const int& xcolbreak) const;
  //!< if xcolbreak >= 0, then use this column for x axis instead of xcol

private:
  int xcol, ycol;  // x & y columns numbers, from 1

  // symbol values:
  // 'o','s','d','^','>','<','p','h','*','$\lambda$',
  // '$\bowtie$', '$\circlearrowleft$', '$\clubsuit$', '$\checkmark$'.
  // These values correspond to:
  // 'Circle','Square','Diamond','Up arrow','Right arrow','Left arrow',
  // 'Pentagon','Hexagon','Star', 'Lambda',
  // 'Bow tie','Circle arrow left','Clubs (as in playing cards)','Tick'. 
  std::string symbol;
  int symbolsize;
  bool symboledge;  // false for no black edge

  TablegraphLineStyle tablegraphlinestyle;

  // Colour: allowed values red, green, blue, yellow, magenta, cyan, black,
  // or r, g, b, y, m, c, k,
  // or any other colour specification understood by matplotlib,
  // e.g. orange, #ff7700.
  std::string colour; 
  std::string label; // for legend
  bool showinlegend; // show line in legend, not used yet

  bool rhaxis;  // true if line belongs to RH axis

  // +1 log file only, 0 [default] both, -1 XML only
  int logorXML;

  // convert string to LineStyle
  static std::string Style(const std::string& style);

};

class TablegraphLine {
  // Just a line in the plot (XML only)
public:
  TablegraphLine(const std::pair<double, double> XY1,
		 const std::pair<double, double> XY2,
		 const int& fw, const int& fd,
		 const std::string& colr="",
		 const std::string& linestyle="",
		 const int& linewidth=-1);

  void init(const std::pair<double, double> XY1,
	    const std::pair<double, double> XY2,
	    const int& fw, const int& fd,
	    const std::string& colr="",
	    const std::string& linestyle="",
	    const int& linewidth=-1);

  //! Return XML format for Pimple
  std::string XMLformat() const;
  
private:
  std::pair<double, double> xy1; // line start
  std::pair<double, double> xy2; // line end

  int fw_, fd_;
  
  TablegraphLineStyle tablegraphlinestyle;
};

//! A plot contains one or more plotlines, etc, and styles
class TableGraphPlot {
public:
  TableGraphPlot();
  //! Construct with graph title
  TableGraphPlot(const std::string& ptitle);

  //! Initialise with graph title
  void init(const std::string& ptitle);

  //! set description text
  void SetDescription(const std::string& Description) {description = Description;}

  //! Define X-axis properties

  //! \param label    for axis, "" to get from data table
  //! \param isinvresolsq true if x axis is 1/d^2
  //! \param range    axis range, null for auto determination
  //! \param integral true if axis values are integral
  void SetXaxis(const std::string& label, const bool& isinvresolsq,
		const scala::Range& range=scala::Range(),
		const bool& integral=false);

  //! Define breaks in the X-axis

  //! \param xcolbr  column number from which to take the broken x value
  //! \param xbreaks list of ranges of breaks (may be empty)
  void SetXbreak(const int& xcolbr, const std::vector<scala::Range>& xbreak,
		 const scala::Range& rangexbreak=scala::Range());

  //! Define Y-axis properties

  //! \param  label    for axis, "" to get from data table
  //! \param  ZeroY    true to run y from zero
  //! \param  range    axis range, null for auto determination
  //! \param  integral true if axis values are integral
  void SetYaxis(const std::string& label, const bool& ZeroY,
		const scala::Range& range=scala::Range(),
		const bool& integral=false);

  //! Define right-hand Y-axis properties

  //! \param  label    for axis, "" to get from data table
  //! \param  ZeroY    true to run y from zero
  //! \param  range    axis range, null for auto determination
  //! \param  integral true if axis values are integral
  //! Must call before AddLine for a RH axis line
  //
  void SetRightYaxis(const std::string& label, const bool& ZeroY,
		     const scala::Range& range=scala::Range(),
		     const bool& integral=false);

  //! Add a line to the plot
  void AddLine(const TableGraphPlotline& pltline);

  //! Add a simple line to the plot
  void AddPlainLine(const TablegraphLine& plainline);

  //! Return XML format for Pimple
  std::string XMLformat() const;
  //! Return format for loggraph
  std::string format(const bool& first) const;

private:
  std::string plottype; // "xy"
  std::string title;
  std::string description;  // optional description
  std::string xlabel, ylabel, ylabel_RH;  // axis labels, if specified
  std::string xscale; // blank or "oneoversqrt" (mostly not needed);
  std::string yscale; // not used
  scala::Range xrange, yrange, yrange_RH; // axis ranges
  bool isRHyaxis;     // true if a RH Y-axis has been specified
  std::vector<scala::Range> xbreaks; // breaks in x axis  
  int xcolbreak;      // x column for breaks
  scala::Range xbreakrange;  // range if using xbreaks

  bool xinvresolsq;
  bool zeroy;  // true if Y axis should start at zero
  bool zeroy_RH;  // true if Y axis should start at zero
  std::vector<scala::Range> ybreaks; // breaks in y-axis, not used
  bool xintegral; // true if x-axis is integral
  bool yintegral; // true if y-axis is integral
  bool yintegral_RH; // true if y-axis is integral
  GraphAxesType axistypes;  // for both x & y

  std::vector<TableGraphPlotline> plotlines;  // plotlines

  std::string formatXbreaks() const;

  std::vector<TablegraphLine> lines; // plain lines
};


//! A class to encapsulate a data table and derived graphs
class TableGraph
//
//! A Table may produce multiple graphs (plots), each of which contains
//! on or more "lines"
//! All lines in a graph share the same x values, but will have different
//! y values, though all with the same scaling (ie x, y1, y2, ...)
//!
//! The assembled plots and data may be returned either as XML data for
//! qloggraph/pimple (XMLformat), or as a plain-text $TABLE suitable for
//! a logfile and loggraph
//!
//! Recommended usage is:-
//!\n  (1) Construct or initialise (\ref TableGraph::init) object with a table title
//!\n  (2) Define each graph (TableGraphPlot) and add it (TableGraph::AddGraph)
//!\n  (3) Store information about all data columns
//!            (TableGraph::StoreColumnFields)
//!\n  (4) Add each line of data (TableGraph::Line)
//!\n  (5) Close table (TableGraph::CloseTable)
//!
//! A graph (TableGraphPlot) is created with attributes defining the axes etc,
//! and a list of lines (TableGraphPlotline), each of which has its own
//! properties
//!
//! The complete formatted table may then be returned as a string,
//! either in XML format (TableGraph::XMLformat)
//!   or as a $TABLE (TableGraph::format)
//!
//! Note that the XML format is a <CCP4Table> block, and will need further
//! wrapping to process in qloggraph/pimple, possibly to encompass multiple
//! <CCP4Table> blocks etc
{
public:
  TableGraph()  {}  // dummy

  //! Store title
  TableGraph(const std::string& Title);
  //! Store title
  void init(const std::string& Title);
  //! Store id string
  void StoreID(const std::string& idstring) {id = idstring;}
  //! Store extra title
  void StoreExtraTitle(const std::string& extraline)
  {extratitle = extraline;}

  //! Add and store a graph (plot)
  void AddGraph(const TableGraphPlot& tgplot);

  //! StoreColumnFields, GetLabels, RawLabels is alternative to ColumnFields
  //! Define format for data table
  // ZeroMark = true to replace zero value with "-"
  void StoreColumnFields(const std::vector<std::string>& Labels,
			 const std::vector<bool>& ZeroMark,
			 const std::string& pformat);

  //! This will probably fail if the number of arguments doesn't match
  //! the format
  std::string Line(const int nc, ...) const;
  //!  Return formatted line, with zeroes by '-' if requested, nc numbers then vector
  std::string Line(const std::vector<double>& val, const int nc, ...) const;
  //!  Return formatted line, with zeroes by '-' if requested
  //! Note that NumberLine, Line & GetLine all append the line
  //!  to internal string array
  // Write nc numbers to returned string, using predefined format
  std::string NumberLine(const int nc, ...) const;

  //! clear line
  void StartLine();
  //! Add integer to next field in line
  void AddToLine(const int& iv);
  //! Add float to next field in line
  void AddToLine(const float& v);
  //! Add double to next field in line
  void AddToLine(const double& v) {AddToLine(float(v));}
  //! replace part of "line" with patch, starting at character firstchar
  void patchLine(const std::string& patch, const int& firstchar);
 //! return line assembled in AddToLine calls (& append to internal array)
  std::string GetLine();


  // Methods to return XML version of table graphs
  //! Return whole table in XML format
  std::string XMLformat() const;
  
  // Methods to generate loggraph formatted $TABLE
  //! Return whole table formatted for loggraph
  std::string format() const;

  //! just format title (DEPRECATED)
  std::string formatTitle() const;

  // Terminate the table
  std::string CloseTable() const;

  //! format label string for loggraph
  std::string GetLabels(const bool& lastmark) const;

  //! DEPRECATED Add a graph, return graph header, loggraph format
  std::string Graph(const std::string& GraphTitle,
		    const std::string& Graphtypestring,
		    const std::vector<int>& columnNumbers);

  //! DEPRECATED Add a graph, return graph header, loggraph format
  std::string Graph(const std::string& GraphTitle,
		    const GraphAxesType& Graphaxestype,
		    const std::vector<int>& columnNumbers);

  //! DEPRECATED Define field widths, & write out labels

  //!   Vector of column labels
  //!  ZeroMark true to replace zeroes by character "-"
  //! pformat is the C-style (printf) format for the table line (with lf)
  //!   This works best with no spaces
  //! Return labels formatted for table output (loggraph format)
  //! if lastmark true [default] add final "$$" after headers
  std::string ColumnFields(const std::vector<std::string>& Labels,
			   const std::vector<bool>& ZeroMark,
			   const std::string& pformat,
			   const bool& lastmark=true); 
  //! Return column labels string
  std::string Labels() const {return labels;} 
  //! Return column labels string, omitting the leader and trailer (if present)
  std::string RawLabels() const;

private:
  class fieldInfo {
  public:
    fieldInfo(const int& fw, const int& dp, const int& tp, const std::string& sf)
      : fieldwidth(fw), dashpos(dp), type(tp), fmt(sf) {}
    int fieldwidth;
    int dashpos;   // position for "-" for missing data, -1 to not use
    int type;      // type: = 0 int; = +1 float
    std::string fmt;
  };

  std::string  title;             // table title
  std::string id;                 // an id string for this graph
  std::string extratitle;         // an optional extra title line

  int ngraphs;                    // number of graphs in table
  std::vector<TableGraphPlot> graphs;     // list of graphs

  // Layout of data table
  int ncolumns;                   // number of columns of data
  std::string prtf_format;        // printf format for data lines
  std::vector<fieldInfo> fields;  // info for each fields
  std::vector<std::string> labelarray; // field labels
  std::string labels;             // all labels
  int lablen; // actual length of real (raw) label string, excluding leader & trailer

  // Formatted data table, from StoreLine entries
  mutable std::vector<std::string> sdatatable;

  std::string  AddInLabel(const std::string& label,
			  const int& ifw, const int& ifd,
			  int& overhang) const;
  mutable std::string line;
  int kfield; // field counter
};

#endif
