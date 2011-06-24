//
//   tablegraph.hh
//
//!  TableGraph class for writing table for loggraph
//!
//! Graph syntax (ccp4 6.1)
//! ============
//! 
//!  $TABLE :table name:
//!  $GRAPHS :graph1 name:graphtype:column_list:
//!          :graph2 name:graphtype:column_list:
//!          :graph 3 ...: ... $$
//!  column1_name column2_name ... $$ any_characters $$
//!   numbers $$
//!
//!  graphtype is
//!  
//!  A[UTO]
//!     for fully automatic scaling (e.g. ... :A:1,2,4,5:)
//!  N[OUGHT]
//!    for automatic y coordinate scaling, where y lowest limit is 0
//!    (e.g. ... :N:1,2,4,5:)
//!  XMIN|XMAXxYMIN|YMAX
//!    for user defined scaling where XMIN ... are axis limits
//!    (e.g. ... :0|100x-1|1:1,2,4,5:)
//!  
//!  

#ifndef TABLEGRAPH_HEADER
#define TABLEGRAPH_HEADER

#define ASSERT assert
#include <assert.h>

#include <string>
#include <vector>
#include "range.hh"

class GraphAxesType {
  //! graph axis type for x and y axes
  //
  //!  Loggraph options are:
  //!  A[UTO]
  //!     for fully automatic scaling of y axis (e.g. ... :A:1,2,4,5:)
  //!  N[OUGHT]
  //!    for automatic y coordinate scaling, where y lowest limit is 0
  //!    (e.g. ... :N:1,2,4,5:)
  //!  XMIN|XMAXxYMIN|YMAX
  //!    for user defined scaling where XMIN ... are axis limits
  //!    (e.g. ... :0|100x-1|1:1,2,4,5:)
  //!
  //!  these may be extended in future
public:
  enum GraphType {AUTO_Y, NOUGHT_Y, XY_SPECIFIED};
  //! construct as default
  GraphAxesType() : graphtype(AUTO_Y){} // default
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
  //! return type
  GraphType Graphtype() const {return graphtype;}
  //! return formatted for graph
  std::string FormatType() const;

private:
  GraphType graphtype;
  scala::Range xrange;
  scala::Range yrange;

  // Fix Yrange to be "sensible"
  // If ZeroY true, then y range should start at 0
  void FixYrange(const bool& ZeroY);
};

class TableGraph
{
public:
  TableGraph()
  {}  // dummy

  //! Store title
  TableGraph(const std::string& Title);
  //! Store title
  void init(const std::string& Title);

  //! Return formatted title
  std::string formatTitle() const;

  //! Add a graph, return graph header
  std::string Graph(const std::string& GraphTitle,
		    const std::string& Graphtypestring,
		    const std::vector<int>& columnNumbers);

  //! Add a graph, return graph header
  std::string Graph(const std::string& GraphTitle,
		    const GraphAxesType& Graphaxestype,
		    const std::vector<int>& columnNumbers);

  //! Define field widths, & write out labels
  //!   Vector of column labels
  //!  ZeroMark true to replace zeroes by character "-"
  //! pformat is the C-style (printf) format for the table line
  //!   This works best with no spaces
  //! Return labels formatted for table output
  //! if lastmark true [default] add final "$$" after headers
  std::string ColumnFields(const std::vector<std::string>& Labels,
			   const std::vector<bool>& ZeroMark,
			   const std::string& pformat,
			   const bool& lastmark=true); 

  //! Column labels string
  std::string Labels() const {return labels;} 
  //! Column labels string, omitting the leader and trailer (if present)
  std::string RawLabels() const;

  //! Write nc numbers to returned string, using predefined format
  //! This will probably fail if the number of arguments doesn't match
  //! the format
  std::string NumberLine(const int nc, ...) const;
  //
  //!  Return formatted line, with zeroes by '-' if requested
  std::string Line(const int nc, ...) const;

  //! clear line
  void StartLine();
  //! Add integer to next field in line
  void AddToLine(const int& iv);
  //! Add float to next field in line
  void AddToLine(const float& v);
  //! return line assembled in AddToLine calls
  std::string GetLine();

  //! Return the table terminating line
  std::string CloseTable() const;

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

  int ngraphs;                    // number of graphs in table
  int ncolumns;                   // number of columns of data

  std::string prtf_format;           // printf format for data lines

  std::vector<fieldInfo> fields;

  std::string labels;
  int lablen; // actual length of real (raw) label string, excluding leader & trailer

  std::string  AddInLabel(const std::string& label,
			  const int& ifw, const int& ifd,
			  int& overhang) const;
  mutable std::string line;
  int kfield; // field counter

  static const std::string LABELLEADER;
  static const std::string LABELTRAILER;
  static const std::string LABELFINAL;
};



#endif
