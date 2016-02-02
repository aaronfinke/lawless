// columnlabels.hh
//
// Classes:
//  column_labels
//  column_select
//  

#ifndef COLUMNLABELS_HEADER
#define COLUMNLABELS_HEADER

#include <string>
#include <map>

// Clipper
#include <clipper/clipper.h>

#include "controls.hh"
#include "hkl_unmerge.hh"
#include "range.hh"

// changed by EK, otherwise clashes with 'OPTIONAL' in Windows
enum col_opt_flag {OF_COMPULSORY,OF_OPTIONAL};

using namespace scala;

namespace MtzIO {  
  typedef std::pair<int,int> IntPair; 
  //======================================================================
  class ColumnData {
    // crystal, dataset, Column label, type
    // xname, dname, label, type
  public:
    ColumnData(){}
    ColumnData(const std::string& Xname, const std::string& Dname,
	       const std::string& Label, const std::string& Type)
      : xname(Xname), dname(Dname), label(Label), type(Type) {}

    // true if two objects have the same xname & dname
    bool SameXDname(const ColumnData& other) const;

    std::string xname;
    std::string dname;
    std::string label; // label in file
    std::string type;  // column type character
  };
  //======================================================================
  class ClipperLabelList {
    // Clipper path & labels for list of items eg /xname/dname/[I(+), SIGI(+), I(-), SIGI(-)]
    // path is "/xname/dname/[label1, label2, ...]"
  public:
    ClipperLabelList(){}
    ClipperLabelList(const std::string& Xname,
		     const std::string& Dname,
		     const std::vector<std::string>& Labels);
    ClipperLabelList(const MtzIO::ColumnData& coldat);

    std::string formatlabels() const;

    std::string xname;
    std::string dname;
    clipper::String path;
    std::vector<std::string> labels;
    bool nosig;
    bool anom;
  };
  //======================================================================
  class ColumnNumberLabel {
  public:
    ColumnNumberLabel(){}
    ColumnNumberLabel(const int& Nc,
		      const std::string& loglab, const std::string& lab,
		      const std::string& typ="")
      :number(Nc), loglabel(loglab) , label(lab), type(typ) {}

    int number;  // column number from 1
    std::string loglabel; // "logical"  label
    std::string label; // actual label in file
    std::string type;  // column type
    Range valuerange;  // range of values, from header
  };
  //======================================================================
  class column_labels {
    // Stores list of columns (label, number) in map container 
    // Column number is 1->Ncolumns
  public:
    column_labels() : nlatticecolumns(0) {setup = false;}   // constructor
    // Store desired column label and COMPULSORY or OPTIONAL flag
    void add (const std::string& loglabel, const col_opt_flag& cflag);

    // True if setup
    bool Setup() const{return setup;}

    // Store labels for (F or I) and its sigma (eg from LABIN) if set
    // Label them as "FI" and "SIGFI"
    void addLabin(const std::string& FIlabel, const std::string& sigFIlabel);

    // return number of columns listed
    int size() const {return columns.size();}
    // return column number (from 0) for requested column
    //  = -1 if not found
    int lookup_col(const std::string& loglabel) const;
    // return actual label, "" if not in range
    std::string Label(const std::string& loglabel) const;
    // return number of lattice columns
    int NlatticeColumns() const {return nlatticecolumns;}

    // Access for setting column numbers
    void start() {pcl = columns.begin(); at_start = true;}  // start iterator
    // Returns next ColumnNumberLabel, or false if beyond end
    bool next(ColumnNumberLabel& CNL);
    // Store updated ColumnNumberLabel at current position
    void Store(ColumnNumberLabel& CNL);
    //! return ColumnNumberLabel for specified column
    ColumnNumberLabel CNL(const std::string& loglabel) const;
    // Store number of latticecolumns
    void SetNlatticeColumns(const int& nlat) {nlatticecolumns = nlat;}

    // format
    std::string format();

  private:
    bool setup;  // true after list is filled from file
    // columns contains a mapping of the "logical" column name to the
    //   actual column number (from 1) and the actual label in the file
    std::map<std::string, ColumnNumberLabel> columns; 
    // iterator for start & next functions
    std::map<std::string, ColumnNumberLabel>::iterator pcl;
    bool at_start;
    int nlatticecolumns;  // number of lattice columns, = 0 if only one
  }; // column_labels
  //===================================================================
  class column_select
  // Column selections for Scala
  {
  public:
    column_select(){}
    column_select(const MtzIO::column_labels& column_list,
		  col_controls& column_selection);

    // Set boolean data flags corresponding to columns present in file
    data_flags DataFlags() const;
   
    int col_h, col_k, col_l, col_misym, col_batch,
      col_I, col_sigI, col_Ipr, col_sigIpr, col_fractioncalc,
      col_Xdet, col_Ydet, col_Rot, col_Width, col_LP, col_Mpart,
      col_ObsFlag, col_BgPkRatio, col_scale, col_sigscale, col_time;
    int col_latnum;
    // extra hkl list for multiple lattice: index to 1st column of group
    std::vector<int> col_lathkl;
    bool scheme2, col_latscale; // true is SCALEn columns present
    int nlatticecolumns;
  };  // column_select
  //======================================================================
  //-----------------------------------------------------------------------------
  class ProcessLabels {
    // This is for merged files
    // A ColLab element is formatted as "/crystal/dataset/label type"
    //    (clipper format)
    //
    // Constructor:
    //  ColLab  clipper column labels from MTZ file
    //  ColumnLabels column labels for "F|I|(+/-)" and "SIGFI" if set on input
    //  IorF true if column is F, false if J (intensity)

  public:
    ProcessLabels(){}

    ProcessLabels(const std::vector<clipper::String>& ColLab,
		  const column_labels& ColumnLabels);

    //  Clipper path & labels for list of items eg F, SIGF
    ClipperLabelList clipperlabellist() const {return clipperlabellist_;}
    // true if columns are F, false if intensity
    bool IorF() const {return IorF_;}
    // true if I+/- or F+/- are present
    bool anom() const {return anom_;}

    bool merged() const {return merged_;}

  private:
    std::vector<ColumnData> ColumnInfo;
    ClipperLabelList clipperlabellist_;
    bool IorF_;
    bool anom_;  // true if anomalous present
    bool merged_;

    //  searches ColumnInfo array for the first column of type "type"
    // returns column number found or -1 if not found
    int FindColumn(const std::string& type) const;

    //  searches ColumnInfo array for column with label
    // returns column number found or -1 if not found
    int FindColumnLabel(const std::string& label) const;

    // return true if ColumnInfo[icol] is of type type 
    bool CheckColumn(const int& icol,
		     const std::string& type) const;

    void failmessage(const std::string& message) const;

  };
} // namespace MtzIO
#endif
