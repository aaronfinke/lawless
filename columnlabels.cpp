// columnlabels.cpp
//
// Classes:
//  column_labels
//  column_select
//


#include "columnlabels.hh"
#include "string_util.hh"
#include "report_errors.hh"

#include <assert.h>
#define ASSERT assert

namespace MtzIO {
  //--------------------------------------------------------------
  bool ColumnData::SameXDname(const ColumnData& other) const
  // true if two objects have the same xname & dname
  {
    return ((xname == other.xname) && (dname == other.dname));
  }
  //--------------------------------------------------------------
  ClipperLabelList::ClipperLabelList(const std::string& Xname,
                     const std::string& Dname,
                     const std::vector<std::string>& Labels)
    : xname(Xname), dname(Dname)
  {
    ASSERT (Labels.size() > 0);
    labels = Labels;
    // if name contains slashes, replace with '*'
    xname = fudgeName(xname);
    dname = fudgeName(dname);

    // Make MTZ path string
    path = "/"+xname+"/"+dname+"/[";
    int n = 0;
    for (size_t k=0; k<Labels.size(); k++) {
      if (Labels[k] != "") {
        if (n++ > 0) {path += ",";}
        path += Labels[k];}
    }
    path += "]";
  }
  //--------------------------------------------------------------
  ClipperLabelList::ClipperLabelList(const MtzIO::ColumnData& coldat)
    : xname(coldat.xname), dname(coldat.dname)
  {
    // if name contains slashes, replace with '*'
    xname = fudgeName(xname);
    dname = fudgeName(dname);
    labels.assign(1, coldat.label);
    // Make MTZ path string
    path = "/"+xname+"/"+dname+"/["+labels[0]+"]";
  }
  //--------------------------------------------------------------
  std::string ClipperLabelList::fudgeName(const std::string& name) const
  // if name contains slashes, replace with '*'
  {
    std::string newname = name;
    if (name.find("/") != std::string::npos) {
      newname = "*";
    }
    return newname;
  }
  //--------------------------------------------------------------
  std::string ClipperLabelList::formatlabels() const
  {
    std::string s;
    for (size_t i=0; i<labels.size(); i++) {
      s += labels[i];
      if (i < labels.size()-1) {s += ", ";}
    }
    return s;
  }
  //--------------------------------------------------------------
  void CheckPairCols(std::string text, const int& col1,const int& col2)
    // Columns col1 & col2 must be both assigned or both unassigned
    // Fail if not
  {
    if ((col1 == 0 && col2 != 0) || (col2 == 0 && col1 != 0))
      ReportErrors::printFatalError
        ( "Column pair must be both present or both absent: "+text+"\n");
  }
  //--------------------------------------------------------------
  void column_labels::add(const std::string& loglabel,
                          const col_opt_flag& cflag)
  // Add a column label entry into list, initialised to
  // column -1 if COMPULSORY [default] or 0 if OPTIONAL
  {
    int cfl;
    if (cflag == OF_OPTIONAL) {  // EK: renamed due to Windows clash
      cfl = 0;
    } else {
      cfl = -1;
    }
    columns[loglabel] =  ColumnNumberLabel(cfl, loglabel, loglabel);
  }
  //--------------------------------------------------------------
  // Store labels for (F or I) and its sigma (eg from LABIN)
  // Label them as "FI" and "SIGFI"
  void column_labels::addLabin(const std::string& FIlabel, const std::string& sigFIlabel)
  {
    if (FIlabel != "") {
      columns["FI"] = ColumnNumberLabel(-1, "FI", FIlabel);
    }
    if (sigFIlabel != "") {
      columns["SIGFI"] = ColumnNumberLabel(-1, "SIGFI", sigFIlabel);
    }
  }
  //--------------------------------------------------------------
  int column_labels::lookup_col(const std::string& loglabel) const
  // return column number (from 0) for requested column, else -1
  {
    if (!setup)
      ReportErrors::printFatalError("column_labels::lookup_col - list not set up");
    std::map<std::string, ColumnNumberLabel>::const_iterator p = columns.find(loglabel);
    if (p == columns.end()) {
      // Label not found
      return -1;
    } else {
      // label found, return index-1 (from 0)
      return (p->second).number - 1;
    }
  }
  //--------------------------------------------------------------
  // return actual label, "" if not in range
  std::string column_labels::Label(const std::string& loglabel) const
  {
    std::map<std::string, ColumnNumberLabel>::const_iterator p = columns.find(loglabel);
    if (p == columns.end()) {
      // Label not found
      return "";
    } else {
      // label found, return actual label
      return (p->second).label;
    }
  }
  //--------------------------------------------------------------
  bool column_labels::next(ColumnNumberLabel& CNL)
  // Returns next ColumnNumberLabel, or false if beyond end
  {
    if (!at_start) {pcl++;}  // increment iterator if not at start
    at_start = false;
    if (pcl == columns.end()) {
      start();
      return false;
    }
    CNL = pcl->second;
    return true;
  }
  //--------------------------------------------------------------
  void column_labels::Store(ColumnNumberLabel& CNL)
  // Store updated ColumnNumberLabel at current position
  {
    pcl->second = CNL;
    setup = true;
  }
  //--------------------------------------------------------------
  //! return ColumnNumberLabel for specified column
  ColumnNumberLabel column_labels::CNL(const std::string& loglabel) const
  {
    std::map<std::string, ColumnNumberLabel>::const_iterator
      p = columns.find(loglabel);
    if (p == columns.end()) {
      // Label not found
      ReportErrors::printFatalError( "Column label not found: "+loglabel+"\n");
    }
    return p->second;
  }
 //--------------------------------------------------------------
  std::string column_labels::format()
  {
    std::string s;
    start();
    while (pcl != columns.end()) {
      s += pcl->second.loglabel + " : " + pcl->second.label +
        clipper::String(pcl->second.number);
      if (pcl->second.number == 0) {
        s += " Absent column";
      } else {
        Range range = pcl->second.valuerange;
        if (range.Valid()) {
          s += " Range: " + clipper::String(range.min())
            + " - " + clipper::String(range.max());
        } else {
          s += " Invalid value range";
        }
      }
      s += "\n";
      pcl++;
    }
    return s;
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  column_select::column_select(const MtzIO::column_labels& column_label_list,
                               col_controls& column_selection)
  {
    // Column assignments
    // compulsory columns
    col_h = column_label_list.lookup_col("H");
    col_k = column_label_list.lookup_col("K");
    col_l = column_label_list.lookup_col("L");
    col_misym = column_label_list.lookup_col("M_ISYM");
    col_batch = column_label_list.lookup_col("BATCH");
    col_I = column_label_list.lookup_col("I");
    col_sigI = column_label_list.lookup_col("SIGI");
    // optional columns
    col_Ipr = column_label_list.lookup_col("IPR");
    col_sigIpr = column_label_list.lookup_col("SIGIPR");
    col_fractioncalc = column_label_list.lookup_col("FRACTIONCALC");
    col_Xdet = column_label_list.lookup_col("XDET");
    col_Ydet = column_label_list.lookup_col("YDET");
    col_Rot = column_label_list.lookup_col("ROT");
    col_Width = column_label_list.lookup_col("WIDTH");
    col_LP = column_label_list.lookup_col("LP");
    col_Mpart = column_label_list.lookup_col("MPART");
    col_ObsFlag = column_label_list.lookup_col("FLAG");
    col_BgPkRatio = column_label_list.lookup_col("BGPKRATIOS");
    col_scale = column_label_list.lookup_col("SCALE");
    col_sigscale = column_label_list.lookup_col("SIGSCALE");
    col_time = column_label_list.lookup_col("TIME");

    // Select profile-fitted (IPR) or integrated (I) column as required
    // Reset column selection if necessary:
    //   if no Ipr column, can't use it
    // and store selection in observation_part class (static method)
    column_selection.SetupColSelection(col_Ipr);
    // Sanity check: I/sigI & Ipr/sigIpr must be both present or both absent
    CheckPairCols("I SIGI",col_I,col_sigI);
    CheckPairCols("IPR SIGIPR",col_Ipr,col_sigIpr);

    // extra columns for multiple lattices, col_latnum, col_lathkl
    col_latnum = -1;
    nlatticecolumns = column_label_list.NlatticeColumns();
    col_lathkl.clear();
    col_latscale = false;
    if (nlatticecolumns > 0) {
      col_latnum = column_label_list.lookup_col("LATTNUM");
      // If this is a scheme 2 file, then a LATTNUM1 must be present
      scheme2 = (column_label_list.lookup_col("LATTNUM1") >= 0);
      std::string label;
      int colnum;
      for (int i=0;i<nlatticecolumns;++i) {
        // Scheme 2, LATTNUMn
        if (scheme2) {
          label = "LATTNUM"+StringUtil::Strip(StringUtil::itos(i+1));
          colnum = column_label_list.lookup_col(label);  // column number from 0
          label = "SCALE"+StringUtil::Strip(StringUtil::itos(i+1));
          if (column_label_list.lookup_col(label) >= 0) {
            col_latscale = true;
          } else if (col_latscale) {
            ReportErrors::printFatalError( "Inconsistent SCALEn columns\n");
          }
        } else {
          // Scheme 2, label should be Hn where n is 1->9
          label = "H"+StringUtil::Strip(StringUtil::itos(i+1));
          colnum = column_label_list.lookup_col(label);  // column number from 0
        }
        if (colnum >= 0) {
          col_lathkl.push_back(colnum);  // column number for "Hn"
        }
      } // end loop nlatticecolumns
      ASSERT (int(col_lathkl.size()) == nlatticecolumns);
    }
  }
  //--------------------------------------------------------------
  data_flags column_select::DataFlags() const
  // Set boolean data flags corresponding to columns present in file
  {
    data_flags flags;  // Constructor sets compulsory columns true, other false
    // Optional
    if (col_Ipr >= 0) flags.is_Ipr = true;
    if (col_sigIpr >= 0) flags.is_sigIpr = true;
    if (col_fractioncalc >= 0) flags.is_fractioncalc = true;
    if (col_Xdet >= 0) flags.is_Xdet = true;
    if (col_Ydet >= 0) flags.is_Ydet = true;
    if (col_Rot >= 0) flags.is_Rot = true;
    if (col_Width >= 0) flags.is_Width = true;
    if (col_LP >= 0) flags.is_LP = true;
    if (col_Mpart >= 0) flags.is_Mpart = true;
    if (col_ObsFlag >= 0) flags.is_ObsFlag = true;
    if (col_BgPkRatio >= 0) flags.is_BgPkRatio = true;
    if (col_scale >= 0) flags.is_scale = true;
    if (col_sigscale >= 0) flags.is_sigscale = true;
    if (col_time >= 0) flags.is_time = true;
    if (col_latnum >= 0) {
      flags.is_latnum = true;
      flags.n_latinfo = nlatticecolumns;
      if (scheme2) {
        // scheme 2 multiple lattices
        if (col_latscale) flags.is_latscale = true;
        flags.is_latinfo = true;
        flags.is_lathkl = false;
      } else { // scheme 1
        if (col_lathkl.size() > 0) flags.is_lathkl = true;
        flags.is_latscale = false;
        flags.is_latinfo = false;
      }
    }
    return flags;
  }
  //--------------------------------------------------------------
  ColumnData ExtractLabelType(const clipper::String& collab)
  // Extract column label & type from clipper column label formatted as
  //    "/crystal/dataset/label type"
  //
  // Returns xname, dname, label, type

  {
    std::vector<clipper::String> substrings = collab.split("/");
    ASSERT (substrings.size() == 3);
    // column label, column type
    std::vector<clipper::String> LabelTypes = substrings.back().split(" ");
    ASSERT (LabelTypes.size() == 2);
    return ColumnData(substrings[0], substrings[1], LabelTypes[0], LabelTypes[1]);
  }
  //--------------------------------------------------------------
  int ProcessLabels::FindColumn(const std::string& type) const
  //  searches ColumnInfo array for the first column of type "type"
  // returns column number found or -1 if not found
  {
    for (size_t i=0;i<ColumnInfo.size();i++) {
      if (ColumnInfo[i].type == type) {
        return i;
      }
    }
    return -1;
  }
  //--------------------------------------------------------------
  int ProcessLabels::FindColumnLabel(const std::string& label) const
  //  searches ColumnInfo array for column with label
  // returns column number found or -1 if not found
  {
    for (size_t i=0;i<ColumnInfo.size();i++) {
      if (ColumnInfo[i].label == label) {
        return i;
      }
    }
    return -1;
  }
  //--------------------------------------------------------------
  bool ProcessLabels::CheckColumn(const int& icol,
                                  const std::string& type) const
  // return true if ColumnInfo[icol] is of type type
  {
    if ((unsigned(icol) < ColumnInfo.size()) &&
        (ColumnInfo[icol].type == type)) {return true;}
    return false;
  }
  //--------------------------------------------------------------
  ProcessLabels::ProcessLabels(const std::vector<clipper::String>& ColLab,
                                     const column_labels& ColumnLabels)
  // This is for merged files
  // A ColLab element is formatted as "/crystal/dataset/label type"
  //    (clipper format)
  //
  //  ColLab  clipper column labels from MTZ file
  //  ColumnLabels column labels for "FI" and "SIGFI" if set on input
  //
  {
    //......................................................
    // Column assignments
    std::vector<int> selectedcolnums;
    IorF_ = false; // true if column is F, false if J (intensity)
    anom_ = false; // no anomalous

    ColumnInfo.resize(ColLab.size()); // data for each column

    clipper::String M_ISYM_label = "M_ISYM"; // a marker for unmerged file

    merged_ = true;
    for (size_t i=0;i<ColLab.size();i++) {
      //  each element of ColumnInfo contains xname, dname, label, type
      ColumnInfo[i] = ExtractLabelType(ColLab[i]);
      if (ColumnInfo[i].label == M_ISYM_label) {
        merged_ = false;
        return;
      }
    }

    std::string xname;
    std::string dname;

    int col1 = -1;

    std::vector<std::string> selectedtypes;
    int ncol;
    bool allowmissing = false;

    if (ColumnLabels.size() > 0) {  //   - - - - Specified columns  FIXME
      // We have column label(s) specified for (IorF) & optionally SIG(IorF)
      col1 = FindColumnLabel(ColumnLabels.Label("FI"));
      if (col1 < 0) { // label not found
        ReportErrors::printFatalError
          ("no column found with name "+ColumnLabels.Label("FI"));
      }
      std::string type = ColumnInfo[col1].type;
      selectedtypes.push_back(type);
      if (type == "K") {
        IorF_ = false;
        anom_ = true;
        ncol = 4;
        // found a I+/- column, the next 3 columns should be of type M,K,M
        selectedtypes.push_back("M");
        selectedtypes.push_back("K");
        selectedtypes.push_back("M");
      } else if (type == "J") {
        IorF_ = false;
        anom_ = false;
        ncol = 2;
        // found a I column, the next column should be of type Q
        selectedtypes.push_back("Q");
      } else if (type == "G") {
        IorF_ = true;
        anom_ = true;
        ncol = 4;
        // F+ and F-, the next 3 columns should be of type L,G,L
        selectedtypes.push_back("L");
        selectedtypes.push_back("G");
        selectedtypes.push_back("L");
      } else if (type == "F") {
        IorF_ = true;
        anom_ = false;
        ncol = 2;
        // F, the next column should be of type Q, but allow missing
        selectedtypes.push_back("Q");
        allowmissing = true;
      } else {
        ReportErrors::printFatalError("chosen column is not intensity or F");
      }
    } else {       // - - - - - No column labels specified on entry
      // Find first viable column of type K, J (intensities), or G, F (amplitudes),
      //  in this order of preference
      //   K   I+/-   sigI type M
      //   J   Imean  sigI type Q
      //   G   F+/-   sigF type L
      //   F   Fmean  sigF type Q
      //  find type J if possible
      if ((col1 = FindColumn("K")) >= 0) {
        // found a I+/- column, the next 3 columns should be of type M,K,M
        ncol = 4;
        selectedtypes.push_back("K");
        selectedtypes.push_back("M");
        selectedtypes.push_back("K");
        selectedtypes.push_back("M");
        IorF_ = false;
        anom_ = true;
      } else if ((col1 = FindColumn("J")) >= 0) {
        ncol = 2;
        // found a I column, the next column should be of type Q
        selectedtypes.push_back("J");
        selectedtypes.push_back("Q");
        IorF_ = false;
        anom_ = false;
      } else if ((col1 = FindColumn("G")) >= 0) {
        ncol = 4;
        // F+ and F-, the next 3 columns should be of type L,G,L
        selectedtypes.push_back("G");
        selectedtypes.push_back("L");
        selectedtypes.push_back("G");
        selectedtypes.push_back("L");
        IorF_ = true;
        anom_ = true;
      } else if ((col1 = FindColumn("F")) >= 0) {
        ncol = 2;
        // F, the next column should be of type Q, but allow missing
        selectedtypes.push_back("F");
        selectedtypes.push_back("Q");
        allowmissing = true;
        IorF_ = true;
        anom_ = false;
      } else {
        ReportErrors::printFatalError("no intensity or F column found");
      }
    }  // end no labels specified

    if (col1 < 0) {
      ReportErrors::printFatalError("no intensity or F column found");
    }

    bool nosig = false;
    for (int i=0;i<ncol;++i) {
      if (col1+i < ColumnInfo.size()) {
        // check types of columns
        if (CheckColumn(col1+i, selectedtypes[i])) {
          // column is of correct type
          selectedcolnums.push_back(col1+i);
        } else {
          if (allowmissing) {
            // column not there, but allowed to be missed
            nosig = true; // missing sigma
          } else {
            failmessage(std::string("Column is of wrong type: ")+
                        ColumnInfo[col1+i].type);
          }
        }
      }
    }

    xname = ColumnInfo[col1].xname;
    dname = ColumnInfo[col1].dname;
    std::vector<std::string> selectedlabels;
    for (size_t k=0; k<selectedcolnums.size(); k++) {
      // Check that selected columns come from same dataset
      if (!ColumnInfo[selectedcolnums[0]].
          SameXDname(ColumnInfo[selectedcolnums[k]])) {
        ReportErrors::printFatalError
                ("Selected columns belong to different datasets");
      }
      selectedlabels.push_back(ColumnInfo[selectedcolnums[k]].label);
    }
    clipperlabellist_ = ClipperLabelList(xname, dname, selectedlabels);
    clipperlabellist_.nosig = nosig;  // flag for no sigma column
    clipperlabellist_.anom = anom_;
  }  // ProcessLabels
  //--------------------------------------------------------------
  void ProcessLabels::failmessage(const std::string& message) const
  {
    ReportErrors::printFatalError("ProcessLabels: "+message);
  }
} // namespace MtzIO
