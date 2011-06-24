// columnlabels.cpp
//
// Classes:
//  column_labels
//  column_select
//  


#include "columnlabels.hh"
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
  ClipperLabelPair::ClipperLabelPair(const std::string& Xname,
				     const std::string& Dname,
				     const std::string& Label1,
				     const std::string& Label2)
    : xname(Xname), dname(Dname), label1(Label1), label2(Label2)
  {
    // Make MTZ path string
    path = "/"+Xname+"/"+Dname+"/["+Label1;
    if (Label2 != "") {path += ","+Label2;}
    path += "]";
  }
  //--------------------------------------------------------------
  void CheckPairCols(std::string text, const int& col1,const int& col2)
    // Columns col1 & col2 must be both assigned or both unassigned
    // Fail if not
  {
    if ((col1 == 0 && col2 != 0) || (col2 == 0 && col1 != 0))
      Message::message(Message_fatal
       ( "Column pair must be both present or both absent: "+text+"\n")); 
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
      Message::message(Message_fatal("column_labels::lookup_col - list not set up")); 
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
      Message::message(Message_fatal
       ( "Column label not found: "+loglabel+"\n")); 
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
	clipper::String(pcl->second.number) + "\n";
      pcl++;
    }
    return s;
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  column_select::column_select(const MtzIO::column_labels& column_list,
                               col_controls& column_selection)
  {
    // Column assignments
    // compulsory columns
    col_h = column_list.lookup_col("H");
    col_k = column_list.lookup_col("K");
    col_l = column_list.lookup_col("L");
    col_misym = column_list.lookup_col("M_ISYM");
    col_batch = column_list.lookup_col("BATCH");
    col_I = column_list.lookup_col("I");
    col_sigI = column_list.lookup_col("SIGI");
    // optional columns
    col_Ipr = column_list.lookup_col("IPR");
    col_sigIpr = column_list.lookup_col("SIGIPR");
    col_fractioncalc = column_list.lookup_col("FRACTIONCALC");
    col_Xdet = column_list.lookup_col("XDET");
    col_Ydet = column_list.lookup_col("YDET");
    col_Rot = column_list.lookup_col("ROT");
    col_Width = column_list.lookup_col("WIDTH");
    col_LP = column_list.lookup_col("LP");
    col_Mpart = column_list.lookup_col("MPART");
    col_ObsFlag = column_list.lookup_col("FLAG");
    col_BgPkRatio = column_list.lookup_col("BGPKRATIOS");
    col_scale = column_list.lookup_col("SCALE");
    col_sigscale = column_list.lookup_col("SIGSCALE");
    col_time = column_list.lookup_col("TIME");

    // Select profile-fitted (IPR) or integrated (I) column as required
    // Reset column selection if necessary:
    //   if no Ipr column, can't use it
    // and store selection in observation_part class (static method)
    column_selection.SetupColSelection(col_Ipr);
    // Sanity check: I/sigI & Ipr/sigIpr must be both present or both absent
    CheckPairCols("I SIGI",col_I,col_sigI);
    CheckPairCols("IPR SIGIPR",col_Ipr,col_sigIpr);
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
  int FindColumn(const std::vector<ColumnData>& ColumnInfo,
		  const std::string& type)
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
  ClipperLabelPair ProcessLabels(const std::vector<clipper::String>& ColLab,
				 const column_labels& ColumnLabels,
				 bool& IorF)
  // This is for merged files
  // A ColLab element is formatted as "/crystal/dataset/label type"
  //    (clipper format)
  //
  // On entry:
  //  ColLab  clipper column labels from MTZ file
  //  ColumnLabels column labels for "FI" and "SIGFI" if set on input
  //  
  // Returns:
  //  ClipperLabelPair   Clipper path & labels for pair of items eg F, SIGF
  //    SIGF omitted from label if not present
  //  IorF true if column is F, false if J (intensity)
  //  
  {
    //......................................................
    // Column assignments
    int colI = -1;
    int colsigI = -1;
    IorF = false; // true if column is F, false if J (intensity)
    int fail = +1;
    
    std::vector<ColumnData> ColumnInfo(ColLab.size()); // data for each column

    clipper::String M_ISYM_label = "M_ISYM"; // a marker for unmerged file

    for (size_t i=0;i<ColLab.size();i++) {
      //  each element of ColumnInfo contains xname, dname, label, type
      ColumnInfo[i] = ExtractLabelType(ColLab[i]);
      if (ColumnInfo[i].label == M_ISYM_label) {
	return ClipperLabelPair();
      }
    }

    if (ColumnLabels.size() > 0) {
      // We have column label(s) specified for (IorF) & optionally SIG(IorF)
      for (size_t i=0;i<ColumnInfo.size();i++) {
	// Look for a given label
	if (colsigI < 0 && ColumnLabels.Label("FI") == ColumnInfo[i].label) {
	  colI = i;
	  if (ColumnInfo[i].type == "J")
	    IorF = false;
	  else if(ColumnInfo[i].type == "F")
	    IorF = true;
	  else
	    Message::message(Message_fatal
			     ("hkl_merged_list: column is not of type J or F"));
	} else if (colsigI < 0 && ColumnLabels.size() > 1) {
	  // Test SIG label
	  if (ColumnLabels.Label("SIGFI") == ColumnInfo[i].label) {
	    colsigI = i;
	    if (ColumnInfo[i].type != "Q") {
	      Message::message(Message_fatal
			       ("hkl_merged_list: SIG column is not of type Q"));
	      fail = +3;
	    }
	  }
	}
      }  // loop columns in file
    } else {
      // No column labels specified on entry
      // Find first viable column of type J or F
      //  find type J if possible
      if ((colI=FindColumn(ColumnInfo, "J")) >= 0) {
	IorF = false;
	// Is the next column type "Q" ie sigma?
	if (unsigned(colI+1) < ColumnInfo.size() &&
	    ColumnInfo[colI+1].type == "Q") {
	  //  Yes    
	  colsigI = colI+1;
	}
      }
      if (colI < 0) {
	// Search for column or type "F"
	if ((colI=FindColumn(ColumnInfo, "F")) >= 0) {
	  IorF = true;
	  // Is the next column type "Q" ie sigma?
	  if (unsigned(colI+1) < ColumnInfo.size() &&
	      ColumnInfo[colI+1].type == "Q") {
	    //  Yes    
	    colsigI = colI+1;
	  }
	}
      }
      if (colI < 0) {
	//  find type K if possible
	if ((colI=FindColumn(ColumnInfo, "K")) >= 0) {
	  IorF = false;
	  // Is the next column type "M" ie sigma?
	  if (unsigned(colI+1) < ColumnInfo.size() &&
	      ColumnInfo[colI+1].type == "M") {
	    //  Yes    
	    colsigI = colI+1;
	  }
	}
      }
    }  // end no labels specified

    if (colI >= 0 && colsigI > 0) {
      fail = 0; 
    } else if (colI >= 0 && colsigI < 0) {
      // Ior F found but not SIG, look for it
      fail = +2;  //no sigI found
      if (unsigned(colI+1) < ColLab.size()) {
	// I or F is not last column
	// Next column might  be sigma
	if (ColumnInfo[colI+1].type == "Q") {
	  colsigI = colI+1;
	  fail = 0;
	} else {
	  fail = +2;  //no sigI found
	}
      }
    }
  
    if (fail == +1) {
      if (ColumnLabels.size() > 0)
	  Message::message(Message_fatal
			 ("no column found with name "+ColumnInfo[colI].label));
	else
	  Message::message(Message_fatal
			 ("no intensity or F column found"));
      }
    if (fail == +3)
      Message::message(Message_fatal
		       ("sigma column is not of type Q"));

    if (fail == +2) {
      // colI found but not colsigI
      // OK
    } else {
      // Check that FI & SIG columns come from same dataset
      if (!ColumnInfo[colI].SameXDname(ColumnInfo[colsigI])) {
	Message::message(Message_fatal
			 ("IorF column belongs to a different dataset from the SIG column"));
      }
    }

    std::string sigLabel;
    if (colsigI >= 0) {sigLabel =  ColumnInfo[colsigI].label;}
    return ClipperLabelPair(ColumnInfo[colI].xname, ColumnInfo[colI].dname,
			    ColumnInfo[colI].label, sigLabel);

  }  // ProcessLabels


} // namespace MtzIO

