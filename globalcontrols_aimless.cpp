
// globalcontrols.cpp

#include "globalcontrols.hh"
#include "scala_util.hh"
#include "string_util.hh"

namespace scala
{
  //------------------------------------------------------------
  void GlobalControls::set_Chiral(const Chirality& ChiralFlag)
  {chiral = ChiralFlag;}

  //------------------------------------------------------------
  void GlobalControls::set_MinIsig(const double& MinIsig)
  {MinIsigRatio =  MinIsig;}
  //------------------------------------------------------------
  OutputControls::OutputControls() {
    // defaults
    mtzoutputtype = MERGED;    // write out merged mtz file
    scaoutputtype = NONE;      // no scalepack-format
    splitmerged = true;    // true to split multiple datasets into separate merged files
    splitunmerged = true;  // true to split multiple datasets into separate unmerged files
    originalhkl = false;  // reduced indices in unmerged output
  }
  //------------------------------------------------------------
  void OutputControls::SetMTZoutputType(const OutputType& outputtype) // set
  {
    if (outputtype == NONE && mtzoutputtype == BOTH) {
      // NONE == NOMERGED, ie turn off merged output
      // if it's not BOTH, then leave as is
      mtzoutputtype = UNMERGED;
    } else if (outputtype == MERGED && mtzoutputtype == UNMERGED) {
      mtzoutputtype = BOTH;
    } else if (outputtype == UNMERGED && mtzoutputtype == MERGED) {
      mtzoutputtype = BOTH;
    } else {
      mtzoutputtype = outputtype;
    }
  }
  //------------------------------------------------------------
  void OutputControls::SetSCAoutputType(const OutputType& outputtype) // set
  {
    if (outputtype == NONE && scaoutputtype == BOTH) {
      // NONE == NOMERGED, ie turn off merged output
      // if it's not BOTH, then leave as is
      scaoutputtype = UNMERGED;
    } else if (outputtype == MERGED && scaoutputtype == UNMERGED) {
      scaoutputtype = BOTH;
    } else if (outputtype == UNMERGED && scaoutputtype == MERGED) {
      scaoutputtype = BOTH;
    } else {
      scaoutputtype = outputtype;
    }
  }
  //------------------------------------------------------------
  //! true if merged output wanted
  bool OutputControls::MTZoutputMerged() const
  {
    return (mtzoutputtype == MERGED || mtzoutputtype == BOTH);
  }
  //------------------------------------------------------------
  //! true if unmerged output wanted
  bool OutputControls::MTZoutputUnmerged() const
  {
    return (mtzoutputtype == UNMERGED || mtzoutputtype == BOTH);
  }
  //------------------------------------------------------------
  //! true if merged output wanted
  bool OutputControls::SCAoutputMerged() const
  {
    return (scaoutputtype == MERGED || scaoutputtype == BOTH);
  }
  //------------------------------------------------------------
  //! true if unmerged output wanted
  bool OutputControls::SCAoutputUnmerged() const
  {
    return (scaoutputtype == UNMERGED || scaoutputtype == BOTH);
  }
  //------------------------------------------------------------
  //! return true if merged files needed for MTZ or SCA
  bool OutputControls::Merged() const
  {
    return (MTZoutputMerged() || SCAoutputMerged());
  }
  //------------------------------------------------------------

  //! return true if unmerged files needed for MTZ or SCA
  bool OutputControls::UnMerged() const
  {
    return (MTZoutputUnmerged() || SCAoutputUnmerged());
  }
  //------------------------------------------------------------
  std::string OutputControls::MakeName(const std::string& name, const std::string& logname) const
  // return name is set, otherwise try to get the enviroment value of logname, else return null
  {
    if (name != "") return name;
    std::string fname = "";
    // try logname environment variable
    if (getenv(logname.c_str()) != NULL) {
        fname = std::string(getenv(logname.c_str()));
    }
    return fname;
  }
  //------------------------------------------------------------

  //! set filenames from here or from environment
  void OutputControls::SetFilenames(const std::string& hkloutname, const std::string& hkloutunmergedname,
                                    const std::string& scaoutname, const std::string& scaoutunmergedname)
  {
    mtzmergedfilename = MakeName(hkloutname, "HKLOUT");

    basefilename = FileNameNoExtension(mtzmergedfilename);  // basename from HKLOUT
    if (basefilename == "") {
      mtzmergedfilename =  "HKLOUT";
      basefilename = mtzmergedfilename;
    }

    std::string ext = FileNameExtension(mtzmergedfilename); // save extension if any
    AddFileExtension(mtzmergedfilename, "mtz");  // if not there already

    // Make unmerged MTZ filename, if not set explicitly, in case it is needed
    mtzunmergedfilename = MakeName(hkloutunmergedname, "HKLOUTUNMERGED");
    if (mtzunmergedfilename == "") {
      // add string to base and put extension back if it was present
      mtzunmergedfilename = basefilename+"_unmerged";
      if (ext != "") {mtzunmergedfilename += "."+ext;}
    }
    AddFileExtension(mtzunmergedfilename, "mtz");  // if not there already

    // Make merged Scalepack filename, if not set explicitly, in case it is needed
    scamergedfilename = MakeName(scaoutname, "SCALEPACK");
    if (scamergedfilename == "") {
      scamergedfilename = basefilename;
    }
    AddFileExtension(scamergedfilename, "sca");  // if not there already

    // Make unmerged Scalepack filename, if not set explicitly, in case it is needed
    scaunmergedfilename = MakeName(scaoutunmergedname, "SCALEPACKUNMERGED");
    if (scaunmergedfilename == "") {
      scaunmergedfilename = FileNameNoExtension(scamergedfilename);
      ext = FileNameExtension(scamergedfilename);
      scaunmergedfilename = basefilename+"_unmerged";
      if (ext != "") {scaunmergedfilename += "."+ext;}
    }
    AddFileExtension(scaunmergedfilename, "sca");  // if not there already
  }
  //------------------------------------------------------------
  //------------------------------------------------------------
  // return output filename, with optional dataset name appended
  std::string OutputControls::FileDatasetName(const std::string& name,
                                              const std::string& datasetname) const
  {
    if (datasetname == "") {return name;}
    std::string ext = FileNameExtension(name);
    return FileNameNoExtension(name)+"_"+datasetname+"."+ext;
  }
  //------------------------------------------------------------
  //------------------------------------------------------------
  void OutputControls::recordToXML(phaser_io::Output& output) const
  // write output  file info to XML
  {
    output.logTab(0,LXML,"<OutputFiles>");
    // enum OutputType {NONE, MERGED, UNMERGED, BOTH};
    output.logTab(0,LXML,
                  StringUtil::MakeXMLtag("OutputType",
                  typestring(mtzoutputtype)));
    output.logTab(0,LXML,
                  StringUtil::MakeXMLtag("SCAOutputType",
                  typestring(scaoutputtype)));
    output.logTab(0,LXML,
                  StringUtil::MakeXMLtag("MTZmergedfilename",
                                         mtzmergedfilename));
    output.logTab(0,LXML,
                  StringUtil::MakeXMLtag("MTZunmergedfilename",
                                         mtzunmergedfilename));
    output.logTab(0,LXML,
                  StringUtil::MakeXMLtag("SCAmergedfilename",
                                         scamergedfilename));
    output.logTab(0,LXML,
                  StringUtil::MakeXMLtag("SCAunmergedfilename",
                                         scaunmergedfilename));
    std::string tf = "False";
    if (splitmerged) {tf = "True";}
    output.logTab(0,LXML, StringUtil::MakeXMLtag("SplitMerged", tf));
    tf = "False";
    if (splitunmerged) {tf = "True";}
    output.logTab(0,LXML, StringUtil::MakeXMLtag("SplitUnmerged", tf));
    tf = "False";
    if (originalhkl) {tf = "True";}
    output.logTab(0,LXML, StringUtil::MakeXMLtag("OriginalHKL", tf));
    output.logTab(0,LXML,"</OutputFiles>");
  }
  //------------------------------------------------------------
  std::string OutputControls::typestring(const OutputType& otype) const
  {
    // enum OutputType {NONE, MERGED, UNMERGED, BOTH};
    std::string s = "";
    if (otype == NONE) {
      s = "NONE";
    } else if (otype == MERGED) {
      s = "MERGED";
    } else if (otype == UNMERGED) {
      s = "UNMERGED";
    } else if (otype == BOTH) {
      s = "BOTH";
    }
    return s;
  }
  //------------------------------------------------------------
  ScoreAccept::ScoreAccept(const double& AcceptanceFraction,
                           const double& MaximumScore)
    : threshold(AcceptanceFraction), scoremax(MaximumScore)
  {}
  //------------------------------------------------------------
  bool ScoreAccept::Accept(const double& score) const
  {
    if (threshold < 0.0) return true;
    return (score >= threshold*scoremax);
  }
}
//--------------------------------------------------------------}
