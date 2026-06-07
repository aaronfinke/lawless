// globalcontrols.hh

#ifndef GLOBAL_AIMLESS_CONTROLS
#define GLOBAL_AIMLESS_CONTROLS

#include "Output.hh"
#include "hkl_datatypes.hh"

using phaser_io::LOGFILE;
using phaser_io::LXML;

namespace scala
{
  //=================================================================
class FlowControl
{
  // Flags to control flow of program
public:
  FlowControl() // Set sensible defaults
    : initialScale(true),
      roughScale(true),
      mainScale(true),
      onlyLambda(false),
      sdcorrectionsinput(false) {}

  void SetOnlyMerge() { // No scaling
    initialScale = false;
    roughScale = false;
    mainScale = false;
  }

  void SetOnlyLambda() { // Wavelength normalization only, no other scaling/outliers
    initialScale = false;  // unity initial scales (ps=bs=ss=ds=1)
    roughScale = false;
    mainScale = false;
    onlyLambda = true;
  }

  void SetAllScales() { // normal scaling
    initialScale = true;
    roughScale = true;
    mainScale = true;
  }

  void SetInitialScale(const bool& flag) {initialScale = flag;}
  void SetRoughScale(const bool& flag) {roughScale = flag;}
  void SetMainScale(const bool& flag) {mainScale = flag;}

  bool OnlyMerge() const { // true if no scaling
    // onlyLambda still builds a (wavelength) scale model, so is not OnlyMerge
    return !initialScale && !roughScale &&  !mainScale && !onlyLambda;
  }
  bool OnlyLambda() const {return onlyLambda;} // true if wavelength-only mode

  bool initialScale;   // true to do initial scaling
  bool roughScale;     // true to do first pass scaling followed by outlier rejection
  bool mainScale;      // true to do main scaling
  bool onlyLambda;     // true for wavelength-normalization-only mode
  bool restore;        // true to restore from dump file
  bool sdoptimise;      // true to optimise SD correction
  bool sdcorrectionsinput;  // true if SD correction parameters explicitly given
};
  //=================================================================

  //--------------------------------------------------------------
  static float ToleranceDefault = 2.0;  // Default value

  class GlobalControls
  {
  public:
    // Null constructor
    GlobalControls() : chiral(CHIRAL), LatTol(ToleranceDefault),
		       ReindexSet(false),
		       assumesameindexing(false) {}
    GlobalControls(const Chirality chiral_in)
      : LatTol(ToleranceDefault), ReindexSet(false),
	assumesameindexing(false)
    {chiral = chiral_in;}

    void set_Chiral(const Chirality& ChiralFlag);
    void set_MinIsig(const double& MinIsig);

    void set_OriginalLattice(const bool& OrigLat) {OriginalLat = OrigLat;}
    void set_Lauegroup(const std::string& LG) {LaueGroup = LG;}
    void set_Spacegroup(const std::string& SG) {SpaceGroup = SG;}
    void set_Reindex(const ReindexOp& H) {Reindx = H; ReindexSet=true;}
    void set_LatticeTolerance(const double& Tol) {LatTol = Tol;}
    void set_SysAbsCheck(const bool& abscheck) {CheckSysAbs = abscheck;}
    void set_TestFirstFile(const bool& TestFirst) {testfirstfile = TestFirst;}
    void set_AssumeSameIndexing(const bool& Assumesameindexing)
      {assumesameindexing = Assumesameindexing;}
    void set_AllowI2(const int& lI2) {allowI2 = lI2;}
    void set_MaxMult(const int& Maxmult) {maxmult = Maxmult;}

    // Access
    Chirality GetChiral() const {return chiral;}
    double GetMinIsig() const {return MinIsigRatio;}

    bool OriginalLattice() const {return OriginalLat;}
    std::string Lauegroup() const {return LaueGroup;}
    std::string Spacegroup() const {return SpaceGroup;}
    ReindexOp Reindex() const {return Reindx;}
    bool IsReindexSet() const {return ReindexSet;}
    double LatticeTolerance() const {return LatTol;}
    bool SysAbsCheck() const {return CheckSysAbs;}
    bool TestFirstFile() const {return testfirstfile;}
    bool AssumeSameIndexing() const {return assumesameindexing;}
    int AllowI2() const {return allowI2;}
    int MaxMult() const {return maxmult;}

  private:
    Chirality chiral;
    double MinIsigRatio;
    bool OriginalLat;    // true to use original lattice
    double LatTol;       // tolerance (degrees) on lattice dimensions
    std::string LaueGroup;   // non-blank to choose Lauegroup
    std::string SpaceGroup;  // non-blank to choose Spacegroup
    ReindexOp Reindx;       //   & reindex operator
    bool ReindexSet;
    bool CheckSysAbs;    // normally true to check systematic absences
    bool testfirstfile;  // test first file for Laue group
    bool assumesameindexing;  // assume all HKLIN files have same indexing
    int  allowI2;        // Allow I2 setting of C2
    int maxmult;         // Maximum multiplicity for scoring
  }; //  GlobalControls
//=================================================================
class OutputControls  
{
public:
  enum OutputType {NONE, MERGED, UNMERGED, BOTH};

  OutputControls();  // defaults

  OutputType  MTZoutputType() const {return mtzoutputtype;} // get
  void SetMTZoutputType(const OutputType& outputtype);  // set
  //! true if merged output wanted
  bool MTZoutputMerged() const;
  //! true if unmerged output wanted
  bool MTZoutputUnmerged() const;

  OutputType  SCAoutputType() const {return scaoutputtype;} // get
  void SetSCAoutputType(const OutputType& outputtype);  // set
  //! true if merged output wanted
  bool SCAoutputMerged() const;
  //! true if unmerged output wanted
  bool SCAoutputUnmerged() const;

  //! return true if merged files needed for MTZ or SCA
  bool Merged() const;

  //! return true if unmerged files needed for MTZ or SCA
  bool UnMerged() const;

  bool  SplitMerged() const {return splitmerged;} //!< return SplitMerged flag
  bool& SplitMerged() {return splitmerged;} //!< set SplitMerged flag

  bool  SplitUnmerged() const {return splitunmerged;} //!< return SplitUnmerged flag
  bool& SplitUnmerged() {return splitunmerged;} //!< set SplitUnmerged flag

  bool  originalHKL() const {return originalhkl;} //!< return originalhkl
  bool& originalHKL() {return originalhkl;} //!< set originalhkl

  //! set filenames from here or from environment
  void SetFilenames(const std::string& hkloutname, const std::string& hkloutunmergedname,
		    const std::string& scaoutname, const std::string& scaoutunmergedname);

  //! return base filename
  std::string Filename() const {return basefilename;}

  //! return output merged MTZ filename, with optional dataset name appended
  std::string Mtzmergedfilename(const std::string& datasetname) const
  {return FileDatasetName(mtzmergedfilename, datasetname);}
  //! return output unmerged MTZ filename
  std::string Mtzunmergedfilename(const std::string& datasetname) const
  {return FileDatasetName(mtzunmergedfilename, datasetname);}
  //! return output merged Scalepackfilename
  std::string Scamergedfilename(const std::string& datasetname) const 
  {return FileDatasetName(scamergedfilename, datasetname);}
  //! return output unmerged Scalepackfilename
  std::string Scaunmergedfilename(const std::string& datasetname) const
  {return FileDatasetName(scaunmergedfilename, datasetname);}
  // write output  file info to XML
  void recordToXML(phaser_io::Output& output) const;

private:
  std::string basefilename; // base file name
  std::string mtzmergedfilename;
  std::string mtzunmergedfilename;
  std::string scamergedfilename;
  std::string scaunmergedfilename;

  // true to split multiple datasets into separate merged MTZ files
  bool splitmerged;
  // true to split multiple datasets into separate unmerged MTZ files
  bool splitunmerged;

  // if originalhkl is true, output original unreduced  hkl,& ISYM = 1
  bool originalhkl;

  OutputType mtzoutputtype;  // MERGED (averaged), UNMERGED for mtz file, or BOTH
  OutputType scaoutputtype;  // MERGED (averaged), UNMERGED for sca file, or BOTH

  // return name is set, otherwise try to get the enviroment value of logname, else return logname
  std::string MakeName(const std::string& name, const std::string& logname) const;

  // return output filename, with optional dataset name appended
  std::string FileDatasetName(const std::string& name,
			      const std::string& datasetname) const;

  std::string typestring(const OutputType& otype) const;

}; // OutputControls
//=================================================================
class ScoreAccept
{
public:
  ScoreAccept() : threshold(-1.0) {};

  ScoreAccept(const double& AcceptanceFraction,
	      const double& MaximumScore);

  bool Accept(const double& score) const;

  bool IfSet() const {return (threshold >= 0.0);}

  double Threshold() const {return threshold*scoremax;}

private:
  double threshold;
  double scoremax;
};


}
#endif
