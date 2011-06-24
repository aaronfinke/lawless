
// globalcontrols.cpp

#include "globalcontrols.hh"

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
    split = true;             // true to split multiple dataset into separate files
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
