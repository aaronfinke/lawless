#include "InputAll.hh"

namespace phaser_io {


  InputAll::InputAll() : OnLine_(false)
  {
    Preprocessor capture(0,0);
    parseCCP4(capture);
  }


  InputAll::InputAll(Preprocessor& capture)  : OnLine_(false)
  { parseCCP4(capture); }

  // Entry to do nothing if "online" (ie if stdin not connected to file)
  InputAll::InputAll(const bool& OnLine)
  {
    Output output;  // dummy
    init(OnLine, false, output);
  }
  // Entry to do nothing if "online" (ie if stdin not connected to file)
  InputAll::InputAll(const bool& OnLine, Output& output)
  {
    init(OnLine, true, output);
  }
  void InputAll::init(const bool& OnLine, const bool& echo, Output& output)
  {
    OnLine_ = OnLine;
      if (! OnLine)
	{
	  Preprocessor capture(0,0);
	  if (echo) {
	    output.logTab(0,phaser_io::LOGFILE, ">>>>> Input command lines <<<<<\n\n");
	    output.logKeywords(phaser_io::LOGFILE, capture.Echo());
	    output.logTab(0,phaser_io::LOGFILE, ">>>>>     End of input    <<<<<\n\n");
	  }
	  parseCCP4(capture);
	}
  }

  InputAll::~InputAll() {}

  //performs the analysis between correlated keyword Input
  void InputAll::Analyse()
  {
    RESO::analyse();
    PARTIALS::analyse();
    EXCLUDE::analyse();
    ANOMALOUS::analyse();
    SCALES::analyse();
    RUNSET::analyse();
    BINS::analyse();
    REJECT::analyse();
    TIE::analyse();
    NAME::analyse();
    REFINE::analyse();
    TITLE::analyse();
    ONLYMERGE::analyse();
    BLANK::analyse();
    SDCORRECTION::analyse();
    INTENSITIES::analyse();
    KEEP::analyse();
    OUTPUT::analyse();
    DUMP::analyse();
    RESTORE::analyse();
    ANALYSIS::analyse();
  }


} // phaser_io
