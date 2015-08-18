//
// report_errors.cpp
//

#include "report_errors.hh"
#include "string_util.hh"

// Clipper
#include <clipper/clipper.h>
#include "clipper/clipper-ccp4.h"
using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;

//--------------------------------------------------------------
phaser_io::Output*  ReportErrors::outputp = NULL;
//--------------------------------------------------------------
// just store pointer to output
ReportErrors::ReportErrors(phaser_io::Output& output)
{
  outputp = &output;
}
//--------------------------------------------------------------
void ReportErrors::printFatalError(const std::string& errormessage)
{
  outputp->logTab(0,LOGFILE, "\n**** ERROR ****\n");
  //outputp->logTab(1,LOGFILE, errormessage);
  //outputp->logTab(0,LOGFILE,   "**** ERROR ****\n\n");
  outputp->logTab(0, LXML, StringUtil::MakeXMLwithclass("FatalErrorMessage", errormessage,
                                                        true, "errormessage"));
  Message::message(Message_fatal(errormessage));
}
//--------------------------------------------------------------
void ReportErrors::printWarning(const std::string& warning,
                                const std::string& xmltag,
                                const bool& text)
// print warning message to logfile as TEXT::WARNING, and if xmltag != "", to XML
{
  if (warning == "") {return;}
  if (text) {
    outputp->logTab(0, LOGFILE,
                    std::string("\n$TEXT:Warning:$$ $$\nWARNING: ")+ warning+"\n$$");
  } else {
    outputp->logTab(0, LOGFILE, warning);
  }
  if (xmltag != "") {
    outputp->logTab(0, LXML, StringUtil::MakeXMLwarning(xmltag, warning));
  }
  //  Message::message(Message_warn(warning));
}
