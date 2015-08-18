//
// report_errors.hh
//

#ifndef REPORT_ERRORS_HEADER
#define REPORT_ERRORS_HEADER


// * Old stuff * // routines from printing.hh
//--------------------------------------------------------------
///void printWarning(const std::string& warning, const std::string& xmltag,
///		  const bool& text, phaser_io::Output& output);
// print warning message to logfile as TEXT::WARNING, and if xmltag != "", to XML
//--------------------------------------------------------------
///void PrintError(const std::string& ErrorMessage,
///		phaser_io::Output& output);

#include <cstring>
#include "Output.hh"
using phaser_io::LOGFILE;
using phaser_io::LXML;

class ReportErrors {
  // report errors or warnings to LOGFILE and maybe LXML
public:
  ReportErrors() {}

  // just store pointer to output
  ReportErrors(phaser_io::Output& output);

  // print fatal error to logfile and to XML as FatalErrorMessage
  static void printFatalError(const std::string& errormessage);

  // print warning to logfile, as TEXT::WARNING if text = true,
  // and to XML if xmltag != ""
  static void printWarning(const std::string& warning,
			   const std::string& xmltag,
			   const bool& text=true);

private:
  static phaser_io::Output*  outputp;  // pointer to Output object

};

#endif
