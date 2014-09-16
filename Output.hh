//  Output.hh
//
// Written by Airlie McCoy for Phaser
//(c) 2000-2005 Cambridge University Technical Services Ltd
//All rights reserved
#ifndef __Output__Class__
#define __Output__Class__
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdio.h>
#include <boost/timer.hpp>
#include <boost/shared_ptr.hpp>
//#include "PhaserError.hh"

#undef PACKAGE

namespace phaser_io {

enum package {
 NO_PACKAGE,PACKAGE_PHENIX,PACKAGE_CCP4,PACKAGE_CIMR
};

enum outStream {
  SUMMARY, LOGFILE, VERBOSE, DEBUG, LXML, RESULT};

  //class Output : public PhaserError
class Output
{
  public:
    Output();
    void setOutput(const Output &);
    Output(const Output &);
    const Output& operator=(const Output&);
  virtual ~Output() throw() { WriteResult(); verboseStream.close(); logfileStream.close();
      summaryStream.close(); xmlStream.close(); }
 
  private:
    package   PACKAGE;
    bool      TOG_CCP4_SUMMARY,VERBOSE_ON,DEBUG_ON,XMLOUT,SILENT;
    std::string FILEROOT,XMLFILE;
    std::ofstream verboseStream,logfileStream,summaryStream,xmlStream;
    std::string verboseString,logfileString,summaryString,resultString;
    int jumpnum,progressCount,progressBars,progressBarSize,progressChunk;
    std::clock_t start_clock;
    double run_time,elapsed_time;
    int max_line_width;

  protected:
  //    boost::shared_ptr<phenix_out> phenix_out_;

  private:
    void        logOutput(outStream,const std::string&);
    std::string formatMessage(const std::string,int,bool,bool=false);
    std::string tab(int);
    std::string getLine(unsigned,char);
    int         Days(double);
    int         Hrs(double);
    int         Mins(double);
    double      Secs(double);
    void        IncrementRunTime();

  public:
    void openOutputStreams(std::string);
    void startClock();

    void logHeader(outStream);
    void logSectionHeader(outStream,std::string);
  //    void logTrailer(PhaserError=PhaserError());
    void logBlank(outStream);
    void logFlush();
    void logUnderLine(outStream,const std::string);
    void logKeywords(outStream,const std::string);
    void logTab(unsigned,outStream,const std::string,bool add_return=true);
    void logTabPrintf(unsigned,outStream,const char*,...);
    void logEllipsisStart(outStream,const std::string);
    void logEllipsisEnd(outStream);
    void logProgressBarStart(outStream, const int);
    void logProgressBarNext(outStream);
    void logProgressBarEnd(outStream);
    void logProgressBarAgain(outStream);
    void logWarning(outStream,const std::string);
    void logElapsedTime(outStream);

    bool isPackageCIMR() { return (PACKAGE == PACKAGE_CIMR); }
    bool isPackageCCP4() { return (PACKAGE == PACKAGE_CCP4); }
    bool isPackagePhenix() { return (PACKAGE == PACKAGE_PHENIX); }
    void setPackageCIMR();
    void setPackageCCP4();
    void setPackagePhenix();
  //    void setPackagePhenix(boost::shared_ptr<phenix_out> const& out);
    void setVerbose(bool b,bool e) { VERBOSE_ON = b; DEBUG_ON = e; }
    void setSilent(bool b) { SILENT = b; }
    bool Verbose() { return VERBOSE_ON; }
    void setFileroot(std::string f) { FILEROOT = f; }
    void setXmlout(std::string f) { XMLOUT = true; XMLFILE = f; } 
    void unsetXmlout() { XMLOUT = false; XMLFILE = ""; }
    bool doXmlout() { return XMLOUT; } 
    // Write out RESULT stream to stdout if CCP4 & anything in it
    void WriteResult();

    std::string Fileroot() { return FILEROOT; }
    std::string XmlFile() { return XMLFILE; } 
    std::string Package();

    std::string verbose() { return verboseString; }
    std::string logfile() { return logfileString; }
    std::string summary() { return summaryString; }
    std::string XML();

    std::string version_date(); //function code in separate Version.cc file
                                //created by SCons so thate version date updated
                                //each time phaser is compiled
    std::string version_number(); //function code in separate Version.cc file
                                  //created by SCons so thate version number 
                                  //corresponds to release or development
  void SetMaxLineWidth(const int& Width) {max_line_width = Width;}
};


} //phaser

#endif
