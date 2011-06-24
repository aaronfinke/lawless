//  Output.cpp
//
// Written by Airlie McCoy for Phaser
//(c) 2000-2005 Cambridge University Technical Services Ltd
//All rights reserved
#include <stdarg.h>
#include <math.h>
#include "Output.hh"
#include "Errors.hh"
#include "version.hh"
#include "jiffy.hh"
#include <cassert>
#include <cctype>

namespace phaser_io {


//----
// master controller - every output goes through this function
//----
void Output::logOutput(outStream where, const std::string& text)
{ 
//Everytime something is output the run_time gets updated.
  IncrementRunTime();

//This is the bit that does the output to the output stream IF THEY ARE OPEN
//if they aren't, the output is just stored in the strings with the code above
//until the streams are opened, and then the strings go to the streams

//SUMMARY   output goes to all three output types
//LOG     output goes to logfileString and summaryString
//VERBOSE output goes to debg and summaryString
//DEBUG   output goes to the debug output if DEBUG_ON is true
//LXML     output to XML file if XMLOUT true
//          CCP4 mode only, not stored in string
//RESULT  (CCP4 only) save & output at end, flagged as SUMMARY & Result
  switch (where) {
    case SUMMARY:
      summaryString += text;
      logfileString += text;
      verboseString += text;
      break;
    case LOGFILE:
      logfileString += text;
      verboseString += text;
      break;
    case VERBOSE:
      verboseString += text;
      break;
    case DEBUG:
      if (DEBUG_ON) 
        verboseString += text;
      break;
    case LXML:
      if (XMLOUT)
	{
	  bool open = xmlStream.is_open();
	  if (!open)
	    {xmlStream.open(XMLFILE.c_str());}
	  xmlStream << text << std::flush;
	}
      break;
    case RESULT:
        resultString += text;
    default:
      break;
  }

  if (isPackageCCP4() && !SILENT)
  {
    switch (where) {
      case SUMMARY:
        //summary file is only file ever opened
        //only store the output in the summary string until the stream in open
        // (at which point the stored string is output)
        if (summaryStream.is_open()) summaryStream << text;
        //logfile is sent directly to standard output
        if (!TOG_CCP4_SUMMARY)
        {
          std::cout << "<!--SUMMARY_BEGIN-->" << std::endl;
          TOG_CCP4_SUMMARY = true;
        }
        std::cout << text;
        break;
      case LOGFILE:
        //logfile is sent directly to standard output
        if (TOG_CCP4_SUMMARY)
        {
          std::cout << "<!--SUMMARY_END-->" << std::endl;
          TOG_CCP4_SUMMARY = false;
        }
        std::cout << text;
        break;
      case VERBOSE: 
        //verbose is sent directly to standard output
        if (VERBOSE_ON)
        {
          if (TOG_CCP4_SUMMARY)
          {
            std::cout << "<!--SUMMARY_END-->" << std::endl;
            TOG_CCP4_SUMMARY = false;
          }
          std::cout << text;
        }
        break;
      case DEBUG: 
        //debug is sent to DEBUG file
        if (DEBUG_ON)
        {
	  if (verboseStream.is_open()) verboseStream << text;
	  if (verboseStream.is_open()) std::cout << text;  //^ & to stdout for now
        }
        break;
      case RESULT: 
	// Don't output yet
	break;
      case LXML: 
	break;
    }
  }
  else if (isPackageCIMR())
  {
    switch (where) {
      case SUMMARY:
        //if the three output files are open, send output to them
        //otherwise store in string and wait for file to be opened, and then flush
        if (!SILENT) std::cout << text;
        if (summaryStream.is_open()) summaryStream << text;
        if (logfileStream.is_open()) logfileStream << text;
        if (verboseStream.is_open()) verboseStream << text;
        break;
      case LOGFILE:
        if (!SILENT) std::cout << text;
        if (logfileStream.is_open()) logfileStream << text;
        if (verboseStream.is_open()) verboseStream << text;
        break;
      case VERBOSE: 
        if (VERBOSE_ON && !SILENT) std::cout << text;
        if (verboseStream.is_open()) verboseStream << text;
        break;
      case DEBUG: 
        if (DEBUG_ON)
        {
          if (!SILENT) std::cout << text;
          if (verboseStream.is_open()) verboseStream << text;
        }
        break;
    case LXML:
      break;
    case RESULT:
      break;
    }
  }
}

void Output::logFlush()
{
  std::cout << std::flush;
  if (summaryStream.is_open()) summaryStream << std::flush;
  if (logfileStream.is_open()) logfileStream << std::flush;
  if (verboseStream.is_open()) verboseStream << std::flush;
}

void Output::openOutputStreams(std::string f)
{
  setFileroot(f);
  //open the files for output
  if (isPackageCCP4())
  {
    if (FILEROOT != "")
    {
      //      std::string sumFilename(FILEROOT + ".sum");
      //      summaryStream.open(const_cast<char*>(sumFilename.c_str()));
      //if there is anything stored in the summaryString, output it now
      //      if (summaryStream.is_open()) summaryStream << summaryString;
      //      std::string dbgFilename(FILEROOT + ".dbg");
      //      verboseStream.open(const_cast<char*>(dbgFilename.c_str()));
      DEBUG_ON = false;
    }
  }
  else if (isPackageCIMR())
  {
    if (FILEROOT != "")
    {
      std::string sumFilename(FILEROOT + ".sum");
      std::string logFilename(FILEROOT + ".log");
      std::string dbgFilename(FILEROOT + ".dbg");
      summaryStream.open(const_cast<char*>(sumFilename.c_str()));
      logfileStream.open(const_cast<char*>(logFilename.c_str()));
      verboseStream.open(const_cast<char*>(dbgFilename.c_str()));
      //if there is anything stored in the summaryString, output it now
      if (summaryStream.is_open()) summaryStream << summaryString;
      if (logfileStream.is_open()) logfileStream << logfileString;
      if (verboseStream.is_open()) verboseStream << verboseString;
    }
  }
  logFlush();
}

Output::Output()
{
  PACKAGE = PACKAGE_CCP4;
  FILEROOT = XMLFILE = "";
  XMLOUT = false;
  SILENT = false;
  start_clock = std::clock();
  run_time = elapsed_time = 0;
  jumpnum = 1;
  max_line_width = 85;
  summaryString = logfileString =  verboseString = "";
  resultString = "";
  if (summaryStream.is_open()) summaryStream.close();
  if (logfileStream.is_open()) logfileStream.close();
  if (verboseStream.is_open()) verboseStream.close();
}

std::string Output::Package()
{
  if (isPackageCCP4()) return "PACKAGE_CCP4";
  if (isPackagePhenix()) return "PACKAGE_PHENIX";
  if (isPackageCIMR()) return "PACKAGE_CIMR";
  return "NO_PACKAGE";
}

void Output::setPackageCIMR()
{
  //The style of output for ccp4 is to markup the summary parts of the standard output
  //with <<!--SUMMARY_BEGIN-->> and to have a flag for verbose output that sends more output to standard out
  //This flag is set when the package is set
  PACKAGE = PACKAGE_CIMR;
  TOG_CCP4_SUMMARY = false;
  VERBOSE_ON = false;
  DEBUG_ON = false;
}

void Output::setPackageCCP4()
{
  //The style of output for ccp4 is to markup the summary parts of the standard output
  //with <<!--SUMMARY_BEGIN-->> and to have a flag for verbose output that sends more output to standard out
  //This flag is set when the package is set
  PACKAGE = PACKAGE_CCP4;
  TOG_CCP4_SUMMARY = false;
  VERBOSE_ON = false;
  DEBUG_ON = false;
}

void Output::setPackagePhenix()
{
  PACKAGE = PACKAGE_PHENIX;
  TOG_CCP4_SUMMARY = false;
  VERBOSE_ON = false;
  DEBUG_ON = false;
}


/*
void Output::setOutput(const Output & init)
{
  setPhaserError(init);
  PACKAGE = init.PACKAGE;
  FILEROOT = init.FILEROOT;
  XMLFILE = init.XMLFILE;
  XMLOUT = init.XMLOUT;
  SILENT = init.SILENT;
  TOG_CCP4_SUMMARY = init.TOG_CCP4_SUMMARY;
  VERBOSE_ON = init.VERBOSE_ON;
  DEBUG_ON = init.DEBUG_ON;
  start_clock = init.start_clock;
  run_time = init.run_time;
  elapsed_time = init.elapsed_time;
  summaryString = init.summaryString;
  logfileString = init.logfileString;
  verboseString = init.verboseString;
  max_line_width = init.max_line_width;
  jumpnum = init.jumpnum;
  phenix_out_ = init.phenix_out_;
}
Output::Output(const Output & init)
{ setOutput(init); }

const Output& Output::operator=(const Output& right)
{
  if (&right != this)   // check for self-assignment
    setOutput(right);
  return *this;
}
*/

//===================
//THE OUTPUT FUNCTION
//===================


void Output::IncrementRunTime()
{
//This stops overflows in the number of ticks and prevents
//negative times being output
//Have to check that the integer std::clock() has increased,
//otherwise just get accumulation of small numbers 
  std::clock_t now_clock = std::clock();
  if (now_clock > start_clock)
  {
    run_time += (now_clock-start_clock)/double(CLOCKS_PER_SEC);
    elapsed_time += (now_clock-start_clock)/double(CLOCKS_PER_SEC);
    start_clock = now_clock;
  }
//overflow - clock() has returned to start of long int
//have to miss out the accumulated time between now_clock and start_clock
  else if (now_clock < start_clock) 
  {
    start_clock = now_clock;
  }
}

int Output::Days(double seconds)
{ return int(seconds/24/60/60); }

int Output::Hrs(double seconds)
{ return int(seconds/60/60-Days(seconds)*24); }

int Output::Mins(double seconds)
{ return int(seconds/60-Days(seconds)*24*60-Hrs(seconds)*60); }

double Output::Secs(double seconds)
{ return seconds-Days(seconds)*24*60*60-Hrs(seconds)*60*60-Mins(seconds)*60; }

void Output::startClock()
{ elapsed_time = 0; }

void Output::logElapsedTime(outStream where)
{
  IncrementRunTime(); 
  logTabPrintf(1,where,"CPU Time: %i hrs %i mins %4.2f secs (%2.2f secs)\n",Hrs(elapsed_time),Mins(elapsed_time),Secs(elapsed_time),elapsed_time);
}


//================
//HELPER FUNCTIONS
//================
std::string Output::tab(int depth)
{ if (depth) return std::string(3*depth,' '); else return ""; }

std::string Output::formatMessage(const std::string text,int depth,bool add_return,bool add_amp)
{
  std::string split_text("");
  int line_width(0);
  int start_line(0);
  int last_space(-1);
  std::string amp(add_amp ? " &" : "");
  for (size_t i = 0; i < text.size(); i++)
  {
    if ((std::isspace)(text[i]))
      last_space = i;
    if (i > 0 && text[i-1] == '\n')
      line_width = 1; //don't reset start_line
    else
      line_width++;
    int tab_max_line_width(max_line_width-tab(depth).size()-amp.size());
    if (line_width > tab_max_line_width)
    {
      if (last_space == start_line-1)
      {
        split_text +=  tab(depth) + text.substr(start_line,tab_max_line_width) + amp + '\n';
        start_line += tab_max_line_width;
        last_space = i = start_line-1;
        line_width = 0;
      }
      else
      {
        int length = last_space - start_line;
        if (length)
          split_text +=  tab(depth) + text.substr(start_line,length) + amp + '\n';
        start_line = i = last_space + 1; //reset the start of the line and the text counter
        line_width = 0;
      }
    }
  }
  int length = text.size() - start_line;
  split_text +=  tab(depth) + text.substr(start_line,length);
  if (add_return && split_text.size() > 0 && split_text[split_text.size()-1] != '\n')
    split_text += '\n';

  return split_text;
}

std::string Output::getLine(unsigned len,char c)
{
  std::string line;
  for (unsigned i = 0; i < std::min(len,static_cast<unsigned>(max_line_width)); i++) line += c;
  return line;
}

void Output::logUnderLine(outStream where, const std::string text)
{ 
  logBlank(where);
  logTab(1,where,text);
  logTab(1,where,getLine(text.size(),'-'));
}

void Output::logBlank(outStream where)
{
  logOutput(where,"\n");
  logFlush();
}

void Output::logKeywords(outStream where, const std::string text)
{ logOutput(where,formatMessage(text,0,true,true)); }

void Output::logTab(unsigned t,outStream where, const std::string text,
		    bool add_return)
{ logOutput(where,formatMessage(text,t,add_return)); }

void Output::logEllipsisStart(outStream where, const std::string text)
{ logOutput(where,formatMessage(text+"...",1,true));  }

void Output::logEllipsisEnd(outStream where)
{ logOutput(where,formatMessage("Done",2,true));  }
  
void Output::logWarning(outStream where,const std::string intext)
{
  std::string text = formatMessage(intext,0,true);
  logBlank(where);
  if (isPackageCCP4()) 
    { logTab(0,SUMMARY,"<B><FONT COLOR=\"#FF8800\">"); }
  int max_len(0),j(1);
  for (size_t i = 0; i < text.size(); i++)
  {
    max_len = std::max(j,max_len);
    (text[i] == '\n') ? j = 1 : j++;
  }
  logTab(0,where,getLine(max_len,'-'));
  logTab(0,where,text);
  logTab(0,where,getLine(max_len,'-'));
  if (isPackageCCP4()) 
    { logTab(0,SUMMARY,"</FONT></B>"); }
  logBlank(where);
}

void Output::logProgressBarStart(outStream where, const int length)
{
  progressCount = 0;
  progressBars = 0;
  float progressBarSizeIdeal(max_line_width-tab(1).size()-7); //7 = "|| DONE"
  progressChunk = std::max(static_cast<int>(1),static_cast<int>(ceil(length/progressBarSizeIdeal)));
  progressBarSize = std::max(static_cast<int>(1),static_cast<int>(floor(length/static_cast<float>(progressChunk))));
  logTabPrintf(1,DEBUG,"One progress bar represents %d steps of function\n",progressChunk);
  logTabPrintf(1,DEBUG,"There will be %d progress bars in total\n",progressBarSize);
  std::string percent_bar = "0%";
  for (int i = 0; i < progressBarSize; i++) percent_bar += " ";
  percent_bar +="100%";
  logTab(1,where,percent_bar);
  logOutput(where,tab(1) + "|");
  logFlush();
}

void Output::logProgressBarNext(outStream where)
{
  //Since alot of time is spent on progress bars we need to increment
  //runtime here each time
  IncrementRunTime();

  progressCount++;
  if (progressCount >= progressChunk) {
    logOutput(where,"="); 
    logFlush();
    progressCount -= progressChunk;
    progressBars++;
  }
}

void Output::logProgressBarAgain(outStream where)
{
  std::string percent_bar = "0%";
  for (int i = 0; i < progressBarSize; i++) percent_bar += " ";
  percent_bar +="100%";
  logTab(1,where,percent_bar);
  logOutput(where,tab(1) + "|");
  logFlush();
  for (int i = 0; i < progressBars; i++) logOutput(where,"=");
}

void Output::logProgressBarEnd(outStream where)
{
  logTab(0,where,"=| DONE");
  logBlank(where); 
  logFlush();
}

void Output::logTabPrintf(unsigned t,outStream where, const char* format, ...)
{
  static const std::size_t temp_size = 8192;
  char temp[temp_size];
  temp[temp_size-1] = '\0';
  va_list arglist;
  va_start(arglist,format);
  vsprintf(temp,format,arglist);
  va_end(arglist);
  assert(temp[temp_size-1] == '\0');
  logOutput(where,formatMessage(temp,t,false));
  logFlush();
}

void Output::logSectionHeader(outStream where,std::string header)
{
  std::string jumpstr;
  if (jumpnum == 1) {
    jumpstr = "Jump to: <a href=\"#jump2\">next section</a>\n";
  }
  else {
    jumpstr = "<a name=\"jump" + itos(jumpnum) + "\">"; //<a name="jumpn">
    jumpstr += "Jump to: <a href=\"#jump"+itos(jumpnum+1)+"\">next section</a>"; //Jump to: <a href="#jump(n+1)">next section</a>
    jumpstr += " <a href=\"#jump"+itos(jumpnum-1)+"\">last section</a>\n";
  }
  int len(0),max_len(0);
  for (size_t i = 0; i < header.size(); i++) 
    if (header[i] != '\n') max_len = std::max(max_len,++len);
    else len = 0;
  logTab(0,where,getLine(max_len,'-'));
  logTab(0,where,header); 
  logTab(0,where,getLine(max_len,'-'));
  logBlank(where);
  jumpnum++;
  logFlush();
}


void Output::logHeader(outStream where)
{
  if (PACKAGE == PACKAGE_CCP4)
  {
    time_t start_time;
    std::time(&start_time);
    // ------ suppressed, using standard CCP4 banner calls ----
    /*
    //this bit works out what the header will look like from define parameters set
    //in Version.h
    //Header will take the form
    //Header will take the form
    //#######################################################################");
    //#######################################################################
    //#######################################################################
    //### CCP4 SUITE: PROGRAM_NAME                PROGRAM_VERSION ###
    //#######################################################################
  
    std::string progtag = "### CCP4 SUITE: " + PROGRAM_NAME;
    std::string progversion =  " " + PROGRAM_VERSION + " ###";
    int max_len = max_line_width;
    int min_len = progtag.size() + progversion.size();
    int spacer_len = std::max(max_len - min_len,0);
    std::string  progheader = progtag + getLine(spacer_len,' ') + progversion;
  
    outStream whereHtml(LOGFILE);
    //    logTab(0,whereHtml,"<pre>");
    //    logTab(0,whereHtml,"<B><FONT COLOR=\"#FF0000\">");
    //    logTab(0,whereHtml,"<html><!-- CCP4 HTML LOGFILE -->");
    logBlank(where);
    logTab(0,where,getLine(progheader.size(),'#'));
    logTab(0,where,getLine(progheader.size(),'#'));
    logTab(0,where,getLine(progheader.size(),'#'));
    logTab(0,where,progheader);
    logTab(0,where,getLine(progheader.size(),'#'));
    if (getenv("USER") != NULL)
    logTab(0,where,"User:         " + std::string(getenv("USER")));
    logTab(0,where,"Run time:     " + std::string(std::ctime(&start_time)));
    logTab(0,where,"Version:      " + std::string(PROGRAM_VERSION));
    */
    if (getenv("OSTYPE") != NULL)
    logTab(0,where,"OS type:      " + std::string(getenv("OSTYPE")));
    logTab(0,where,"Release Date: " + std::string(PROGRAM_DATE));
    logBlank(where);
    //    logTab(0,where,"Please reference: Collaborative Computational Project, Number 4. 1994. \"The CCP4 Suite: Programs for Protein Crystallography\". Acta Cryst. D50, 760-763.");
    //    logBlank(where);
    //    logTab(0,whereHtml,"<!--END--></FONT></B>");
    //    logBlank(where);
    // Write XML stuff if activated
    if (XMLOUT)
      {
	std::string time = std::string(std::ctime(&start_time));
	time.resize(time.length()-1); // remove trailing \n
	logTab(0, LXML,
	"<"+PROGRAM_NAME+" version=\""+PROGRAM_VERSION+"\" RunTime=\""+
		time+"\">");
      }
  }
  logFlush();
}
  // ---------------------------------------------------------------------------
  void Output::WriteResult()
  // Write out RESULT stream to stdout if CCP4 & anything in it
  {
    if (PACKAGE == PACKAGE_CCP4) {
      if (resultString.size() > 0) {
	std::cout << "\n<!--SUMMARY_BEGIN--> $TEXT:Result: $$ $$"
		  << std::endl;
	std::cout << resultString;
	std::cout << "\n$$ <!--SUMMARY_END-->" << std::endl;
	resultString = ""; // Clear string
      }
    }
  }


 /*
void Output::logTrailer(PhaserError err)
{
  outStream where(LOGFILE);
  setPhaserError(err);
  if (Failure())
    logWarning(SUMMARY,ErrorName() + " ERROR: " + ErrorMessage());
  std::string trailer = Failure() ? "EXIT STATUS: FAILURE" : "EXIT STATUS: SUCCESS";
  time_t finish_time;
  std::time(&finish_time);
  std::string jumpstr = "<a name=\"jump" + itos(jumpnum)+"\">Jump to: <a href=\"#jump"+ itos(jumpnum-1)+"\">last section</a>\n";
  outStream whereHtml(LOGFILE);
  logTab(0,where,getLine(trailer.size(),'-'));
  logTab(0,where,trailer);
  logTab(0,where,getLine(trailer.size(),'-'));
  logBlank(where);
  IncrementRunTime();
  logTabPrintf(0,where,"CPU Time: %i days %i hrs %i mins %4.2f secs (%2.2f secs)\n",Days(run_time),Hrs(run_time),Mins(run_time),Secs(run_time),run_time);
  logTab(0,where,"Finished: " + std::string(std::ctime(&finish_time))); 
  logBlank(where);
  if (isPackageCCP4())
  {
    logTab(0,whereHtml,"</pre>");
    logTab(0,whereHtml,"</html>");
  }
  openOutputStreams(FILEROOT);
  if (PhaserError::Failure() && XMLOUT)
  {
    FILE* outfile;
    std::string filename = XMLFILE;
    if (outfile = fopen(const_cast<char*>(filename.c_str()), "w"))
    {
      fprintf(outfile,"<phaser version=\"%s\" ostype=\"%s\">\n",Output::version_number().c_str(),getenv("OSTYPE"));
      fprintf(outfile,"%s",Output::XML().c_str());
      fprintf(outfile,"</phaser>\n");
      fclose(outfile);
    }
  }
}
 */
}//namespace phaser
