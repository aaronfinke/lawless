//(c) 2000-2007 Cambridge University Technical Services Ltd
//All rights reserved
#include "Errors.h"
#include <cstdio>

namespace phaser {

  error::error(std::string const& msg) throw()
  {
    msg_ = std::string("Program Error: ") + msg;
  }

  error::error(const char* file, long line, std::string const& msg,
               bool internal) throw()
  {
    const char *s = "";
    if (internal) s = " internal";
    char buf[64];
    std::string sfile(file);
    std::string::size_type i = sfile.rfind("/"); //unix file separator
    sfile.erase(sfile.begin(),sfile.begin()+i+1);
    std::string::size_type j = sfile.rfind("\\"); //microsoft file separator
    sfile.erase(sfile.begin(),sfile.begin()+j+1);
    std::sprintf(buf, "%ld", line);
    msg_ =   std::string("Program") + s + " error in source file "
              + sfile + " (line " + buf + ")\n";
    if (msg.size()) msg_ += "*** "+ msg + " ***";
    if (internal) msg_ += "\nPlease email this log file to ccp4@ccp4.ac.uk";
  }

  error::~error() throw() {}

  const char* error::what() const throw()
  {
     return msg_.c_str();
  }

} //phaser
