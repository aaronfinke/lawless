// Written by Airlie McCoy for Phaser
//(c) 2000-2005 Cambridge University Technical Services Ltd
//All rights reserved
#include "Errors.hh"
extern "C" {
  #include <stdio.h>
}

namespace phaser_io {

  InputError::InputError(std::string e,std::string m): message(m),echo(e) { }
std::string InputError::Message() { return message;}
std::string InputError::Echo() { return echo; }

  SyntaxError::SyntaxError(std::string e,std::string m): message(m),echo(e) {  }
std::string SyntaxError::Message() { return message; }
std::string SyntaxError::Echo() { return echo; }

WriteError::WriteError(std::string f,std::string m): filename(f),message(m) {  }
std::string WriteError::Message() { return "File " + filename + ": " + message; }
std::string WriteError::Filename() { return filename; }

PreprocessorError::PreprocessorError(std::string e,std::vector<std::string> m): echo(e)
{
  files_not_found.resize(m.size()); 
  for (size_t i = 0; i < files_not_found.size(); i++)
    files_not_found[i] = m[i];
  if (files_not_found.size() == 1)
    message = "File \"" + files_not_found[0] + "\" could not be opened for reading";
  if (files_not_found.size() > 1)
  {
    message += "Files ";
    int i = 0;
    for (i = 0; i < int(files_not_found.size())-1; i++) message += "\"" +files_not_found[i] + "\" ";
    message += "and \""  + files_not_found[++i] + "\" could not be opened for reading";
  }
}
std::string PreprocessorError::Message() { return message;}
std::string PreprocessorError::partialEcho() { return echo; }

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
    sprintf(buf, "%ld", line);
    msg_ =   std::string("Program") + s + " Error: "
              + file + "(" + buf + ")";
    if (msg.size()) msg_ += std::string(": ") + msg;
  }

  error::~error() throw() {}

  const char* error::what() const throw()
  {
     return msg_.c_str();
  }

} //phaser_io
