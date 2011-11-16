// file_util.cpp

#include <stdlib.h>
#include "file_util.hh"

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;


FILE* OpenFile(const std::string& Filename, const bool& Write,
	       const bool Binary)
{
  std::string filename = Filename;
  if (getenv(filename.c_str()) != NULL)
    {filename = std::string(getenv(Filename.c_str()));}
  std::string rw = "r";
  if (Write) {rw = "w";}
  if (Binary) rw+= "b";
  FILE* fp = fopen(filename.c_str(), rw.c_str());
  if (fp == NULL) {
    Message::message(Message_fatal
		     ("OpenFile: cannot open file "+filename));
  }
  return fp;
}

void WriteToFile(const std::string& filename,
			  const std::string& text)
// write text to filename
{
  FILE* dfile = OpenFile(filename, true);
  if (dfile == NULL) {
    clipper::Message::message(Message_fatal
			      ("Failed to open file "+filename));
  }
  fprintf(dfile, "%s", text.c_str());
  fclose(dfile);
}


