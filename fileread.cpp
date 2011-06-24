// fileread.cpp

#include "fileread.hh"

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;

//--------------------------------------------------------------
Fileread::Fileread(std::ifstream& File,
		   const std::string& Filename,
		   const std::string& Label)
  : label(Label), file(File) 
{
  if (!file) {
    filename = Filename;
    clipper::Message::message(Message_fatal
      ("Failed to open "+label+" file "+filename));
  }
}
//--------------------------------------------------------------
std::string Fileread::GetTag() const
// read tag
{
  file >> s;
  if (file.eof()) {EOFerror("");}
  return s;
}
//--------------------------------------------------------------
void Fileread::ReadTag(const std::string& tag) const
// throws exception if next field != tag
{
  GetTag();
  if (s == tag) return;
  clipper::Message::message(Message_fatal
     (label+" tag error: "+s+" != "+tag));
}
//--------------------------------------------------------------
void Fileread::Skip() const
//! skip to after next "{"
{
  while (GetTag() != "{") {}  
}
//--------------------------------------------------------------
bool Fileread::CheckEnd() const
//! True if next tag is "}", position after that
{
  std::string tag = GetTag();
  if (tag == "}") return true; // "}" found
  while (GetTag() != "}") {}  // skip
  return false;
}
//--------------------------------------------------------------
int Fileread::SkipSection(const int& level) const
//! skip over "{}" block, nested if necessary
{
  std::string tag;
  int lev = level;
  while (true) {
    tag = GetTag();
    if (tag == "{") {
      // new inner block found
      lev++;
      if (lev > 1) {
	lev = SkipSection(lev);
      }      
    } else if (tag == "}") {
      return lev-1;
    } 
  }
}
//--------------------------------------------------------------
int Fileread::Int() const
// read one integer
{
  int i;
  file >> i;
  if (file.eof()) {EOFerror("");}
  return i;
}
//--------------------------------------------------------------
double Fileread::Double() const
// read one double
{
  double d;
  file >> d;
  if (file.eof()) {EOFerror("");}
  return d;
}
//--------------------------------------------------------------
//! read integer vector length N
std::vector<int> Fileread::IntVec(const int& N) const
{
  std::vector<int> v;
  for (int i=0;i<N;++i) {
    v.push_back(Int());
  }
  return v;
}
//--------------------------------------------------------------
//! read double vector length N
std::vector<double> Fileread::DoubleVec(const int& N) const
{
  std::vector<double> v;
  for (int i=0;i<N;++i) {
    v.push_back(Double());
  }
  return v;
}
//--------------------------------------------------------------
void Fileread::EOFerror(const std::string& tag) const
{
  clipper::Message::message(Message_fatal
    ("FILEREAD error:"+label+" end of file when looking for "+tag));
}
//--------------------------------------------------------------
