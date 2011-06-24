// fileread.hh
//
// Class for read file and decode it, for eg RESTORE command

#ifndef FILEREAD_HEADER
#define FILEREAD_HEADER

#include <iostream>
#include <fstream>
#include <vector>

class Fileread {
public:
  //! Constructor from input stream object & label for errors
  Fileread(std::ifstream& File, 
	   const std::string& Filename,
	   const std::string& Label);

  //! return filename
  std::string Filename() const {return filename;}

  //! read one tag, leave in internal string
  std::string GetTag() const;

  //! return tag from internal string
  std::string Tag() const {return s;}

  //! throws exception if next field != tag
  void ReadTag(const std::string& tag) const;

  //! skip to after next "{"
  void Skip() const;
  
  //! True if next tag is "}", position after that
  bool CheckEnd() const;

  //! skip over "{}" block, nested if necessary, level 0 is exit
  int SkipSection(const int& level) const;

  //! read one integer
  int Int() const;

  // read one double
  double Double() const;

  //! read integer vector length N
  std::vector<int> IntVec(const int& N) const;

  //! read double vector length N
  std::vector<double> DoubleVec(const int& N) const;


private:
  std::string label;     // label for any error messages
  std::string filename;
  std::ifstream& file;
  mutable std::string s;

  void EOFerror(const std::string& tag) const;

};
#endif
