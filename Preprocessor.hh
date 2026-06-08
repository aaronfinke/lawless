//  Preprocessor.hh
//
// Written by Airlie McCoy for Phaser
//(c) 2000-2005 Cambridge University Technical Services Ltd
//All rights reserved
#ifndef __Preprocessor__Class__
#define __Preprocessor__Class__

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include "Errors.hh"

namespace phaser_io {

  class Preprocessor
  {
  private:
    std::string echo;
 
  public:
    Preprocessor();
    Preprocessor(std::string);
    // Usually call with argc, argv
    // Optional 3rd argument = false to not read stdin
    Preprocessor(int, char*[],bool ReadInput=true);
    ~Preprocessor();
    std::string Echo();
    std::string stoup(const std::string & str);
    void addLine(std::string);
    void deleteLine(std::string);
    void addFile(std::string); 
    std::vector<std::string> getEndKeys(); 
  };

} //phaser_io

#endif
