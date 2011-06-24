//  Preprocessor.cpp
//
// Written by Airlie McCoy for Phaser
//(c) 2000-2005 Cambridge University Technical Services Ltd
//All rights reserved
#include "Preprocessor.hh"
#include <cctype>
#include <stdio.h>

namespace phaser_io {

Preprocessor::Preprocessor() 
{}

Preprocessor::Preprocessor(std::string e) 
{ echo = e; }
    
Preprocessor::Preprocessor(int argc, char* argv[], bool ReadInput)
{
  std::istream& input_stream = std::cin;
  std::vector<std::string> files_not_found;
  echo = "";
  //turns the command line into keyworded input
  // assuming KEYWORD VALUE pairs on command line
  // arguments beginning with "-" are assumed to be unpaired
  for(int i=1; i<argc; i++)
  {
    echo += std::string(argv[i]) + " " ;
    if (argv[i][0] != '-') 
      // Add second of pair
      if (++i < argc) echo += std::string(argv[i]);
    echo += "\n";
  }
  // Exit here if not reading anything other than command line
  if (! ReadInput) return;

  std::string new_line;
  while (!input_stream.eof()) 
  {
    std::getline(input_stream, new_line,'\n');
    echo += new_line + '\n';
  //  std::istringstream line_stream(const_cast<char*>(new_line.c_str()),new_line.size());
    std::istringstream line_stream(const_cast<char*>(new_line.c_str()));
    std::string keyword("");
    line_stream >> keyword; //remove leading white spaces
    std::vector<std::string> end_keys = getEndKeys();
    bool exit_input_stream(false);
    for (size_t i = 0; i < end_keys.size(); i++)
      if (!stoup(keyword).find(end_keys[i])) 
        exit_input_stream = true;
    if (exit_input_stream) break;
    if (!stoup(keyword).find("@"))  
    {
      std::string filename("");
      if (keyword.size() == 1) //just the @ character
        line_stream >> filename; //take next white space separated string
      else
        filename = keyword.substr(1,keyword.size()-1); //take off the leading @
      std::ifstream file_stream(const_cast<char*>(filename.c_str()));
      if (file_stream)
      {
        char ch;
        while (file_stream.get(ch)) echo += ch;
        file_stream.close();
      }
      else
      {
        //do not throw an error immediately 
        //capture the rest of the input and then decide what to do with the unopened files later
        //user can either catch exception and forget it, or choose to exit
        files_not_found.push_back(filename);
      }
    }
  }
  if (files_not_found.size())
  throw PreprocessorError(echo,files_not_found);
}
  
Preprocessor::~Preprocessor() { }

std::string Preprocessor::Echo() 
{ return echo; }

void  Preprocessor::addLine(std::string line) 
{
  std::vector<std::string> end_keys = getEndKeys();
  for (size_t i = 0; i < end_keys.size(); i++)
    deleteLine(end_keys[i]);
  echo += line + "\n"; 
}

void  Preprocessor::deleteLine(std::string key)
{ 
 // std::istringstream echo_stream(const_cast<char*>(echo.c_str()),echo.size());
 // echo = "";
  size_t pos(0);
  while (pos < echo.size()) 
  {
    size_t pos_return = echo.find('\n',pos);
    if (pos_return != std::string::npos);
    {
      std::string new_line = echo.substr(pos,pos_return-pos);
      int key_pos(0); //first position that is interesting (not a space)
      while (key_pos < int(new_line.size()) && (std::isspace)(new_line[key_pos])) key_pos++;
      if (int(new_line.find(key)) == key_pos) echo.erase(pos,pos_return-pos+1);
      else pos = pos_return+1;
    }
  }
}

void Preprocessor::addFile(std::string filename) 
{
  std::vector<std::string> end_keys = getEndKeys();
  for (size_t i = 0; i < end_keys.size(); i++)
    deleteLine(end_keys[i]);
  std::ifstream file_stream(const_cast<char*>(filename.c_str()));
  if (file_stream)
  {
    char ch;
    while (file_stream.get(ch)) echo += ch;
    file_stream.close();
  }
  else
  {
    std::vector<std::string> files_not_found(1,filename);
    throw PreprocessorError(echo,files_not_found);
  }
}

std::string Preprocessor::stoup(const std::string & str)
{
  char ch;
  std::stringstream s;
  std::string upper_str="";
  s << str;
  while (s.get(ch))
    upper_str += std::toupper(ch); //convert to uppercase
  return upper_str;
}

std::vector<std::string> Preprocessor::getEndKeys() 
{
  std::vector<std::string> end_keys;
  end_keys.push_back("END");
  end_keys.push_back("QUIT");
  end_keys.push_back("STOP");
  end_keys.push_back("KILL");
  end_keys.push_back("EXIT");
  end_keys.push_back("GO");
  /// No!  end_keys.push_back("RUN");
  end_keys.push_back("START");
  return end_keys;
}

} // phaser_io

