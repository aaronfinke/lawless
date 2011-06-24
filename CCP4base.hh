// CCP4base.hh
//
// Written by Airlie McCoy for Phaser
//(c) 2000-2005 Cambridge University Technical Services Ltd
//All rights reserved
#ifndef __CCP4base__Class__
#define __CCP4base__Class__
#include <vector>
#include <set>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include "Errors.hh"
#include "Preprocessor.hh"
#include "InputBase.hh"

// Clipper
#include <clipper/clipper.h>
typedef clipper::Vec3<float> FVect3;
typedef clipper::Vec3<double> DVect3;
typedef clipper::Mat33<double> DMat33;

typedef double floatType;

namespace phaser_io {

  class CCP4base 
  {
  protected: //members
    std::string  keywords;
    floatType       number_value;
    std::string  string_value;   
    Token_value  curr_tok;
    std::vector<std::string> possible_keys;
    std::vector<inputPtr> possible_fns;

  private:
    bool annotate_keywords;

  public: //member functions
    CCP4base() 
    {
      string_value = "";
      number_value = 0;
      curr_tok = END ;
      //don't clear the arrays of pointers to member functions or else the
      //class will not build up the list as functions are added
    
      //clear the keywords so that instantiating a new input object in the scope
      //of another clears the keyword list
      keywords = "";
      annotate_keywords = true; //with the Warnings
    }
    virtual ~CCP4base() { }

    // NOTE:
    //  a "token" is the type of token NAME< NUMBER, ENDLINE, END
    //  a "key" is a keyword

    // push keyword into list of possible_keys
    void Add_Key(std::string card);
    // Skip to end of line, warn of unexpected tokens, returns END or ENDLINE
    Token_value skip_line(std::istringstream& input_stream); 
    // Get next token from input
    //  Returns: (possible values of token)
    //   NUMBER   number in number_value (& in string_value)
    //            adds to keywords string
    //   NAME     non-numeric token in string_value
    //            adds to keywords string
    //   ASSIGN   '=' found
    //   ENDLINE  end of line
    //   END      end of input
    Token_value get_token(std::istringstream& input_stream);
    // Check line for keyword (must be at start of line)
    Token_value get_key(std::istringstream& input_stream);
    // look for "ON" or "OFF"
    bool getBoolean(std::istringstream& input_stream);
    // Returns rest of line
    std::string getLine(std::istringstream& input_stream);
    // Returns next token as filename (allows ,#,!)
    std::string getFileName(std::istringstream& input_stream);
    // Checks that next token is "=", throws error if not
    void getAssign(std::istringstream& input_stream);
    // Returns next token as string
    std::string getString(std::istringstream& input_stream);
    // Returns next token as number (error if invalid)
    floatType get1num(std::istringstream& input_stream);
    // Returns next 3 tokens as vector (error if invalid)
    DVect3 get3nums(std::istringstream& input_stream);
    // Returns next 9 token as matrix (error if invalid)
    DMat33 get9nums(std::istringstream& input_stream);
    // Tests current key (1st 4 characters, or less for specials)
    bool keyIs(std::string key);
    // Check token against list
    bool tokenIs(std::istringstream& input_stream,int va_len, ...);
    bool tokenIs(int va_len, ...);
    // Check current or next key against list of possibles
    // Return true if OK
    bool compulsoryKey(int va_len, ...);
    bool compulsoryKey(std::istringstream& input_stream,int va_len, ...);
    bool optionalKey(std::istringstream& input_stream,int va_len, ...);
    bool optionalKey(int va_len, ...);

    void parseCCP4(Preprocessor capture);
    floatType isperc(floatType percent);
    std::string Keywords();
};

} // phaser_io

#endif
