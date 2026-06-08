// Written by Airlie McCoy for Phaser
//(c) 2000-2005 Cambridge University Technical Services Ltd
//All rights reserved
#ifndef __Phaser_Io__Error__Classes__
#define __Phaser_Io__Error__Classes__
//@#include <phaser/Phaser.h>
#include <string>
#include <vector>
#include <exception>

namespace phaser_io {


#define PROGRAM_INTERNAL_ERROR() ::phaser_io::error(__FILE__, __LINE__)
#define PROGRAM_NOT_IMPLEMENTED() ::phaser_io::error(__FILE__, __LINE__, \
              "Not implemented.")
#define PROGRAM_ASSERT(bool) \
  if (!(bool)) throw ::phaser_io::error(__FILE__, __LINE__,\
    "Consistency check (" # bool ") failure.")

  class error : public std::exception
  {
    public:
      explicit
      error(std::string const& msg) throw();

      error(const char* file, long line, std::string const& msg = "",
            bool internal = true) throw();

      virtual ~error() throw();

      virtual const char* what() const throw();

    protected:
      std::string msg_;
  };

class InputError : public std::exception
{
  private:
    std::string message;
    std::string echo;
  public:
    virtual const char* what() const throw() { return message.c_str(); }
    virtual ~InputError() throw() {}
    InputError(std::string e,std::string m);
    std::string Message();
    std::string Echo();
};

class SyntaxError : public std::exception
{
  private:
    std::string message;
    std::string echo;
  public:
    virtual const char* what() const throw() { return message.c_str(); }
    virtual ~SyntaxError() throw() {}
    std::string Message();
    SyntaxError(std::string,std::string);
    std::string Echo();
};

class PreprocessorError : public std::exception
{
  private:
    std::string echo;
    std::vector<std::string> files_not_found;
    std::string message;
  public:
    virtual const char* what() const throw() { return message.c_str(); }
    virtual ~PreprocessorError() throw() {}
    std::string Message();
    PreprocessorError(std::string,std::vector<std::string>);
    std::string partialEcho();
};

class WriteError : public std::exception
{
  private:
    std::string filename;
    std::string message;
  public:
    virtual const char* what() const throw() { return message.c_str(); }
    virtual ~WriteError() throw() {}
    std::string Message();
    std::string Filename();
    WriteError(std::string,std::string);
};

} // phaser_io

#endif
