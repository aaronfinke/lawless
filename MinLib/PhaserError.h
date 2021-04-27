//(c) 2000-2015 Cambridge University Technical Services Ltd
//All rights reserved
#ifndef __PhaserErrorClass__
#define __PhaserErrorClass__
#include <exception>
#include <string>
#include <cstdlib>

namespace phaser {

// these error codes are offset with 64 as in
// /usr/include/sysexists.h
enum errorType {
 NULL_ERROR, //0
 SYNTAX,     //64
 INPUT,      //65
 FILEOPEN,   //66
 MEMORY,     //67
 KILLFILE,   //68
 KILLTIME,   //69
 FATAL,      //70
 UNHANDLED,  //71
 UNKNOWN     //72
};

class PhaserError : public std::exception
{
  protected:
    errorType   type;
    std::string message;
    std::string echo;
    static const int EXIT_BASE = 63;
  public:
    PhaserError(): type(NULL_ERROR),message(""),echo("") { }
    PhaserError(errorType t): type(t),message(""),echo("") { }
    PhaserError(errorType t,std::string m): type(t),message(m),echo("") { }
    PhaserError(errorType t,const char* m): type(t),message(m),echo("") { }
    PhaserError(errorType t,std::string e,std::string m): type(t),message(m),echo(e) { }

    void setPhaserError(PhaserError err) { type=err.ErrorType(); message=err.ErrorMessage(); echo=err.Cards(); }
    void setPhaserError(errorType t,std::string m) { type=t; message=m; }
    errorType   ErrorType() const  { return type; }
    std::string ErrorMessage() const { return message; }
    std::string ErrorName() const;
    bool        Success() const { return (type == NULL_ERROR); }
    bool        Failure() const { return (type != NULL_ERROR); }
    int         ExitCode() const { return Success() ? EXIT_SUCCESS: (type + EXIT_BASE); }
    std::string Cards() const   { return echo; }
//virtual
// Commented out to avoid compiler warning
//    virtual const char* what() const throw() { return ErrorName().c_str(); } //more useful than message
    virtual ~PhaserError() throw() {}
};

} //phaser

#endif
