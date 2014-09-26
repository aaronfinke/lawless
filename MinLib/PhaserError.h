//(c) 2000-2007 Cambridge University Technical Services Ltd
//All rights reserved
#ifndef __PhaserErrorClass__
#define __PhaserErrorClass__
#include <exception>
#include <string>

namespace phaser {

enum errorType {
 NO_ERROR,
 SYNTAX,INPUT,
 FILEOPEN,FATAL,MEMORY,UNHANDLED,UNKNOWN
};

class PhaserError : public std::exception
{
  protected:
    errorType   type;
    std::string message;
  public:
    PhaserError(): type(NO_ERROR),message("") { }
    PhaserError(errorType t): type(t),message("") { }
    PhaserError(errorType t,std::string m): type(t),message(m) { }
    PhaserError(errorType t,const char* m): type(t),message(m) { }
    void setPhaserError(PhaserError err) { type=err.ErrorType(); message=err.ErrorMessage(); }
    void setPhaserError(errorType t,std::string m) { type=t; message=m; }
    errorType   ErrorType() const  { return type; }
    std::string ErrorMessage() const { return message; }
    std::string ErrorName() const;
    std::string XML() const;
    bool        Success() const { return (type == NO_ERROR); }
    bool        Failure() const { return (type != NO_ERROR); }
//virtual
    virtual const char* what() const throw() { return ErrorName().c_str(); } //more useful than message
    virtual ~PhaserError() throw() {}
};

} //phaser

#endif
