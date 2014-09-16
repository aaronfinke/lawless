//(c) 2000-2007 Cambridge University Technical Services Ltd
//All rights reserved
#ifndef __ProtocolBaseClass__
#define __ProtocolBaseClass__

#include "phaser_types.h"
//#include <phaser/include/Phaser.h>

namespace phaser {

//abstract - reminder! not objects of an abstract
//base class can be instantiated
class ProtocolBase  //abstract 
{
  public:
    ProtocolBase() {}
    virtual ~ProtocolBase() {}
  
//these functions are called either in callRefine or in setProtocol
//in the base object
    virtual bool           is_default() const;
    virtual bool           is_off() const;
    virtual void           setDEFAULT();
    virtual bool           getFIX(int) const;
    virtual std::string    getMINIMIZER() const;
    virtual unsigned       getNCYC() const;
    virtual std::string    getTXT(int) const;
    virtual floatType      getNUM(int) const;
    virtual std::string    unparse() const;
    virtual std::string    logfile() const;
};

inline bool           ProtocolBase::is_default() const { return false; }
inline bool           ProtocolBase::is_off() const { return false; }
inline void           ProtocolBase::setDEFAULT() { }
inline bool           ProtocolBase::getFIX(int) const { return false; }
inline std::string    ProtocolBase::getMINIMIZER() const { return ""; }
inline unsigned       ProtocolBase::getNCYC() const { return 0; }
inline std::string    ProtocolBase::getTXT(int i) const { return ""; }
inline floatType      ProtocolBase::getNUM(int i) const { return 0; }
inline std::string    ProtocolBase::unparse() const { return ""; }
inline std::string    ProtocolBase::logfile() const { return ""; }

typedef boost::shared_ptr<ProtocolBase> protocolPtr;

} //phaser

#endif
