// Written by Airlie McCoy for Phaser
//(c) 2000-2005 Cambridge University Technical Services Ltd
//All rights reserved
#ifndef __InputBase__Class__
#define __InputBase__Class__
#include <sstream>

namespace phaser_io {

enum Token_value {
  NAME, NUMBER, PLUS, MINUS, END, ENDLINE, ASSIGN
};

//reminder! no objects of an abstract base class can be instantiated
class InputBase  //abstract 
{
  public:
    InputBase() {}
    virtual ~InputBase() { }
    virtual Token_value parse(std::istringstream&) = 0;
    virtual void        analyse(void) = 0;
};

typedef InputBase* inputPtr;

} // phaser_io

#endif
