//(c) 2000-2007 Cambridge University Technical Services Ltd
//All rights reserved
#ifndef __Phaser__Error__Classes__
#define __Phaser__Error__Classes__
#include "PhaserError.h"
#include <exception>

namespace phaser {

#define PHASER_INTERNAL_ERROR() ::phaser::error(__FILE__, __LINE__)
#define PHASER_NOT_IMPLEMENTED() ::phaser::error(__FILE__, __LINE__, \
              "Not implemented.")
#define PHASER_ASSERT(bool) \
  if (!(bool)) throw ::phaser::error(__FILE__, __LINE__,\
    "Consistency check (" # bool ") failed.")

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

} //phaser

#endif
