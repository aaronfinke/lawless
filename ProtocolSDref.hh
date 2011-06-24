// ProtocolSdref.hh


#ifndef PROTOCOLSDREF_HEADER
#define PROTOCOLSDREF_HEADER

#include "phaser_types.hh"
#include "ProtocolBase.h"
#include "jiffy.h"

namespace phaser {


class ProtocolSdref : public ProtocolBase
{
  private:
    unsigned    NCYC;
    std::string MINIMIZER;

    void init()
    {
      NCYC = 50;
      MINIMIZER = "BFGS";
    }

  public: 
    ProtocolSdref() 
    { init(); }

    ProtocolSdref(const int& ncyc, const std::string& minimizer="BFGS")
    {
      init();
      if (ncyc >= 0) NCYC = ncyc;
      if (minimizer == "BFGS" || minimizer == "NEWTON" || minimizer == "DESCENT")
        MINIMIZER = minimizer;
    }

  //-------------------------
  //concrete member functions
  //-------------------------

    std::string getMINIMIZER() const
    { return MINIMIZER; }

    unsigned getNCYC() const
    { return NCYC; }
   

    std::string logfile() const
    {
      return "";
      ///      std::string tab(3,' ');
	///      return tab + "Refinement protocol for this macrocycle:  FREE\n" ;
    }
};

} //phaser

#endif
