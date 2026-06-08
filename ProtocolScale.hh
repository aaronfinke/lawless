// ProtocolScale.hh


#ifndef PROTOCOLSCALE_HEADER
#define PROTOCOLSCALE_HEADER

#include "phaser_types.hh"
#include "ProtocolBase.h"
#include "jiffy.h"

namespace phaser {


class ProtocolScale : public ProtocolBase
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
    ProtocolScale() 
    { init(); }

    ProtocolScale(const int& ncyc, const std::string& minimizer="BFGS")
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
