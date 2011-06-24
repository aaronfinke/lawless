// interpretcommandline.hh
//

#ifndef INTERPRETCOMMANDLINE_HEADER
#define INTERPRETCOMMANDLINE_HEADER

#include "CCP4base.hh"

namespace phaser_io {

  class InterpretCommandLine : public InputBase, virtual public CCP4base
    // Read file names from command line
    //  HKLIN, XDSIN, HKLREF, XYZIN, HKLOUT, XMLOUT 
    // Although these may also be read from commands, there is
    // separate code here since the syntax may vary slightly,
    // specifically [HKLIN] may be omitted if it is the only file 
  {
  public:
    InterpretCommandLine(int argc, char* argv[]);
    InterpretCommandLine(Preprocessor& CommandLine);
    Token_value parse(std::istringstream&) {return END;}
    void analyse(void) {}

    std::string getHKLIN1();   // return 1st one or ""
    std::vector<std::string> getHKLIN();
    std::string getXDSIN();
    std::string getSCAIN();
    std::string getHKLREF();
    std::string getHKLOUT();
    std::string getXMLOUT();
    std::string getXYZIN();
    //  copyFlag true to just copy file
    bool CopyFlag() const {return copy;}
  private:
    std::vector<std::string> HklinNames;
    std::string XDSinName;
    std::string SCAinName;
    std::string HklrefName;
    std::string HkloutName;
    std::string XmloutName;
    std::string XyzinName;

    bool copy;  // true from option "-c[opy]", just copy file

    void initialise(Preprocessor& CommandLine);
  };

} // phaser_io
#endif
