// interpretcommandline.hh
//

#ifndef INTERPRETCOMMANDLINE_HEADER
#define INTERPRETCOMMANDLINE_HEADER

#include "CCP4base.hh"
#include "Output.hh"

namespace phaser_io {

  class InterpretCommandLine : public InputBase, virtual public CCP4base
    // Read file names from command line
    //  HKLIN, XDSIN, HKLREF, XYZIN, HKLOUT, XMLOUT 
    // Although these may also be read from commands, there is
    // separate code here since the syntax may vary slightly,
    // specifically [HKLIN] may be omitted if it is the only file 
  {
  public:
    InterpretCommandLine(int argc, char* argv[], phaser_io::Output& output);
    InterpretCommandLine(Preprocessor& CommandLine, phaser_io::Output& output);
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
    std::string getHKLOUTUNMERGED();
    std::string getSCAOUT();
    std::string getSCAOUTUNMERGED();

    //  copyFlag true to just copy file
    bool CopyFlag() const {return copy;}
  private:
    std::vector<std::string> HklinNames;
    std::string XDSinName;
    std::string SCAinName;
    std::string HklrefName;
    std::string HkloutName;
    std::string HkloutUnmergedName;
    std::string ScaoutName;
    std::string ScaoutUnmergedName;
    std::string XmloutName;
    std::string XyzinName;

    bool copy;  // true from option "-c[opy]", just copy file

    void initialise(Preprocessor& CommandLine, phaser_io::Output& output);

    // return true if field is one of the recognised "logical" file names
    bool otherFiles(const std::string& field) const;

  };

} // phaser_io
#endif
