//  interpretcommandline.cpp

#include "interpretcommandline.hh"
#include "scala_util.hh"
#include "jiffy.hh"
#include "string_util.hh"
#include "version.hh"

using phaser_io::LOGFILE;
using phaser_io::stoup;

namespace phaser_io {
  //--------------------------------------------------------------
  InterpretCommandLine::InterpretCommandLine(int argc, char* argv[],
       phaser_io::Output& output)
    : HklrefName(""), HkloutName(""), XmloutName(""), XyzinName("")
  {
    Preprocessor CommandLine(argc,argv,false);
    initialise(CommandLine, output);
  }
  //--------------------------------------------------------------
InterpretCommandLine::InterpretCommandLine(Preprocessor& CommandLine,
     phaser_io::Output& output)
    : XDSinName(""), HklrefName(""), HkloutName(""), XmloutName(""), XyzinName("")
  {
    initialise(CommandLine, output);
  }
  //--------------------------------------------------------------
  void InterpretCommandLine::initialise(Preprocessor& CommandLine,
     phaser_io::Output& output)
  {
    HklinNames.clear();
    parseCCP4(CommandLine);
    // command line has been divide into lines each containing
    // either switch (string beginning "-" or pairs of tokens
    // BUT we cannot assume things come in pairs, if there were wild-cards
    // so reparse

    // Split line either at " " or "\n"
    std::vector<std::string> fields =
      StringUtil::split(CommandLine.Echo(), " ", "\n");

    copy = false;
    noinput = false;
    run = true;

    int ifld = 0;

    commandlineArguments = "";  // for echoing command line arguments

    while (ifld < int(fields.size())) {
      if (fields[ifld][0] == '-')       {
        // Switch, ie string beginning with '-'
        //                  std::cout << "Command line switch found: "
        //                            << string_value << "\n";
        std::string option = fields[ifld++];
        if (option == "--help") {
          std::string s = "Aimless " + PROGRAM_VERSION + "\n";
          s += "Usage: aimless [options] hklin <filein> hklout <fileout> etc\n";
          s += " File assignments may be done on the command line or as input commands\n";
          s += "Options:\n";
          s += " --no-input   run immediately without waiting for input\n";
          s += "\nSee aimless.html for full documentation\n";
          output.logTab(0,LOGFILE, s);
          run = false;
          return;
        }
        if ((option == "--version") || (option == "-v")) {
          std::string s = "Aimless " + PROGRAM_VERSION + "\n";
          output.logTab(0,LOGFILE, s);
          run = false;
          return;
        }
        if (option.substr(0,2) == "-c") {
          copy = true;
          commandlineArguments += "-copy\n";
        }
        if (option == "--no-input") {
          noinput = true;
          commandlineArguments += "--no-input\n";
        }
      } else  {
        bool fieldpair = true;
        if (stoup(fields[ifld]) == "HKLIN") {
          HklinNames.push_back(fields[++ifld]);
        }
        else if (stoup(fields[ifld]) == "XDSIN")  {
          XDSinName = fields[++ifld];
        }
        else if (stoup(fields[ifld]) == "SCAIN")  {
          SCAinName = fields[++ifld];
        }
        else if (stoup(fields[ifld]) == "HKLREF") {
          HklrefName = fields[++ifld];
        }
        else if (stoup(fields[ifld]) == "HKLOUT") {
          HkloutName = fields[++ifld];
        }
        else if ((stoup(fields[ifld]) == "UNMERGEDOUT")|
                 (stoup(fields[ifld]) == "HKLOUTUNMERGED")) {
          HkloutUnmergedName = fields[++ifld];
        }
        else if (stoup(fields[ifld]) == "SCALEPACK") {
          ScaoutName = fields[++ifld];
        }
        else if (stoup(fields[ifld]) == "SCALEPACKUNMERGED") {
          ScaoutUnmergedName = fields[++ifld];
        }
        else if (stoup(fields[ifld]) == "XMLOUT") {
          XmloutName = fields[++ifld];
        }
        else if (stoup(fields[ifld]) == "XYZIN") {
          XyzinName = fields[++ifld];
        }
        else if (otherFiles(stoup(fields[ifld]))) {
          // Other file name, ignore
          ifld++;
        }
        else {
          HklinNames.push_back(fields[ifld]);
          fieldpair = false;
        }
        if (fieldpair && ifld > 0) {
          commandlineArguments +=  fields[ifld-1] + " ";
        }
        commandlineArguments += fields[ifld] + "\n";
        ifld++;
      }
    }
    // Add .mtz if needed, ie non-blank and no extension already
    for (size_t i=0;i<HklinNames.size();i++) {
      scala::AddFileExtension(HklinNames[i],"mtz");
    }
    scala::AddFileExtension(HklrefName,"mtz");
    scala::AddFileExtension(HkloutName,"mtz");
    scala::AddFileExtension(XmloutName,"xml");
    scala::AddFileExtension(XyzinName,"pdb");
  }
  //--------------------------------------------------------------
  void InterpretCommandLine::printCommandLine(phaser_io::Output& output) const
  {
    if (commandlineArguments != "") {
      output.logTab(0,LOGFILE,
                    "==== Command line arguments ====\n" +
                    commandlineArguments + "\n");
    }
  }
  //--------------------------------------------------------------
  std::string InterpretCommandLine::getHKLIN1()
  {if (HklinNames.size() > 0) {
      return HklinNames[0];
    } else {return "";}
  }
  //--------------------------------------------------------------
  std::vector<std::string> InterpretCommandLine::getHKLIN()
  {return HklinNames;}
  //--------------------------------------------------------------
  std::string InterpretCommandLine::getXDSIN()
  {return XDSinName;}
  //--------------------------------------------------------------
  std::string InterpretCommandLine::getSCAIN()
  {return SCAinName;}
  //--------------------------------------------------------------
  std::string InterpretCommandLine::getHKLREF()
  {return HklrefName;}
  //--------------------------------------------------------------
  std::string InterpretCommandLine::getHKLOUT()
  {return HkloutName;}
  //--------------------------------------------------------------
  std::string InterpretCommandLine::getXMLOUT()
  {return XmloutName;}
  //--------------------------------------------------------------
  std::string InterpretCommandLine::getXYZIN()
  {return XyzinName;}
  //--------------------------------------------------------------
  std::string InterpretCommandLine::getHKLOUTUNMERGED()
  {return HkloutUnmergedName;}
  //--------------------------------------------------------------
  std::string InterpretCommandLine::getSCAOUT()
  {return ScaoutName;}
  //--------------------------------------------------------------
  std::string InterpretCommandLine::getSCAOUTUNMERGED()
  {return ScaoutUnmergedName;}
  //--------------------------------------------------------------
  bool InterpretCommandLine::otherFiles
  (const std::string& field) const
  // return true if field is one of the recognised "logical" file names
  {
    std::string files[] = {
      "NORMPLOT", "ANOMPLOT", "ROGUES", "ROGUEPLOT", "CORRELPLOT",
      "SCALES", "TILEIMAGE"};
    const int N = 7;  // number of filenames
    std::vector<std::string> names(files, files+N);
    for (size_t k=0; k<names.size(); k++) {
      if (field == names[k]) {
        return true;
      }
    }
    return false;
  }
  //--------------------------------------------------------------
}  // phaser_io
