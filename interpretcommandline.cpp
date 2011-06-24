//  interpretcommandline.cpp

#include "interpretcommandline.hh"
#include "scala_util.hh"
#include "jiffy.hh"
#include "string_util.hh"

using phaser_io::stoup;

namespace phaser_io {
  //--------------------------------------------------------------
  InterpretCommandLine::InterpretCommandLine(int argc, char* argv[])
    : HklrefName(""), HkloutName(""), XmloutName(""), XyzinName("")
  {
    Preprocessor CommandLine(argc,argv,false);
    initialise(CommandLine);
  }
  //--------------------------------------------------------------
  InterpretCommandLine::InterpretCommandLine(Preprocessor& CommandLine)
    : XDSinName(""), HklrefName(""), HkloutName(""), XmloutName(""), XyzinName("")
  {
    initialise(CommandLine);
  }
  //--------------------------------------------------------------
  void InterpretCommandLine::initialise(Preprocessor& CommandLine)
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
    
    int ifld = 0;

    while (ifld < int(fields.size())) {
      if (fields[ifld][0] == '-')	{
	// Switch, ie string beginning with '-'	
	//		    std::cout << "Command line switch found: " 
	//			      << string_value << "\n";
	if (fields[ifld++].substr(0,2) == "-c") {
	  copy = true;
	}
      } else  {
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
	else if (stoup(fields[ifld]) == "XMLOUT") {
	  XmloutName = fields[++ifld];
	}
	else if (stoup(fields[ifld]) == "XYZIN") {
	  XyzinName = fields[++ifld];
	}
	else {
	  HklinNames.push_back(fields[ifld]);
	}
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
  //--------------------------------------------------------------
}  // phaser_io

