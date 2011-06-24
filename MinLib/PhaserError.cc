//(c) 2000-2007 Cambridge University Technical Services Ltd
//All rights reserved
#include "PhaserError.h"

namespace phaser {

std::string PhaserError::XML() const
{
  std::string xml;
  std::string status(Failure() ? "error" : "ok");
  xml += "  <status>" + status + "</status>\n";
  if (Failure())
  {
    xml += "  <error>\n";
    xml += "    <name>" +  ErrorName() + "</name>\n";
    xml += "    <message>"  +  ErrorMessage()  +  "</message>\n";
    xml += "  </error>\n";
  }
  return xml;
}

std::string PhaserError::ErrorName() const
{
  if (type == NO_ERROR) return "NO";
  else if (type == SYNTAX) return "SYNTAX";
  else if (type == INPUT) return "INPUT";
  else if (type == FILEOPEN) return "FILE OPENING";
  else if (type == FATAL) return "FATAL RUNTIME";
  else if (type == MEMORY) return "OUT OF MEMORY";
  else if (type == UNHANDLED) return "UNHANDLED";
  else if (type == UNKNOWN) return "UNKNOWN";
  return "";
}

} //phaser
