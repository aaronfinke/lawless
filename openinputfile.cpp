// openinputfile.cpp

//  setup_columns

#include "openinputfile.hh"
#include "columnlabels.hh"

namespace MtzIO
{
//--------------------------------------------------------------

column_labels setup_columns()
// Create list of required columns
{
  column_labels column_list;
  column_list.add("H",OF_COMPULSORY);
  column_list.add("K",OF_COMPULSORY);
  column_list.add("L",OF_COMPULSORY);
  column_list.add("M_ISYM",OF_COMPULSORY); // cmtzlib translates "/" to "_"
  column_list.add("BATCH",OF_COMPULSORY);
  column_list.add("I",OF_COMPULSORY);
  column_list.add("SIGI",OF_COMPULSORY);
  column_list.add("IPR",OF_OPTIONAL);
  column_list.add("SIGIPR",OF_OPTIONAL);
  column_list.add("SCALE",OF_OPTIONAL);
  column_list.add("SIGSCALE",OF_OPTIONAL);
  column_list.add("TIME",OF_OPTIONAL);
  column_list.add("XDET",OF_OPTIONAL);
  column_list.add("YDET",OF_OPTIONAL);
  column_list.add("ROT",OF_OPTIONAL);
  column_list.add("WIDTH",OF_OPTIONAL);
  column_list.add("MPART",OF_OPTIONAL);
  column_list.add("FRACTIONCALC",OF_OPTIONAL);
  column_list.add("LP",OF_OPTIONAL);
  column_list.add("FLAG",OF_OPTIONAL);
  column_list.add("BGPKRATIOS",OF_OPTIONAL);

  return column_list;
}
} // MtzIO
