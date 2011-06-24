// openinputfile.hh

//  setup_columns

#ifndef OPENINPUTFILE_HEADER
#define OPENINPUTFILE_HEADER

#include "mtz_unmerge_io.hh"

namespace MtzIO
{
  //--------------------------------------------------------------
  const int MAXNCOLUMNS = 21;
  // Create list of required columns
  column_labels setup_columns();
} // MtzIO

#endif
