// file_util.hh

#ifndef FILE_UTIL_HEADER
#define FILE_UTIL_HEADER

#include <stdio.h>
#include <string>

FILE* OpenFile(const std::string& Filename, const bool& Write,
	       const bool Binary=false);

// write text to filename
void WriteToFile(const std::string& filename,
		 const std::string& text);

#endif
