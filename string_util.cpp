// string_util.cpp

#include <stdarg.h>
#include <cstdio>

#include "string_util.hh"
#define ASSERT assert
#include <assert.h>
#include <clipper/clipper.h>


//--------------------------------------------------------------
std::string StringUtil::Strip(const std::string& s, const char& x)
// Strip out space and x characters from string
// x defaults to null
{
  std::string ss;
  for (size_t i=0;i<s.size();i++) {
    if (s[i] != ' ' && s[i] != x) ss.push_back(s[i]);
  }
  return ss;
}
//--------------------------------------------------------------
std::string StringUtil::StripNull(const std::string& s)
// Strip out null characters from string
{
  std::string ss;
  for (size_t i=0;i<s.size();i++)
    if (s[i] != '\0') ss.push_back(s[i]);
  return ss;
}
//--------------------------------------------------------------
std::string StringUtil::Trim(const std::string& s)
// Trim off leading & trailing spaces from string
{
  std::string ss;
  if (s.size() > 0) {
      bool instring = false;  // first non-space character
      for (size_t i=0;i<s.size();i++) {
	if (instring || s[i] != ' ') {
	  ss.push_back(s[i]);
	  instring = true;
	}
      }
      // now remove trailing spaces
      int j = ss.size()-1;
      while (j >=0) 
	{if (ss[j--] != ' ') break;}
      ss.resize(j+2);
    }
  return ss;
}
//--------------------------------------------------------------
// Return string of length <fieldwidth> with text centred on position
// cenpos (numbered from 0)
std::string StringUtil::CentreString(const std::string& text, const int& fieldwidth,
			 const int& cenposition)
{
  int cenpos = cenposition;
  if (cenpos == 0) cenpos = fieldwidth/2;
  //  ASSERT (cenpos < fieldwidth);
  int len = text.size();
  if (len >= fieldwidth) return text;
  // number of trailing spaces
  int p2 = clipper::Util::max(0,clipper::Util::min
			      (fieldwidth-1-(cenpos+len/2), fieldwidth-1)); 
  int p1 = fieldwidth - len - p2; // number of leading spaces
  std::string pad1;
  if (p1 > 0) pad1 = std::string(p1, ' ');
  std::string pad2;
  if (p2 > 0) {
    pad2 = std::string(p2, ' ');
  }
  //    std::cout << "CS " << fieldwidth << " "
  //	      << pad1.size() << " " << text.size()<< " " << pad2.size() << "\n"; //
  return pad1+text+pad2;
}
//--------------------------------------------------------------
// Return string of length Max(fieldwidth, Length(text))
std::string StringUtil::PadString(const std::string& text, const int& fieldwidth)
{
  int len = text.size();
  if (len >= fieldwidth) return text;
  std::string pad1((fieldwidth-len), ' ');
  return text+pad1;
}
//--------------------------------------------------------------
// Return right-justified string
std::string StringUtil::RightString(const std::string& text, const int& fieldwidth)
{
  int len = text.size();
  if (len >= fieldwidth) return text;
  std::string pad1((fieldwidth-len), ' ');
  return pad1+text;
}
//--------------------------------------------------------------
// Return left-justified string
std::string StringUtil::LeftString(const std::string& text, const int& fieldwidth)
{
  int len = text.size();
  if (len >= fieldwidth) return text;
  std::string pad1((fieldwidth-len), ' ');
  return text+pad1;
}
//--------------------------------------------------------------
std::vector<std::string> StringUtil::split(const std::string& str,
			       const std::string& sep)
// Return substrings split at string "sep" (excluded)
// modified from clipper_types.cpp
{
  std::vector<std::string> splitstr;
  size_t tokbeg = 0, tokend = 0;
  while (1) {
    tokbeg = str.find_first_not_of(sep, tokend);
    if (tokbeg == std::string::npos) return splitstr;
    tokend = str.find_first_of(sep, tokbeg);
    splitstr.push_back(str.substr(tokbeg, tokend-tokbeg) );
    if (tokend == std::string::npos) return splitstr;
  }
}
//--------------------------------------------------------------
std::vector<std::string> StringUtil::split(const std::string& str,
			       const std::string& sep1, const std::string& sep2)
// Return substrings split at string "sep1" or "sep2" (excluded)
// modified from clipper_types.cpp
{
  std::vector<std::string> splitstr;
  size_t tokbeg = 0, tokend = 0;
  while (1) {
    tokbeg = clipper::Util::max(str.find_first_not_of(sep1, tokend),
		 str.find_first_not_of(sep2, tokend));
    if (tokbeg == std::string::npos) return splitstr;
    tokend = clipper::Util::min(str.find_first_of(sep1, tokbeg),
		 str.find_first_of(sep2, tokbeg));
    if (tokend-tokbeg > 0) {
      splitstr.push_back(str.substr(tokbeg, tokend-tokbeg) );
    }
    if (tokend == std::string::npos) return splitstr;
  }
}
//--------------------------------------------------------------
// <tag><data</tag>
std::string StringUtil::MakeXMLtag(const std::string& tag, const std::string& data)
{
  return "<"+tag+">"+data+"</"+tag+">";
}
//--------------------------------------------------------------
std::string StringUtil::itos(const int f, const int w)
{ std::ostringstream s; s.width(w); s.setf(std::ios::fixed);s << f; return s.str(); }
//--------------------------------------------------------------
std::string StringUtil::ftos(const float f, const int w, const int d)
{ std::ostringstream s; s.width(w); s.setf(std::ios::fixed); s.precision(d);s << f; return s.str(); }
//--------------------------------------------------------------
std::string StringUtil::etos(const float f, const int w, const int d)
{ std::ostringstream s; s.width(w); s.precision(d);s << f; return s.str(); }
//--------------------------------------------------------------
std::string StringUtil::ftos(const double f, const int w, const int d)
{ std::ostringstream s; s.width(w); s.setf(std::ios::fixed); s.precision(d);s << f; return s.str(); }
//--------------------------------------------------------------
std::string StringUtil::etos(const double f, const int w, const int d)
{ std::ostringstream s; s.width(w); s.precision(d);s << f; return s.str(); }
//--------------------------------------------------------------
std::string StringUtil::WrapLine(const std::string& line,
	 const int& pagewidth, const int& nindent, const std::string& sepc)
// wrap after field terminated by character sepc (default " ")
// if sepc = " ", exclude it from field
{
  if (int(line.size()) <= pagewidth) return line; // nothing to do
  std::string str;
  std::string sep = sepc;
  if (sep == "") sep = " ";
  // make list all end of field positions
  std::vector<int> eon;
  size_t tokbeg = 0, tokend = 0;
  tokbeg = line.find_first_not_of(sep, tokend); // start of first field
  if (tokbeg == std::string::npos) return line; // return everything
  while (true) {
    // look for end of field
    tokend = line.find_first_of(sep, tokbeg);
    if (tokend == std::string::npos) break;   // end of line
    if (sep != " ") tokend++;  // include non-space separator
    eon.push_back(tokend);
    //^    std::cout << eon.back() << "eon\n"; //^
    // look for start of next field
    tokbeg = line.find_first_not_of(sep, tokend);
    if (tokbeg == std::string::npos) break; // end of line
  }
  // we have now in eon a list of character positions one beyond each field
  // except for the last
  int i1 = 0;  // position of start of current line
  int pgw = pagewidth;// current page width
  // Indent lines after first by nindent fields, and try to line up next field
  int indent=0;
  int nextfldw = 0;
  if (nindent > 0 && nindent < int(eon.size())) {
    indent = eon[nindent-1];
    nextfldw = eon[nindent] - eon[nindent-1]; // width of next field
  }

  for (size_t ifd=0;ifd<eon.size()-2;++ifd) { // allow two trailing fields
    if (eon[ifd] > pgw+i1) {
      int indt = indent + nextfldw - (eon[ifd+1]-eon[ifd]);
      if (indt < 0) indt = 0;
      //^      std::cout <<"wrap "<<ifd<<" "<<eon[ifd]<<" "<<i1<<" "<<pgw<<" "<<indt<<"\n";
      std::string indentstring(indt,' ');
      str += line.substr(i1, eon[ifd]-i1)+"\n"+indentstring;
      i1 =  eon[ifd];
      if (pgw == pagewidth) pgw -= indent; // indent after first line
    }
  }
  // last bit
  if (i1 < int(line.size())) {
    str += line.substr(i1, line.size()-i1);
  }
  return str;
}
//--------------------------------------------------------------
std::string StringUtil::formatFraction(const double& fr, const int& width)
//!< Format a real number as a fraction ie "n/m"
/*! \param  fr  number of format
    \param  width  field width for decimal version, if fail to find suitable fraction
*/
{
  std::string s;
  const int MAXDEN = 24; // maximum denominator to try
  double f = fr;
  if (fr < 0.0) { // negative, add sign
    f = -fr;
    s += "-";
  }
  double tol = 0.5/double(MAXDEN);  // tolerance for nearest integer
  int den = -1;
  for (int i=1;i<=MAXDEN;++i) { // try integer divisors up to MAXDEN
    if ((f*double(i) - floor(f*double(i)+tol)) < tol) {
      // found a suitable denominator i
      den = i;
      break;
    }
  }
  if (den < 0) {
    // not found, just format as decimal number
    int w = width;
    int d = w-3;  // number of decimal digits
    if (d<3) {d = 3;w = d+3;}
    if (f > 1.0) {
      int l = int(log10(f+0.0001));
      d = w - (l+3);
      if (d<0) {d=0; w=l+3;}
    }
    s = clipper::String(f, w, d);
  } else {
    // format fraction
    int num = clipper::Util::intr(f*double(den)); // numerator
    s += clipper::String(num);
    if (den > 1) s+= "/" + clipper::String(den);
  }
  return StringUtil::Strip(s);
}
//--------------------------------------------------------------
std::string StringUtil::ToUpper(const std::string& s)
//!< returns uppercase version of string
{
  std::string ss = s;
  for (size_t i=0;i<s.size();++i) {
    ss[i] = std::toupper(ss[i]);
  }
  return ss;
}
//--------------------------------------------------------------
std::string StringUtil::BuftoLine(const std::string buf)
// Extract line from buffer, removing any trailing Cr or Lf characters
{
  std::string line(buf);
  size_t ll = line.size();
  while (ll > 0) {
    if (line[ll-1] != '\n' && line[ll-1] != '\r') {break;}
    ll--;
  }
  return line.substr(0,ll);
}
//--------------------------------------------------------------
std::string StringUtil::FormatSaveVector(const std::vector<int> ivec)
// format integer vector for dump/save
{
  std::string s = "";
  std::string line = "";
  const unsigned int MAXLINE = 100;
  for (size_t i=0;i<ivec.size();++i) {
    line += " "+clipper::String(ivec[i]);
    if (line.size() > MAXLINE) {
      s += line+"\n";
      line = "";
    }
  }
  return s+line+"\n";
}
//--------------------------------------------------------------
std::string StringUtil::FormatSaveVector(const std::vector<double> vec)
// format double vector for dump/save
{
  std::string s = "";
  std::string line = "";
  const unsigned int MAXLINE = 100;
  for (size_t i=0;i<vec.size();++i) {
    line += " "+clipper::String(vec[i]);
    if (line.size() > MAXLINE) {
      s += line+"\n";
      line = "";
    }
  }
  return s+line+"\n";
}
//======================================================================
//! just add leading tabs to string and newline if not there already
std::string FormatOutput::logTab(const int& tab, const std::string& text)
{
  std::string nl;
  if (text[text.size()-1] != '\n') nl = "\n";
  return std::string(3*tab, ' ') + text + nl;
}
//--------------------------------------------------------------
//! format using vsstringf
std::string FormatOutput::logTabPrintf(const int& tab,
				  const char* formattext,...)
{
  static const std::size_t temp_size = 8192;
  char temp[temp_size];
  temp[temp_size-1] = '\0';
  va_list arglist;
  va_start(arglist,formattext);
  vsprintf(temp,formattext,arglist);
  va_end(arglist);
  assert(temp[temp_size-1] == '\0');
  return std::string(3*tab, ' ') + std::string(temp);
}
//--------------------------------------------------------------
//! just add newline if not there already
std::string FormatOutput::logWarning(const std::string& text)
{
  std::string nl;
  if (text[text.size()-1] != '\n') nl = "\n";
  return "\nWARNING! " + text + nl;
}

//======================================================================
//! constructor from type, maximum value, minimum field width, and precision
Numberfield::Numberfield(const bool& IntType, const float& MaxValue,
			 const int& MinWidth, const int& Precision)
  : width(-1)
{
  init(IntType, MaxValue, MinWidth, Precision);
}
//--------------------------------------------------------------
//! (re)initialise from type, minimum field width, and precision
  /*! \param IntType    true if integer type
      \param MaxValue   maximum ||value|| for field
      \param MinWidth   minimum field width
      \param Precision  number of significant figures
   */
void Numberfield::init(const bool& IntType, const float& MaxValue,
		       const int& MinWidth, const int& Precision)
{
  type = (IntType) ? +1 : 0; ;
  if (width < 0) {  // unknown width
    // Number of digits for integer part + 1 for luck
    int l = int(log10(MaxValue)+1.001);
    int prec = clipper::Util::max(Precision, l+1);
    dec =  clipper::Util::max(0, prec-l-1);  // number after decimal point
    width = clipper::Util::max(MinWidth, l+dec+3);
    //^    std::cout <<"Numberfield "<<MinWidth <<" "<<width <<" "
    //^	      <<Precision<<" "<<prec
    //^	      <<" "<<l<<" "<<dec <<" "<<label1<<" "<<label2<<"\n";
  } else {
    // width from previous construction, reset to minimum if necessary
    if (width < MinWidth) {
      width = MinWidth;
      dec =  clipper::Util::max(0, width-Precision-1);  // number after decimal point
    }
  }
}
//--------------------------------------------------------------
//! (re)initialise from type, field width, and number of characters after decimal point
void Numberfield::init(const int& Type, const int& Width, const int& Dec,
	    const std::string& Label1, const std::string& Label2)
{type = Type; width = Width; dec = Dec; label1 = Label1; label2 = Label2;}
//--------------------------------------------------------------

//--------------------------------------------------------------



