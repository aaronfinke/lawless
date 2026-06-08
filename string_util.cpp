// string_util.cpp

#include <stdarg.h>
#include <cstdio>
#include <limits>

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
std::string StringUtil::Unquote(const std::string& sin)
// Trim off spaces and leading & trailing quotes from string
{
  std::string s = Trim(sin);
  std::string ss;
  if (s.size() > 0) {
    bool instring = false;  // first non-quote character
    for (size_t i=0;i<s.size();i++) {
      if (instring || ((s[i] != '\'') && (s[i] != '"'))) {
        ss.push_back(s[i]);
        instring = true;
      }
    }
    // now remove trailing spaces
    int j = ss.size()-1;
    while (j >=0) {
      if ((ss[j] != '\'') && (ss[j] != '"')) {
        break;
      }
      j--;
    }
    ss.resize(j+1);
  }
  return ss;
}
//--------------------------------------------------------------
bool StringUtil::isquoted(const std::string& s)
// true if a quoted string
{
  if (s.size()== 0) {return false;}

  bool instring = false;
  int k = -1;
  for (size_t i=0;i<s.size();i++) {
    if ((s[i] == '\'') || (s[i] == '"')) {
      k = i;
      break;
    }
  }
  if (k < 0) {return false;}

  // search for trailing quote
  int j = s.size()-1;
  bool found = false;
  while (j > k) {
    if ((s[j] == '\'') || (s[j] == '"')) {
      found = true;
      break;
    }
    j--;
  }
  return found;
}
//--------------------------------------------------------------
int StringUtil::isanumber(const std::string& s)
// return +1 if s is a valid integer
//        -1 if s is a valid floating point number
//         0 if s is a string (not a number)
// NB not foolproof
{
  if (isquoted(s)) {return 0;} // quoted string
  std::string ss = StringUtil::Trim(s);  // strip leading and trailing spaces
  int type = +1;
  for (size_t i=0;i<ss.size();i++) {
    if (isspace(s[i])) {
      return 0;  // no embedded whitespace in a number
    }
    if (!std::isdigit(s[i])) {
      // Allow +, -, e, E, .
      if (s[i] == '.' || s[i] == 'e' || s[i] == 'E') {
        type = -1;
      } else if (!(s[i] == '+' || s[i] == '-')) {
        return 0; // non-numeric character
      }
    }
  }
  return type;
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
  //          << pad1.size() << " " << text.size()<< " " << pad2.size() << "\n"; //
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
std::string StringUtil::onespace(const std::string& s)
// Reduce spaces in string to single spaces
{
  std::vector<std::string> splitstring = StringUtil::split(s, " ");
  std::string sss;
  for (size_t k=0; k<splitstring.size(); k++) {
    sss += splitstring[k];
    if (k < splitstring.size()-1) {
      sss += " ";
    }
  }
  return sss;
}
//--------------------------------------------------------------
std::string StringUtil::removespaces(const std::string& s)
// Replace spaces in string with "_"
{
  std::vector<std::string> splitstring = StringUtil::split(s, " ");
  std::string sss;
  for (size_t k=0; k<splitstring.size(); k++) {
    sss += splitstring[k];
    if (k < splitstring.size()-1) {
      sss += "_";
    }
  }
  return sss;
}
//--------------------------------------------------------------
std::string StringUtil::valueSD(const double& v, const double& sd,
                                const int& totalfw,
                                const int& fw, const int& fd)
// Return value v with sd in brackets (converted to integer)
// totalfw total field width (if <=0, calculate),
// fw field width for value, fd number of decimals
{
  // scale to convert sd to integer
  double sdscale = 1.0;
  if (fd > 0) {
    for (int i=0;i<fd;++i) {sdscale *= 10.0;}
  }
  int nd = int(log10(sd*sdscale))+1;   // number of digits in sd
  bool decinsd = false;
  if (nd > fd) {
    nd += 1;
    decinsd = true;  // decimal point in sd
  }
  int tfw = std::max(totalfw, fw + 2 + nd); // set total width

  std::string s;
  if (decinsd) { // real sd
    s = StringUtil::Strip(StringUtil::ftos(sd, nd, fd));
  } else { // integer sd
    s = StringUtil::Strip(StringUtil::itos(Nint(sd*sdscale), nd));
  }
  s = "("+s+")";
  std::string ss = StringUtil::PadString(StringUtil::ftos(v, fw, fd)+s, tfw);
  return ss;
}
//--------------------------------------------------------------
// <tag><data</tag>
std::string StringUtil::MakeXMLtag(const std::string& tag, const std::string& data,
                                   const bool& edit)
{
  if (edit) {
    return "<"+tag+">"+XMLstring(data)+"</"+tag+">"; // edited to remove "<" characters etc
  } else {
    return "<"+tag+">"+data+"</"+tag+">"; // unedited
  }
}
//--------------------------------------------------------------
// <tag class="messageclass"><data</tag>
std::string StringUtil::MakeXMLwithclass(const std::string& tag, const std::string& data,
                                         const bool& edit, const std::string& messageclass)
{
  std::string content = data;
  if (edit) {
    content = XMLstring(data); // edited to remove "<" characters etc
  }
  return "<"+tag+" class=\""+messageclass+"\">\n"+content+"</"+tag+">";
}
//--------------------------------------------------------------
// <tag class="warningmessage"><data</tag>
std::string StringUtil::MakeXMLwarning(const std::string& tag, const std::string& data,
                                       const bool& edit)
{
  return MakeXMLwithclass(tag, data, edit, "warningmessage");
}
//--------------------------------------------------------------
//! make XML tag <tag>value</tag>
std::string StringUtil::MakeXMLtag(const std::string& tag, const int& value,
                                   const int& w)
{
  std::string s = StringUtil::Strip(StringUtil::itos(value, w));
  return MakeXMLtag(tag, s);
}
//--------------------------------------------------------------
//! make XML tag <tag>value</tag>
std::string StringUtil::MakeXMLtag(const std::string& tag,
                                   const double& value,
                                   const int& w, const int& d)
{
  std::string s = StringUtil::Strip(StringUtil::ftos(value, w,d));
  return MakeXMLtag(tag, s);
}
//--------------------------------------------------------------
std::string StringUtil::XMLstring(const std::string& s0)
//! return string modified to replace non-X/HTML characters &<>, append NL if long
{
  std::string s;
  for (size_t i=0;i<s0.size();++i) {
    if (s0[i] == '&') {s += "&amp;";}
    else if (s0[i] == '<') {s += "&lt;";}
    else if (s0[i] == '>') {s += "&gt;";}
    else {s += s0[i];}
  }
  const size_t MAXLENGTH = 100;
  if (s0.size() > MAXLENGTH) {
    s += "\n";
  }
  return s;
}
//--------------------------------------------------------------
std::string StringUtil::itos(const int f, const int w)
{ std::ostringstream s; s.width(w); s.setf(std::ios::fixed);s << f; return s.str(); }
//--------------------------------------------------------------
std::string StringUtil::itos(const int f)
{ std::ostringstream s;s << f; return s.str(); }
//--------------------------------------------------------------
std::string StringUtil::ftos(const float f, const int w, const int d)
{ std::ostringstream s; s.width(w); s.setf(std::ios::fixed); s.precision(d);s << f; return s.str(); }
//--------------------------------------------------------------
std::string StringUtil::ftos(const float f)
{ std::ostringstream s; s << f; return s.str(); }
//--------------------------------------------------------------
std::string StringUtil::etos(const float f, const int w, const int d)
{ std::ostringstream s; s.width(w); s.precision(d);s << f; return s.str(); }
//--------------------------------------------------------------
std::string StringUtil::ftos(const double f, const int w, const int d)
{ double ff = f;
  if (std::abs(ff) < std::numeric_limits<double>::min()) {
    ff = 0.0;
  }
  std::ostringstream s; s.width(w); s.setf(std::ios::fixed); s.precision(d);
  s << ff; return s.str();
}
//--------------------------------------------------------------
std::string StringUtil::ftos(const double f)
{ double ff = f;
  if (std::abs(ff) < std::numeric_limits<double>::min()) {
    ff = 0.0;
  }
  std::ostringstream s; s << ff; return s.str();
}
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
// Extract line from buffer, removing any leading or trailing
// Cr or Lf characters
{
  std::string line(buf);
  size_t ll = line.size();
  while (ll > 0) {
    if (line[ll-1] != '\n' && line[ll-1] != '\r') {break;}
    ll--;
  }
  size_t l0 = 0;
  while (l0 < ll) {
    if (line[l0] != '\n' && line[l0] != '\r') {break;}
    l0++;
  }
  return line.substr(l0,ll);
}
//--------------------------------------------------------------
std::string StringUtil::FormatSaveVector(const std::vector<int> ivec,
                                         const bool& finalNL)
// format integer vector for dump/save
// finalNL true to add newline at end
{
  std::string s = "";
  std::string line = "";
  const unsigned int MAXLINE = 100;
  for (size_t i=0;i<ivec.size();++i) {
    line += " "+clipper::String(ivec[i]);
    if (line.size() > MAXLINE && i < (ivec.size()-1)) {
      s += line+"\n";
      line = "";
    }
  }
  if (finalNL) {line += "\n";}
  return s+line;
}
//--------------------------------------------------------------
//! format real array for dump/save, each row delimited by "{}"

std::string StringUtil::FormatSaveArray(const clipper::Array2d<double>& VC)
{
  int nrows = VC.rows();
  int ncols = VC.cols();
  std::string s = "";

  for (int ir=0;ir<nrows;++ir) { // loop rows
    std::vector<double> v(ncols);
    for (int ic=0;ic<ncols;++ic) { // loop columns
      v[ic] = VC(ir,ic);
    }
    s += "{"+FormatSaveVector(v)+"}\n";
  }
  return s;
}
//--------------------------------------------------------------
std::string StringUtil::FormatSaveVector(const std::vector<double> vec,
                                         const bool& finalNL)
// format integer vector for dump/save
// finalNL true to add newline at end
{
  std::string s = "";
  std::string line = "";
  const unsigned int MAXLINE = 100;
  for (size_t i=0;i<vec.size();++i) {
    line += " "+clipper::String(vec[i]);
    if (line.size() > MAXLINE && i < (vec.size()-1)) {
      s += line+"\n";
      line = "";
    }
  }
  if (finalNL) {line += "\n";}
  return s+line;
}
//--------------------------------------------------------------
std::string StringUtil::FormatXMLcrossTable(const std::string& elementid, const std::string& tableid,
                                            const std::vector<std::string>& names,
                                            const std::string& valTag,
                                            const std::vector<std::pair<double,int> >& valCount)
// elementid      name for XML element
// tableid        id=tableid
// names          for each column/row
// valTag         XML tag string for value
// valCount       value and count
//  valCount array is in order:-
//    ab, ac, ad, ...
//        bc, bd, ...
//            cd, ...
{
  std::string s = "\n<"+Strip(elementid)+" id=\""+tableid+"\">\n";
  int nval = names.size();
  int npairs = nval*(nval-1)/2;
  ASSERT (int(valCount.size()) == npairs);

  // column headers, names from 1
  s += "<columnheaders>";
  for (int i=1;i<nval;++i) {
    s += StringUtil::MakeXMLtag("label", names[i]);
  }
  s += "\n</columnheaders>\n";
  std::string v;

  int k = 0;  // index to valCount
  for (int j=0;j<nval-1;++j) { // loop rows from 0 -> nval-2
    // row header
    s += "<row><label>"+names[j]+"</label>";
    for (int i=1;i<nval;++i) { // loop columns
      if (i <= j) {
        // blank (redundant) entry
        v = StringUtil::MakeXMLtag(valTag, " ");
        v += StringUtil::MakeXMLtag("Number", " ");
      } else {
        //      std::cout <<"i,j,k " <<i<<" "<<j<<" "<<k<<"\n";
        // values
        v = StringUtil::MakeXMLtag
          (valTag, StringUtil::ftos(valCount[k].first,7,3));
        v += StringUtil::MakeXMLtag("Number", StringUtil::itos(valCount[k].second,5));
        k++;
      }
      s += v;
    } // end column loop
    s += "</row>\n";
  } // end row loop
  s += "</"+Strip(elementid)+">";
  return s;
}
//======================================================================
//! just add leading tabs to string and newline if not there already
std::string FormatOutput::logTab(const int& tab, const std::string& text,
                                 const bool& add_return)
{
  std::string nl;
  if (add_return && (text[text.size()-1] != '\n')) nl = "\n";
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
  vsnprintf(temp,8192,formattext,arglist);
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
    //^       <<Precision<<" "<<prec
    //^       <<" "<<l<<" "<<dec <<" "<<label1<<" "<<label2<<"\n";
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
//======================================================================
//! construct with one citation
Citation::Citation(const std::string& citation,
                   const std::string& link)
{
  citations.assign(1, citation);
  links.assign(1,link);
}
//--------------------------------------------------------------
//! add a citation
void Citation::AddCitation(const std::string& citation,
                           const std::string& link)
{
  citations.push_back(citation);
  links.push_back(link);
}
//--------------------------------------------------------------
//! make citation string for log file, with html link
std::string Citation::MakeLogCitation() const
{
  std::string s = "$TEXT:Reference: $$ Please cite $$\n";
  for (size_t i=0; i<citations.size(); i++) {
    s += citations[i]+"\n";
    if (links[i] != "") {
      s += "<a href=\""+links[i]+"\">\n";
      s += "<b>PDF</b></a>\n";
    }
  }
  s += "$$\n";
  return s;
}
//--------------------------------------------------------------
//! make citation string for XML file, with html link  FIXME
// std::string Citation::MakeXMLCitation() const;
