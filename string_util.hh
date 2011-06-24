#ifndef STRING_UTIL
#define STRING_UTIL

#include <string>
#include <vector>
#include "range.hh"

//======================================================================
class StringUtil {
  //! String utilities
public:
//--------------------------------------------------------------
//! Strip out space and x characters from string
// x defaults to null
static std::string Strip(const std::string& s, const char& x = '\0');
//--------------------------------------------------------------
//! Strip out null characters from string
static std::string StripNull(const std::string& s);
//--------------------------------------------------------------
//! Trim off leading & trailing spaces from string
static std::string Trim(const std::string& s);
//--------------------------------------------------------------
//! Return string of length <fieldwidth> with text centred on position cenpos in it
/*! If string length is > fieldwidth, return full string */
static std::string CentreString(const std::string& text, const int& fieldwidth,
			 const int& cenpos=0);
//--------------------------------------------------------------
//! Return string of length Max(fieldwidth, Length(text))
static std::string PadString(const std::string& text, const int& fieldwidth);
//--------------------------------------------------------------
//! Return right-justified string
static std::string RightString(const std::string& text, const int& fieldwidth);
//--------------------------------------------------------------
//! Return left-justified string
static std::string LeftString(const std::string& text, const int& fieldwidth);
//--------------------------------------------------------------
//! Return substrings split at string "sep" (excluded)
static std::vector<std::string> split(const std::string& str,
			       const std::string& sep);
//--------------------------------------------------------------
//! Return substrings split at string "sep1" or "sep2" (excluded)
static std::vector<std::string> split(const std::string& str,
			   const std::string& sep1, const std::string& sep2);
//--------------------------------------------------------------
//! make XML tag <tag><data</tag>
  static std::string MakeXMLtag(const std::string& tag, const std::string& data);
//--------------------------------------------------------------
//! int to string conversion (%wd)
  static std::string itos(const int f, const int w);
//--------------------------------------------------------------
//! fixed point float to string conversion (%w.df)
   static std::string ftos(const float f, const int w, const int d);
//--------------------------------------------------------------
//! floating point float to string conversion (%w.df)
   static std::string etos(const float f, const int w, const int d);
//--------------------------------------------------------------
//! fixed point double to string conversion (%w.df)
   static std::string ftos(const double f, const int w, const int d);
//--------------------------------------------------------------
//! floating point double to string conversion (%w.df)
   static std::string etos(const double f, const int w, const int d);
//--------------------------------------------------------------
//! returns wrapped line with indent
/*! 
 \param line       line to wrap
 \param pagewidth  number of characters in wrapped line (may run over at end)
 \param indent     indent lines after first by this number of items in 1st line
 \param sepc       terminator for item, default space, included if not space
*/
  static std::string WrapLine(const std::string& line,
    const int& pagewidth, const int& nindent, const std::string& sepc=" ");

//--------------------------------------------------------------
  static std::string formatFraction(const double& fr, const int& width);
  //!< Format a real number as a fraction ie "n/m"
  /*! \param  fr  number of format
    \param  width  field width for decimal version, if fail to find suitable fraction
  */
  //--------------------------------------------------------------
  static std::string ToUpper(const std::string& s);
  //!< returns uppercase version of string
  //--------------------------------------------------------------
  //! Extract line from buffer, removing any trailing Cr or Lf characters
  static std::string BuftoLine(const std::string buf);
  //--------------------------------------------------------------
  //! format integer vector for dump/save
  static std::string FormatSaveVector(const std::vector<int> ivec);
  //--------------------------------------------------------------
  //! format real vector for dump/save
    static std::string FormatSaveVector(const std::vector<double> vec);
};
//======================================================================
//! Format output in a similar way to the phaser Output class, but just return as a string
class FormatOutput {
public:
  //! just add leading tabs to string and newline if not there already
  /*! usually redundant */
  static std::string logTab(const int& tab, const std::string& text);
  //! format using sstringf
  static std::string logTabPrintf(const int& tab,
				  const char* formattext,...);
  static std::string logWarning(const std::string& text);
};
//======================================================================
//! Information for formatting a number field
class Numberfield {
public:
  Numberfield():type(-1),width(-1) {}
  //! constructor from type, field width, and number of characters after decimal point
  /*! \param Type   = 0 float, +1 int, +2 string, -1 item type (float or int)
      \param Width  field width
      \param Dec    number of characters after decimal point
      \param label1 column title for 1st line of header
      \param label2 column title for 2nd line of header (default blank)
   */
  Numberfield(const int& Type, const int& Width, const int& Dec,
	      const std::string& Label1="", const std::string& Label2="")
    :type(Type),width(Width),dec(Dec),label1(Label1),label2(Label2){}

  //! (re)initialise from type, field width, and number of characters after decimal point
  void init(const int& Type, const int& Width, const int& Dec,
	    const std::string& Label1="", const std::string& Label2="");

  //! constructor from type, minimum field width, and precision
  /*! \param IntType    true if integer type
      \param MaxValue   maximum ||value|| for field
      \param MinWidth   minimum field width
      \param Precision  number of significant figures
   */
Numberfield(const bool& IntType, const float& MaxValue,
	    const int& MinWidth, const int& Precision);

  //! (re)initialise from type, minimum field width, and precision
  /*! \param IntType    true if integer type
      \param MaxValue   maximum ||value|| for field
      \param MinWidth   minimum field width
      \param Precision  number of significant figures
   */
void init(const bool& IntType, const float& MaxValue,
	  const int& MinWidth, const int& Precision);

  int type;  // = 0 float, +1 int, +2 string, -1 item type
  int width; // field width
  int dec;   // digits after point (type 0)
  scala::Range range;
  std::string label1; // labels for column heading rows
  std::string label2;
};


#endif
