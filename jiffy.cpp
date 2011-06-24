#include "jiffy.hh"
#include <complex>
#include <string>
#include <cstdlib>
#include <ctype.h>

namespace phaser_io {

//-----------------------------------------------------------------------
// itos - int to std::string
//-----------------------------------------------------------------------
std::string itos(const int & i)
{
  std::string str;
  std::stringstream s;
  s << i;
  s >> str;
  return str;
}

//-----------------------------------------------------------------------
// ftos - float to std::string
//-----------------------------------------------------------------------
std::string ftos(const float & f)
{
  std::string str;
  std::stringstream s;
  s << f;
  s >> str;
  return str;
}

//-----------------------------------------------------------------------
// dtos - double/float to std::string
//-----------------------------------------------------------------------
std::string dtos(const float & f)
{
  std::string str;
  std::stringstream s;
  s << f;
  s >> str;
  return str;
}

std::string dtos(const double & f)
{
  std::string str;
  std::stringstream s;
  s << f;
  s >> str;
  return str;
}

//-----------------------------------------------------------------------
// ctos - complex to std::string
//-----------------------------------------------------------------------
std::string ctos(const std::complex<float> & c)
{ return "(" + dtos(std::real(c)) + "," +  dtos(std::imag(c)) + ")"; }

std::string ctos(const std::complex<double> & c)
{ return "(" + dtos(std::real(c)) + "," +  dtos(std::imag(c)) + ")"; }

//-----------------------------------------------------------------------
// btos - bool to std::string
//-----------------------------------------------------------------------
std::string btos(const bool& i)
{
  if (i) return "true ";
  return "false";
}

//-----------------------------------------------------------------------
// stoa - std::string to alpha
//-----------------------------------------------------------------------
std::string stoa(const std::string & str)
{ 
  std::stringstream s;
  std::string t;
  char ch;
  s << str;
  do { if (!s.get(ch)) return "";
     } while (!isalpha(ch)); //get a-z,A-Z
  t = ch;
  while(s.get(ch) && isalpha(ch)) t += ch;
  return t;
} 

//-----------------------------------------------------------------------
// stod - std::string to digit
//-----------------------------------------------------------------------
int stod(const std::string & str)
{ 
  std::stringstream s;
  char ch;
  std::string t;
  s << str;
  do { if (!s.get(ch)) return 0;
     } while (!isdigit(ch)); //get digits
  t = ch;
  while(s.get(ch) && isdigit(ch)) t += ch;
  return std::atoi(t.c_str());
}

//-----------------------------------------------------------------------
// vtos - vector to std::string
//-----------------------------------------------------------------------
std::string ivtos(const scitbx::vec3<int> & f)
{ return itos(f[0]) + " " + itos(f[1]) + " " + itos(f[2]) ; }

std::string dvtos(const scitbx::vec3<double> & f)
{ return dtos(f[0]) + " " + dtos(f[1]) + " " + dtos(f[2]) ; }

std::string dvtos(const scitbx::vec3<float> & f)
{ return dtos(f[0]) + " " + dtos(f[1]) + " " + dtos(f[2]) ; }

//-----------------------------------------------------------------------
// mtos - matrix to std::string
//-----------------------------------------------------------------------
std::string imtos(const scitbx::mat3<int> & f)
{
  return itos(f(0,0)) + " " + itos(f(0,1)) + " " + itos(f(0,2)) + " "
       + itos(f(1,0)) + " " + itos(f(1,1)) + " " + itos(f(1,2)) + " "
       + itos(f(2,0)) + " " + itos(f(2,1)) + " " + itos(f(2,2));
}
std::string dmtos(const scitbx::mat3<float> & f)
{
  return dtos(f(0,0)) + " " + dtos(f(0,1)) + " " + dtos(f(0,2)) + " "
       + dtos(f(1,0)) + " " + dtos(f(1,1)) + " " + dtos(f(1,2)) + " "
       + dtos(f(2,0)) + " " + dtos(f(2,1)) + " " + dtos(f(2,2));
}
std::string dmtos(const scitbx::mat3<double> & f)
{
  return dtos(f(0,0)) + " " + dtos(f(0,1)) + " " + dtos(f(0,2)) + " "
       + dtos(f(1,0)) + " " + dtos(f(1,1)) + " " + dtos(f(1,2)) + " "
       + dtos(f(2,0)) + " " + dtos(f(2,1)) + " " + dtos(f(2,2));
}

//-----------------------------------------------------------------------
// itoaniso - anisotropic integer to  std::string
//-----------------------------------------------------------------------
std::string itoaniso(const int & i)
{
  if (i == 0) return " HH ";
  if (i == 1) return " KK ";
  if (i == 2) return " LL ";
  if (i == 3) return " HK ";
  if (i == 4) return " HL ";
  if (i == 5) return " KL ";
  return " -- ";
}

//-----------------------------------------------------------------------
// stoup - std::string to uppercase
//-----------------------------------------------------------------------
std::string stoup(const std::string & str)
{ 
  char ch;
  std::stringstream s;
  std::string t="";
  s << str;
  while (s.get(ch))
    t += toupper(ch); //convert to uppercase
  return t;
} 

//-----------------------------------------------------------------------
// isfloat - is the std::string a float?
//-----------------------------------------------------------------------
bool isfloat(std::string ch)
{
  std::string separators(" \t");
  std::size_t start = ch.find_first_not_of(separators);
  std::size_t stop = ch.find_last_not_of(separators);
  ch = ch.substr(start,stop-start);

  //first character is different, can't start with exponent 
  if (!ch.size()) return false;
  if (!((isdigit)(ch[0]) || ch[0]=='.' || ch[0]=='+'|| ch[0]=='-' ))
    return false;
//rest of string is all digits or . or e or E or + or -(exponent)
  for (size_t i = 0; i < ch.size(); i++)
    if (!((isdigit)(ch[i]) || ch[i]=='.' || ch[i]=='e' || ch[i]=='E' || ch[i]=='+'|| ch[i]=='-'))
      return false;
  return true;
}

//-----------------------------------------------------------------------
// fmod_pos - return modulus if pos or neg
//-----------------------------------------------------------------------
//floatType fmod_pos(floatType x, floatType y)
//{
//  if (x < 0) return fmod(x+1,y) + y - 1;
//  return fmod(x,y);
//} 

}
