#ifndef __PHASER_JIFFY__
#define __PHASER_JIFFY__
//#include <phaser/Phaser.h>
#include <sstream> 
#include <string>

#include "phaser_types.hh"

using namespace phaser;

namespace phaser_io {

std::string dtos(const float&);
std::string dtos(const double&);
std::string dvtos(const scitbx::vec3<float>&);
std::string dvtos(const scitbx::vec3<double>&);
std::string ivtos(const scitbx::vec3<int>&);
std::string dmtos(const scitbx::mat3<float>&);
std::string dmtos(const scitbx::mat3<double>&);
std::string imtos(const scitbx::mat3<int>&);
std::string ftos(const float&);
std::string itos(const int&);
std::string itoaniso(const int&);
std::string btos(const bool&);
std::string ctos(const std::complex<float>&);
std::string ctos(const std::complex<double>&);
std::string stoup(const std::string&);
int         stod(const std::string&);
std::string stoa(const std::string&);

bool isposi(const double &);
bool isfloat(std::string);
  //floatType fmod_pos(floatType, floatType);

}//end namespace phaser

#endif

