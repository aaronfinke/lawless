// imagearray.cpp

#define ASSERT assert
#include <assert.h>
#include <string.h>
#include "imagearray.hh"
#include "file_util.hh"
#include "string_util.hh"

using clipper::Message;
using clipper::Message_fatal;

//--------------------------------------------------------------
Imagearray::Imagearray()
//! construct empty
  : imgsize1(-1), scale(-1.0), dataread(0)
{
  lenheader = 512;
}
//--------------------------------------------------------------
Imagearray::Imagearray(const clipper::Array2d<double>& array)
//! construct from 2D array
  : imgsize1(-1), scale(-1.0), dataread(0)
{
  init(array);
}
//--------------------------------------------------------------
void Imagearray::init(const clipper::Array2d<double>& array)
//! initialise from 2D array
{
  lenheader = 512;
  title = "ADSC format image";
  size1 = array.rows();
  size2 = array.cols();
  if (imgsize1 < 0) {
    imgsize1 = size1; // by default
    imgsize2 = size2;
  }
  if (scale < 0.0) {scale = 1000.0;}
  data.resize(array.size());
  dataread = +1;  // data contains scale*value

  for (int i=0;i<size1;++i) { // loop x
    for (int j=0;j<size2;++j) { // loop y
      ASSERT (array(i,j) >= 0.0);
      uint16_t d = uint16(scale*array(i,j));
      size_t k = i*size2 + j;
      data[k] = d;
    }}
}
//--------------------------------------------------------------
void Imagearray::init(const clipper::Array2d<scala::MeanSD>& array, const bool& Mean)
//! initialise from 2D array of MeanSD, use mean or SD, if mean = true or false
{
  lenheader = 512;
  title = "ADSC format image, mean deviation";
  size1 = array.rows();
  size2 = array.cols();
  if (imgsize1 < 0) {
    imgsize1 = size1; // by default
    imgsize2 = size2;
  }
  if (scale < 0.0) {scale = 1000.0;}
  data.resize(array.size());
  dataread = +1;  // data contains scale*value
  double value;

  for (int i=0;i<size1;++i) { // loop x
    for (int j=0;j<size2;++j) { // loop y
      if (Mean) {
	value = array(i,j).Mean();
      } else {
	value = array(i,j).SD();
      }
      uint16_t d = uint16_t(scale*value);
      size_t k = i*size2 + j;
      data[k] = d;
    }}
}
//--------------------------------------------------------------
// write out
void Imagearray::Write(const std::string& filename) const
{
  // output binary file
  FILE* opfile = OpenFile(filename, true, true);
  std::string header = MakeHeader();
  char hd[lenheader];
  strcpy(hd, std::string(lenheader,' ').c_str());
  strcpy(hd, header.c_str());
  int n = fwrite(hd, sizeof(char), lenheader, opfile);
  if (n < header.size()) {
    Message::message(Message_fatal
		     ("Imagearray::Write failed to write header"));
  }
  n = fwrite(&*data.begin(), sizeof(uint16_t), data.size(), opfile);
  if (n < data.size()) {
    Message::message(Message_fatal
		     ("Imagearray::Write failed to write data"));
  }
  fclose(opfile);
}
//--------------------------------------------------------------
  std::string Imagearray::MakeHeader() const {
  // ADSC header
  //   {
  //   HEADER_BYTES=  512;
  //   COMMENT=Written by xds2cor;
  //   DIM=2;
  //   SIZE1=384;
  //   SIZE2=384;
  //   DETECTOR_SN=918;
  //   DATE=Tue Aug 15 18:01:00 2006
  //   TYPE=unsigned_short;
  //   BYTE_ORDER=little_endian;
  //   SCALE=4096;
  //   BINNING=8;
  //   IMG_PIXEL_SIZE=0.051294;
  //   IMG_SIZE1=6144;
  //   IMG_SIZE2=6144;
  //   }^L
  std::string header = "{\n";
  header += HeaderLineInt("HEADER_BYTES",lenheader,4);
  header += "COMMENT="+title+";\n";
  header += "DIM=2;\n";
  header += HeaderLineInt("SIZE1", size1, 6);
  header += HeaderLineInt("SIZE2", size2, 6);
  header += std::string("DETECTOR_SN=1;\n")+
    "TYPE=unsigned_short;\n"+
    "BYTE_ORDER=little_endian;\n";
  header += HeaderLineInt("SCALE",int(scale),10);
  header += std::string("BINNING=1;\n")+
    "IMG_PIXEL_SIZE=0.10;\n";
  header += HeaderLineInt("IMG_SIZE1", imgsize1, 6);
  header += HeaderLineInt("IMG_SIZE2", imgsize2, 6);
  header += "DISTANCE=100.;\n";
  header += "}\f";
  return header;
}
//--------------------------------------------------------------
std::string Imagearray::HeaderLineReal(const std::string& label,
				      const double& value,
				      const int& fw,const int& dw) const
{
  // format real number line
  return StringUtil::Strip(label+"="+StringUtil::ftos(value,fw,dw)+";\n");
}
//--------------------------------------------------------------
std::string Imagearray::HeaderLineInt(const std::string& label,
				     const int& value,
				     const int& fw) const
{
  // format integer number line
  return StringUtil::Strip(label+"="+StringUtil::itos(value,fw)+";\n");
}
//--------------------------------------------------------------
void Imagearray::ReadFile(const std::string& filename)
//! Read ADSC image file
{
  FILE* infile = OpenFile(filename, false, true); // binary read
  char cheader[lenheader+1];  // buffer for header
  int nread = fread(cheader, sizeof(char), lenheader, infile);
  ASSERT (nread == lenheader);
  cheader[lenheader] = '\0';

  clipper::String header(cheader);
  //  std::cout << header <<"\n";
  // Split on "\n"
  std::vector<clipper::String> lines = header.split("\n");
  for (size_t l=0;l<lines.size();++l) {
    // Token lines contain "="
    if (lines[l].find("=") != std::string::npos) {
      //^      std::cout << "\n" << lines[l] << "\n";
      std::vector<clipper::String> tokens = lines[l].split("=");
      //^      for (size_t i=0;i<tokens.size();++i) {std::cout <<" "<<tokens[i];}
      if (tokens[0] == "SIZE1") {
	size1 = tokens[1].i();
      } else if (tokens[0] == "SIZE2") {
	size2 = tokens[1].i();
      } else if (tokens[0] == "SIZE2") {
	size2 = tokens[1].i();
      } else if (tokens[0] == "IMG_SIZE1") {
	imgsize1 = tokens[1].i();
      } else if (tokens[0] == "IMG_SIZE2") {
	imgsize2 = tokens[1].i();
      } else if (tokens[0] == "SCALE") {
	  scale = tokens[1].f64();
      }
    }
  } // end loop lines

  // Now read the image
  int npxl = size1*size2;
  data.resize(npxl);
  nread = fread(&*data.begin(), sizeof(uint16_t), npxl, infile);
  ASSERT (nread == npxl);
  dataread = -1; //  = -1 data contains scale/value
}
//--------------------------------------------------------------
clipper::Array2d<double> Imagearray::GetImage(const bool& invert) const
//! Return data array as
//  (a) value from scale/value if invert = false
//  (b) 1/value from scale/value if invert = true
{
  clipper::Array2d<double> array;
  if (dataread == 0) {
    // no data
    return array;
  }
  array.resize(size1,size2);
  for (int i=0;i<size1;++i) { // loop x
    for (int j=0;j<size2;++j) { // loop y
      size_t k = i*size1 + j;
      if (!invert && dataread < 0) { // data =scale/value
	array(i,j) = scale/data[k];
      } else { // data = scale*value
	array(i,j) = data[k]/scale;
      }
    }}
  return array;
}
//--------------------------------------------------------------
