// Class for reading or writing image array
//  at present just ADSC image format

#ifndef IMAGEARRAY_HEADER
#define IMAGEARRAY_HEADER

#include <stdint.h>

// Clipper
#include <clipper/clipper.h>

#include "scala_util.hh"

class Imagearray {
public:
  Imagearray(); // construct as empty
  //! construct from 2D array
  Imagearray(const clipper::Array2d<double>& array);
  //! initialise from 2D array
  void init(const clipper::Array2d<double>& array);
  //! initialise from 2D array of MeanSD, use mean or SD, if mean = true or false
  void init(const clipper::Array2d<scala::MeanSD>& array, const bool& Mean);

  //! Set image size
  void SetImageSize(const int& sz1, const int& sz2)
  {imgsize1=sz1;imgsize2=sz2;}

  //! Set scale
  void SetScale(const double& scl)
  {scale=scl;}

  // write out
  void Write(const std::string& filename) const;

  //! Read ADSC image file
  void ReadFile(const std::string& filename);

  //! Return data array as
  //  (a) value from scale/value if invert = false
  //  (b) 1/value from scale/value if invert = true
  clipper::Array2d<double> GetImage(const bool& invert) const;

private:
  // Header information
  int lenheader;
  std::string title;
  // Array size (may be smaller than image size)
  int size1, size2;
  // Image size (may be larger than array size)
  int imgsize1, imgsize2;
  // Scale factor for values val = scale/data
  double scale;

  std::vector<uint16_t> data; // the data array as short unsigned integers
  int dataread;  // = 0 empty, = -1 data contains scale/value
		 // = +1 data contains scale * value

  // format header lines

  std::string MakeHeader() const;

  std::string HeaderLineReal(const std::string& label,
			     const double& value,
			     const int& fw,const int& dw) const;
  
  std::string HeaderLineInt(const std::string& label,
			    const int& value,
			    const int& fw) const;

};


#endif
