//
// initialscales.hh
//
// Initial rough scaling by making average intensities equal
//

#ifndef INITIALSCALES_HEADER
#define INITIALSCALES_HEADER

#include "hkl_unmerge.hh"
#include "scalemodel.hh"
#include "controls.hh"
#include "Output.hh"


typedef std::pair<float,float> FPair;
typedef std::pair<double,double> DPair;

namespace scala {

class InitialData 
// Class to store initial scaling data, for FoxHolmes scaling
{
public:
  InitialData();
  InitialData(const clipper::Array2d<double>& AvI);

  int Npar() const {return np;}    // N valid rotation ranges 
  int NparAll() const {return npall;} // total N rotation ranges

  // return for each rotation range true if there are data, else false
  std::vector<bool> validRanges() const {return validranges;}

  // Data are in 2D array AvI(rotation, resolution)
  // For each resolution bin, we want to make all the <Irot> equal over
  //   all rotation bins  ie <I(rot,reso)>/scales(rot)  constant for all "rot"

  // return one observation, I, sigma pair, length np
  // Return false if end of data
  bool ObsArray(std::vector<DPair>& obs) const;

private:
  int npall;  // Number of parameters = nrotranges
  int np;     // Number of active parameters <= nrotranges
  // Data is a 2D array(nrotranges, nresbins)
  const clipper::Array2d<double>*  avi;
  std::vector<bool> validranges;
  mutable int next;
};


  // Get initial estimates of primary scales, from making intensity
  // averages equal
  void InitialScales(hkl_unmerge_list& hkl_list, ScaleModel& AllScales,
  		     const all_controls& controls,
		     phaser_io::Output& output);

}

#endif
