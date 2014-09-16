// intensitybin.hh
//
// Intensity bin class

#ifndef INTENSITYBIN_HEADER
#define INTENSITYBIN_HEADER

#include "range.hh"
///#include "score_datatypes.hh"
#include "scala_util.hh"

namespace scala
{
  // -----------------------------------------------------------
  class IntensityBin
  {
    // Intensity binning based on exponential Wilson distribution
    /*
      Wilson distribution
        f(I) = p(I) dI = (1/pi.S) exp ( -I/S)
	 where S = SigmaN

        then f(0) = 1/pi.S

      Dividing f(I) into N divisions we have for the j'th limit
        (j = 0,N-1)

      p(I)j = (1/pi.S)[1 - (j+1)/N] =  (1/pi.S) exp ( -I/S)
      ie   ln(1-(j+1)/N) = -I/S

      If Ik is the value of I at j = jk (top of reference bin)
      S = - Ik / (ln(1 - (jk+1)/N))
      & Ij = - S ln(1-(j+1)/N)  (top of j'th bin, j = 0,Nbin-1)
    */
  public:
    IntensityBin(){}

    // Constructor
    //  Nbin    number of bins
    //  jk      reference bin number (0 < jk < Nbin-1)
    //  Ik      intensity at top of reference bin jk
    //  Imax    maximum possible intensity (top of (Nbin-1)'th bin)
    IntensityBin(const int& Nbin, const int& jk, const float& Ik,
		 const float& Imax);

    // Return bin number, update ranges
    int bin(const float& I);

    // Return number of bins
    int NumberBins() const {return nbin;}

    // Middle of bin
    float middle(const int& bin) const;
    // Mean
    float mean(const int& bin) const;
    // limits of bin (actual limits if any data, else range)
    RPair bounds(const int& bin) const;
    // Number in bin
    int Count(const int& bin) const;
    // Clear counts
    void Clear();


  private:
    int nbin;
    std::vector<float> bintop;
    std::vector<Range> binrange;  // actual ranges in each bin
    std::vector<MeanSD> binmean;
    std::vector<int> count;
  };

}
#endif
