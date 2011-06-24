// intensitybin.cpp

#include "intensitybin.hh"
// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;


namespace scala 
{
  //--------------------------------------------------------------
    // Constructor
    //  Nbin    number of bins
    //  jk      reference bin number (0 < jk < Nbin-1)
    //  Ik      intensity at top of reference bin jk
    //  Imax    maximum possible intensity (top of (Nbin-1)'th bin)
  IntensityBin::IntensityBin(const int& Nbin, const int& jk, const float& Ik,
			     const float& Imax)
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
  {    
    nbin = Nbin;
    float S = -Ik/log(1.-float(jk+1)/float(Nbin));
    bintop.resize(nbin);
    binrange.resize(nbin);
    binmean.resize(nbin);
    count.resize(nbin);

    for (int i=0;i<nbin-1;i++) {
      bintop[i] = -S * log(1.0 - float(i+1)/float(nbin));
      binrange[i].clear();
      binmean[i].clear();
      count[i] = 0;
      //^	std::cout << "i, Itop" << i << " " << bintop[i] << "\n";
    }
    bintop[nbin-1] = Imax;
    //^    std::cout << "i, Ibin " << nbin-1 << " " << bintop[nbin-1] << "\n";
    //^-
  }
  //--------------------------------------------------------------
  // Return bin number, update ranges
  int IntensityBin::bin(const float& I)
  {
    int j = -1;
    for (int i=0;i<nbin;i++)  {
      if (I < bintop[i]) {
	j = i;
	break;
      }
    }
    if (j < 0) j = nbin-1;
    binrange[j].update(I);   // update min & max
    binmean[j].Add(I);
    count[j]++;  // count
    return j;
    }
  //--------------------------------------------------------------
  float IntensityBin::middle(const int& bin) const
  {
    if (count.at(bin) > 0) {
      return 0.5*(binrange[bin].min() + binrange[bin].max());
    }
    // Range limits
    float lower = 0.0;
    if (bin > 0) lower = bintop[bin-1];
    return 0.5*(lower + bintop[bin]);
  }
  //--------------------------------------------------------------
  float IntensityBin::mean(const int& bin) const
  {
    if (count.at(bin) > 0) {
      return binmean[bin].Mean();
    }
    // Range limits
    float lower = 0.0;
    if (bin > 0) lower = bintop[bin-1];
    return 0.5*(lower + bintop[bin]);
  }
  //--------------------------------------------------------------
  // limits of bin (actual limits if any data, else range)
  RPair IntensityBin::bounds(const int& bin) const
  {
    if (count.at(bin) > 0) {
      // real data limits
      return RPair(binrange[bin].min(), binrange[bin].max());
    }
    // Range limits
    float lower = 0.0;
    if (bin > 0) lower = bintop[bin-1];
    return RPair(lower, bintop[bin]);
  }
  //--------------------------------------------------------------
  // Number in bin
  int IntensityBin::Count(const int& bin) const
  {
    return count.at(bin);
  }
  //--------------------------------------------------------------
  // Clear counts
  void IntensityBin::Clear()
  {
    for (size_t i=0;i<count.size();i++) {
      count[i] = 0;
      binrange[i].clear();
    }
  }
  //--------------------------------------------------------------
}
