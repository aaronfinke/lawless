// eprob.cpp
//
// Class to deal with probabilities of normalised structure factor E
//
// p(E) = exp(-E^2)        acentric
//      = erfc(E/sqrt(2))  centric
//

#include <iostream>
#include <cmath>

#include "eprob.hh"
#include "util.hh"
#include "string_util.hh"

// ------------------------------------------------------------
const float EProb::MAXIMUM_EMAX = 15.0;
const float EProb::MAXIMUM_EMAX_CENTRIC = 22.0;
const float EProb::DEFAULT_EMAX = 10.0;
// for large negative E, Min(E) = -NEGATIVE_EMAX_RATIO * Emax
const float NEGATIVE_EMAX_RATIO = 0.5;
const float NEGATIVE_EMAX2_RATIO = NEGATIVE_EMAX_RATIO * NEGATIVE_EMAX_RATIO;
// ------------------------------------------------------------
EProb::EProb()
{
  init(DEFAULT_EMAX);
}
// ------------------------------------------------------------
//! initialise from acentric Emax
 EProb::EProb(const float& Emax)
{
  init(Emax);
}
// ------------------------------------------------------------
//! initialise from acentric Emax
void EProb::init(const float& Emax)
{
  if (Emax <= 0.0) {
    Clear();
  } else {
    emaxacen = Min(Emax, MAXIMUM_EMAX);
    emaxcentric = CentricEmax(emaxacen); // centric Emax from acentric
    emaxacen2 = emaxacen*emaxacen;
    emaxcentric2 = emaxcentric*emaxcentric;
  }
}
// ------------------------------------------------------------
//! clear, ie flag as no check
void EProb::Clear()
{
  emaxacen = -1.0; emaxcentric = -1.0;
  emaxacen2 = -1.0; emaxcentric2 = -1.0;
}
// ------------------------------------------------------------
float EProb::CentricEmax(const float& emaxacen) const
// Get centric Emax from acentric value,
// using binary search (cf fortran routine fndepb in Scala)
{
  double eacen = Min(emaxacen, MAXIMUM_EMAX); // starting value
  double stope = eacen*eacen; // upper bracket
  double prob = exp(-stope);  // target probability
  double start = 0.0;         // lower bracket
  double step  = 1.0;         // step
  double p = 0.0;
  double rt2 = sqrt(2.0);
  const double TOLERANCE = 0.00001;
  const double FAC = 3.0;
  double e = start;

  while (step > TOLERANCE) {
    p = erfc(e/rt2);
    if (p < prob) {
      // reduce step
      start = e - step;
      stope = e;
      step  = step/FAC;
      e = start;
    } else {
      // increase step
      e = e + step;
    }
  }
  // If we fall out of the end, this answer is near enough
  double ecen = Min(e - 0.5*step, MAXIMUM_EMAX_CENTRIC);
  return float(ecen);
}
// ------------------------------------------------------------
//! return true if E > limits
bool EProb::TooBig(const float& E2, const bool& Centric) const
{
  if (E2 > 0.0) { // normal positive E^2
    if (Centric) {
      return (E2 > emaxcentric2);
    } else {
      return (E2 > emaxacen2);
    }
  } else { // test for very negative E^2
    if (Centric) {
      return (E2 < -NEGATIVE_EMAX2_RATIO*emaxcentric2);
    } else {
      return (E2 < -NEGATIVE_EMAX2_RATIO*emaxacen2);
    }
  }
}
// ------------------------------------------------------------
//! return true if E > limits, multiplied by factor
bool EProb::TooBig(const float& E2, const bool& Centric,
                   const float& factor) const
{
  if (E2 > 0.0) { // normal positive E^2
    if (Centric) {
      return (E2 > emaxcentric2*factor);
    } else {
      return (E2 > emaxacen2*factor);
    }
  } else { // test for very negative E^2
    if (Centric) {
      return (E2 < -NEGATIVE_EMAX2_RATIO*emaxcentric2*factor);
    } else {
      return (E2 < -NEGATIVE_EMAX2_RATIO*emaxacen2*factor);
    }
  }
}
// ------------------------------------------------------------
std::string EProb::format() const
{
  if (emaxacen > 0.0) {
    std::string s ="Reflections judged implausibly large will be rejected\n";
    s += "     Maximum and minimum normalised F (ie E) for acentric reflection"+
      StringUtil::ftos(emaxacen,10,2)+", "+
      StringUtil::ftos(-NEGATIVE_EMAX_RATIO*emaxacen,10,2)+"\n";
    s += "     Maximum and minimum normalised F (ie E) for centric reflection "+
      StringUtil::ftos(emaxcentric,10,2)+", "+
      StringUtil::ftos(-NEGATIVE_EMAX_RATIO*emaxcentric,10,2)+"\n";
    //    float p = exp(-emaxacen*emaxacen);
    //    s += "     Minimum probability before reflection is rejected  "+
    //      StringUtil::etos(p,10,3)+"\n";
    return s;
  } else {
    return "No maximum E test\n";
  }
}
// ------------------------------------------------------------
std::string EProb::formatXML() const
{
  std::string s = "\n<EmaxTest>";
  if (emaxacen > 0.0) {
    s += StringUtil::MakeXMLtag("EmaxAcentric", StringUtil::ftos(emaxacen,7,2));
    s += StringUtil::MakeXMLtag("EmaxCentric", StringUtil::ftos(emaxcentric,7,2));
  } else {
    s += StringUtil::MakeXMLtag("EmaxAcentric", -1.0);
  }
  s += "\n</EmaxTest>\n";
  return s;
}
// ------------------------------------------------------------
