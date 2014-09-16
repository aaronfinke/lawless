// eprob.hh
//
// Class to deal with probabilities of normalised structure factor E
//
// p(E) = exp(-E^2)        acentric
//      = erfc(E/sqrt(2))  centric
//

#ifndef EPROB_HEADER
#define EPROB_HEADER

#include <string>

class EProb
{
public:
  EProb();
  EProb(const float& Emax);

  //! initialise from acentric Emax
  void init(const float& Emax);

  //! clear, ie flag as no check
  void Clear();

  //! return true if no test
  bool Null() const {return (emaxacen <= 0.0);}

  //! return true if E > limits
  bool TooBig(const float& E2, const bool& Centric) const;

  std::string format() const;

private:
  float emaxacen;  // acentric
  float emaxcentric;
  float emaxacen2;  // squares
  float emaxcentric2;

  static const float MAXIMUM_EMAX;
  static const float MAXIMUM_EMAX_CENTRIC;
  static const float DEFAULT_EMAX;

  // Get centric Emax from acentric value,
  // using binary search (cf fortran routine fndepb in Scala)
  float CentricEmax(const float& emaxacen) const;

};
#endif
