//
// tie.hh
//
// Tie things


#ifndef TIE_HEADER
#define TIE_HEADER

#include "phaser_types.hh"

class TieGradient
// contribution to gradient
{
public:
  TieGradient(){}
  TieGradient(const int& i, const floatType& g) : index(i), Grad(g) {}
  int index;      // index into global parameter list
  floatType Grad; // gradient contribution
};

class TieHessian
// contribution to Hessian
{
  // NB indices here are from 0, but are used from 1
  // in actual Hessian
public:
  TieHessian(){}
  TieHessian(const int& i1, const int& i2, const floatType& H)
    : index1(i1), index2(i2), Hij(H) {}
  int index1;  // Indices for Hessian contribution
  int index2;  //   ... to dR/dp12
  floatType Hij;
};

class Tie
{
public:
  Tie(){}
  // Constructor for restraint to target
  Tie(const int& index, const double& Target, const double& Weight)
    : target(Target), weight(Weight) {kpidx.assign(1, index);}
  // Constructor for tie between 2
  Tie(const int& kindex1, const int& kindex2, const double& Weight)
    : target(0.0), weight(Weight)  
  {kpidx.resize(2); kpidx[0] = kindex1; kpidx[1] = kindex2;}
  // Constructor for tie between 2 or more parameters
  Tie(const std::vector<int>& kindex, const double& Weight)
    : target(0.0), weight(Weight)  {kpidx = kindex;}

  // Returns contribution to restraint residual
  floatType R(const std::vector<double>& params);
  // Returns list of contributions to gradient
  // NB must call R() first
  std::vector<TieGradient> Gradient(const std::vector<double>& params);
  // Returns list of contributions to Hessian
  // NB must call R() first
  std::vector<TieHessian> Hessian(const std::vector<double>& params);

  std::string format() const;

private:
  std::vector<int> kpidx; // list of global parameter indices 
  double target;            // target value if nparam == 1
  double weight;            // weight = 1/sigma^2
  double d;
};

#endif
