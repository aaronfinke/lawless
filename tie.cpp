//
//  tie.cpp
//

#include "tie.hh"
#include "string_util.hh"

//--------------------------------------------------------------
floatType Tie::R(const std::vector<double>& params)
// Returns contribution to restraint residual
{
  if (kpidx.size() == 1) {
    d = params[kpidx[0]] - target;
    return weight * d * d;
  } else if (kpidx.size() == 2) {
    d = params[kpidx[0]] - params[kpidx[1]];
    return weight * d * d;
  } else {
    d = 0.0;
    for (size_t i=0;i<kpidx.size();++i) {
      d += params[kpidx[i]];}
    d /= double(kpidx.size());  // average parameter
    double sd = 0.0;
    double s;
    for (size_t i=0;i<kpidx.size();++i) {
      s = (params[kpidx[i]] - d);
      sd += s*s;
    }
    return weight * sd;
  }
}
//--------------------------------------------------------------
std::vector<TieGradient> Tie::Gradient(const std::vector<double>& params)
// Returns list of contributions to gradient
// Assumes d has been calculated before by call to R()
{
  std::vector<TieGradient> TG;
  if (kpidx.size() == 1) {
    // d is p(i) - target
    TG.push_back(TieGradient(kpidx[0], weight*d));
  } else if (kpidx.size() == 2) {
    // d is p(i) - p(j)
    TG.push_back(TieGradient(kpidx[0], +weight*d));
    TG.push_back(TieGradient(kpidx[1], -weight*d));
  } else {
    // d is <p>
    for (size_t i=0;i<kpidx.size();++i) {
      TG.push_back(TieGradient(kpidx[i],
               weight*(params[kpidx[i]] - d)/double(kpidx.size())));
    }
  }
    return TG;
}
//--------------------------------------------------------------
std::vector<TieHessian> Tie::Hessian(const std::vector<double>& params)
// Returns list of contributions to Hessian
{
  std::vector<TieHessian> TH;
  if (kpidx.size() == 1) {
    TH.push_back(TieHessian(kpidx[0], kpidx[0], weight));
  } else if (kpidx.size() == 2) {
    TH.push_back(TieHessian(kpidx[0], kpidx[0], weight));
    TH.push_back(TieHessian(kpidx[1], kpidx[1], weight));
    TH.push_back(TieHessian(kpidx[0], kpidx[1], -weight));
  } else {
    double wn = 1./double(kpidx.size()*kpidx.size());
    for (size_t i=0;i<kpidx.size();++i) {
      for (size_t j=i;j<kpidx.size();++j) {
        TH.push_back(TieHessian(kpidx[i], kpidx[j], weight*wn));
      }}
  }
  return TH;
}
//--------------------------------------------------------------
std::string Tie::format() const
{
  std::string s;
  if (kpidx.size() == 1) {
    // restraint to target
    s += "Tie index "+StringUtil::itos(kpidx[0],5)+
      " to target "+StringUtil::ftos(target,8,4)+
      " weight "+StringUtil::ftos(weight,8,4);
  } else if (kpidx.size() == 2) {
    // tie between two
    s += "Tie index "+StringUtil::itos(kpidx[0],5)+
      " to index "+StringUtil::itos(kpidx[1],5)+
      " weight "+StringUtil::ftos(weight,8,4);
  } else if (kpidx.size() > 2) {
    // tie between > two
    s += "Tie togetherindices ";
    for (size_t i=0; i<kpidx.size(); i++) {
      s += StringUtil::itos(kpidx[i],5);
      if (i<kpidx.size()-1) {s += ", ";}
    }
    " weight "+StringUtil::ftos(weight,8,4);
  } else {
    s = "No tie";
  }
  return s;
}
