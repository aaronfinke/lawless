//
//  anisotropicmodel.hh
//

#ifndef ANISOTROPICMODEL_HEADER
#define ANISOTROPICMODEL_HEADER

// Clipper
#include <clipper/clipper.h>
//#include "clipper/clipper-contrib.h"
//#include "intensity_scale.h"

#include "anisotropy.hh"  // for class OrthogonalAnisotropy etc
#include "crystaltype.hh"

using namespace scala;

namespace normalise {
  //--------------------------------------------------------------
  class AnisotropicModel {
  public:
    AnisotropicModel() : nparams(0) {}

    AnisotropicModel(const AnisotropicAnalysis& anisotropicanalysis);

    void init(const AnisotropicAnalysis& anisotropicanalysis);

    void init(const clipper::Cell& cCell,
	      const CrystalSystem& crysSys,
	      const clipper::U_aniso_orth& U_aniso_orth);

    void clear();

    //! return [Uorth]
    //   where [Uorth] is the anisotropic B-factor (no isotropic part)
    clipper::U_aniso_orth U_aniso_orth() const {return u_aniso_orth;}

    //! number of variable parameters
    int Nparameters() const {return nparams;}

    //! set parameters
    void SetParameters(const std::vector<double>& params);

    //! Get vector of variable parameters
    std::vector<double> GetParameters() const {return params;}

    // get lower bound for parameter, depending on type: return false if unbounded
    bool GetLowerBound(const int& Ipar, double& Lower) const;
    // get upper bound for parameter, depending on type: return false if unbounded
    bool GetUpperBound(const int& Ipar, double& Upper) const;
    // get "large shift" value for parameter, depending on type
    double GetLargeShift(const int& Ipar) const;

    std::string format() const;

    //! return multiplying scale for index hkl
    //  exp(-2pi^2 d*T [Uorth] d*)
    //   where [Uorth] is the anisotropic B-factor (no isotropic part)
    //   d* is the orthogonalised reciprocal space coordinate [B] h
    double scale(const Hkl& hkl) const;

    // Return scale and derivative vector for index hkl
    double fderiv(const bool& deriv,
		  const scala::Hkl& hkl,
		  std::vector<double>& dkdp) const;

    // Setting
    //! Set u_aniso_orth, umat components and refinable parameters params
    //  remove isotropic part    
    void setUmat(const std::vector<double> pars);

    //! Set u_aniso_orth, umat components and refinable parameters params
    void setUmat(const clipper::U_aniso_orth& U_aniso_orth);

    //! isotropic U^2 subtracted from anisotropic matrix
    double IsotropicPart() const {return isotropicpart/clipper::Util::twopi2();}

    //! modify umat components to make trace = 0.0
    //  store isotropicpart, update umat and params
    void forcePureAnisotropic();

  private:
    clipper::Cell ccell;
    clipper::U_aniso_orth u_aniso_orth;
    CrystalSystem cryssys;

    int nparams;  // number of refinable parameters
    // 6 parameters for [Uorth], not all always refinable
    // in order 00, 11, 22, 01, 02, 12
    std::vector<double> umat;
    std::vector<double> params;  // refinable parameters
    // index into refinable parameters for each umat element
    std::vector<int> lpumat;
    double isotropicpart;

    // set parameters from umat components, allowing for symmetry
    void setParams();

    void U_aniso_orth_from_Umat();

    // set umat components from parameters, allowing for symmetry
    void setUmatFromParams(const std::vector<double> pars);

  };

} // namespace NormaliseNameSpace
#endif
