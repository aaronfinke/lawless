//
// referencescalemodel.hh
//
// scale to apply to reference intensities to scale to test set
//   k(h) = koverall * exp(-2 B s^2) * A(h)
//     where s = (sin theta/lambda)
//           A(h) = exp(-2pi^2 d*T [Uorth] d*)
//             where [Uorth] is the anisotropic B-factor (no isotropic part)
//           d* is the orthogonalised reciprocal space coordinate [B] h
//
//


#ifndef REFERENCESCALEMODEL_HEADER
#define REFERENCESCALEMODEL_HEADER

#include "anisotropicmodel.hh"
#include "hkl_merge.hh"
#include "mergedlist.hh"

namespace scala {
  //--------------------------------------------------------------
  class ReferenceScaleModel {
  public:

    // types of parameter for refinement (isoB not refined
    enum ScaleParameterType {NONE, SCALE, ANISO};

    ReferenceScaleModel() : status(-99) {}

    ReferenceScaleModel(const MergedList& mergedobslist,
			const hkl_merge& hklmergelist,
			const double& toleranceratio);

    void init(const MergedList& mergedobslist,
	      const hkl_merge& hklmergelist,
	      const double& toleranceratio);

    void init();

    //! return true if OK
    bool isOK() const {return (status == 0);}

    //! format reason for any error
    std::string formatError() const;

    // toleranceratio = 1.0 for difference > maximum resolution,
    //    larger tolerance is more lax
    // set status = -1 if the two lists have different symmetry (point group),
    // or +1 if cell is too different
    // return false if status != 0
    // see also private function
    bool checkCompatible(const hkl_unmerge_list& hkl_list,
			 const hkl_merge& hklmergelist,
			 const double& toleranceratio);

    //! return scale for index hkl
    double scale(const Hkl& hkl) const;

    normalise::AnisotropicModel anisotropicModel() const {return anisomodel;}

    double fderiv(const bool& deriv, const scala::Hkl& hkl,
		  std::vector<double>& dvdp) const;

    // Return restraint target, and optionally gradient & Hessian contributions
    double TieValues(const bool& DoGradient, const bool& DoHessian,
		     const std::vector<double>& params,
		     std::vector<double>& dRdpi,
		     std::vector<TieHessian>& Htie);

    //! number of parameters
    int Nparameters() const {return nparams;}

    //! set parameters
    void SetParameters(const std::vector<double>& params);

    //! final fixup of parameters:
    //   transfer isotropic B from anisomodel to radial component
    void fixUp();

    //! Get vector of parameters
    std::vector<double> GetParameters() const;

    // get lower bound for parameter, depending on type: return false if unbounded
    bool GetLowerBound(const int& Ipar, double& Lower) const;
    // get upper bound for parameter, depending on type: return false if unbounded
    bool GetUpperBound(const int& Ipar, double& Upper) const;
    // get "large shift" value for parameter, depending on type
    double GetLargeShift(const int& Ipar) const;

    std::string format() const;

    std::pair<ScaleParameterType, int>
      GetParameterType(const int& Ipar) const;

    //! return anisotropic part for index hkl
    double anisopart(const Hkl& hkl) const;

  private:
    int status;
    double kscale;
    double kscale0;
    double isoB;   // isotropic Bfactor
    normalise::AnisotropicModel anisomodel;
    int nparams;
    int idx0aniso; // index of 1st aniso parameter
    Scell scell;

    int nrefcommon;

    double celldiff;

    // Ties
    int nties;
    std::vector<Tie> ties;

    void setNparameters();

    // toleranceratio = 1.0 for difference > maximum resolution,
    //    larger tolerance is more lax
    // sets status -1 if the two lists have different symmetry (point group),
    // or +1 if cell is too different
    bool checkCompatible(const MergedList& mergedobslist,
			 const hkl_merge& hklmergelist,
			 const double& toleranceratio);

    void getWilsonScale(const MergedList& mergedobslist,
			const hkl_merge& hklmergelist);

  };
} // namespace Normalise

#endif
