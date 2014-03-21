// anisotropy.hh


#ifndef ANISO_HEADER
#define ANISO_HEADER

// Clipper
#include <clipper/clipper.h>
#include "clipper/clipper-contrib.h"
#include "intensity_scale.h"

#include "hkl_unmerge.hh"
#include "sdmodel.hh"

typedef clipper::Vec3<double> DVect3;
typedef clipper::Mat33<double> DMat33;

namespace scala {
//----------------------------------------------------------------------------
class OrthogonalAnisotropy {
  //! Generate and store anisotropy tensor from intensity data
public:
  OrthogonalAnisotropy(){}

  //! construct from intensity list
  OrthogonalAnisotropy(clipper::HKL_data<clipper::data32::I_sigI>& isigi);

  //! initialise from intensity list
  void init(clipper::HKL_data<clipper::data32::I_sigI>& isigi);

  clipper::U_aniso_orth u_aniso_orth() const {return uanorth;}

  clipper::U_aniso_frac u_aniso_frac(const clipper::Cell& cell)
  {return sfscl.u_aniso_orth().u_aniso_frac(cell);}

  //! return eigenvalues (orthogonal frame), sorted as closest to a*,b*,c*
  //  U values applied to intensities
  DVect3  EigenValuesOrth() const {return eigval;}

  //! return eigenvectors (orthogonal frame), sorted as closest to a*,b*,c*
  std::vector<DVect3> EigenVectorsOrth() const {return eigvec;}

private:
  clipper::Iscale_aniso<float> sfscl;
  clipper::U_aniso_orth uanorth;  // U tensor, orthogonal frame
  DVect3 eigval;   // eigenvalues
  std::vector<DVect3> eigvec; // eigenvectors (orthogonal frame)

  //! Sort eigenvectors as closest to a*, b*, c*, and eigenvalues
  void SortEigenVectorsOrth();
};
//----------------------------------------------------------------------------
  class AnisotropicAnalysis {
  public:
    AnisotropicAnalysis(){}

    //! construct from unmerged list
    AnisotropicAnalysis(const hkl_unmerge_list& hkl_list,
			const int& datasetindex,
			const SDmodel& SDM);

    //! initialise from unmerged list
    void init(const hkl_unmerge_list& hkl_list,
	      const int& datasetindex,
	      const SDmodel& SDM);

    //! initialise (for testing)
    void init(const hkl_symmetry& ssymmetry,
	      const Scell& cscell);

    //! initialise from intensity list (NB returned scaled)
    void init(const hkl_symmetry& ssymmetry,
	      const Scell& cscell,
	      clipper::HKL_data<clipper::data32::I_sigI>& isigi);

    //! set cone angle in degrees
    void SetConeAngle(const double& angledegrees);
    //! return cone angle in degrees
    double ConeAngle() const;

    // return nearest axis index, if within cone around that axis, -1 if outside, and weight
    std::pair<int,double> Axis(const Hkl& hkl, const Rtype& invresolsq) const;

    // return projections on to principal directions
    // 1) if (!abplane), return 3 projections
    // 2) if (abplane), return 1 projection on to axis if doplane false in element [1]
    //                  return 2 projections on to plane and axis if doplane true
    DVect3 Projection(const Hkl& hkl, const bool& doplane=false) const;

    //! return OrthogonalAnisotropy
    OrthogonalAnisotropy GetOrthogonalAnisotropy() const
    {return orthogonalanisotropy;}

    //! return eigenvectors (orthogonal frame), sorted as closest to a*,b*,c*
    std::vector<DVect3> PrincipalAxes() const {return principalaxes;}

    //! return true if principal axes are general (triclinic, monoclinic)
    bool AreGeneralAxes() const {return lowsymmetry;}

    //! return eigenvalues for B(amplitude) (orthogonal frame), sorted as closest to a*,b*,c*
    DVect3  EigenValuesOrth() const;

    //! return difference between maximum & minimum B-factor
    double BfactorDifference() const;

    //! return number of reflections used in fit, for lowsymmetry (else = 0)
    int NreflUsed() const {return nreflused;}

    //! return true if first direction is plane perpendicular to 3rd direction
    bool IsPlane() const {return abplane;}

    //! return true if cubic
    bool IsCubic() const {return cubic;}

    //! return cell
    // unit cell in clipper convention
    clipper::Cell Cell() const {return ccell;}


    //! format types of analyses done
    std::string formattype() const;

    //! return labels for three axes
    std::vector<std::string> Axesformat() const;

    //! return crystal system
    CrystalSystem crysSys() const {return cryssys;}

  private:
    // principal directions for analysis, along a*, b*, c* for higher symmetry
    // along principal components for monoclinic & triclinic
    // Orthogonal space in Clipper convention, unit vectors
    std::vector<DVect3> principalaxes;
    OrthogonalAnisotropy orthogonalanisotropy;  // anisotropy

    CrystalSystem cryssys;
    bool lowsymmetry;  // true if monoclinic or triclinic, get axes from fit
    int nreflused;     // number of reflections used to fit axes, if lowsymmetry
    bool cubic;   // true if cubic, no analysis
    // true if analysis is relative to ab plane and c* (trigonal, hexagonal, tetragonal)
    bool abplane;
    bool rlattice;  // true for rhombohedral in R-setting

    double coneangle; // radians
    double cosconeangle;
    double sinconeangle;
    hkl_symmetry symmetry;
    clipper::Cell ccell;  // unit cell in clipper convention

    //! return formatted normalised principal axes in hkl space (fractional)
    std::string formatHvector(const clipper::Coord_reci_frac& crdrf) const;
        // Normalise to a maximum value of 1.0
    clipper::Coord_reci_frac NormaliseVector(const clipper::Coord_reci_frac& crdrf) const;

    void SetPrincipalDirectionsGeneral
    (const hkl_unmerge_list& hkl_list,
     const int& datasetindex,
     const SDmodel& SDM);



  };
}
#endif
