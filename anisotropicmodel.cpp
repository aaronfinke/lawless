//
//  anisotropicmodel.cpp
//


#include "anisotropicmodel.hh"
#include "string_util.hh"

namespace normalise {
  //--------------------------------------------------------------
   AnisotropicModel::AnisotropicModel
   (const AnisotropicAnalysis& anisotropicanalysis)
   {
     init(anisotropicanalysis);
   }
  //--------------------------------------------------------------
  void AnisotropicModel::init(const AnisotropicAnalysis& anisotropicanalysis)
  {
    OrthogonalAnisotropy orthogonalanisotropy =
      anisotropicanalysis.GetOrthogonalAnisotropy();
    ccell = anisotropicanalysis.Cell();
    cryssys = anisotropicanalysis.crysSys();
    u_aniso_orth = orthogonalanisotropy.u_aniso_orth();
    isotropicpart = 0.0;
    setUmat(u_aniso_orth);
    setParams();  // from u_aniso_orth
  }
  //--------------------------------------------------------------
  void AnisotropicModel::init(const clipper::Cell& cCell,
			      const CrystalSystem& crysSys,
			      const clipper::U_aniso_orth& U_aniso_orth)
  {
    ccell = cCell;
    cryssys = crysSys;
    u_aniso_orth = U_aniso_orth;
    isotropicpart = 0.0;
    setUmat(u_aniso_orth);
    setParams();  // from u_aniso_orth
  }
  //--------------------------------------------------------------
  void AnisotropicModel::clear()
  {
    u_aniso_orth = clipper::U_aniso_orth(0.0);
    isotropicpart = 0.0;
  }
  //--------------------------------------------------------------
  void AnisotropicModel::setUmat(const std::vector<double> pars)
  // set umat components from parameters, allowing for symmetry
  {
    ASSERT (int(pars.size()) == nparams);
    setUmatFromParams(pars);
    setParams();       // set parameters from umat
    U_aniso_orth_from_Umat();
  }
  //--------------------------------------------------------------
  void AnisotropicModel::setUmat(const clipper::U_aniso_orth& U_aniso_orth)
  // Store u_aniso_orth, umat components and refinable parameters params
  // note that umat components are 2 pi^2 * U components
  //   cf clipper/core/resol_basisfn.cpp::BasisFn_aniso_gaussian
  //     (NB not quite the same, 6 parameters not 7!)
  {
    umat.assign(6, 0.0);
    umat[0] = U_aniso_orth.mat00() * clipper::Util::twopi2();
    umat[1] = U_aniso_orth.mat11() * clipper::Util::twopi2();
    umat[2] = U_aniso_orth.mat22() * clipper::Util::twopi2();
    umat[3] = U_aniso_orth.mat01() * clipper::Util::twopi2();
    umat[4] = U_aniso_orth.mat02() * clipper::Util::twopi2();
    umat[5] = U_aniso_orth.mat12() * clipper::Util::twopi2();
    setParams();       // reset parameters from umat
  }
  //--------------------------------------------------------------
  //! set parameters
  void AnisotropicModel::SetParameters(const std::vector<double>& params)
  {
    if (params.size() > 0) {
      setUmat(params);
    }
  }
  //--------------------------------------------------------------
  void AnisotropicModel::setParams()
  // set parameters from umat components, allowing for symmetry
  // set indices lpumat to point from umat components to params
  {
    lpumat.resize(6);
    for (int j=0;j<6;++j) {lpumat[j] = j;}
    if (cryssys == TRICLINIC) {
      nparams = 6;
      params = umat;  // just copy all
    } else if (cryssys == MONOCLINIC) {
      nparams = 4;
      params.resize(nparams);
      params[0] = umat[0];  // U00
      params[1] = umat[1];
      params[2] = umat[2];
      params[3] = umat[4];  // U02  ie a c
      lpumat[3] = -1;
      lpumat[4] = 3;
      lpumat[5] = -1;
    } else if (cryssys == ORTHORHOMBIC) {
      nparams = 3;
      params.resize(nparams);
      params[0] = umat[0];  // U00
      params[1] = umat[1];
      params[2] = umat[2];
      lpumat[3] = -1;
      lpumat[4] = -1;
      lpumat[5] = -1;
    } else if (cryssys == TETRAGONAL) {
      nparams = 2;
      params.resize(nparams);
      params[0] = umat[0];  // U00 = U11
      params[1] = umat[2];  // U22
      lpumat[0] = 0;
      lpumat[1] = 0;
      lpumat[2] = 1;
      lpumat[3] = -1;
      lpumat[4] = -1;
      lpumat[5] = -1;
    } else if (cryssys == TRIGONAL &&
	       scala::RhombohedralAxes(Scell(ccell).UnitCell())) {
      nparams = 2;
      params.resize(nparams);
      params[0] = umat[0];  // U00 = U11 = U22
      params[1] = umat[3];  // U01 = U02 = U12
      lpumat[0] = 0;
      lpumat[1] = 0;
      lpumat[2] = 0;
      lpumat[3] = 1;
      lpumat[4] = 1;
      lpumat[5] = 1;
    } else if (cryssys == TRIGONAL || cryssys == HEXAGONAL) {
      nparams = 3;
      params.resize(nparams);
      params[0] = umat[0];  // U00 = U11
      params[1] = umat[2];  // U22
      params[2] = umat[3];  // U01
      lpumat[0] = 0;
      lpumat[1] = 0;
      lpumat[2] = 1;
      lpumat[3] = 2;
      lpumat[4] = -1;
      lpumat[5] = -1;
    } else if (cryssys == CUBIC) {
      nparams = 0;
      params.clear();
      lpumat.assign(6,-1);
    } else { // shouldn't happen
      Message::message(Message_fatal
		       ("AnisotropicModel: undefined Bravais lattice\n"));
    }
  }
  //--------------------------------------------------------------
  void AnisotropicModel::setUmatFromParams(const std::vector<double> pars)
  // set umat components from parameters, allowing for symmetry
  {
    umat.assign(6, 0.0);
    if (cryssys == TRICLINIC) {
      umat = pars;  // just copy all
    } else if (cryssys == MONOCLINIC) {
      umat[0] = pars[0];  // U00
      umat[1] = pars[1];
      umat[2] = pars[2];
      umat[4] = pars[3];  // U02  ie ac
    } else if (cryssys == ORTHORHOMBIC) {
      umat[0] = pars[0];  // U00
      umat[1] = pars[1];
      umat[2] = pars[2];
    } else if (cryssys == TETRAGONAL) {
      umat[0] = pars[0];  // U00 = U1
      umat[1] = pars[0];  // U00 = U1
      umat[2] = pars[1];  // U22
    } else if (cryssys == TRIGONAL &&
	       scala::RhombohedralAxes(Scell(ccell).UnitCell())) {
      umat[0] = pars[0];  // U00 = U11 = U22
      umat[1] = pars[0];  // U00 = U11 = U22
      umat[2] = pars[0];  // U00 = U11 = U22
      umat[3] = pars[1];  // U01 = U02 = U12
      umat[4] = pars[1];  // U01 = U02 = U12
      umat[5] = pars[1];  // U01 = U02 = U12
    } else if (cryssys == TRIGONAL || cryssys == HEXAGONAL) {
      umat[0] = pars[0];  // U00 = U11
      umat[1] = pars[0];  // U00 = U11
      umat[2] = pars[1];  // U00 = U11
      umat[3] = pars[2];  // U01
    } else if (cryssys == CUBIC) {
    } else { // shouldn't happen
      Message::message(Message_fatal
		       ("AnisotropicModel: undefined Bravais lattice\n"));
    }
  }
  //--------------------------------------------------------------
  void AnisotropicModel::U_aniso_orth_from_Umat()
  // Set u_aniso_orth from umat components
  //   cf clipper/core/resol_basisfn.cpp::BasisFn_aniso_gaussian
  //     (NB not quite the same!)
  {
    u_aniso_orth = clipper::U_aniso_orth
      (umat[0]/clipper::Util::twopi2(), umat[1]/clipper::Util::twopi2(),
       umat[2]/clipper::Util::twopi2(), umat[3]/clipper::Util::twopi2(),
       umat[4]/clipper::Util::twopi2(), umat[5]/clipper::Util::twopi2());
  }
  //--------------------------------------------------------------
  void AnisotropicModel::forcePureAnisotropic()
  // modify umat components to make trace = 0.0
  {
    if (nparams > 0) {
      isotropicpart = (umat[0] + umat[1] + umat[2])/3.0;
      umat[0] -= isotropicpart;
      umat[1] -= isotropicpart;
      umat[2] -= isotropicpart;
      U_aniso_orth_from_Umat();
      setParams();  // reset params
    }
  }
  //--------------------------------------------------------------
  //! return multiplying scale for index hkl
  double AnisotropicModel::scale(const scala::Hkl& hkl) const
  {
    if (nparams <= 0) {return 1.0;} // null model
    std::vector<double> dkdp;
    double scale = fderiv(false, hkl, dkdp);
    return scale;
  }
  //--------------------------------------------------------------
  // Return scale and derivative vector for index hkl
  double AnisotropicModel::fderiv(const bool& deriv,
				  const scala::Hkl& hkl,
				  std::vector<double>& dkdp) const
  //  exp(-2pi^2 d*T [Uorth] d*)
  //   where [Uorth] is the anisotropic B-factor (no isotropic part)
  //   d* is the orthogonalised reciprocal space coordinate [B] h
  // umat components are already multiplied by +2pi^2
  {
    double scale = 1.0;
    if (nparams <= 0) { // null model
      dkdp.clear();
      return scale;
    }

    clipper::HKL HKL = hkl.HKL();
    clipper::Coord_reci_orth xs =
      clipper::Coord_reci_frac(HKL).coord_reci_orth(ccell);
    double c[6];
    c[0] = -xs[0]*xs[0];	  
    c[1] = -xs[1]*xs[1];	  
    c[2] = -xs[2]*xs[2];	  
    c[3] = -2.0*xs[0]*xs[1];
    c[4] = -2.0*xs[0]*xs[2];
    c[5] = -2.0*xs[1]*xs[2];
    scale = exp(c[0]*umat[0] + c[1]*umat[1] + c[2]*umat[2] +
		c[3]*umat[3] + c[4]*umat[4] + c[5]*umat[5]);

    if (deriv) {
      // dscale/dp for active parameters
      dkdp.assign(nparams, 0.0);
      // for each of the 6 umat[j] components, lpumat[j] indexes
      //  the active parameter, or = -1 if fixed
      for (int j=0;j<6;++j) {
	int k = lpumat[j];
	if (k >= 0) {
	  dkdp[k] += c[j];
	}
      }
      for (int k=0;k<nparams;++k) {
	dkdp[k] *= scale;
      }
    }
    return scale;
  }
  //--------------------------------------------------------------
  // get lower bound for parameter, depending on type: return false if unbounded
  bool AnisotropicModel::GetLowerBound(const int& Ipar, double& Lower) const
  {
    return false;
  }
  //--------------------------------------------------------------
  // get upper bound for parameter, depending on type: return false if unbounded
  bool AnisotropicModel::GetUpperBound(const int& Ipar, double& Upper) const
  {
      return false;
  }
  //--------------------------------------------------------------
  // get "large shift" value for parameter, depending on type
  double AnisotropicModel::GetLargeShift(const int& Ipar) const
  {
    return 0.1;
  }
  //--------------------------------------------------------------
  std::string AnisotropicModel::format() const
  {
    std::string s = std::string("Anisotropic scale model:")+
      " a = exp ( xT [Uo] x)\n"+
      "where x = [B] h, and [Uo] are on orthogonal axes\n";
    s += std::string("            (")+
      " "+StringUtil::ftos(u_aniso_orth.mat00(), 9, 5)+
      " "+StringUtil::ftos(u_aniso_orth.mat01(), 9, 5)+
      " "+StringUtil::ftos(u_aniso_orth.mat02(), 9, 5)+" )\n";
    s += std::string("   [Uo]  =  (")+
      " "+StringUtil::ftos(u_aniso_orth.mat01(), 9, 5)+
      " "+StringUtil::ftos(u_aniso_orth.mat11(), 9, 5)+
      " "+StringUtil::ftos(u_aniso_orth.mat12(), 9, 5)+" )\n";
    s += std::string("            (")+
      " "+StringUtil::ftos(u_aniso_orth.mat02(), 9, 5)+
      " "+StringUtil::ftos(u_aniso_orth.mat12(), 9, 5)+
      " "+StringUtil::ftos(u_aniso_orth.mat22(), 9, 5)+" )\n";
    return s;
  }
  //--------------------------------------------------------------

}
