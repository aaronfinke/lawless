//
// referencescalemodel.cpp
// scale to apply to reference intensities to scale to test set
//   k(h) = koverall * exp(-2 B s^2) * A(h)
//     where s = (sin theta/lambda)
//           A(h) = exp(-2pi^2 d*T [Uorth] d*)
//             where [Uorth] is the anisotropic B-factor (no isotropic part)
//           d* is the orthogonalised reciprocal space coordinate [B] h
//
//

#include "referencescalemodel.hh"
#include "scala_util.hh"
#include "string_util.hh"
#include "score_datatypes.hh"

namespace scala {
  //--------------------------------------------------------------
  ReferenceScaleModel::ReferenceScaleModel(const MergedList& mergedobslist,
					   const hkl_merge& hklmergelist,
					   const double& toleranceratio)
  {
    init(mergedobslist, hklmergelist, toleranceratio);
  }
  //--------------------------------------------------------------
  void ReferenceScaleModel::init()
  {
    status = 0;
    kscale = 1.0;
    kscale0 = 1.0;
    isoB = 0.0;   // isotropic Bfactor
    anisomodel.clear();
  }
  //--------------------------------------------------------------
  void ReferenceScaleModel::init(const MergedList& mergedobslist,
				 const hkl_merge& hklmergelist,
				 const double& toleranceratio)
  //  initialise from unmerged list and a merged reference list
  {
    status = 0;
    checkCompatible(mergedobslist, hklmergelist, toleranceratio); // sets status
    if (status != 0) {
      return;
    }
    
    scell = Scell(mergedobslist.Cell());
    nties = 0;
    
    // isotropic scaling on mean intensities
    // sets status = +2 if failed
    getWilsonScale(mergedobslist, hklmergelist);
    
    anisomodel.init(scell.ClipperCell(),
		    hkl_symmetry(mergedobslist.spacegroup()).CrysSys(),
		    clipper::U_aniso_orth(0.0));
		    
    setNparameters();
  }
  //--------------------------------------------------------------
  bool ReferenceScaleModel::checkCompatible(const MergedList& mergedobslist,
					    const hkl_merge& hklmergelist,
					    const double& toleranceratio)
  // toleranceratio = 1.0 for difference > maximum resolution,
  //    larger tolerance is more lax
  // set status = -1 if the two lists have different symmetry (point group),
  // or +1 if cell is too different
  // return false if status != 0
  {
    SpaceGroup SGtest = SpaceGroup(mergedobslist.spacegroup()).PointGroup();
    SpaceGroup  SGref = SpaceGroup(hklmergelist.spacegroup()).PointGroup();
    if (SGtest != SGref) {
      status = -1;
    }
    clipper::Cell celltest = mergedobslist.Cell();
    clipper::Cell cellref  = hklmergelist.Cell().ClipperCell();
    // minimum high resolution  of the two datasets
    double minreshigh = std::max(mergedobslist.resHigh(),
				 hklmergelist.resHigh());
    // "difference" in A
    if (!celltest.equals(cellref, toleranceratio*minreshigh)) {
      status = +1;
    }
    return (status == 0);
  }
  //--------------------------------------------------------------
  bool ReferenceScaleModel::checkCompatible(const hkl_unmerge_list& hkl_list,
					    const hkl_merge& hklmergelist,
					    const double& toleranceratio)
  // toleranceratio = 1.0 for difference > maximum resolution,
  //    larger tolerance is more lax
  // set status = -1 if the two lists have different symmetry (point group),
  // or +1 if cell is too different
  // return false if status != 0
  {
    status = 0;
    SpaceGroup SGtest = hkl_list.symmetry().GetSpaceGroup().PointGroup();
    SpaceGroup  SGref = SpaceGroup(hklmergelist.spacegroup()).PointGroup();
    if (SGtest != SGref) {
      status = -1;
    }
    clipper::Cell celltest = hkl_list.Cell().ClipperCell();
    clipper::Cell cellref  = hklmergelist.Cell().ClipperCell();
    // minimum high resolution  of the two datasets
    double minreshigh = std::max(hkl_list.ResRange().ResHigh(),
				 hklmergelist.resHigh());
    // "difference" in A
    if (!celltest.equals(cellref, toleranceratio*minreshigh)) {
      status = +1;
    }
    return (status == 0);
  }
  //--------------------------------------------------------------
  std::string ReferenceScaleModel::format() const
  {
    std::string s;
    s =
      "Scale = k exp(-2B (sin theta/lambda)^2) * Aniso(h):\n";
    s += "  k = "+
      StringUtil::ftos(kscale, 9, 4);
    s += "  relative Wilson B-factor: "+
      StringUtil::ftos(isoB, 6, 1)+"\n";
    s += anisomodel.format();
    return s;
  }
  //--------------------------------------------------------------
  //! format reason for any error
  std::string ReferenceScaleModel::formatError() const
  {
    std::string s = "";
    if (status == 0) {return s;}
    if (status == -99) {
      s = "Not initialised";
    } else if (status == -1) {
      s = "Reference and test lists have different point groups";
    } else if (status == +1) {
      s = "Reference and test lists have different unit cells, difference "+
	StringUtil::ftos(celldiff,7,3) + "A";
    } else if (status == +2) {
      s = "Wilson scaling between reference and test lists failed";
    } else {
      s = "Unknown error";
    }
    return s;
  }
  //--------------------------------------------------------------
  void ReferenceScaleModel::getWilsonScale(const MergedList& mergedobslist,
					   const hkl_merge& hklmergelist)
  // sets status = +2 if failed
  {
    // minimum high resolution  of the two datasets
    double minreshigh = std::max(mergedobslist.resHigh(),
				 hklmergelist.resHigh());
    ResoRange resrange(100000.0, minreshigh, mergedobslist.num_reflections());
    int nbins = resrange.Nbins();

    std::vector<MeanValue> Imeantest(nbins);
    std::vector<MeanValue> Imeanref(nbins);

    // loop all test reflections
    nrefcommon = 0;
  
    typedef clipper::HKL_data_base::HKL_reference_index HRI;
    // Observed data
    mergedobslist.start(0); // dataset 0, should be only one

    IsigI Iobs;
    while (mergedobslist.next(Iobs)) {
      clipper::HKL hkl = mergedobslist.hkl().HKL();
      // find equivalent in reference list, if present
      IsigI Isref = hklmergelist.Isig(hkl);
      if (Isref.sigI() > 0.0) {
	ftype  I = Iobs.I();     // Iobs
	ftype sd = Iobs.sigI();
	if (sd > 0.0) {

	  ftype invresolsq = mergedobslist.invresolsq();
	  // Resolution bin
	  int mres = resrange.bin(invresolsq);
	  nrefcommon++;
	  Imeantest[mres].Add(I);
	  Imeanref[mres].Add(Isref.I());
	}
      }
    }

    // Wilson scale
    // <Itest>[j] ~= k exp(-2 B sj^2) <Iref>[j] for j'th resolution bin
    // ie ln(<Itest>[j]/<Iref>[j]) ~= ln(k) - 2 B sj^2
    //  plot vs sj^2, slope = -2 B, intercept = ln(k)

    LinearFit line;

    for (int mres=0;mres<nbins;++mres) {
      if (Imeanref[mres].Mean() != 0.0) {
	double ratio = Imeantest[mres].Mean()/Imeanref[mres].Mean();
	double s2 = 0.25*resrange.middle(mres); // 1/(4d^2)
	double w = 1.0;
	// downweight 1st & last bins
	if (mres == 0 || mres == nbins-1) {w = 0.5;}
	line.add(s2, log(ratio), w);
      }
    }
    if (line.Number() > 0) {
      RPair result = line.result();
      kscale = exp(result.second);
      kscale0 = kscale;
      isoB = - 0.5 * result.first;
    } else {
      status = +2;
    }
  }
  //--------------------------------------------------------------
  void ReferenceScaleModel::setNparameters()
  {
    nparams = 1;
    idx0aniso = nparams;
    nparams += anisomodel.Nparameters();
  }
  //--------------------------------------------------------------
  //! return scale for index hkl
  double ReferenceScaleModel::scale(const Hkl& hkl) const
  {
    double radial = kscale * exp (-0.5 * isoB * hkl.invresolsq(scell));
    double aniso  = anisomodel.scale(hkl);
    return radial * aniso;
  }
  //--------------------------------------------------------------
  //! return anisotropic part for index hkl
  double ReferenceScaleModel::anisopart(const Hkl& hkl) const
  {
    return anisomodel.scale(hkl);
  }
  //--------------------------------------------------------------
  double ReferenceScaleModel::fderiv(const bool& deriv, const scala::Hkl& hkl,
			     std::vector<double>& dvdp) const
  {
    dvdp.resize(nparams);
    std::vector<double> dadp; // anisotropy

    double biso = exp (-0.5 * isoB * hkl.invresolsq(scell));
    double radial = kscale * biso;
    double drdp = biso;  // scale
    
    double aniso = anisomodel.fderiv(true, hkl, dadp);

    // radial dvdp = dradial/dp * aniso
    dvdp[0] = drdp * aniso;

    // aniso dvdp = radial * daniso/dp
    for (size_t i=0; i<dadp.size(); i++) { 
      dvdp[idx0aniso+i] = radial * dadp[i];
    }

    return radial * aniso;
  }
  //--------------------------------------------------------------
  // Return restraint target, and optionally gradient & Hessian contributions
  double ReferenceScaleModel::TieValues(const bool& DoGradient,
				const bool& DoHessian,
				const std::vector<double>& params,
				std::vector<double>& dRdpi,
				std::vector<TieHessian>& Htie)
  // Return restraint target, and optionally gradient & Hessian contributions
  //
  // On entry:
  //  DoGradient  true to calculate gradient dRdpi
  //  DoHessian   true to calculate Hessian terms Htie
  //
  // On exit:
  //  dRdpi(nparameters)  gradient contribution for each parameter
  //  Htie                list of indexed Hessian contributions
  //  
  {
    if (DoGradient) {
      dRdpi.assign(nparams, 0.0); // clear derivatives
      if (DoHessian) {
	Htie.clear();                   // and Hessian
      }}

    double R = 0.0;
    for (int itie=0;itie<nties;++itie) { // loop ties
      R += ties[itie].R(params);
      if (DoGradient) {
	std::vector<TieGradient> grad = ties[itie].Gradient(params);
	for (size_t i=0;i<grad.size();++i) {
	  dRdpi[grad[i].index] = grad[i].Grad;
	}
	if (DoHessian) {
	  std::vector<TieHessian> hessian = ties[itie].Hessian(params);
	  Htie.insert(Htie.end(), hessian.begin(), hessian.end());
	}
      }
    }
    return R;
  }
  //--------------------------------------------------------------
  //! set parameters
  void ReferenceScaleModel::SetParameters(const std::vector<double>& params)
  {
    ASSERT (int(params.size()) == nparams);
    std::vector<double>::const_iterator pos1 = params.begin();  // start of range
    std::vector<double>::const_iterator pos2;                   // end of range
    // Radial model, 1 parameter
    pos2 = pos1 + 1;
    kscale = params[0];
    pos1 = pos2;
    // Anisotropic model
    pos2 = pos1 + anisomodel.Nparameters();
    anisomodel.SetParameters(std::vector<double>(pos1, pos2));
  }
  //--------------------------------------------------------------
  void ReferenceScaleModel::fixUp()
  // final fixup of parameters:
  //   transfer isotropic B from anisomodel to radial model
  {
    // Apply isotropic part removed from anisotropic B-factor
    anisomodel.forcePureAnisotropic();
    isoB += anisomodel.IsotropicPart();
  }
  //--------------------------------------------------------------
  //! Get vector of parameters
  std::vector<double> ReferenceScaleModel::GetParameters() const
  {
    std::vector<double> params;
    std::vector<double> p;
    params.push_back(kscale);

    p = anisomodel.GetParameters();
    params.insert(params.end(), p.begin(), p.end());
    return params;
  }
  //--------------------------------------------------------------
  std::pair<ReferenceScaleModel::ScaleParameterType, int>
     ReferenceScaleModel::GetParameterType(const int& Ipar) const
  {
    ASSERT (Ipar < nparams);
    ScaleParameterType type;
    int idx = -1;
    if (Ipar == 0) {
      type = SCALE;
      idx = Ipar;
    } else {
      type = ANISO;
      idx = Ipar - 1;
    }
    ASSERT (idx >= 0);
    return std::pair<ScaleParameterType, int>(type, idx);
  }
  //--------------------------------------------------------------
  // get lower bound for parameter, depending on type: return false if unbounded
  bool ReferenceScaleModel::GetLowerBound(const int& Ipar, double& Lower) const
  {
    std::pair<ScaleParameterType, int>  partypeidx = GetParameterType(Ipar);
    switch (partypeidx.first) {
    case ReferenceScaleModel::SCALE:
      return 0.001;
    case ReferenceScaleModel::ANISO:
      return anisomodel.GetLowerBound(partypeidx.second, Lower);
    }
    return false;
  }
  //--------------------------------------------------------------
  // get upper bound for parameter, depending on type: return false if unbounded
  bool ReferenceScaleModel::GetUpperBound(const int& Ipar, double& Upper) const
  {
    std::pair<ScaleParameterType, int>  partypeidx = GetParameterType(Ipar);
    switch (partypeidx.first) {
    case ReferenceScaleModel::SCALE:
      return 1000.;
    case ReferenceScaleModel::ANISO:
      return anisomodel.GetUpperBound(partypeidx.second, Upper);
    }
    return false;
  }
  //--------------------------------------------------------------
  // get large shift for parameter, depending on type
  double ReferenceScaleModel::GetLargeShift(const int& Ipar) const
  {
    std::pair<ScaleParameterType, int>  partypeidx = GetParameterType(Ipar);
    switch (partypeidx.first) {
    case ReferenceScaleModel::SCALE:
      return 0.5*kscale0; // use initial value to get large shift
    case ReferenceScaleModel::ANISO:
      return anisomodel.GetLargeShift(partypeidx.second);
    }
    return 0.0;
  }
  //--------------------------------------------------------------


} // namespace scala
