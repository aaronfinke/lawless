// aniso.cpp

#include "anisotropy.hh"
#include "mergedlist.hh"
#include "string_util.hh"

using clipper::Message;
using clipper::Message_fatal;



namespace scala {
//-------------------------------------------------------------------------------
  OrthogonalAnisotropy::OrthogonalAnisotropy
  (clipper::HKL_data<clipper::data32::I_sigI>& isigi)
  //! construct from intensity list
{
  init(isigi);
}
//-------------------------------------------------------------------------------
void OrthogonalAnisotropy::init(clipper::HKL_data<clipper::data32::I_sigI>& isigi)
  //! initialise from intensity list
{
  sfscl = clipper::Iscale_aniso<float>(3.0);
  sfscl(isigi);

  uanorth = sfscl.u_aniso_orth();
  clipper::Matrix<double> Uorth(3,3); 
  for (int j=0;j<3;++j) {
    for (int i=0;i<3;++i) {
      Uorth(i,j) = uanorth(i,j);
    }}
  eigvec.resize(3);
  std::vector<double> eo = Uorth.eigen();
  for (int j=0;j<3;++j) {
    for (int i=0;i<3;++i) {
      eigvec[j][i] = Uorth(i,j);  // store eigenvectors (orthogonal frame)
    }
    eigval[j] = eo[j];  // store eigenvalues
  }
  // sort on closest to a*b*c*
  SortEigenVectorsOrth();
}
//--------------------------------------------------------------------------
void OrthogonalAnisotropy::SortEigenVectorsOrth()
//! Sort eigenvectorsas closest to a*, b*, c*, and eigenvalues
{
  std::vector<DVect3> eigvecsrt(3);
  DVect3 eigvalsrt;
  std::vector<int> close(3, -1);
  
  for (int i=0;i<3;++i) { // loop a*, b*, c*
    double evcmax = 0.0;
    for (int j=0;j<3;++j) { // loop vectors
      if (std::abs(eigvec[j][i]) > evcmax) {
	evcmax = std::abs(eigvec[j][i]);
	close[i] = j;
      }
    }
  }
  // close[i] for i -> a*, b*, c* is closest vector
  for (int i=0;i<3;++i) { // loop a*, b*, c*
    if (close[i] < 0) { // fail
      Message::message(Message_fatal("OrthogonalAnisotropy::SortEigenVectorsOrth fail\n"));
    }
    eigvecsrt[i] = eigvec[close[i]];
    eigvalsrt[i] = eigval[close[i]];
  }
  // replace
  for (int i=0;i<3;++i) {
    eigvec[i] = eigvecsrt[i];
  }
  eigval = eigvalsrt;
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
  //! construct from unmerged list
  AnisotropicAnalysis::AnisotropicAnalysis(const hkl_unmerge_list& hkl_list,
					   const int& datasetindex,
					   const SDmodel& SDM)
  {
    init(hkl_list,datasetindex,SDM);
  }
//--------------------------------------------------------------------------
  //! initialise from unmerged list
  void AnisotropicAnalysis::init(const hkl_unmerge_list& hkl_list,
				 const int& datasetindex,
				 const SDmodel& SDM)
  {
    // Anisotropy analysis is based on 3 directions, sometimes 2
    // Cases
    //  1) triclinic, monoclinic: directions come from the eigenvectors of
    //     the orthogonal anisotropic Uorth tensor. For monoclinc, the 2nd one
    //     is along b* (k)
    //  2) orthorhombic: directions are along a*, b*, c*
    //  3) trigonal, hexagonal, tetragonal, rhombohedral (H setting):
    //     one direction is along c*, the other "direction" is the a* b* plane
    //  4) cubic: all directions are equivalent, no anisotropy

    symmetry = hkl_list.symmetry(); // store symmetry
    cryssys = symmetry.CrysSys();
    ccell = hkl_list.Cell().ClipperCell();

    // Defaults
    lowsymmetry = false;
    cubic = false;
    abplane = false;
    rlattice = false;

    // case (1)
    if (cryssys == TRICLINIC || cryssys == MONOCLINIC) {
      // Low symmetry, get principal axes from anisotropic U tensor
      lowsymmetry = true;
      SetPrincipalDirectionsGeneral(hkl_list, datasetindex, SDM);
    // case (2)
    } else if (cryssys == ORTHORHOMBIC){  // orthorhombic, just set to a*, b*, c*
      principalaxes.assign(3, DVect3(0.0,0.0,0.0));
      for (int i=0;i<3;++i) {
	principalaxes[i][i] = 1.0;
      }
    // case (3)
    } else if (cryssys == TETRAGONAL || cryssys == TRIGONAL
	       || cryssys == HEXAGONAL) {
      abplane = true;  // analysis against plane perpendicular to c*
      principalaxes.assign(3, DVect3(0.0,0.0,1.0));  // all along c*

      // Is it rhombohedral in R setting?
      if (cryssys == TRIGONAL && RhombohedralAxes(hkl_list.Cell().UnitCell())) {
	rlattice = true;
	DVect3 diagonal = clipper::Coord_reci_frac(1.0,1.0,1.0).coord_reci_orth(ccell);
	principalaxes.assign(3, diagonal.unit()); // all along diagonal, unit vector
      }
      // case (4)
    } else if (cryssys == CUBIC) {
      cubic = true;
      principalaxes.assign(3, DVect3(0.0,0.0,0.0));  // dummy
    } else { // shouldn't happen
      Message::message(Message_fatal
		       ("AnisotropicAnalysis: undefined Bravais lattice\n"));
    }
  }
//--------------------------------------------------------------------------
  //! initialise (for testing)
void AnisotropicAnalysis::init(const hkl_symmetry& ssymmetry,
			       const Scell& cscell)
  {
    // Anisotropy analysis is based on 3 directions, sometimes 2
    // Cases
    //  1) triclinic, monoclinic: directions come from the eigenvectors of
    //     the orthogonal anisotropic Uorth tensor. For monoclinc, the 2nd one
    //     is along b* (k)
    //  2) orthorhombic: driections are along a*, b*, c*
    //  3) trigonal, hexagonal, tetragonal, rhombohedral (H setting):
    //     one direction is along c*, the other "direction" is the a* b* plane
    //  4) cubic: all directions are equivalent, no anisotropy

    symmetry = ssymmetry; // store symmetry
    cryssys = symmetry.CrysSys();
    ccell = cscell.ClipperCell();

    // Defaults
    lowsymmetry = false;
    cubic = false;
    abplane = false;
    rlattice = false;

    // case (1), (2)  (not for real use)
    if (cryssys == TRICLINIC || cryssys == MONOCLINIC ||
	cryssys == ORTHORHOMBIC){  // just set to a*, b*, c*
      principalaxes.assign(3, DVect3(0.0,0.0,0.0));
      for (int i=0;i<3;++i) {
	principalaxes[i][i] = 1.0;
      }
    // case (3)
    } else if (cryssys == TETRAGONAL || cryssys == TRIGONAL
	       || cryssys == HEXAGONAL) {
      abplane = true;  // analysis against plane perpendicular to c*
      principalaxes.assign(3, DVect3(0.0,0.0,1.0));  // all along c*

      // Is it rhombohedral in R setting?
      if (cryssys == TRIGONAL && RhombohedralAxes(cscell.UnitCell())) {
	rlattice = true;
	DVect3 diagonal = clipper::Coord_reci_frac(1.0,1.0,1.0).coord_reci_orth(ccell);
	principalaxes.assign(3, diagonal.unit()); // all along diagonal, unit vector
      }
      // case (4)
    } else if (cryssys == CUBIC) {
      cubic = true;
      principalaxes.assign(3, DVect3(0.0,0.0,0.0));  // dummy
    } else { // shouldn't happen
      Message::message(Message_fatal
		       ("AnisotropicAnalysis: undefined Bravais lattice\n"));
    }
  }
//--------------------------------------------------------------------------
//! initialise from intensity list
void AnisotropicAnalysis::init(const hkl_symmetry& ssymmetry,
			       const Scell& cscell,
			       clipper::HKL_data<clipper::data32::I_sigI>& isigi)
{
  init(ssymmetry, cscell);  // initialise symmetry etc
  if (cryssys == TRICLINIC || cryssys == MONOCLINIC) {
    // Low symmetry, get principal axes from anisotropic U tensor
    lowsymmetry = true;
    // Get anisotropy
    OrthogonalAnisotropy orthogonalanisotropy(isigi);
    principalaxes = orthogonalanisotropy.EigenVectorsOrth(); // store directions
  }
}
//--------------------------------------------------------------------------
  void AnisotropicAnalysis::SetPrincipalDirectionsGeneral
  (const hkl_unmerge_list& hkl_list,
   const int& datasetindex,
   const SDmodel& SDM)
  // set principalaxes from data
  {
    // merged list for given dataset
    MergedList mergedlist(hkl_list, SDM);
    clipper::HKL_data<clipper::data32::I_sigI>& isigi =
      mergedlist.ImeanForDataset(datasetindex);
    // Get anisotropy
    OrthogonalAnisotropy orthogonalanisotropy(isigi);
    principalaxes = orthogonalanisotropy.EigenVectorsOrth(); // store directions
  }
//--------------------------------------------------------------------------
  clipper::Coord_reci_frac AnisotropicAnalysis::NormaliseVector(const clipper::Coord_reci_frac& crdrf) const
  // Normalise to a modulus of 1.0
{
  return clipper::Coord_reci_frac(Vec3<>(crdrf).unit());
  /*
  double maxval = 0.0;
  for (int i=0;i<3;++i) { // get maximum value
    if (std::abs(crdrf[i]) > maxval) {maxval = std::abs(crdrf[i]);}
  }
  DVect3  v3(crdrf);
  v3 = (1.0/maxval) * v3;
  return clipper::Coord_reci_frac(v3); // scaled vector
  */
}
  // ------------------------------------------------------------
  //! set AnisotropicAxis angle in degrees
  void AnisotropicAnalysis::SetConeAngle(const double& angledegrees)
  {
    coneangle = clipper::Util::d2rad(angledegrees);
    cosconeangle = cos(coneangle);
    sinconeangle = sin(coneangle);
  }  
  // ------------------------------------------------------------
  //! return cone angle in degrees
  double AnisotropicAnalysis::ConeAngle() const
  {
    return clipper::Util::rad2d(coneangle);
  }
// ------------------------------------------------------------
std::pair<int,double> AnisotropicAnalysis::Axis(const Hkl& hkl, const Rtype& invresolsq) const
  // return nearest principal axis index, if within cone around that axis,
  // -1 if outside
  {
    if (cubic) {
      return std::pair<int,double>(0,0.0);   // cubic, isotropic so always return index 0
    }
    double dstar = sqrt(invresolsq);  // 1/d
    int jaxis = -1;  // nearest axis, if with cone

    DVect3 projections = Projection(hkl); // 1 or 3 projections
    double cosang;
    double wt = 0.0;   // weight = (cosang - cosconeangle)/(1-cosconeangle)

    if (abplane) {
    //   abplane, just get angle to axis and test against 0 and 90deg
    //     projections[0] is projection on to axis
      cosang = Min(1.0,std::abs(projections[2])/dstar);
      //^
      //      std::cout << hkl.format() <<" " <<dstar<<" "<<projections[2]/dstar
      //		<<" " <<ang <<" angles\n"; 
      //^-
      if (cosang > cosconeangle) { // within cone around axis?
	jaxis = 2;
      } else if (Close(cosang, 0.0, sinconeangle)) {
	// within coneangle of 90 degrees from axis
	jaxis = 0;
	cosang = sqrt(1.0 - cosang*cosang);
	//^	std::cout << hkl.format() <<" " << clipper::Util::rad2d(ang) <<" plane\n"; //^
      }
    } else {
      double maxcosang = 0.0; // cos angle
      for (int i=0;i<3;++i) { //angles from each axis
	// dot products of x with principal axis vectors (unit length)
	// dot products, projections, have length = d*
	cosang = Min(1.0,std::abs(projections[i])/dstar);
	if (cosang > cosconeangle && cosang > maxcosang) {
	  jaxis = i;
	  maxcosang = Max(maxcosang, cosang);
	}
      }
      cosang = maxcosang;
    }
    if (jaxis >= 0) {
	wt = (cosang - cosconeangle)/(1-cosconeangle);
	//^
	//	std::cout <<"coneaxis wt " <<hkl.format()<<" "<<clipper::Util::rad2d(acos(cosang))
	//		  <<" "<<wt<<"\n"; //^-
    }
    return std::pair<int,double>(jaxis, wt);
  }
  // ------------------------------------------------------------
  DVect3 AnisotropicAnalysis::Projection(const Hkl& hkl, const bool& doplane) const
  // return projections on to principal directions
  // 1) if (!abplane), return 3 projections
  // 2) if (abplane),
  //   if doplane false, return projection on to axis in element [2]
  //   if doplane true, return also projection on to plane in elements [0] and [1]
  {
    // orthogonalised coordinate, Clipper convention
    clipper::Coord_reci_orth horth =  hkl.HKL().coord_reci_orth(ccell);
    DVect3 projections(0.,0.,0.);

    if (abplane) { // For testing against axis & plane, abplane = true
      projections[2] = horth * principalaxes[2];
      if (doplane) { // radius of projection on to plane
	projections[0] = sqrt(Max(0.0,horth*horth - projections[2]*projections[2]));
	projections[1] = projections[0];
      }
    } else {
      // 3 directions to test
      for (int i=0;i<3;++i) { //angles from each axis
	// dot products of x with principal axis vectors (unit length)
	projections[i] = horth * principalaxes[i];
      }
    }
    return projections;
  }
  //--------------------------------------------------------------------------
  std::string AnisotropicAnalysis::formattype() const
  //! format types of analyses done
  {
    if (cubic) {
      return "isotropically in cubic symmetry";
    }
    std::string s =
      "within an maxangle of "+StringUtil::Strip(clipper::String(ConeAngle()))+" degrees";
     
    if (abplane) {
      if (rlattice) {
	s+= " perpendicular to and along the (111) axis";
      } else {
	s += " of h k plane and of l axis";
      }
    } else {
      s+= " of principal axes of anisotropy";
    }
    s += "\nweighted according to angle, w = [cos(angle) - cos(maxangle)]/[1 - cos(maxangle)],\n";
    return s;
  }
//--------------------------------------------------------------------------
//! return labels for three axes
std::vector<std::string> AnisotropicAnalysis::Axesformat() const
{
  std::vector<std::string> labels(3);
  // a*, b*, c*
  labels[0] = "h axis";
  labels[1] = "k axis";
  labels[2] = "l axis";
  if (abplane) {
    labels[1] = "";
    if (rlattice) {
      labels[2] = "(1 1 1) axis, diagonal 3-fold";
      // plane perpendicular to (111)
      labels[0] = "plane perpendicular to (111) axis";
    } else {
      // plane perpendicular to l, set 1st label
      labels[0] = "h k plane";
    }
  } else if (lowsymmetry) {
    double tol = 0.001;
    for (int j=0;j<3;++j) {
      clipper::Coord_reci_orth crdro(principalaxes[j]);
      clipper::Coord_reci_frac crdrf = crdro.coord_reci_frac(ccell);
      //^
      //      std::cout << "principalaxis " << principalaxes[j].format() <<"\n";
      //      std::cout << "crdro " << crdro.format() <<"\n";
      //      std::cout << "crdrf " << crdrf.format() <<"\n\n";
      //^-
      int jz = 0;  // count "zeroes"
      for (int i=0;i<3;++i) {
	if (std::abs(crdrf[i]) < tol) {jz++;}
      }
      // jz == 2 if along reciprocal cell axis, leave alone
      if (jz < 2) {
	// not along axis
	labels[j] = formatHvector(crdrf);
      }
    } // end loop axes
  }
  return labels;
}
//--------------------------------------------------------------------------
  std::string AnisotropicAnalysis::formatHvector(const clipper::Coord_reci_frac& crdrf) const
  //! return formatted normalised vector in hkl space (fractional)
{
  clipper::Coord_reci_frac nc = NormaliseVector(crdrf);
  std::string s;
  const double TOL = 0.0001; 
  bool started = false;
  char chkl[] = {'h','k','l'};
  for (int i=0;i<3;++i) {
    std::string sign = "";
    if (std::abs(nc[i]) > TOL) {
      double c = nc[i];
      if (started) {
	sign = "+";
	if (nc[i] < 0.0) {
	  sign = "-";
	  c = std::abs(c);
	}
      }
      if (sign != "") {s += " "+sign+" ";}
      s += StringUtil::Strip(StringUtil::ftos(c,8,2))+" "+chkl[i];
      started = true;
    }}
  return s;
}
}
