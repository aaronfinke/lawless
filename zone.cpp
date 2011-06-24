// zone.cpp
#include <iostream>

#include "zone.hh"
#include "scala_util.hh"
#include "probfunctions.hh"
#include "score_datatypes.hh"
#include "string_util.hh"

#include <assert.h>
#define ASSERT assert

namespace scala
{  //--------------------------------------------------------------
  double IscoreVal(const IsigI& Isig)
  // Score function value from I,sig
  //  = I/sigI 
  {
    if (Isig.I() < 0.0) {return 0.0;}
    return Isig.I()/(Isig.sigI());
    // I
    //    return Isig.I();
  }
  //--------------------------------------------------------------
  OneDFourier::OneDFourier(const std::vector<int>& Ngrid)
  {
    ngrid = Ngrid;
    xgrid.resize(Ngrid.size());
    Sumx.resize(Ngrid.size());

    for (size_t is=0;is<ngrid.size();is++)
      {
	Sumx[is] = 0.0;
	if (ngrid[is] > 1)
	  xgrid[is] = 1.0/double(ngrid[is]);
	else
	  xgrid[is] = 0.0;
      }
    nobs = 0;
  }
  //--------------------------------------------------------------
  OneDFourier::OneDFourier(const std::vector<double>& Xgrid)
  {
    xgrid = Xgrid;
    ngrid.resize(xgrid.size());
    Sumx.resize(xgrid.size());

    for (size_t is=0;is<ngrid.size();is++)
      {
	Sumx[is] = 0.0;
	if (xgrid[is] > 1)
	  ngrid[is] = Nint(1.0/xgrid[is]);
	else
	  ngrid[is] = 0;
      }
    nobs = 0;
  }
  //--------------------------------------------------------------
  void OneDFourier::AddRef(const int& j, const double& Val)
  //  Fourier component = 2 pi j . x

  //  j is 1-D index
  {
    double aj = j;
    double twopij = clipper::Util::twopi() * aj;

    for (size_t is=0;is<ngrid.size();is++)
      {
	Sumx[is] += Val * cos(twopij * xgrid[is]);
      }
    nobs++;
  }
  //--------------------------------------------------------------
    // Return number of contributions
  int OneDFourier::Nobs() const {return nobs;}
  //--------------------------------------------------------------
  // Return Fourier values at each grid point
  // normalised by first point (at 0)
  std::vector<double> OneDFourier::FourierVal() const
  {
    std::vector<double> F(ngrid.size(),0.0);
    if (nobs > 0)
      {
	double F0;
	for (size_t is=0;is<ngrid.size();is++)
	  {
	    F[is] = Sumx[is]/double(nobs);
	    if (is == 0) F0 = F[is];   // normalise to zero'th term
	    if (F0 > 0.0) {
	      F[is] = F[is]/F0;
	    } else {
	      F[is] = 0.0;
	    }
	    //^
	    //^	    std::cout << "1DFV " << is << " " << F[is] << "\n";
	  }
      }
    return F;
  }
  //--------------------------------------------------------------
  // static constants
  float Zone::neighbourFraction = 0.02;   // default
  double Zone::MinControlSD = 0.05;
  //--------------------------------------------------------------
  Zone::Zone(const std::string& Axis, const int& Nfold,
	     const hkl_symmetry& LGsym)
  // Initialise from Axis (a|b|c) & order (2,3,4,6)
  // Nfold negated to test only the 1/Nfold point
  // (sets validpoint array)
  {
    axis = true;
    diagonal = false;
    direction = Axis;
    glide = "";
    order = std::abs(Nfold);
    bool disable = false;
    if (Nfold < 0) disable = true;

    // singleprob = true if only the last probability (Fourier)
    // point is useful. This is the case for all glides
    // and 2- & 3-fold screws
    if (order > 3)
      {singleprob = false;}
    else
      {singleprob = true;}

    // Vector to convert input hkl to Fourier index
    //  set here to 001 (ie 00l), permuted in "init"
    cond = clipper::Vec3<int>(0,0,1);  


    if (Axis == "a") {
      permute = ReindexOp("k,l,h");
    }
    else if (Axis == "b") {
      permute = ReindexOp("l,h,k");
    }
    else if (Axis == "c") {
      permute = ReindexOp("h,k,l");
    }
    else {
      clipper::Message::message(clipper::Message_fatal
				("Illegal argument to Zone Axis: "+Axis));
    }
    lgsymm = LGsym;  // Laue group symmetry
    init(disable);
  }
  //--------------------------------------------------------------
  Zone::Zone(const std::string& GlidePlane, const std::string& Glide,
	     const hkl_symmetry& LGsym)
  // Initialise from Glide Plane & translation direction
  {
    axis = false;
    diagonal = false;
    direction = GlidePlane;
    glide = Glide;
    order = 2;
    bool disable = false;
    // singleprob = true if only the last probability (Fourier)
    // point is useful. This is the case for all glides
    // and 2- & 3-fold screws
    singleprob = true;
    
    if (GlidePlane =="a") {
      permute = ReindexOp("k,l,h");
    }
    else if (GlidePlane =="b") {
      permute = ReindexOp("l,h,k");
    }
    else if (GlidePlane =="c") {
      permute = ReindexOp("h,k,l");
    }
    else if (GlidePlane =="110") {
      permute = ReindexOp("h,k,l");
      diagonal = true;
    }
    else if (GlidePlane =="101") {
      permute = ReindexOp("l,h,k");
      diagonal = true;
    }
    else if (GlidePlane =="011") {
      permute = ReindexOp("k,l,h");
      diagonal = true;
    }
    else {
      clipper::Message::message(clipper::Message_fatal
				("Illegal argument to Zone Plane: "+GlidePlane));
    }

    // cond is vector to convert input hkl to Fourier index
    //  set here, later permuted in "init"
    // Cases:
    //   Principle zone  (hk0)
    //     a,b,c   010   k=2n
    //     n       110   h+k=2n
    //     d       110   h+k=4n
    //   Diagonal zone (hhl)
    //     c       001   l=2n
    //     n       001   l=2n
    //     d       111   h+k+l=4n

    if (Glide =="a" || Glide =="b" || Glide =="c") 
      {
	if (Glide =="a") {
	  cond = clipper::Vec3<int>(1,0,0);  // Glide shift x+1/2
	}
	else if (Glide =="b") {
	  cond = clipper::Vec3<int>(0,1,0);  // Glide shift y+1/2
	}
	else if (Glide =="c") {
	  cond = clipper::Vec3<int>(0,0,1);  // Glide shift z+1/2
	}
	// Transform cond vector to internal standard frame
	clipper::Vec3<double> v = clipper::Vec3<double>(cond) * permute;
	cond = clipper::Vec3<int>(Nint(v[0]), Nint(v[1]), Nint(v[2]));
      }
    else if (Glide =="n") {
      if (diagonal) {
	cond = clipper::Vec3<int>(0,0,1);  // n(110) Glide shift z+1/2
      }
      else {
	cond = clipper::Vec3<int>(1,1,0);  // n(c) Glide shift x+1/2,y+1/2
      }
    }
    else if (Glide =="d") {
      order = 4;
      disable = true;  // only 4n not 2n
      if (diagonal) {
	cond = clipper::Vec3<int>(1,1,1);  // d(110) Glide shift x,y,z+1/4
      }
      else {
	cond = clipper::Vec3<int>(1,1,0);  // d(c) Glide shift x+1/4,y+1/4
      }
    }
    else {
      clipper::Message::message(clipper::Message_fatal
				("Illegal argument to Zone Glide: "+GlidePlane));
    }
    lgsymm = LGsym;  // Laue group symmetry
    init(disable);
  }
  //--------------------------------------------------------------
  void Zone::init(const bool& disable)
  {
    minIndx =  1000000;
    maxIndx = -1000000;
    // Set up grid intervals
    for (int i=1;i<=order;i++) {
      if (i == 1) {
	ngrid.push_back(i);
      } else if (order%i == 0) {
	ngrid.push_back(i);
      }
    }
    npoint = ngrid.size();

    if (disable) {
      // Nfold negative, validate only 0'th & last point
      validpoint = std::vector<bool>(npoint, false);
      validpoint[0] = true;
      validpoint[npoint-1] = true;
    } else {
      validpoint = std::vector<bool>(npoint, true);
    }

    fsum = OneDFourier(ngrid);  // 1D Fourier
    controlsd.assign(npoint, -1.0);
    controlmean.assign(npoint, 0.0);

    valid = false;
    results = false;
    // Transform cond vector to external frame  [P] q
    clipper::Vec3<double> v = permute * clipper::Vec3<double>(cond);
    cond = clipper::Vec3<int>(Nint(v[0]), Nint(v[1]), Nint(v[2]));
    for (int i=0;i<3;i++) {qcond[i] = cond[i];}

    permuteIndex = permute;
    NonIdentityReindex = false; // flag identity reindex

    prob_yes = -1.0;
  }
  //--------------------------------------------------------------
  bool Zone::InZone(const Hkl& hkl) const 
  // returns true if hkl is in defined zone
  //  hkl is the original index in reference (lattice) frame 
  {
    // Put into the Laue group frame
    Hkl hklL = hkl.change_basis(reindexmat);
    // Reduce to asymmetric unit in Laue group
    // This is needed since eg in high symmetry Laue groups some axes
    // are equivalent, while in lower symmetry groups they are distinct
    //%    if (axis) {
      int isym;
      hklL = lgsymm.put_in_asu(hklL, isym);
      //%    }
    // Transform from Laue group frame
    // Permute to internal standard setting (eg axis is always 00l)
    Hkl h = hklL.change_basis(permute);
    if (axis) {
      // Axis, internally test axis is 00l
      if (h[0] == 0 && h[1] == 0) {
	return true;
      }
    } else if (diagonal) {
      // Glide on diagonal, internally hhl
      if (h[0] == h[1]) return true;
    } else {
      // Glide on principle zone, internally hk0
      if (h[2] == 0) return true;
    }
    return false;
  }
  //--------------------------------------------------------------
  int AxisDirection(const clipper::Vec3<double> vaxis) 
  //  Returns:
  //   jc   axis direction (0,1,2) if principle axis, -1 if not
  {
    int nc=0;
    int jc=-1;
    // Count non-zero elements
    for (int i=0;i<3;i++)
      {
	if (std::abs(vaxis[i]) > 0.02) {
	  nc++;
	  jc = i;
	}
      }
    
    // Not principle	  
    if (nc != 1) jc = -1;
    return jc;
  }
  //--------------------------------------------------------------
  clipper::Vec3<double> Zone::DirectionRefFrame() const
  // Get zone direction in reference frame
  //  vector returned is axis direction or zone normal (for glide plane)
  {
    // Axis in reference frame is (001)[P]^-1[H]^-1 = (001)[HP]^-1
    // Glide plane normal in reference frame is same
    clipper::Vec3<double> laxis(0.0,0.0,1.0);
    clipper::Vec3<double> vaxis = laxis * permuteIndex.inverse();
    return vaxis;
  }
  //--------------------------------------------------------------
  std::string Zone::formatRefFrame(const int Sub) const
  // Format in reference frame
  {
    return formatNewFrame(ReindexOp(), Sub);
  }
  //--------------------------------------------------------------
  std::string Zone::formatLGFrame(const int Sub) const
  // Format in Laue group frame
  {
    return formatNewFrame(reindexmat, Sub);
  }
  //--------------------------------------------------------------
  clipper::Vec3<int> IntVector(const clipper::Vec3<double> v)
  {
    const double tol=0.001;
    // Make integer version of vector by multiplication
    if (v.is_null())
      {return clipper::Vec3<int>(0,0,0);}
    // Unit vector
    clipper::Vec3<double> v1 = v.unit();
    // Find smallest non-zero element
    double small = 1.0;
    for (int i=0;i<3;i++)
      {if (std::abs(v1[i]) > 0.001) {small = Min(small, std::abs(v1[i]));}}
    double sc = 1./small;
    for (int i=0;i<3;i++)  {v1[i] *= sc;}
    bool integral = true;
    for (int i=0;i<3;i++)  
      {if (!Close(v1[i], double(Nint(v1[i])), tol)) integral = false;}
    if (!integral)
      {clipper::Message::message(clipper::Message_fatal
				 ("IntVector: non-integral vector "+v1.format()));}
    clipper::Vec3<int> iv;
    for (int i=0;i<3;i++) {iv[i] = Nint(v1[i]);}
    return iv;
  }
  //--------------------------------------------------------------
  std::string DirectionFormat(const int& jc, const clipper::Vec3<double>& vd,
			      const bool& diagonal)
  // jc = 0,1,2 for principle direction
  // vd   direction
  // diagonal   true if diagonal
  {
    std::vector<char> sabc(3);
    sabc[0] = 'a';
    sabc[1] = 'b';
    sabc[2] = 'c';
    // Diagonals
    std::vector<std::string> sdiag(3);
    sdiag[0] = "011";
    sdiag[1] = "101";
    sdiag[2] = "110";
    std::string dir;
    if (jc >= 0) {
      // Principal direction
      if (diagonal)
	{dir = sdiag[jc];}
      else
	{dir = sabc[jc];}
    } else {
      // Try to make integer version of direction
      clipper::Vec3<int> iv = IntVector(vd);
      for (int i=0;i<3;i++) {
	dir += clipper::String(iv[i],4);
	//	  if (i<2) dir +=",";
      }
    }
    return StringUtil::Strip(dir);
  }
  //--------------------------------------------------------------
  std::string Zone::formatNewFrame(const ReindexOp& RefToNew,
				   const int Sub) const
  // Format in "new" frame
  // RefToNew is reindex operator from "reference" (lattice) frame
  // (stored with StoreReindex function) to some new frame in
  // which to format the zone
  // if Sub > 1, the screw component is order/Sub
  {
    // Direction of axis or glide in new frame
    clipper::Vec3<double> vd = DirectionRefFrame() * RefToNew;
    int jc = AxisDirection(vd);
    std::string dir = DirectionFormat(jc, vd, diagonal);
    std::string dref = "";
    if (!RefToNew.IsIdentity()) {
      // print lattice direction as well
      dref = DirectionFormat(AxisDirection(DirectionRefFrame()),
			   DirectionRefFrame(), diagonal);
    }      

    if (axis) {
      // Axis
      std::string l = "screw axis "+
	clipper::String(order,1)+"("+clipper::String(order/Sub,1)+
	") ["+dir+"]";
      if (dref == "") {
	l += "    ";
      } else {
	l += " ("+dref+")";
      }
      return l;
    } else {
      // Glide
      // Direction of axis or glide in new frame
      std::string gld = glide;
      clipper::Vec3<double> vg = qcond * RefToNew;
	int jg = AxisDirection(vg);
	if (jg >= 0)
	  {gld = DirectionFormat(jg, vg, false);}
	
	return "glide plane "+gld+"("+dir+")";
    }
  }
  //--------------------------------------------------------------
  template<class T>  std::string Cond(const clipper::Vec3<T>& cond,
				      const int mult)
  // format condition as eg h+k=2n
  {
    std::vector<char> shkl(3);
    shkl[0] = 'h';
    shkl[1] = 'k';
    shkl[2] = 'l';
    std::string sc;
    bool first = true;
    // Force all elements positive
    clipper::Vec3<T> cv(cond);
    for (int i=0;i<3;i++) cv[i] = std::abs(cond[i]);

    for (int i=0;i<3;i++)
      {
	if (!Close(cv[i], T(0), T(0.01)))
	  {
	    if (!first) sc += "+";
	    if (!Close(cv[i], T(1), T(0.01))) 
	      {
		sc += clipper::String(cv[i]);
	      }
	    sc += shkl[i];
	    first = false;
	  }
      }
    sc += "="+clipper::String(mult,1)+"n";
    return StringUtil::Strip(sc);
  }
  //--------------------------------------------------------------
  // Format condition in reference frame
  std::string Zone::FormatConditionRefFrame(const int Sub) const
  {
    return FormatConditionNewFrame(ReindexOp(), Sub);
  }
  //--------------------------------------------------------------
  // Format condition in Laue group frame
  std::string Zone::FormatConditionLGFrame(const int Sub) const
  {
    return FormatConditionNewFrame(reindexmat, Sub);
  }
  /*?

  std::string Zone::FormatCondition(const int Sub) const
  // Format condition in constructor frame, eg "0kl: h+k = 2n"
  {
  if (axis && Sub > 0 && !validpoint[Sub])
  {
  return "";  // Obscured point by lattice absences
  }

  std::vector<char> shkl(3);
  shkl[0] = 'h';
  shkl[1] = 'k';
  shkl[2] = 'l';
  std::string sc;

  Hkl hz;
  if (axis)
  // Axis: direction hz is axis
  {
  hz = Hkl(0,0,1).change_basis(permute.inverse());
  }
  else if (diagonal)
  // Diagonal glide: direction hz is normal
  {
  hz = Hkl(1,1,0).change_basis(permute.inverse());
  }
  else
  // Principle glide: direction hz is zone
  {
  hz = Hkl(1,1,0).change_basis(permute.inverse());
  }

  if (axis || !diagonal)
  {
  for (int i=0;i<3;i++)
  {
  if (hz[i] == 0)
  sc += "0";
  else
  sc += shkl[i];
  }
  }
  else
  // Diagonal glide, just list them
  {
  if (direction == "110")
  {sc = "hhl";}
  else if (direction == "101")
  {sc = "hkh";}
  else if (direction == "011")
  {sc = "hkk";}
  }

  sc += ": ";

  int mult = order;
  if (axis && Sub > 0)
  {
  mult = ngrid[Sub];
  }
  sc += Cond<int>(cond, Sub, mult);
  return sc;
  }?*/
  //--------------------------------------------------------------
  std::string Zone::FormatConditionNewFrame(const ReindexOp& RefToNew,
					    const int Sub) const
  // Format condition in newframe
  // RefToNew is reindex operator from "reference" (lattice) frame
  // (stored with StoreReindex function) to some new frame in
  // which to format the zone
  {
    if (axis && Sub > 0 && !validpoint[Sub])
      {
	return "";  // Obscured point by lattice absences
      }

    std::vector<char> shkl(3);
    shkl[0] = 'h';
    shkl[1] = 'k';
    shkl[2] = 'l';
    std::string sc;
    // Axis in reference frame is (001)[P]^-1[H]^-1 = (001)[HP]^-1
    // Glide plane normal in reference frame is same
    clipper::Vec3<double> vaxis = DirectionRefFrame() * RefToNew;
    int jc = AxisDirection(vaxis);

    std::string sax;
    std::string sgl;
    char zg = ' ';

    if (jc >= 0)
      {
	for (int i=0;i<3;i++)
	  {
	    vaxis[i] = std::abs(vaxis[i]);  // force positive
	    if (vaxis[i] > 0.01)
	      {
		sax += shkl[i];
		if (diagonal)
		  {sgl += shkl[i];}
		else
		  {sgl += '0';}
	      }
	    else
	      {
		sax += '0';
		if (diagonal)
		  {if (zg == ' ') zg = shkl[i];
		    sgl += zg;
		  }
		else
		  {sgl += shkl[i];}
	      }
	  }
      }
    else
      // Not principle	  
      {
	sax = "[";
	// Try to make integer version of direction
	clipper::Vec3<int> iv = IntVector(vaxis);
	for (int i=0;i<3;i++) {
	  sax += clipper::String(iv[i],3);
	  // if (i<2) sax +=",";
	}
	sax += "]";
	sax = StringUtil::Strip(sax);
	sgl = sax;
      }

    if (axis)
      {
	sc += sax;
      }
    else
      {
	sc = sgl;
      }
    
    int mult = order;
    if (axis && Sub > 0)
      {
	mult = ngrid[Sub];
      }

    sc += ": " + Cond<double>(RefToNew.inverse()*qcond, mult);

    return sc;
  }
  //--------------------------------------------------------------
  // Axis direction in reference frame
  std::string Zone::Direction() const
  {
    // Direction of axis or glide in new frame
    clipper::Vec3<double> vd = DirectionRefFrame();
    int jc = AxisDirection(vd);
    return DirectionFormat(jc, vd, diagonal);
  }
  //--------------------------------------------------------------
  bool Zone::SameZone(const Zone& other) const
    // Return true if two zones apply to the same hkl zone
    // Glides only, eg b(a), c(a), n(a) all apply to 0kl
    {
      // It is sufficient to compare directions, but after reindexing
      bool same = true;
      clipper::Vec3<double> vd1 = DirectionRefFrame();
      clipper::Vec3<double> vd2 = other.DirectionRefFrame();
      for (int i=0;i<3;i++)
	{if (std::abs(vd1[i]) != std::abs(vd2[i]))
	    {same = false;}
	}
      return same;
    }
  //--------------------------------------------------------------
  // Add into totals of Fourier components
  // For glide planes, do this immediately
  // For axes, just store I, sigI. Fourier synthesis is done in function
  //  GetFourierValues, to allow for neighbour corrections
  void Zone::AddRef(const Hkl& hkl, const IsigI& Isig, const float& sSqr)
  {
    // Put into the Laue group frame
    Hkl hklL = hkl.change_basis(reindexmat);
    // If axis, reduce to asymmetric unit in Laue group
    // This is needed since eg in high symmetry Laue groups some axes
    // are equivalent, while in lower symmetry groups they are distinct
    // For Glide planes, leave as is - NO don't
    //%    if (axis) {
    int isym;
    hklL = lgsymm.put_in_asu(hklL, isym);
    //%    }
    // Transform back to lattice frame
    Hkl h = hklL.change_basis(reindexmat.inverse());
    
    InvResRange.update(sSqr);
    // Fourier index j = [h].[qcond]
    //  Fourier component = 2 pi j . x
    int j = Nint(clipper::Vec3<double>(h) * qcond);
    if (axis) {
      // Store IsigI values for axis only
      IndxIsigI.push_back(IndexIsigI(j,Isig));
    } else {
      // Glide plane, add in immediately
      fsum.AddRef(j, IscoreVal(Isig));
    }
    minIndx = Min(minIndx, std::abs(j));
    maxIndx = Max(maxIndx, std::abs(j));
  }
  //--------------------------------------------------------------
  // Store SD of control reflections for Z-score
  void Zone::StoreMeanSD(const std::vector<double>& Mean,
			 const std::vector<double>& SD)
  {
    ASSERT (int(Mean.size()) == npoint && int(SD.size()) == npoint);
    ZeroControlSD = false;
    UnitControlMean.assign(npoint, false);

    for (int i=0;i<npoint;++i) {
      if (SD[i] < 0.00001) {
	// Set flag to indicate a zero SD
	ZeroControlSD = true;
	controlsd[i] = Max(MinControlSD, SD[i]);
	if (i > 0 && std::abs(Mean[i] - 1.0) < 0.001) {
	  UnitControlMean[i] = true;
	}
      } 
      controlsd[i] = Max(MinControlSD, SD[i]);
    }
    controlmean = Mean;
  }
  //--------------------------------------------------------------
  std::vector<IsigI> Zone::AverageI(const std::vector<IndexIsigI>& IdxIs) const
  // Average replicant (& inverse) axial observations
  // On entry:
  //  IdxIs      vector of index,I,sigma for each observation
  //
  // Returns:
  //  Is   vector of average IsigI indexed by axial index, (0,0) if no value
  {
      std::vector<double> mnI(maxIndx+1, 0.0);
      std::vector<double> wt(maxIndx+1, 0.0);
      std::vector<IsigI> Is(maxIndx+1);

      for (size_t i=0;i<IdxIs.size();i++) {
	IsigI Isig = IdxIs[i].Isig;
	int j = std::abs(IdxIs[i].index);
	double w = 1./(Isig.sigI()*Isig.sigI());
	mnI.at(j) += w * Isig.I();
	wt.at(j) += w;
      }
      for (int j=0;j<maxIndx+1;j++) {
	if (wt.at(j) > 0.0) {
	  Is.at(j) = IsigI(mnI.at(j)/wt.at(j), 1./sqrt(wt.at(j)));
	} else {
	  Is.at(j) = IsigI(0.0, 0.0);
	}
      }
      return Is;
  }
  //--------------------------------------------------------------
  std::vector<double> Zone::GetFourierValues() const
  {
    // Return normalised Fourier values
    // For axes, add in Fourier coefficients first: already done for glide planes

    // For axes, "correct" each intensity by subtracting a fraction of neighbouring
    // intensities, if present, add into Fourier sums
    if (axis) {
      // get average I for each index
      averageIsigI = AverageI(IndxIsigI);
      std::vector<IndexIsigI> adjustedIndxIsigI(IndxIsigI.size());
      // Correct intensities
      for (size_t i=0;i<IndxIsigI.size();i++) {
	IsigI Isig = IndxIsigI[i].Isig;
	float I = Isig.I();
	//^
	//^	float I0 = I;
	//^	double v0 =  IscoreVal(Isig);
	//^-
	int j = std::abs(IndxIsigI[i].index);
	if (I > 0.0) {
	  // No correction if already negative
	  if (j > 0) {
	    if (averageIsigI[j-1].I() > 0.0) {
	      I -= neighbourFraction * averageIsigI[j-1].I();
	      I = Max(0.0, I); // don't allow to go negative
	    }
	  }
	  if (j < maxIndx) {
	    if (averageIsigI[j+1].I() > 0.0) {
	      I -= neighbourFraction * averageIsigI[j+1].I();
	      I = Max(0.0, I); // don't allow to go negative
	    }
	  }
	  Isig.I() = I;
	} else {
	  // I <= 0
	  Isig.I() = 0.0;
	}
	adjustedIndxIsigI[i] = IndexIsigI(j, Isig);  // store adjusted values
	fsum.AddRef(j, IscoreVal(Isig));
	//^
	//^	std::cout << "I, I' " << j << " " << I0 << " " << I << " " 
	//^		  << v0 << " " << IscoreVal(Isig) << "\n";
      }
      // Average adjusted values for printing
      adjustedIsigI = AverageI(adjustedIndxIsigI);
    }  // end if axis
    return fsum.FourierVal();
  }
  //--------------------------------------------------------------
  // Store reindex from some standard reference frame, for comparison of zones
  void Zone::StoreReindex(const ReindexOp& reindexin) 
  {
    reindexmat = reindexin;
    NonIdentityReindex = ! reindexmat.IsIdentity(); // flag non-identity
    // update qcond, ql = [H] qc
    clipper::Vec3<double> v(cond);
    qcond = reindexmat * v;
    permuteIndex = reindexmat * permute;
  }
  //--------------------------------------------------------------
  // Return Fourier values at each grid point
  std::vector<double> Zone::FourierVal() const
  {
    return GetFourierValues();
  }
  //--------------------------------------------------------------
  void Zone::TestValidIndices() const
  {
    if ((axis && IndxIsigI.size() < 2) ||
	(maxIndx <= minIndx))
      {
	valid = false;
      }
    else
      valid = true;
  }
  //--------------------------------------------------------------
  // Internal function to calculate scores
  void Zone::CalcResults() const
  // Prerequisites:
  //  1) SD(Fouriervalue) controlsd
  //  2) Control Mean  controlmean
  //
  // This function sets:
  //  1) valid flag (in TestValidIndices)
  //  2) Pfor[i], i=1,npoint
  //  3) If singleprob, sets prob_yes
  //  4) results = true
  //
  // Algorithm:
  //  We have a calculated value of the Fourier coefficient(s),
  //  expressed as a fraction of the origin, and also an
  //  estimate of the standard deviation from the control sample.
  //  
  //  We want to calculate the probability that the reflection
  //  condition is true, ie the reflections are systematically
  //  absent, as opposed to not true. In this calculation, we
  //  need to allow for the fact that the Fourier coefficients
  //  are not independent.
  //
  //  Write the putative screw axis as M(q), with a pure rotation
  //  represented with q = M  (eg 6(6) == 6(0)). Glide planes
  //  behave in the same way as 2(1) axes (even for the 4n
  //  d glides, where we only analyse Fourier point at 1/4,
  //  not at 1/2).
  //
  //  Fourier analysis is done at N points corresponding to the
  //  factors of M, eg for M=6, at points 1/2, 1/3, and 1/6.
  //  Write the Fourier peak heights (relative to the origin) as
  //  vector v(j), j=1,N, with corresponding sample orders
  //  p(j) eg (2,3,6) for M=6.
  //
  //  For M>3, the values v(j) are not independent, so we need to
  //  test them simultaneously. This can be done by considering
  //  the "distance" between the score vector v and its ideal value
  //  for each possible screw component q, a vector e(q) = (e(q)(j))
  //
  //  The ideal values are either 0 or 1, according to whether the
  //  Fourier order corresponds to the screw order.
  //  e(q)(j) = E(v(j) | q) =  1  if M/(q p(j)) is integral
  //                        =  0  if M/(q p(j)) is not integral
  //  thus, for a 4-fold
  //
  //                        j    = 1          2 
  //                        p(j) = 2          4
  //   q   M/q  Condition         0.5        0.25   Fourier position
  //   4    1      1n             1/2(0)     1/4(0)     M/qp(E(v))
  //   2    2      2n             2/2(1)     2/4(0)
  //   1    4      4n             4/2(1)     4/4(1)
  //  
  //  Then the "distance" between measured & ideal values is
  //  d = || v - e(q) ||            dmax = Sqrt(N)
  //
  //  For M = 2 or 3, N = 1 and the vectors reduce to a single value.
  //
  // 6-folds are a bit more complicated
  //   M = 6
  //                      j    = 1          2         3
  //                      p(j) = 2          3         6
  //   q   M/q  Condition      0.5       0.333      0.1666      Fourier position
  //   6    1      1n         1/2(0)     1/3(0)     1/6(0)      M/q p(E(v))
  //   3    2      2n         2/2(1)     2/3(=1/6)  2/6(=1/3)   point at 1/3 = point at 1/6
  //   2    3      3n         3/2(=1/6)  3/3(1)     3/6(=1/2)   point at 1/2 = point at 1/6
  //   1    6      6n         6/2(1)     6/3(1)     6/6(1)
  //
  // ie for 6(2) and 6(3) we need to test for equality between 2 of the Fourier points
  //    q               e1  e2  e3   target values at 1/2, 1/3, 1/6
  //         p(j)        2   3   6
  //    6  6(0) l=1n     0   0   0
  //    3  6(3) l=2n     1  v3  v2   v2 is value at j=2 etc
  //    2  6(2) l=3n     v3  1  v1
  //    1  6(1) l=6n     1   1   1
  // thus eg for q=3, 6(3),  d^2 = (v1-1)^2 + 2*(v2 - v3)^2  (factor of 2 to get v2-v3 & v3-v2)
  //
  //  Given d, a probability distribution can be constructed to give
  //  P(Aq | d). Here a Lorentzian distribution around 0 is used with a
  //  "width" parameter = SD estimate from random reflections.
  //
  //  Try the following
  // *  P(!A | d) is calculated assuming the "ideal" value of the Fourier
  // *  value v is not necessarily = 0 if the systematic absence condition
  // *  is not true, ie allowing for pseudosymmetry
  // *  But we can integrate over all possible values of the ideal
  // *  v = m, 0 <= m < dmax, with some model for the prior probability of m
  // *  with a high value = 1 = P(m=dmax) and a low value = 0 = P(m=0)
  //
  // *  p(!A|m) = N(m,sigma)
  // *  p(!A|v) = Integral{p(!A|v,m) p(m) dm}
  //
  // *  We can model p(m) =  eg. Sqrt(1 - (1-m)^2) and numerically integrate
  // *  Then normalise to make P(A|v) + P!A|v) = 1
  // *  The values are then renormalised by making the sum over all values
  // *  of q
  //    
  {
    if (results) return;   // already calculated
    TestValidIndices();
    std::vector<double> V(1,0.0);
    prunedData.assign(npoint, false);

    // Normalised Fourier values (V[0] = 1.0 unless no valid data)
    if (valid) {
      V = GetFourierValues();
      ASSERT (int(V.size()) == npoint);
    }

    bool DEBUG = false;
    double average_controlsd = 0.0;
    if (valid) {
      int ninvalid = 0;
      for (int i=1;i<npoint;i++) {
	if (controlsd[i] < -0.0001) {
	  clipper::Message::message(clipper::Message_fatal
				    ("Zone: "+
				     formatRefFrame()+": results requested with SD unset"));
	}
	// Check for strange values of the control mean, which probably indicate a systematic
	// non-random sample of indices, eg all odd ones already eliminated
	if (DEBUG) {
	  std::cout << "ZoneCalc "
		    << controlsd[i] << " " << controlmean[i]
		    << "  V[i] = " << V[i] << "\n";
	}
	if (UnitControlMean[i]) {
	  prunedData[i] = true;
	  ninvalid++;
	}
	average_controlsd += controlsd[i];
      }
      // Zone is invalid if all points are invalid
      if (ninvalid == npoint) {valid = false;}
      average_controlsd /= double(npoint-1);
    }

    if (!valid || V[0] <= 0.0)  {
      // No data or invalid data
      Pfor.resize(npoint);
      for (int i=0;i<npoint;i++)	{
	Pfor[i] = -1.0;
      }
      prob_yes = -1.0;
      return;
    }

    V[0] = 0.0; // not interested in this!
    double Ptot = 0.0;
    Pfor.resize(npoint);
    int Nf = npoint-1;
    double dmax = sqrt(double(Nf));
    // Offset Lorentzian probability so that P(dmax) = 0 
    double Poffset = TruncatedLorentzianProb(dmax, 0.0, average_controlsd, 0.0, dmax);
    //^ 
    //^    Poffset = 0.0; //^!!!!

    // Prior distribution of E(v|!A)
    //DM_2sqrt DM(0.0, dmax);  // Low probability at dmax, high at 0
    //    DM_1minusmSq DM(0.0, dmax);  // Low probability at dmax, high at 0
    DM_1minusmCu DM(0.0, dmax);  // Low probability at dmax, high at 0
    IntgrtProb IP(DM);      // object to do integrated probability
    if (DEBUG) {std::cout << "Zone calc: npoint " << npoint << "\n";}

    // We have Nf = (npoint-1) non-independent points in Fourier space to consider
    // writing these values as v1, v2, ... v(Nf), and considering these as axes
    // running from 0 to 1, then for each possible screw value (it only matters
    // for 4- & 6-fold screws), there is an optimum point {E(v1), E(v2), ...},
    // Get "distance" from optimum point for each possible screw (glide) value i
    //
    // Special for 6-folds:
    //  for 6(3) E(v2) = E(v3) so use v2-v3 instead of v2-E(v2), etc
    //  for 6(2) E(v1) = E(v3) so use v1-v3 instead of v1-E(v1), etc
    //
    for (int i=0;i<npoint;i++) { // Loop M(q) screw or glide translation
      int q = order/ngrid[i];
      if (DEBUG) {std::cout << "\ni = " << i << " q = " << q <<"\n";}
      double d = 0.0;
      int np = 0;
      if (order == 6) { // special for 6-fold   -------------------> 6
	if (q == 1) { // q = 1, link v1, v2 and v3, ie j = 1,2,3
	  double av = 0.0;
	  for (int j=1;j<npoint;j++) {
	    av += V[j];
	    np++;
	  }
	  av /= double(np);
	  // for 6(1) we have two sorts of targets
	  //  (1) v1=v2=v3
	  //  (2) v1=v2=v3=1
	  // weight (1) higher, total weight = 1.0 so put on same scale as other 6-folds
	  double w1 = 2./3.;
	  double w2 = 1. - w1;
	  for (int j=1;j<npoint;j++) {
	    d += w1*(V[j] - av)*(V[j] - av);  // v1,2,3 should be equal
	    d += w2*(1.- V[j])*(1.- V[j]);    //  ... and = 1
	  }
	} else if (q == 2) { // q = 2, link v1 and v3, ie j = 1 and 3
	  d += 2.*(V[1] - V[3])*(V[1] - V[3]);
	  d += (1.- V[2])*(1.- V[2]);  // v2 = 1
	  np +=3;
	} else if (q == 3) { // q = 3, link v2 and v3, ie j = 2 and 3
	  d += 2.*(V[2] - V[3])*(V[2] - V[3]);
	  d += (1.- V[1])*(1.- V[1]);  // v1 = 1
	  np +=3;
	} else {
	  for (int j=1;j<npoint;j++) {
	    d += V[j]*V[j];
	    np++;
	  }
	}
	//                                     -------------------> 6
      } else {
	// Loop each non-origin Fourier point
	for (int j=1;j<npoint;j++) {
	  if (validpoint[j]) {	
	    // Screw condition   M/q, q = M/ngrid[i]
	    // Fourier order          p = ngrid[j]
	    int Movq = ngrid[i];
	    if (Movq%ngrid[j] ==  0) {  // (M/q)%p
	      // M/q is integral multiple of p(j), this Fourier point should be present
	      //  "distance" from expectation value 1 = (1 - v)
	      d += (1.- V[j])*(1.- V[j]);
	      np++;
	    } else {
	      // M/q is not an integral multiple of p(j), this Fourier point should be absent
	      //  "distance" from expectation value 0 = v
	      if (V[j] > 0.0) {
		d += V[j]*V[j];
	      }
	      np++;
	    } // valid point
	  }
	}
      }
      if (np > 0) {
	// Sqrt to get "distance"
	d = Min(dmax,sqrt(d));
	if (i == 0) {
	  // First point q = 0, ie condition absent
	  // Do something different for this one,
	  //     to allow for possibility of pseudosymmetry
	  //   Integrated Gaussian probability around possible values
	  Pfor[i] = IP.LorentzProb(d, average_controlsd, 0, dmax);
	  //		Pfor[i] = IP.Prob(d, controlsd[i], 0, 0);
	  if (DEBUG) {std::cout << "Zone: qA " << i << " q = " << q
				<< " " << d << " " << average_controlsd << " " << Pfor[i] << "\n";}
	} else {
	  // Probability offset to make P(dmax) = 0
	  Pfor[i] = Max(0.0, TruncatedLorentzianProb(d, 0.0, controlsd[i], 0.0, dmax) - Poffset);
	  if (DEBUG) {std::cout << "Zone: qB " << i << " q = " << q
				<< " " << d << " " << controlsd[i] << " " << Pfor[i] << "\n";}
	}
      } else {Pfor[i] = 1.0;}
      if (validpoint[i]) {	
	Ptot += Pfor[i];
      }
    }  // end loop q
    if (Ptot > 0.0) {
      for (int i=0;i<npoint;i++) {  // Loop O(i) screw or glide translation
	Pfor[i] /= Ptot;
      }
    }
    // Store single value for probability if appropriate
    if (prob_yes < 0.0) {
      if (singleprob) {
	prob_yes = Pfor.back();
      } else {
	prob_yes = -1.0;
      }
    }

    results = true;
  }
  //--------------------------------------------------------------
  // Store probability
  void Zone::StoreProb(const double& ProbYes)
  {
    if (!singleprob)
      clipper::Message::message(clipper::Message_fatal
				("Zone::StoreProb called for case of no single value "));
    prob_yes = ProbYes;
  }
  //--------------------------------------------------------------
  // Return probability 
  double Zone::Prob() const
  {
    if (!singleprob)
      clipper::Message::message(clipper::Message_fatal
				("Zone::Prob called for case of no single value "));
    CalcResults();
    return prob_yes;
  }
  //--------------------------------------------------------------
  int Zone::Nobs() const 
  // Return number of contributions
  {
    if (axis) {
      return IndxIsigI.size();
    }
    return fsum.Nobs();
  }
  //--------------------------------------------------------------
  // Return probability at each grid point
  std::vector<double> Zone::p() const
  {
    CalcResults();
    return Pfor;
  }
  //--------------------------------------------------------------
  // Return list of indices
  std::vector<int> Zone::Indices() const
  {
    std::vector<int> IndxList;
    TestValidIndices();
    if (! valid) return IndxList;
    if (axis)
      {
	if (IndxIsigI.size() > 1)
	  {
	    for (size_t i=0;i<IndxIsigI.size();i++)
	      {
		IndxList.push_back(IndxIsigI[i].index);
	      }
	  }
      }
    else
      {
	if (maxIndx > minIndx)
	  // only if more than one index
	  {
	    for (int i=minIndx;i<=maxIndx;i++)
	      {
		IndxList.push_back(i);
	      }
	  }
      }
    return IndxList;
  }
  //--------------------------------------------------------------
  void Zone::dump() const
  {
    std::cout << "\n" << formatNewFrame(ReindexOp(), 1) << "\n"
	      << "  Order " << order << "  Npoint " << npoint
	      << "\n"
	      << "  Zone  " << direction
	      << "     Condition: " << FormatConditionNewFrame(ReindexOp(), 0) << "\n";

    std::cout << permute.format() << "\n";
    std::cout << "  Cond  ";
    for (int i=0;i<3;i++) std::cout << " " << cond[i];
    std::cout << "  QCond  ";
    for (int i=0;i<3;i++) std::cout << " " << qcond[i];
    std::cout << "  Validpoint: ";
    for (int i=0;i<npoint;i++) std::cout << " " << validpoint[i];
    std::cout << "\n";
    std::cout << "Reindex matrix: \n" << reindexmat.format() << "\n\n";
    std::cout << "Permute Index matrix: \n" << permuteIndex.format() << "\n\n";
    CalcResults();
    if (singleprob)
      {
	std::cout << "Single probability: " << prob_yes << "\n";
      }
    std::cout << "Probability values: ";
    for (int i=1;i<npoint;i++) {std::cout << " " << Pfor[i];}
    std::cout << "\n";
  }
  //--------------------------------------------------------------
  // Return two test hkl which will be systematically absent
  // if in this zone, for glide planes
  // RefToNew is reindex operator from "reference" (lattice) frame
  // (stored with StoreReindex function) to some new frame in
  // which to return the indices
  std::vector<Hkl> Zone::GlideTestHkl(const ReindexOp& RefToNew) const
  {
    ASSERT (!axis);  // Glides only
    std::vector<Hkl> Chkl;
    Hkl hkl(0,0,0);
    // Values for index
    int odd1 =  7;
    int odd2 = 11;
    int even =  2;
    //    int idx = odd1;

    if (!diagonal)
      // Principle zone glide
      {
	if (glide == "d")
	  // d(a),(b),(c)
	  //    reflections eg u u' 0, 2g 2u 0  (for hk0, d(c))
	  {
	    hkl[2] = 0;
	    hkl[0] = odd1;
	    hkl[1] = odd2;
	    Chkl.push_back(hkl);
	    hkl[0] = 2*odd1;
	    hkl[1] = 2*even;
	    Chkl.push_back(hkl);
	  }
	else if (glide == "n")
	  // n(a),(b),(c)
	  //    reflections eg  g u 0, u g 0 (for hk0, n(c))
	  {
	    hkl[2] = 0;
	    hkl[0] = odd1;
	    hkl[1] = even;
	    Chkl.push_back(hkl);
	    hkl[0] = even;
	    hkl[1] = odd1;
	    Chkl.push_back(hkl);
	  }
	else if (glide == "a" || glide == "b" || glide == "c")
	  // a, b, c glide
	  //    
	  //    reflections u u 0, u g 0 (for hk0, b(c))
	  //                u u 0, g u 0 (for hk0, a(c))
	  {
	    hkl[0] = odd1;
	    hkl[1] = odd2;
	    hkl[2] = 0;
	    Chkl.push_back(hkl);
	    // Transform cond vector to internal standard frame
	    clipper::Vec3<double> v = clipper::Vec3<double>(cond) * permute;
	    // Which is the non-zero element?	    
	    int jc = -1;
	    for (int i=0;i<3;i++)
	      {if (Nint(v[i]) !=0) jc = i;}
	    ASSERT (jc<2);
	    hkl[jc] = odd1;        // a or b glide in internal frame (hk0)
	    hkl[1-jc] = even;      // jc = 0 or 1
	    Chkl.push_back(hkl);
	  }
	else
	  clipper::Message::message(clipper::Message_fatal
				    ("Zone::GlideTestHkl: shouldn't happen"));
      }
    else
      // 110 etc  treat all as 110, hhl zone
      {
	if (glide == "d")
	  // d(110) 
	  //    reflections eg u u u', 2u 2u 2u'
	  {
	    hkl[0] = odd1;
	    hkl[1] = odd1;
	    hkl[2] = odd2;
	    Chkl.push_back(hkl);
	    hkl[0] = 2*odd1;
	    hkl[1] = 2*odd1;
	    hkl[2] = 2*odd2;
	    Chkl.push_back(hkl);
	  }
	else if (glide == "a" || glide == "b" || glide == "c" || glide == "n")
	  // c or n(110) 
	  //    reflections eg u u u', g g u
	  {
	    hkl[0] = odd1;
	    hkl[1] = odd1;
	    hkl[2] = odd2;
	    Chkl.push_back(hkl);
	    hkl[0] = even;
	    hkl[1] = even;
	    hkl[2] = odd2;
	    Chkl.push_back(hkl);
	  }
	else
	  clipper::Message::message(clipper::Message_fatal
				    ("Zone::GlideTestHkl: shouldn't happen"));
      }

    // Unpermute to original form, convert to "New" frame
    //  using RefToNew operator [Hr]]
    // generated hkl is h'
    // then h(new) = h' [P]^-1 [H]^-1 [Hr] = [HP]^-1 [Hr]
    for (size_t i=0;i<Chkl.size();i++)
      {
	Chkl[i] = Chkl[i].change_basis(permuteIndex.inverse() * RefToNew);
      }

    return Chkl;
  }
  //--------------------------------------------------------------
  std::vector<int> PairIntVec(const int& i1, const int& i2)
  // pair of integers as vector
  {
    std::vector<int> pint(2);
    pint[0] = i1;
    pint[1] = i2;
    return pint;
  }
  //--------------------------------------------------------------
  // Return test hkl's for systematically absent test
  // test hkl are index 2 & 3 (* cond)
  // Also return pair of flags for each screw component,
  // for index 2 & 3 respectively
  // If these flags match systematic absence condition
  // (+1 absent, 0 present) then that screw is present
  // 
  // for principal axes only, not 110)

  // RefToNew is reindex operator from "reference" (lattice) frame
  // (stored with StoreReindex function) to some new frame in
  // which to return the indices

  std::vector<Hkl> Zone::AxisTestHkl(const ReindexOp& RefToNew,
				     std::vector<std::vector<int> >& screwabsence) const
  {
    ASSERT (axis);  // Axes only
    std::vector<Hkl> Chkl;

    if (direction == "a" || direction == "b" || direction == "c")
      {
	// cond is 100, 010 or 001
	Hkl hkl;
	int idx = 2;
	for (int i=0;i<3;i++) hkl[i] = idx * cond[i];
	Chkl.push_back(hkl);
	idx = 3;
	for (int i=0;i<3;i++) hkl[i] = idx * cond[i];
	Chkl.push_back(hkl);

	screwabsence.resize(npoint);
	screwabsence[0] = PairIntVec(0,0);  // first point is always present
	// Just tabulate everything (not clever)
	if (order == 2)
	  {
	    screwabsence[1] = PairIntVec(0,+1);  // 2(1)
	  }
	else if (order == 3)
	  {
	    screwabsence[1] = PairIntVec(+1,0);  // 3(1)
	  }
	else if (order == 4)
	  {
	    if (validpoint[1])
	      screwabsence[1] = PairIntVec(0,+1);  // 4(2) only if valid
	    else
	      screwabsence[1] = PairIntVec(-1,-1);   // 4(2) only if valid
	    screwabsence[2] = PairIntVec(+1,+1);   // 4(1)
	  }
	else if (order == 6)
	  {
	    screwabsence[1] = PairIntVec(0,+1);  // 6(3)
	    screwabsence[2] = PairIntVec(+1,0);  // 6(2)
	    screwabsence[3] = PairIntVec(+1,+1);   // 6(1)
	  }
      }

    // convert to "New" frame
    //  using RefToNew operator [Hr]]
    // generated hkl is h
    // then h(new) = h [H]^-1 [Hr]
    for (size_t i=0;i<Chkl.size();i++)
      {
	Chkl[i] = Chkl[i].change_basis(reindexmat.inverse() * RefToNew);
      }

    return Chkl;
  }
  //--------------------------------------------------------------
  bool Zone::Compare(const Zone& other) const
  {
    if (axis != other.axis) return false;
    if (order != other.order) return false;
    // Use reindex operator to check that both cond is the same
    // ignoring sign changes
    bool same = true;
    for (int i=0;i<3;i++)
      {
	if (std::abs(qcond[i]) != std::abs(other.qcond[i]))
	  same = false;
      }
    if (same)
      // Check direction
      {
	clipper::Vec3<double> vd1 = DirectionRefFrame();
	clipper::Vec3<double> vd2 = other.DirectionRefFrame();
	for (int i=0;i<3;i++)
	  {if (std::abs(vd1[i]) != std::abs(vd2[i]))
	      {same = false;}
	  }
      }
    return same;
  }
  //--------------------------------------------------------------
  bool SysAbsScore::IsZoneInGroup(const int& iz) const
  // Return true if zone iz is in group
  {
    return (std::find(zonesingroup.begin(), zonesingroup.end(), iz) !=
	    zonesingroup.end());
  }
  //--------------------------------------------------------------

}
