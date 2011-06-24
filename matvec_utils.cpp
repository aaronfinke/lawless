// matvec_utils.cpp
// Utility routines

#include "matvec_utils.hh"
#include "util.hh"
#include "scala_util.hh"
#include "string_util.hh"

namespace MVutil
{
  // * *  * *  * *  * *  * *  * *  * *  * *  * *  * *  * *  * * 
  // Setting clipper matrices
  //--------------------------------------------------------------
  clipper::Mat33<double> SetCMat33(const float m[3][3])
  //                     ^^^^^^^
  {
    return clipper::Mat33<double>(m[0][0], m[0][1], m[0][2],
				  m[1][0], m[1][1], m[1][2],
				  m[2][0], m[2][1], m[2][2]);
  }
  //--------------------------------------------------------------
  clipper::Mat33<double> SetCMat33(const float m[9])
  //                     ^^^^^^^
  {
    return clipper::Mat33<double>(m[0], m[3], m[6],
				  m[1], m[4], m[7],
				  m[2], m[5], m[8]);
  }
  //--------------------------------------------------------------
  clipper::Mat33<double> SetCMat33(const double m[3][3])
  //                     ^^^^^^^
  {
    return clipper::Mat33<double>(m[0][0], m[0][1], m[0][2],
				  m[1][0], m[1][1], m[1][2],
				  m[2][0], m[2][1], m[2][2]);
  }
  //--------------------------------------------------------------
  clipper::Mat33<double> SetCMat33(const std::vector<double>& m)
    //                 ^^^^^^^
  {
    return clipper::Mat33<double>(m.at(0), m.at(1), m.at(2),
 				   m.at(3), m.at(4), m.at(5),
				   m.at(6), m.at(7), m.at(8));

  }
  //--------------------------------------------------------------
  clipper::Mat33<double> SetCMat33(CSym::ccp4_symop& op)
    //                 ^^^^^^^
  {
    return clipper::Mat33<double>(op.rot[0][0], op.rot[0][1], op.rot[0][2],
				  op.rot[1][0], op.rot[1][1], op.rot[1][2],
				  op.rot[2][0], op.rot[2][1], op.rot[2][2]);
  }
  //--------------------------------------------------------------
  // Diagonal matrix
  clipper::Mat33<double> SetDiagCMat33(const clipper::Vec3<double>& v)
  {
    return clipper::Mat33<double>(v[0],0.,0.,
				  0.,v[1],0.,
				  0.,0.,v[2]);
  }
  // * *  * *  * *  * *  * *  * *  * *  * *  * *  * *  * *  * * 
  // Setting Vector
  //--------------------------------------------------------------
  std::vector<double> SetVMat33(const clipper::Mat33<double>& R)
  // Convert Mat33 to vector type
  {
    std::vector<double> Mv(9);
    int k=0;
    for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
	Mv[k++] = R(i,j);
    return Mv;
  }
  //--------------------------------------------------------------
  std::vector<double> SetVMat33(const float m[3][3])
  // Convert 3x3 array to vector type
  {
    std::vector<double> Mv(9);
    int k=0;
    for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
	Mv[k++] = m[i][j];
    return Mv;
  }
  // Tests
  //--------------------------------------------------------------
  bool is_mat33_ident(const clipper::Mat33<float>& R)
    // Returns true if 3x3 matrix R is identity
  {
    return R.equals(Mat33<float>::identity(), 1.0e-4);
  }

  //--------------------------------------------------------------
  bool is_mat33_ident(const clipper::Mat33<double>& R)
    // Returns true if 3x3 matrix R is identity
  {
    return R.equals(Mat33<double>::identity(), 1.0e-4);
  }
  //--------------------------------------------------------------
  bool is_rtop_ident(const clipper::RTop<double>& R)
    // Returns true if 3x3 matrix + translation R is identity
  {
    return R.equals(RTop<double>::identity(), 1.0e-4);
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  std::string FormatSymop_as_hkl(const clipper::Symop& op,
				 const std::string& brackets)
  {
    //^    std::cout <<"FormatSymop_as_hkl\n" << op.format() << "\n"; //^
    // transpose, convert to symop and format
    std::string s = clipper::Symop(clipper::RTop<>
		 (clipper::RTop<>(op).rot().transpose())).format();
    // replace xyz by hkl
    size_t p;
    while ((p = s.find("x")) != std::string::npos) {
      s.replace(p,1,"h");}
    while ((p = s.find("y")) != std::string::npos) {
      s.replace(p,1,"k");}
    while ((p = s.find("z")) != std::string::npos) {
      s.replace(p,1,"l");}
    return StringUtil::Strip(brackets[0]+s+brackets[1]);
  }
  //--------------------------------------------------------------
  std::string FormatReindex_as_hkl(const clipper::RTop<double>& op,
				   const std::string& brackets)
  // a reindex operator H post-multiplies a row vector h
  // h'T = hT [H]
  {
    //^    std::cout <<"FormatReindex_as_hkl\n" << op.format() << "\n"; //^
    char chkl[] = {'h','k','l'};
    std::string s;
    int width = 7; // fall-back field width
    const double TOL = 0.00001;
    for (int i=0;i<3;++i) { // loop columns
      bool first = true;
      for (int j=0;j<3;++j) { // loop rows
	double r = op.rot()(j,i);
	if (std::abs(r) > TOL) { // not zero
	  if (std::abs(r+1.0) < TOL) { // -1
	    s += "-";
	  } else {
	    if (!first && r > 0.0) s += "+";
	    if (std::abs(r-1.0) > TOL) { // not +1.0
	      s += StringUtil::formatFraction(op.rot()(j,i), width);
	    }
	  }
	  s += chkl[j];
	  first = false;
	}
      }
      if (i<2) s += ",";
    }
    return StringUtil::Strip(brackets[0]+s+brackets[1]);
  }
  // * *  * *  * *  * *  * *  * *  * *  * *  * *  * *  * *  * * 
  // Other stuff
  //--------------------------------------------------------------
clipper::Vec3<int> IntVec(const clipper::Vec3<double>& vector)
  // make integral version of vector by scaling up, within tolerance
  // leave argument vector unchanged
  {
    const double tolerance = 0.05;
    double dmin = 1.0;
    double d;

    // Unit version
    clipper::Vec3<double> v = vector.unit();

    // Find smallest value greater than tolerance
    for (int i=0;i<3;i++)
      {
	d = std::abs(v[i]);
	if (d > tolerance && d < (1.0-tolerance))
	  {
	    dmin = Min(d, dmin);
	  }
      }
    // multiplier is inverse of minimum value
    v = (1./dmin) * v;
    return clipper::Vec3<int>(Nint(v[0]),Nint(v[1]), Nint(v[2]));
  }
  // * *  * *  * *  * *  * *  * *  * *  * *  * *  * *  * *  * * 
  //--------------------------------------------------------------
  std::vector<double> Eigen33(const DMat33& M, DMat33& E)
  // Returns eigenvalues & eigenvectors (columns in E) of real symmetric matrix
  {
    clipper::Matrix<double> SM(3,3);
    for (int i=0;i<3;++i)
      for (int j=0;j<3;++j)
	SM(i,j) = M(i,j);
    std::vector<double> EV = SM.eigen(true);  // sorted with smallest 1st
    E = DMat33(SM(0,0),SM(0,1),SM(0,2),
	     SM(1,0),SM(1,1),SM(1,2),
	     SM(2,0),SM(2,1),SM(2,2));
    return EV;
  }
  //--------------------------------------------------------------
  DMat33 LsqTransform(const std::vector<DVect3>& x, const std::vector<DVect3>& y,
		      double& Resid, const bool& AvoidPlanar)
  {
    // if the 3x3 matrix [U] transforms x -> y, with a possible error e
    // then e = [U] x - y
    //
    //  Return [U] & mean residual Resid
    //
    // Minimising E = 1/2N Sum(e.e)  ie least-squares gives
    //   [U] = [R] [S]^-1
    // where
    //   [S]ij  =  Sum(xi xj)    i, j = 1,3 for 3-vector
    //   [R]ij  =  Sum(yi xj)
    //
    // On entry:
    //  x, y    corresponding vector lists
    //  AvoidPlanar  if true, generate a few orthogonal vectors if needed to avoid the ambiguities
    //          introduced in all x vectors are coplanar (not included in residual E)
    //          This should only arise in an orthogonal coordinate frame
    //          when it's relevant
    //   
    ASSERT (x.size() == y.size());
    DMat33 S = SetCMat33(std::vector<double>(9,0.0));  // clear matrices
    DMat33 R = SetCMat33(std::vector<double>(9,0.0));

    const double THRESHOLD = 1.0e-3;       

    for (size_t n=0;n<x.size();++n) {
      for (int i=0;i<3;++i) {
	for (int j=0;j<3;++j) {
	  S(i,j) += x[n][i] * x[n][j];
	  R(i,j) += y[n][i] * x[n][j];
	}}
    } // end loop observations
    std::cout << "LsqT det[S] " << S.det() << "\n";
    std::cout << "LsqT det[R] " << R.det() << "\n";
    
    if (AvoidPlanar) {
      if (S.det() < THRESHOLD || R.det() < THRESHOLD) {
	// Add in a few cross-product vectors
	const int KSAMPLE = 20;   // 5%
	int nadd = Max(1,int(x.size())/KSAMPLE);
	std::cout << "Adding cross-vectors " << nadd <<"\n";
	int ninc = x.size()/(2*nadd);  // increment for index
	int n = 0; // 1st of pair
	int m = x.size()-1; // 2nd of pair
	
	for (int k=0;k<nadd;++k) {
	  //  a pair of corresponding vector orthogonal to one x,y pair
	  DVect3 xc = DVect3::cross(x[n], x[m]);
	  DVect3 yc = DVect3::cross(y[n], y[m]);
	  for (int i=0;i<3;++i) {  // add in to sums
	    for (int j=0;j<3;++j) {
	      S(i,j) += xc[i] * xc[j];
	      R(i,j) += yc[i] * xc[j];
	    }}
	  n += ninc;
	  m -= ninc; 
	}
      }
    }

    DMat33 U = R * S.inverse();
    
    Resid = 0.0; // residual 1/2N Sum(e.e)
    for (size_t n=0;n<x.size();++n) {
      DVect3 e = U * x[n]  - y[n];
      Resid += e * e;
    }
    Resid /= double(x.size());
    
    return U;
  }
  // * *  * *  * *  * *  * *  * *  * *  * *  * *  * *  * *  * * 
  //--------------------------------------------------------------
  clipper::Vec3<double> AxisDirection(const clipper::Mat33<double>& R)
  {
    // Get axis direction for rotation matrix R
    // Returns normalised axis direction
    //
    // Algorithm (from Kevin Cowtan)
    // If we rotate three non-colinear points (eg axis vector
    // 100, 010, 001) then the axis is perpendicular to the two
    // largest difference vectors,
    // ie to the cross-product with largest magnitude
    //
    // For matrix R, get cross-product of all columns of (R - I)
    // choose the one with the largest magnitude, & normalise
    //
    // This works even for rotations in non-orthogonal space, as long
    // as they are crystallographic (ie non-orthogonality is 
    // preserved by this operator)
    
    // The identity matrix
    clipper::Mat33<double> Imat;
    Imat = clipper::Mat33<double>::identity();
    
    // R - I
    clipper::Mat33<double> RmI;
    for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
	RmI(i,j) = R(i,j) - Imat(i,j);
    
    // Extract columns
    std::vector<clipper::Vec3<double> > cols(3);
    
    for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
	(cols[i])[j] = RmI(j,i);
    
    // Cross products for all 3 pairs
    std::vector<clipper::Vec3<double> > cp(3);
    std::vector<double> mag(3);
    
    double maxv = -1.0;
    int maxc = 0;
    
    for (int i=0;i<3;i++)
      {
	int j = (i+1)%3;
	int k = (j+1)%3;
	cp[i] = clipper::Vec3<double>::cross(cols[j], cols[k]);
	mag[i] = MVutil::Modulus(cp[i]);
	if (mag[i] > maxv)
	  {
	    // largest cross-product
	    maxv = mag[i];
	    maxc = i;
	  }
      }
    
    return cp[maxc].unit();
  }
  //--------------------------------------------------------------
}
