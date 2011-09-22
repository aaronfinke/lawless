// \file hkl_datatypes.cpp

#include <vector>
#include <algorithm>

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;

#include "csymlib.h"    // CCP4 symmetry stuff

#include "rotation.hh"
#include "hkl_datatypes.hh"
#include "scala_util.hh"
#include "range.hh"
#include "matvec_utils.hh"
#include "string_util.hh"
#include "hkl_symmetry.hh"

namespace scala
{
  //--------------------------------------------------------------
  //------------------------------------------------------------
  std::string Chiral_as_string(const Chirality& chiral)
  {
    if (chiral == UNKNOWN) return "UNKNOWN";
    if (chiral == CHIRAL) return "CHIRAL";
    if (chiral == NONCHIRAL) return "NONCHIRAL";
    if (chiral == CENTROSYMMETRIC) return "CENTROSYMMETRIC";
    return "";
  }
  //--------------------------------------------------------------
  // Orientation things
  //
  DVect3 SpindleToPrincipleAxes(const DMat33& DUB)
  // Returns angles in degrees between principle reciprocal axes a*, b*, c*
  // and the rotation spindle (along z). Angles in range 0 - 90 degress
  // On entry:
  //  DUB is setting matrix [D][U][B]
  {
    DVect3 angles;
    DVect3 spindle(0.0,0.0,1.0);  // spindle along z
    for (int i=0;i<3;++i) {
      DVect3 axis(0.0,0.0,0.0);
      axis[i] = 1.0;  // a*, b*, c*
      angles[i] = clipper::Util::rad2d(acos
	    (std::abs(clipper::Vec3<>::dot(spindle, (DUB * axis).unit()))));
    }
    return angles;
  }
  //--------------------------------------------------------------
  int SmallestComponent(const DVect3& v)
  // Returns index of smallest component of vector
  {
    int ismall = 0; // EK: renamed due to Windows clash
    double min = +1.0e30;
    for (int i=0;i<3;++i) {
      if (v[i] < min) {
        ismall = i;
	min = v[i];
      }
    }
    return ismall;
  }
  //--------------------------------------------------------------
  std::string FormatSpindleToPrincipleAxes(const DMat33& DUB)
  // Format closest reciprocal axis to spindle
  {
    DVect3 angles = SpindleToPrincipleAxes(DUB);
    int ismall = SmallestComponent(angles); // EK: renamed due to Windows clash
    std::string axis;
    if (ismall == 0) {axis = "a*";}
    else if (ismall == 1) {axis = "b*";}
    else if (ismall == 2) {axis = "c*";}
    return axis+" (angle "+clipper::String(angles[ismall],6,3)+
      " degrees)";
  }
  //--------------------------------------------------------------
  int ClosestPrincipleAxisToSpindle(const DMat33& DUB)
  // 0,1,2 for h,k,l
  {
    DVect3 angles = SpindleToPrincipleAxes(DUB);
    return SmallestComponent(angles);
  }
  //--------------------------------------------------------------
  ReindexOp::ReindexOp(const std::string& Operator)
  // Construct from string eg "2h+k,k,l"
  {
    float rt44[10][4][4];  // Allocate excess space in case of error
    int Noper = CSym::symfr_driver(Operator.c_str(), rt44);
    if (Noper != 1)
      Message::message(Message_fatal
		       ("ReindexOp: syntax error in operator; "+Operator));

    Mat33<double> H;
    Vec3<double> v;
    for (int i=0;i<3;i++)
      {
	for (int j=0;j<3;j++)
	  // We want the transpose of the matrix since reindex operator H 
	  // applies to index h such that h'T = hT H
	  {H(j,i) = rt44[0][i][j];}
	v[i] = rt44[0][i][3];
      }
    rot() = H;
    trn() = v;
  }
  //--------------------------------------------------------------
  std::string ReindexOp::as_hkl() const
  {
    std::string st =  MVutil::FormatReindex_as_hkl(*this);
    return st;
  }
  //--------------------------------------------------------------
  std::string ReindexOp::as_hkl_XML() const
  {
    std::string st =  MVutil::FormatReindex_as_hkl(*this);
    return "<ReindexOperator>"+st+"</ReindexOperator>";
  }
  //--------------------------------------------------------------
  std::string ReindexOp::as_matrix() const
  // Rotation only!
  {
    std::string fmat = "   h'   = ( h k l ) (";
    for (int i=0;i<3;i++)
      {
	if (i > 0) fmat += "                    (";
	for (int j=0;j<3;j++)
	  {
	    fmat += " "+String(rot()(i,j),7,4);
	  }
	fmat += " )\n";
      }
    return fmat;
  }
  //--------------------------------------------------------------
  std::string ReindexOp::as_XML() const
  {
    std::string fmat = "<ReindexMatrix>";
    for (int i=0;i<3;i++)
      {
	if (i > 0) fmat += "                  ";
	for (int j=0;j<3;j++)
	  {
	    fmat += " "+String(rot()(i,j),6,4);
	  }
	fmat += "\n";
      }
    fmat += "    </ReindexMatrix>\n";
    return fmat;
  }
  //--------------------------------------------------------------
  bool operator < (const ReindexOp& a,const ReindexOp& b)
  // Sort on deviation, putting identity before anything else
  {
    if (MVutil::is_rtop_ident(a)) return true;
    if (MVutil::is_rtop_ident(b)) return false;
    return a.deviation < b.deviation;
  }
  //--------------------------------------------------------------
  bool ReindexOp::equals(const ReindexOp& b, const double& tol) const
  {
    return (rot().equals(b.rot(), tol) && trn().equals(b.trn(), tol));
  }
  //--------------------------------------------------------------
  void ReindexOp::FindSimplest(const SpaceGroup& symm)
  // Apply symmetry operators to reindex operator H and choose the simplest
  // Criteria for "simplest"
  // 1) identity
  // 2) smallest sum of absolute value of all elements
  // 3) minimum number of negatives
  //  cf. AlternativeBases::SimplestOp
  //
  //   h'T = hT [H]  reindex      hiT = hT [Ri] symmetry
  //   [H'] = [Ri] [H]   where [Ri] is the reciprocal space symmetry operator
  {
    Mat33<double> H = rot();  // initial reindex
    Mat33<double> Hbest = H;  // initial reindex
    Mat33<double> Imat =  Mat33<double>::identity();  // identity matrix
    const double TOL = 0.001;  // tolerance

    double besttot = 10000000.; // sum of absolute value of all elements
    int minneg = 100000; // minimum number of negatives
    int kbest = -1;    

    for (int isym=0;isym<symm.num_primops();++isym) {
      Mat33<double> Ht = H * symm.symop(isym).rot().inverse();
      if (Ht.equals(Imat, TOL)) {
	rot() = Imat;
	//^
	//	std::cout << "*** FindSimplest: identity with symop " << isym << "\n"
	//		  << symm.symop(isym).rot().inverse().format() << "\n"; //^-
	return;  // always accept identity operation if found
      }
      double total = 0.0; // total
      int nneg = 0; // number of negatives
      for (int i=0;i<3;i++)
	for (int j=0;j<3;j++) {
	  total += std::abs(Ht(i,j));
	  if (Ht(i,j) < 0.0) nneg += 1;
	}
      if (total < besttot-TOL) {
	besttot = total;
	minneg = nneg;
	Hbest = Ht;
	kbest = isym;
      } else if (Close<float>(total, besttot, TOL)) {
	if (nneg < minneg) {
	  besttot = total;
	  minneg = nneg;
	  Hbest = Ht;
	  kbest = isym;
	}
      }
      //^
      //      std::cout << "\n\nFindSimplest: symop " << isym << "\n"
      //		  << symm.symop(isym).rot().inverse().format()
      //		<< "\n\nReindex:\n" << Ht.format()
      //		<< " total " << total << " nneg " << nneg
      //		<< "\n"; 
      //^-
    }
    rot() = Hbest;
    //^
    //    std::cout << "\n***FindSimplest: best reindex " << as_hkl() << " with symop " << kbest << "\n"
    //	      << symm.symop(kbest).rot().inverse().format() << "\n\n"; //^-
  }
  //--------------------------------------------------------------
  //! change basis of symmetry operator
  clipper::Symop ReindexOp::Symop(const clipper::Symop& symop) const
  {
    // S' = H^-1 S H
    // t' = H^-1 t
    clipper::Mat33<double> HR = rot();
    clipper::RTop<double> Sp(HR.inverse()*symop.rot()*HR, HR.inverse() * symop.trn());
    return clipper::Symop(Sp);
  }
  //--------------------------------------------------------------
  bool operator == (const ReindexOp& a,const ReindexOp& b)
  {
    return a.equals(b);
  }
  //--------------------------------------------------------------
  bool operator != (const ReindexOp& a,const ReindexOp& b)
  {
    return !(a.equals(b));
  }
  //--------------------------------------------------------------
   ReindexOp operator * (const ReindexOp& a,const ReindexOp& b)
   {
     return ReindexOp(RTop<double>(a) * RTop<double>(b));
   }
  //--------------------------------------------------------------
   ReindexOp operator * (const RTop<double>& a,const ReindexOp& b)
   {
     return ReindexOp(a * b);
   }
  //--------------------------------------------------------------
   ReindexOp operator * (const ReindexOp& a,const RTop<double>& b)
   {
     return ReindexOp(a * b);
   }
  //--------------------------------------------------------------
  Vec3<double> operator * (const Vec3<double> v, const ReindexOp& R)
  // vT . R
  {
    return v*R.rot() + R.trn();
  }

  //--------------------------------------------------------------
  //--------------------------------------------------------------
  /*! Construct and initialise a metric tensor, given a set of real or
    reciprocal cell parameters.
    \param a Length of \b a axis in Angstroms or reciprocal Angstroms.
    \param b Length of \b b axis in Angstroms or reciprocal Angstroms.
    \param c Length of \b c axis in Angstroms or reciprocal Angstroms.
    \param alph Angle between \b b and \b c in degrees.
    \param beta Angle between \b a and \b c in degrees.
    \param gamm Angle between \b a and \b b in degrees.
  */
  // from Clipper with extensions
  //  Note that the off-diagonal elements are stored doubled so that they can be used
  //  directly to multiply vectors
  MetricTensor::MetricTensor( const ftype& a, const ftype& b, const ftype& c, const ftype& alph, const ftype& beta, const ftype& gamm )
  {
    double degtorad = atan(1.0)/45.0;
    m00 = a*a;
    m11 = b*b;
    m22 = c*c;
    m01 = 2.0*a*b*cos(gamm*degtorad);
    m02 = 2.0*a*c*cos(beta*degtorad);
    m12 = 2.0*b*c*cos(alph*degtorad);
  }
  //--------------------------------------------------------------
  MetricTensor::MetricTensor(const std::vector<Dtype> cell)
  {
    double degtorad = atan(1.0)/45.0;
    m00 = cell[0]*cell[0];
    m11 = cell[1]*cell[1];
    m22 = cell[2]*cell[2];
    m01 = 2.0*cell[0]*cell[1]*cos(cell[5]*degtorad);
    m02 = 2.0*cell[0]*cell[2]*cos(cell[4]*degtorad);
    m12 = 2.0*cell[1]*cell[2]*cos(cell[3]*degtorad);
  }
  //--------------------------------------------------------------
  MetricTensor::MetricTensor(const Mat33<Dtype> MetricMatrix)
  // Constructor from full matrix - no check on equality of
  // off-diagonal elements, just add them together to store
  // 2 * a.b etc
  {
    m00 = MetricMatrix(0,0);
    m11 = MetricMatrix(1,1);
    m22 = MetricMatrix(2,2);
    m01 = MetricMatrix(0,1) + MetricMatrix(1,0);
    m02 = MetricMatrix(0,2) + MetricMatrix(2,0);
    m12 = MetricMatrix(1,2) + MetricMatrix(2,1);
  }
  //--------------------------------------------------------------
  std::vector<Dtype> MetricTensor::cell() const
  // Return unit cell dimensions (real or reciprocal)
  // Angles in degrees
  {
    double degtorad = atan(1.0)/45.0;
    std::vector<Dtype> cell(6);
    cell[0] = sqrt(m00);
    cell[1] = sqrt(m11);
    cell[2] = sqrt(m22);
    cell[3] = acos(0.5 * m12 / (cell[1] * cell[2]))/degtorad;
    cell[4] = acos(0.5 * m02 / (cell[2] * cell[0]))/degtorad;
    cell[5] = acos(0.5 * m01 / (cell[0] * cell[1]))/degtorad;
    return cell;
  }

  //--------------------------------------------------------------
  MetricTensor MetricTensor::inverse() const
  {
    //      Mat33<Dtype> MetricMatrix = matrix();
    return MetricTensor((*this).matrix().inverse());
  }
  //--------------------------------------------------------------
  Mat33<Dtype> MetricTensor::matrix() const
  {
    return Mat33<Dtype>(m00,     0.5*m01, 0.5*m02,
			0.5*m01, m11,     0.5*m12,
			0.5*m02, 0.5*m12,  m22);
  }

  //--------------------------------------------------------------
  String MetricTensor::format() const
  {
    return "m00=" + String(m00) + " m11=" + String(m11) + " m22=" + String(m22) + " m01=" + String(m01) + " m02=" + String(m02) + " m12=" + String(m12);
  }


  //--------------------------------------------------------------
  Scell::Scell(const std::vector<Dtype>& real_cell)
  // Construct Scell object from vector cell
  {
    init(real_cell);
  }

  //--------------------------------------------------------------
  Scell::Scell(const float* cell)
  // Construct Scell object from float array cell
  {
    std::vector<Dtype> real_cell(cell, &cell[6]);
    init(real_cell);
  }
  //--------------------------------------------------------------
  Scell::Scell(const clipper::Cell& ccell)
  // Construct Scell object from clipper::Cell
  {
    std::vector<Dtype> real_cell(6);
    real_cell[0] = ccell.a();
    real_cell[1] = ccell.b();
    real_cell[2] = ccell.c();
    real_cell[3] = ccell.alpha_deg();
    real_cell[4] = ccell.beta_deg();
    real_cell[5] = ccell.gamma_deg();
    init(real_cell);
  }
  //--------------------------------------------------------------
  Scell::Scell(const double& a,const double& b,const double& c,
	       const double& alpha, const double& beta, const double& gamma)
  // Construct Scell object from real cell
  {
    std::vector<Dtype> real_cell(6);
    real_cell[0] = a;
    real_cell[1] = b;
    real_cell[2] = c;
    real_cell[3] = alpha;
    real_cell[4] = beta;
    real_cell[5] = gamma;
    init(real_cell);
  }
  //--------------------------------------------------------------
  Scell::Scell(const MetricTensor& metric_tensor, const bool real)
  // Construct from real or reciprocal metric tensor
  {
    if (real) {
      // real cell
      init(metric_tensor.cell());
    } else {
      // reciprocal cell
      init(metric_tensor.inverse().cell());
    }
  }
  //--------------------------------------------------------------
  // Construct cell class from real cell
  //  derive reciprocal cell and orthogonalisation matrix B
  void Scell::init(const std::vector<Dtype>& real_cell)
  {

    double degtorad = atan(1.0)/45.0;

    cell_.resize(6);
    rec_cell_.resize(6);

    for (int i=0;i<6;i++) cell_[i] = real_cell[i];

    // Real metric tensor
    metric_tensor_ = MetricTensor(cell_);
    // Reciprocal metric tensor
    recip_metric_tensor_ = metric_tensor_.inverse();
    // Reciprocal cell
    rec_cell_ = recip_metric_tensor_.cell();

    // Orthogonalisation matrix (A**-1)
    // Reciprocal cell
    double as = rec_cell_[0];  // a*
    double bs = rec_cell_[1];  // b*
    double cs = rec_cell_[2];  // c*
    //    double cas = cos(rec_cell_[3]*degtorad);  // cos(alpha*)
    double cbs = cos(rec_cell_[4]*degtorad);  // cos(beta*)
    double cgs = cos(rec_cell_[5]*degtorad);  // cos(gamma*)
    //    double sas = sin(rec_cell_[3]*degtorad);          // sin(alpha*)
    double sbs = sin(rec_cell_[4]*degtorad);          // sin(beta*)
    double sgs = sin(rec_cell_[5]*degtorad);          // sin(gamma*)

    double cc = cell_[2];
    double ca = cos(cell_[3]*degtorad);

    Bmat_ = Mat33<Dtype>(as, bs*cgs, cs*cbs,
			 0., bs*sgs, -cs*sbs*ca,
			 0.,     0., 1.0/cc);
  }
  //--------------------------------------------------------------
  std::string Scell::formatPrint(const bool& newline) const
  {
    std::string s;
    for (int i = 0; i < 6; ++i)
      s += FormatOutput::logTabPrintf(0, "%7.2f ", cell_[i]);
    if (newline) s += "\n";
    return s;
  }
  //--------------------------------------------------------------
  void Scell::dump() const
  {
    std::cout << "\n\n      Real cell: ";
    for (int i = 0; i < 6; ++i)
      std::cout << cell_[i] << "  ";
    std::cout << std::endl;
    std::cout << "Reciprocal cell: ";
    for (int i = 0; i < 6; ++i)
      std::cout << rec_cell_[i] << "  ";
    std::cout << std::endl;

    std::cout << "RealMetricTensor:\n"
	      << metric_tensor_.format() << "\n";
    std::cout << "RecipMetricTensor:\n"
	      << recip_metric_tensor_.format() << "\n";

    std::cout << "Bmat:\n" 
	      << Bmat_.format() << "\n";

  }
  //--------------------------------------------------------------
  Scell Scell::change_basis(const ReindexOp& reindex_op) const
  {
    // H reindex operator (ignores translations)
    // B current orthogonalisation matrix
    // new (UB)' = B H(T)^-1
    Mat33<Dtype> UBp = Bmat_ * reindex_op.transpose().inverse();
    // Reciprocal metric tensor = (UB)'(T) (UB)'
    // construct new cell object from reciprocal metric tensor
    return Scell(UBp.transpose() * UBp, false);
  }
  //--------------------------------------------------------------
  std::string Scell::format(const int w, const int p) const
  {
    return FormatCell(cell_,w,p);
  }
  //--------------------------------------------------------------
  std::string Scell::xml() const
  {
    std::string line = "<cell>\n";
    line += "   <a>"+clipper::String(cell_[0],7,4)+"</a>\n";
    line += "   <b>"+clipper::String(cell_[1],7,4)+"</b>\n";
    line += "   <c>"+clipper::String(cell_[2],7,4)+"</c>\n";
    line += "   <alpha>"+clipper::String(cell_[3],7,4)+"</alpha>\n";
    line += "   <beta>"+clipper::String(cell_[4],7,4)+"</beta>\n";
    line += "   <gamma>"+clipper::String(cell_[5],7,4)+"</gamma>\n";
    line += "</cell>\n";
    return line;
  }
  //--------------------------------------------------------------
  double Scell::Volume() const
  {
    double degtorad = atan(1.0)/45.0;
    double a = cell_[0];
    double b = cell_[1];
    double c = cell_[2];
    double alpha = cell_[3] * degtorad;
    double beta  = cell_[4] * degtorad;
    double gamma = cell_[5] * degtorad;
    double vol = a * b * c *
      sqrt( 2.0*cos(alpha)*cos(beta)*cos(gamma)
	    - cos(alpha)*cos(alpha)
	    - cos( beta)*cos( beta)
	    - cos(gamma)*cos(gamma) + 1.0 );
    return vol;
  }
  //--------------------------------------------------------------
  // Return true if all angles are within Tol of given test values
  bool Scell::AngleTest(const double& AlphaTest,
			const double& BetaTest,
			const double& GammaTest,
			const double Tol) const
  {
    if (Close<double,double>(cell_[3], AlphaTest, Tol))
      if (Close<double,double>(cell_[4], BetaTest, Tol))
	if (Close<double,double>(cell_[5], GammaTest, Tol))
	  return true;
    return false;
  }
  //--------------------------------------------------------------
  // Returns true if cells agree within "tolerance"
  bool Scell::equals(const Scell& other, const double& tolA, const double& tolD) const
  {
    bool OK = true;
    for (int i=0;i<3;i++)
      {if (!Close(cell_[i], other.cell_[i], tolA)) OK = false;}
    for (int i=3;i<6;i++)
      {if (!Close(cell_[i], other.cell_[i], tolD)) OK = false;}
    return OK;
  }
  //--------------------------------------------------------------
  // Returns true if cells agree within "tolerance"
  bool Scell::equalsTol(const Scell& other, const double& AngularTolerance) const
  {
    // Very rough: convert angular tolerance to "length" using
    // average cell edge
    //	double edge = Max(Max(test_cell[0], test_cell[1]), test_cell[2]);
    double edge = 0.333 * (cell_[0] + cell_[1] + cell_[2]);
    double max_diff = clipper::Util::d2rad(AngularTolerance) * edge;
    return equals(other, max_diff, AngularTolerance);
  }
  //--------------------------------------------------------------
  bool Scell::null() const
  // Returns true if any value is close to zero
  {
    Dtype tol = 0.1;
    for (int i=0;i<6;i++) {
      if (std::abs(cell_[i]) < tol) return true;
    }
    return false;
  }
  //--------------------------------------------------------------
  /*! If the "difference" between cells is greater than half the resolution,
    then a reflection may be mapped on to the next reflection in the other cell
   */
  double Scell::Difference(const Scell& other) const
  // "difference" between two cells in A, maximum allowed distance
  // algorithm from clipper::Cell::equals
  // Even Kevin doesn't understand this, though presumably V^1/3 is an
  //  average cell length
  {
    double s = 0.0;
    for (int j=0;j<3;++j) {
      for (int i=0;i<3;++i) {
	s += (Bmat_(i,j) - other.Bmat_(i,j))*(Bmat_(i,j) - other.Bmat_(i,j));
      }}  // sum of squares of orthogonalisation matrix differences
    // Average cell volume
    double vav = (Volume() + other.Volume())*0.5;
    double v43 = pow(vav, 1.333333); // V^4/3
    //    std::cout << "\nScell::Difference " << s
    //	      << " " << Volume()<< " " << other.Volume() << "\n"; //^
    return sqrt(v43 * s);
  }
  //--------------------------------------------------------------
  //! Construct from list of cells
  UnitCellSet::UnitCellSet(const std::vector<Scell>& Cells) :cells(Cells)
  {
    averagecell = Average(); // average all cells
  }
  //--------------------------------------------------------------
  //! Initialise from list of cells
  void UnitCellSet::init(const std::vector<Scell>& Cells)
  {
    cells = Cells;
    averagecell = Average(); // average all cells
  }
  //--------------------------------------------------------------
  void UnitCellSet::AddCell(const Scell& Cell) 
  {
    cells.push_back(Cell); 
    averagecell = Average(); // average all cells
  }
  //--------------------------------------------------------------
  Scell UnitCellSet::Average(const int& idxexclude) const
  // Average list of cells
  // if idxexclude >= 0, exclude entry with this index, < 0 include all
  {
    int nc = cells.size();
    ASSERT (idxexclude < nc);
    if (idxexclude >= 0) nc--;
    if (nc <= 0) {return Scell();}
    else if (nc == 1) {
      int j = (idxexclude+1)%2; // 1 or 0
      return cells[j];
    }
    // We have 2 or more,average
    std::vector<Dtype> sumcell(6,0.0);
    for (size_t k=0; k<cells.size(); k++) {
      if (idxexclude < 0 || int(k) != idxexclude) {
	for (int i=0; i<6; i++) {sumcell[i] += cells[k].UnitCell()[i];}
      }
    }
    for (int i=0; i<6; i++) {sumcell[i] /= nc;}
    return Scell(sumcell);
  }
  //--------------------------------------------------------------
  //! return list of deviations (A) from average of other cells
  std::vector<double> UnitCellSet::Deviations() const
  {
    std::vector<double> dv(Number(), 0.0);
    if (Number() <= 1) return dv;
    Scell cell;
    for (int i=0;i<Number();++i) {
      cell = Average(i); // average cell excluding i'th
      dv[i] = cell.Difference(cells[i]); // difference between cells
    }
    return dv;
  }
  //--------------------------------------------------------------
  //! return worst deviation (A)
  double UnitCellSet::WorstDeviation() const
  {
    double d = -10000.;
    std::vector<double> dv = Deviations();
    for (size_t i=0;i<dv.size();++i) {
      d = Max(d, dv[i]);
    }
    return d;
  }
  //--------------------------------------------------------------
  //! return list of rms deviations for each cell parameter
  std::vector<double> UnitCellSet::RmsD() const
  {
    std::vector<double> dv(6, 0.0);
    if (Number() <= 1) return dv;
    std::vector<MeanSD> rms(6);
    for (int i=0;i<Number();++i) {
      std::vector<Dtype> cv =  cells[i].UnitCell();
      for (int j=0;j<6;++j) {
	rms[j].Add(cv[j]);
      }
    }
    for (int j=0;j<6;++j) {
      dv[j] = rms[j].SD();
    }
    return dv;
  }
  //--------------------------------------------------------------
  // Change basis: new h' = h * reindex_op 
  // return false if indices are non-integral
  // reindex_op may include translations
  bool Hkl::change_basis(Hkl& Newhkl, const ReindexOp& reindex_op) const
  {
    Vec3<double> vh(*this);
    Vec3<double> v = vh * reindex_op;
    Newhkl = Hkl(Nint(v[0]), Nint(v[1]), Nint(v[2]));
    bool OK = true;
    for (int i=0;i<3;i++) {
      if (std::abs(v[i]-double(Newhkl[i])) > 0.05) OK = false;
    }
    return OK; 
  }
  //--------------------------------------------------------------
  Hkl Hkl::change_basis(const ReindexOp& reindex_op) const
  // reindex_op may include translations
  {
    Vec3<double> vh(*this);
    Vec3<double> v = vh * reindex_op;
    return Hkl(Nint(v[0]), Nint(v[1]), Nint(v[2]));
  }

  //--------------------------------------------------------------
  std::string Hkl::format() const
  { return "HKL = ("+String(h())+","+String(k())+","+String(l())+")"; }

  //--------------------------------------------------------------
  int Hkl::code() const
  // encode hkl triple as integer (not guaranteed unique)
  {
    return (((*this)[0] & 0x3FF) << 20) |
      (((*this)[1] & 0x3FF) << 10) |
      (((*this)[2] & 0x3FF));
  }
  //--------------------------------------------------------------
  Hkl HklDecode(const int& kode)
  // decode integer as hkl triple
  {
    int h = (kode << 2) >> 22;
    int k = (kode << 12) >> 22;
    int l = (kode << 22) >> 22;
    return Hkl(h,k,l);
  }
  //--------------------------------------------------------------
  bool Hkl::IsGeneral() const
  {
    if (h()==0 || k()==0 || l()==0) return false;
    if ((h()==k()) || (k()==l()) || (l()==h())) return false;
    return true;
  }
  //--------------------------------------------------------------
  void PxdName::update(const PxdName& NewPxd)
  // Replace elements by any non-blank elements in NewPxd
  {
    if (NewPxd.pname() != "") {
      pname_ = NewPxd.pname();
    }
    if (NewPxd.xname() != "") {
      xname_ = NewPxd.xname();
    }
    if (NewPxd.dname() != "") {
      dname_ = NewPxd.dname();
    }
  }
  //--------------------------------------------------------------
  std::string PxdName::formatPrint() const
  {
    int lm = Max(pname_.size(), xname_.size());
    lm = Max(lm, int(dname_.size()));
    lm++;
    int lpn = Max(lm, int(pname_.size()+1));
    int lxn = Max(lm, int(xname_.size()+1));
    int ldn = Max(lm, int(dname_.size()+1));
    return FormatOutput::logTab(2,
				StringUtil::PadString("Project: "+pname_, lpn)+
				StringUtil::PadString(" Crystal: "+xname_, lxn)+
				StringUtil::PadString(" Dataset: "+dname_, ldn));
  }
  //--------------------------------------------------------------
  std::string PxdName::format() const
  {
    return  pname_ + "/" + xname_ + "/" + dname_;
  }
  //--------------------------------------------------------------
  Xdataset::Xdataset(const PxdName& pxdname, const Scell& cell,
		     const float& wavel, const int& setid)
    : pxdname_(pxdname), setid_(setid), cell_(cell), wavel_(wavel)
  {
    allcells_.AddCell(cell);
    allwavel_.push_back(wavel);
  }
  //--------------------------------------------------------------
  void Xdataset::add_batch(const int& batch_num)
  {
    batches.push_back(batch_num);
  }
  //--------------------------------------------------------------
  void Xdataset::AddRunIndex(const int& RunIndex)
  {
    run_index_list.push_back(RunIndex);
  }
  //--------------------------------------------------------------
  void Xdataset::AddCellWavelength(const Scell& newcell, const float& wavel)
  //! add in another cell and wavelength
  {
    allcells_.AddCell(newcell);
    allwavel_.push_back(wavel);
    AverageCellWavelength();
  }
  //--------------------------------------------------------------
  std::string Xdataset::formatPrint() const
  {
    std::string s = FormatOutput::logTab(1,"\n * Dataset information *\n");
    s += pxdname_.formatPrint();
    s += FormatOutput::logTabPrintf(1,"Unit cell:  ");
    s += cell_.formatPrint();
    s += FormatOutput::logTabPrintf(1,
				    "Wavelength: %8.5f\n", wavel_);
    s += FormatOutput::logTabPrintf(1,"Runs: ");
    for (size_t i=0;i<run_index_list.size();i++)
      {s += FormatOutput::logTabPrintf(1," %3d", run_index_list[i]+1);}
    return s+"\n";
  }
  //--------------------------------------------------------------
  //! average multiple cells and wavelengths, return false if they differ by more than tolerance
  void Xdataset::AverageCellWavelength()
  // Tolerance in A
  {
    ASSERT (unsigned(allcells_.Number()) == allwavel_.size());
    if (allcells_.Number() == 0) return;
    if (allcells_.Number() == 1) { // just one
      cell_ = allcells_.Cells()[0];
      wavel_ = allwavel_[0];
      return;
    }
    // allcells_ is an UnitCellSet object
    cell_ = allcells_.Average();
    wavel_ = MeanSD(allwavel_).Mean();
  }
  //--------------------------------------------------------------
  double Xdataset::WorstDeviation() const
  {
    return allcells_.WorstDeviation();
  }
  //--------------------------------------------------------------
  std::string Xdataset::formatAllCells() const
  //!< format cell & wavelength list if more than one
  {
    std::string s;
    if (allcells_.Number() <= 1) return s;
    std::string blank(22,' ');
    double devmax = WorstDeviation();
    const double TOL = 0.001;
    if (devmax < TOL) {
      s += "         Files for dataset contain "+clipper::String(allcells_.Number(), 3)+
	" near identical cells, "+
	" maximum deviation "+StringUtil::ftos(devmax, 6,4)+"\n";
      return s;
    }
    s += "\n"+std::string(10,' ')+
      "'deviation' is estimate of the difference of the cell (in A)\n";
    s += std::string(12,' ')+
      "from the mean of the others\n";
    s += std::string(10,' ')+
      "At worst this should be less the half the maximum resolution of the data\n\n";
    s += blank +
      "    a       b       c     alpha    beta   gamma  deviation lambda\n";
    if (allcells_.Number() <= 1) return s;
    std::vector<Scell> cells = allcells_.Cells();  // cells
    // differences from mean of others
    std::vector<double> dv = allcells_.Deviations();

    for (int i=0;i<allcells_.Number();++i) {
      s += blank;
      s += cells[i].formatPrint(false)+
	StringUtil::ftos(dv[i], 7,2)+
	StringUtil::ftos(allwavel_[i],9,4)+
	"\n";
    }
    s += std::string(7,' ')+"RMS deviation: ";
    dv = allcells_.RmsD();
    for (int i=0;i<6;++i) {
      s += StringUtil::ftos(dv[i], 7,2)+" ";
    }
    return s;
  }
  //--------------------------------------------------------------
  void Xdataset::change_basis(const ReindexOp& reindex_op)
  {
    cell_ = cell_.change_basis(reindex_op);
  }
  //--------------------------------------------------------------
  bool operator == (const Xdataset& a,const Xdataset& b)
  // Equality just tests pxdname 
  {
    return (a.pxdname_ == b.pxdname_);
  }
  //--------------------------------------------------------------
  DMat33 MakePermutationMatrix(const DVect3& s0, const DVect3& e0)
  // Permutation matrix to convert header frame to Cambridge frame
  //  Cambridge frame has
  //    source vector s0 (anti-parallel to beam) approximately along -x
  //      ie x is along beam (if perpendicular to rotation axis)
  //    principal rotation axis e0 exactly along z
  //  x(Cambridge) = [Q] x(header)
  {
    // xq = first column of [Q] = -s0 (for now)
    DVect3 xq = (-1.0 * s0).unit();
    // zq = e0
    DVect3 zq = e0.unit();
    // yq = zq x xq
    DVect3 yq = DVect3::cross(zq, xq).unit();
    // xq = yq x zq
    xq = DVect3::cross(yq, zq).unit();

    return DMat33(xq[0], yq[0], zq[0],
		  xq[1], yq[1], zq[1],
		  xq[2], yq[2], zq[2]);
  }
  //--------------------------------------------------------------
  Batch::Batch()
  // A dummy batch
  {
    batchinfo.num = 1;
    Xdataset_index = 0;
    run_index = -1;
    accepted = true;
    valid_cell = false;
    valid_Umat = false;  
    valid_time = 0;
    valid_phi = false;
    phioffset = 0;
    offset = 0;
    file_num = 1;
    phiscan = false;
    pole = 0;
    initBatchInfo();
  } // Batch constructor
  //--------------------------------------------------------------
  void Batch::initBatchInfo()
  // initialise batchinfo to zeroes
  {
    int intbuf[NBATCHINTEGERS];
    float fltbuf[NBATCHREALS];
    for (int i=0;i<NBATCHINTEGERS;++i) {intbuf[i] = 0;}
    for (int i=0;i<NBATCHREALS;++i) {fltbuf[i]=0.0;}
    int status = CMtz::MtzArrayToBatch(intbuf, fltbuf, &batchinfo);
    if (!status) {
      Message::message(Message_fatal
		       ("Batch::initBatchInfo fail"));
    }
    strcpy(batchinfo.title, "");                /**< batch title */	      
    strcpy(batchinfo.gonlab[0], "        ");    /**< names of the three axes */
    strcpy(batchinfo.gonlab[1], "        ");    // 8-characters only! 
    strcpy(batchinfo.gonlab[2], "        ");
  }
  //--------------------------------------------------------------
  Batch::Batch(const CMtz::MTZBAT& batch,
	       const bool& accept, const int& idataset)
    :   batchinfo(batch), Xdataset_index(idataset), accepted(accept)
  {
    offset = 0;
    file_num = 1;
    run_index = -1;
    phioffset = 0.0;
    phiscan = false;
    pole = 0;
    PDUinv = DMat33::identity();
    E1E2 = DMat33::identity();

    const float tolerance = 0.001;
    valid_cell = true;
    // Check validity of cell & orientation matrix
    if (batchinfo.cell[0]*batchinfo.cell[1]*batchinfo.cell[2] == 0.0)
      {valid_cell = false;}
    if (valid_cell) {bcell = Scell(batchinfo.cell);}

    U = MVutil::SetCMat33(batchinfo.umat); // [U]
    valid_Umat = true;
    if (U.det() <= 0.9)   valid_Umat = false;  
    // If Umat == Identity, it is probably invalid
    if (MVutil::is_mat33_ident(U)) {
      valid_Umat = false;
    }
    if (!valid_cell) valid_Umat = false;

    if (valid_Umat)  {
      // Generate derived matrices etc in "Cambridge" frame
      //  Cambridge frame has
      //    source vector s0 (anti-parallel to beam) along -x
      //      ie x is along beam
      //    principal rotation axis e0 along z
      // Rotation axes
      DVect3 e1(batchinfo.e1[0],batchinfo.e1[1],batchinfo.e1[2]);
      DVect3 e2(batchinfo.e2[0],batchinfo.e2[1],batchinfo.e2[2]);
      DVect3 e3(batchinfo.e3[0],batchinfo.e3[1],batchinfo.e3[2]);
      // Principal rotation axis e0
      DVect3 e0 = e1;
      // Source vector from batch header, anti-parallel to beam
      DVect3 ss0(batchinfo.source[0],batchinfo.source[1],batchinfo.source[2]);
      // Permutation matrix to convert header frame to Cambridge frame
      Q = MakePermutationMatrix(ss0, e0);
      // Transform relevent vectors
      s0 = Q * ss0;  // source vector, stored
      // Datum matrix
      e1 = Q * e1;
      e2 = Q * e2;
      e3 = Q * e3;
      DVect3 dtm(batchinfo.datum[0],batchinfo.datum[1],batchinfo.datum[2]);
      // angles = 0 for non-existent axes
      if (batchinfo.ngonax == 1) {dtm[1] = 0.0;}
      if (batchinfo.ngonax < 3)  {dtm[2] = 0.0;}
      // Datum matrix rotation
      if (batchinfo.ngonax == 2) {
	// Two axes, construct e3 orthogonal to e1 & e2
	e3 = clipper::Vec3<>::cross(e1, e2).unit();
      }
      spindle = e1;  // phi/omega scan
      scala::Rotation D;
      if (batchinfo.ngonax == 1) {
	// Single axis, e1
	D = scala::Rotation(e1, dtm[0]);
      } else {
	if (batchinfo.jsaxs == 1) {
	  // Omega scan
	  D = scala::Rotation(scala::Euler_explicit
				(e1, clipper::Util::d2rad(dtm[0]),
				 e2, clipper::Util::d2rad(dtm[1]),
				 e3, clipper::Util::d2rad(dtm[2])));
	  phiscan = false;
	} else if (batchinfo.jsaxs == 3) {
	  // Phi scan
	  // datum [D] = [Phi0]
	  D = scala::Rotation(scala::Euler_explicit
				(e1, 0.0,
				 e2, 0.0,
				 e3, clipper::Util::d2rad(dtm[2])));
	  phiscan = true;
	  spindle = e3;  // phi/omega scan
	  // E1E2 = [Omega][Chi/Kappa]
	  E1E2 = scala::Rotation(scala::Euler_explicit
				   (e1, clipper::Util::d2rad(dtm[0]),
				    e2, clipper::Util::d2rad(dtm[1]),
				    e3, 0.0)).matrix();
	}
      }
      // "Missetting" angles [A] = [PhiXYZ][U]
      DMat33 A = U;
      if (batchinfo.misflg > 0)	{
	DVect3 x(1.0,0.0,0.0);
	DVect3 y(0.0,1.0,0.0);
	DVect3 z(0.0,0.0,1.0);
	scala::Rotation PhiXYZ(scala::Euler_explicit
				 (x, clipper::Util::d2rad(batchinfo.phixyz[0][0]),
				  y, clipper::Util::d2rad(batchinfo.phixyz[0][1]),
				  z, clipper::Util::d2rad(batchinfo.phixyz[0][2])));
	A = PhiXYZ.matrix() * U;
	U = A;
      }
      A = Q * A;  // transform to frame
      // Store [D][PhiXYZ][U]
      DU = D.matrix() * A;
      // [D][U][B] ([B] in is A^-1)
      DUB = DU * bcell.Bmat();
    } // end if valid_Umat
    
    float phirange = batchinfo.phiend  - batchinfo.phistt;
    valid_phi = true;
    if (std::abs(phirange) < tolerance) valid_phi = false;

    float timerange = batchinfo.time2  - batchinfo.time1;
    valid_time = +1;
    if (std::abs(timerange) < tolerance) {
      valid_time = 0;
      // Substitute phi
      if (valid_phi) {
	timerange = phirange;
	batchinfo.time1 =  batchinfo.phistt;
	batchinfo.time2 =  batchinfo.phiend;
	valid_time = -1;
      }
    }

    if (batchinfo.phirange > tolerance) {
      // phirange is explicitly in file, check phi2 modulo 360
      if (std::abs(phirange - batchinfo.phirange) > tolerance) {
	double dif = batchinfo.phistt+batchinfo.phirange-batchinfo.phiend;
	if (std::abs(fmod(dif,360.0)) > tolerance) {
	  Message::message(Message_warn(
		"Inconsistent rotation information in header for batch "+
					clipper::String(batchinfo.num)+
		"\n  Phi1: "+clipper::String(batchinfo.phistt)+
		" Phi2: "+clipper::String(batchinfo.phiend)+
		" DelPhi: "+clipper::String(batchinfo.phirange)));
	} else
	  // Fix up phi end if off by multiple of 360
	  {batchinfo.phiend = batchinfo.phistt + batchinfo.phirange;}
      }
    }

  } // Batch constructor
  //--------------------------------------------------------------
  void Batch::change_basis(const ReindexOp& reindex_op)
  // change cell and orientation matrix for change of basis
  {
    // Check for valid cell
    if (valid_cell) {
      Scell bcell(batchinfo.cell);
      // New cell
      Scell newcell = bcell.change_basis(reindex_op);
      for (int i=0;i<6;i++) batchinfo.cell[i] = newcell[i];
      // Check for valid orientation
      if (valid_Umat) {
	// Rotate crystal orientation
	DMat33 U = MVutil::SetCMat33(batchinfo.umat);
	// U' = U * B * H(T)^-1 * B'^-1
	DMat33 RU = U *
	  bcell.Bmat() * reindex_op.transpose().inverse() * newcell.Bmat().inverse();
	for (int i=0;i<3;i++) {
	  for(int j=0;j<3;j++) {
	    batchinfo.umat[j*3+i] = RU(i,j);}}
	// Update DU' = DU * U^-1 * U'
	DU = DU * U.inverse() * RU;
	// [D][U][B]
	DUB = DU * bcell.Bmat();
      }
    }
  }
  //--------------------------------------------------------------
  void Batch::print() const
  // Can't do this to "Output" object
  {
    CMtz::MtzPrintBatchHeader(&batchinfo);
  }
  //--------------------------------------------------------------
  void Batch::SetCellConstraint(const std::vector<int>& lbcell)
  {
    for (int i=0;i<6;i++) batchinfo.lbcell[i] = lbcell[i];
  }
  //--------------------------------------------------------------
  void Batch::SetCell(const Scell& newcell)
  {
    for (int i=0;i<6;i++) batchinfo.cell[i] = newcell[i];
    valid_cell = true;
    bcell = newcell;
  }
  //--------------------------------------------------------------
  void Batch::OffsetPhi(const float& Phioffset)
  // offset stored phi1, phi2 to allow to 360deg wrap-around
  {
    phioffset = Phioffset;
    batchinfo.phistt += phioffset;
    batchinfo.phiend += phioffset;
  }
  //--------------------------------------------------------------
  void Batch::OffsetTime(const float& Timeoffset)
  // offset stored time1, time2 to allow to 360deg wrap-around
  // if copied from Phi
  {
    timeoffset = Timeoffset;
    batchinfo.time1 += timeoffset;
    batchinfo.time2 += timeoffset;
  }
  //--------------------------------------------------------------
  void Batch::SetRunIndex(const int& RunIndex)
  {
    run_index= RunIndex;
  }
  //--------------------------------------------------------------
  // offset batch number
  void Batch::OffsetNum(const int& Offset) {
    batchinfo.num += Offset;
    offset = Offset;  // and store it
  }
  //--------------------------------------------------------------
  // Orientation information:
  //   x = [R] [D] [U] [B] h
  //  where h   = (h k l)
  //       [B]  orthogonalisation matrix from unit cell
  //       [U]  orientation matrix
  // Two cases for rotation:
  //   1) rotation is about the outer goniostat axis e1
  //      (a) single axis Phi around e1
  //          [R] = [Phi]      [D] = [Phi0], often = [I]
  //      (b) 3-axis goniostat, spindle rotation by omega around e1
  //          [R] = [Omega]    [D] = [Chi/Kappa][Phi]
  //   2) 3-axis goniostat, spindle rotation by phi around e3
  //          [R] = [Omega][Chi/Kappa][Phi]     [D] = [Phi0]
  //              = [E1E2] [Phi]  (phi varying)
  //              
  // Transformation operations:
  //   four useful reference frames
  //    (a) reciprocal lattice indices h
  //    (b) camera frame (Ewald sphere frame)
  //        s(r)  = x = [R][D][U][B]h
  //    (c) zero rotation angle frame (phi/omega = 0)
  //        s(r0) = [D][U][B]h = [R]^-1 s(r)
  //    (d) polar orthogonalised crystal frame to put a defined
  //      reciprocal lattice axis ~ along z
  //        s(pole) = [P][B]h
  //                = [P][DU]^-1 s(r0)
  //       [P] is a cyclic permutation matrix (so [P]^-1 = [P]transpose)
  //
  // NB in the following FtoF rountines, all vectors (apart from hkl) are
  // in dimensionless reciprocal lattice units (ir s(r), s(r0), s(pole)
  //--------------------------------------------------------------
  DVect3 Batch::HtoSr(const Hkl& hkl, const float& phi) const
  // reciprocal lattice index h -> s(r) = [R][D][U][B]h camera frame
  {
    DVect3 hl;  // lambda * hkl
    for (int i=0;i<3;++i) {
      hl[i] = batchinfo.alambd * hkl[i];
    }
    // spindle rotation [Phi]
    //   for single axis or omega scan [R] = [Phi]
    DMat33 R = scala::Rotation
      (spindle, clipper::Util::d2rad(phi)).matrix();
    if (phiscan) {
      // 3-axis Phi scan, [R] = [E1E2][Phi]
      R = E1E2 * R;
    }
    return R * DUB * hl;
  }
  //--------------------------------------------------------------
  DVect3 Batch::HtoSr0(const Hkl& hkl) const
  // reciprocal lattice index h -> s(r0) = [D][U][B]h zero rotation angle frame
  {
    DVect3 hl;  // lambda * hkl
    for (int i=0;i<3;++i) {
      hl[i] = batchinfo.alambd * hkl[i];
    }
    return DUB * hl;
  }
  //--------------------------------------------------------------
  DVect3 Batch::SrtoSr0(const DVect3& sr, const float& phi) const
  // camera frame s(r) -> s(r0) = [R]^-1 s(r) zero rotation angle frame
  {
    // spindle rotation [Phi]
    //   for single axis or omega scan [R] = [Phi]
    DMat33 R = scala::Rotation
      (spindle, clipper::Util::d2rad(phi)).matrix();
    if (phiscan) {
      // 3-axis Phi scan, [R] = [E1E2][Phi]
      R = E1E2 * R;
    }
    return R.transpose() * sr;
  }
  //--------------------------------------------------------------
  DVect3 Batch::Sr0toP(const DVect3& sr0) const
  // zero rotation angle frame s(r0) -> polar orthogonalised crystal frame
  //    s(pole) = [P] [DU]^-1 s(r0)
  {
    return PDUinv * sr0;
  }
  //--------------------------------------------------------------
  void Batch::SetPole(const int& Pole)
  {
    pole = Pole;
    if (pole < 0) {
      //  Unspecified ABSORPTION pole, choose closest reciprocal axis
      //    numbered 1,2,3
      pole =
	SmallestComponent(SpindleToPrincipleAxes(DUB)) + 1;
    }
    // Generate permutation matrix P for this batch to put pole along z
    // polar orthogonal crystal frame  x(polar) = [P] [B] h
    if (pole != 0) {
      DMat33 P(0,0,0,0,0,0,0,0,0);
      int l = pole-1;  // pole = 1,2,3
      P(2, l) = 1.0;    // z axis
      P(0,++l%3) = 1.0; // x axis
      P(1,++l%3) = 1.0; // y axis
      PDUinv = (DU * P).inverse();
    }
  }
  //--------------------------------------------------------------
  //! Return range of detector pixel coordinates
  std::vector<std::vector<float> > Batch::DetectorCoordinateRange() const
  {
    std::vector<float> lims(2);
    std::vector<std::vector<float> > xylims(2);
    lims[0] = batchinfo.detlm[0][0][0]; // XDET
    lims[1] = batchinfo.detlm[0][0][1];
    xylims[0] = lims;
    lims[0] = batchinfo.detlm[0][1][0]; // YDET
    lims[1] = batchinfo.detlm[0][1][1];
    xylims[1] = lims;
    return xylims;
  }
  //--------------------------------------------------------------
  std::string Batch::SpindleToPrincipleAxis() const
  {
    return FormatSpindleToPrincipleAxes(DUB);
  }
  //--------------------------------------------------------------
  std::string Batch::format() const
  // formatted version
  // (most of this copied from cmtzlib.c CMtz::MtzPrintBatchHeader)
  {
    std::string s = "Orientation data for batch "+clipper::String(num())+"   ";
    if (batchinfo.ldtype == 1) {
      s += "oscillation (2D) data\n\n";
    } else if (batchinfo.ldtype == 2) {
      s += "area detector (3D) data\n\n";
    } else if (batchinfo.ldtype == 3) {
      s += "Laue data\n\n";
    } else {
      s += "unknown data type\n\n";
    }
    s += "   Crystal number ..................."+clipper::String(batchinfo.ncryst)+
       "\n   Associated dataset ID ............"+clipper::String(batchinfo.nbsetid)+
       "\n   Cell dimensions .................."+FormatCell(cell())+"\n";
    s += "   Cell fix flags ...................";
    for (int i=0;i<6;++i) {
      s += clipper::String(batchinfo.lbcell[i],7);
    }
    std::string s2;
    if (!batchinfo.misflg) {
      s += "\n   Orientation matrix U .............";
      s2 = "\n       (including setting angles)    ";
    } else {
      s += "\n   Standard orientation matrix [U] ..";
      s2 = "\n                                     ";
    }
    for (int i=0;i<3;++i) {for (int j=0;j<3;++j) {
	s += FormatOutput::logTabPrintf(0,"%10.4f",U(i,j));
      }
      if (i<2) {s += s2;}
      s2 = "\n                                     ";
    }
    s += "\n";
    if (batchinfo.misflg == 1) {
      s += FormatOutput::logTabPrintf(0,"   %s %6.2f %6.2f %6.2f\n",
	 "Missetting angles PhiX PhiY PhiZ..",
         batchinfo.phixyz[0][0],batchinfo.phixyz[0][1],batchinfo.phixyz[0][2]);
    } else if (batchinfo.misflg > 1) {
      s += FormatOutput::logTabPrintf(0,"   %s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
         "Missetting angles PhiX PhiY PhiZ..",
         batchinfo.phixyz[0][0],batchinfo.phixyz[0][1],batchinfo.phixyz[0][2],
         batchinfo.phixyz[1][0],batchinfo.phixyz[1][1],batchinfo.phixyz[1][2]);
    }
    std::string axis("none");
    if (batchinfo.jumpax == 1) axis = "a*";
    if (batchinfo.jumpax == 2) axis = "b*";
    if (batchinfo.jumpax == 3) axis = "c*";
    s += FormatOutput::logTabPrintf(0,"   %s%s%s   %s\n",
	    "Reciprocal axis nearest ",batchinfo.gonlab[0],"..",axis.c_str());
    if (!batchinfo.lcrflg) {
    s += FormatOutput::logTabPrintf(0,"   %s %6.3f \n",
	   "Mosaicity ........................",batchinfo.crydat[0]);
    } else {
      s += FormatOutput::logTabPrintf(0,"   %s %6.3f %6.3f \n",
	   "Mosaicity (horizontal, vertical)..",batchinfo.crydat[0],batchinfo.crydat[1]);
    }
    s += FormatOutput::logTabPrintf(0,"   Datum goniostat angles (degrees)..");
    for (int i = 0; i < batchinfo.ngonax; ++i) 
      s += FormatOutput::logTabPrintf(0," %8.3f",batchinfo.datum[i]);
    s += FormatOutput::logTabPrintf(0,"\n");

    if (batchinfo.jsaxs > 0 && batchinfo.jsaxs <= batchinfo.ngonax) 
      s += FormatOutput::logTabPrintf(0,"   %s  %s \n",
      "Scan axis ........................",batchinfo.gonlab[batchinfo.jsaxs-1]);
    s += FormatOutput::logTabPrintf(0,"   %s %8.3f %8.3f \n   %s %8.3f \n   %s %8.2f %8.2f \n",
	 "Start & stop Phi angles (degrees).",
	   batchinfo.phistt,batchinfo.phiend,
	 "Range of Phi angles (degrees).....",batchinfo.phirange,
         "Start & stop time (minutes).......",batchinfo.time1,batchinfo.time2);

    if (batchinfo.nbscal == 4) {
      s += FormatOutput::logTabPrintf(0,"   %s %9.4f %9.4f \n   %s %9.4f %9.4f \n",
      "   Batch scale & SD .................",batchinfo.bscale,batchinfo.sdbscale,
      "   Batch B-factor & SD ..............",batchinfo.bbfac,batchinfo.sdbfac);
    }

    s += FormatOutput::logTabPrintf(0,"   %s  \n   %s %7d \n   %s %s %s %9.4f %9.4f %9.4f \n   %s %s %s %9.4f %9.4f %9.4f \n   %s %s %s %9.4f %9.4f %9.4f \n",
         " Crystal goniostat information :-",
	 "   Number of goniostat axes..........",batchinfo.ngonax,
	 "   Goniostat vectors.....",
	    batchinfo.gonlab[0],"....",batchinfo.e1[0],batchinfo.e1[1],batchinfo.e1[2],
	 "                    .....",batchinfo.gonlab[1],"....",batchinfo.e2[0],batchinfo.e2[1],batchinfo.e2[2],
	 "                    .....",batchinfo.gonlab[2],"....",batchinfo.e3[0],batchinfo.e3[1],batchinfo.e3[2]);

    s += FormatOutput::logTabPrintf(0,
      "   %s \n   %s%9.4f %9.4f %9.4f \n   %s%9.4f %9.4f %9.4f \n",
	 " Beam information :-",
	 "   Idealized X-ray beam vector.......",
	 batchinfo.source[0],batchinfo.source[1],batchinfo.source[2],
	 "   X-ray beam vector with tilts......",
	 batchinfo.so[0],batchinfo.so[1],batchinfo.so[2]);

    if (batchinfo.lbmflg == 0) {
      s += FormatOutput::logTabPrintf(0,"   %s %9.5f %9.5f \n",
           "   Wavelength and dispersion ........",batchinfo.alambd,batchinfo.delamb);
    } else if (batchinfo.lbmflg == 1) {
      s += FormatOutput::logTabPrintf(0,"   %s %9.5f %9.5f %9.5f \n   %s %7.3f %7.3f \n",
      "   Wavelength and dispersion ........",batchinfo.alambd,batchinfo.delamb,batchinfo.delcor,
      "   Divergence .......................",batchinfo.divhd,batchinfo.divvd);
    }
    
    s += FormatOutput::logTabPrintf(0," Detector information :-\n   Number of detectors...............%7d \n",batchinfo.ndet);
    s += FormatOutput::logTabPrintf(0,"   %s%9.3f\n%s%9.3f\n%s%7.1f%7.1f%7.1f%7.1f\n",
    "   Crystal to Detector distance (mm).",batchinfo.dx[0],
    "   Detector swing angle..............",batchinfo.theta[0],
    "   Pixel limits on detector..........",
    batchinfo.detlm[0][0][0],batchinfo.detlm[0][0][1],batchinfo.detlm[0][1][0],batchinfo.detlm[0][1][1]);
    return s;
  }
  //--------------------------------------------------------------
    // For sorting on batch number
  bool operator < (const Batch& a,const Batch& b)
  {
    return (a.batchinfo.num < b.batchinfo.num);
  }
  //--------------------------------------------------------------
  bool in_datasets(const int& setid,
		   const std::vector<Xdataset>& datasets,
		   int& idataset)
  // Return true if dataset setid is in datasets list
  //  & return dataset index idataset (-1 if not)
  // If setid == 0, assign to first dataset 
  {
    if (setid <= 0) {
	idataset = 0;
	return true;
    }
    for (size_t k = 0; k < datasets.size(); ++k) {
      if (setid == datasets[k].setid()) {
	idataset = k;
	return true;
      }
    }
    idataset = -1;
    return false;
  }
  //--------------------------------------------------------------
  BatchSelection::BatchSelection()
  {clear();}
  //--------------------------------------------------------------
  void BatchSelection::clear()
  {
    batchlist.clear();
    fileseries_list.clear();
    batchranges.clear();
    fileseries_range.clear();
    flaglist.clear();
  }
  //--------------------------------------------------------------
  void BatchSelection::AddBatch(const int& batch, const int& fileSeriesList)
  // Add in one batch to selection
  {
    batchlist.push_back(batch);
    fileseries_list.push_back(fileSeriesList);
  }
  //--------------------------------------------------------------
  void BatchSelection::AddRange(const int& batch1, const int& batch2,
				const int& fileSeriesList,
				const int& flag)
  // Add in batch range to selection
  {
    batchranges.push_back(IntRange(batch1,batch2));
    fileseries_range.push_back(fileSeriesList);
    flaglist.push_back(flag);
  }
  //--------------------------------------------------------------
  bool BatchSelection::InSelection(const int& batch,
				   const int& fileSeriesTest) const
  {
    return (FindInSelection(batch, fileSeriesTest) >= 0);
  }
  //--------------------------------------------------------------
  int BatchSelection::FindInSelection(const int& batch,
				      const int& fileSeriesTest) const
  // Return index in list if in selection, for given file series, else -1
  // Two possibilities for selection:
  //  1) Specified selection on final numbering, ie from sole file or
  //  after renumbering of 2nd or subsequent file, fileSeriesList == 0,
  //  only test if fileSeriesTest == 0
  //
  //  2) Specified selection on original file numbering from 2nd or
  //  subsequent file, fileSeriesList > 0, only test if
  //  fileSeriesTest > 0
  {
    for (size_t i=0;i<batchlist.size();i++) {
      if (fileseries_list[i] == 0 && fileSeriesTest == 0) {
	// final numbering
	if (batch == batchlist[i]) return int(i);
      } else if (fileseries_list[i] > 0 && fileseries_list[i] == fileSeriesTest) {
	// original numbering
	if (batch == batchlist[i]) return int(i);
      }
    }
    for (size_t i=0;i<batchranges.size();i++) {
      if (fileseries_range[i] == 0 && fileSeriesTest == 0) {
	// final numbering
	if (batchranges[i].InRange(batch)) return int(i);
      } else if (fileseries_range[i] > 0 && fileseries_range[i] == fileSeriesTest) {
	// original numbering
	if (batchranges[i].InRange(batch)) return int(i);
      }
    }
    // not in range or not tested
    return -1;
  }
  //--------------------------------------------------------------
    //! return  numerical flag for batch, = -1 if not in list
  int BatchSelection::FlagNumber(const int& batch)
  {
    int i;
    if (batchlist.size() > 0) {
      Message::message(Message_fatal
       ("BatchSelection::FlagNumber not valid for batch list, only for batch range"));
    }
    if ((i = FindInSelection(batch, 0)) >= 0) {
      return flaglist[i];
    } else {
      return -1;
    }
  }
  //--------------------------------------------------------------
  //! return list of unique numerical flags specified
  std::vector<int> BatchSelection::UniqueFlags() const
  {
    std::vector<int> uniqueflags;
    for (size_t i=0;i<flaglist.size();++i) {
      if (uniqueflags.size() > 0) {
	if (std::find(uniqueflags.begin(), uniqueflags.end(), flaglist[i])
	    != uniqueflags.end()) {continue;} // skip if already in unique list
      }
      uniqueflags.push_back(flaglist[i]);
    }
    return uniqueflags;
  }
  //--------------------------------------------------------------
  //! return list of batch ranges for specified flag, all ranges if flag < 0
  std::vector<IntRange> BatchSelection::BatchRanges(const int& flag)
  {
    ASSERT (batchranges.size() == flaglist.size());
    if (flag < 0) return batchranges;
    std::vector<IntRange> bsel;
    for (size_t i=0;i<batchranges.size();++i) {
      if (flaglist[i] == flag) {
	bsel.push_back(batchranges[i]);
      }
    }
    return bsel;
  }
  //--------------------------------------------------------------
  void BatchSelection::CheckSeries(const int& NumFileSeries) const
  // Check that fileSeries specified on selection commands match
  // specified files. Fail here if not
  {
    for (size_t i=0;i<batchlist.size();i++) {
      if (fileseries_list[i] > NumFileSeries) {
	Message::message(Message_fatal
			 ("Batch selection specifies file[series] "+
			  clipper::String(fileseries_list[i],3)+
			  " which does not exist"));
      }
    }
    for (size_t i=0;i<batchranges.size();i++) {
      if (fileseries_range[i] > NumFileSeries) {
	Message::message(Message_fatal
			 ("Batch selection specifies file[series] "+
			  clipper::String(fileseries_range[i],3)+
			  " which does not exist"));
      }
    }
  }
  //--------------------------------------------------------------

} // end namespace scala

// end hkl_datatypes.cpp
