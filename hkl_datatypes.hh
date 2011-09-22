// hkl_datatypes.hh
//
// Data type definitions
// Just reflection datatypes


#ifndef HKL_DATATYPES_HEADER
#define HKL_DATATYPES_HEADER

// Clipper
#include <clipper/clipper.h>
#include "clipper/core/clipper_precision.h"
using clipper::ftype;
using clipper::Vec3;
using clipper::Mat33;
using clipper::String;
using clipper::Metric_tensor;
typedef clipper::Vec3<double> DVect3;
typedef clipper::Mat33<double> DMat33;
typedef clipper::Vec3<float> FVect3;
typedef clipper::Mat33<float> FMat33;
typedef clipper::Vec3<int> IVect3;

#include "cmtzlib.h"    // CCP4 MTZlib headers (namespace CMtz)
#include "csymlib.h"    // CCP4 symmetry stuff
#include "matvec_utils.hh"  // Matrix & vector utilities
#include "util.hh"

typedef float  Rtype;
typedef double Dtype;

typedef std::pair<int,int> IndexPair;
typedef std::pair<float,float> RPair;
typedef std::pair<double,double> DPair;


namespace scala
{
  enum Chirality {UNKNOWN, CHIRAL, NONCHIRAL, CENTROSYMMETRIC};
  std::string Chiral_as_string(const Chirality& chiral);

  class SpaceGroup; // forward definition

  enum AnomalousClass {ALL, BOTH, IPLUS, IMINUS};

  //--------------------------------------------------------------
  // Orientation things
  //
  // Returns angles in degrees between principle reciprocal axes a*, b*, c*
  // and the rotation spindle (along z)
  // On entry:
  //  DUB is setting matrix [D][U][B]
  DVect3 SpindleToPrincipleAxes(const DMat33& DUB);
  //--------------------------------------------------------------
  // Returns index of smallest component of vector
  int SmallestComponent(const DVect3& v);
  //--------------------------------------------------------------
  // Format closet reciprocal axis to spindle
  std::string FormatSpindleToPrincipleAxes(const DMat33& DUB);
  //--------------------------------------------------------------
  class IntRange;  // forward definition
  //================================================================
  typedef clipper::datatypes::I_sigI<float> ClipperI_sigI;
  class IsigI
  //! I, sigI pair, based on clipper::datatypes::I_sigI but with some differences
  {
  public:
    IsigI() {clipper::Util::set_null(I_); clipper::Util::set_null(sigI_); }
    IsigI(const float& I, const float& sigI) {I_=I; sigI_=sigI;}
    IsigI(const IsigI& Is) {
      if (this != &Is) {I_=Is.I(); sigI_=Is.sigI();}
    }
    // Construct from clipper I_sigI
    IsigI(const ClipperI_sigI& cIsig) {I_=cIsig.I(); sigI_=cIsig.sigI();}

    // Accessors
    const float&    I() const {return  I_;}
    const float& sigI() const {return sigI_;}
    // write access
    float&    I() {return    I_;}
    float& sigI() {return sigI_;}

    void scale(const float& a) {I_ *= a; sigI_ *=a;}
    void scale(const double& a) {I_ *= float(a); sigI_ *=float(a);}

    IsigI scaleIs(const float& a) const {return IsigI(I_*a, sigI_*a);}
    IsigI scaleIs(const double& a) const {return IsigI(I_*a, sigI_*a);}

    bool missing() const { return (clipper::Util::is_nan(I_) || clipper::Util::is_nan(sigI_)); }
    bool missingI() const { return (clipper::Util::is_nan(I_)); }
    bool missingsigI() const { return (clipper::Util::is_nan(sigI_)); }


  private:
    float I_, sigI_;
  };

  //======================================================================
  //--------------------------------------------------------------
  class ReindexOp : public RTop<double>
    //! Reindex operator H  (plus optional translation, usually 0,0,0)
    /*! Applies to index h such that h'T = hT H
      thus H matrices concatenate in order H' = H1 H2 */
  {
  public:
    ReindexOp() : RTop<double>(RTop<double>::identity()),
      strict(false), deviation(0.0) {}
    //! constructor from matrix (and optional vector)
    ReindexOp(const Mat33<double>& Hin,
	      const Vec3<double>& Vin = Vec3<double>(0.0,0.0,0.0))
      : RTop<double>(Hin, Vin),
      strict(false), deviation(0.0) {}
    //! constructor from RTop
    ReindexOp(const RTop<double>& RTin)
      : RTop<double>(RTin),
      strict(false), deviation(0.0) {}
    //! copy constructor
    ReindexOp(const ReindexOp& Rdxin)
      : RTop<double>(Rdxin)    {
      strict = Rdxin.strict;
      deviation = Rdxin.deviation;
    }

    //! constructor from string eg "2h+k,k,l"
    ReindexOp(const std::string& Operator);

    //! store "strict" flag for later return      
    void SetStrict(const bool& Strict) {strict = Strict;}
    //! return stored "strict" flag
    bool Strict() const {return strict;}
    //! store cell deviation for later return, and for sorting on minimum deviation
    void SetDeviation(const double& Deviation) {deviation = Deviation;}
    //! return stored deviation
    double Deviation() const {return deviation;}
    
    //! transpose (of rotation part)
    Mat33<double> transpose() const {return rot().transpose();}
    //! inverse
    ReindexOp inverse() const {return ReindexOp(RTop<double>::inverse());}
    
    //! access as vector 
    std::vector<double> as_Dvec() const {return MVutil::SetVMat33(rot());}
    
    //! returns true if identity    
    bool IsIdentity() const {return MVutil::is_rtop_ident(*this);}
    //! returns true if determinant is positive (desirable!)
    bool Positive() const {
      return (rot().det() > 0.0);}
    //! Returns true if there is a translation component
    bool IsTranslation() const {
      return trn() != clipper::Vec3<double>(0.0,0.0,0.0);}

    //! change basis of symmetry operator
    clipper::Symop Symop(const clipper::Symop& symop) const;

    //! format as string [h,k,l]
    std::string as_hkl() const;
    //! format as matrix
    std::string as_matrix() const;
    //! format as XML matrix
    std::string as_XML() const;
    //! format XML operator
    std::string as_hkl_XML() const;
    //! true if two reindex operators are equal within tolerance
    bool equals(const ReindexOp& b,
		const double& tol = 1.0e-6) const;
    
    //! Apply symmetry operators to reindex operator H and choose the simplest
    void FindSimplest(const SpaceGroup& symm);

    //! for sort on deviation, leaving identity at beginning
    friend bool operator < (const ReindexOp& a,const ReindexOp& b);
    
    friend bool operator == (const ReindexOp& a,const ReindexOp& b);
    friend bool operator != (const ReindexOp& a,const ReindexOp& b);
    
    friend ReindexOp operator * (const ReindexOp& a,const ReindexOp& b);
    friend Vec3<double> operator * (const Vec3<double> v, const ReindexOp& R);
    friend ReindexOp operator * (const RTop<double>& a,const ReindexOp& b);
    friend ReindexOp operator * (const ReindexOp& a,const RTop<double>& b);
    
     
  private:
    bool strict;
    double deviation;
  };
  //--------------------------------------------------------------
  //! Metric tensor  (from Clipper with extensions)
  /*! The metric tensor is used to determine a distance in real or
    reciprocal space using fraction coordinates or Miller indices. It
    is symmetrical, so only the upper triangle is stored with the
    off-diagonal elements doubled. 

    Note that the metric tensor is independent of the
    orthogonalisation convention, although the reciprocal metric tensor
    [M] = [B]T [B]  where [B] is a reciprocal orthogonalisation matrix
  */
  class MetricTensor
  {
  public:
    //! null constructor
    inline MetricTensor() {} //!< Null constructor
    //! constructor: takes parameters of normal or inverse cell
    //  angles in degrees
    MetricTensor( const ftype& a, const ftype& b, const ftype& c, const ftype& alph, const ftype& beta, const ftype& gamm );
    //! constructor from cell dimension vector
    MetricTensor(const std::vector<Dtype> cell);
    //! constructor from matrix
    MetricTensor(const Mat33<Dtype> MetricMatrix);
    //! apply metric to vector
    inline ftype lengthsq( const Vec3<>& v ) const
    { return ( v[0]*(v[0]*m00 + v[1]*m01 + v[2]*m02) +
	       v[1]*(v[1]*m11 + v[2]*m12) + v[2]*(v[2]*m22) ); }
    //! apply metric to int vector
    inline ftype lengthsq( const Vec3<int>& v ) const
    {	ftype h = ftype(v[0]); ftype k = ftype(v[1]); ftype l = ftype(v[2]);
      return h*(h*m00 + k*m01 + l*m02) + k*(k*m11 + l*m12) + l*(l*m22); }

    //! extract cell dimensions (real or reciprocal)
    std::vector<Dtype> cell() const;

    //! invert real <-> reciprocal
    MetricTensor inverse() const;

    //! return as matrix
    Mat33<Dtype> matrix() const;

    String format() const;  //!< return formatted String representation
  private:
    ftype m00, m11, m22, m01, m02, m12;
  };

  //======================================================================
  class Scell
  //! Unit cell class
  /*! This has a different orthogonalisation convention to that in
    clipper::Cell. This one follows the "Cambridge" convention used in
    Mosflm, Scala etc
   */
  {
  public:
    Scell() {cell_.assign(6,0.0);}
    //! constructor from real cell vector
    //  derive reciprocal cell and orthogonalisation matrix B
    Scell(const std::vector<Dtype>& rcell);
    //! constructor from real cell vector
    Scell(const float* cell);
    //! construct from real [default] or reciprocal metric tensor
    Scell(const MetricTensor& metric_tensor, const bool real=true);
    //! construct from clipper::Cell
    Scell(const clipper::Cell& ccell);
    //! construct from real cell
    Scell(const double& a,const double& b,const double& c,
	  const double& alpha, const double& beta, const double& gamma);
    //! initialise from cell vector
    void init(const std::vector<Dtype>& rcell);

    std::string formatPrint(const bool& newline=true) const;  //!< simple print
    void dump() const;   //!< fuller debug print

    //! return orthogonalisation matrix [B]
    inline Mat33<Dtype> Bmat() const {return Bmat_;}
    //! return reciprocal metric tensor
    inline MetricTensor rec_metric_tensor() const
    {return recip_metric_tensor_;}

    //! return real cell as a vector
    inline std::vector<Dtype>  UnitCell() const {return cell_;}
    //! return reciprocal cell as a vector
    inline std::vector<Dtype>  ReciprocalCell() const {return rec_cell_;}
    //! return cell element
    inline const Dtype& operator[](const int& i) const {return cell_[i];}
    //! change basis: apply reindex operator
    /*! if [H] is reindex operator (ignoring translations)
      [B] current orthogonalisation matrix
      new [UB]' = [B] [H]T^-1
      new reciprocal metric tensor [B']T[B] = [UB]'T [UB]' -> new cell
    */
    Scell change_basis(const ReindexOp& reindex_op) const;

    //! format cell
    std::string format(const int w=7, const int p=4) const;
    //! return XML representation
    std::string xml() const;
    //! return clipper::Cell
    clipper::Cell ClipperCell() const
    {return clipper::Cell(clipper::Cell_descr
	  (cell_[0],cell_[1],cell_[2],cell_[3],cell_[4],cell_[5]));}
    //! return cell volume
    double Volume() const;

    //! return true if all angles are within Tol of given test values
    bool AngleTest(const double& AlphaTest,
		   const double& BetaTest,
		   const double& GammaTest,
		   const double Tol=0.2) const;
    //! returns true if cells agree within "tolerances" in A & degrees
    bool equals(const Scell& other, const double& tolA=0.02, const double& tolD=0.02) const;
    //! returns true if cells agree within "tolerances" in degrees (lengths & angles)
    bool equalsTol(const Scell& other, const double& AngularTolerance) const;

    //! "difference" between two cells in A, maximum allowed distance
    double Difference(const Scell& other) const;

    //! returns true if any value is close to zero
    bool null() const;

  private:
    std::vector<Dtype> cell_;
    std::vector<Dtype> rec_cell_;
    Mat33<Dtype> Bmat_;
    MetricTensor  metric_tensor_;
    MetricTensor  recip_metric_tensor_;
  };
  //======================================================================
  class UnitCellSet {
    //! A group of unit cells, to averages and deviations etc
  public:
    UnitCellSet(){}
    //! Construct from list of cells
    UnitCellSet(const std::vector<Scell>& Cells);
    //! Initialise from list of cells
    void init(const std::vector<Scell>& Cells);

    //! Add in a cell
    void AddCell(const Scell& Cell);

    //! Number of cells stored
    int Number() const {return cells.size();}

    //! return all cells stored
    std::vector<Scell> Cells() const {return cells;}

    // return average cell
    Scell AverageCell() const {return averagecell;}

    //! return list of deviations (A) from average of other cells
    std::vector<double> Deviations() const;
    //! return worst deviation (A)
    double WorstDeviation() const;
    //! return list of rms deviations for each cell parameter
    std::vector<double> RmsD() const;

    //! Average list of cells, if idxexclude >= 0, exclude entry with this index
    Scell Average(const int& idxexclude=-1) const;

  private:
    std::vector<Scell> cells;
    Scell averagecell;

  }; // UnitCellSet
  //======================================================================
  //! reflection 'Miller' index
  /*! Copied from Clipper::HKL with some simplifications */
  class Hkl : public Vec3<int>
  {
  public:
    inline Hkl() {}                //!< null constructor
    inline explicit Hkl( const Vec3<int>& v ) :
      Vec3<int>( v ) {}            //!< constructor: copy/convert
    inline Hkl( const int& h, const int& k, const int& l ) :
      Vec3<int>( h, k, l ) {}      //!< constructor: from H,K,L
    inline Hkl( const Vec3<double>& d ) :
      Vec3<int>(Nint(d[0]), Nint(d[1]), Nint(d[2])) {} //!< constructor from double
    inline Hkl( const clipper::HKL& chkl) :
      Vec3<int>(chkl.h(),chkl.k(),chkl.l()) {} //!< constructor from clipper::HKL
    inline const int& h() const { return (*this)[0]; }  //!< get h
    inline const int& k() const { return (*this)[1]; }  //!< get k
    inline const int& l() const { return (*this)[2]; }  //!< get l
    inline int& h() { return (*this)[0]; }  //!< set h
    inline int& k() { return (*this)[1]; }  //!< set k
    inline int& l() { return (*this)[2]; }  //!< set l
    bool IsGeneral() const; //!< true if not in any potential zero level
    //! return clipper HKL
    clipper::HKL HKL() const {return clipper::HKL(*this);}

    //! return real version of hkl
    inline Vec3<Dtype> real() const
    {return Vec3<Dtype>(Dtype((*this)[0]),Dtype((*this)[1]),Dtype((*this)[2]));}
    //! return inverse resolution squared for this reflection in given cell
    inline Dtype invresolsq( const Scell& cell ) const
    {return cell.rec_metric_tensor().lengthsq(real());}
    //! change basis: new h' = h * reindex_op 
    Hkl change_basis(const ReindexOp& reindex_op) const;
    //! change basis: new h' = h * reindex_op. return false if indices are non-integral
    bool change_basis(Hkl& Newhkl, const ReindexOp& reindex_op) const;
    //! return orthogonalised vector [B] h
    Vec3<Dtype> orth(const Scell& cell) const
    {return (cell.Bmat() * real());}
    std::string format() const;  //!< return formatted String representation
    //! returned packed form as index code
    //  Note that this uses 10-bit packing, so is not guaranteed to
    //  produce a unique code, but it is good enough for some purposes
    int code() const;

    friend inline Hkl operator -(const Hkl& h1)
    { return Hkl( -h1.h(), -h1.k(), -h1.l() ); }
    friend inline Hkl operator +(const Hkl& h1, const Hkl& h2)
    { return Hkl( h1.h()+h2.h(), h1.k()+h2.k(), h1.l()+h2.l() ); }
    friend inline Hkl operator -(const Hkl& h1, const Hkl& h2)
    { return Hkl( h1.h()-h2.h(), h1.k()-h2.k(), h1.l()-h2.l() ); }
    friend inline Hkl operator *(const int& s, const Hkl& h1)
    { return Hkl( s*h1.h(), s*h1.k(), s*h1.l() ); }
  };
  //--------------------------------------------------------------
  Hkl HklDecode(const int& kode);

  //======================================================================
  class PxdName
  //! Project name / Crystal name / Dataset name
  {
  public:
    PxdName() : pname_(""), xname_(""), dname_("") {}
    //! constructor from names
    PxdName(const std::string& pname_in,
	    const std::string& xname_in,
	    const std::string& dname_in)
      :  pname_(pname_in), xname_(xname_in), dname_(dname_in) {}
  
  
    const std::string& pname() const {return pname_;} //!< get project name
    const std::string& xname() const {return xname_;} //!< get crystalname
    const std::string& dname() const {return dname_;} //!< get dataset name
    
    std::string& pname() {return pname_;} //!< set project name
    std::string& xname() {return xname_;} //!< set crystalname 
    std::string& dname() {return dname_;} //!< set dataset name

    //! Replace elements by any non-blank elements in NewPxd
    void update(const PxdName& NewPxd);
  
    std::string formatPrint() const; //!< format for printing
    std::string format() const; //!< simple format pname/xname/dname

    //! equality
    friend inline bool operator ==(const PxdName& pxd1, const PxdName& pxd2)
    {
      return (pxd1.pname()==pxd2.pname() &&
	      pxd1.xname()==pxd2.xname() &&
	      pxd1.dname()==pxd2.dname());
    }
    //! true if blank
    bool is_blank() const
    {
      return (pname_=="" &&
	      xname_=="" &&
	      dname_=="");
    }

  private:
    std::string pname_, xname_, dname_;
  };
  //======================================================================
  class Xdataset
  //! Crystal/Dataset object
  /*!  for these purposes we are not really interested in the crystal
    level in the hierarchy, just datasets */
  {
  public:
    Xdataset(){}
    //! constructor from names, cell, wavelength, ID index
    Xdataset(const PxdName& pxdname, const Scell& cell,
	     const float& wavel, const int& setid);
  
    //! Add batch number to list for this dataset
    void add_batch(const int& batch_num);
    void ClearBatchList() {batches.clear();} //!< clear batch list
    int num_batches() const {return batches.size();} //!< number of batches
    //! change basis: reindex to get new cell
    void change_basis(const ReindexOp& reindex_op);

    void AddRunIndex(const int& RunIndex); //!< add to run index list
    void ClearRunList() {run_index_list.clear();} //!< clear run index list
    std::vector<int> RunIndexList() const {return run_index_list;} //!< return run index list
  
    PxdName pxdname() const {return pxdname_;} //!< return PXD names
    int setid() const {return setid_;} //!< get set ID index
    int& setid() {return setid_;}  //!< set set ID index

    Scell cell() const {return cell_;} //!< return cell
    Scell& cell() {return cell_;} //!< set cell
    float wavelength() const {return wavel_;} //!< return wavelength
    float& wavelength() {return wavel_;} //!< set wavelength

    float Mosaicity() const {return av_mosaic;} //!< return average mosaicity
    float& Mosaicity() {return av_mosaic;}  //!< set average mosaicity
   
    std::string formatPrint() const; //!< format

    //! add in another cell and wavelength
    void AddCellWavelength(const Scell& newcell, const float& wavel);

    //! average multiple cells and wavelengths
    void AverageCellWavelength(); 
    //! return number of cells/wavelengths
    int NumberofCells() const {return allcells_.Number();}

    std::string formatAllCells() const; //!< format cell & wavelength list if more than one
    //! return worst deviation (A), = 0 if only one
    double WorstDeviation() const;

    //! equality, just tests pxdname
    friend bool operator == (const Xdataset& a,const Xdataset& b);
    
  private:
    PxdName pxdname_;
    int setid_;
    Scell cell_;
    float wavel_;
    float av_mosaic;
    // List of all batch numbers belonging to this dataset
    std::vector<int> batches;
    // List of runs (indices)
    std::vector<int> run_index_list;
    // List of unit cells if multiple runs
    UnitCellSet allcells_;
    // List of wavelengths if multiple runs
    std::vector<float> allwavel_;
  };
  //======================================================================
  class BatchNumber
  //! A batch number and its original version (before offset, if any)
  {
  public:
    BatchNumber() : number(0), original_number(0) {}
    //! constructor from number == original number
    BatchNumber(const int& n) : number(n), original_number(n) {}
    //! constructor from number & original number
    BatchNumber(const int& n, const int& norig):  number(n), original_number(norig) {}

    int Number() const {return number;}  //!< return batch number
    int OriginalNumber() const {return original_number;} //!< return original number

  private:
    int number;
    int original_number;
  };
  //======================================================================
  class Batch
  //! Batch object
  /*! This class encapsulates all the batch information from the MTZ header*/
  {
  public:
    Batch();  // a dummy batch
    //! constructor from MTZ batch, accept flag, dataset index
    Batch(const CMtz::MTZBAT& batch, const bool& accept, const int& idataset);
    //! change basis, reindex cell and orientation
    void change_basis(const ReindexOp& reindex_op);
    //! store cell constraint flags
    void SetCellConstraint(const std::vector<int>& lbcell);
    //! store cell
    void SetCell(const Scell& newcell);
    //! set file serial number  [default = 1]
    int& FileNumber() {return file_num;}
    //! get file serial number
    int  FileNumber() const {return file_num;}

    //! offset the stored phi1 & phi2, to allow to 360deg wrap-around
    void OffsetPhi(const float& Phioffset);
    float PhiOffset() const {return phioffset;} //!< return phi offset
    void SetRunIndex(const int& RunIndex); //!< store run index
    int RunIndex() const {return run_index;} //!< return run index

    int num() const {return batchinfo.num;} //!< get batch number
    int& num() {return batchinfo.num;}  //!< set batch number

    //! offset batch number & store offset [default = 0]
    void OffsetNum(const int& Offset);
    //! retrieve batch number offset [default = 0]
    int  BatchNumberOffset() const {return offset;}  // get

    int index() const {return Xdataset_index;} //!< get dataset index
    int& index() {return Xdataset_index;} //!< set dataset index
    int DatasetID() const {return batchinfo.nbsetid;} //!< return dataset ID
    int& DatasetID() {return batchinfo.nbsetid;}  //!< set dataset ID

    //! return start phi1
    float Phi1() const {return valid_phi ? batchinfo.phistt : 0.0f;}
    //! return end phi2
    float Phi2() const {return valid_phi ? batchinfo.phiend : 0.0f;}
    //! return mid-phi 
    float MidPhi() const {return valid_phi ?  0.5*(batchinfo.phistt+batchinfo.phiend) : 0.0f;}
    bool ValidPhi() const {return valid_phi;} //!< true if valid phi
    float PhiRange() const
    {return Phi2() - Phi1();} //!< return phi range = phi2 - phi1

    // valid_time: > 0 valid time information (time1, time2)
    //             < 0 time inferred from phi
    //              -1 set = phi
    //              -2 inverted from descending phi
    //             = 0 no time information
    //! return time1 if valid else 0
    float Time1() const {return (valid_time!=0) ? batchinfo.time1 : 0.0f;}
    //! return time2 if valid else 0
    float Time2() const {return (valid_time!=0) ? batchinfo.time2 : 0.0f;}
    //! return mid-time if valid else 0
    float MidTime() const {return (valid_time!=0) ? 0.5*(batchinfo.time1+batchinfo.time2) : 0.0f;}
    //! true if either time column present or phi as proxy
    bool IsValidTime() const {return (valid_time!=0);}
    //! true if phi is proxy for time
    bool IsTimePhi() const {return (valid_time<0);}
    //! store flag for time column present in input file
    void SetTimeColumnPresent() {valid_time = +1;}
    //! store flag for time column absent from input file
    void SetTimeColumnAbsent() {valid_time = 0;}
    //! store flag for phi as proxy for time
    void SetPhiAsProxyTime() {valid_time = -1;}
    //! store flag for time value taken from phi has been reset to make ascending 
    void SetTimeReversedFromPhi() {valid_time = -2;}
    //! offset stored time1, time2 to allow to 360deg wrap-around if time is copied from Phi
    void OffsetTime(const float& Timeoffset);
    //! Store time range limits
    void StoreTimeRange(const float& Time1, const float& Time2) {
      batchinfo.time1 =  Time1;
      batchinfo.time2 =  Time2;
    }


    bool Accepted() const {return accepted;}  //!< return accepted flag
    void SetAccept(const bool& accept) {accepted = accept;} //!< set accepted flag

    Scell cell() const {return bcell;} //!< return cell
    float Mosaicity() const {return batchinfo.crydat[0];} //!< return mosaicity
    float Wavelength() const {return batchinfo.alambd;} //!< return wavelength

    //! true if we have valid orientation information
    bool ValidOrientation() const {return valid_Umat;}
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

    //! reciprocal lattice index h -> s(r) = [R][D][U][B]h camera frame
    DVect3 HtoSr(const Hkl& hkl, const float& phi) const;
    //! reciprocal lattice index h -> s(r0) = [D][U][B]h zero rotation angle frame
    DVect3 HtoSr0(const Hkl& hkl) const;
    //! camera frame s(r) -> s(r0) = [R]^-1 s(r) zero rotation angle frame
    DVect3 SrtoSr0(const DVect3& sr, const float& phi) const;
    //! zero rotation angle frame s(r0) -> polar orthogonalised crystal frame
    /*!    s(pole) = [P] [DU]^-1 s(r0) */
    DVect3 Sr0toP(const DVect3& sr0) const;

    DMat33 Umat() const {return U;} //!< return orientation matrix [U]

    //! return reciprocal axis (name) closest to spindle (a*, b*, c*)
    std::string SpindleToPrincipleAxis() const;

    //! return source vector s0
    DVect3 Source() const {return s0;}
    //! true if 3-axis phi scan (ie not omega scan)
    bool PhiScan() const {return phiscan;}

    //! set absorption pole permutation, 1,2,3 for a*,b*,c*, <0 automatic
    /*! = 0 no pole
      < 0 choose closest to spindle */
    void SetPole(const int& Pole);
    //! return absorption pole, = 0 unset
    int Pole() const {return pole;}

    //! Return range of detector pixel coordinates
    std::vector<std::vector<float> > DetectorCoordinateRange() const;

    //! return number of detectors
    int Ndet() const {return batchinfo.ndet;}

    //! return MTZ batch: only use this for writing MTZ file
    /*! otherwise information retrieval should be through explicit calls */
    CMtz::MTZBAT batchdata() const {return batchinfo;}

    //! return formatted version as from CCP4 mtzlib
    std::string format() const;
    //! print using CCP4 mtzlib routine
    void print() const;

    //! for sorting on batch number
    friend bool operator < (const Batch& a,const Batch& b);

  private:
    // initialise batchinfo to zeroes
    void initBatchInfo();

    CMtz::MTZBAT batchinfo;   // All in original frame
    int Xdataset_index;      // index into crystal/dataset list (-1 if rejected)
    int run_index;           // index into run list
    bool accepted;
    bool valid_cell;
    bool valid_Umat;
    // valid_time: > 0 valid time information (time1, time2)
    //             < 0 time inferred from phi
    //              -1 set = phi
    //              -2 inverted from descending phi
    //             = 0 no time information
    int valid_time;
    bool valid_phi;  // valid phi information (phi1, phi2)
    float phioffset;   // Offset in phi from original
    float timeoffset;   // Offset in time from original
    int offset;      // number offset from original
    int file_num;    // serial number of file from which this batch came

    DVect3 spindle; // goniostat vector for principle rotation axis (= e1 or e3)

    // Stored information, in Cambridge frame
    //  Cambridge frame has
    //    source vector s0 (anti-parallel to beam) approximately along -x
    //      ie x is along beam (if perpendicular to rotation axis)
    //    principal rotation axis e0 exactly along z
    DMat33 Q;  // conversion matrix from input frame to Cambridge
    DVect3 s0;  // source vector
    DMat33 U;  // [U]     [Orientation]
    DMat33 DU; // [D][U]  [Datum][Orientation]
    DMat33 DUB; // [D][U][B]  [Datum][Orientation][Orthogonalisation]
		   

    int pole;      // absorption pole, 0 if unset
    DMat33 PDUinv; // [P][[D][U]]^-1 for crystal frame absorption pole
    bool phiscan;  // true for Phi scan in 3-axis system, else false
    DMat33 E1E2;   // [[Omega][Chi/Kappa]] matrix for Phi scan in
		    // 3-axis system, else [I]

    Scell bcell;
  };
  //======================================================================
  // Return true if dataset setid is in datasets list
  //  & return dataset index idataset (-1 if not)
  bool in_datasets(const int& setid,
		   const std::vector<Xdataset>& datasets,
		   int& idataset);
  //======================================================================
  class BatchSelection
  //! A selection of batches, ie a collection of batch ranges & a batch list
  /*! Batches selections may be flagged with a "file series number":
    if this is > 0, then the numbers apply to the original numbers in the
    file, before any offsets. If == 0 then it applies to the offset numbers.
    They may also be flagged with an  numerical flag eg a run number indicating
    what "class" they belong to, > 0 if valid, = 0 if unset
 */
  {
  public:
    BatchSelection();
    void clear(); //!< clear list

    //! return true if no entries
    bool Null() const {return (batchlist.size()==0 && flaglist.size()==0);}

    //! Add in one batch to selection, qualified by fileSeriesList
    /*! fileSeriesList = 0 means final batch numbering (after any
      renumbering for file series > 1) */
    void AddBatch(const int& batch, const int& fileSeriesList);
    //! Add in batch range to selection, qualified by fileSeriesList
    /*! fileSeriesList = 0 to mean final batch numbering (after any
      renumbering for file series > 1) */
    void AddRange(const int& batch1, const int& batch2,
		  const int& fileSeriesList, const int& flag=-1);

    //! Return index in list if in selection, for given file series, else -1
    int FindInSelection(const int& batch,
			const int& fileSeriesTest) const;

    //! Return true if in selection, for given file series
    bool InSelection(const int& batch, const int& fileSeriesTest) const;

    //!< Two possibilities for selection:
    //!<  1) Specified selection on final numbering, ie from sole file or
    //!<  after renumbering of 2nd or subsequent file, fileSeriesList == 0,
    //!<  only test if fileSeriesTest == 0
    //<!
    //<!  2) Specified selection on original file numbering from 2nd or
    //<!  subsequent file, fileSeriesList > 0, only test if
    //<!  fileSeriesTest > 0

    //! return numerical flag for batch, = -1 if not in list
    int FlagNumber(const int& batch);

    //! return list of unique numerical flags specified
    std::vector<int> UniqueFlags() const;

    //! return list of batch ranges for specified flag, all ranges if flag < 0
    std::vector<IntRange> BatchRanges(const int& flag);

    //! Check that fileSeries specified on selection commands match specified files. Fail here if not
    void CheckSeries(const int& NumFileSeries) const;

  private:
    std::vector<int> batchlist;         // list of batch number
    std::vector<int> fileseries_list;   // file series numbers for list 
    std::vector<IntRange> batchranges;  // batch ranges
    std::vector<int> fileseries_range;  // file series numbers for ranges
    std::vector<int> flaglist;          // a numerical flag eg a run number
  };
//--------------------------------------------------------------
  class FileRead
  //! Status of a file, was it read succesfully etc?
  {
  public:
    FileRead() :read(false), nexcluded(0) {}
    //! constructor from flags
    /*! \param Opened  true if successfully opened
        \param Read    true if successfully read
        \param SymmetrySet true if file contained symmetry information
        \param Nexcluded   number of observation(parts) excluded
    */
    FileRead(const bool& Opened, const bool& Read, 
	     const bool& SymmetrySet, const int& Nexcluded)
      : opened(Opened), read(Read),
	nexcluded(Nexcluded), symmetrySet(SymmetrySet)	 {}
    
    //! return true if successfully opened
    bool Opened() const {return opened;}
    //! return true if file was read				    
    bool Read() const {return read;}
    //! return true if symmetry set from file
    bool SymmetrySet() const {return symmetrySet;}
    //! return number excluded
    int Nexcluded() const {return nexcluded;}

  private:
    bool opened;        // true if successfully opened
    bool read;          // true if file was read
    int  nexcluded;     // number of observation(parts) excluded
    bool symmetrySet;   // true if the file contained symmetry information
			// (eg lattice type)
  };
}  // namespace scala
#endif
// end  hkl_datatypes.hh
