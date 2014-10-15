// Tile correction


#ifndef TILE_HEADER
#define TILE_HEADER

// Clipper
#include <clipper/clipper.h>

#include "range.hh"
#include "fileread.hh"
#include "hkl_datatypes.hh"
#include "tie.hh"
#include "scala_util.hh"

namespace scala {
  class TileBase;
  //--------------------------------------------------------------
  class DetectorType {
    //! Detector type and characteristics
  public:
    enum Type {UNKNOWN, CCD1, CCD2x2, CCD3x3, PILATUS6M, PILATUS2M, PIXEL};

    DetectorType() :dettype(UNKNOWN), typestr("Unknown"),
		    ndet(0), ntilex(0), ntiley(0) {}

    //! construct from Batch object
    DetectorType(const Batch& batch);

    //! construct from arguments
    DetectorType(const Type& Dtype, const std::string& TypeLabel,
    		 const std::vector<std::vector<float> >& Detrange);

    //! construct from arguments
    DetectorType(const std::string& TypeLabel,
		 const std::vector<std::vector<float> >& Detrange);

    //! valid: non-zero detector coordinate range
    bool Valid() const;

    //! Detector range on Xdet
    Range XdetRange() const;
    //! Detector range on Ydet
    Range YdetRange() const;

    //! Number of tiles deduced from type
    int NtileX() const {return ntilex;}

    //! Number of tiles deduced from type
    int NtileY() const {return ntiley;}

    Type type() const {return dettype;}

    //! return a string corresponding to the detector type
    std::string TypeLabel() const;

    //! return a string corresponding to the detector type
    std::string TypeLabel(const DetectorType::Type& dtype) const;

    //! return a detector type corresponding to the string
    Type TypeFromLabel(const std::string& typelabel) const;

    //! test for equality on type and detector range (not on number of tiles)
    bool equals(const DetectorType& b) const;

    friend bool operator == (const DetectorType& a,const DetectorType& b);
    friend bool operator != (const DetectorType& a,const DetectorType& b);

  private:
    Type dettype;
    std::string typestr; // detector type (if we can deduce it!)
    int ndet;            // number of detectors included in this class
    std::vector<std::vector<float> > detrange;  // coordinate range
    int ntilex, ntiley;  // number of tiles on Xdet, Ydet
    bool pixelcoords;    // true if coordinates seem to be pixels
  };
  //--------------------------------------------------------------
  class DetectorStatistics;
  class hkl_unmerge_list;

  class DetectorAnalysis {
    //! Analysis of detector scales etc, for multiple detectors if necessary
  public:
    DetectorAnalysis(){}
    //! Construct to match reflection list
    DetectorAnalysis(const hkl_unmerge_list& hkl_list);
    //! Initialise to match reflection list
    void init(const hkl_unmerge_list& hkl_list);

    //! Construct explicitly for testing
    DetectorAnalysis(const std::vector<std::vector<float> >& Detrange);

    void AddStats(const float& I, const float& AvI,
		  const int& runidx,
		  const int& xdet, const int& ydet);

    void WriteImages(const std::string& fname) const;

  private:
    int ndet;  // number of different detectors
    std::vector<DetectorStatistics> detectorstatistics; // for each detector
    std::vector<int> idxrun; // index into detectorstatistics array for each run
  };
  //--------------------------------------------------------------
  class DetectorStatistics {
    //! Analysis of detector scales etc, for one detector
  public:
    DetectorStatistics(){}
    //! construct from detector type object
    DetectorStatistics(const DetectorType& Dtype);
    //! initialise from detector type object
    void init(const DetectorType& Dtype);

    //! Set pixel binning, default = 8
    void SetBinning(const int& Ngpxlx, const int& Ngpxly);

    void AddStats(const float& I, const float& AvI,
		  const int& xdet, const int& ydet);

    // Scale factors for each pixel group
    void WriteImageScales(const std::string& fname) const;
    // Deviations for each pixel group
    void WriteImageDeviations(const std::string& fname) const;

  private:
    DetectorType detectortype;
    int ngpxlX, ngpxlY; // pixel binning
    Range xdetrange;
    Range ydetrange;
    int nbx, nby;       // number of bins along x, y

    clipper::Array2d<double> sumwIxy;     // Sum w Ixy for position xy
    clipper::Array2d<double> sumwIothers; // Sum w Iothers for position xy
    clipper::Array2d<MeanSD> avdelta;     // (I-AvI) for position xy

  };
  //--------------------------------------------------------------
  class DetectorScale
  {
    // Scale factor for each position on the detector
    // Parameterised by tile
  public:
    // scale type for each tile
    enum DetectorScaleType {NONE, FLAT, CCD1, CCD2, CCD3,
			    PIXEL, AUTOMATIC};
    DetectorScale() : detectorscaletype(NONE), ntilex(0), ntiley(0) {}
    DetectorScale(const DetectorScaleType& DetScaleType,
		  const int& nTileX, const int& nTileY,
		  const DetectorType& dettype);

    void init(const DetectorScaleType& DetScaleType,
	      const int& nTileX, const int& nTileY,
	      const DetectorType& dettype);

    void init(const DetectorScaleType& DetScaleType,
	      const int& nTileX, const int& nTileY,
	      const Range& Xrange, const Range& Yrange);

    void init();

    // reassign parameters and tilescales objects after Restore
    void reinit();

    void setSymmetric(const bool& symmetric=false);

    ~DetectorScale();

    // Copy & copy constructor should fail or be done properly due to pointers
    DetectorScale(const DetectorScale& detscale);
    DetectorScale& operator= (const DetectorScale& detscale);

    DetectorScaleType Type() const {return detectorscaletype;}

    //! Return true if Valid
    bool Valid() const;

    // Store parameter vector (length nparams)
    void StoreParameters(const std::vector<double>& parameters);
    // Retrieve parameter vector (length nparams)
    std::vector<double> Parameters() const;

    //! lower bounds for ipar'th parameter
    bool LowerBound(const int& ipar, double& Lower) const;
    //! upper bounds for ipar'th parameter
    bool UpperBound(const int& ipar, double& Upper) const;
    //! "large shift" for ipar'th parameter
    double LargeShift(const int& ipar) const;

    // Store Nobs vector (length nparams)
    void StoreNobservations(const std::vector<int>& Nobs);
    // Retrieve number of contributions (length nparams)
    std::vector<int> Nobservations() const {return nobsPar;}

    // clear all counts, eg tile corners
    void clearCounts();

    //! number of parameters
    int Number() const {return nparams;}

    //! return number of parameters/tile
    int NparamsTile() const;

    //! Return list of ties: sdties are sds for weight, idx0 is index to first global
    // parameter for setting ties, since they refer to the global parameter index
    //
    // For CCD tiles:
    //  sdties[0] for r
    //  sdties[1] for w
    //  sdties[2] for A
    //  sdties[3] for x0, y0
    //  sdties[4] for Fourier coefficients
    std::vector<Tie> Ties(const std::vector<double> sdties,
			  const int& idx0) const;

    // Return scale for detector coordinates Xdet, Ydet
    double Scale(const std::pair<float,float>& XYdet) const;

    // Return scale & derivatives for detector coordinates Xdet, Ydet
    double ScaleDeriv(const std::pair<float,float>& XYdet,
		      std::vector<double>& dgdp) const;

    // Return scale & derivatives for  detector coordinates Xdet, Ydet
    void ScaleDeriv(const std::pair<float,float>& XYdet,
		    double& scale, std::vector<double>& dgdp) const;

    std::string format() const;

    std::string formatparameters(const std::vector<double>& sds) const;

    std::string formatTies() const;

    //! write output image of correction factors
    // as ADSC format image
    void WriteImage(const std::string& imagefilename) const;

    // Format all information into a labelled save format for later restoration
    std::string FormatSave() const;

    // restore
    void Restore(Fileread& FR);

    static std::string formatType(const  DetectorScaleType& type);
    static DetectorScaleType Type(const std::string& scaletypelabel);


  private:
    DetectorType type;
    DetectorScaleType detectorscaletype;  // type of correction
    int ntilex, ntiley;  // number of tiles on Xdet, Ydet
    int nparams;         // total number of parameters
    Range xdrange, ydrange;       // detector coordinate ranges

    clipper::Array2d<TileBase*> tilescales; // the tiles
    clipper::Array2d<int> idx_tile;         // index to first parameter for each tile
    std::vector<int> nobsPar;  //  number of observations for each parameter

    // Ties: note that some parameters between different tiles may be
    // restrained together

    //! find which tile ipar'th parameter belongs to
    std::pair<int,int> WhichTile(const int& ipar) const;

    // Return scale & derivatives for  detector coordinates Xdet, Ydet
    void ScaleDeriv(const bool& Deriv,
		    const std::pair<float,float>& XYdet,
		    double& scale, std::vector<double>& dgdp) const;
  }; // class DetectorScale
  //--------------------------------------------------------------
  class RadialFunctionBase {
  public:
    RadialFunctionBase(){};

    //! calculate the value of the function and store intermediates 
    virtual double value(const double& z, const double& A) = 0;

    //! derivatives, must follow a call to value() 
    virtual std::vector<double> deriv(const int& nparams,
				      const double& w) const = 0;

  private:
  };
  //--------------------------------------------------------------
  class RadialFunctionGompertzCDF : public RadialFunctionBase{
    // suggested by Rob Nicolls, but not apparently better than erfc
  public:
    RadialFunctionGompertzCDF() : A_(-1000.0) {} 

    //! calculate the value of the function and store intermediates 
    double value(const double& z, const double& A);

    //! derivatives, must follow a call to value() 
    std::vector<double> deriv(const int& nparams,
			      const double& w) const;

  private:
    double z_;
    double A_;
    // scale = A f(z) + 1 - A
    // f(z) = 1 - exp(g(z))
    // g(z) = -exp(-z)-1
    double expmz;
    double gz;
    double fz;
  };
  //--------------------------------------------------------------
  class RadialFunctionErfc : public RadialFunctionBase{
  public:
    RadialFunctionErfc() : A_(-1.0) {} 

    //! calculate the value of the function and store intermediates 
    double value(const double& z, const double& A);

    //! derivatives, must follow a call to value() 
    std::vector<double> deriv(const int& nparams,
			      const double& w) const;

  private:
    double z_;
    double A_;
    // scale = A f(z) + 1 - A
    // f(z) = 0.5 * erfc(z)
    double gz;
    double fz;
  };
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  class TileBase {  // base class for tile scale
  public:
    TileBase(){}
    TileBase(const double& Xmax, const double& Ymax) {}

    virtual ~TileBase() {}

    virtual void init(const double& Xmax, const double& Ymax) = 0;

    virtual void clearCounts() = 0;

    //! symmetric = false to allow A to vary around the tile
    virtual void setSymmetric(const bool& symmetric) {} // dummy

    // Pure virtual functions
    virtual int Nparams() const = 0; //!< number of parameters
    //! Store parameters
    virtual void StoreParameters(const std::vector<double>& parameters) = 0;
    //! Retrieve parameter vector (length nparams)
    virtual std::vector<double> Parameters() const = 0;
    //! return vector of ties, given SDs and 1st global parameter index
    virtual std::vector<Tie> Ties(const std::vector<double> sdties,
				   const int& idx0)
    {return std::vector<Tie>();} // default return empty vector

    //! vector of indices and weights for each parameter to be restrained across tiles 
    virtual std::pair<std::vector<int>, std::vector<double> > 
    TiedParameters(const std::vector<double> sdties,
		   const int& idx0)
    {return std::pair<std::vector<int>, std::vector<double> >();}


    //! lower bounds for ipar'th parameter
    virtual double LowerBound(const int& ipar) const = 0;
    //! upper bounds for ipar'th parameter
    virtual double UpperBound(const int& ipar) const = 0;
    //! "large shift" for ipar'th parameter
    virtual double LargeShift(const int& ipar) const = 0;

    //! Return edge radius
    virtual double Rad0() const {return 0.0;}
    //! Return Xcentre
    virtual double Xcentre() const {return 0.0;};
    //! Return Ycentre
    virtual double Ycentre() const {return 0.0;};

    // Return scale & derivatives for tile coordinates Xt, Yt
    virtual void ScaleDeriv(const bool& Deriv,
			    const double& Xt, const double& Yt,
			    double& scale, std::vector<double>& dgdp) const = 0;

    // Format layout for printing
    virtual std::string format() const = 0;
    // Format formulation for printing
    virtual std::string formattype() const = 0;

    //! format parameters for printing
    virtual std::vector<std::string> formatparameters(const std::vector<double>& sds) const
    {return std::vector<std::string>(1,"");}

    virtual std::string formatTies() const
    {return "";}

    // Format all information into a labelled save format for later restoration
    virtual std::string FormatSave() const {return "";}

    // restore
    virtual void Restore(Fileread& FR) {}

    //! number of smoothing parameters for each of r,w,A, CCD only
    virtual int NparamsSmooth() const {return 0;}

    //! record grid coordinates
    void SetGridCoordinates(const int& kx, const int& ky)
    {ix = kx; iy=ky;}

    //! record grid coordinates
    std::pair<int,int> GridCoordinates() const
    {return std::pair<int,int>(ix, iy);}


  protected:
    int ix, iy;  // grid coordinates of this tile
  };

//==================================================================
class FourierSmooth {
  // A four- or five parameter Fourier class,
  //    like Hendrickson-Lattman coefficients
  // No constant term if 4 parameters
  // for angle p,
  //   v = A cos(p) + B sin(p) + C cos(2p) + D sin(2p)
  // ie the 1st two complex Fourier coefficients
  //
public:
  FourierSmooth():nparams(-1){}
  //! construct or initialise from constant value
  FourierSmooth(const double& flatlevel, const bool& isconstant=true);
  //! construct or initialise from parameters
  //  FourierSmooth(const std::vector<double> parameters)
  //  {setParameters(parameters);}

  //! set either: true for 4 parameters; false no constant term E, 4 params
  void setIsConstant(const bool& isconstant);

  //! set 4 or 5 parameters
  void setParameters(const std::vector<double> parameters);

  //! set a constant level, ie set E, A=B=C=D=0
  void setLevel(const double& flatlevel=0.0);

  //! return parameters
  std::vector<double> GetParameters() const {return parameters;}

  //! get value at angle phi
  double Value(const double& phi) const;

  //! get value at angle phi and dvdp its derivatives wrt parameters
  double ValueDerivatives(const double& phi,
			  std::vector<double>& dvdp) const;

  std::string format() const;
  std::string dump() const;

  int NumberParameters() {return nparams;}

private:
  int nparams;
  std::vector<double> parameters;
};
  // Note there are 3 CCDTile classes, default is CCD2
  //  CCD1  r, w, A all radially symmetric
  //  CCD2  just A varies with azimuthal angle around tile centre
  //         (unless circularlysymmetric is set, default off)
  //  CCD3  r,w, and A all vary with azimuthal angle around tile centre
  //         (unless circularlysymmetric is set, default off)
  //--------------------------------------------------------------
  class CCDTile3 : public TileBase {
    //! A CCD tile, correct for fall-off in the corners
    // Coordinates within the tile are defined from 0->xmax, 0->ymax

    // The model:
    //  distance from centre d = sqrt((X-x0)^2 + Y-y0)^2)
    //  3 parameters defining the radial fall-off:
    //    r  radius for start of fall-off
    //    w  half-width of fall off
    //    A  amplitude of fall-off
    // r,w,A vary with azimuthal angle around tile centre
    // r,w stored as fractions of rad0 = tile half-width in pixels
    //
    //  then inverse scale g = (A/2) erfc [ (2/w)(d - r - w) ] + 1 - A
    //  d, r, w all in same units
    //
  public:
    CCDTile3(){}
    CCDTile3(const double& Xmax, const double& Ymax);
    void init(const double& Xmax, const double& Ymax);

    void clearCounts();

    //! symmetric = false to allow A to vary around the tile
    void setSymmetric(const bool& symmetric);

    int Nparams() const {return nparams;} //!< number of parameters
    //! number of smoothing parameters for each of r,w,A
    int NparamsSmooth() const {return nparams_smooth;}

    // Parameter order: r,w,A, each one with 5 parameters (ABCDE)
    // r,w in pixels
    // Store parameters
    void StoreParameters(const std::vector<double>& parameters);
    // Retrieve parameter vector (length nparams)
    std::vector<double> Parameters() const;

    //! return vector of internal ties, given SDs and 1st global parameter index
    std::vector<Tie> Ties(const std::vector<double> sdties,
			  const int& idx0);
    //! vector of indices and weights for each parameter to be restrained across tiles 
    std::pair<std::vector<int>, std::vector<double> > 
    TiedParameters(const std::vector<double> sdties,
		   const int& idx0);

    //! lower bounds for ipar'th parameter
    double LowerBound(const int& ipar) const;
    //! upper bounds for ipar'th parameter
    double UpperBound(const int& ipar) const;
    //! "large shift" for ipar'th parameter
    double LargeShift(const int& ipar) const;

    //! Return edge radius, pixels
    double Rad0() const {return rad0;}
    //! Return Xcentre, pixels
    double Xcentre() const {return x0*rad0;}
    //! Return Ycentre
    double Ycentre() const {return y0*rad0;}

    // Return scale & derivatives for tile coordinates Xt, Yt
    void ScaleDeriv(const bool& Deriv,
		    const double& Xt, const double& Yt,
		    double& scale, std::vector<double>& dgdp) const;

    // Format layout for printing
    std::string format() const;
    // Format type for printing
    std::string formattype() const;

    //! format parameters for printing
    std::vector<std::string> formatparameters(const std::vector<double>& sds) const;

    // Format all information into a labelled save format for later restoration
    std::string FormatSave() const;

    // restore
    void Restore(Fileread& FR);

  private:
    double xmax, ymax; // tile coordinates are 0->xmax, 0->ymax
    double rad0;     // radius to edge, pixels
    // Parameters
    bool circularlysymmetric; // false if A varies with polar coordinate
    int nparams;  // == 15
    int nparams_smooth;   // == 5, for each of r,w,A
    double xc0, yc0;   //  centre of coordinates, pixels
    double x0, y0;     // effective centre of taper, fraction of rad0
    
    double r0, w0, A0;  // constant part
    FourierSmooth rfs, wfs, Afs; // fall-off parameters, azimuth dependent

    double twooverrootpi;  // 2/sqrt(pi)

    // number of observations in each corner
    // corner is defined as
    // d > dcrnmin = sqrt(1/2((xmax/2)^2+(ymax/2)^2))
    mutable clipper::Array2d<int> ncorners;
    double dcrnmin;
    mutable RadialFunctionErfc radfunc;

    std::vector<std::string>
    format5(const bool& hasSd, const double& v0, const double& sd,
	    const FourierSmooth& vfs, const std::vector<double>& sdfs,
	    const int& width, const std::string& label) const;

  }; // end class CCDTile3
  //--------------------------------------------------------------
  class CCDTile1 : public TileBase {
    //! A CCD tile, correct for fall-off in the corners
    // Coordinates within the tile are defined from 0->xmax, 0->ymax

    // The model:
    //  Circularly symmetric around the point (x0,y0)
    //  distance from centre d = sqrt((X-x0)^2 + Y-y0)^2)
    //  3 parameters defining the radial fall-off:
    //    r  radius for start of fall-off
    //    w  half-width of fall off
    //    A  amplitude of fall-off
    //
    //  then inverse scale g = (A/2) erfc [ (2/w)(d - r - w) ] + 1 - A
    //

  public:
    CCDTile1(){}
    CCDTile1(const double& Xmax, const double& Ymax);
    void init(const double& Xmax, const double& Ymax);

    void clearCounts();

    int Nparams() const {return nparams;} //!< number of parameters

    // Parameter order: r,w,A,x0,y0
    // Store parameters
    void StoreParameters(const std::vector<double>& parameters);
    // Retrieve parameter vector (length nparams)
    std::vector<double> Parameters() const;

    //! return vector of internal ties, given SDs and 1st global parameter index
    std::vector<Tie> Ties(const std::vector<double> sdties,
			  const int& idx0);

    //! vector of indices and weights for each parameter to be restrained across tiles 
    std::pair<std::vector<int>, std::vector<double> > 
    TiedParameters(const std::vector<double> sdties,
		   const int& idx0);

    //! lower bounds for ipar'th parameter
    double LowerBound(const int& ipar) const;
    //! upper bounds for ipar'th parameter
    double UpperBound(const int& ipar) const;
    //! "large shift" for ipar'th parameter
    double LargeShift(const int& ipar) const;

    //! Return edge radius
    double Rad0() const {return rad0;}
    //! Return Xcentre
    double Xcentre() const {return x0*rad0;}
    //! Return Ycentre
    double Ycentre() const {return y0*rad0;}

    // Return scale & derivatives for tile coordinates Xt, Yt
    void ScaleDeriv(const bool& Deriv,
		    const double& Xt, const double& Yt,
		    double& scale, std::vector<double>& dgdp) const;

    // Format layout for printing
    std::string format() const;
    // Format type for printing
    std::string formattype() const;

    //! format parameters for printing
    std::vector<std::string> formatparameters(const std::vector<double>& sds) const;

    // Format all information into a labelled save format for later restoration
    std::string FormatSave() const;

    // restore
    void Restore(Fileread& FR);

  private:
    double xmax, ymax; // tile coordinates are 0->xmax, 0->ymax
    double rad0;     // radius to edge
    // Parameters
    int nparams;  // == 5
    double xc0, yc0;   //  centre of coordinates, pixels
    double x0, y0;     // effective centre of taper, fraction of rad0
    double r, w, A;    // fall-off parameters
    double twooverrootpi;  // 2/sqrt(pi)

    // number of observations in each corner
    // corner is defined as
    // d > dcrnmin = sqrt(1/2((xmax/2)^2+(ymax/2)^2))
    mutable clipper::Array2d<int> ncorners;
    double dcrnmin;
    mutable RadialFunctionErfc radfunc;
  }; // end class CCDTile1
  //--------------------------------------------------------------
  class CCDTile2 : public TileBase {
    //! A CCD tile, correct for fall-off in the corners
    // Coordinates within the tile are defined from 0->xmax, 0->ymax

    // The model:
    //  Circularly symmetric around the point (x0,y0) except for amplitude A
    //  distance from centre d = sqrt((X-x0)^2 + Y-y0)^2)
    //  3 parameters defining the radial fall-off:
    //    r  radius for start of fall-off
    //    w  half-width of fall off
    //    A  amplitude of fall-off
    //  A optionally varies with azimuthal angle around tile centre as
    //     A(constant) * (1 + Afs) where Afs is a 4-parameter Fourier series
    //
    //  then inverse scale g = (A/2) erfc [ (2/w)(d - r - w) ] + 1 - A
    //  d, r, w all in same units
    //
  public:
    CCDTile2(){}
    CCDTile2(const double& Xmax, const double& Ymax);
    void init(const double& Xmax, const double& Ymax);

    void clearCounts();

    //! symmetric = false to allow A to vary around the tile
    void setSymmetric(const bool& symmetric);

    int Nparams() const {return nparams;} //!< number of parameters
    //! number of smoothing parameters for each of r,w,A
    int NparamsSmooth() const {return nparams_smooth;}

    // Parameter order: r,w,A, each one with 5 parameters (ABCDE)
    // r,w in pixels
    // Store parameters
    void StoreParameters(const std::vector<double>& parameters);
    // Retrieve parameter vector (length nparams)
    std::vector<double> Parameters() const;

    //! return vector of internal ties, given SDs and 1st global parameter index
    std::vector<Tie> Ties(const std::vector<double> sdties,
			  const int& idx0);
    //! vector of indices and weights for each parameter to be restrained across tiles 
    std::pair<std::vector<int>, std::vector<double> > 
    TiedParameters(const std::vector<double> sdties,
		   const int& idx0);

    //! lower bounds for ipar'th parameter
    double LowerBound(const int& ipar) const;
    //! upper bounds for ipar'th parameter
    double UpperBound(const int& ipar) const;
    //! "large shift" for ipar'th parameter
    double LargeShift(const int& ipar) const;

    //! Return edge radius, pixels
    double Rad0() const {return rad0;}
    //! Return Xcentre, pixels
    double Xcentre() const {return x0*rad0;}
    //! Return Ycentre
    double Ycentre() const {return y0*rad0;}

    // Return scale & derivatives for tile coordinates Xt, Yt
    void ScaleDeriv(const bool& Deriv,
		    const double& Xt, const double& Yt,
		    double& scale, std::vector<double>& dgdp) const;

    // Format layout for printing
    std::string format() const;
    // Format type for printing
    std::string formattype() const;

    std::string formatTies() const;

    //! format parameters for printing
    std::vector<std::string> formatparameters(const std::vector<double>& sds) const;

    // Format all information into a labelled save format for later restoration
    std::string FormatSave() const;

    // restore
    void Restore(Fileread& FR);

  private:
    double xmax, ymax; // tile coordinates are 0->xmax, 0->ymax
    double rad0;     // radius to edge
    // Parameters
    bool circularlysymmetric; // false if A varies with polar coordinate
    int nparams;  // == 5 or 9
    int nparams_smooth;   // == 5, for each of A
    double xc0, yc0;   //  centre of coordinates, pixels
    double x0, y0;     // effective centre of taper, fraction of rad0
    double r;
    double w;
    double A0;

    std::vector<double> ties_;  // stored for printing

    FourierSmooth Afs; // fall-off parameters, azimuth dependent

    double twooverrootpi;  // 2/sqrt(pi)

    // number of observations in each corner
    // corner is defined as
    // d > dcrnmin = sqrt(1/2((xmax/2)^2+(ymax/2)^2))
    mutable clipper::Array2d<int> ncorners;
    double dcrnmin;
    mutable RadialFunctionErfc radfunc;
  }; // end class CCDTile2
   //--------------------------------------------------------------
  class FlatTile : public TileBase {
    //! flat tile, a single scale for each tile
  public:
    FlatTile():scale(1.0){}
    FlatTile(const double& Xmax, const double& Ymax):scale(1.0){nparams=1;}
    void init(const double& Xmax, const double& Ymax) {nparams=1;}

    void clearCounts(){}

    int Nparams() const {return nparams;} //!< number of parameters

    // Store parameters
    void StoreParameters(const std::vector<double>& parameters);
    // Retrieve parameter vector (length nparams)
    std::vector<double> Parameters() const;

    //! lower bounds for ipar'th parameter
    double LowerBound(const int& ipar) const;
    //! upper bounds for ipar'th parameter
    double UpperBound(const int& ipar) const;
    //! "large shift" for ipar'th parameter
    double LargeShift(const int& ipar) const;

    // Return scale & derivatives for tile coordinates Xt, Yt
    void ScaleDeriv(const bool& Deriv,
		    const double& Xt, const double& Yt,
		    double& scale, std::vector<double>& dgdp) const;

    // Format layout for printing
    std::string format() const;
    // Format type for printing
    std::string formattype() const;

    //! format parameters for printing
    std::vector<std::string> formatparameters(const std::vector<double>& sds) const;

    // Format all information into a labelled save format for later restoration
    std::string FormatSave() const;

    // restore
    void Restore(Fileread& FR);

  private:
    int nparams;
    double scale;
  };  // class FlatTile
  //--------------------------------------------------------------
  class TilePixel : public TileBase {
    //! A pixellated tile, one independent scale for each ngpxlX x ngpxlY pixel group
    // Coordinates within the tile are defined from 0->xmax, 0->ymax

    //  group coordinates jx, jy = X/ngpxlX, Y/ngpxlY
    //  then inverse scale g = gscale(jx,jy)
    //

  public:
    TilePixel(){}
    TilePixel(const double& Xmax, const double& Ymax);
    void init(const double& Xmax, const double& Ymax);

    void clearCounts(){}

    int Nparams() const {return nparams;} //!< number of parameters

    // Parameter order: jx = 0->njx-1, jy = 0->njy-1
    // Store parameters
    void StoreParameters(const std::vector<double>& parameters);
    // Retrieve parameter vector (length nparams)
    std::vector<double> Parameters() const;

    //! lower bounds for ipar'th parameter
    double LowerBound(const int& ipar) const;
    //! upper bounds for ipar'th parameter
    double UpperBound(const int& ipar) const;
    //! "large shift" for ipar'th parameter
    double LargeShift(const int& ipar) const;

    // Return scale & derivatives for tile coordinates Xt, Yt
    void ScaleDeriv(const bool& Deriv,
		    const double& Xt, const double& Yt,
		    double& scale, std::vector<double>& dgdp) const;

    // Format layout for printing
    std::string format() const;
    // Format type for printing
    std::string formattype() const;

    //! format parameters for printing
    std::vector<std::string> formatparameters(const std::vector<double>& sds) const;

    // Format all information into a labelled save format for later restoration
    std::string FormatSave() const;

    // restore
    void Restore(Fileread& FR);

  private:
    double xmax, ymax; // tile coordinates are 0->xmax, 0->ymax
    // Parameters
    int nparams;  // ==  njx * njy
    clipper::Array2d<double> scalexy;

    int ngpxlX, ngpxlY; // pixel binning
    int njx, njy;       // number of bins in each direction, = xmax/ngpxlX etc

  }; // end class TilePixel
}
#endif
