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

    DetectorType() :type(UNKNOWN), typestr("Unknown"),
		    ndet(0), ntilex(0), ntiley(0) {}

    //! construct from Batch object
    DetectorType(const Batch& batch);

    //! construct from arguments
    DetectorType(const Type& Dtype, const std::string& TypeLabel,
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

    //! return a string corresponding to the detector type
    std::string TypeLabel(const DetectorType::Type& dtype);

    //! test for equality on type and detector range (not on number of tiles)
    bool equals(const DetectorType& b) const;

    friend bool operator == (const DetectorType& a,const DetectorType& b);
    friend bool operator != (const DetectorType& a,const DetectorType& b);

  private:
    Type type;
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
    enum DetectorScaleType {NONE, FLAT, CCD, PIXEL, AUTOMATIC}; // scale type for each tile
    DetectorScale() : detectorscaletype(NONE), ntilex(0), ntiley(0) {}
    DetectorScale(const DetectorScaleType DetScaleType,
		  const int& nTileX, const int& nTileY,
		  const DetectorType& dettype);

    void init(const DetectorScaleType DetScaleType,
	      const int& nTileX, const int& nTileY,
	      const DetectorType& dettype);

    void init(const DetectorScaleType DetScaleType,
	      const int& nTileX, const int& nTileY,
	      const Range& Xrange, const Range& Yrange);

    void init();

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

    int Number() const {return nparams;}

    //! Return list of ties: sdties are sds for weight, idx0 is index to first global
    // parameter for setting ties, since they refer to the global parameter index
    //
    // For CCD tiles:
    //  sdties[0] for r
    //  sdties[1] for w
    //  sdties[2] for A
    //  sdties[3] for x0, y0
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

    std::string formatparameters() const;

    //! write output image of correction factors
    // as ADSC format image
    void WriteImage(const std::string& imagefilename) const;

    // Format all information into a labelled save format for later restoration
    std::string FormatSave() const;

    // restore
    void Restore(Fileread& FR);

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
  //--------------------------------------------------------------
  class TileBase {  // base class for tile scale
  public:
    TileBase(){}
    TileBase(const double& Xmax, const double& Ymax) {}

    virtual void init(const double& Xmax, const double& Ymax) = 0;

    // Pure virtual functions
    virtual int Nparams() const = 0; //!< number of parameters
    //! Store parameters
    virtual void StoreParameters(const std::vector<double>& parameters) = 0;
    //! Retrieve parameter vector (length nparams)
    virtual std::vector<double> Parameters() const = 0;

    //! lower bounds for ipar'th parameter
    virtual double LowerBound(const int& ipar) const = 0;
    //! upper bounds for ipar'th parameter
    virtual double UpperBound(const int& ipar) const = 0;
    //! "large shift" for ipar'th parameter
    virtual double LargeShift(const int& ipar) const = 0;

    //! Return maximum radius
    virtual double Radmax() const = 0;
    //! Return Xcentre
    virtual double Xcentre() const = 0;
    //! Return Ycentre
    virtual double Ycentre() const = 0;

    // Return scale & derivatives for tile coordinates Xt, Yt
    virtual void ScaleDeriv(const bool& Deriv,
			    const double& Xt, const double& Yt,
			    double& scale, std::vector<double>& dgdp) const = 0;

    // Format layout for printing
    virtual std::string format() const = 0;
    // Format formulation for printing
    virtual std::string formattype() const = 0;

    //! format parameters for printing
    virtual std::vector<std::string> formatparameters() const {}

    // Format all information into a labelled save format for later restoration
    std::string FormatSave() const;

    // restore
    void Restore(Fileread& FR);

  protected:
  };
  //--------------------------------------------------------------
  class CCDTile : public TileBase {
    //! A CCD tile, correct for fall-off in the corners
    // Coordinates within the tile are defined from 0->xmax, 0->ymax

    // The model:
    //  Radially symmetric around the point (x0,y0)
    //  distance from centre d = sqrt((X-x0)^2 + Y-y0)^2)
    //  3 parameters defining the radial fall-off:
    //    r  radius for start of fall-off
    //    w  half-width of fall off
    //    A  amplitude of fall-off
    //
    //  then inverse scale g = (A/2) erfc [ (2/w)(d - r - w) ] + 1 - A
    //

  public:
    CCDTile(){}
    CCDTile(const double& Xmax, const double& Ymax);
    void init(const double& Xmax, const double& Ymax);

    int Nparams() const {return nparams;} //!< number of parameters

    // Parameter order: r,w,A,x0,y0
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

    //! Return maximum radius
    double Radmax() const {return radmax;}
    //! Return Xcentre
    double Xcentre() const {return xc0;}
    //! Return Ycentre
    double Ycentre() const {return yc0;}

    // Return scale & derivatives for tile coordinates Xt, Yt
    void ScaleDeriv(const bool& Deriv,
		    const double& Xt, const double& Yt,
		    double& scale, std::vector<double>& dgdp) const;

    // Format layout for printing
    std::string format() const;
    // Format type for printing
    std::string formattype() const;

    //! format parameters for printing
    std::vector<std::string> formatparameters() const;

    // Format all information into a labelled save format for later restoration
    std::string FormatSave() const;

    // restore
    void Restore(Fileread& FR);

  private:
    double xmax, ymax; // tile coordinates are 0->xmax, 0->ymax
    double radmax;     // radius into farthest corner
    // Parameters
    int nparams;  // == 5
    double xc0, yc0;   //  centre of coordinates
    double x0, y0;     // effective centre of taper
    double r, w, A;    // fall-off parameters
    double twooverrootpi;  // 2/sqrt(pi)

    // number of obseravtions in each corner
    // corner is defined as
    // d > dcrnmin = sqrt(1/2((xmax/2)^2+(ymax/2)^2))
    mutable clipper::Array2d<int> ncorners;
    double dcrnmin;
  }; // end class CCDTile
  //--------------------------------------------------------------
  class FlatTile : public TileBase {
    //! flat tile, a single scale for each tile
  public:
    FlatTile():scale(1.0){}
    FlatTile(const double& Xmax, const double& Ymax):scale(1.0){nparams=1;}
    void init(const double& Xmax, const double& Ymax) {nparams=1;}

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

    //! Return maximum radius
    double Radmax() const {return 0.0;}
    //! Return Xcentre
    double Xcentre() const {return 0.0;}
    //! Return Ycentre
    double Ycentre() const {return 0.0;}

    // Return scale & derivatives for tile coordinates Xt, Yt
    void ScaleDeriv(const bool& Deriv,
		    const double& Xt, const double& Yt,
		    double& scale, std::vector<double>& dgdp) const;

    // Format layout for printing
    std::string format() const;
    // Format type for printing
    std::string formattype() const;

    //! format parameters for printing
    std::vector<std::string> formatparameters() const;

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

    // Dummies:-
    //! Return maximum radius
    double Radmax() const {return 0.0;}
    //! Return Xcentre
    double Xcentre() const {return 0.0;}
    //! Return Ycentre
    double Ycentre() const {return 0.0;}

    // Return scale & derivatives for tile coordinates Xt, Yt
    void ScaleDeriv(const bool& Deriv,
		    const double& Xt, const double& Yt,
		    double& scale, std::vector<double>& dgdp) const;

    // Format layout for printing
    std::string format() const;
    // Format type for printing
    std::string formattype() const;

    //! format parameters for printing
    std::vector<std::string> formatparameters() const;

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
