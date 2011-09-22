
#define ASSERT assert
#include <assert.h>

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;

#include "tile.hh"
#include "string_util.hh"
#include "imagearray.hh"
#include "hkl_unmerge.hh"

namespace scala {
  //--------------------------------------------------------------
  //! construct from Batch object
  DetectorType::DetectorType(const Batch& batch)
    :type(UNKNOWN), typestr("Unknown"),
     ndet(1), ntilex(0), ntiley(0)
  {
    if (batch.Ndet() > 1) {
      Message::message(Message_fatal
   ("DetectorType: cannot cope with more than one detector, update program"));
    }
    // At present, the only detector information we have is the
    // range of detector pixel coordinates
    detrange = batch.DetectorCoordinateRange();
    // Usually this contains pixel coordinates from 0
    // Check if this is so: if not, leave type as Unknown
    pixelcoords = true;
    if (Nint(detrange[0][0]) != 0 || Nint(detrange[1][0]) != 0) {
      pixelcoords = false;
    }
    if (((detrange[0][1] - detrange[0][0]) < 127.) ||
      ((detrange[1][1] - detrange[1][0]) < 127.)) {
      pixelcoords = false;
    }
    if (pixelcoords) {
      // Try to deduce detector type
    }
  }
  //--------------------------------------------------------------
  //! construct from arguments
  DetectorType::DetectorType(const Type& Dtype, const std::string& TypeLabel,
			     const std::vector<std::vector<float> >& Detrange)
    :type(Dtype), typestr(TypeLabel),
     ndet(1)
  {
    if (ndet > 1) {
      Message::message(Message_fatal
   ("DetectorType: cannot cope with more than one detector, update program"));
    }
    ntilex = ntiley = 1;
    if (type == CCD2x2) {ntilex = ntiley = 2;}
    if (type == CCD3x3) {ntilex = ntiley = 3;}

    // At present, the only detector information we have is the
    // range of detector pixel coordinates
    detrange = Detrange;
    // Usually this contains pixel coordinates from 0
    // Check if this is so: if not, leave type as Unknown
    pixelcoords = true;
    if (Nint(detrange[0][0]) != 0 || Nint(detrange[1][0]) != 0) {
      pixelcoords = false;
    }
    if (((detrange[0][1] - detrange[0][0]) < 127.) ||
      ((detrange[1][1] - detrange[1][0]) < 127.)) {
      pixelcoords = false;
    }
    if (pixelcoords) {
      // Try to deduce detector type
    }
  }
  //--------------------------------------------------------------
  //! valid: non-zero detector coordinate range
  bool DetectorType::Valid() const
  {
    bool OK = true;
    if (std::abs(detrange[1][0] - detrange[0][0]) < 0.001) {OK = false;}
    if (std::abs(detrange[1][1] - detrange[0][1]) < 0.001) {OK = false;}
    return OK; 
  }
  //--------------------------------------------------------------
  Range DetectorType::XdetRange() const
  //! Detector range on Xdet
  {
    return Range(detrange[0][0], detrange[0][1]);
  }
  //--------------------------------------------------------------
  Range DetectorType::YdetRange() const
  //! Detector range on Ydet
  {
    return Range(detrange[1][0], detrange[1][1]);
  }
  //--------------------------------------------------------------
  std::string DetectorType::TypeLabel(const DetectorType::Type& dtype)
  //! return a string corresponding to the detector type
  {
    std::string s;
    switch (dtype) {
    case UNKNOWN:
      s = "Unknown";
      break;
    case CCD1:  // one tile CCD
      s = "CCD 1-tile";
      break;
    case CCD2x2:
      s = "CCD 2x2-tile";
      break;
    case CCD3x3:
      s = "CCD 3x3-tile";
      break;
    case PILATUS6M:
      s = "Pilatus6M";
      break;
    case PILATUS2M:
      s = "Pilatus2M";
      break;
    }
    return s;
  }
  //--------------------------------------------------------------
  bool DetectorType::equals(const DetectorType& b) const
  //! test for equality on type and detector range (not on number of tiles)
  {
    if (type != b.type) return false;
    if (ndet != b.ndet) return false;
    for (int i=0;i<2;++i) {for (int j=0;j<2;++j) {
	if (std::abs(detrange[i][j]-b.detrange[i][j]) > 0.001) {
	  return false;
	}
      }}
    return true;
  }
  //--------------------------------------------------------------
  bool operator == (const DetectorType& a,const DetectorType& b)
  {
    return a.equals(b);
  }
  //--------------------------------------------------------------
  bool operator != (const DetectorType& a,const DetectorType& b)
  {
    return !(a.equals(b));
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  //! Construct to match reflection list
  DetectorAnalysis::DetectorAnalysis(const hkl_unmerge_list& hkl_list)
  {
    init(hkl_list);
  }
  //--------------------------------------------------------------
  //! Initialise to match reflection list
  void DetectorAnalysis::init(const hkl_unmerge_list& hkl_list)
  {
    int nrun = hkl_list.num_runs();
    std::vector<DetectorType> detectortypes(nrun);
    idxrun.assign(nrun, -1);
    int k = -1;
    detectorstatistics.clear();
    
    for (int irun=0;irun<nrun;irun++) {
      int b0 = hkl_list.RunList()[irun].BatchSerial0();  // 1st batch serial 
      Batch bat0 = hkl_list.Batches()[b0];    // first batch in run
      detectortypes[irun] = DetectorType(bat0);
      if (irun == 0) { // 1st irun
	k++;
	idxrun[irun] = k;  // index for run
	// Create new statistics object
	detectorstatistics.push_back(DetectorStatistics(detectortypes[irun]));
      } else { // run > 0
	// Check against earlier runs
	bool found = false;
	for (int jrun=0;jrun<irun;jrun++) {
	  if (detectortypes[jrun] == detectortypes[0]) {
	    // irun is same detector as jrun
	    idxrun[irun] = idxrun[jrun];
	    found = true;
	    break;
	  }
	}
	if (!found) {
	  // irun is new dataset
	  idxrun[irun] = ++k;
	  detectorstatistics.push_back
	    (DetectorStatistics(detectortypes[irun]));
	}
      }
    }  // end loop runs
    ndet = k+1;
    ASSERT (ndet == int(detectorstatistics.size()));
  }
  //--------------------------------------------------------------
  //! Construct explicitly for testing
  DetectorAnalysis::DetectorAnalysis(const std::vector<std::vector<float> >& Detrange)
  {
    int nrun = 1;
    DetectorType detectortypes;
    idxrun.assign(nrun, -1);
    detectorstatistics.clear();
    
    DetectorType::Type Dtype(DetectorType::UNKNOWN);
    std::string Typelabel("Unknown");

    detectortypes = DetectorType(Dtype, Typelabel, Detrange);
    idxrun[0] = 0;  // index for run
    // Create new statistics object
    detectorstatistics.push_back(DetectorStatistics(detectortypes));
    ndet = 1;
    ASSERT (ndet == int(detectorstatistics.size()));
  }
  //--------------------------------------------------------------
  void DetectorAnalysis::AddStats(const float& I, const float& AvI,
				  const int& runidx,
				  const int& xdet, const int& ydet)
  {
    detectorstatistics[idxrun[runidx]].AddStats(I, AvI, xdet, ydet);
  }
  //--------------------------------------------------------------
  void DetectorAnalysis::WriteImages(const std::string& fname) const
  {
    std::string imagefilename = fname;
    if (imagefilename == "") imagefilename = "DETECTORIMAGE";
    if (getenv(imagefilename.c_str()) != NULL) { // it's an environment variable
      imagefilename = std::string(getenv(imagefilename.c_str()));
    }

    for (int idsc=0;idsc<ndet;++idsc) {
      std::string basename = FileNameNoExtension(imagefilename);
      std::string ext = FileNameExtension(imagefilename);
      if (ext == "") {ext = "img";}
      std::string trailer = "_"+StringUtil::Strip(StringUtil::itos(idsc+1,4))+"."+ext;
      std::string name = basename+"_scales"+trailer;
      detectorstatistics[idsc].WriteImageScales(name);
      name = basename+"_deviations"+trailer;
      detectorstatistics[idsc].WriteImageDeviations(name);
    }
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  DetectorStatistics::DetectorStatistics(const DetectorType& Dtype)
  //! construct from detector type object
  {
    init(Dtype);
  }
  //--------------------------------------------------------------
  void DetectorStatistics::init(const DetectorType& Dtype)
  //! initialise from detector type object
  {
    detectortype = Dtype;
    // Pixel binning
    const int kbin = 32;
    SetBinning(kbin, kbin);
  }    
  //--------------------------------------------------------------
  void DetectorStatistics::SetBinning(const int& Ngpxlx, const int& Ngpxly)
  {
    ngpxlX = Ngpxlx;
    ngpxlY = Ngpxly;
    xdetrange = detectortype.XdetRange();
    ydetrange = detectortype.YdetRange();
    nbx = Nint(xdetrange.AbsRange()/double(ngpxlX));
    nby = Nint(ydetrange.AbsRange()/double(ngpxlY));
    xdetrange.SetNbin(nbx);
    ydetrange.SetNbin(nby);
    sumwIxy.resize(nbx, nby, 0.0);
    sumwIothers.resize(nbx, nby, 0.0);
    avdelta.resize(nbx, nby, MeanSD());
  }
  //--------------------------------------------------------------
  void DetectorStatistics::AddStats(const float& I, const float& AvI,
				    const int& xdet, const int& ydet)
  {
    int jx = xdetrange.tbin(xdet);
    int jy = ydetrange.tbin(ydet);
    ASSERT ((jx >= 0) && (jy >= 0));
    //^
    if (jx >= sumwIxy.rows() || jy >= sumwIxy.cols() || 
	jx >= sumwIothers.rows() || jy >= sumwIothers.cols() ||
	jx >= avdelta.rows() || jy >= avdelta.cols()) {
      std::cout << "Fail\n";
    }
    float w = 1.0;  // unit weights
    sumwIxy(jx,jy) += w * I;
    sumwIothers(jx,jy) += w * AvI;
    float d = (I-AvI);
    avdelta(jx,jy).Add(d);    // deviations
  }
  //--------------------------------------------------------------
  // Scale factors for each pixel group
  void DetectorStatistics::WriteImageScales(const std::string& fname) const
  {
    //^
    //    std::cout <<"Sizes: "
    //	      <<sumwIxy.rows()<<" "<< sumwIxy.cols() <<" "
    //	      <<sumwIothers.rows()<<" "<< sumwIothers.cols() <<" "
    //	      <<avdelta.rows()<<" "<< avdelta.cols() <<"\n";
    //^-
    // Calculate scales = Sum w I(xy) / Sum w Iothers(xy)
    clipper::Array2d<double> scale(sumwIxy.rows(), sumwIxy.cols());
    for (int ix=0;ix<sumwIxy.rows();++ix) {
      for (int iy=0;iy<sumwIxy.cols();++iy) {
	if (sumwIothers(ix,iy) != 0.0) {
	  scale(ix,iy) /= sumwIothers(ix,iy);
	} else {
	  scale(ix,iy) = 1.0;
	}
	if (scale(ix,iy) < 0.0) {
	  scale(ix,iy) = 0.0;
	}
      }}
    Imagearray imagearray;
    imagearray.SetScale(1000.);
    imagearray.init(scale);
    imagearray.Write(fname);
  }
  //--------------------------------------------------------------
  // Deviations for each pixel group
  void DetectorStatistics::WriteImageDeviations(const std::string& fname) const
  {
    // write out signed deviations (I - <I>)
    Imagearray imagearray;
    imagearray.SetScale(1000.);
    imagearray.init(avdelta, true);
    imagearray.Write(fname);
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  DetectorScale::DetectorScale(const DetectorScaleType DetScaleType,
			       const int& nTileX, const int& nTileY,
			       const DetectorType& dettype)
  {
    init(DetScaleType, nTileX, nTileY, dettype);
  }
  //--------------------------------------------------------------
  void DetectorScale::init(const DetectorScaleType DetScaleType,
			   const int& nTileX, const int& nTileY,
			   const DetectorType& dettype)
  {
    detectorscaletype = DetScaleType;
    ntilex = nTileX;
    ntiley = nTileY;
    if (ntilex <= 0) {ntilex = dettype.NtileX();}
    if (ntiley <= 0) {ntiley = dettype.NtileY();}

    tilescales.resize(ntilex, ntiley);
    idx_tile.resize(ntilex, ntiley);
    // Store ranges and number of tiles
    xdrange = dettype.XdetRange();
    ydrange = dettype.YdetRange();
    init();
  }
  //--------------------------------------------------------------
  void DetectorScale::init(const DetectorScaleType DetScaleType,
			   const int& nTileX, const int& nTileY,
			   const Range& Xrange, const Range& Yrange)
  {
    detectorscaletype = DetScaleType;
    ntilex = nTileX;
    ntiley = nTileY;
    tilescales.resize(ntilex, ntiley);
    idx_tile.resize(ntilex, ntiley);

    // Store ranges and number of tiles
    xdrange = Xrange;
    ydrange = Yrange;
    init();
  }
  //--------------------------------------------------------------
  void DetectorScale::init()
  {
    xdrange.SetNbin(ntilex);
    double sizex = xdrange.AbsRange()/double(ntilex);
    ydrange.SetNbin(ntiley);
    double sizey = ydrange.AbsRange()/double(ntiley);
    nparams = 0;
    for (int i=0;i<ntilex;++i) { // loop x
      for (int j=0;j<ntiley;++j) { // loop y
	if (detectorscaletype == FLAT) {
	  tilescales(i,j) = new FlatTile;
	} else if (detectorscaletype == CCD) {
	  tilescales(i,j) = new CCDTile;
	} else if (detectorscaletype == PIXEL) {
	  tilescales(i,j) = new TilePixel;
	}
	tilescales(i,j)->init(sizex, sizey);
	idx_tile(i,j) = nparams;
	nparams += tilescales(i,j)->Nparams();
      }}
  }
  //--------------------------------------------------------------
  DetectorScale::~DetectorScale() {
    if (ntilex > 0 && ntiley > 0) {
      for (int i=0;i<ntilex;++i) { // loop x
	for (int j=0;j<ntiley;++j) { // loop y
	  delete tilescales(i,j);
	}}
    }
  }
  //--------------------------------------------------------------
  // Copy & copy constructor should fail or be done properly due to pointers
  DetectorScale::DetectorScale(const DetectorScale& detscale) 
    : detectorscaletype(NONE), ntilex(0), ntiley(0)
  {
    if (detscale.detectorscaletype != NONE) {
      Message::message(Message_fatal
		       ("DetectorScale: illegal copy constructor"));
    }
  }
  //--------------------------------------------------------------
  DetectorScale& DetectorScale::operator= (const DetectorScale& detscale) {
    if (detscale.detectorscaletype != NONE) {
      Message::message(Message_fatal
		       ("DetectorScale: illegal copy operator"));
    }
    init(NONE, 0, 0, Range(), Range());
    return *this;
  }
  //--------------------------------------------------------------
  //! Return true if Valid
  bool DetectorScale::Valid() const
  {
    if (detectorscaletype == NONE) {
      return true;
    } else if (detectorscaletype == FLAT) {
      // Only sensible if more than one tile
      if (ntilex <= 0 || ntiley <= 0) {return false;}
    } else if (detectorscaletype == CCD) {
      // At least one tile
      if (ntilex <= 1 || ntiley <= 1) {return false;}
      if (nparams <= 0) {return false;}
    }
    return true;
  }
  //--------------------------------------------------------------
  // Store parameter vector (length nparams)
  void DetectorScale::StoreParameters(const std::vector<double>& parameters)
  {
    ASSERT (int(parameters.size()) == nparams);
    std::vector<double>::const_iterator pos1 = parameters.begin();  // start of range
    std::vector<double>::const_iterator pos2;                   // end of range
    
    for (int i=0;i<ntilex;++i) { // loop x
      for (int j=0;j<ntiley;++j) { // loop y
	pos2 = pos1 + tilescales(i,j)->Nparams();
	tilescales(i,j)->StoreParameters(std::vector<double>(pos1, pos2));
	//^
	//	std::cout << "Tile "<<i<<" "<<j
	//		  <<" "<<tilescales(i,j)->formatparameters()<<"\n";
	//^-
	pos1 = pos2;
      }}
  }
  //--------------------------------------------------------------
  // Retrieve parameter vector (length nparams)
  std::vector<double> DetectorScale::Parameters() const
  {
    std::vector<double> params;
    for (int i=0;i<ntilex;++i) { // loop x
      for (int j=0;j<ntiley;++j) { // loop y
	std::vector<double> pars = tilescales(i,j)->Parameters();
	params.insert(params.end(), pars.begin(), pars.end());
      }}
    return params;
  }
  //--------------------------------------------------------------
  void DetectorScale::StoreNobservations(const std::vector<int>& Nobs)
  // Store Nobs vector (length nparams)
  {
    ASSERT (int(Nobs.size()) == nparams);
    nobsPar = Nobs;
  }
  //--------------------------------------------------------------
  //! find which tile ipar'th parameter belongs to
  std::pair<int,int> DetectorScale::WhichTile(const int& ipar) const
  {
    int i,j;
    bool found = false;
    for (i=0;i<ntilex;++i) { // loop x
      for (j=0;j<ntiley;++j) { // loop y
	if (ipar < idx_tile(i,j)) {
	  found = true; j--; break;}
      }
      if (found) {break;}
    }
    if (!found) {
      i=ntilex-1; j=ntiley-1;
    }
    return  std::pair<int,int>(i,j);
  }
  //--------------------------------------------------------------
  //! lower bounds for ipar'th parameter
  bool DetectorScale::LowerBound(const int& ipar, double& Lower) const
  {
    std::pair<int,int> idxtile = WhichTile(ipar);
    Lower = tilescales(idxtile.first, idxtile.second)->
      LowerBound(ipar-idx_tile(idxtile.first, idxtile.second));
    return true;
  }
  //--------------------------------------------------------------
  //! upper bounds for ipar'th parameter
  bool DetectorScale::UpperBound(const int& ipar, double& Upper) const
  {
    std::pair<int,int> idxtile = WhichTile(ipar);
    Upper = tilescales(idxtile.first, idxtile.second)->
      UpperBound(ipar-idx_tile(idxtile.first, idxtile.second));
    return true;
  }
  //--------------------------------------------------------------
  //! "large shift" for ipar'th parameter
  double DetectorScale::LargeShift(const int& ipar) const
  {
    std::pair<int,int> idxtile = WhichTile(ipar);
    return tilescales(idxtile.first, idxtile.second)->
      LargeShift(ipar-idx_tile(idxtile.first, idxtile.second));
  }
  //--------------------------------------------------------------
  std::vector<Tie> DetectorScale::Ties(const std::vector<double> sdties,
				       const int& idx0) const
  // Return list of ties: sdties are sds for weight, idx0 is index to first global
  // parameter for setting ties, since they refer to the global parameter index
  //
  // For CCD tiles:
  //  sdties[0] for r
  //  sdties[1] for w
  //  sdties[2] for A
  //  sdties[3] for x0, y0
  {
    std::vector<Tie> ties;
    if (detectorscaletype != CCD) {return ties;}
    // Only ties for CCD tiled detector (at present)

    ASSERT (sdties.size() == 4);
    // x0, y0 tie to centre position
    // sdties[3] is relative to tile size
    int idx = idx0 + 3;  // first x0 parameter
    for (int i=0;i<ntilex;++i) { // loop tile x
      for (int j=0;j<ntiley;++j) { // loop tile y
	double weight = sdties[3] * tilescales(i,j)->Radmax();
	weight = 1./(weight*weight);
	// tie x0 to centre
	ties.push_back(Tie(idx, tilescales(i,j)->Xcentre(), weight));
	// tie y0 to centre
	ties.push_back(Tie(idx+1, tilescales(i,j)->Ycentre(), weight));
	idx += tilescales(i,j)->Nparams(); // increment for next tile tie
      }}

    // Tie r,w,A parameters together for all tiles
    for (int k=0;k<3;++k) { // loop parameters 0,1,2 = r,w,A
      // Use Radmax from central tile to scale weight for r & w
      double weight = sdties[k];
      if (k<2) {
	weight *= tilescales(ntilex/2,ntiley/2)->Radmax();
      }
      weight = 1./(weight*weight);

      std::vector<int> kindex;  // parameter indices for this group of parameters
      idx = idx0 + k;   // starting index
      for (int i=0;i<ntilex;++i) { // loop tile x
	for (int j=0;j<ntiley;++j) { // loop tile y
	  kindex.push_back(idx);
	  idx += tilescales(i,j)->Nparams(); // increment for next tile tie
	}}
      ties.push_back(Tie(kindex, weight));  // add tie to list
    } // end loop r,w,A

    return ties;
  }
  //--------------------------------------------------------------
  // Return scale for detector coordinates Xdet, Ydet
  double DetectorScale::Scale(const std::pair<float,float>& XYdet) const
  {
    double scale;
    std::vector<double> dgdp;
    ScaleDeriv(false, XYdet, scale, dgdp);
    return scale;
  }
  //--------------------------------------------------------------
  // Return scale & derivatives for detector coordinates Xdet, Ydet
  double DetectorScale::ScaleDeriv(const std::pair<float,float>& XYdet,
				   std::vector<double>& dgdp) const
  {
    double scale;
    ScaleDeriv(true, XYdet, scale, dgdp);
    return scale;
  }
  //--------------------------------------------------------------
  // Return scale & derivatives for  detector coordinates Xdet, Ydet
  void DetectorScale::ScaleDeriv(const std::pair<float,float>& XYdet,
				 double& scale, std::vector<double>& dgdp) const
  {
    ScaleDeriv(true, XYdet, scale, dgdp);
  }
  //--------------------------------------------------------------
    // Return scale & derivatives for  detector coordinates Xdet, Ydet
  void DetectorScale::ScaleDeriv(const bool& Deriv,
				 const std::pair<float,float>& XYdet,
				 double& scale, std::vector<double>& dgdp) const
  {
    // Which tile?
    int xtile = xdrange.bin(XYdet.first);
    int ytile = ydrange.bin(XYdet.second);
    // Tile coordinates
    double Xt = XYdet.first - xdrange.bounds(xtile).first;
    double Yt = XYdet.second - ydrange.bounds(ytile).first;
    //^
    //    std::cout << "xtl,ytl,Xt,Yt " << xtile <<" "<<ytile<<" "<<Xt<<" "<<Yt<<"\n";

    // get scales & optional derivatives
    std::vector<double> dgdt;
    tilescales(xtile,ytile)->ScaleDeriv(Deriv, Xt, Yt, scale, dgdt);
    dgdp.assign(nparams,0.0);
    std::copy(dgdt.begin(), dgdt.end(), dgdp.begin()+idx_tile(xtile,ytile));
    //^
    //    std::cout << "DetectorScale::ScaleDeriv "
    //	      << Xt <<" "<<Yt<<"\n";
    //    for (size_t i=0;i<dgdp.size();++i) {std::cout <<" "<<dgdp[i];}
    //    std::cout <<"\n";
    //^-
  }
  //--------------------------------------------------------------
  std::string DetectorScale::format() const
  //! Format for printing scale type
  {
    std::string s;
    if (detectorscaletype == NONE) {
      s = "No detector scale";
    } else {
      // mid-tile
      s = tilescales(ntilex/2, ntiley/2)->format();
      s += ", number of tiles: "+StringUtil::StringUtil::itos(ntilex,2)+
	" x "+StringUtil::StringUtil::itos(ntiley,2);
    }
    return s;
  }
  //--------------------------------------------------------------
  std::string DetectorScale::formatparameters() const
  //! format parameters for printing
  {
    std::string s = "\n";
    if (detectorscaletype == NONE) {return s;}

    s += "Detector scales for"+
    StringUtil::itos(ntilex,2)+" x"+StringUtil::itos(ntiley,2)+
      " tiles\n";
    s += tilescales(ntilex/2, ntiley/2)->formattype()+"\n";

    s += "    Tile parameters arranged with Xdet across and Ydet up\n";
    s += "    Numbers in the corners are the number of observations contributing to the scales\n";

    // tile scale
    for (int j=ntiley-1;j>=0;j--) { // loop y backwards
      // Get formatted squares for x tiles
      // these should all have the same number of lines
      std::vector<std::vector<std::string> > sxtiles(ntilex);
      for (int i=0;i<ntilex;++i) { // loop x
	  sxtiles[i] = tilescales(i,j)->formatparameters();
	  if (i>0) {
	    ASSERT (sxtiles[i].size() == sxtiles[0].size());
	  }
	}
      // Print them out across the page (to string)
      for (size_t k=0;k<sxtiles[0].size();++k) { // loop lines
	for (int i=0;i<ntilex;++i) { // loop x tiles
	  s += sxtiles[i][k];
	}
	s += "\n";
      }
    } // end loop y tiles
    return s;
  }
  //--------------------------------------------------------------
  void DetectorScale::WriteImage(const std::string& imagefilename) const
  //! write output image of correction factors
  // as ADSC format image
  {
    int nxp = Nint(xdrange.AbsRange());
    int nyp = Nint(xdrange.AbsRange());

    clipper::Array2d<double> image(nxp, nyp);

    for (int y=0;y<nyp;++y) {
      float Ydet = y;
      for (int x=0;x<nxp;++x) {
	float Xdet = x;
	image(x,y) = Scale(std::pair<float,float>(Xdet,Ydet));
      }}
    Imagearray imagearray;
    imagearray.SetScale(1000.);
    imagearray.init(image);
    imagearray.Write(imagefilename);
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  CCDTile::CCDTile(const double& Xmax, const double& Ymax)
  {
    init(Xmax, Ymax);
  }
  //--------------------------------------------------------------
  void CCDTile::init(const double& Xmax, const double& Ymax)
  // coordinates in range 0->Xmax, 0->Ymax
  {
    xmax = Xmax;
    ymax = Ymax;
    nparams = 5;
    twooverrootpi = 2.0/sqrt(clipper::Util::pi()); // 2/sqrt(pi)
    // Centre of coordinate system (for ties)
    xc0 = xmax/2.0;
    yc0 = ymax/2.0;
    // Set default (initial) values
    x0 = xc0;
    y0 = yc0;
    double xymax = Max(x0, y0);  // larger radius
    radmax = sqrt(2.*xymax*xymax); // to corner
    r = radmax * 0.7;
    w = radmax * 0.1;
    A = 0.1;
    ncorners.resize(2,2,0); // counts in corners, initialise to 0
    // dcrnmin = sqrt(1/2((xmax/2)^2+(ymax/2)^2)) limit for corner
    dcrnmin = sqrt(0.5*0.25*(xmax*xmax + ymax*ymax));
  }
  //--------------------------------------------------------------
  // Parameter order: r,w,A,x0,y0,
  // Store parameters
  void CCDTile::StoreParameters(const std::vector<double>& parameters)
  {
    ASSERT (int(parameters.size()) == nparams);
    r = parameters[0];
    w = parameters[1];
    A = parameters[2];
    x0 = parameters[3];
    y0 = parameters[4];
    //^
    std::cout << "r,w,A,x0,y0 ";
    for (int i=0;i<5;++i) {std::cout <<" "<<parameters[i];}
    std::cout <<"\n";  //^-
  }
  //--------------------------------------------------------------
  // Retrieve parameter vector (length nparams)
  std::vector<double> CCDTile::Parameters() const
  {
    std::vector<double> par(nparams);
    par[0] = r;
    par[1] = w;
    par[2] = A;
    par[3] = x0;
    par[4] = y0;
    return par;
  }
  //--------------------------------------------------------------
  //! lower bounds for ipar'th parameter
  double CCDTile::LowerBound(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    if (ipar <= 1) {
      return 0.1 * radmax; // r or w
    } else if (ipar == 2 ) {
      return 0.0; // A
    }
    return -0.5 * radmax;  // x0 or y0
  }
  //--------------------------------------------------------------
  //! upper bounds for ipar'th parameter
  double CCDTile::UpperBound(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    if (ipar <= 1) {
      return radmax; // r or w
    } else if (ipar == 2 ) {
      return 0.5; // A
    }
    return +0.5 * radmax;  // x0 or y0
  }
  //--------------------------------------------------------------
  //! "large shift" for ipar'th parameter
  double CCDTile::LargeShift(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    if (ipar <= 1) {
      return 0.3 * radmax; // r or w
    } else if (ipar == 2 ) {
      return 0.2; // A
    }
    return 0.4 * radmax;  // x0 or y0
  }
  //--------------------------------------------------------------
  // Return scale & derivatives for tile coordinates Xt, Yt
  void CCDTile::ScaleDeriv(const bool& Deriv,
				 const double& Xt, const double& Yt,
				 double& scale,
				 std::vector<double>& dgdp) const
  {
    double d = sqrt((Xt-x0)*(Xt-x0) + (Yt-y0)*(Yt-y0));
    double z = 2.0*(d-r-w)/w;
    double erfz = erfc(z);
    scale = 0.5 * A * erfz + 1.0 - A;
    // Count observations in each corner
    if (d > dcrnmin) { // in a corner
      int i = (Xt-xc0) < 0.0 ? 0 : 1; // 0 or 1 if left or right
      int j = (Yt-yc0) < 0.0 ? 0 : 1; // 0 or 1 if bottom or top
      ncorners(i,j)++;
    }

    //^
    //    std::cout <<"\nXt,x0,Yt,y0 "<<Xt<<" "<<x0<<" "<<Yt<<" "<<y0<<"\n";
    //    std::cout <<"d,z,erfz "<<d<<" "<<z<<" "<<erfz<<"\n"; //^-

    if (Deriv) {
      dgdp.resize(nparams);
      dgdp[0] = (twooverrootpi * A / w) * exp(-z*z);  // dgdr = -dgdd
      dgdp[1] = dgdp[0] * (1. + 0.5*z);               // dgdw 
      dgdp[2] = 0.5 * erfz - 1.0;                     // dgdA
      if (d == 0.0) {
	dgdp[3] = 0.0;
	dgdp[4] = 0.0;
      } else {
	dgdp[3] = - dgdp[0] * (x0 - Xt)/d;  // dgdx0 = dg/dd * dd/dx0
	dgdp[4] = - dgdp[0] * (y0 - Yt)/d;  // dgdy0 = dg/dd * dd/dy0
      }
    }
  }
  //--------------------------------------------------------------
  std::string CCDTile::format() const
  //! format scale type for printing
  {
    std::string text = "Tile correction for CCD detector";
    return text;
  }
  //--------------------------------------------------------------
  std::string CCDTile::formattype() const
  //! format scale type for printing
  {
    std::string text = format()+"\n\n";
    text += "   Scale up the corners of the tiles to allow for fall-off in the taper\n";
    text += "   Inverse scale g calculated from coordinate within the tile (xd,yd) as\n";
    text += "     g = (A/2) erfc(z) + 1 - A\n";
    text += "      where A is amplitude and z = 2(d - r - w)/w\n";
    text += "       d is the distance of the pixel (xd, yd) from the effective tile centre (x0, y0)\n";
    text += "       r is the radius of the fall-off\n";
    text += "       w is the width (steepness) of the fall-off\n";
    return text;
  }
  //--------------------------------------------------------------
  std::vector<std::string> CCDTile::formatparameters() const
  //! format parameters for printing
  {
    std::vector<std::string> s;
    std::string line;
    int width = 32; // width of window
    line = "";
    for (int i=0;i<width;++i) {line += "-";};
    s.push_back(line);
    // counts for top corners
    line = "| "+
      StringUtil::LeftString(StringUtil::Strip(StringUtil::itos(ncorners(0,1),6)),6)+
      StringUtil::PadString(" ",width-12-4)+
      StringUtil::RightString(StringUtil::Strip(StringUtil::itos(ncorners(1,1),6)),6)+
      +" |";
    s.push_back(line);
    line = "| "+ StringUtil::PadString(" ",width-4)+ " |"; // "blank" line
    s.push_back(line);
    // r, w, A
    line = "| "+
      StringUtil::CentreString(("r="+StringUtil::ftos(r/radmax,6,2)+
				", w="+StringUtil::ftos(w/radmax,6,2)+
				", A="+StringUtil::ftos(A,6,2)), width-4)+
      " |";
    s.push_back(line);
    line = "| "+
      StringUtil::CentreString
      (("Tile centre: "+StringUtil::Strip(StringUtil::ftos(x0,7,1))+
	", "+StringUtil::Strip(StringUtil::ftos(y0,7,1))), width-4)+
      " |";
    s.push_back(line);
    line = "| "+ StringUtil::PadString(" ",width-4)+ " |"; // "blank" line
    s.push_back(line);
    // counts for bottom corners
    line = "| "+
      StringUtil::LeftString(StringUtil::Strip(StringUtil::itos(ncorners(0,0),6)),6)+
      StringUtil::PadString(" ",width-12-4)+
      StringUtil::RightString(StringUtil::Strip(StringUtil::itos(ncorners(1,0),6)),6)+
      +" |";
    s.push_back(line);
    line = "";
    for (int i=0;i<width;++i) {line += "-";};
    s.push_back(line);
    return s;
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  // Store parameters (just one)
  void FlatTile::StoreParameters(const std::vector<double>& parameters)
  {
    ASSERT (int(parameters.size()) == nparams);
    scale = parameters[0];
  }
  //--------------------------------------------------------------
  // Retrieve parameter vector (length nparams)
  std::vector<double> FlatTile::Parameters() const
  {
    std::vector<double> par(nparams);
    par[0] = scale;
    return par;
  }
  //--------------------------------------------------------------
  //! lower bounds for ipar'th parameter
  double FlatTile::LowerBound(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    return 0.01;
  }
  //--------------------------------------------------------------
  //! upper bounds for ipar'th parameter
  double FlatTile::UpperBound(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    return 10.0;
  }
  //--------------------------------------------------------------
  //! "large shift" for ipar'th parameter
  double FlatTile::LargeShift(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    return 0.1;
  }
  //--------------------------------------------------------------
  // Return scale & derivatives for tile coordinates Xt, Yt
  void FlatTile::ScaleDeriv(const bool& Deriv,
			    const double& Xt, const double& Yt,
			    double& scale,
			    std::vector<double>& dgdp) const
  {
    dgdp.assign(1,1.0);
    scale = 1.0;
  }
  //--------------------------------------------------------------
  std::string FlatTile::format() const
  //! format scale type for printing
  {
    std::string text = "Flat constant tile correction";
    return text;
  }
  //--------------------------------------------------------------
  std::string FlatTile::formattype() const
  //! format scale type for printing
  {
    std::string text = format()+"\n";
    text += "  One single scale/tile\n";
    return text;
  }
  //--------------------------------------------------------------
  std::vector<std::string> FlatTile::formatparameters() const
  //! format parameters for printing
  {
    std::vector<std::string> s(1,"  Tile scale: "+StringUtil::ftos(scale,8,3));
    return s;
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  TilePixel::TilePixel(const double& Xmax, const double& Ymax)
  {
    init(Xmax, Ymax);
  }
  //--------------------------------------------------------------
  void TilePixel::init(const double& Xmax, const double& Ymax)
  // coordinates in range 0->Xmax, 0->Ymax
  {
    xmax = Xmax;
    ymax = Ymax;
    ngpxlX = ngpxlY = 64; // pixel binning
    njx = xmax/ngpxlX;
    njy = ymax/ngpxlY;

    nparams = njx * njy;

    scalexy.resize(njx, njy, 1.0);  // set all scales to 1.0
  }
  //--------------------------------------------------------------
  // Parameter order: r,w,A,x0,y0,
  // Store parameters
  void TilePixel::StoreParameters(const std::vector<double>& parameters)
  {
    ASSERT (int(parameters.size()) == nparams);
    int k = -1;
    for (int jx=0;jx<njx;++jx) {
      for (int jy=0;jy<njy;++jy) {
	scalexy(jx,jy) = parameters[++k];
      }}
  }
  //--------------------------------------------------------------
  // Retrieve parameter vector (length nparams)
  std::vector<double> TilePixel::Parameters() const
  {
    std::vector<double> par(nparams);
    int k = -1;
    for (int jx=0;jx<njx;++jx) {
      for (int jy=0;jy<njy;++jy) {
	par[++k] = scalexy(jx,jy);
      }}
    return par;
  }
  //--------------------------------------------------------------
  //! lower bounds for ipar'th parameter
  double TilePixel::LowerBound(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    return 0.0;  // x0 or y0
  }
  //--------------------------------------------------------------
  //! upper bounds for ipar'th parameter
  double TilePixel::UpperBound(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    return 5.0;  // x0 or y0
  }
  //--------------------------------------------------------------
  //! "large shift" for ipar'th parameter
  double TilePixel::LargeShift(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    return 0.2;  // x0 or y0
  }
  //--------------------------------------------------------------
  // Return scale & derivatives for tile coordinates Xt, Yt
  void TilePixel::ScaleDeriv(const bool& Deriv,
				 const double& Xt, const double& Yt,
				 double& scale,
				 std::vector<double>& dgdp) const
  {
    int jx = Xt/ngpxlX;
    int jy = Yt/ngpxlY;
    scale = scalexy(jx,jy);

    if (Deriv) {
      dgdp.assign(nparams,0.0);
      dgdp[jx*njx+jy] = 1.0;
    }
  }
  //--------------------------------------------------------------
  std::string TilePixel::format() const
  //! format scale type for printing
  {
    std::string text = "Pixel by pixel tile correction";
    return text;
  }
  //--------------------------------------------------------------
  std::string TilePixel::formattype() const
  //! format scale type for printing
  {
    std::string text = format()+"\n\n";
    text += "   Determine independent scale for each pixel\n";
    return text;
  }
  //--------------------------------------------------------------
  std::vector<std::string> TilePixel::formatparameters() const
  //! format parameters for printing
  {
    std::vector<std::string> s;
    return s;
  }
  //--------------------------------------------------------------
}
