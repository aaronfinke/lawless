
#define ASSERT assert
#include <assert.h>

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;

#include "tile.hh"
#include "string_util.hh"
#include "imagearray.hh"
#include "hkl_unmerge.hh"
#include "jiffy.hh"

using phaser_io::itos;
using phaser_io::dtos;
using phaser_io::ftos;

namespace scala {
  //--------------------------------------------------------------
  //! construct from Batch object
  DetectorType::DetectorType(const Batch& batch)
    :dettype(UNKNOWN), typestr("Unknown"),
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
      int xdetrange = Nint(detrange[0][1] - detrange[0][0]);
      int ydetrange = Nint(detrange[1][1] - detrange[1][0]);
      if (xdetrange == 3072 && ydetrange == 3072) {
        dettype = CCD3x3;
        ntilex = 3;
        ntiley = 3;
      }
    }
  }
  //--------------------------------------------------------------
  //! construct from arguments
  DetectorType::DetectorType(const Type& Dtype, const std::string& TypeLabel,
                             const std::vector<std::vector<float> >& Detrange)
    :dettype(Dtype), typestr(TypeLabel),
     ndet(1)
  {
    if (ndet > 1) {
      Message::message(Message_fatal
   ("DetectorType: cannot cope with more than one detector, update program"));
    }
    ntilex = ntiley = 1;
    if (dettype == CCD2x2) {ntilex = ntiley = 2;}
    if (dettype == CCD3x3) {ntilex = ntiley = 3;}

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
  //! construct from arguments
  DetectorType::DetectorType(const std::string& TypeLabel,
                             const std::vector<std::vector<float> >& Detrange)
    : typestr(TypeLabel), ndet(1)
  {
    if (ndet > 1) {
      Message::message(Message_fatal
   ("DetectorType: cannot cope with more than one detector, update program"));
    }
    //
    dettype = TypeFromLabel(typestr);
    ntilex = ntiley = 1;
    if (dettype == CCD2x2) {ntilex = ntiley = 2;}
    if (dettype == CCD3x3) {ntilex = ntiley = 3;}

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
  std::string DetectorType::TypeLabel() const
  //! return a string corresponding to the detector type
  {
    return TypeLabel(dettype);
  }
  //--------------------------------------------------------------
  std::string DetectorType::TypeLabel(const DetectorType::Type& dtype) const
  //! return a string corresponding to the detector type
  {
    std::string s;
    switch (dtype) {
    case UNKNOWN:
      s = "Unknown";
      break;
    case CCD1:  // one tile CCD
      s = "CCD-1-tile";
      break;
    case CCD2x2:
      s = "CCD-2x2-tile";
      break;
    case CCD3x3:
      s = "CCD-3x3-tile";
      break;
    case PILATUS6M:
      s = "Pilatus6M";
      break;
    case PILATUS2M:
      s = "Pilatus2M";
      break;
    default:
      s = "Unknown";
    }
    return s;
  }
  //--------------------------------------------------------------
  DetectorType::Type DetectorType::TypeFromLabel(const std::string& typelabel) const
  //! return a detector type corresponding to the string
  {
    DetectorType::Type t = UNKNOWN;

    if (typelabel == "Unknown") {
      t = UNKNOWN;
    } else if (typelabel == "CCD-1-tile") {  // one tile CCD
      t = CCD1;
    } else if (typelabel == "CCD-2x2-tile") {
      t = CCD2x2;
    } else if (typelabel == "CCD 3x3-tile") {
      t = CCD3x3;
    } else if (typelabel == "Pilatus6M") {
      t = PILATUS6M;
    } else if (typelabel == "Pilatus2M") {
      t = PILATUS2M;
    }
    return t;
  }
  //--------------------------------------------------------------
  bool DetectorType::equals(const DetectorType& b) const
  //! test for equality on type and detector range (not on number of tiles)
  {
    if (dettype != b.dettype) return false;
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
    //        <<sumwIxy.rows()<<" "<< sumwIxy.cols() <<" "
    //        <<sumwIothers.rows()<<" "<< sumwIothers.cols() <<" "
    //        <<avdelta.rows()<<" "<< avdelta.cols() <<"\n";
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
  DetectorScale::DetectorScale(const DetectorScaleType& DetScaleType,
                               const int& nTileX, const int& nTileY,
                               const DetectorType& dettype)
  {
    init(DetScaleType, nTileX, nTileY, dettype);
  }
  //--------------------------------------------------------------
  void DetectorScale::init(const DetectorScaleType& DetScaleType,
                           const int& nTileX, const int& nTileY,
                           const DetectorType& dettype)
  {
    detectorscaletype = DetScaleType;
    type = dettype;
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
  void DetectorScale::init(const DetectorScaleType& DetScaleType,
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
        } else if (detectorscaletype == CCD3) {
          tilescales(i,j) = new CCDTile3;
        } else if (detectorscaletype == CCD1) {
          tilescales(i,j) = new CCDTile1;
        } else if (detectorscaletype == CCD2) {
          tilescales(i,j) = new CCDTile2;
        } else if (detectorscaletype == PIXEL) {
          tilescales(i,j) = new TilePixel;
        }
        tilescales(i,j)->init(sizex, sizey);
        tilescales(i,j)->SetGridCoordinates(i,j);
      }}
    setSymmetric(false);
  }
  //--------------------------------------------------------------
  void DetectorScale::reinit()
  // reassign parameters and tilescales objects after Restore
  {
    for (int i=0;i<ntilex;++i) { // loop x
      for (int j=0;j<ntiley;++j) { // loop y
        if (tilescales(i,j) != NULL) {
          delete tilescales(i,j);
        }
      }
    }
    init();
  }
  //--------------------------------------------------------------
  void DetectorScale::clearCounts()
  {
    for (int i=0;i<ntilex;++i) { // loop x
      for (int j=0;j<ntiley;++j) { // loop y
        tilescales(i,j)->clearCounts();
      }}
  }
  //--------------------------------------------------------------
  void DetectorScale::setSymmetric(const bool& symmetric)
  {
    nparams = 0;
    for (int i=0;i<ntilex;++i) { // loop x
      for (int j=0;j<ntiley;++j) { // loop y
        tilescales(i,j)->setSymmetric(symmetric);
        idx_tile(i,j) = nparams;
        nparams += tilescales(i,j)->Nparams();
      }}
    //    std::cout << "DetectorScale::nparams " <<nparams <<"\n";  //^-
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
    } else if ((detectorscaletype == CCD1) ||
               (detectorscaletype == CCD2) ||
               (detectorscaletype == CCD3)) {
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
  std::vector<Tie> DetectorScale::Ties(const std::vector<double> tie_tile,
                                       const int& idx0) const
  // Return list of ties: tie_tile are sds for weight, idx0 is index to first global
  // parameter for setting ties, since they refer to the global parameter index
  //
  // For CCD tiles:
  //  tie_tile[0] for r
  //  tie_tile[1] for w
  //  tie_tile[2] for A
  //  tie_tile[3] for x0, y0
  //  tie_tile[4] for Fourier coefficients (if needed)
  {
    std::vector<Tie> ties;
    if (! ((detectorscaletype == CCD1) ||
           (detectorscaletype == CCD2) ||
           (detectorscaletype == CCD3)))
      {return ties;}
    // Only ties for CCD tiled detector (at present)
    int idx = idx0;
    // for ties between tiles:
    //  for each tile, a list of parameter indices and a list of weights
    std::vector<std::pair<std::vector<int>, std::vector<double> > > kindexwt;
     int kcentral = -1;
    if ((ntilex%2 != 0) && (ntiley%2 != 0)) {
      // there is a central tile if both nx & ny are odd
      kcentral = (ntilex/2)*ntiley+ntiley/2;
      if (kcentral == 0) kcentral = -1;
    }
    for (int i=0;i<ntilex;++i) { // loop tile x
      for (int j=0;j<ntiley;++j) { // loop tile y
        std::vector<Tie> tileties = tilescales(i,j)->Ties(tie_tile, idx);
        // ties within tile
        ties.insert(ties.end(), tileties.begin(), tileties.end());
        // Ties between tiles
        kindexwt.push_back(tilescales(i,j)->TiedParameters(tie_tile, idx));
        idx += tilescales(i,j)->Nparams(); // point to 1st parameter of next tile
      }}  // end tile loop

    int ntiles = ntilex*ntiley;
    ASSERT (int(kindexwt.size()) == ntiles);
    int npars = kindexwt[0].first.size();  // number of tied parameters/tile
    for (int j=0;j<npars;++j) { // loop parameters
      std::vector<int> kindex; // indices
      for (size_t k=0; k<kindexwt.size(); k++) {  // loop tiles
        kindex.push_back(kindexwt[k].first[j]);  // add in parameter index for tile
        if (int(k) == kcentral) { // overweight central tile
          kindex.push_back(kindexwt[k].first[j]);  // add in again
          kindex.push_back(kindexwt[k].first[j]);
        }
      }
      if (kindex.size() > 1) {ties.push_back(Tie(kindex, kindexwt[0].second[j]));}
    }
    //^
    //    std::cout << "Ties:\n";
    //    for (size_t j=0; j<ties.size(); j++) {
    //      std::cout << ties[j].format() <<"\n";
    //    }
    //^-

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
    //    std::cout << "xtl,ytl,Xt,Yt " << xtile <<" "<<ytile<<" "<<Xt<<" "<<Yt<<std::endl;
    // get scales & optional derivatives
    std::vector<double> dgdt;
    tilescales(xtile,ytile)->ScaleDeriv(Deriv, Xt, Yt, scale, dgdt);
    dgdp.assign(nparams,0.0);
    std::copy(dgdt.begin(), dgdt.end(), dgdp.begin()+idx_tile(xtile,ytile));
    //^
    //    if (DEBUG) {
    //      std::cout << "DetectorScale::ScaleDeriv "
    //          << Xt <<" "<<Yt<<"\n";
    //      for (size_t i=0;i<dgdp.size();++i) {std::cout <<" "<<dgdp[i];}
    //      std::cout <<"\n";
    //    }
    //^-
  }
  //--------------------------------------------------------------
    //! return number of parameters/tile
  int DetectorScale::NparamsTile() const
  {return tilescales(0,0)->Nparams();}
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
  std::string DetectorScale::formatparameters(const std::vector<double>& sds) const
  //! format parameters for printing
  {
    std::string s = "\n";
    if (detectorscaletype == NONE) {return s;}

    s += "Detector scales for"+
    StringUtil::itos(ntilex,2)+" x"+StringUtil::itos(ntiley,2)+
      " tiles\n";
    s += tilescales(ntilex/2, ntiley/2)->formattype()+"\n";

    if (detectorscaletype != PIXEL) {
      s += "    Tile parameters arranged with Xdet across and Ydet down\n";
      if (detectorscaletype == CCD1 ||detectorscaletype == CCD2 ||detectorscaletype == CCD3) {
        s += "    Numbers in the corners are the number of observations contributing to the scales\n";
      }
    }
    int npartile = NparamsTile() ;

    // tile scale
    for (int j=0;j<ntiley;++j) { // loop y
      // Get formatted squares for x tiles
      // these should all have the same number of lines
      std::vector<std::vector<std::string> > sxtiles(ntilex);
      for (int i=0;i<ntilex;++i) { // loop x
        // SDs for this tile: note that in parameter list tile loops are x, then y (fast)
        //  so tile number is i*ntilex + j
        int idx1 = (i*ntilex + j) * npartile;
        std::vector<double> sdstile(sds.begin()+idx1, sds.begin()+idx1+npartile);
        sxtiles[i] = tilescales(i,j)->formatparameters(sdstile);
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
  std::string DetectorScale::formatTies() const
  {
    return tilescales(ntilex/2, ntiley/2)->formatTies();
  }
  //--------------------------------------------------------------
  std::string DetectorScale::FormatSave() const
  // return formatted version for save and restore
  {
    std::string dump = "DetectorScale V1 {\n";
    dump += "DetectorType " + type.TypeLabel() + "\n";
    dump += "ScaleType " + formatType(detectorscaletype) + "\n";
    dump += "Ntilex "+itos(ntilex)+"\n";
    dump += "Ntiley "+itos(ntiley)+"\n";
    dump += "XDrange "+ftos(xdrange.min())+" "+ftos(xdrange.max()) + "\n";
    dump += "YDrange "+ftos(ydrange.min())+" "+ftos(ydrange.max()) + "\n";
    for (int j=0;j<ntiley;++j) { // loop y
      for (int i=0;i<ntilex;++i) { // loop x
        dump += "TileXY " + itos(i) + " " + itos(j) + "\n";
        dump += tilescales(i,j)->FormatSave();
        dump += "ParameterIndex " +itos(idx_tile(i,j)) + "\n";
      } // x
    } // y
    return dump+"}\n";
  }
  //--------------------------------------------------------------
  // restore
  void DetectorScale::Restore(Fileread& FR)
  {
    FR.ReadTag("DetectorScale"); // fails if tag does not match
    if (FR.GetTag() != "V1") {  // version check
      clipper::Message::message(Message_fatal
        ("DetectorScale::Restore incompatible version in "+FR.Filename()));
    }
    FR.Skip();
    std::string typelabel;
    FR.ReadTag("DetectorType"); typelabel = FR.GetTag();
    std::string detectorscaletypelabel;
    FR.ReadTag("ScaleType"); detectorscaletypelabel = FR.GetTag();
    detectorscaletype = Type(detectorscaletypelabel);
    FR.ReadTag("Ntilex"); ntilex = FR.Int();
    FR.ReadTag("Ntiley"); ntiley = FR.Int();
    double a1, a2;
    FR.ReadTag("XDrange"); a1 = FR.Double(); a2 = FR.Double();
    xdrange = Range(a1, a2);
    FR.ReadTag("YDrange"); a1 = FR.Double(); a2 = FR.Double();
    ydrange = Range(a1, a2);
    tilescales.resize(ntilex, ntiley);
    idx_tile.resize(ntilex, ntiley);
    reinit(); // recreate tilescales objects
    int ix, jy;
    for (int j=0;j<ntiley;++j) { // loop y
      for (int i=0;i<ntilex;++i) { // loop x
        FR.ReadTag("TileXY"); ix = FR.Int(); jy = FR.Int();
        ASSERT ((ix == i) && (jy == j));
        tilescales(i,j)->Restore(FR);
        FR.ReadTag("ParameterIndex"); idx_tile(i,j) = FR.Int();
      } // x
    } // y
    if (!FR.CheckEnd()) {
      clipper::Message::message(Message_warn
        ("DetectorScale::Restore unexpected tag "+FR.Tag()));
    }
    std::vector<std::vector<float> > detrange(2);
    detrange[0].resize(2);
    detrange[1].resize(2);
    detrange[0][0] = xdrange.min();
    detrange[0][1] = xdrange.max();
    detrange[1][0] = ydrange.min();
    detrange[1][1] = ydrange.max();
    type = DetectorType(DetectorType::UNKNOWN, typelabel, detrange);
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
  std::string DetectorScale::formatType(const  DetectorScaleType& type)
  {
    if (type == NONE) {return "NONE";}
    if (type == FLAT) {return "FLAT";}
    if (type == CCD1) {return "CCD1";}
    if (type == CCD2) {return "CCD2";}
    if (type == CCD3) {return "CCD3";}
    if (type == PIXEL) {return "PIXEL";}
    if (type == AUTOMATIC) {return "AUTOMATIC";}
    return "";
  }
  //--------------------------------------------------------------
  DetectorScale::DetectorScaleType DetectorScale::Type(const std::string& scaletypelabel)
  {
    if (scaletypelabel == "NONE") {return NONE;}
    if (scaletypelabel == "FLAT") {return FLAT;}
    if (scaletypelabel == "CCD1") {return CCD1;}
    if (scaletypelabel == "CCD2") {return CCD2;}
    if (scaletypelabel == "CCD3") {return CCD3;}
    if (scaletypelabel == "PIXEL") {return PIXEL;}
    if (scaletypelabel == "AUTOMATIC") {return AUTOMATIC;}
    return NONE;
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  CCDTile3::CCDTile3(const double& Xmax, const double& Ymax)
  {
    init(Xmax, Ymax);
  }
  //--------------------------------------------------------------
  void CCDTile3::init(const double& Xmax, const double& Ymax)
  // coordinates in range 0->Xmax, 0->Ymax
  {
    xmax = Xmax;
    ymax = Ymax;
    nparams = 15;
    twooverrootpi = 2.0/sqrt(clipper::Util::pi()); // 2/sqrt(pi)

    // Centre of coordinate system (for ties)
    xc0 = xmax/2.0;
    yc0 = ymax/2.0;
    rad0 = Max(xc0, yc0);  // larger edge radius, pixels
    // Set default (initial) values
    x0 = 0.0;
    y0 = 0.0;

    rfs.setIsConstant(false);
    r0 = 0.7;
    rfs.setLevel(0.0);  // relative to rad0
    wfs.setIsConstant(false);
    w0 = 0.2;
    wfs.setLevel(0.0);
    Afs.setIsConstant(false);
    A0 = 0.1;
    Afs.setLevel(0.0);
    nparams_smooth = rfs.NumberParameters(); // for each of r,w,A

    ncorners.resize(2,2,0); // counts in corners, initialise to 0
    // dcrnmin = sqrt(1/2((xmax/2)^2+(ymax/2)^2)) limit for corner
    dcrnmin = (sqrt(0.5*0.25*(xmax*xmax + ymax*ymax)))/rad0;
  }
  //--------------------------------------------------------------
  void CCDTile3::clearCounts()
  {
    ncorners.resize(2,2,0); // counts in corners, initialise to 0
  }
  //--------------------------------------------------------------
  //! symmetric = false to allow A to vary around the tile
  void CCDTile3::setSymmetric(const bool& symmetric)
  {
    circularlysymmetric = symmetric;
    rfs.setIsConstant(false); // no constant term in Fourier smoothing
    rfs.setLevel(0.0);
    wfs.setIsConstant(false); // no constant term in Fourier smoothing
    wfs.setLevel(0.0);
    Afs.setIsConstant(false); // no constant term in Fourier smoothing
    Afs.setLevel(0.0);
    if (circularlysymmetric) {
      nparams_smooth = 0; // for r,w,A
    } else {
      nparams_smooth = Afs.NumberParameters(); // for r,w,A
    }
    nparams = 3*nparams_smooth + 3;  // r,w,A,Asmooth(4)
  }
  //--------------------------------------------------------------
  // Parameter order: r0, w0, A0, rfs(4), wfs(4), Afs(4)
  // Store parameters
  void CCDTile3::StoreParameters(const std::vector<double>& parameters)
  {
    ASSERT (int(parameters.size()) == nparams);

    r0 = parameters[0];
    w0 = parameters[1];
    A0 = parameters[2];
    if (!circularlysymmetric) {
      std::vector<double>::const_iterator pp = parameters.begin()+3;
      std::vector<double> r_params(pp, pp+nparams_smooth);
      pp += nparams_smooth;
      std::vector<double> w_params(pp, pp+nparams_smooth);
      pp += nparams_smooth;
      std::vector<double> A_params(pp, pp+nparams_smooth);
      rfs.setParameters(r_params);
      wfs.setParameters(w_params);
      Afs.setParameters(A_params);
    }
    //^
    //    std::cout <<"\nr : " <<r0<<" "<< rfs.format() <<"\n";
    //    std::cout <<"w : " <<w0<<" "<< wfs.format() <<"\n";
    //    std::cout <<"A : " <<A0<<" "<< Afs.format() <<"\n";
    //    std::cout <<"\n";  //^-
  }
  //--------------------------------------------------------------
  // Retrieve parameter vector (length nparams)
  // Parameter order: r0, w0, A0, rfs(4), wfs(4), Afs(4)
  std::vector<double> CCDTile3::Parameters() const
  {
    std::vector<double> par;
    par.push_back(r0);
    par.push_back(w0);
    par.push_back(A0);
    if (!circularlysymmetric) {
      std::vector<double> r_params = rfs.GetParameters();
      std::vector<double> w_params = wfs.GetParameters();
      std::vector<double> A_params = Afs.GetParameters();
      par.insert(par.end(), r_params.begin(), r_params.end());
      par.insert(par.end(), w_params.begin(), w_params.end());
      par.insert(par.end(), A_params.begin(), A_params.end());
    }
    ASSERT (int(par.size()) == nparams);
    return par;
  }
  //--------------------------------------------------------------
  //! return vector of internal ties, given SDs and 1st global parameter index
  std::vector<Tie> CCDTile3::Ties(const std::vector<double> tie_tile,
                                   const int& idx0)
  {
    // 12 parameters expressed as Fourier series, r,w relative to rad0
    // For all parameter types r,w,A, tie all Fourier coefficients to 0.0
    // no ties if SD = 0
    std::vector<Tie> ties;
    if (!circularlysymmetric && tie_tile[4] > 0.0) {
      int nparams_tile = Nparams();
      double SD = tie_tile[4];
      double weight = 1.0/(SD*SD);
      int idx = idx0+3;  // skip r0, w0, A0

      for (int k=3;k<nparams_tile;++k) { // loop parameters / tile (=15-3)
        // tie ABCD to 0.0
        ties.push_back(Tie(idx, 0.0, weight));
        idx++;
      }
    }
    return ties;
  }
  //--------------------------------------------------------------
  //! vector of indices and weights for each parameter to be restrained across tiles
  std::pair<std::vector<int>, std::vector<double> >
  CCDTile3::TiedParameters(const std::vector<double> tie_tile,
                 const int& idx0)
  {
    // tie r0,w0,A0 tiles
    std::pair<std::vector<int>, std::vector<double> > kindexwt;

    int idx = idx0;

    for (int k=0;k<3;++k) { // loop parameters 0,1,2 = r0,w0,A0
      if (tie_tile[k] > 0.0) {
        double weight = 1.0/(tie_tile[k]*tie_tile[k]);
        kindexwt.first.push_back(idx);
        kindexwt.second.push_back(weight);
      }
      idx++;
    }
    return kindexwt;
  }
  //--------------------------------------------------------------
  //! lower bounds for ipar'th parameter
  double CCDTile3::LowerBound(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    if (ipar == 0) { // r0
      return 0.1;
    } else if (ipar == 1) { // w0
      return 0.1;
    } else if (ipar == 2) { // A0
      return -0.001;
    }
    return -1.0;  // Fourier terms
  }
  //--------------------------------------------------------------
  //! upper bounds for ipar'th parameter
  double CCDTile3::UpperBound(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    if (ipar == 0) { // r0
      return 1.41;
    } else if (ipar == 1) { // w0
      return 1.0;
    } else if (ipar == 2) { // A0
      return 0.5;
    }
    return +1.0;  // Fourier terms
  }
  //--------------------------------------------------------------
  //! "large shift" for ipar'th parameter
  double CCDTile3::LargeShift(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    if (ipar == 0) { // r0
      return 0.3;
    } else if (ipar == 1) { // w0
      return 0.2;
    } else if (ipar == 2) { // A0
      return 0.2;
    }
    return 0.05;  // Fourier terms
  }
  //--------------------------------------------------------------
  // Return scale & derivatives for tile coordinates Xt, Yt
  void CCDTile3::ScaleDeriv(const bool& Deriv,
                            const double& Xt, const double& Yt,
                            double& scale,
                            std::vector<double>& dgdp) const
  {
    // Xt, Yt in pixels
    double x = (Xt-xc0)/rad0;
    double y = (Yt-yc0)/rad0;

    double d = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)); // relative to rad0
    double phi = atan2((y-y0), (x-x0));

    // For Fourier coefficients
    std::vector<double> drdp;
    std::vector<double> dwdp;
    std::vector<double> dAdp;
    double fr = 0.0;
    double fw = 0.0;
    double fA = 0.0;
    if (!circularlysymmetric) {
      fr = rfs.ValueDerivatives(phi, drdp); // rfs and dr/dp vector
      fw = wfs.ValueDerivatives(phi, dwdp); // w and dw/dp vector
      fA = Afs.ValueDerivatives(phi, dAdp); // A and dA/dp vector
    }
    double r = r0*(1.0+fr);
    double w = w0*(1.0+fw);
    double A = A0*(1.0+fA);

    double z = 2.0*(d-r-w)/w;
    scale = radfunc.value(z, A);

    // Count observations in each corner
    if (d > dcrnmin) { // in a corner
      int i = (x-x0) < 0.0 ? 0 : 1; // 0 or 1 if left or right
      int j = (y-y0) < 0.0 ? 0 : 1; // 0 or 1 if bottom or top
      ncorners(i,j)++;
    }

    if (Deriv) {
      const int NPARAMBASE = 3;
      dgdp = radfunc.deriv(NPARAMBASE, w);  // fills 1st 3 slots in dgdp
      if (!circularlysymmetric) {
        // d/d(Fourier coefficients)
        for (size_t i=0; i<drdp.size(); i++) {
          // dg/dr(ABCD) = dg/dr * r0 * df(phi)/dr(ABCD)
          drdp[i] *= r0*dgdp[0];
          dwdp[i] *= w0*dgdp[1];  //  etc
          dAdp[i] *= A0*dgdp[2];
        }
        // dg/dq0 = dg/dq dq/dq0 = dg/dq (1+fq)  for q=r,w,A
        dgdp[0] *= (1.0+fr); // dg/dr0
        dgdp[1] *= (1.0+fw); // dg/dw0
        dgdp[2] *= (1.0+fA); // dg/dA0
        dgdp.insert(dgdp.end(), drdp.begin(), drdp.end());
        dgdp.insert(dgdp.end(), dwdp.begin(), dwdp.end());
        dgdp.insert(dgdp.end(), dAdp.begin(), dAdp.end());
      }
      ASSERT (int(dgdp.size()) == nparams);
    }
  }
  //--------------------------------------------------------------
  std::string CCDTile3::format() const
  //! format scale type for printing
  {
    std::string text = "Tile correction for CCD detector";
    return text;
  }
  //--------------------------------------------------------------
  std::string CCDTile3::formattype() const
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
    text += "   r,w and A are each parameterised as a constant and a 4-parameter Fourier series\n";
    text += "     eg r = r0 * (1 + f(phi))\n";
    return text;
  }
  //--------------------------------------------------------------
  std::vector<std::string>
  CCDTile3::format5(const bool& hasSd, const double& v0, const double& sd,
                    const FourierSmooth& vfs, const std::vector<double>& sdfs,
                    const int& width, const std::string& label) const
  // format v0{vfs} with optional sds on 2nd line
  {
    std::vector<std::string> lines;
    std::string line;
    if (hasSd) {
      line = StringUtil::Strip(StringUtil::ftos(v0,6,2));
    } else {
      line = StringUtil::Strip(StringUtil::valueSD(v0,sd,0,6,2));
    }
    line += vfs.format();
    line = "| "+ StringUtil::CentreString(label+"="+line, width-4) + " |";
    lines.push_back(line);
    if (hasSd) {
      std::string s;
      int nsmooth = sdfs.size();
      for (int i=0;i<nsmooth;++i) {
        s += StringUtil::ftos(sdfs[i],6,2);
        if (i < nsmooth-1) {s += ",";}
      }
      line = "| "+ StringUtil::CentreString("SDs("+StringUtil::Strip(s)+")", width-4) + " |";
      lines.push_back(line);
    }

    return lines;
  }
  //--------------------------------------------------------------
  std::vector<std::string> CCDTile3::formatparameters(const std::vector<double>& sds) const
  //! format parameters for printing
  // Parameter order: r0, w0, A0, rfs(4), wfs(4), Afs(4)
  {
    int nsds = sds.size(); // number of SDs given, if any
    if (nsds > 0 ) {
      ASSERT (nsds == nparams);
    }
    std::vector<std::string> s;
    std::string line;
    int width = 37; // width of window
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

    std::vector<std::string> lines;
    // r, w, A
    int idxfs = 3;  // index for first radial variation parameter (rfs)
    int nfs = 4;    // number of fs parameters
    bool hasSd = (nsds > 0);
    double sd = 0.0;
    std::vector<double> sdfs;
    if (hasSd) {
      sd = sds[0];  // r0
      sdfs.assign(sds.begin()+idxfs, sds.begin()+idxfs+nfs);
    }
    lines = format5(hasSd, r0, sd, rfs, sdfs, width, "r");
    s.insert(s.end(),lines.begin(), lines.end());

    if (hasSd) {
      sd = sds[1]; // w0
      idxfs += nfs;
      sdfs.assign(sds.begin()+idxfs, sds.begin()+idxfs+nfs);
    }
    lines = format5(hasSd, r0, sd, wfs, sdfs, width, "w");
    s.insert(s.end(),lines.begin(), lines.end());

    if (hasSd) {
      sd = sds[2]; // A0
      idxfs += nfs;
      sdfs.assign(sds.begin()+idxfs, sds.begin()+idxfs+nfs);
    }
    lines = format5(hasSd, A0, sd, Afs, sdfs, width, "A");
    s.insert(s.end(),lines.begin(), lines.end());

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
  // Format all information into a labelled save format for later restoration
  std::string CCDTile3::FormatSave() const
  {
    std::string dump = "CCDTile3 V1 {\n"; // with version number
    dump += "XYmax " + ftos(xmax) + " " + ftos(ymax) + "\n";
    dump += "Rad0 " + ftos(rad0) + "\n";
    if (circularlysymmetric) {
      dump += "Symmetric\n";
    } else {
      dump += "Nonsymmetric\n";
    }
    dump += "Nparams " + itos(nparams) + " " + itos(nparams_smooth) + "\n";
    dump += "XYc0 " + ftos(xc0) + " " + ftos(yc0) + "\n";
    dump += "XY0 " + ftos(x0) + " " + ftos(y0) + "\n";
    dump += "Parameters\n" + StringUtil::FormatSaveVector(Parameters());
    dump += "Ncorners " + itos(ncorners(0,0)) + " " + itos(ncorners(0,1))
      + " " + itos(ncorners(1,0)) + " " + itos(ncorners(1,1)) + "\n";
    return dump+"}\n";
  }
  //--------------------------------------------------------------
  void CCDTile3::Restore(Fileread& FR)
  // restore
  {
    FR.ReadTag("CCDTile3"); // fails if tag does not match
    if (FR.GetTag() != "V1") {  // version check
      clipper::Message::message(Message_fatal
        ("CCDTile3::Restore incompatible version in "+FR.Filename()));
    }
    FR.Skip();
    FR.ReadTag("XYmax"); xmax = FR.Double(); ymax = FR.Double();
    FR.ReadTag("Rad0"); rad0 = FR.Double();
    std::string symm = FR.GetTag();
    circularlysymmetric = true;
    if (symm == "Nonsymmetric") {
      circularlysymmetric = false;
    }
    FR.ReadTag("Nparams"); nparams = FR.Int(); nparams_smooth = FR.Int();
    FR.ReadTag("XYc0"); xc0 = FR.Double(); yc0 = FR.Double();
    FR.ReadTag("XY0"); x0 = FR.Double(); y0 = FR.Double();
    FR.ReadTag("Parameters");
    std::vector<double> parameters = FR.DoubleVec(nparams);
    StoreParameters(parameters);
    FR.ReadTag("Ncorners");
    ncorners(0,0) = FR.Int(); ncorners(0,1) = FR.Int();
    ncorners(1,0) = FR.Int(); ncorners(1,1) = FR.Int();
    if (!FR.CheckEnd()) {
      clipper::Message::message(Message_warn
        ("CCDTile3::Restore unexpected tag "+FR.Tag()));
    }
  }
  //--------------------------------------------------------------
  CCDTile1::CCDTile1(const double& Xmax, const double& Ymax)
  {
    init(Xmax, Ymax);
  }
  //--------------------------------------------------------------
  void CCDTile1::init(const double& Xmax, const double& Ymax)
  // coordinates in range 0->Xmax, 0->Ymax
  // symmetric correction, function of x & y separately
  {
    xmax = Xmax;
    ymax = Ymax;
    nparams = 5;
    twooverrootpi = 2.0/sqrt(clipper::Util::pi()); // 2/sqrt(pi)

    // Centre of coordinate system (for ties)
    xc0 = xmax/2.0;
    yc0 = ymax/2.0;
    rad0 = Max(xc0, yc0);  // larger edge radius, pixels
    // Set default (initial) values
    x0 = 0.0;
    y0 = 0.0;

    r = 0.6;
    w = 0.4;
    A = 0.2;
    ncorners.resize(2,2,0); // counts in corners, initialise to 0
    // dcrnmin = sqrt(1/2((xmax/2)^2+(ymax/2)^2)) limit for corner
    dcrnmin = 0.5*0.25*(xmax*xmax + ymax*ymax)/(rad0*rad0);
  }
  //--------------------------------------------------------------
  void CCDTile1::clearCounts()
  {
    ncorners.resize(2,2,0); // counts in corners, initialise to 0
  }
  //--------------------------------------------------------------
  // Parameter order: r,w,A,x0,y0,
  // Store parameters
  void CCDTile1::StoreParameters(const std::vector<double>& parameters)
  {
    ASSERT (int(parameters.size()) == nparams);
    r = parameters[0];
    w = parameters[1];
    A = parameters[2];
    x0 = parameters[3];
    y0 = parameters[4];
    //^
    //    std::cout << "r,w,A,x0,y0 ";
    //    for (int i=0;i<5;++i) {std::cout <<" "<<parameters[i];}
    //    std::cout <<"\n";  //^-
  }
  //--------------------------------------------------------------
  // Retrieve parameter vector (length nparams)
  std::vector<double> CCDTile1::Parameters() const
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
  //! return vector of ties, given SDs and 1st global parameter index
  std::vector<Tie> CCDTile1::Ties(const std::vector<double> tie_tile,
                                   const int& idx0)
  {
    // Note that r,w,x0,y0 etc are in fractions of rad0
    std::vector<Tie> ties;
    if (tie_tile[3] > 0.0) {
      // x0, y0 tie to centre position
      // tie_tile[3] is relative to tile size
      int idx = idx0 + 3;  // first x0 parameter
      double weight = tie_tile[3];
      weight = 1./(weight*weight);
      // tie x0 to centre
      ties.push_back(Tie(idx, 0.0, weight));
      // tie y0 to centre
      ties.push_back(Tie(idx+1, 0.0, weight));
    }
    return ties;
  }
  //--------------------------------------------------------------
  //! vector of indices and weights for each parameter to be restrained across tiles
  std::pair<std::vector<int>, std::vector<double> >
  CCDTile1::TiedParameters(const std::vector<double> tie_tile,
                 const int& idx0)
  {
    std::pair<std::vector<int>, std::vector<double> > kindexwt;
    // Tie r,w,A parameters together for all tiles
    int idx = idx0;   // starting global parameter index
    for (int k=0;k<3;++k) { // loop parameters 0,1,2 = r,w,A
      if (tie_tile[k] > 0.0) {
        double weight = tie_tile[k];
        weight = 1./(weight*weight);
        kindexwt.second.push_back(weight);
        std::cout << k <<" tie for r,w,A = 0,1,2\n";
        std::vector<int> kindex;  // parameter indices for this group of parameters
        kindexwt.first.push_back(idx);
      }
      idx++;
    } // end loop r,w,A
    return kindexwt;
  }
  //--------------------------------------------------------------
  //! lower bounds for ipar'th parameter
  double CCDTile1::LowerBound(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    if (ipar <= 1) {
      return 0.1; // r or w
    } else if (ipar == 2 ) {
      return -0.001; // A
    }
    return -0.5;  // x0 or y0
  }
  //--------------------------------------------------------------
  //! upper bounds for ipar'th parameter
  double CCDTile1::UpperBound(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    if (ipar == 0) { // r
      return 1.41;
    } else if (ipar == 1) { // w
      return 1.0;
    } else if (ipar == 2 ) {
      return 0.5; // A
    }
    return +0.5;  // x0 or y0
  }
  //--------------------------------------------------------------
  //! "large shift" for ipar'th parameter
  double CCDTile1::LargeShift(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    if (ipar == 0) { // r
      return 0.3;
    } else if (ipar == 1) { // w
      return 0.2;
    } else if (ipar == 2) { // A
      return 0.2;
    }
    return 0.4;  // x0 or y0
  }
  //--------------------------------------------------------------
  // Return scale & derivatives for tile coordinates Xt, Yt
  void CCDTile1::ScaleDeriv(const bool& Deriv,
                            const double& Xt, const double& Yt,
                            double& scale,
                            std::vector<double>& dgdp) const
  {
    // symmetric correction, function of x & y separately
    // Xt, Yt in pixels
    double x = (Xt-xc0)/rad0;  // relative to centre, scaled by rad0
    double y = (Yt-yc0)/rad0;  // ie in range -1 to +1

    double dx = x-x0;
    double dy = y-y0;
    double d2 = dx*dx + dy*dy;
    double d = sqrt(d2);
    double z = 2.0*(d-r-w)/w;
    scale = radfunc.value(z, A);
    /*
    double expmz = exp(-z); // also = dg/dz
    // scale = A f(z) + 1 - A
    // f(z) = 1 - exp(g(z))
    // g(z) = -exp(-z)-1
    double gz = -expmz - 1.0;
    double fz = 1.0 - exp(gz);
    scale = A * fz + 1.0 - A;
    */

    // Count observations in each corner
    if (d2 > dcrnmin) { // in a corner
      int i = dx < 0.0 ? 0 : 1; // 0 or 1 if left or right
      int j = dy < 0.0 ? 0 : 1; // 0 or 1 if bottom or top
      ncorners(i,j)++;
    }

    if (Deriv) {
      dgdp = radfunc.deriv(nparams, w);  // fills 1st 3 slots in dgdp
      /*
      // df/dz = exp(gz)exp(-z); ds/dz = A df/dz
      double dsdz = - A * exp(gz) * expmz;
      dgdp[0] = - dsdz * 2.0/w;          // ds/dr = -ds/dd
      dgdp[1] = dgdp[0] *(1.0 + 0.5*z);  // ds/dw
      dgdp[2] = fz - 1.0;                // ds/dA
      */
      // ds/dx0 = ds/dd dd/dx0
      // dd/dx0 =
      if (d == 0.0) {
        dgdp[3] = dgdp[0]; // eg if y=y0, d=x-x0, dd/dx0 = -1
        dgdp[4] = dgdp[0]; //  dd/dy0 = -1
      } else {
        dgdp[3] = - dgdp[0] * (x0 - x) / d;  // dgdx0 = dg/dd * dd/dx0
        dgdp[4] = - dgdp[0] * (y0 - y) / d;  // dgdy0 = dg/dd * dd/dy0
      }
    }
  }
  //--------------------------------------------------------------
  std::string CCDTile1::format() const
  //! format scale type for printing
  {
    std::string text = "Tile correction for CCD detector";
    return text;
  }
  //--------------------------------------------------------------
  std::string CCDTile1::formattype() const
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
  std::vector<std::string> CCDTile1::formatparameters(const std::vector<double>& sds) const
  //! format parameters for printing
  // Parameter order r,w,A,x0,y0
  {
    int nsds = sds.size(); // number of SDs given, if any
    if (nsds > 0 ) {
      ASSERT (nsds == nparams);
    }
    std::vector<std::string> s;
    std::string line;
    int width = 34; // width of window
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
    if (nsds == 0) {
      line = "r="+StringUtil::ftos(r,6,2)+
        ", w="+StringUtil::ftos(w,6,2)+
        ", A="+StringUtil::ftos(A,6,2);
    } else {
      line = "r="+StringUtil::valueSD(r,sds[0],0,6,2)+
        ", w="+StringUtil::valueSD(w,sds[1],0,6,2)+
        ", A="+StringUtil::valueSD(A,sds[2],0,6,2);
    }
    line = "| "+ StringUtil::CentreString(StringUtil::Strip(line), width-4)+ " |";
    s.push_back(line);
    if (nsds == 0) {
      line ="Tile centre: "+StringUtil::Strip(StringUtil::ftos(x0,7,2))+
        ", "+StringUtil::Strip(StringUtil::ftos(y0,7,2));
    } else {
      line ="Tile centre: "+StringUtil::Strip(StringUtil::valueSD(x0,sds[3],0,7,2))+
        ", "+StringUtil::Strip(StringUtil::valueSD(y0,sds[4],0,7,2));
    }
    line = "| "+ StringUtil::CentreString(line, width-4)+ " |";
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
  // Format all information into a labelled save format for later restoration
  std::string CCDTile1::FormatSave() const
  {
    std::string dump = "CCDTile1 V1 {\n"; // with version number
    dump += "XYmax " + ftos(xmax) + " " + ftos(ymax) + "\n";
    dump += "Rad0 " + ftos(rad0) + "\n";
    dump += "Nparams " + itos(nparams) + "\n";
    dump += "XYc0 " + ftos(xc0) + " " + ftos(yc0) + "\n";
    dump += "XY0 " + ftos(x0) + " " + ftos(y0) + "\n";
    dump += "Parameters\n" + StringUtil::FormatSaveVector(Parameters());
    dump += "Ncorners " + itos(ncorners(0,0)) + " " + itos(ncorners(0,1))
      + " " + itos(ncorners(1,0)) + " " + itos(ncorners(1,1)) + "\n";
    return dump+"}\n";
  }
  //--------------------------------------------------------------
  void CCDTile1::Restore(Fileread& FR)
  // restore
  {
    FR.ReadTag("CCDTile1"); // fails if tag does not match
    if (FR.GetTag() != "V1") {  // version check
      clipper::Message::message(Message_fatal
        ("CCDTile1::Restore incompatible version in "+FR.Filename()));
    }
    FR.Skip();
    FR.ReadTag("XYmax"); xmax = FR.Double(); ymax = FR.Double();
    FR.ReadTag("Rad0"); rad0 = FR.Double();
    FR.ReadTag("Nparams"); nparams = FR.Int();
    FR.ReadTag("XYc0"); xc0 = FR.Double(); yc0 = FR.Double();
    FR.ReadTag("XY0"); x0 = FR.Double(); y0 = FR.Double();
    FR.ReadTag("Parameters");
    std::vector<double> parameters = FR.DoubleVec(nparams);
    StoreParameters(parameters);
    FR.ReadTag("Ncorners");
    ncorners(0,0) = FR.Int(); ncorners(0,1) = FR.Int();
    ncorners(1,0) = FR.Int(); ncorners(1,1) = FR.Int();
    if (!FR.CheckEnd()) {
      clipper::Message::message(Message_warn
        ("CCDTile1::Restore unexpected tag "+FR.Tag()));
    }
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  CCDTile2::CCDTile2(const double& Xmax, const double& Ymax)
  {
    init(Xmax, Ymax);
  }
  //--------------------------------------------------------------
  void CCDTile2::init(const double& Xmax, const double& Ymax)
  // coordinates in range 0->Xmax, 0->Ymax
  {
    xmax = Xmax;
    ymax = Ymax;
    twooverrootpi = 2.0/sqrt(clipper::Util::pi()); // 2/sqrt(pi)

    // Centre of coordinate system (for ties)
    xc0 = xmax/2.0;
    yc0 = ymax/2.0;
    rad0 = Max(xc0, yc0);  // larger edge radius, pixels
    // Set default (initial) values
    x0 = 0.0;
    y0 = 0.0;

    r = 0.7;  // relative to rad0
    w = 0.4;
    A0 = 0.2;
    setSymmetric(false);  // set default to allow A to vary
    ties_.assign(7,1.0);

    ncorners.resize(2,2,0); // counts in corners, initialise to 0
    // dcrnmin = sqrt(1/2((xmax/2)^2+(ymax/2)^2)) limit for corner
    dcrnmin = (sqrt(0.5*0.25*(xmax*xmax + ymax*ymax)))/rad0;
  }
  //--------------------------------------------------------------
  void CCDTile2::clearCounts()
  {
    ncorners.resize(2,2,0); // counts in corners, initialise to 0
  }
  //--------------------------------------------------------------
  //! symmetric = false to allow A to vary around the tile
  void CCDTile2::setSymmetric(const bool& symmetric)
  {
    circularlysymmetric = symmetric;
    Afs.setIsConstant(false); // no constant term in Fourier smoothing
    Afs.setLevel(0.0);
    if (circularlysymmetric) {
      nparams_smooth = 0; // for A
    } else {
      nparams_smooth = Afs.NumberParameters(); // for A
    }
    nparams = nparams_smooth + 5;  // r,w,x0,y0,A,Asmooth(4)
  }
  //--------------------------------------------------------------
  // Parameter order: r,w,A0,x0,y0
  // Store parameters
  void CCDTile2::StoreParameters(const std::vector<double>& parameters)
  {
    ASSERT (int(parameters.size()) == nparams);
    r = parameters[0];
    w = parameters[1];
    A0 = parameters[2];
    x0 = parameters[3];
    y0 = parameters[4];

    if (!circularlysymmetric) {
      std::vector<double>::const_iterator pp = parameters.begin() + 5;
      std::vector<double> A_params(pp, pp+nparams_smooth);
      Afs.setParameters(A_params);
    }
    //^
    //    std::cout <<"Tile "<<ix<<","<<iy<<", r : " << r;
    //    std::cout <<"  w : " << w;
    //    std::cout <<"  A0 : " <<A0<<" "<< Afs.format() <<"\n";
    //    std::cout << "x0, y0: " << x0 <<" " << y0<<"\n";
    //    std::cout <<"\n";  //^-
  }
  //--------------------------------------------------------------
  // Retrieve parameter vector (length nparams)
  std::vector<double> CCDTile2::Parameters() const
  {
    std::vector<double> par(nparams);
    par[0] = r;
    par[1] = w;
    par[2] = A0;
    par[3] = x0;
    par[4] = y0;
    if (!circularlysymmetric) {
      std::vector<double> A_params = Afs.GetParameters();
      for (size_t j=0; j<A_params.size(); j++) {
        par[j+5] = A_params[j];
      }
    }
    ASSERT (int(par.size()) == nparams);
    return par;
  }
  //--------------------------------------------------------------
  //! return vector of internal ties, given SDs and 1st global parameter index
  std::vector<Tie> CCDTile2::Ties(const std::vector<double> tie_tile,
                                   const int& idx0)
  {
    std::vector<Tie> ties;
    ties_ = tie_tile;
    int idx = idx0;  // don't skip r, w
    //    int idx = idx0 + 2;  // skip r, w
    double weight;

    // ties for r and w to target values
    double rtarget = tie_tile[5];
    if (rtarget > 0.0) {
      weight = 1.0/(4.0*tie_tile[0]*tie_tile[0]);  // r weight
      ties.push_back(Tie(idx, rtarget, weight));
    }
    idx++;
    double wtarget = tie_tile[6];
    if (wtarget > 0.0) {
      weight = 1.0/(4.0*tie_tile[1]*tie_tile[1]);  // w weight
      ties.push_back(Tie(idx, wtarget, weight));
    }
    idx++;

    // Tie A0 to 0.0
    if (tie_tile[2] > 0.0) {
      weight = 1.0/(4.0*tie_tile[2]*tie_tile[2]);  // sd*2 for tie to zero
      ties.push_back(Tie(idx, 0.0, weight));
    }
    idx++;

    // x0, y0 tie to centre position
    // tie_tile[3] is relative to tile size as are x0, y0 parameters
    if (tie_tile[3] > 0.0) {
      weight = tie_tile[3]; // rad0 in pixels
      weight = 1./(weight*weight);
      // tie x0 to centre
      ties.push_back(Tie(idx++, 0.0, weight));
      // tie y0 to centre
      ties.push_back(Tie(idx++, 0.0, weight));
    } else {
      idx += 2;
    }
    if (!circularlysymmetric && tie_tile[4] > 0.0) {
      weight = 1.0/(tie_tile[4]*tie_tile[4]);
      // For parameter type A, tie all Fourier coefficients to 0.0
      for (int j=0;j<nparams_smooth;++j) { // ABCD
        // tie ABCD to 0.0
        ties.push_back(Tie(idx++, 0.0, weight));
      }
    }
    return ties;
  }
  //--------------------------------------------------------------
  //! vector of indices and weights for each parameter to be restrained across tiles
  std::pair<std::vector<int>, std::vector<double> >
  CCDTile2::TiedParameters(const std::vector<double> tie_tile,
                 const int& idx0)
  {
    //  r,w relative to rad0
    // For parameter types r,w,A, tie across tiles
    std::pair<std::vector<int>, std::vector<double> > kindexwt;
    ties_ = tie_tile;

    int idx = idx0;
    double weight;

    // r
    if (tie_tile[0] > 0.0) {
      weight = 1./(tie_tile[0]*tie_tile[0]);
      kindexwt.second.push_back(weight);
      kindexwt.first.push_back(idx);
    }
    idx++;
    // w
    if (tie_tile[1] > 0.0) {
      weight = 1./(tie_tile[1]*tie_tile[1]);
      kindexwt.second.push_back(weight);
      kindexwt.first.push_back(idx);
    }
    idx++;
    // A0
    if (tie_tile[2] > 0.0) {
      weight = 1./(tie_tile[2]*tie_tile[2]);
      kindexwt.second.push_back(weight);
      kindexwt.first.push_back(idx);
    }
    return kindexwt;
  }
  //--------------------------------------------------------------
  //! lower bounds for ipar'th parameter
  double CCDTile2::LowerBound(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    if (ipar == 0) {
      return 0.1; // r
    } else if (ipar == 1) {
      return 0.1; //  w
    } else if (ipar == 2) { // A0
      return 0.002;
    } else if (ipar <= 4 ) {
      return -0.5;  // x0 or y0
    }
    return -0.8; // ABCD for A
  }
  //--------------------------------------------------------------
  //! upper bounds for ipar'th parameter
  double CCDTile2::UpperBound(const int& ipar) const
  {
    if (ipar == 0) {
      return 1.41; // r
    } else if (ipar == 1) {
      return 1.0; //  w
    } else if (ipar == 3) { // A0
      return 1.0;
    } else if (ipar <= 4 ) {
      return +0.5;  // x0 or y0
    }
    return +0.8; // ABCD for A
  }
  //--------------------------------------------------------------
  //! "large shift" for ipar'th parameter
  double CCDTile2::LargeShift(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    if (ipar == 0) {
      return 0.3; // r
    } else if (ipar == 1) {
      return 0.3; //  w
    } else if (ipar == 3) { // A0
      return 0.3;
    } else if (ipar <= 4 ) {
      return +0.05;  // x0 or y0
    }
    return +0.05; // ABCD for A
  }
  //--------------------------------------------------------------
  // Return scale & derivatives for tile coordinates Xt, Yt
  void CCDTile2::ScaleDeriv(const bool& Deriv,
                                 const double& Xt, const double& Yt,
                                 double& scale,
                                 std::vector<double>& dgdp) const
  {
    // Xt, Yt in pixels
    double x = (Xt-xc0)/rad0;
    double y = (Yt-yc0)/rad0;

    double d = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)); // relative to rad0
    double phi = atan2((y-y0), (x-x0));

    std::vector<double> dAdp;
    double fA = 0.0;
    if (!circularlysymmetric) {
      fA = Afs.ValueDerivatives(phi, dAdp); // ABCD and dA/dp vector
    }
    double A = A0*(1.0 + fA); // A


    double z = 2.0*(d-r-w)/w;

    scale = radfunc.value(z, A);

    //*    scale = 0.5 * A * erfz + 1.0 - A;
    //^
    //    std::cout <<"\nXt,x0,Yt,y0 "<<Xt<<" "<<x0<<" "<<Yt<<" "<<y0<<"\n";
    //    std::cout <<"d,z,erfz "<<d<<" "<<z<<" "<<erfz<<"\n"; //^-

    if (Deriv) {
      // Count observations in each corner, only if calculating derivatives
      if (d > dcrnmin) { // in a corner
        int i = (x-x0) < 0.0 ? 0 : 1; // 0 or 1 if left or right
        int j = (y-y0) < 0.0 ? 0 : 1; // 0 or 1 if bottom or top
        ncorners(i,j)++;
      }
      dgdp = radfunc.deriv(nparams, w);  // fills 1st 3 slots in dgdp
      /*
        dgdp.resize(nparams);
        dgdp[0] = (twooverrootpi * A / w) * exp(-z*z);  // dgdr = -dgdd
        dgdp[1] = dgdp[0] * (1. + 0.5*z);               // dgdw
        double dgdA = 0.5 * erfz - 1.0;         // dgd(Atotal)
        dgdp[2] = dgdA *(1.0 + fA);  // dg/d(Aconstant) = dg/dA dA/d(Aconstant)
      */
      // dg/dA0 = dg/dA dA/dA0 = dg/dA (1 + fA)
      dgdp[2] = dgdp[2] *(1.0 + fA); // dg/d(Aconstant) = dg/dA dA/d(Aconstant)
      // x0, y0
      if (d == 0.0) {
        dgdp[3] = - dgdp[0]; // eg if y=y0, d=x-x0, dd/dx0 = -1
        dgdp[4] = - dgdp[0]; //  dd/dy0 = -1
      } else {
        dgdp[3] = - dgdp[0] * (x0 - x) / d;  // dgdx0 = dg/dd * dd/dx0
        dgdp[4] = - dgdp[0] * (y0 - y) / d;  // dgdy0 = dg/dd * dd/dy0
      }
      if (!circularlysymmetric) {
        for (size_t i=0; i<dAdp.size(); i++) {
          dgdp[i+5] = dAdp[i] * dgdp[2] * A0;
        }
      }
    }
  }
  //--------------------------------------------------------------
  std::string CCDTile2::format() const
  //! format scale type for printing
  {
    std::string text = "Tile correction for CCD detector";
    return text;
  }
  //--------------------------------------------------------------
  std::string CCDTile2::formattype() const
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
    if (circularlysymmetric) {
      text += "   r,w,A0 and x0,y0 are refined, all circularly symmetric\n";
    } else {
      text += "   r,w,A and x0,y0 are refined\n";
      text += "   A is parameterised as a constant (A0) + a 4-parameter Fourier series\n";
    }
    return text;
  }
  //--------------------------------------------------------------
  std::vector<std::string> CCDTile2::formatparameters(const std::vector<double>& sds) const
  //! format parameters for printing
  // Parameter order: r,w,A0,x0,y0,A(fourier)
  {
    int nsds = sds.size(); // number of SDs given, if any
    if (nsds > 0 ) {
      ASSERT (nsds == nparams);
    }
    std::vector<std::string> s;
    std::string line;
    int width = 37; // width of window
    line = "";
    for (int i=0;i<width;++i) {line += "-";};
    s.push_back(line);
    // counts for top corners
    line = "| "+
      StringUtil::LeftString(StringUtil::Strip(StringUtil::itos(ncorners(0,0),6)),6)+
      StringUtil::PadString(" ",width-12-4)+
      StringUtil::RightString(StringUtil::Strip(StringUtil::itos(ncorners(1,0),6)),6)+
      +" |";
    s.push_back(line);
    line = "| "+ StringUtil::PadString(" ",width-4)+ " |"; // "blank" line
    s.push_back(line);
    // r, w
    if (nsds == 0) {
      line = "| "+
        StringUtil::CentreString("r="+StringUtil::ftos(r,6,2)+
                               ", w="+StringUtil::ftos(w,6,2), width-4)+" |";
    } else { // with sds
      line = "| "+
        StringUtil::CentreString("r="+StringUtil::valueSD(r,sds[0],0,6,2)+
                                 ", w="+StringUtil::valueSD(w,sds[1],0,6,2), width-4)+" |";
    }
    s.push_back(line);
    // A0 (parameter 2), 4 Fourier terms (parameters 5-8)
    if (nsds == 0) {
      line = StringUtil::Strip(StringUtil::ftos(A0,6,2));
    } else {
      line = StringUtil::Strip(StringUtil::valueSD(A0,sds[2],0,6,2));
    }
    if (!circularlysymmetric) {
      line += Afs.format();
    }
    line = "| "+ StringUtil::CentreString("A="+line, width-4) + " |";
    s.push_back(line);
    if (nsds > 0 && !circularlysymmetric) {
      line = "";
      for (int i=0;i<4;++i) { // parameters 5-8
        line += StringUtil::ftos(sds[i+5],6,2);
        if (i < 3) {line += ",";}
      }
      line = "| "+ StringUtil::CentreString("SDs("+StringUtil::Strip(line)+")", width-4) + " |";
      s.push_back(line);
    }

    if (nsds == 0) {
      line = StringUtil::Strip(StringUtil::ftos(x0,7,2))+
        ", "+StringUtil::Strip(StringUtil::ftos(y0,7,2));
    } else {
      line = StringUtil::Strip(StringUtil::valueSD(x0,sds[3],0,7,2))+
        ", "+StringUtil::Strip(StringUtil::valueSD(y0,sds[4],0,7,2));
    }
    line = "| "+
      StringUtil::CentreString
      ("Tile centre: "+line, width-4)+
      " |";
    s.push_back(line);
    line = "| "+ StringUtil::PadString(" ",width-4)+ " |"; // "blank" line
    s.push_back(line);
    // counts for bottom corners
    line = "| "+
      StringUtil::LeftString(StringUtil::Strip(StringUtil::itos(ncorners(0,1),6)),6)+
      StringUtil::PadString(" ",width-12-4)+
      StringUtil::RightString(StringUtil::Strip(StringUtil::itos(ncorners(1,1),6)),6)+
      +" |";
    s.push_back(line);
    line = "";
    for (int i=0;i<width;++i) {line += "-";};
    s.push_back(line);
    return s;
  }
  //--------------------------------------------------------------
  std::string CCDTile2::formatTies() const
  {
    std::string s =
      std::string("Detector parameters r,w,A0 will be TIED across the tiles,")+
      " with SDs ";
    s += StringUtil::ftos(ties_[0],7,3)+","+StringUtil::ftos(ties_[1],7,3)+","
      +StringUtil::ftos(ties_[2],7,4)+"\n";
    if (ties_.size() > 5 && ties_[5] > 0.0) {
      s += "  radius parameter r will be tied to target "+
        StringUtil::ftos(ties_[5],7,4)+" with SD"+
        StringUtil::ftos(ties_[0],7,4)+"\n";
    }
    if (ties_.size() > 6 && ties_[6] > 0.0) {
      s += "  width parameter w will be tied to target "+
        StringUtil::ftos(ties_[6],7,4)+" with SD"+
        StringUtil::ftos(ties_[1],7,4)+"\n";
    }
    s += std::string
      ("  amplitude parameter A0 will be tied to zero with SD")+
      StringUtil::ftos(ties_[2],7,4)+"\n";
    s += std::string
      ("  and tile centre positions (x0,y0) will be tied to the true centre with SD")+
      StringUtil::ftos(ties_[3],7,3)+"\n";
    if (!circularlysymmetric) {
      s += std::string(
        "   Fourier coefficients of variation of A will be tied to zero with SD")+
        StringUtil::ftos(ties_[4],7,3)+"\n";
    }
    return s;
  }
  //--------------------------------------------------------------
  // Format all information into a labelled save format for later restoration
  std::string CCDTile2::FormatSave() const
  {
    std::string dump = "CCDTile2 V1 {\n"; // with version number
    dump += "XYmax " + ftos(xmax) + " " + ftos(ymax) + "\n";
    dump += "Rad0 " + ftos(rad0) + "\n";
    if (circularlysymmetric) {
      dump += "Symmetric\n";
    } else {
      dump += "Nonsymmetric\n";
    }
    dump += "Nparams " + itos(nparams) + " " + itos(nparams_smooth) + "\n";
    dump += "XYc0 " + ftos(xc0) + " " + ftos(yc0) + "\n";
    dump += "XY0 " + ftos(x0) + " " + ftos(y0) + "\n";
    dump += "Parameters\n" + StringUtil::FormatSaveVector(Parameters());
    dump += "Ncorners " + itos(ncorners(0,0)) + " " + itos(ncorners(0,1))
      + " " + itos(ncorners(1,0)) + " " + itos(ncorners(1,1)) + "\n";
    return dump+"}\n";
  }
  //--------------------------------------------------------------
  void CCDTile2::Restore(Fileread& FR)
  // restore
  {
    FR.ReadTag("CCDTile2"); // fails if tag does not match
    if (FR.GetTag() != "V1") {  // version check
      clipper::Message::message(Message_fatal
        ("CCDTile2::Restore incompatible version in "+FR.Filename()));
    }
    FR.Skip();
    FR.ReadTag("XYmax"); xmax = FR.Double(); ymax = FR.Double();
    FR.ReadTag("Rad0"); rad0 = FR.Double();
    std::string symm = FR.GetTag();
    circularlysymmetric = true;
    if (symm == "Nonsymmetric") {
      circularlysymmetric = false;
    }
    FR.ReadTag("Nparams"); nparams = FR.Int(); nparams_smooth = FR.Int();
    FR.ReadTag("XYc0"); xc0 = FR.Double(); yc0 = FR.Double();
    FR.ReadTag("XY0"); x0 = FR.Double(); y0 = FR.Double();
    FR.ReadTag("Parameters");
    std::vector<double> parameters = FR.DoubleVec(nparams);
    StoreParameters(parameters);
    FR.ReadTag("Ncorners");
    ncorners(0,0) = FR.Int(); ncorners(0,1) = FR.Int();
    ncorners(1,0) = FR.Int(); ncorners(1,1) = FR.Int();
    if (!FR.CheckEnd()) {
      clipper::Message::message(Message_warn
        ("CCDTile2::Restore unexpected tag "+FR.Tag()));
    }
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  double RadialFunctionGompertzCDF::value(const double& z, const double& A)
  //! calculate the value of the function and store intermediates
  {
    // Gompertz distribution CDF
    z_ = z;
    A_ = A;
    expmz = exp(-z); // also = dg/dz
    // scale = A f(z) + 1 - A
    // f(z) = 1 - exp(g(z))
    // g(z) = -exp(-z)
    gz = -expmz;
    fz = 1.0 - exp(gz);
    double scale = A * fz + 1.0 - A;
    return scale;
  }
  //--------------------------------------------------------------
  //! derivatives ds/dp0, for r,w,A, must follow a call to value
  std::vector<double> RadialFunctionGompertzCDF::deriv(const int& nparams,
                                                       const double& w) const
  {
    ASSERT (A_ >= -999.);
    std::vector<double> dsdp(nparams); // nparams >= 3
    // df/dz = exp(gz)exp(-z); ds/dz = A df/dz
    double dsdz = - A_ * exp(gz) * expmz;
    dsdp[0] = - dsdz * 2.0/w;          // ds/dr = -ds/dd
    dsdp[1] = dsdp[0] *(1.0 + 0.5*z_);  // ds/dw
    dsdp[2] = fz - 1.0;                // ds/dA
    return dsdp;
  }
  //--------------------------------------------------------------
  double RadialFunctionErfc::value(const double& z, const double& A)
  //! calculate the value of the function and store intermediates
  {
    // erfc
    z_ = z;
    A_ = A;
    // scale = A f(z) + 1 - A
    // f(z) = 0.5*erfc(z)
    fz = 0.5*erfc(z);
    double scale = A * fz + 1.0 - A;
    return scale;
  }
  //--------------------------------------------------------------
  //! derivatives ds/dp0, for r,w,A, must follow a call to value
  std::vector<double> RadialFunctionErfc::deriv(const int& nparams,
                                                const double& w) const
  {
    std::vector<double> dsdp(nparams); // nparams >= 3
    // df/dz = -(1/sqrt(pi))exp(-z^2); ds/dz = A df/dz
    double twooverrootpi = 2.0/sqrt(clipper::Util::pi()); // 2/sqrt(pi)
    double dsdz = - A_ * twooverrootpi * exp(-z_*z_);
    dsdp[0] = - dsdz/w;          // ds/dr = -ds/dd
    dsdp[1] = dsdp[0] *(1.0 + 0.5*z_);  // ds/dw
    dsdp[2] = fz - 1.0;                // ds/dA
    return dsdp;
  }
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
  std::vector<std::string> FlatTile::formatparameters(const std::vector<double>& sds) const
  //! format parameters for printing
  {
    if (sds.size() == 0) {
      std::vector<std::string> s(1,"  Tile scale: "+StringUtil::ftos(scale,8,3));
      return s;
    } else {
      std::vector<std::string> s(1,"  Tile scale: "+
                                 StringUtil::valueSD(scale,sds[0],0,8,3));
      return s;
    }
  }
  //--------------------------------------------------------------
  // Format all information into a labelled save format for later restoration
  std::string FlatTile::FormatSave() const
  {
    std::string dump = "FlatTile V1 {\n"; // with version number
    dump += "Nparams " + itos(nparams) + "\n";
    dump += "Scale " + ftos(scale) + "\n";
    return dump+"}\n";
  }
  //--------------------------------------------------------------
  void FlatTile::Restore(Fileread& FR)
  // restore
  {
    FR.ReadTag("FlatTile"); // fails if tag does not match
    if (FR.GetTag() != "V1") {  // version check
      clipper::Message::message(Message_fatal
        ("FlatTile::Restore incompatible version in "+FR.Filename()));
    }
    FR.Skip();
    FR.ReadTag("Nparams"); nparams = FR.Int();
    FR.ReadTag("Scale"); scale = FR.Double();
    if (!FR.CheckEnd()) {
      clipper::Message::message(Message_warn
        ("FlatTile::Restore unexpected tag "+FR.Tag()));
    }
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
    //^
    //    std::cout << "TilePixel::init " << Xmax<<" "<<Ymax<<" "
    //        <<njx <<" "<<njy<<std::endl;
    //^-

    scalexy.resize(njx, njy, 1.0);  // set all scales to 1.0
  }
  //--------------------------------------------------------------
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
    return 0.0;
  }
  //--------------------------------------------------------------
  //! upper bounds for ipar'th parameter
  double TilePixel::UpperBound(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    return 5.0;
  }
  //--------------------------------------------------------------
  //! "large shift" for ipar'th parameter
  double TilePixel::LargeShift(const int& ipar) const
  {
    ASSERT (ipar >= 0 && ipar < nparams);
    return 0.2;
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
    text += "   Independent scale determined for each group of "+
      StringUtil::Strip(StringUtil::itos(ngpxlX, 5))+" x "+
      StringUtil::Strip(StringUtil::itos(ngpxlY, 5))+ " pixels\n";
    return text;
  }
  //--------------------------------------------------------------
  std::vector<std::string> TilePixel::formatparameters(const std::vector<double>& sds) const
  //! format parameters for printing
  {
    std::vector<std::string> s;
    return s;
  }
  //--------------------------------------------------------------
  // Format all information into a labelled save format for later restoration
  std::string TilePixel::FormatSave() const
  {
    std::string dump = "TilePixel V1 {\n"; // with version number
    dump += "XYmax " + ftos(xmax) + " " + ftos(ymax) + "\n";
    dump += "Nparams " + itos(nparams) + "\n";
    dump += "NgpxlXY " + itos(ngpxlX) + " " + itos(ngpxlY) + "\n";
    dump += "NjXY " + itos(njx) + " " + itos(njy) + "\n";
    dump += "Parameters\n" + StringUtil::FormatSaveVector(Parameters());
    return dump+"}\n";
  }
  //--------------------------------------------------------------
  void TilePixel::Restore(Fileread& FR)
  // restore
  {
    FR.ReadTag("TilePixel"); // fails if tag does not match
    if (FR.GetTag() != "V1") {  // version check
      clipper::Message::message(Message_fatal
        ("TilePixel::Restore incompatible version in "+FR.Filename()));
    }
    FR.Skip();
    FR.ReadTag("XYmax"); xmax = FR.Double(); ymax = FR.Double();
    FR.ReadTag("Nparams"); nparams = FR.Int();
    FR.ReadTag("NgpxlXY"); ngpxlX = FR.Int(); ngpxlY = FR.Int();
    FR.ReadTag("NjXY"); njx = FR.Int(); njy = FR.Int();
    FR.ReadTag("Parameters");
    std::vector<double> parameters = FR.DoubleVec(nparams);
    StoreParameters(parameters);
    if (!FR.CheckEnd()) {
      clipper::Message::message(Message_warn
        ("TilePixel::Restore unexpected tag "+FR.Tag()));
    }
  }
  //--------------------------------------------------------------
  //=======================================================================
  //! construct or initialise from constant value
  FourierSmooth::FourierSmooth(const double& flatlevel,
                               const bool& isconstant)
  {
    setIsConstant(isconstant);
    setLevel(flatlevel);
  }
  //-------------------------------------------------------------------------
  //! set either: true for 4 parameters; false no constant term E, 4 params
  void FourierSmooth::setIsConstant(const bool& isconstant)
  {
    nparams = 4;
    if (isconstant) {nparams = 5;} // constant term E present
    parameters.resize(nparams);
  }
  //------------------------------------------------------------------------
  void FourierSmooth::setParameters(const std::vector<double> params)
  {
    ASSERT (int(params.size()) == nparams);
    parameters = params;
  }
  //--------------------------------------------------------------------------
  //! set a constant level, ie set E, A=B=C=D=0
  void FourierSmooth::setLevel(const double& flatlevel)
  {
    parameters.assign(nparams,0.0);
    if (nparams == 5 ){parameters[4] = flatlevel;}
  }
  /*
//----------------------------------------------------------------------------
//! set parameters from 5 values, ideally at pi/5 + n pi/2
void FourierSmooth::determineParameters(const std::vector<double>& phivalues,
                                        const std::vector<double>& values)
{
  if (int(phivalues.size()) != 5) {
    clipper::Message::message(Message_fatal
        ("FourierSmooth::determineParameters: must have 5 values"));
  }
  ASSERT (phivalues.size() == values.size());
  // We have 5 observational equations, i=1,5, values v[i], angles p[i]
  //   v[i] = A cos(p[i]) + B sin(p[i]) + C cos(2p[i]) + D sin(2p[i])
  // so we solve the equations
  //  v = [P] parameters
  //    where the row [P]i. = A cos(p[i]) + B sin(p[i]) + C cos(2p[i]) + D sin(2p[i])
  clipper::Matrix<double> P(5,5);
  for (size_t j=0; j<phivalues.size(); j++) {  // build the matrix
    P(j, 0) = cos(phivalues[j]);
    P(j, 1) = sin(phivalues[j]);
    P(j, 2) = cos(2.0*phivalues[j]);
    P(j, 3) = sin(2.0*phivalues[j]);
    P(j, 4) = 1.0;
    //^
    std::cout << "P["<<j<<" : " << phivalues[j] <<" : "<<
      P(j, 0) <<" " << P(j, 1) <<" " << P(j, 2) <<" " << P(j, 3) <<"\n";
  }
  // and solve it
  parameters = P.solve(values);
  //^
  std::cout << "determineParameters: new values ";
  //^-
}
  */
//-----------------------------------------------------------------------
//! get value at angle phi
double FourierSmooth::Value(const double& phi) const
{
  double cp = cos(phi);
  double sp = sin(phi);
  double c2p = cos(2.0*phi);
  double s2p = sin(2.0*phi);
  double value = parameters[0] * cp + parameters[1] * sp +
    parameters[2] * c2p + parameters[3] * s2p;
  if (nparams == 5) value += parameters[4];
  return value;
}
//-----------------------------------------------------------------------
//! get value at angle phi and dvdp its derivatives wrt parameters
double FourierSmooth::ValueDerivatives(const double& phi,
                                       std::vector<double>& dvdp) const
{
  double cp = cos(phi);
  double sp = sin(phi);
  double c2p = cos(2.0*phi);
  double s2p = sin(2.0*phi);
  double value = parameters[0] * cp + parameters[1] * sp +
    parameters[2] * c2p + parameters[3] * s2p;
  if (nparams == 5) value += parameters[4];
  dvdp.resize(nparams);
  dvdp[0] = cp;  // dv/dA = cos(phi)
  dvdp[1] = sp;  //  etc
  dvdp[2] = c2p;
  dvdp[3] = s2p;
  if (nparams == 5) {dvdp[4] = 1.0;}
  return value;
}
//-----------------------------------------------------------------------
std::string FourierSmooth::dump() const
{
  std::string s = "FourierSmooth: ABCD[E] = ";
  for (size_t j=0; j<parameters.size(); j++) {
    s += " "+StringUtil::ftos(parameters[j], 9, 4);
  }
  return s+"\n";
}
//-----------------------------------------------------------------------
std::string FourierSmooth::format() const
{
  std::string s = "";
  size_t last = parameters.size();
  if (nparams == 5) {
    last = parameters.size()-1;
    s += StringUtil::ftos(parameters.back(), 5, 2);
  }
  s += " {";
  for (size_t j=0; j<last; j++) {
    if (j != 0) {s += ",";}
    s += StringUtil::ftos(parameters[j], 5, 2);
  }
  s += "}";
  return StringUtil::Strip(s);
}
//-----------------------------------------------------------------------


}
