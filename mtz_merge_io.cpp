// mtz_merge_io.cpp
//
// Phil Evans 2009
//
// MTZ io for merged files (read only)
//

#include "mtz_merge_io.hh"
#include "hkl_datatypes.hh"
#include "string_util.hh"
#include "report_errors.hh"

namespace MtzIO {
  //--------------------------------------------------------------
  std::string FixRhombohedralSymbolMess(const std::string& SGsymbol,
                                        const Scell& cell)
  {
    std::string Symbol = SGsymbol;
    int kR = Symbol.find("R");
    if (kR < int(Symbol.size())) {
      // Symbol contains "R"
      if (cell.AngleTest(90.,90.,120.)) {
        // Cell is hexagonal (angles 90,90,120), reset lattice symbol to "H"
        Symbol[kR] = 'H';
      }
    }
    return Symbol;
  }
  //--------------------------------------------------------------
  //! Constructor: does nothing
  MtzMrgFile::MtzMrgFile()
    : fileopen(false)
  {}
  //--------------------------------------------------------------
  MtzMrgFile::~MtzMrgFile()
  {}
  //--------------------------------------------------------------
  bool MtzMrgFile::open_read(const std::string& filename_in)
  // Open file for reading
  // returns false if fails
  {
    merged = true;
    bool valid = false;
    if (fileopen)
      ReportErrors::printFatalError("MtzMrgFile: open_read - File already open");
    if ( filename_in == "")
      ReportErrors::printFatalError("MtzMrgFile: open_read - no filename given");

    // store filename
    if (getenv(filename_in.c_str()) != NULL) {
      filenamein = std::string(getenv(filename_in.c_str()));
    } else {
      filenamein = filename_in;
    }

    // open file
    clipper::CCP4MTZfile mtzin;
    try {
      mtzin.open_read(filenamein);
    }
    catch (Message_fatal) {
      // Failed to open file, missing or incomplete
      return valid;
    }
    fileopen = true;

    mtzfile_resolution = mtzin.resolution();
    ResMax = mtzfile_resolution.limit();

    mcell = Scell(mtzin.cell());

    spacegroup.init(mtzin.spacegroup());
    spacegroupsymbol = spacegroup.symbol_xHM();
    spg_status = mtzin.spacegroup_confidence();

    // Clipper seems to return spacegroup R3 as "R3" even on hexagonal axes
    //  so change "R" to "H" if cell is hexagonal
    ///    spacegroupsymbol = FixRhombohedralSymbolMess
    ///      (mtzin.spacegroup().symbol_hm(), mcell);

    // Check for column label "M_ISYM" as marker for unmerged file
    merged = true;
    std::vector<clipper::String> ColLab = mtzin.column_labels();
    std::vector<std::vector<clipper::String> > LabelTypes(ColLab.size());
    ProcessLabels processlabels(ColLab, column_labels());
    merged = processlabels.merged();
    mtzin.close_read();
    valid = true;
    return valid;
  }
  //--------------------------------------------------------------
  FileRead MtzMrgFile::MakeHklList(const std::string& mtzname,
                                   file_select& file_sel,
                                   col_controls& column_selection,
                                   MtzIO::column_labels& column_list,
                                   const scala::PxdName& InputPxdName,
                                   const scala::Scell& cell,
                                   std::string& output,
                                   const int& verbose,
                                   hkl_unmerge_list& hkl_list)
  // Fill an unmerged hkl_list from a merged file
  // On entry:
  //  mtzname            name of MTZ file (or logical name)
  //  file_sel           flags for general selection
  //                     - dataset selection
  //                     - batch selection
  //                     - resolution limits
  //                     - detector coordinate rejection ranges
  //  column_list        list of column names wanted by the program
  //  InputPxdName       PXD name to override dataset information from
  //                     MTZ file
  //  output             output string for printing
  //  verbose            set verbosity level
  //                      = 0 silent, = +1 usual summary
  //                      >= +2 debug
  //
  // On exit:
  //  hkl_list  has been filled and closed
  //
  // Returns FileRead.Opened() true if the MTZ file has been successfully opened
  // Returns FileRead.Read() true if the MTZ file has been read, false if has failed
  {
    if (!merged) {  // file must be merged for this function
      ReportErrors::printFatalError("MtzMrgFile::MakeHklList: not a merged file");
    }
    if (!hkl_list.IsEmpty()) {
      ReportErrors::printFatalError("MtzMrgFile::MakeHklList: hkl_list is not empty");
    }

    clipper::CCP4MTZfile mtzin;
    try {
      mtzin.open_read(filenamein);
    }
    catch (Message_fatal) {
      // Failed to open file, missing or incomplete
      return FileRead(false, false, false, 0);
    }
    output = "";
    if (verbose > 0) {
      output += FormatOutput::logTabPrintf(0,
          "\nReflection list generated from merged file: %s\n",filenamein.c_str());
      output += FormatOutput::logTabPrintf(0,
                   "\nTitle: %s\n\n", mtzin.title().c_str());
      output += FormatOutput::logTabPrintf(0,
              "   Space group from HKLIN file : %s\n",
                   mtzin.spacegroup().symbol_hm().c_str());
      output += FormatOutput::logTabPrintf(0, "   Cell: ");
      for (int i=0;i<6;i++) output += FormatOutput::logTabPrintf(0,"%8.2f",
                                 mcell[i]);
      output += FormatOutput::logTab(0,"\n");
      output += FormatOutput::logTabPrintf(0,
                          "   Maximum resolution in file:  %8.2f\n",
                           mtzin.resolution().limit());
    }
    mtzfile_resolution = mtzin.resolution();
    // Read all data into Clipper objects
    clipper::ftype ResoLimit = file_sel.reslimits().ResHigh();
    clipper::HKL_info hkl_info_list;
    clipper::HKL_data<clipper::data32::I_sigI> IsigData;
    clipper::HKL_data<clipper::data32::I_sigI_ano> IsigDataAnom;

    bool anom;
    ClipperLabelList labelthings = ReadData(mtzin, ResoLimit, column_list,
                                            hkl_info_list, anom, IsigData, IsigDataAnom,
                                            mtzdataset, (verbose>0), output);
    // IsigDataAnom populated if anom = true, IsigData always populated
    bool NoSigI = labelthings.nosig;  // true if no sigma column
    // Now construct hkl_list
    std::string title = mtzin.title();
    int Nref = hkl_info_list.num_reflections();

    // Construct dataset (clipper doesn't give us a project)
    std::vector<scala::Dataset> DataSets;
    PxdName pxdname("", labelthings.xname, labelthings.dname);
    if (!InputPxdName.is_blank()) {
      pxdname = InputPxdName;
    }
    DataSets.push_back(scala::Dataset(scala::Xdataset(pxdname,
                      Scell(mtzin.cell()), mtzdataset.wavelength(), 1)));
    // One batch
    std::vector<Batch> Batches(1);
    Batches[0].PXDname() = pxdname;
    int setid = 1;
    Batches[0].DatasetID() = setid;
    Batches[0].SetWavelength(float(mtzdataset.wavelength()));

    hkl_list.init(title, Nref,
                  hkl_symmetry(spacegroup), all_controls(),
                  DataSets, Batches);

    hkl_list.SetSpaceGroupStatus(spg_status);

    std::vector<clipper::String> chistory = mtzin.history();
    std::vector<std::string> history;
    for (size_t i=0; i<chistory.size(); i++) {
      history.push_back(chistory[i]); // clipper:String to std::string
    }
    hkl_list.addHistory(history);

    int isym = 1;
    int batch = 1;
    Rtype Xdet = 0.0;  Rtype Ydet = 0.0;
    Rtype phi = 0.0; Rtype time = 0.0;
    Rtype fraction_calc = 1.0;
    Rtype width = 0.0;
    Rtype LP = 0.0;
    int Npart = 1;
    int Ipart = 1;
    ObservationFlag ObsFlag;

    // Store all reflections
    // For resolution range found in file
    Range InvResRange;

    // Set start
    clipper::HKL_info::HKL_reference_index hkl_index = IsigData.first();
    at_start = true;
    clipper::data32::I_sigI Isig;
    clipper::data32::I_sigI_ano IsigAnom;
    //IsigI Is;
    Rtype I, Ipr, Ip, Im;
    Rtype sigI = 1.0;   // Dummy sigma = 1
    Rtype sigIpr = 1.0;
    Rtype sigIp = 1.0;   // Dummy sigma = 1
    Rtype sigIm = 1.0;
    bool haveIp, haveIm;
    scala::Hkl hred, hredm;
    int isymm;
    int nI = 0;
    int nIp = 0;
    int nIm = 0;
    int nIpacen = 0;
    int nImacen = 0;

    while (next(hkl_index)) {  // increments index if not at_start
      scala::Hkl hkl(hkl_index.hkl());  // hkl of current reflection
      hred  = hkl_list.symmetry().put_in_asu(hkl, isym);
      hredm = hkl_list.symmetry().put_in_asu(-hkl, isymm);  // -hkl
      bool written = false;
      if (anom) {
        haveIp = false;
        haveIm = false;
        bool centric = spacegroup.hkl_class(hkl.HKL()).centric();
        IsigAnom = IsigDataAnom[hkl_index];
        //      std::cout << "makelist " << hkl.format() <<" "<<IsigDataAnom[hkl_index].I()
        //                << " " << IsigData[hkl_index].I()<<"\n"; //^^
        Ip = 0.0;
        Im = 0.0;
        sigIp = 0.0;
        sigIm = 0.0;
        if (!clipper::Util::is_null(IsigAnom.I_pl())) { // I+ not null
          haveIp = true;
          Ip = IsigAnom.I_pl();
          if (!NoSigI) {
            sigIp = IsigAnom.sigI_pl();
          } else {
            sigIp = 1.0;
          }
        }
        if (!clipper::Util::is_null(IsigAnom.I_mi())) { // I- not null
          haveIm = true;
          Im = IsigAnom.I_mi();
          if (!NoSigI) {
            sigIm = IsigAnom.sigI_mi();
          } else {
            sigIm = 1.0;
          }
        }
        // If centric, check that I+ and I- are the same
        if (centric) {
          sigI = 0.0;
          if (haveIp && haveIm) {
            // both present, always take the unweighted average, even if same
            I = 0.5*(Ip + Im);
            sigI = 0.5*(sigIp + sigIm);
          } else {
            // just one
            if (haveIp) {
              I = Ip;
              sigI = sigIp;
            } else if (haveIm) {
              I = Im;
              sigI = sigIm;
            }
          }
          Ip = I;
          sigIp = sigI;
          Im = I;
          sigIm = sigI;
        } // end centric
        // still anom, centric or acentric
        // Store these observations, but not if sigI <= 0
        if (sigIp > 0.0) {
          written = true;
          Ipr = Ip;
          sigIpr = sigIp;
          hkl_list.store_part(hred, isym, batch, Ip, sigIp, Ipr, sigIpr,
                              Xdet, Ydet, phi, time,
                              fraction_calc, width, LP,
                              Npart, Ipart, ObsFlag);
          nIp++;
          if (!centric) {nIpacen++;}
        }
        if (!centric && sigIm > 0.0) { // no I- for centric
          Ipr = Im;
          sigIpr = sigIm;
          hkl_list.store_part(hredm, isymm, batch, Im, sigIm, Ipr, sigIpr,
                              Xdet, Ydet, phi, time,
                              fraction_calc, width, LP,
                              Npart, Ipart, ObsFlag);
          nIm++;
          if (!centric) {nImacen++;}
        }
      } else { // no anomalous
        Isig = IsigData[hkl_index];
        if (!clipper::Util::is_null(Isig.I_pl())) { // I+ not null
          scala::Hkl hred = hkl_list.symmetry().put_in_asu(hkl, isym);
          I = Isig.I_pl();
          Ipr = I;
          if (!NoSigI) {
            sigI = Isig.sigI_pl();
            sigIpr = sigI;
          } else {
            sigI = 1.0;
            sigIpr = sigI;
          }
          // Store this observation, but not if sigI <= 0
          if (sigI > 0.0) {
            written = true;
            hkl_list.store_part(hred, isym, batch, I, sigI, Ipr, sigIpr,
                                Xdet, Ydet, phi, time,
                                fraction_calc, width, LP,
                                Npart, Ipart, ObsFlag);
            nI++;
          }
        }

      }
      if (written) {
        InvResRange.update( hkl_index.invresolsq());  //smin, smax
      }
    }
    // end data read
    bool sorted = true;
    hkl_list.close_part_list(ResoRange(InvResRange), sorted);

    // Set into column_list the actual column numbers for columns requested by program
    //    (in column_list)
    //^!    column_list.get_col_lookup(*this);
    //^!    column_select col_select(column_list);

    // Store data presence flags
    //^!    hkl_list.StoreDataFlags(col_select.DataFlags());
    // Store file name
    hkl_list.AppendFileName(filenamein);

    std::string msg = "";
    if (nI + nIp + nIm == 0) {
      // no data at all
      ReportErrors::printFatalError("No valid data recognised in file");
    } else if (anom) {
      // check that we have both I+ and I-
      if (nIp == 0 || nIm == 0) {
        if (nIp == 0) {
          msg = "No valid data in I+ column";
        } else if (nIm == 0) {
          msg = "No valid data in I- column";
        }
      } else if (nIpacen != nImacen) {
        // unequal numbers of acentric I+ and I-
        output += FormatOutput::logTab
          (1,"NB: some acentric I+ or I- values missing");
        output += FormatOutput::logTabPrintf
          (2,"Number of valid I+:%9d, I-%9d\n", nIpacen, nImacen);
        std::string s = "Some acentric I+ or I- values missing";
        ReportErrors::printText(s, "DataMissing", false);
        s = StringUtil::Strip(StringUtil::itos(nIpacen, 9));
        ReportErrors::printText(s, "N_Iplus", false);
        s = StringUtil::Strip(StringUtil::itos(nImacen, 9));
        ReportErrors::printText(s, "N_Iminus", false);
      }
    } else {
      if (nI == 0) {
        msg = "No valid data in I column";
      }

    }
    if (msg != "") {
      ReportErrors::printWarning(msg, "DataMissing", true);
    }
    return FileRead(true, true, true, 0);
  }
  //--------------------------------------------------------------
  ClipperLabelList MtzMrgFile::ReadData(clipper::CCP4MTZfile& mtzin,
                   const double& ResoLimit,
                   const MtzIO::column_labels& column_list,
                   clipper::HKL_info& hkl_info_list,
                   bool& anom,
                   clipper::HKL_data<clipper::data32::I_sigI>& IsigData,
                   clipper::HKL_data<clipper::data32::I_sigI_ano>& IsigDataAnom,
                   clipper::MTZdataset& mtzdataset,
                   const bool& verbose,
                   std::string& output)
  // Read all selected data from MTZ file into clipper objects
  // hkl_info_list,  mtzdataset,
  // data: always fill IsigData (no anomalous), set missing sigmas to 1.0
  //    if anomalous data present in file, also fill IsoDataAnom
  // sets setanom = true if anomalous present
  // Returns label list
  {
    ResMax = mtzin.resolution().limit();
    if (ResoLimit > 0.0)
      ResMax = Max(ResoLimit, ResMax);

    // Set hkl list to desired resolution
    // reflections outside limits will be discarded
    hkl_info_list =
      clipper::HKL_info(spacegroup, mtzin.cell(),
                        clipper::Resolution(ResMax));

    if (verbose) {
      if (ResoLimit > 0.0 && ResMax > mtzfile_resolution.limit()+0.001) {
        output += FormatOutput::logTabPrintf(0, "Maximum resolution in file %s: %8.3f",
                                             filenamein.c_str(), mtzfile_resolution.limit());
        output += FormatOutput::logTabPrintf(0,"  restricted to %8.3f", ResMax);
        output += "\n";
      }
    }

    ProcessLabels processlabels(mtzin.column_labels(), column_list);
    IorF = processlabels.IorF();
    anom =  processlabels.anom();

    ClipperLabelList labelthings = processlabels.clipperlabellist();
    columnlabelsused = labelthings.labels;

    bool NoSigI = labelthings.nosig;  // true if there is no sigI column

    if (verbose) {
      if (IorF) {
        // column found is F
        output += FormatOutput::logTab(1, "Columns for amplitudes F (squared to I): "+
                                       labelthings.formatlabels() +"\n");
      } else {
        // column found is I
        output += FormatOutput::logTab(1, "Columns for intensities I: "+
                                       labelthings.formatlabels()+"\n");
      }
    }

    // Clipper MTZ dataset
    mtzin.import_dataset(mtzdataset, labelthings.path);

    //......................................................
    // Read header info & hkl list
    mtzin.import_hkl_info(hkl_info_list, false);
    IsigData.init(hkl_info_list, hkl_info_list.cell());
    if (anom) { // initialise if anomalous data is present
      IsigDataAnom.init(hkl_info_list, hkl_info_list.cell());
    }
    clipper::HKL_info::HKL_reference_index ih;
    if (IorF) {
      // File contains F, square all values
      // I = F^2
      // sigI = 2 F sigF  + sigF^2
      // temporary for F
      clipper::HKL_data_base* FsigData = NULL;  // pointer to the data, with
      if (anom) {
        FsigData = new clipper::HKL_data<clipper::data32::F_sigF_ano> (hkl_info_list);
      } else {
        FsigData = new clipper::HKL_data<clipper::data32::F_sigF> (hkl_info_list);
      }
      // Read in data
      mtzin.import_hkl_data(*FsigData, labelthings.path);
      mtzin.close_read();

      clipper::data32::I_sigI Isig;
      clipper::data32::I_sigI_ano IsigAnom;
      double F;
      const double iscale = 0.1;  // scale down F^2, by iscale^2

      for (ih = hkl_info_list.first(); !ih.last(); ih.next()) {
        Isig.set_null();
        IsigAnom.set_null();
        double sigF = 1.0;
        double sigI = 1.0/(iscale*iscale);
        bool OK = false;
        // OK if  1. NoSigI && F OK, or  2. F OK
        if (anom) {
          // F+, sigF+, F-, sigF-
          clipper::data32::F_sigF_ano fsig =
            (*dynamic_cast<clipper::HKL_data<clipper::data32::F_sigF_ano>*>(FsigData))[ih];
          if (!clipper::Util::is_null(fsig.f_pl())) { // F not null
            F = fsig.f_pl();
            IsigAnom.I_pl() = F * F;
            if (!NoSigI) {
              if (!clipper::Util::is_null(fsig.sigf_pl())) { // F not null
                sigF = fsig.sigf_pl();
                sigI = 2.*F*sigF + sigF*sigF;
              }
            }
            IsigAnom.sigI_pl() = sigI;
            OK = true;
          }
          // set I-
          if (!clipper::Util::is_null(fsig.f_mi())) { // F not null
            F = fsig.f_mi();
            IsigAnom.I_mi() = F * F;
            if (!NoSigI) {
              if (!clipper::Util::is_null(fsig.sigf_mi())) { // F not null
                sigF = fsig.sigf_mi();
                sigI = 2.*F*sigF + sigF*sigF;
              }
            }
            IsigAnom.sigI_mi() = sigI;
            OK = true;
          }
          if (OK) {
            IsigAnom.scale(iscale);
            IsigDataAnom[ih] = IsigAnom;
            // set IsigData as well
            IsigData[ih] = clipper::data32::I_sigI(IsigAnom.I(), IsigAnom.sigI());
          }
        } else { // no anomalous
          clipper::data32::F_sigF fsig =
            (*dynamic_cast<clipper::HKL_data<clipper::data32::F_sigF>*>(FsigData))[ih];
          if (!clipper::Util::is_null(fsig.f())) { // F not null
            F = fsig.f();
            Isig.I() = F * F;
            if (!NoSigI) {
              if (!clipper::Util::is_null(fsig.sigf())) { // F not null
                sigF = fsig.sigf();
                sigI = 2.*F*sigF + sigF*sigF;
              }
            }
            Isig.sigI() = sigI;
            Isig.scale(iscale);
            IsigData[ih] = Isig;
          }
        }
      }
    } else {
      // Read in data
      if (anom) {
        mtzin.import_hkl_data(IsigDataAnom, labelthings.path);
        mtzin.close_read();
        // Populate IsigData array as well
        for (ih = IsigDataAnom.first(); !ih.last(); ih.next()) {
          if (!clipper::Util::is_null(IsigDataAnom[ih].I())) { // I not null
            IsigData[ih] = clipper::data32::I_sigI(IsigDataAnom[ih].I(), IsigDataAnom[ih].sigI());
          }
        }
      } else {
        mtzin.import_hkl_data(IsigData, labelthings.path);
        mtzin.close_read();
      }
    }
    return labelthings;
  }
  //--------------------------------------------------------------
//--------------------------------------------------------------
  //--------------------------------------------------------------
  // Next IsigI, returns false if end of list
  bool MtzMrgFile::next(clipper::HKL_info::HKL_reference_index& hkl_index)
  {
    if (!at_start) {
      // increment index
      hkl_index.next();
    }
    at_start = false;
    if (hkl_index.last()) return false;
    return true;
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
} // namespace MtzIO
