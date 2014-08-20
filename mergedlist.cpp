// mergedlist.cpp
//
// Merged data for one or more datasets, in clipper classes

#include "mergedlist.hh"
#include "file_util.hh"
#include "string_util.hh"
#include "cellgroup.hh"
#include "mtz_utils.hh"
#include "scala_util.hh"

using clipper::Message;
using clipper::Message_fatal;

namespace scala {
  // ---------------------------------------------------------
  void MergedDatasetIntensities::init(const clipper::HKL_info& hkl_info_list,
                                      const clipper::Cell& ccell)
  {
    Imean.init(hkl_info_list, ccell);
    Ipm.init(hkl_info_list, ccell);
  }
// ---------------------------------------------------------
  MergedList::MergedList(const hkl_unmerge_list& hkl_list, const SDmodel& SDM,
                         const std::string& Title, const int& datasetindex)
  {
    init(hkl_list, SDM, Title, datasetindex);
  }
  // ---------------------------------------------------------
  void MergedList::init(const hkl_unmerge_list& hkl_list, const SDmodel& SDM,
                        const std::string& Title, const int& datasetindex)
  // extract merged (averaged) data from hkl_list (with SDs corrected by SDM)
  // and store by dataset
  // If datasetindex >=0, store only that dataset
  {
    // Add new type to clipper registry
    clipper::CCP4MTZ_type_registry::add_group( "J_sigJ_ano", "IANO" );
    dataset_index = datasetindex;
    ndatasets = hkl_list.num_datasets();
    if (dataset_index >= 0) {
      ASSERT (dataset_index < ndatasets);
      ndatasets = 1;
    } else if (dataset_index == -2) {
      ndatasets = 1;  // combine datasets
    }

    maxintensity = -1000.;
    MeanValue meanI;
    title = Title;

    historylines = hkl_list.getHistory();

    // First construct the hkl list for all unique reflections
    reflection this_refl;
    std::vector<clipper::HKL> hkls;
    hkl_list.rewind();
    while (hkl_list.next_reflection(this_refl) >= 0)  {
      hkls.push_back(this_refl.hkl().HKL());
    }
    // Impose lattice symmetry constraints on cell
    CCtbxSym::CellGroup CG(hkl_list.symmetry().GetSpaceGroup());
    scala::Scell HKLcell = CG.constrain(hkl_list.cell());  // constrained cell
    clipper::Cell ccell = HKLcell.ClipperCell();

    // Initialise clipper::HKL_info list
    //   // space group, cell, resolution
    hkl_info_list.init(hkl_list.symmetry().ClipperGroup(), ccell,
                       clipper::Resolution(hkl_list.ResRange().ResHigh()),
                       false);  // no generation
    hkl_info_list.add_hkl_list(hkls);   // add hkl list
    spg_status = hkl_list.MtzSym().spg_confidence; // status of space group

    if (dataset_index < 0) {
      datasets = hkl_list.AllDatasets();
    } else { // just selected dataset
      datasets.assign(1, hkl_list.dataset(dataset_index));
    }

    // Initialise data objects for each dataset
    nrefdts.assign(ndatasets, 0);
    resmaxdts.assign(ndatasets, 0.0);

    datasetdata.resize(ndatasets);
    std::string simean = "IMEAN,SIGIMEAN";
    std::string sipm   = "I(+), SIGI(+), I(-), SIGI(-)";
    int jdts;

    for (int idts=0;idts<ndatasets;++idts) {
      // idts is local index to dataset stored, which may be just one
      // Individual cells for each dataset
      HKLcell = CG.constrain(datasets[idts].cell());  // constrained cell
      clipper::Cell dcell = HKLcell.ClipperCell();
      datasetdata[idts].init(hkl_info_list, dcell);
      // construct clipper-style mtzpaths
      PxdName pxdname = datasets[idts].pxdname();  // consensus name, all crystals for this dataset
      clipper::String mtzpath = pxdname.xname()+"/"+pxdname.dname()+"/";
      datasetdata[idts].mtzpath = mtzpath+"["+simean+", "+sipm+"]";;
      datasetdata[idts].mtzpathImean = mtzpath+"["+simean+"]";
      datasetdata[idts].mtzpathIpm = mtzpath+"["+sipm+"]";
      datasetdata[idts].cset =
        clipper::MTZdataset(pxdname.dname(), datasets[idts].wavelength());
      datasetdata[idts].cxtl =
        clipper::MTZcrystal(pxdname.xname(), pxdname.pname(), dcell);
    }

    // Read the data
    SelectedObservations allobs;     // all I+ and I-
    SelectedObservations obsplus;    // just I+
    SelectedObservations obsminus;   // just I-

    clipper::xtype data[10];  // 10 in case

    hkl_list.rewind();
    while (hkl_list.next_reflection(this_refl) >= 0)  {
      bool Centric = hkl_list.symmetry().is_centric(this_refl.hkl());
      Rtype invrsq = this_refl.invresolsq();
      //  Apply current SD correction to reflection (all observations)
      SDM.CorrectReflection(this_refl);
      for (int idts=0;idts<ndatasets;++idts) {
        // jdts is global index
        if (dataset_index == -2)
          {jdts = -1;}  // all datasets read together
        else if (dataset_index < 0)
          {jdts = idts;} // all datasets read separately
        else
          {jdts = dataset_index;} // selected dataset
        // mean I
        allobs.init(this_refl, jdts, ALL); // all data for selected dataset
        if (allobs.Number() > 0) {
          nrefdts[idts]++;
          resmaxdts[idts] = Max(resmaxdts[idts], invrsq);
          IsigI avI = allobs.Average();
          data[0] = avI.I();
          data[1] = avI.sigI();
          datasetdata[idts].Imean.data_import(this_refl.hkl().HKL(), data);
          maxintensity = Max(maxintensity, avI.I());
          meanI.Add(avI.I());
          data[2] = 0.0;
          data[3] = 0.0;
          data[4] = 0.0;
          if (!Centric) {
            data[0] = 0.0;
            data[1] = 0.0;
            obsplus.init(this_refl, jdts, IPLUS); // I+ data for selected dataset
            if (obsplus.Number() > 0) {
              avI = obsplus.Average();
              data[0] = avI.I();  // I+
              data[1] = avI.sigI();
            }
            obsminus.init(this_refl, jdts, IMINUS); // I- data for selected dataset
            if (obsminus.Number() > 0) {
              avI = obsminus.Average();
              data[2] = avI.I();  // I-
              data[3] = avI.sigI();
            }
          } else { // centric, I+ = I- = <I>
            data[2] =avI.I();
            data[3] =avI.sigI();
          }
          datasetdata[idts].Ipm.data_import(this_refl.hkl().HKL(), data);
        }
      } // end loop datasets
    } // end loop reflections
    meanintensity = meanI.Mean();
  }
  // ---------------------------------------------------------
  int MergedList::WriteDatasetToMTZ(const std::string& outfilename,
                                    const int& datasetIndex) const
  // Write data for datasetIndex to MTZ file
  // Return number of reflections written
  {
    int idx = InternalDTSindex(datasetIndex); // allow for one or all datasets stored
    ASSERT (datasetIndex < ndatasets);
    if (nrefdts.at(idx) <= 0) {
      Message::message(Message_fatal("MergedList::WriteMTZ no data for dataset"+
                                     clipper::String(datasetIndex)));
    }
    clipper::CCP4MTZfile mtzout;
    mtzout.open_write(outfilename);

    mtzout.set_title(title);

    // Update history
    std::vector<clipper::String>  history =
      MtzIO::addToHistory(historylines);

    mtzout.set_history(history);

    mtzout.set_spacegroup_confidence(spg_status);

    mtzout.export_crystal(datasetdata[idx].cxtl,
                          datasetdata[idx].mtzpath);
    mtzout.export_dataset(datasetdata[idx].cset,
                          datasetdata[idx].mtzpath);

    mtzout.export_hkl_info(hkl_info_list);
    mtzout.export_hkl_data(datasetdata[idx].Imean,
                           datasetdata[idx].mtzpathImean);
    mtzout.export_hkl_data(datasetdata[idx].Ipm,
                           datasetdata[idx].mtzpathIpm); ///!!

    mtzout.close_write();
    return nrefdts[idx];
  }
  // ---------------------------------------------------------
  bool MergedList::CheckNullImean(const clipper::data32::I_sigI& MIsig) const
  // returns false if I or sigI are Nan or sig = 0
  {
    if (clipper::Util::is_nan(MIsig.I()) ||
        clipper::Util::is_nan(MIsig.sigI()) ||
        MIsig.sigI() == 0.0) {
      return false;
    }
    return true;
  }
  // ---------------------------------------------------------
  int MergedList::CheckNullIano(const clipper::data32::J_sigJ_ano& MIsig) const
  // returns 0 if OK, -1 if both missing, +1 if I+ missing +2 if I- missing
  {
    int stat = 0;
    if (clipper::Util::is_nan(MIsig.I_pl()) ||
        clipper::Util::is_nan(MIsig.sigI_pl()) ||
        MIsig.sigI_pl() == 0.0) {
      stat = +1;
    }
    if (clipper::Util::is_nan(MIsig.I_mi()) ||
        clipper::Util::is_nan(MIsig.sigI_mi()) ||
        MIsig.sigI_mi() == 0.0) {
      if (stat == 0) {stat = +2;}
      else if (stat == +1) {stat = -1;}
    }
    return stat;
  }
  // ---------------------------------------------------------
  int MergedList::WriteDatasetToSCA(const std::string& outfilename,
                                    const int& datasetIndex) const
  // Write data for datasetIndex to SCA file
  // Return number of reflections written
  {
    ASSERT (datasetIndex < ndatasets);
    int idx = InternalDTSindex(datasetIndex); // allow for one or all datasets stored
    if (nrefdts.at(idx) <= 0) {
      Message::message(Message_fatal("MergedList::WriteSCA no data for dataset"+
                                     clipper::String(datasetIndex)));
    }

    FILE* scafile = OpenFile(outfilename, true);
    if (scafile == NULL) {
      Message::message(Message_fatal("Can't write file "+outfilename));
    }

    fprintf(scafile, "    1\n -987\n");

    float scale = 999900.0;  // keep scaled intensity in format %8.1f
    if (maxintensity > scale) {
      scale = scale/maxintensity;
    } else {
      scale = 1.0;
    }

    // cell
    // Impose lattice symmetry constraints on cell
    CCtbxSym::CellGroup CG(hkl_info_list.spacegroup());
    scala::Scell HKLcell = CG.constrain(datasets.at(idx).cell());  // constrained cell
    std::vector<Dtype> scell = HKLcell.UnitCell();
    for (int i=0;i<6;++i) {
      fprintf(scafile, "%10.3f", scell[i]);
    }
    std::string sgname = StringUtil::Strip(hkl_info_list.spacegroup().symbol_hm());
    char HorR = 'H'; // default H setting
    if (sgname[0] == 'H' || sgname[0] == 'R') {
      if (RhombohedralAxes(scell)) { // true if not H
        HorR = 'R';
      }}

    sgname = SGnameHtoR(sgname, HorR);  // R -> H
    // Space group name
    fprintf(scafile, " %s\n", sgname.c_str());

    clipper::HKL_info::HKL_reference_index ih;
    for (ih = hkl_info_list.first();!ih.last(); ih.next()) { // loop reflections
      bool centric = hkl_info_list.spacegroup().hkl_class(ih.hkl()).centric();
      if (centric) {
        if (CheckNullImean(datasetdata[idx].Imean[ih])) {
          float I = scale * datasetdata[idx].Imean[ih].I();
          float sigI = scale * datasetdata[idx].Imean[ih].sigI();
          fprintf(scafile, "%4d%4d%4d%8.1f%8.1f\n",
                  ih.hkl().h(), ih.hkl().k(), ih.hkl().l(), I, sigI);
        }
      } else { // I+ and I-
        int stat = CheckNullIano(datasetdata[idx].Ipm[ih]);
        if (stat >= 0) {
          float Ip = 0.0;
          float sigIp = 0.0;
          float Im = 0.0;
          float sigIm = 0.0;
          if (stat == 0 || stat == +1) {
            Ip = scale * datasetdata[idx].Ipm[ih].I_pl();
            sigIp = scale * datasetdata[idx].Ipm[ih].sigI_pl();
          }
          if (stat == 0 || stat == +2) {
            Im = scale * datasetdata[idx].Ipm[ih].I_mi();
            sigIm = scale * datasetdata[idx].Ipm[ih].sigI_mi();
          }
          if (sigIp <= 0.0) sigIp = -1;
          if (sigIm <= 0.0) sigIm = -1;
          fprintf(scafile, "%4d%4d%4d%8.1f%8.1f%8.1f%8.1f\n",
                  ih.hkl().h(), ih.hkl().k(), ih.hkl().l(), Ip, sigIp, Im, sigIm);
        }
      }
    } // end loop reflections
    return nrefdts[idx];
  }
  // ---------------------------------------------------------
  float MergedList::InvResMax(const int& datasetIndex) const
  // max(1/d^2) for given dataset, = 0 if unset
  {
    int idx = InternalDTSindex(datasetIndex); // allow for one or all datasets stored
    return resmaxdts.at(idx);
  }
  // ---------------------------------------------------------
  int MergedList::InternalDTSindex(const int& datasetIndex) const
  // return internal index to dataset datasetIndex
  // If dataset_index >=0, then only this dataset has been stored, so return 0
  // If dataset_index == -2, then all datasets have been stored together, return 0
  // If dataset_index <0, then all datasets have been stored, so return datasetIndex
  {
    int idx = datasetIndex;
    if (dataset_index >= 0) {  // only one dataset stored
      ASSERT (dataset_index == datasetIndex); // check that it is this one
      idx = 0; // index to only dataset
    } else if (dataset_index == -2) {  // combined datasets stored
      ASSERT (datasetIndex == 0); // check that it is the only one
      idx = 0; // index to only dataset
    } else { // all datasets stored, check request is for a valid one
      ASSERT ((datasetIndex >= 0) && (datasetIndex < ndatasets));
    }
    return idx;
  }
  // ---------------------------------------------------------
  clipper::Spacegroup MergedList::spacegroup() const
  {
    return hkl_info_list.spacegroup();
  }
  // ---------------------------------------------------------
  clipper::Cell MergedList::Cell() const
  {
    return hkl_info_list.cell();
  }
  // ---------------------------------------------------------
  double MergedList::resHigh() const
  {
    return hkl_info_list.resolution().limit();
  }
  // ---------------------------------------------------------
  // Reset current reflection pointer to first reflection, for given dataset
  void MergedList::start(const int& datasetIndex) const
  {
    current_dataset_index = InternalDTSindex(datasetIndex);
    hkl_index = datasetdata[current_dataset_index].Imean.first();
    at_start = true;
  }
  // ---------------------------------------------------------
  // Next IsigI, returns false if end of list
  bool MergedList::next(IsigI& Is) const
  {
    if (!at_start) {
      // increment index
      hkl_index.next();
    }
    at_start = false;
    if (hkl_index.last()) return false;
    Is = datasetdata[current_dataset_index].Imean[hkl_index];
    return true;
  }
  // ---------------------------------------------------------
  // get hkl for current reflection
  Hkl  MergedList::hkl() const
  {
    return Hkl(hkl_index.hkl());
  }
  // ---------------------------------------------------------
  // get clipper hkl for current reflection
  clipper::HKL MergedList::HKL() const
  {
    return hkl_index.hkl();
  }
  //--------------------------------------------------------------
  // resolution of current reflection
  double MergedList::invresolsq() const
  {
    return hkl_index.invresolsq();
  }
  // ---------------------------------------------------------
  //! return reference to Imean data for given dataset
  const clipper::HKL_data<clipper::data32::I_sigI>&
  MergedList::ImeanForDataset(const int& datasetIndex) const
  {
    int idx = InternalDTSindex(datasetIndex); // allow for one or all datasets stored
    return datasetdata.at(idx).Imean;
  }
  // ---------------------------------------------------------
  //! return reference to Imean data for given dataset
  clipper::HKL_data<clipper::data32::I_sigI>&
  MergedList::ImeanForDataset(const int& datasetIndex)
  {
    int idx = InternalDTSindex(datasetIndex); // allow for one or all datasets stored
    return datasetdata.at(idx).Imean;
  }
  // ---------------------------------------------------------
}
