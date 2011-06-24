// mergedlist.cpp
//
// Merged data for one or more datasets, in clipper classes

#include "mergedlist.hh"
#include "file_util.hh"
#include "string_util.hh"

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
  MergedList::MergedList(const hkl_unmerge_list& hkl_list, const SDmodel& SDM)
  {
    init(hkl_list, SDM);
  }
  // ---------------------------------------------------------
  void MergedList::init(const hkl_unmerge_list& hkl_list, const SDmodel& SDM)
  // extract merged (averaged) data from hkl_list (with SDs corrected by SDM)
  // and store by dataset
  {
    // Add new type to clipper registry
    clipper::CCP4MTZ_type_registry::add_group( "J_sigJ_ano", "IANO" );
    ndatasets = hkl_list.num_datasets();
    maxintensity = -1000.;
    // First construct the hkl list for all unique reflections
    reflection this_refl;
    std::vector<clipper::HKL> hkls;
    hkl_list.rewind();
    while (hkl_list.next_reflection(this_refl) >= 0)  {
      hkls.push_back(this_refl.hkl().HKL());
    }
    clipper::Cell ccell = hkl_list.Cell().ClipperCell();

    // Initialise clipper::HKL_info list
    //   // space group, cell, resolution
    hkl_info_list.init(hkl_list.symmetry().ClipperGroup(), ccell,
		       clipper::Resolution(hkl_list.ResRange().ResHigh()),
		       false);  // no generation
    hkl_info_list.add_hkl_list(hkls);   // add hkl list

    xdatasets = hkl_list.AllXdatasets();

    // Initialise data objects for each dataset
    nrefdts.assign(ndatasets, 0);
    resmaxdts.assign(ndatasets, 0.0);

    datasetdata.resize(ndatasets);
    std::string simean = "IMEAN,SIGIMEAN";
    std::string sipm   = "I(+), SIGI(+), I(-), SIGI(-)";

    for (int idts=0;idts<ndatasets;++idts) {
      datasetdata[idts].init(hkl_info_list, ccell);
      // construct clipper-style mtzpaths
      PxdName pxdname = xdatasets[idts].pxdname();
      clipper::String mtzpath = pxdname.xname()+"/"+pxdname.dname()+"/";
      datasetdata[idts].mtzpath = mtzpath+"["+simean+", "+sipm+"]";;
      datasetdata[idts].mtzpathImean = mtzpath+"["+simean+"]";
      datasetdata[idts].mtzpathIpm = mtzpath+"["+sipm+"]";
      datasetdata[idts].cset =
	clipper::MTZdataset(pxdname.dname(), xdatasets[idts].wavelength());
      datasetdata[idts].cxtl =
	clipper::MTZcrystal(pxdname.xname(), pxdname.pname(), ccell);
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
	// mean I
      	allobs.init(this_refl, idts, ALL); // all data for selected dataset
	if (allobs.Number() > 0) {
	  nrefdts[idts]++;
	  resmaxdts[idts] = Max(resmaxdts[idts], invrsq);
	  IsigI avI = allobs.Average();
	  data[0] = avI.I();
	  data[1] = avI.sigI();
	  datasetdata[idts].Imean.data_import(this_refl.hkl().HKL(), data);
	  maxintensity = Max(maxintensity, avI.I());
	  data[2] = 0.0;
	  data[3] = 0.0;
	  data[4] = 0.0;
	  if (!Centric) {
	    data[0] = 0.0;
	    data[1] = 0.0;
	    obsplus.init(this_refl, idts, IPLUS); // I+ data for selected dataset
	    if (obsplus.Number() > 0) {
	      avI = obsplus.Average();
	      data[0] = avI.I();  // I+
	      data[1] = avI.sigI();
	    }
	    obsminus.init(this_refl, idts, IMINUS); // I- data for selected dataset
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
  }
  // ---------------------------------------------------------
  int MergedList::WriteDatasetToMTZ(const std::string& outfilename,
				    const int& datasetIndex) const
  // Write data for datasetIndex to MTZ file
  // Return number of refelcetions written
  {
    ASSERT (datasetIndex < ndatasets);
    if (nrefdts.at(datasetIndex) <= 0) {
      Message::message(Message_fatal("MergedList::WriteMTZ no data for dataset"+
				     clipper::String(datasetIndex)));
    }
    clipper::CCP4MTZfile mtzout;
    mtzout.open_write(outfilename);

    mtzout.set_title(title);
    std::vector<clipper::String> history(1,"from Aimless");
    mtzout.set_history(history);

    mtzout.export_crystal(datasetdata[datasetIndex].cxtl,
			  datasetdata[datasetIndex].mtzpath);
    mtzout.export_dataset(datasetdata[datasetIndex].cset,
			  datasetdata[datasetIndex].mtzpath);

    mtzout.export_hkl_info(hkl_info_list);
    mtzout.export_hkl_data(datasetdata[datasetIndex].Imean,
			   datasetdata[datasetIndex].mtzpathImean);
    mtzout.export_hkl_data(datasetdata[datasetIndex].Ipm,
			   datasetdata[datasetIndex].mtzpathIpm); ///!!

    mtzout.close_write();
    return nrefdts[datasetIndex];
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
    if (nrefdts.at(datasetIndex) <= 0) {
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
    std::vector<Dtype> scell = xdatasets.at(datasetIndex).cell().UnitCell();
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
	if (CheckNullImean(datasetdata[datasetIndex].Imean[ih])) {
	  float I = scale * datasetdata[datasetIndex].Imean[ih].I();
	  float sigI = scale * datasetdata[datasetIndex].Imean[ih].sigI();
	  fprintf(scafile, "%4d%4d%4d%8.1f%8.1f\n",
		  ih.hkl().h(), ih.hkl().k(), ih.hkl().l(), I, sigI);
	}
      } else { // I+ and I-
	int stat = CheckNullIano(datasetdata[datasetIndex].Ipm[ih]);
	if (stat >= 0) {
	  float Ip = 0.0;
	  float sigIp = 0.0;
	  float Im = 0.0;
	  float sigIm = 0.0;
	  if (stat == 0 || stat == +1) {
	    Ip = scale * datasetdata[datasetIndex].Ipm[ih].I_pl();
	    sigIp = scale * datasetdata[datasetIndex].Ipm[ih].sigI_pl();
	  }
	  if (stat == 0 || stat == +2) {
	    Im = scale * datasetdata[datasetIndex].Ipm[ih].I_mi();
	    sigIm = scale * datasetdata[datasetIndex].Ipm[ih].sigI_mi();
	  }
	  if (sigIp <= 0.0) sigIp = -1;
	  if (sigIm <= 0.0) sigIm = -1;
	  fprintf(scafile, "%4d%4d%4d%8.1f%8.1f%8.1f%8.1f\n",
		  ih.hkl().h(), ih.hkl().k(), ih.hkl().l(), Ip, sigIp, Im, sigIm);
	}
      }
    } // end loop reflections
    return nrefdts[datasetIndex];
  }
  // ---------------------------------------------------------
  float MergedList::InvResMax(const int& datasetIndex) const
  // max(1/d^2) for given dataset, = 0 if unset
  {
    ASSERT (datasetIndex >= 0);
    ASSERT (datasetIndex < ndatasets);
    return resmaxdts.at(datasetIndex);
  }
  // ---------------------------------------------------------
}
