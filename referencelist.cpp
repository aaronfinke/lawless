// referencelist.cpp
//

#include "referencelist.hh"
#include "columnlabels.hh"
#include "refinetargets.hh"
#include "refinereferencescale.hh"
#include "string_util.hh"
#include "timer.hh"
#include "report_errors.hh"

using phaser_io::LOGFILE;

namespace scala {
  //------------------------------------------------------
  ReferenceList::ReferenceList(const std::string& filename,
                               const std::string& labI,
                               const std::string& labsigI,
                               const double& resoLimit,
                               const bool& verbose,
                               phaser_io::Output& output)
  {
    init(filename, labI, labsigI, resoLimit, verbose, output);
  }
  //------------------------------------------------------
  void ReferenceList::init(const std::string& filename,
                           const std::string& labI,
                           const std::string& labsigI,
                           const double& resoLimit,
                           const bool& verbose,
                           phaser_io::Output& output)
  {
    hklrefname = filename;
    labI_ = labI;
    labsigI_ = labsigI;
    resolimit = resoLimit;

    inputtype = "HKLREF";
    bulksolvent = false;

    hklmergelist.init(filename, verbose, output);
    MtzIO::column_labels column_list;
    column_list.addLabin(labI_, labsigI_);
    hklmergelist.read(column_list, resolimit, verbose, output);
    output.logTabPrintf(1,LOGFILE,
                     "Number of reflections read %8d\n",
                     hklmergelist.num_obs());
    referencescalemodel.init();
  }
  //------------------------------------------------------
  void ReferenceList::init(const std::string& xyzin,
                           const double& resoLimit,
                           const bool& verbose,
                           phaser_io::Output& output)
  // Initialise from coordinate list
  {
    hklrefname = xyzin;
    resolimit = resoLimit;
    inputtype = "XYZIN";
    bulksolvent = false;

    hklmergelist.CreateFromAtoms(hklrefname, "", Scell(),
                                 resolimit, verbose, output);
    output.logTabPrintf(1,LOGFILE,
                     "Number of reflections calculated %8d\n",
                     hklmergelist.num_obs());
    referencescalemodel.init();
  }
  //------------------------------------------------------
  // Return scaled I sigI for given hkl
  IsigI ReferenceList::Isig(const Hkl& h) const
  {
    double scale = referencescalemodel.scale(h);
    return hklmergelist.Isig(h).scaleIs(scale);
  }
  //------------------------------------------------------
  // Return unscaled I sigI for given hkl
  IsigI ReferenceList::Isig0(const Hkl& h) const
  {
    return hklmergelist.Isig(h);
  }
  //------------------------------------------------------
  // column labels used
  std::vector<std::string> ReferenceList::columnLabels() const
  {
    return hklmergelist.columnLabels();
  }
  //------------------------------------------------------
  void ReferenceList::recordScaleReference(phaser_io::Output& output) const
  {
    std::string s = "Scaling against reference file: "+hklrefname;
    std::string s1,slabels,source;
    bool Fs = false;
    if (inputtype == "HKLREF") {
      ASSERT  (hklmergelist.fromReflectionData());
      if (hklmergelist.Amplitudes()) {
        s1 += "  amplitudes F squared to intensities";
        source = "squared F";
        Fs = true;
      } else {
        s1 += "  intensities";
        source = "intensities";
      }
      std::vector<std::string> columnlabels = columnLabels();
      if (columnlabels.size() > 0) {
        s1 += ", columns: ";
        for (size_t k=0; k<columnlabels.size(); k++) {
          slabels += " " + columnlabels[k];
        }
      }
    } else if (inputtype == "XYZIN") {
      s1 += "  calculated Fc^2 from coordinates";
      source = "coordinates";
    }
    output.logTab(0,LOGFILE, s);
    output.logTab(0,LOGFILE, s1+slabels);

    // XML
    output.logTab(0,LXML,"<ScalingType>");
    output.logTab(0,LXML,
                  StringUtil::MakeXMLtag("SourceType", source));
    output.logTab(0,LXML,
                  StringUtil::MakeXMLtag("Filename", hklrefname));
    if (slabels != "") {
      output.logTab(0,LXML,
                    StringUtil::MakeXMLtag("ColumnLabels", slabels));
    }
    if (Fs) {
      ReportErrors::printWarning(
           "It is better to use intensities than squared Fs",
           "ReferenceWarning");
    }
    output.logTab(0,LXML,"</ScalingType>");
  }
  //------------------------------------------------------
  bool ReferenceList::checkCompatible(const hkl_unmerge_list& hkl_list,
                                      const double& toleranceratio)
  // toleranceratio = 1.0 for difference > maximum resolution,
  //    larger tolerance is more lax
  // set status = -1 if the two lists have different symmetry (point group),
  // or +1 if cell is too different
  // return false if status != 0
  // see also private function
  {
    return referencescalemodel.checkCompatible(hkl_list, hklmergelist, toleranceratio);
  }
  //------------------------------------------------------
  bool ReferenceList::scaleToObserved(const hkl_unmerge_list& hkl_list,
                                      const int& datasetindex,
                                      const SDmodel& SDM,
                                      const double& toleranceratio,
                                      phaser_io::Output& output)
  // Scale reference list to observed, ie generate ReferenceScaleModel
  //  hkl_list        the observed data
  //  datasetindex    = -1 for all datasets
  //  SDM             sd model
  //  toleranceratio  = 1.0 for cell difference > maximum resolution,
  //                  larger tolerance is more lax
  //
  // Return false if datasets are incompatible
  {
    RefineTargets::REFINETARGETTYPES targettype(RefineTargets::LNCOSH);
    std::string s = "\nScaling reference data "+
      inputtype+" to observed data";
    s += " with minimisation target type "+
      RefineTargets::format(targettype)+"\n";
    output.logTab(0,LOGFILE,s);

    // MergedList list from unmerged
    MergedList mergedobslist(hkl_list, SDM, "", datasetindex);

    // Initialise model
    referencescalemodel.init(mergedobslist,
                             hklmergelist, toleranceratio);
    if (!referencescalemodel.isOK()) {
      return false;  // fail
    }

    int nprocs = 1;
    Timer timer;
    RefineReferenceScale refinereferencescale(mergedobslist.ImeanForDataset(0),
                                              hklmergelist,
                                              referencescalemodel,
                                              mergedobslist.meanIntensity(),
                                              nprocs);

    refinereferencescale.setTargetType(targettype);  // ln cosh

    int Ncycles = 10;
    // default protocols
    phaser::protocolPtr cPtr(new phaser::ProtocolScale(Ncycles));

    phaser::Minimizer Min;
    //^    output.setVerbose(true, true);
    Min.run(refinereferencescale, cPtr, output);    // run minimiser
    referencescalemodel.fixUp();  // remove isotropic part of anisotropic B

    output.logTabPrintf(0,LOGFILE,"\nTime for refinement: %8.3f\n",
                        timer.Stop());
    output.logTabPrintf(0,LOGFILE,"Number of reflections used: %8d\n",
                        refinereferencescale.Nobservations());
    output.logTab(0,LOGFILE,
                  "\n"+referencescalemodel.format());
    //^    output.setVerbose(false, false);

    return true;
  }
  //------------------------------------------------------
  bool ReferenceList::SFcalcScaleToObserved(const std::string& xyzin,
                                            const hkl_unmerge_list& hkl_list,
                                            const int& datasetindex,
                                            const SDmodel& SDM,
                                            const double& toleranceratio,
                                            const bool& verbose,
                                            phaser_io::Output& output)
  // Make reference list from atoms, scale to observed with bulk solvent (CLipper),
  //  and generate ReferenceScaleModel
  //  xyzin         file name for coordinate file
  //  hkl_list        the observed data
  //  datasetindex    = -1 for all datasets
  //  SDM             sd model
  //  toleranceratio  = 1.0 for cell difference > maximum resolution,
  //                  larger tolerance is more lax
  //
  //
  // Return false if datasets are incompatible
  {
    hklrefname = xyzin;
    inputtype = "XYZIN";
    bulksolvent = true;
    // MergedList list from unmerged
    MergedList mergedobslist(hkl_list, SDM, "", datasetindex);

    hklmergelist.createFromAtomsBulk(xyzin, mergedobslist,
                                     verbose, output);
    output.logTabPrintf(0,LOGFILE,
                        "%9d reflections generated and scaled\n",
                        hklmergelist.num_obs());
    referencescalemodel.init();  // unit scales
    return true;
  }

} // namespace scala
