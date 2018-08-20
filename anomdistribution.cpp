//
// anomdistribution.cpp
//

#include "anomdistribution.hh"
#include "selectedobservations.hh"
#include "tablegraph.hh"
#include "jiffy.hh"
#include "string_util.hh"

using phaser_io::LOGFILE;
using phaser_io::LXML;
using phaser_io::itos;


namespace scala {
  AllAnomDistributions::AllAnomDistributions(const hkl_unmerge_list& hkl_list,
                                             const SDmodel& SDM,
                                             const all_controls& controls,
                                             const AnalyseAnom& analysanom,
                                             const ResoRange& ResRange)
  // Analyse distribution of anomalous differences to get estimate
  // of maximum likely values, for all datasets
  //
  {
    reflection this_refl;

    // Resolution ranges
    int nresbin =  ResRange.Nbins();
    resrange = ResRange;
    // Number of datasets
    ndatasets = hkl_list.num_datasets();
    bool correlAnom;

    // rms DelAnom by dataset & resolution
    anomdistributions.resize(ndatasets);
    // Set number of resolution bins
    for (int id=0;id<ndatasets;++id) {
      anomdistributions[id].init(nresbin,
                                 hkl_list.dataset(id).pxdname(),
                                 analysanom.RmsDelAnom().at(id));
    }
    //  datasets
    std::vector<Dataset> datasets = hkl_list.AllDatasets();
    pxdnames.resize(ndatasets);
    wavelengths.resize(ndatasets);
    dnames.resize(ndatasets);
    for (int id=0;id<ndatasets;id++) {
      pxdnames[id] = datasets[id].pxdname();
      wavelengths[id] = datasets[id].wavelength();
      dnames[id] = datasets[id].Dname();
    }

    basedataset = controls.datasetcontrol.BaseDataset();
    int ncorrel = ndatasets*(ndatasets - 1)/2;    // number of anomalous cc
    int ndispcc = (ndatasets-2)*(ndatasets - 1)/2; // number of dispersive cc
    if (ncorrel > 0) {
      cca.resize(ncorrel);  // CC between anomalous differences
      ccadtsindex.resize(ncorrel);
      for (int i=0;i<ncorrel;++i) {cca[i].resize(nresbin);}
      int k=0;
      for (int j=0;j<ndatasets-1;++j) {
        for (int i=j+1;i<ndatasets;++i) {
          // for each CC, store pair of dataset indices
          ccadtsindex[k++] = std::pair<int,int>(j,i);
        }}
      ASSERT (k == ncorrel);
      if (ndispcc > 0) {
        ccd.resize(ndispcc); // CC between dispersive differences
        ccddtsindex.resize(ndispcc);
        for (int i=0;i<ndispcc;++i) {ccd[i].resize(nresbin);}
        // Choose base dataset as the one with the shortest wavelength
        // unless specified
        if (basedataset < 0) {
          double wvl = 10000.;
          for (int id=0;id<ndatasets;id++) {
            if (wavelengths[id] < wvl) {
              wvl = wavelengths[id];
              basedataset = id;
            }
          }
          if (wvl < 0.001) basedataset = 0;
        }
        // Store dataset index pairs for each CC
        k = 0;
        for (int j=0;j<ndatasets-1;++j) {
          if (j != basedataset) {
            for (int i=j+1;i<ndatasets;++i) {
              if (i != basedataset) {
                ccddtsindex[k++] = std::pair<int,int>(j,i);
              }
            }
          }
        }
        ASSERT (k == ndispcc);
      }
    }  // end ncorrel > 0

    double danom;
    std::vector<double> danomdts(ndatasets); // DelAnom for each dataset
    std::vector<double> Imeandts(ndatasets); // <I> for each dataset, for dispersive values
    SelectedObservations obsall, obsplus, obsminus;
    int nacc = 0;

    while (hkl_list.next_reflection(this_refl) >= 0)  {
      Rtype invresolsq = this_refl.invresolsq();
      bool Centric = hkl_list.symmetry().is_centric(this_refl.hkl());
      // Resolution bin
      int mres = ResRange.bin(invresolsq);
      //  Apply current SD correction to reflection (all observations)
      SDM.CorrectReflection(this_refl);
      correlAnom = true;  // flag for at least 2 I+ & 2 I- observations

      // loop datasets
      nacc = 0;
      danomdts.assign(ndatasets,0.0);
      for (int id=0;id<ndatasets;id++) {
        obsall.init(this_refl, id, ALL);
        if (obsall.Number() > 0) {
          Imeandts[id] = obsall.Average().I();
          if (!Centric) {
            obsplus.init(this_refl, id, IPLUS);
            if (obsplus.Number() > 0) {
              IsigI Iplus = obsplus.Average();
              if (obsplus.Number() <= 1) {correlAnom = false;}
              obsminus.init(this_refl, id, IMINUS);
              if (obsminus.Number() > 0) {
                IsigI Iminus = obsminus.Average();
                if (obsminus.Number() <= 1) {correlAnom = false;}
                // DelAnom
                double sig = Iplus.sigI()*Iplus.sigI() +
                  Iminus.sigI()*Iminus.sigI();
                if (sig > 0.0) {
                  // Store delAnom, & count reflections used for
                  // half-dataset correlations (ie with n+ & n- > 1)
                  danom = (Iplus.I()-Iminus.I());
                  anomdistributions[id].Add(mres, danom, correlAnom, obsplus, obsminus);
                  danomdts[id] = danom;
                  nacc++;
                }
              }
            }
          }
        }
      } // end loop datasets
      if (nacc > 1) {
        // correlations across datasets
        AddCorrelations(danomdts, Imeandts, mres);
      }
    } // end loop reflections
  }
  // ------------------------------------------------------------
  // Store slopes of normal probability anomplot for each dataset into Anomdistribution
  void AllAnomDistributions::SetSlope(const std::vector<double>& slope)
  {
    ASSERT (slope.size() == anomdistributions.size());
    for (size_t id=0;id<slope.size();++id) {
      anomdistributions[id].SetSlope(slope[id]);
    }
  }
  // ------------------------------------------------------------
  bool AllAnomDistributions::IsAnomalous(const all_controls& controls) const
  // return true if it appears that any dataset has significant anomalous
  // At present, anomalous scattering is considered to be present if any one of
  // the following is true (defaults in brackets):
  //  1) Anomplot slope > anomslopethreshold (1.3)
  //  2) CCanom > anomCCthreshold (0.3) in more than anomNbinthreshold bins (2), or overall,
  //     or interdataset CCanom
  //  3) RCRanom > anomRCRthreshold (1.3) in more than anomNbinthreshold bins (2), or overall
  //
  // It should be possible to estimate probabilities, but this will do for now
  {
    bool isanomalous = false;

    for (size_t k=0; k<cca.size(); k++) {  // loop cross-correlation CCanom
      // Interdataset CCs
      int nccanom  = 0; // ... CCanom
      correl_coeff cc;
      for (size_t mres=0;mres<cca[k].size();++mres) {
        if (cca[k][mres].result().val >
            controls.anomalouscontrol.anomCCthreshold) {
          nccanom++;
        }
        cc += cca[k][mres];
      }    // resolution bin loop
      if (nccanom > controls.anomalouscontrol.anomNbinthreshold ||
          cc.result().val > controls.anomalouscontrol.anomCCthreshold) {
          isanomalous = true;   // overall value above threshold
      }
    }   // end loop cross terms

    for (int id=0;id<ndatasets;id++) { // loop datasets
      if (anomdistributions[id].Slope() >
          controls.anomalouscontrol.anomslopethreshold) {
        isanomalous = true;
      }
      // count resolution bins above threshold for ...
      int nccanom  = 0; // ... CCanom
      int nrcranom = 0; // ... RCRanom
      int nbin = anomdistributions[id].Halfdataset().NresBin();
      if (anomdistributions[id].Halfdataset().CCanom().result().val >
            controls.anomalouscontrol.anomCCthreshold) {
        isanomalous = true;   // overall value above threshold
      }
      if (anomdistributions[id].Halfdataset().RMScorrelRatio() >
          controls.anomalouscontrol.anomRCRthreshold) {
        isanomalous = true;   // overall value above threshold
      }

      for (int mres=0;mres<nbin;++mres) {
        if (anomdistributions[id].Halfdataset().CCanom(mres).result().val >
            controls.anomalouscontrol.anomCCthreshold) {
          nccanom++;
        }
        if (anomdistributions[id].Halfdataset().RMScorrelRatio(mres) >
            controls.anomalouscontrol.anomRCRthreshold) {
          nrcranom++;
        }
      }
      if (nccanom > controls.anomalouscontrol.anomNbinthreshold) {
        isanomalous = true;
      }
    } // end loop datasets
    return isanomalous;
  }
  // ------------------------------------------------------------
  void AllAnomDistributions::AddCorrelations
  (const std::vector<double>& danomdts, const std::vector<double>& Imeandts,
   const int& mres)
  // Add in to correlation sums, resolution bin mres
  {
    ASSERT (int(danomdts.size()) == ndatasets);
    // Anomalous differences
    int k = 0;
    for (int j=0;j<ndatasets-1;++j) {
      for (int i=j+1;i<ndatasets;++i) {
        if (danomdts[i]!=0.0 && danomdts[j]!=0.0) {
          cca[k][mres].add(danomdts[i], danomdts[j]);
        }
        k++;
      }}
    if (ndatasets > 2) { // no correlation if < 3 datasets
      k = 0;
      for (int j=0;j<ndatasets-1;++j) {
        if (j != basedataset) {
          double dj = Imeandts[j] - Imeandts[basedataset];
          for (int i=j+1;i<ndatasets;++i) {
            if (i != basedataset) {
              if (dj!=0.0 && (Imeandts[i] - Imeandts[basedataset])!=0.0) {
                ccd[k][mres].add(Imeandts[i] - Imeandts[basedataset], dj);
              }
              k++;
            }
          }
        }
      }
    }
  }
  // ------------------------------------------------------------
  void AllAnomDistributions::Print(phaser_io::Output& output) const
  // Print correlation tables
  {
    if (ndatasets < 2) return;
    output.logTab(0, LOGFILE,
  "\nCorrelation coefficients for anomalous & dispersive differences between different datasets");
    output.logTab(0, LOGFILE,
                  "==========================================================================================\n");
    output.logTab(0, LOGFILE,"\nDatasets and wavelengths:\n");
    for (int id=0;id<ndatasets;++id) {
      std::string s("    ");
      if (id == basedataset) s = "base";
      output.logTabPrintf(0, LOGFILE,"%4d %4s %8.5f %s\n",
                          id+1, s.c_str(), wavelengths[id],
                          pxdnames[id].format().c_str());
    }

    std::string title = "=== Correlation of Anomalous Differences between datasets";
    std::string graphtitle = "Anom CCs v resln -";
    for (int id=0;id<ndatasets;++id) {graphtitle += " "+dnames[id];}
    std::string cl1 = "1st dataset         ";
    std::string cl2 = "2nd dataset         ";
    std::vector<correl_coeff> allcc;  // totals
    FormatTable(title, graphtitle, cl1, cl2, cca,
                ccadtsindex, false, allcc, output);

    // Format cross-correlation table
    title = "\nOverall correlation of Anomalous Differences between datasets\n";
    title += "      (Numbers in brackets)\n\n";
    CrossCorrelation(title, "DatasetAnomalousCorrelation",allcc, false, output);

    if (ndatasets < 3) return;

    // Dispersive differences
    title = "Correlation of Dispersive Differences between datasets";
    graphtitle = "Dispersive CCs v resln -";
    for (int id=0;id<ndatasets;++id) {graphtitle += " "+dnames[id];}
    cl1 = "1st difference      ";
    cl2 = "2nd difference      ";
    FormatTable(title, graphtitle, cl1, cl2, ccd,
                ccddtsindex, true, allcc, output);
    // Format cross-correlation table
    title = "\nCorrelation between datasets of Dispersive Differences from base set\n";
    title += "      (Numbers in brackets)\n\n";
    CrossCorrelation(title, "DatasetDispersiveCorrelation", allcc, true, output);
  }
  // ------------------------------------------------------------
  void AllAnomDistributions::FormatTable(const std::string& title,
                                         const std::string& graphtitle,
                                         const std::string& ccl1,
                                         const std::string& ccl2,
                                         const std::vector<std::vector<correl_coeff> >& cc,
                                         const std::vector<std::pair<int,int> >& ccidx,
                                         const bool& diff,
                                         std::vector<correl_coeff>& allcc,
                                         phaser_io::Output& output) const
  // diff = true for dispersive differences
  // private
  {
    TableGraph table(title);
    std::string id = "Graph-";
    if (diff) {
      id += "DispersiveDifferences";
    } else {
      id += "AnomalousDifferences";
    }
    table.StoreID(id);

    int ng = cc.size();     // number of graphs

    TableGraphPlot graph(graphtitle);
    for (int id=0;id<ng;++id) { // loop graphs for each dataset pair
      graph.AddLine(TableGraphPlotline(2,id*2+4));  // default colour
    }
    graph.SetXaxis("", true);  // x axis is 1/d^2
    Range yrange(0.0, 1.0);   // Y from 0 to 1
    graph.SetYaxis("", true, yrange);
    table.AddGraph(graph);

    // Column labels, zero-field flags & format
    std::vector<std::string> collabels;
    std::vector<bool> Zero;
    std::string fmt = "%3d%8.4f%6.2f ";  // leader
    std::string fmtn  = "%8.3f%8d";      // each part
    collabels.push_back("N");         // 1
    collabels.push_back("1/d^2");     // 2
    collabels.push_back("Dmid");      // 3
    Zero.push_back(false);
    Zero.push_back(false);
    Zero.push_back(false);
    std::string cl1 = ccl1; // header lines
    std::string cl2 = ccl2; // header lines
    for (size_t i=0;i<cc.size();++i) { // loop anomalous CCs for dataset pairs
      int idx1 = ccidx[i].first;  // dataset indices for this CC
      int idx2 = ccidx[i].second;
      std::string cl;
      if (diff) {
        // dispersive difference
        std::string bd = itos(basedataset+1);
        cl = StringUtil::Strip(itos(idx1+1)+bd+"-"+itos(idx2+1)+bd);
        std::string dl1 = StringUtil::Strip((dnames[idx1]+"-"+dnames[basedataset]));
        std::string dl2 = StringUtil::Strip((dnames[idx2]+"-"+dnames[basedataset]));
        cl1 += StringUtil::CentreString(dl1,16);
        cl2 += StringUtil::CentreString(dl2,16);
      } else {
        // anomalous
        cl = StringUtil::Strip(itos(idx1+1)+"-"+itos(idx2+1));
        cl1 += StringUtil::CentreString(dnames[idx1],16);
        cl2 += StringUtil::CentreString(dnames[idx2],16);
      }
      collabels.push_back("CC"+cl);  // eg CC1-2  (CC)
      collabels.push_back("N"+cl);   // eg N1-2   (number)
      Zero.push_back(true);
      Zero.push_back(true);
      fmt += fmtn;
    } // end loop dataset pairs
    fmt += "\n";

    table.StoreColumnFields(collabels, Zero, fmt);

    allcc.assign(cc.size(), correl_coeff());  // totals
    for (int ir=0;ir<resrange.Nbins();++ir) {
      int nonzero = 0; // number of non-zero item
      for (size_t i=0;i<cc.size();++i) {
        if (cc[i][ir].result().count > 0) {nonzero++;}
      }
      if (nonzero > 0) {
        table.StartLine();
        table.AddToLine(ir+1);
        table.AddToLine(resrange.middle(ir));
        table.AddToLine(resrange.middleA(ir));
        for (size_t i=0;i<cc.size();++i) {
          table.AddToLine(cc[i][ir].result().val);
          table.AddToLine(cc[i][ir].result().count);
          allcc[i] += cc[i][ir];
        }
        table.GetLine();
      }
    }  // end loop res bins
    table.CloseTable();

    // Output table to log file and XML
    output.logTab(0,LOGFILE, "\n"+table.format());
    output.logTab(0,LXML,table.XMLformat());

    std::string line = "Overall           ";
    char buf[256];
    for (size_t i=0;i<cc.size();++i) {
      sprintf(buf, fmtn.c_str(),
              allcc[i].result().val, allcc[i].result().count);
      line += std::string(buf);
    }
    output.logTab(0,LOGFILE, line);
  }
  // ------------------------------------------------------------
  void AllAnomDistributions::CrossCorrelation(const std::string& title,
                                              const std::string& tableid,
                                              const std::vector<correl_coeff>& allcc,
                                              const bool& diff,
                                              phaser_io::Output& output) const
  // format CC of anomalous differences between datasets as table
  // private
  {
    std::string outstring = title;
    std::string line1 = "            ";
    std::string line2 = line1;
    std::vector<std::string> ldf;
    for (int id=0;id<ndatasets;++id) { // loop datasets
      if (diff) {
        if(id != basedataset) {
          ldf.push_back(StringUtil::Strip((dnames[id]+"-"+dnames[basedataset])));
        }
      } else {
        ldf.push_back(dnames[id]);
      }
    }
    int ncc = ldf.size();

    for (int id=1;id<ncc;++id) { // loop datasets or differences from 2nd
      line1 += StringUtil::CentreString(ldf.at(id),8);
      line2 += StringUtil::CentreString("*",8);
    }
    outstring += line1+"\n"+line2+"\n";
    int k=0;
    for (int j=0;j<ncc-1;++j) {  // loop lines
      line1 = StringUtil::LeftString(ldf[j], 10)+"*";
      line2 = "           ";
      for (int i=1;i<ncc;++i) {
        if (i < j+1) {
          line1 += "        ";
          line2 += "        ";
        } else {
          line1 += StringUtil::ftos(allcc[k].result().val, 8, 3);
          line2 += StringUtil::CentreString
            (StringUtil::Strip("("+itos(allcc[k].result().count)+")"), 8);
          k++;
        }
      }
      outstring += line1+"\n"+line2+"\n";
    }  // end loop lines
    output.logTab(0,LOGFILE, "\n"+outstring);
    // XML output
    std::vector<std::pair<double,int> > valCount(allcc.size());
    for (size_t k=0; k<allcc.size(); k++) {
      valCount[k] = std::pair<double,int>(allcc[k].result().val, allcc[k].result().count);
    }
    output.logTab(0,LXML,
                  StringUtil::FormatXMLcrossTable("crosstable",
                                                  tableid, ldf, "CC", valCount));
  }
  // ------------------------------------------------------------
  // ------------------------------------------------------------
  void AnomDistribution::init(const int& Nresbin,
                              const PxdName& Dataset_pxd,
                              std::vector<MeanSD>& RMSdelanom)
  // Store number of resolution bins & clear arrays
  {
    nresbin = Nresbin;
    rmsDelAnom.assign(nresbin, MeanSD());
    nDelAnom.assign(nresbin,0);
    halfdataset.init(Nresbin, Dataset_pxd);
    halfdataset.StoreRMS(RMSdelanom);  // RMS values for this dataset, resolution bins
  }
  // ------------------------------------------------------------
  void AnomDistribution::Add(const int& mres,
                             const double& delAnom, const bool& correlAnom,
                             SelectedObservations& obsplus,
                             SelectedObservations& obsminus)
  // Store delAnom, & count reflections used for
  // half-dataset correlations (ie with n+ & n- > 1, correlAnom true),
  // and the halfdataset correlations
  //  for resolution range mres
  {
    rmsDelAnom[mres].Add(delAnom);
    // count reflections used for half-dataset correlations
    if (correlAnom) {
      nDelAnom[mres]++;
      halfdataset.AddAnom(mres, obsplus, obsminus);
    }
  }
  // ------------------------------------------------------------
  std::string AnomDistribution::formatStatus(const AnomalousStatus::anomalousStatus& anomalousstatus,
                                             const ResolutionLimit& anomresolimitCC)
  // static
  // (1) anomalousstatus indicates:-
  //   (a) explicit user input of anomalous ON or OFF, ANOMALOUS_ON_*  or ANOMALOUS_OFF_*,
  //       or if not specified just  ANOMALOUS_*
  //   (f) whether an initial estimate (from QQ-plot or CCanom) suggests an anomalous signal
  //       _* == _FOUND or _ABSENT
  //
  // (2) anomresolimitCC is the resolution limit for a "good" anomalous signal
  //
  // These two flags may be inconsistent
  {
    bool reslimitvalid = anomresolimitCC.valid(); // false if unset
    bool found = AnomalousStatus::isFound(anomalousstatus);       // true if FOUND
    bool consistent = true;

    std::string s;
    if (anomalousstatus == AnomalousStatus::ANOMALOUS_ON_FOUND) {
      s = "Anomalous flag switched ON in input, strong anomalous signal found";
    } else if (anomalousstatus == AnomalousStatus::ANOMALOUS_OFF_FOUND) {
      s = std::string("WARNING WARNING\n")+
        "Anomalous flag switched OFF in input but there appears to be a significant anomalous signal";
    } else if (anomalousstatus == AnomalousStatus::ANOMALOUS_FOUND) {
      s = "There appears to be a significant anomalous signal so anomalous flag was switched ON";
    } else if (anomalousstatus == AnomalousStatus::ANOMALOUS_ON_ABSENT) {
      s = "Anomalous flag switched ON in input but the anomalous signal is weak";
    } else if (anomalousstatus == AnomalousStatus::ANOMALOUS_OFF_ABSENT) {
      s = "Anomalous flag switched OFF in input, anomalous signal is weak";
    } else if (anomalousstatus == AnomalousStatus::ANOMALOUS_ABSENT) {
      s = "The anomalous signal appears to be weak so anomalous flag was left OFF";
    }

    std::string s2 = "";
    if (reslimitvalid) {
      // check for consistency between alternative metrics
      if (found && (anomresolimitCC.Status() < 0)) {
        consistent = false;
        s2 = ", but no anomalous resolution limit could be determined";
      } else if (!found && (anomresolimitCC.Status() >= 0)) {
        consistent = false;
        s2 = ", but a resolution limit could still be determined:\n"+
          anomresolimitCC.formatbrief(true);
      }
      if (found and consistent) {
        s2 = "\n" + anomresolimitCC.formatbrief(true);
      }
    }
    s += s2;

    return s;
  }
  // ------------------------------------------------------------
} // namespace scala
