// runcorrelations.cpp
//
// Matrix of correlations between runs
//
// Could do correlations on I or E^2, but there is no significant difference if
// they are done in resolution bins, then the bin CCs averaged to get an overall
// value
//

#include <vector>

#include "runcorrelations.hh"
#include "scala_util.hh"
#include "Output.hh"
#include "string_util.hh"

#include "tablegraph.hh"

using phaser_io::LOGFILE;
using phaser_io::LXML;


namespace scala {

  // ------------------------------------------------------------
  RunCorrelations::RunCorrelations(const hkl_unmerge_list& hkl_list,
                                   const SDmodel& SDM,
                                   const Normalise& NormRes,
                                   const int& nresbins)
  {
    init(hkl_list, SDM, NormRes, nresbins);
  }
  // ------------------------------------------------------------
  void RunCorrelations::init(const hkl_unmerge_list& hkl_list,
                             const SDmodel& SDM,
                             const Normalise& NormRes,
                             const int& nresbins)
  // Accumulate CC on I and E^2 for each run pair in resolution bins,
  // then average over resolution bins. This gives a better average
  // over resolution than calculating CC over all data
  {
    nruns = hkl_list.num_runs();
    if (nruns < 2) {return;}
    hkl_list_p = &hkl_list;
    normres = &NormRes;
    int nmatrix = (nruns * (nruns-1)) / 2; // unique pairs
    mnCC_E2.assign(nmatrix, 0.0);
    nmeanCC.assign(nmatrix, 0);

    std::vector<Run> runlist = hkl_list_p->RunList();
    runnumbers.resize(nruns);
    // Get actual run numbers for each run index
    for (int irun=0;irun<nruns;++irun) {
      runnumbers[irun] = runlist[irun].RunNumber();
    }

    resrange = hkl_list.ResRange();
    resrange.SetNbins(nresbins);

    // for each run pair, for each resolution bin
    //    std::vector<std::vector<correl_coeff> > CC_E2_res;  // on E^2

    CC_E2_res.resize(nmatrix);
    for (int k=0; k<nmatrix; k++) {
      CC_E2_res[k].resize(nresbins);
    }

    const bool VARIANCEWEIGHT = true; // FIXME

    observation this_obs;
    int iobs;

    hkl_list.rewind();
    reflection this_refl;
    while (hkl_list.next_reflection(this_refl) >= 0)  {
      SDM.CorrectReflection(this_refl);

      std::vector<MeanVariance> meanIrun(nruns);

      // loop accepted observations to get weighted mean I by run
      while ((iobs = this_refl.next_observation(this_obs)) >= 0) {
        int irun = this_obs.run();
        float weight = 1.0/(this_obs.ksigI()*this_obs.ksigI());
        meanIrun[irun].Add(this_obs.kI(), weight);
      } // end loop observations

      // Normalisation factor (multiplying)
      float normscale = normres->CorrAvg(this_refl.invresolsq());
      // resolution bin
      int mres = resrange.bin(this_refl.invresolsq());

      int k = 0;
      for (int j=0;j<nruns-1;++j) {
        if (meanIrun[j].Count() > 0) {
          for (int i=j+1;i<nruns;++i) {
            if (meanIrun[i].Count() > 0) {
              double w = 1.0;
              if (VARIANCEWEIGHT) {
                double mv = normscale*normscale*(meanIrun[i].VarianceFromWeights() +
                                                 meanIrun[j].VarianceFromWeights());
                if (mv <= 0.0) {
                  w = 0.0;
                } else {
                  w = 1.0/mv;
                }
              }
              if (w > 0.0) {
                k = j*nruns - j*(j+1)/2 + i - j - 1;
                //^std::cout << i<<" "<<j<<" "<<k<<" i,j,k refl\n";
                CC_E2_res[k][mres].add(normscale*meanIrun[i].Mean(),
                                       normscale*meanIrun[j].Mean(), w);
              }
            }
          } // i
        }
      } // j
    }  // end loop reflections

    // Average over resolution bins
    resweighttype = +2;  // 0 = unit, +1 Var(CC), +2 Number
    int k = 0;
    for (int j=0;j<nruns-1;++j) {
      for (int i=j+1;i<nruns;++i) {
        MeanValue m2;
        int n = 0;
        for (int l=0;l<nresbins;++l) {
          double w = 1.0;
          if (resweighttype == +1) {
            double sdCC = CC_E2_res[k][l].SD();
            w = 0.0;
            if (sdCC > 0.0) {
              w = 1.0/(sdCC*sdCC);
            }
          } else if (resweighttype == +2) {
            w = CC_E2_res[k][l].Number();
          }
          if (w > 0.0) {
            m2.Add(CC_E2_res[k][l].CC(), w);
            n += CC_E2_res[k][l].Number();
          }
        }
        //^std::cout << i<<" "<<j<<" "<<k<<" i,j,k adding\n";
        mnCC_E2[k] = m2.Mean();
        nmeanCC[k] += n;
        k++;
      } // i
    } // j

    //^
    //    k = 0;
    //    for (int i=0;i<nresbins;++i) {
    //      double CC = CC_E2_res[k][i].CC();
    //      double sdCC = CC_E2_res[k][i].SD();
    //      std::cout << "Rbin, CC, sdCC " <<i<<" "<<CC<<" "<<sdCC<<"\n";
    //    } //^-
  }
  // ------------------------------------------------------------
  void RunCorrelations::formatTable(phaser_io::Output& output) const
  {
    if (nruns < 2) {return;}

    output.logTab(0, LOGFILE,
                  "\nMatrix of correlations of E^2 between runs");
    output.logTab(0, LOGFILE,
                    "==========================================");

    output.logTab(0, LOGFILE,
                  std::string("\nWeighted correlation coefficients are calculated for each")+
                  " run pair in resolution ranges,\n  then averaged over bins\n\n");

    std::string s;
    if (resweighttype == 0) {
      s = "Bins in averages over resolution are equally weighted";
    } else if (resweighttype == +1) {
      s = "Bins in averages over resolution are weighted by 1/Var(CC)";
    } else if (resweighttype == +2) {
      s = "Bins in averages over resolution are weighted by number of reflections";
    }
    if (s != "") {
      output.logTab(0, LOGFILE, s);
      output.logTab(0, LOGFILE, "\n");
    }

    // List run information
    std::vector<Run> runlist = hkl_list_p->RunList();
    ASSERT (int(runlist.size()) == nruns);
    for (int irun=0;irun<nruns;++irun) {
      std::string sr = "Run "+StringUtil::itos(runnumbers[irun]);
      sr += " maximum resolution "+
        StringUtil::ftos(runlist[irun].GetResoRange().ResHigh(), 7, 3);
      sr += " Dataset: "+
        hkl_list_p->dataset(runlist[irun].DatasetIndex()).pxdname().format();
      output.logTab(0, LOGFILE, sr);
    }

    const int MAXCOLS = 12;
    const int MAXROWS = 20;
    bool clipped = false;

    std::string line1 = "\n         Run  ";
    for (int i=1;i<nruns;++i) {
      if (i >= MAXCOLS) {
        line1 += " ...";
        break;
      } else {
        line1 += StringUtil::itos(runnumbers[i], 7);
      }
    }
    output.logTab(0, LOGFILE, line1);

    int k = 0;
    for (int j=0;j<nruns-1;++j) {
      line1 = "Run"+StringUtil::itos(runnumbers[j], 4) + " CC(E^2) ";
      std::string line2 = "          N     ";

      for (int i=1;i<nruns;++i) {
        if (i <= j) {
          line1 += "       ";
          line2 += "       ";
        } else {
          if (i >= MAXCOLS) {
            if (!clipped) {line1 += " ...";}
            clipped = true;
          } else {
            //^std::cout << i<<" "<<j<<" "<<k<<" i,j,k\n";
            line1 += StringUtil::ftos(mnCC_E2[k], 7, 3);
            line2 += StringUtil::itos(nmeanCC[k], 7);
          }
          k++;
        }
      }
      if (j <= MAXROWS) {
        output.logTab(0, LOGFILE, line1);
        output.logTab(0, LOGFILE, line2);
      }
    }

    if (clipped) {
      output.logTab(0, LOGFILE," ...");
      output.logTab(0, LOGFILE,
            std::string("\nSome matrix entries have been suppressed to save paper:")+
                    " all entries are written to XML");
    }


    // Write all CCs (by resolution) to XML
    writeCCtoXML(output);

    const int MAXGRAPHRUNS = 8;
    int nrungraph = Min(nruns, MAXGRAPHRUNS); // not too many runs to graph
    resolutiongraph(nrungraph, output);
  }
  // ------------------------------------------------------------
  void RunCorrelations::resolutiongraph(const int& nrungraph,
                                        phaser_io::Output& output) const
  // graph pairs up to run index nrungraph
  {
    std::string s = "Run pair correlations by resolution";
    TableGraph table(" ==== "+s);
    table.StoreID("Graph-RunCorrelationsVsResolution");

    std::vector<std::string> labels;
    int k = 0;
    for (int j=0;j<nruns-1;++j) {
      for (int i=j+1;i<nruns;++i) {
        std::string s = "CC"+StringUtil::itos(runnumbers[i],3) + "-" +
          StringUtil::itos(runnumbers[j],3);
        labels.push_back(StringUtil::Strip(s));
        k++;
      }
    }
    int ncc = labels.size();

    Range xrange = resrange; // x axis range to full resolution limit
    xrange.first() = 0.0;    // from 0

    TableGraphPlot graph(s);
    std::string description =
      "Correlations between runs may point out bad runs; see also the full matrix table";
    graph.SetDescription(description);

    for (int k=0;k<ncc;++k) {
      graph.AddLine(TableGraphPlotline(2,k+4)); // column numbers for x,y
    }
    graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
    graph.SetYaxis("", true);  // y axis from 0 to maximum
    table.AddGraph(graph);

    // Define column labels
    std::vector<std::string> collabels;
    collabels.push_back("N");         // 1
    collabels.push_back("1/d^2");     // 2
    collabels.push_back("Dmid");      // 3
    for (int k=0;k<ncc;++k) {
      collabels.push_back(labels[k]);     // 4
    }
    int nc = collabels.size();
    // which columns should have "0.0" replaced by "-"
    bool z[] =
      {false, false, false};
    std::vector<bool> Zero(z, z+3);
    std::string fmt;
    for (int k=0;k<ncc;++k) {
      Zero.push_back(true);
      fmt+= "%7.3f"; // excluding 1st 3 columns
    }
    fmt += "\n";
    // store labels, zero flags and format
    table.StoreColumnFields(collabels, Zero, "%3d%8.4f%7.2f"+fmt);

    std::vector<MeanValue> mnCC_E2(ncc);

    int n=1;
    int nc3 = 3;
    for (int i=0;i<resrange.Nbins();++i) {
      bool write = false;
      std::vector<double> val(ncc, 0.0);
      for (int k=0;k<ncc;++k) {
        if (CC_E2_res[k][i].Number() > 0) {
          write = true;
          val[k] = CC_E2_res[k][i].CC();
          mnCC_E2[k].Add(CC_E2_res[k][i].CC());
        }
      }
      if (write) {
        // Store each table line
        table.Line(val, nc3, n++, resrange.middle(i),
                   resrange.middleA(i));
      }
    }
    table.CloseTable();
    output.logTab(0,LOGFILE, "\n"+table.format());
    output.logTab(0,LXML,table.XMLformat());

    s = "Overall:          ";
    for (int k=0;k<ncc;++k) {
      s += StringUtil::ftos(mnCC_E2[k].Mean(), 7, 3);
    }
    output.logTab(0,LOGFILE, s);
  }
  // ------------------------------------------------------------
  void RunCorrelations::writeCCtoXML(phaser_io::Output& output) const
  // output all CC by resolution to XML
  {
    // XML output
    // 1) overall values as a table for ccp4i2
    int k = 0;
    std::vector<std::string> names(nruns);
    std::vector<std::pair<double,int> > valCount(mnCC_E2.size());
    for (int i=0;i<nruns;++i) {
      names[i] = "Run"+StringUtil::Strip(StringUtil::itos(runnumbers[i]));
    }
    for (int j=0;j<nruns-1;++j) {
      for (int i=j+1;i<nruns;++i) {
        valCount[k] = std::pair<double,int>(mnCC_E2[k], nmeanCC[k]);
        k++;
      }
    }
    output.logTab(0,LXML,
                  StringUtil::FormatXMLcrossTable
                  ("crosstable", "RunCorrelation",
                   names, "CC", valCount));

    // 2) the full gory details, by resolution
    output.logTab(0, LXML, "<runcorrelations>");
    output.logTab(0, LXML,
                  StringUtil::MakeXMLtag("Nruns", StringUtil::itos(nruns,3)));

    output.logTab(0,LXML,
          StringUtil::MakeXMLtag("Nresolutionbins",
                                 StringUtil::itos(resrange.Nbins(), 3)));
    std::string line;
    for (int i=0;i<resrange.Nbins();++i) {
      line += StringUtil::ftos(resrange.middle(i), 8, 4);
    }
    // middle of resolution bins as 1/d^2
    output.logTab(0,LXML,
          StringUtil::MakeXMLtag("invresolsqBins",
                                 StringUtil::onespace(line)));
    k = 0;
    for (int i=1;i<nruns;++i) {
      std::string line1 = StringUtil::MakeXMLtag("runi",
                                     StringUtil::itos(runnumbers[i], 3));
      for (int j=0;j<i;++j) {
        std::string line = "<correlation>"+line1+
          StringUtil::MakeXMLtag("runj",
                                 StringUtil::itos(runnumbers[j], 3));
        line += StringUtil::MakeXMLtag("CC",
                                    StringUtil::ftos(mnCC_E2[k], 6, 3));
        line += StringUtil::MakeXMLtag("Number",
                                    StringUtil::itos(nmeanCC[k], 7));
        line += "\n";

        std::string s;
        for (int l=0;l<resrange.Nbins();++l) {
          s += " "+StringUtil::ftos(CC_E2_res[k][l].CC(),7,3);
        }
        line += "  "+StringUtil::MakeXMLtag("CCbyresolution", s, false);
        line += "</correlation>\n";
        output.logTab(0,LXML, StringUtil::onespace(line));
        k++;
      }
    }

    // close XML block
    output.logTab(0, LXML, "</runcorrelations>");

  }
  // ------------------------------------------------------------

};
