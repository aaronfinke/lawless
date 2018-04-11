//
// secondaryscalestats.cpp
//
// Statistics about the applied and calculated secondary scales
//

#include "secondaryscalestats.hh"
#include "string_util.hh"
#include "tablegraph.hh"

namespace scala {
  //--------------------------------------------------------------
  void SecondaryScaleStats::init(const ScaleModel& allScales)
  {
    nsecscales = allScales.numbersecondaryscales();
    if (nsecscales == 0) {return;}  // no secondary scales

    angleinterval = 10.0;  // degrees
    allscales = &allScales;  // pointer to scale model

    // Phi values
    nphi = Nint(360.0/angleinterval);
    phivalues.resize(nphi);
    for (size_t k=0; k<phivalues.size(); k++) {
      phivalues[k] = k * angleinterval;
    }
    // Theta values
    ntheta = Nint(180.0/angleinterval);
    thetavalues.resize(ntheta);
    for (size_t k=0; k<thetavalues.size(); k++) {
      thetavalues[k] = k * angleinterval;
    }

    secscales.resize(nsecscales);
    for (size_t k=0; k<secscales.size(); k++) {
      secscales[k].resize(ntheta, nphi);
    }

    // Make table of calculated corrections, and get range for analysis
    secscalerange.clear();
    secscalerange.SetNbin(1);
    for (size_t ks=0; ks<nsecscales;ks++) {
      for (size_t kt=0; kt<thetavalues.size(); kt++) {
        for (size_t kp=0; kp<phivalues.size(); kp++) {
          double c =
            allScales.secscale(ks,
                               clipper::Util::d2rad(thetavalues[kt]),
                               clipper::Util::d2rad(phivalues[kp]));
          if (c != 0.0) {
            c = 1.0/c;
            secscalerange.update(c);
          }
          secscales[ks](kt,kp) = c;
        }
      }
    }

    setupHistogram();  // initialise histogram from secscalerange
  }
  //--------------------------------------------------------------
  void SecondaryScaleStats::setupHistogram()
  // initialise histogram from secscalerange
  {
    const double HISTOBINWIDTH = 0.02;
    // Try to make a nice-looking range
    double lower = secscalerange.min() * 0.8;
    // multiple of bin width
    lower = Nint(lower/HISTOBINWIDTH) * HISTOBINWIDTH - 0.5*HISTOBINWIDTH;
    double upper = secscalerange.max() * 1.2;
    // multiple of bin width
    upper = Nint(upper/HISTOBINWIDTH) * HISTOBINWIDTH + 0.5*HISTOBINWIDTH;
    int nbins = Nint((upper - lower)/HISTOBINWIDTH);
    Range historange(lower, upper, true, nbins);
    secscalehisto.init(historange);
  }
  //--------------------------------------------------------------
  void SecondaryScaleStats::getScaleStats(hkl_unmerge_list& hkl_list,
                                          const bool& onlyUseSingletons)
  // go through list to get secondary scale statistics
  {
    if (nsecscales == 0) {return;}
    reflection this_refl;
    observation this_obs;
    double thetap, phip;
    double radinterval = clipper::Util::d2rad(angleinterval);
    double pi = clipper::Util::pi();
    int ks, kt, kp;

    sscount.resize(nsecscales);
    for (size_t k=0; k<secscales.size(); k++) {
      sscount[k].resize(ntheta, nphi, 0);
    }

    hkl_list.rewind();

    // Loop accepted reflections
    while (hkl_list.next_reflection(this_refl) >= 0) {
      while (this_refl.next_observation(this_obs) >= 0) {
        // get theta (0->pi), phi (-pi->+pi)
        this_obs.GetS2(thetap, phip);
        kt = std::min(std::max(int(thetap/radinterval), 0), ntheta-1);
        kp = std::min(std::max(int((phip + pi)/radinterval), 0), nphi-1);
        ks = allscales->secondaryScaleIndex(this_obs.run());
        double c = allscales->secscale(ks, thetap, phip);
        //      std::cout << this_obs.hkl_original().format() <<" "<<ks
        //                <<" "<<thetap<<" "<<kt<<"  "<<phip<<" "<<kp<<" ***\n";
        sscount[ks](kt,kp)++;
        secscalehisto.add(c);  // add to histogram
      }
    }
  }
  //--------------------------------------------------------------
  std::string SecondaryScaleStats::secscaletype(const std::string& sectype) const
  {
    std::string s;
    size_t j;
    if (sectype.find("SECONDARY") != std::string::npos) {
      s += "in goniometer frame (SECONDARY)";
    } else if ((j = sectype.find("ABSORPTION")) != std::string::npos) {
      s += "in crystal frame (ABSORPTION), pole " + sectype.substr(j+10);
    }
    return s;
  }
  //--------------------------------------------------------------
  void SecondaryScaleStats::PrintSecondaryCorrections(phaser_io::Output& output) const
  // Print secondary corrections as 2D array
  {
    if (nsecscales == 0) {return;}
    // Secondary
    std::string s = "\nSecondary scale corrections";
    std::vector<std::string> sectypes = allscales->secondaryscaletypes();
    if (sectypes.size() == 1) {
      s += ", "+secscaletype(sectypes[0]);
      output.logTab(0,LOGFILE, s);
    } else if (sectypes.size() > 1) {
      output.logTab(0,LOGFILE, s);
      for (size_t k=0; k<sectypes.size(); k++) {
        s = StringUtil::itos(k+1,2)+". "+secscaletype(sectypes[0]);
        output.logTab(1,LOGFILE, s);
      }
    }

    output.logTab(0,LOGFILE,
                  "\nCalculated for polar angles of theta (colatitude from 0 at N pole, 180 at S pole) and phi (longitude)\n");
    output.logTab(0,LOGFILE,
                  "Printed only for angular ranges containing data");

    output.logTabPrintf(0,LOGFILE,
          "\nRange of secondary corrections: %5.3f - %5.3f\n",
                        secscalerange.min(), secscalerange.max());

    // Print to log file, theta across
    std::string line;
    for (int j=0;j<nsecscales;++j) {
      output.logTabPrintf(0,LOGFILE,"\nSecondary scale number %3d\n", j+1);
      line = " Theta";
      for (size_t kt=0; kt<thetavalues.size(); kt++) {
        line += StringUtil::ftos(thetavalues[kt], 4, 0)+" ";
      }
      output.logTab(0,LOGFILE, line);
      output.logTab(0,LOGFILE, "  Phi");

      for (size_t kp=0; kp<phivalues.size(); kp++) {
        line = StringUtil::ftos(phivalues[kp], 5, 0)+" ";
        for (size_t kt=0; kt<thetavalues.size(); kt++) {
          std::string sv = StringUtil::ftos(secscales[j](kt,kp),5,2);
          //          std::cout << "j,kt,kp,ssc "<<j<<" "<<kt<<" "<<kp<<" "<< sscount[j](kt,kp)<<"\n";
          if (sscount.size() > 0) {
            if (sscount[j](kt,kp) == 0) {
              sv = "   - ";
            }
          }
          line += sv;
        }
        output.logTab(0,LOGFILE, line);
        /*  don't print numbers, for now anyway
        line = "     ";
        if (sscount.size() > 0) {
        for (size_t kt=0; kt<thetavalues.size(); kt++) {
        std::string sv = StringUtil::itos(sscount[j](kt,kp),5);
        //            std::cout << "j,kt,kp,ssc "<<j<<" "<<kt<<" "<<kp<<" "<< sscount[j](kt,kp)<<"\n";
        if (sscount[j](kt,kp) == 0) {
        sv = "    -";
        }
        line += sv;
        }
        output.logTab(0,LOGFILE, line);
        }
        */
      }
    }
    // XML theta slow, phi fast
    for (int j=0;j<nsecscales;++j) {
      output.logTab(0,LXML,"<SecondaryCorrection>");
      output.logTab(0,LXML,
                    StringUtil::MakeXMLtag("ScaleNumber", j+1));
      output.logTab(0,LXML,
                    StringUtil::MakeXMLtag("PhiValues",
                           StringUtil::FormatSaveVector(phivalues, false)));

      for (size_t kt=0; kt<thetavalues.size(); kt++) {
        line = StringUtil::ftos(thetavalues[kt], 5, 0)+" ";
        output.logTab(0,LXML,"<secscales>");
        output.logTab(0,LXML,
                      StringUtil::MakeXMLtag("Theta", thetavalues[kt]));
        std::vector<double> corrphi(nphi);
        std::vector<int> countphi(nphi);

        for (size_t kp=0; kp<phivalues.size(); kp++) {
          corrphi[kp] = secscales[j](kt,kp);
          if (sscount.size() > 0) {
            countphi[kp] = sscount[j](kt,kp);
          }
        }
        output.logTab(0,LXML,
                      StringUtil::MakeXMLtag("corrections",
                            StringUtil::FormatSaveVector(corrphi, false)));
        if (sscount.size() > 0) {
          output.logTab(0,LXML,
                        StringUtil::MakeXMLtag("counts",
                            StringUtil::FormatSaveVector(countphi, false)));
        }
        output.logTab(0,LXML,"</secscales>");

      }
      output.logTab(0,LXML,"</SecondaryCorrection>");
    }
    PrintHistogram(output);
  }
  //--------------------------------------------------------------
  void SecondaryScaleStats::PrintHistogram(phaser_io::Output& output) const
  {
    // $TABLE  start
    TableGraph table("Histogram of secondary corrections");
    table.StoreID("Graph-SecondaryCorrectionHistogram");

    TableGraphPlot graph("Histogram of secondary corrections");
    graph.AddLine(TableGraphPlotline(2,3,"red"));
    graph.SetYaxis("", true);  // Y from zero
    table.AddGraph(graph);

    std::vector<std::string> collabels;
    collabels.push_back("N");          // 1
    collabels.push_back("SecScale");   // 2
    collabels.push_back("Number");     // 3
    int nc = collabels.size();

    std::vector<bool> Zero(nc, false);
    Zero[2] = true;
    table.StoreColumnFields(collabels, Zero, "%5d%10.2f%9d\n");

    // Histogram
    std::vector<int> histocounts = secscalehisto.Counts();
    Range hrange = secscalehisto.HRange();
    int nbins = hrange.Nbins();
    ASSERT (nbins == int(histocounts.size()));

    for (int i=0;i<nbins;++i) {
      if (histocounts[i] > 0) {
        table.Line(nc, i+1,
                   hrange.middle(i), histocounts[i]);
      }
    }
    table.CloseTable();
    output.logTab(0,LOGFILE, "\n"+table.format());
    output.logTab(0,LXML,table.XMLformat());
  }
}
