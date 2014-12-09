// applyscales.cpp
//
// calculate & store all scales
//

#include "scalemodel.hh"
#include "hkl_unmerge.hh"
#include "scala_util.hh"
#include "applyscales.hh"
#include "tablegraph.hh"
#include "string_util.hh"
#include "printing.hh"

using phaser_io::LOGFILE;
using phaser_io::LXML;

namespace scala {
  //----------------------------------------------------------------
  ApplyScales::ApplyScales(const ScaleModel& AllScales,
                           hkl_unmerge_list& hkl_list,
                           const bool& onlyUseSingletons)
  {
    scale(AllScales, hkl_list, onlyUseSingletons);
  }
  //----------------------------------------------------------------
  void ApplyScales::scale(const ScaleModel& AllScales,
                          hkl_unmerge_list& hkl_list,
                          const bool& onlyUseSingletons)
  //! Apply scales to all data, return mean(I) within resolution limits
  // if onlyUseSingletons true, do not attempt to apply scales to overlaps
  {
    reflection this_refl;
    observation this_obs;
    hkl_list.rewind();
    resrange = hkl_list.ResLimRange();  // resolution limits
    meani.clear();  // mean scaled I within resolution limits

    int nresbin =  resrange.Nbins();
    meansdk.assign(nresbin, MeanValue());  // Mean sd(1/g)
    meansdkbatch.assign(hkl_list.num_batches(), MeanValue());  // Mean sd(1/g) by batch
    meanrelsddiff.assign(nresbin, MeanValue());  // Mean sigIcorrected/sigI
    maxrelsddiff = -1000.0;
    maxsdk = -1000.0;
    bool usesdparameter = (AllScales.parameterSDusage() != scala::ScaleSpecification::NONE);

    // * * * * Loop reflections
    // loop all reflections unconditionally
    for (int jref=0;jref<hkl_list.num_reflections();++jref) {
      this_refl = hkl_list.get_reflection(jref);
      Rtype invresolsq = this_refl.invresolsq();
      int mres = resrange.tbin(invresolsq);
      bool inrange = (mres >= 0); // true if in reso limits

      //  Loop all observations
      for (int i=0;i<this_refl.num_observations();++i) {
        this_obs = this_refl.get_observation(i);
        // apply scale to observation
        AllScales.ScaleObs(this_obs, invresolsq, onlyUseSingletons);
        this_refl.replace_observation(this_obs);
        if (inrange) {
          meani.Add(this_obs.kI());

          if (usesdparameter) {
            // Statistics on effect of parameter variance
            Rtype sigI0 = this_obs.sigI(); // sigI before scaling
            Rtype sigI  = this_obs.ksigI(); // sigI after scaling & correction
            ASSERT (sigI > 0.0);
            Rtype gscale = this_obs.Gscale(); // inverse scale g
            Rtype varg = this_obs.varGscale(); // inverse scale g

            // sd(1/g) = sd(g)/g^2; sdk = sd(k)/k = sd(k) * g = sd(g)/g
            Rtype sdk = sqrt(varg)/gscale;
            meansdk[mres].Add(sdk);
            if (sdk > maxsdk) {
              maxsdk = sdk;
              maxsdkobs.init(this_obs, invresolsq);
            }
            int batchn = this_obs.Batch();  // batch number
            int jbatch = hkl_list.batch_serial(batchn); // batch serial
            meansdkbatch[jbatch].Add(sdk);

            // relative change of sd(I)
            Rtype relsddiff =std::max(0.0f,(sigI - sigI0/gscale))/sigI;
            maxrelsddiff = std::max(maxrelsddiff, relsddiff);
            // ratio of sds after and before allowance for sd(1/g)
            meanrelsddiff[mres].Add(relsddiff);
          }
        } // in range
      }
      hkl_list.replace_reflection(this_refl); // store updated reflection
    }
  }
  //----------------------------------------------------------------
  void ApplyScales::print(const std::vector<Batch>& batches,
                          phaser_io::Output& output) const
  {
    printResolution(output);
    //printBatch(batches, output);
  }
  //----------------------------------------------------------------
  void ApplyScales::printResolution(phaser_io::Output& output) const
  {
    if (maxrelsddiff < -999.0) {return;}
    std::string s =
      "\nEffect of allowing for parameter variance in estimation of sig(I)\n";
    s +=
      "================================================================\n\n";
    s += "Corrected sd'(I') = Sqrt[sd(I')^2 + I'^2 (sd(k)/k)^2];  I' = kI\n";

    s += FormatOutput::logTabPrintf(0, "\n%s %8.3f\n",
         "Maximum relative correction relSD = (sd'(I') - sd(I'))/sd'(I') = ",
                                    maxrelsddiff);
    s += FormatOutput::logTabPrintf(0, "\n%s %8.4f\n",
         "Maximum relative sd(scale), sd(k)/k = ",
                                    maxsdk);
    s += FormatOutput::logTabPrintf(1, "%s\n     %s\n",
                                    " for observation:",
                                    maxsdkobs.format().c_str());
    s += "\nAnalysis by resolution\n\n";

    output.logTab(0,LXML,StringUtil::MakeXMLtag
                  ("MaximumRelSDcorrection", maxrelsddiff, 8, 3));
    output.logTab(0,LXML,StringUtil::MakeXMLtag
                  ("MaximumRelSDscale", maxsdk, 8, 4));

    TableGraph table("Effect of parameter variance on sd(I)");
    table.StoreID("Graph-ParameterVariance");

    TableGraphPlot graph("Mean sd(k)/k and relative Delta(sd(I)) vs. resolution");

    Range xrange = resrange; // x axis range to full resolution limit
    xrange.first() = 0.0;    // from 0

    graph.AddLine(TableGraphPlotline(2,4));
    graph.AddLine(TableGraphPlotline(2,5));
    graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
    graph.SetYaxis("", true);  // y axis from 0 to maximum
    table.AddGraph(graph);

    std::vector<std::string> collabels;
    collabels.push_back("N");         // 1
    collabels.push_back("1/d^2");     // 2
    collabels.push_back("Dmid");      // 3
    collabels.push_back("sd(k)/k");    // 4
    collabels.push_back("relSD");     // 5
    collabels.push_back("Number");   // 6
    bool z[] = {false, false, false, true, true, true};
    std::vector<bool> Zero(z, z+6);
    std::string fmt = "%8.4f %8.3f %8d\n"; // excluding 1st 3 columns
    table.StoreColumnFields(collabels, Zero, "%3d%8.4f%7.2f"+fmt);
    int nc = collabels.size();

    int nbins = meansdk.size();
    MeanValue meansdkall, meanrelsddiffall;
    int ntotal = 0;
    int n = 1;
    for (int mres=0;mres<nbins;++mres) {
      table.Line(nc, n++, resrange.middle(mres), resrange.middleA(mres),
                 meansdk[mres].Mean(), meanrelsddiff[mres].Mean(),
                 meansdk[mres].Count());
      meansdkall += meansdk[mres];
      meanrelsddiffall += meanrelsddiff[mres];
      ntotal += meansdk[mres].Count();
    }
    s += table.format();
    output.logTab(0,LXML, table.XMLformat());

    std::string leader = "Overall:          ";
    fmt = leader+fmt;
    s += FormatOutput::logTabPrintf(0, fmt.c_str(),
                                    meansdkall.Mean(),
                                    meanrelsddiffall.Mean(),
                                    ntotal);

    output.logTab(0,LOGFILE, s);
  }
  //----------------------------------------------------------------
  void ApplyScales::printBatch(const std::vector<Batch>& batches,
                               phaser_io::Output& output) const
  {
    if (maxrelsddiff < -999.0) {return;}
    std::string s = "\n\nAnalysis by batch\n\n";

    TableGraph table("Effect of parameter variance on sd(I), by batch");
    table.StoreID("Graph-ParameterVarianceByBatch");

    TableGraphPlot graph("Mean sd(k)/k and relative Delta(sd(I)) vs. batch");

    // Breaks in X axis
    int xcolbr = 2;  // column for real batch number
    Xbreaks xbreaks(batches, -1);
    std::vector<Range> xbreaklist = xbreaks.get_breaks();
    // overall batch number range, for XML plot
    Range xrange(xbreaks.get_batchnumberrange());
    Range xnrange; // dummy for $TABLE range
    graph.SetXbreak(xcolbr, xbreaklist, xrange);
    graph.SetXaxis("", false, xnrange, true);
    graph.SetYaxis("", true);  // y axis from 0 to maximum
    graph.AddLine(TableGraphPlotline(1,3));
    table.AddGraph(graph);

    std::vector<std::string> collabels;
    collabels.push_back("N");         // 1
    collabels.push_back("Batch");     // 2
    collabels.push_back("sd(k)/k");    // 3
    collabels.push_back("Number");   // 4
    int nc = collabels.size();
    bool z[] = {false, false, true, false};
    std::vector<bool> Zero(z, z+nc);
    std::string fmt = "%5d %7d %8.4f %8d\n";
    table.StoreColumnFields(collabels, Zero, fmt);

    int n = 1;
    int nbatches = batches.size();
    for (int i=0;i<nbatches;++i) {
      table.Line(nc, n++, batches[i].num(),
                 meansdkbatch[i].Mean(), meansdkbatch[i].Count());
    }
    s += table.format();
    output.logTab(0,LXML, table.XMLformat());
    output.logTab(0,LOGFILE, s);
  }
  //----------------------------------------------------------------
  UnusualObservation::UnusualObservation(const observation& this_obs, const Rtype& Invresolsq)
  {
    init(this_obs, Invresolsq);
  }
  //----------------------------------------------------------------
  void UnusualObservation::init(const observation& this_obs, const Rtype& Invresolsq)
  {
    hkl = this_obs.hkl_original();
    invresolsq= Invresolsq;
    batch = this_obs.Batch();
    I = this_obs.kI();       // scaled I
    sigI0 = this_obs.sigI(); // sigI before scaling
    sigI  = this_obs.ksigI(); // sigI after scaling & correction
    gscale = this_obs.Gscale(); // inverse scale g
    varg = this_obs.varGscale(); // inverse scale g
  }
  //----------------------------------------------------------------
  std::string UnusualObservation::format() const
  {
    if (batch <= 0) {return "";}
    std::string s = "HKL "+StringUtil::itos(hkl[0],4)+","+
      StringUtil::itos(hkl[1],4)+","+StringUtil::itos(hkl[2],4);
    if (invresolsq > 0.0) {
      s += ", resolution(A): "+StringUtil::ftos(1.0/sqrt(invresolsq),7, 2);
    }
    s += ", batch: "+StringUtil::itos(batch,7);
    std::string s1;
    if (gscale > 0.0) {
      s += ", scale: "+StringUtil::ftos(1.0/gscale, 8, 3);
      if (varg >= 0.0) {
        s += ", sd(scale): "+StringUtil::ftos(sqrt(varg)/(gscale*gscale),8,3);
      }
      s1 += "I': "+StringUtil::ftos(I, 8, 1);
      s1 += ", sd(I'): "+StringUtil::ftos(sigI0/gscale, 8, 1);
      s1 += ", corrected sd'(I'): "+StringUtil::ftos(sigI, 8, 1);
    }
    return StringUtil::onespace(s)+"\n         "+StringUtil::onespace(s1);
  }
  //----------------------------------------------------------------
}
