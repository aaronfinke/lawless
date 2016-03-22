//
// comparesds.cpp
//
// Get average sd(I)/samplesd(I) by resolution and intensity
// Only valid if sampleSD is used

#include "comparesds.hh"
#include "util.hh"
#include "tablegraph.hh"

using phaser_io::LOGFILE;
using phaser_io::LXML;

namespace scala {
  //---------------------------------------------------------------
  CompareSDs::CompareSDs(const ResoRange& ResRange,
                         const IntensityBin& Irange)
  {
    init(ResRange, Irange);
  }
  //---------------------------------------------------------------
  void CompareSDs::init(const ResoRange& ResRange,
                        const IntensityBin& Irange)
  {
    resrange = ResRange;
    int nresbins =  resrange.Nbins();
    meanSDratioreso.resize(nresbins);
    int nintbins = Irange.NumberBins();
    meanSDratioInt.resize(nintbins);
  }
  //---------------------------------------------------------------
  void CompareSDs::add(SelectedObservations& selobs,
                       const int& mres, const int& mint)
  // add in contributions
  {
    if (selobs.Number() == 0 || !selobs.SampleSDused()) {
      // no used observations or too few for sampleSD
      return;
    }

    // sample SDs for each observation
    std::vector<double> samplesds = selobs.sampleSDI();

    observation obs;
    selobs.reset_next();

    double ratio = selobs.SDvariance()/selobs.SDsample();
    meanSDratioreso[mres].Add(ratio); // mean(SD/sampleSD) by resolution
    meanSDratioInt[mint].Add(ratio);  // mean(SD/sampleSD) by intensity
  }
  //---------------------------------------------------------------
  void CompareSDs::printByResolution(phaser_io::Output& output) const
  // Also prints header messages, so call before byIntensity
  {
    output.logTab(0,LOGFILE,
                  std::string("\n\nComparison of sample SD and estimated SD\n")+
                  "========================================\n");
    output.logTab(0,LOGFILE,
                  std::string("\nWe have two estimates of sigma(Imean):")+
                  "SDsample from the sample variance of each reflection;\n"+
                  " and SD(Imean) from the individual SDs\n");

    output.logTab(0,LOGFILE,
  "Here the average ratio SDsample/SD(Imean) is analysed against resolution and intensity\n");

    TableGraph table(" Compare SD estimates by resolution");
    table.StoreID("Graph-CompareSDsVsResolution");
    Range xrange = resrange; // x axis range to full resolution limit
    xrange.first() = 0.0;    // from 0

    TableGraphPlot graph("Average SDsample/SD");
    std::string description =
      std::string("SDsample from sample variance");
    graph.SetDescription(description);

    graph.AddLine(TableGraphPlotline(2,4)); // column numbers for x,y
    graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
    graph.SetYaxis("", true);  // y axis from 0 to maximum
    table.AddGraph(graph);

    // Define column labels
    std::vector<std::string> collabels;
    collabels.push_back("N");         // 1
    collabels.push_back("1/d^2");     // 2
    collabels.push_back("Dmid");      // 3
    collabels.push_back("SDsample/SD");      // 4
    int nc = collabels.size();
    // which columns should have "0.0" replaced by "-"
    bool z[] =
      {false, false, false, true};
    std::vector<bool> Zero(z, z+4);
    std::string fmt = "%12.3f\n";
    // store labels, zero flags and format
    table.StoreColumnFields(collabels, Zero, "%3d%8.4f%7.2f"+fmt);

    int n=1;
    MeanValue overall;
    for (int i=0;i<resrange.Nbins();++i) {
      table.Line(nc, n++, resrange.middle(i),
                 resrange.middleA(i), meanSDratioreso[i].Mean());
      overall += meanSDratioreso[i];
    }
    table.CloseTable();

    // Output table to log file and XML
    output.logTab(0,LOGFILE, "\n"+table.format());
    output.logTab(0,LXML,table.XMLformat());

    fmt = "Overall:          "+fmt;
    output.logTabPrintf(0,LOGFILE,fmt.c_str(), overall.Mean());
  }
  //---------------------------------------------------------------
  void CompareSDs::printByIntensity(const IntensityBin& irange,
                                    phaser_io::Output& output) const
  {
    TableGraph table(" Compare SD estimates by intensity");
    table.StoreID("Graph-CompareSDsVsIntensity");

    TableGraphPlot graph("Average SDsample/SD");
    std::string description =
      std::string("SDsample from sample variance");
    graph.SetDescription(description);
    graph.AddLine(TableGraphPlotline(1,2));
    graph.SetYaxis("", true);  // y axis from 0 to maximum
    table.AddGraph(graph);

    std::vector<std::string> collabels;
    collabels.push_back("Imax");      // 1
    collabels.push_back("SDsample/SD");      // 2
    int nc = collabels.size();

    bool z[] = {false, true};
    std::vector<bool> Zero(z, z+nc);
    std::string fmt = "%12.3f\n";
    table.StoreColumnFields(collabels, Zero, "%10.0f"+fmt);

    MeanValue overall;
    for (int i=0;i<irange.NumberBins();++i) {
      table.Line(nc, irange.bounds(i).second,
                 meanSDratioInt[i].Mean());
      overall +=  meanSDratioInt[i];
    }

    table.CloseTable();

    // Output table to log file and XML
    output.logTab(0,LOGFILE, "\n"+table.format());
    output.logTab(0,LXML,table.XMLformat());

    fmt = "Overall:  "+fmt;
    output.logTabPrintf(0,LOGFILE,fmt.c_str(), overall.Mean());
  }
  //---------------------------------------------------------------
}
