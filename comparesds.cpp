//
// comparesds.cpp
//
// Get average sd(I)/samplesd(I) by resolution and intensity
// Only valid if sampleSD is used

#include "comparesds.hh"
#include "util.hh"
#include "tablegraph.hh"
#include "string_util.hh"
#include "printing.hh"

using phaser_io::LOGFILE;
using phaser_io::LXML;

namespace scala {
  //---------------------------------------------------------------
  CompareSDs::CompareSDs(const ResoRange& ResRange,
                         const IntensityBin& Irange,
                         const int& Minimumsample,
                         const Batchgroup& batchgroup)
  {
    init(ResRange, Irange, Minimumsample, batchgroup);
  }
  //---------------------------------------------------------------
  void CompareSDs::init(const ResoRange& ResRange,
                        const IntensityBin& Irange,
                        const int& Minimumsample,
                        const Batchgroup& batchgroup)
  {
    resrange = ResRange;
    int nresbins =  resrange.Nbins();
    meanSDratioreso.resize(nresbins);
    CCSDreso.resize(nresbins);
    int nintbins = Irange.NumberBins();
    meanSDratioInt.resize(nintbins);
    CCSDint.resize(nintbins);

    int nbtgp = batchgroup.numberofgroups(); // number of batch groups
    meanChiSqBtGp.resize(nbtgp);
    meanChiSqSampleBtGp.resize(nbtgp);

    pbatchgroup = &batchgroup;

    minimumsample = Minimumsample;

    meanChiSqreso.resize(nresbins); // mean(ChiSq-Variance) by resolution
    meanChiSqInt.resize(nintbins);  // mean(ChiSq-Variance) by intensity
    meanChiSqSamplereso.resize(nresbins); // mean(ChiSq-Sample) by resolution
    meanChiSqSampleInt.resize(nintbins);  // mean(ChiSq-Sample) by intensity

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

    //    // sample SDs for each observation
    //    std::vector<double> samplesds = selobs.sampleSDI();
    //    observation obs;
    //    selobs.reset_next();

    double ratio = selobs.SDvariance()/selobs.SDsample();
    meanSDratioreso[mres].Add(ratio); // mean(SD/sampleSD) by resolution
    meanSDratioInt[mint].Add(ratio);  // mean(SD/sampleSD) by intensity

    CCSDreso[mres].add(selobs.SDvariance(), selobs.SDsample());
    CCSDint[mint].add(selobs.SDvariance(), selobs.SDsample());
    //^^
    //    if (mint == int(CCSDint.size())-1) {
    //    std::cout <<"SDs "<<mint<<"  "
    //        << selobs.SDvariance()<< "  "<<selobs.SDsample()<<"\n";
    //    }

    double chiSqvariance = selobs.chiSq(false);
    double chiSqsample   = selobs.chiSq(true);

    //    std::cout << "chisq "<<chiSqvariance<<" "<<chiSqsample<<" "<<mres<<" "<<mint<<std::endl; //c^
    if (chiSqvariance != 0.0) {
      meanChiSqreso[mres].Add(chiSqvariance); // mean(ChiSq-Variance) by resolution
      meanChiSqInt[mint].Add(chiSqvariance);  // mean(ChiSq-Variance) by intensity
    }
    if (selobs.SampleSDused() && (chiSqsample != 0.0)) {
      meanChiSqSamplereso[mres].Add(chiSqsample); // mean(ChiSq-Sample) by resolution
      meanChiSqSampleInt[mint].Add(chiSqsample);  // mean(ChiSq-Sample) by intensity
    }

    // Analysis by Batchgroup
    std::vector<float> del2 = selobs.Delta2(); // individual SDs
    std::vector<float> del3 = selobs.Delta3(); // partitioned sample SDs
    int idx;
    observation this_obs;
    while ((idx=selobs.next_observation(this_obs)) >= 0) {  // loop all valid observations
      // idx points to valid observation in list, same as in delta2 array
      if (del2[idx] != 0.0) {
        int batchn = this_obs.Batch();  // batch number
        int jbatchgroup = pbatchgroup->batchgroup(batchn);
        meanChiSqBtGp[jbatchgroup].Add(del2[idx]*del2[idx]);
        if (del3[idx] != 0.0) {
          meanChiSqSampleBtGp[jbatchgroup].Add(del3[idx]*del3[idx]);
        }
      }
    }

    /*
    //^   various sanity checks
    std::vector<float> del1 = selobs.Deviations();
    std::vector<float> del2 = selobs.Delta2();
    std::vector<float> del3 = selobs.Delta3();
    ASSERT (del2.size() == del1.size());
    ASSERT (del3.size() == del1.size());

    std::vector<double> weights = selobs.weights(); // = 0.0 if not used
    //^
    int large = -1;
    const float limit = 6.0;
    for (size_t k=0; k<del1.size(); k++) {
      if (weights[k] > 0.0) {
        if (std::abs(del2[k]) > limit) {
          large = k;
        }
        csq1.Add(del1[k]*del1[k]);
        csq2.Add(del2[k]*del2[k]);
        csq3.Add(del3[k]*del3[k]);
        //      std::cout << del1[k]<<" "<<del2[k]<<" "<<del3[k]<<" del*\n";
      }
    }
    //^
    if (large >= 0) {
      reflection refl = selobs.Reflection();
      IsigI avIsigI = selobs.Average();

      std::vector<double> sampleSDI = selobs.sampleSDI();
      std::vector<float> delI = selobs.DelI();

      std::cout << "Large del2 "<<" "<<selobs.hkl().format() << " "<<large
                <<" "<<del2[large]
                <<" "<<avIsigI.I()<<" "<<avIsigI.sigI()<<std::endl;
      double fac = sqrt(double(del1.size())/double(del1.size()-1));

      for (size_t k=0; k<del1.size(); k++) {
        if (weights[k] > 0.0) {
          observation obs = refl.get_observation(k);
          std::string sm = "  ";
          if (std::abs(del2[k]) > limit) {
            sm = "$ ";}

          double d2 = fac*delI[k]/obs.ksigI();
          double d3 = fac*delI[k]/sampleSDI[k];
          printf("%s %3d %8.3f%8.3f%8.3f   %8.3f %8.1f%8.1f%8.1f %8.2f %8.3f %8.3f\n",
                 sm.c_str(),int(k), del1[k],del2[k],del3[k],
                 obs.Gscale(), obs.kI(),obs.ksigI(), sampleSDI[k],
                 delI[k], d2, d3);

      //      std::cout << "del* "<<del1[k]<<" "<<del2[k]<<" "<<del3[k]<<std::endl;
        }
      }
      std::cout <<"<<\n";
    }
    //^-
    */


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
    output.logTab(0,LOGFILE,
                  std::string("Also correlation coefficients between SDsample and SD(Imean),\n")+
                  "  and goodness of fit Mean(chi^2) estimated using both SD values");

    output.logTab(0,LOGFILE,
                  "\nSample SDs are used for reflections with more than "+
                  StringUtil::itos(minimumsample,2)+" observations\n");

    //^
    //      std::cout << "Msq " << csq1.Mean()<<" : "<< csq2.Mean()<<" : "
    //                << csq3.Mean()<<" : N " << csq1.Count()<<  std::endl;

    TableGraph table(" Compare SD estimates by resolution");
    table.StoreID("Graph-CompareSDsVsResolution");
    Range xrange = resrange; // x axis range to full resolution limit
    xrange.first() = 0.0;    // from 0

    TableGraphPlot graph("Average SD/SDsample");
    std::string description =
      std::string("SDsample from sample variance");
    graph.SetDescription(description);

    graph.AddLine(TableGraphPlotline(2,4)); // column numbers for x,y
    graph.AddLine(TableGraphPlotline(2,5)); // column numbers for x,y
    graph.AddLine(TableGraphPlotline(2,6)); // column numbers for x,y
    graph.AddLine(TableGraphPlotline(2,7)); // column numbers for x,y
    graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
    graph.SetYaxis("", true);  // y axis from 0 to maximum
    table.AddGraph(graph);

    // Define column labels
    std::vector<std::string> collabels;
    collabels.push_back("N");         // 1
    collabels.push_back("1/d^2");     // 2
    collabels.push_back("Dmid");      // 3
    collabels.push_back("SD/SDsample");      // 4
    collabels.push_back("ChiSqVarn");      // 5
    collabels.push_back("ChiSqSmpl");      // 6
    collabels.push_back("CC(SDs)");      // 7
    collabels.push_back("Number");      // 8
    int nc = collabels.size();
    // which columns should have "0.0" replaced by "-"
    bool z[] =
      {false, false, false, true, true, true, true, true};
    std::vector<bool> Zero(z, z+nc);
    std::string fmt = "%12.3f%12.3f%12.3f%12.3f%9d\n";
    // store labels, zero flags and format
    table.StoreColumnFields(collabels, Zero, "%3d%8.4f%7.2f"+fmt);

    int n=1;
    MeanValue overall;
    MeanValue meanChiSqoverall;
    MeanValue meanChiSqSampleoverall;
    correl_coeff CCSDresooverall;

    for (int i=0;i<resrange.Nbins();++i) {
      table.Line(nc, n++, resrange.middle(i), resrange.middleA(i),
                 meanSDratioreso[i].Mean(),
                 meanChiSqreso[i].Mean(), meanChiSqSamplereso[i].Mean(),
                 CCSDreso[i].CC(),
                 meanSDratioreso[i].Count());
      overall += meanSDratioreso[i];
      meanChiSqoverall += meanChiSqreso[i];
      meanChiSqSampleoverall += meanChiSqSamplereso[i];
      CCSDresooverall += CCSDreso[i];
    }
    table.CloseTable();

    // Output table to log file and XML
    output.logTab(0,LOGFILE, "\n"+table.format());
    output.logTab(0,LXML,table.XMLformat());

    fmt = "Overall:          "+fmt;
    output.logTabPrintf(0,LOGFILE,fmt.c_str(), overall.Mean(),
                        meanChiSqoverall.Mean(),
                        meanChiSqSampleoverall.Mean(),
                        CCSDresooverall.CC(), overall.Count());
  }
  //---------------------------------------------------------------
  void CompareSDs::printByIntensity(const IntensityBin& irange,
                                    phaser_io::Output& output) const
  {
    TableGraph table(" Compare SD estimates by intensity");
    table.StoreID("Graph-CompareSDsVsIntensity");

    TableGraphPlot graph("Average SD/SDsample and correlation");
    std::string description =
      std::string("SDsample from sample variance");
    graph.SetDescription(description);
    graph.AddLine(TableGraphPlotline(1,2));
    graph.AddLine(TableGraphPlotline(1,3));
    graph.SetYaxis("", true);  // y axis from 0 to maximum
    table.AddGraph(graph);

    std::vector<std::string> collabels;
    collabels.push_back("Imax");      // 1
    collabels.push_back("SD/SDsample");      // 2
    collabels.push_back("CC(SDs)");      // 3
    collabels.push_back("Number");      // 3
    int nc = collabels.size();

    bool z[] = {false, true, true, true};
    std::vector<bool> Zero(z, z+nc);
    std::string fmt = "%12.3f%12.3f%9d\n";
    table.StoreColumnFields(collabels, Zero, "%10.0f"+fmt);

    MeanValue overall;
    correl_coeff ccoverall;
    for (int i=0;i<irange.NumberBins();++i) {
      table.Line(nc, irange.mean(i),
                 meanSDratioInt[i].Mean(), CCSDint[i].CC(),
                 meanSDratioInt[i].Count());
      overall +=  meanSDratioInt[i];
      ccoverall += CCSDint[i];
    }

    table.CloseTable();

    // Output table to log file and XML
    output.logTab(0,LOGFILE, "\n"+table.format());
    output.logTab(0,LXML,table.XMLformat());

    fmt = "Overall:  "+fmt;
    output.logTabPrintf(0,LOGFILE,fmt.c_str(), overall.Mean(), ccoverall.CC(), overall.Count());
  }
  //---------------------------------------------------------------
  void CompareSDs::printByBatch(const int& datasetIndex,
                                phaser_io::Output& output) const
  {
    TableGraph table(" Compare SD estimates by batch");
    table.StoreID("Graph-CompareSDsVsBatch");

    size_t ngroups = pbatchgroup->numberofgroups();
    ASSERT (meanChiSqBtGp.size() == ngroups);

    // no symbol if too many points
    int symbolsize = GetSymbolSizeforNpoints(ngroups);
    TableGraphPlot graph("Chi^2 v. batch, for observation SDs and sample SDs");

    graph.AddLine(TableGraphPlotline(1,5,"red","",symbolsize));   // ChiSq-Variance
    graph.AddLine(TableGraphPlotline(1,6,"blue","",symbolsize));  // ChiSq-Sample
    graph.SetYaxis("", true);  // Y from zero


    std::vector<Batch> batches = pbatchgroup->batches(); // batch list
    // Breaks in X axis
    int xcolbr = 4;  // column for real batch number
    Xbreaks xbreaks(batches, datasetIndex);
    std::vector<Range> xbreaklist = xbreaks.get_breaks();
    graph.SetXbreak(xcolbr, xbreaklist);

    Range xrange(xbreaks.get_batchnumberrange());  // overall batch number range
    graph.SetXaxis("", false, xrange, true);
    table.AddGraph(graph);

    std::vector<std::string> collabels;
    collabels.push_back("N");         // 1
    collabels.push_back("Run");       // 2
    collabels.push_back("Phi");       // 3
    collabels.push_back("Batch");     // 4
    collabels.push_back("ChiSqVarn");      // 5
    collabels.push_back("ChiSqSmpl");      // 6
    collabels.push_back("Number");    // 7
    int nc = collabels.size();

    std::vector<bool> Zero(nc, false);
    Zero[4] = true;
    Zero[5] = true;
    table.StoreColumnFields(collabels, Zero,
                          "%5d%5d%8.2f%8d%10.3f%10.3f%10d\n");

    std::vector<Run> RunList = pbatchgroup->runlist();

    int n=1;
    for (size_t i=0;i<ngroups;++i) {  // print even batches that have no reflections
      int ib = pbatchgroup->batchserial(i);
      //  batch in this dataset
      if (batches[ib].datasetindex() == datasetIndex && batches[ib].Accepted()) {
        // ... but not rejected batches
        if (meanChiSqBtGp[i].Count() > 0) {
          table.Line(nc, ib+1,
                     RunList[batches[ib].RunIndex()].RunNumber(),
                     batches[ib].MidPhi(), batches[ib].num(),
                     meanChiSqBtGp[i].Mean(), meanChiSqSampleBtGp[i].Mean(),
                     meanChiSqBtGp[i].Count());
        }
        n++;
      }
    }
    table.CloseTable();
    output.logTab(0,LOGFILE, "\n"+table.format());
    output.logTab(0,LOGFILE,table.RawLabels());
    output.logTab(0,LXML,table.XMLformat());
  }
  //---------------------------------------------------------------
}
