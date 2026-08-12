// printing.cpp

#include "printing.hh"
#include "version.hh"
#include "tablegraph.hh"
#include "numbercomplete.hh"
#include "halfdataset.hh"
#include "string_util.hh"
#include "anisotropy.hh"

using phaser_io::LOGFILE;
using phaser_io::LXML;

//--------------------------------------------------------------
float FractionN(const float& scale, const int& num, const int& den)
// return scale * (num/den), = 0 if den = 0  (eg scale = 1 or 100)
{
  if (den == 0) {
    return 0.0f;
  } else {
    return scale * (float(num)/float(den));
  }
}
//--------------------------------------------------------------
float FractionN(const float& scale, const float& num, const int& den)
// return scale * (num/den), = 0 if den = 0
{
  if (den == 0) {
    return 0.0f;
  } else {
    return scale * (num/float(den));
  }
}
//--------------------------------------------------------------
float FractionN(const float& scale, const double& num, const int& den)
// return scale * (num/den), = 0 if den = 0
{
  if (den == 0) {
    return 0.0f;
  } else {
    return scale * (num/float(den));
  }
}
//--------------------------------------------------------------
void PrintTitle( phaser_io::Output& output)
{
  const int sp = 52; // interior width of the box (chars between the '*' borders)
  // centre a string within the box interior, return the full bordered line
  auto boxline = [&](const std::string& text) {
    int len = int(text.size());
    if (len > sp) len = sp;
    int n1 = (sp - len)/2;
    int n2 = sp - len - n1;
    return std::string("        *") + std::string(n1,' ')
         + text.substr(0,len) + std::string(n2,' ') + "*\n";
  };
  std::string border = "        " + std::string(sp+2,'*') + "\n";

  output.logTab(0,LOGFILE, "\n"+border);
  output.logTab(0,LOGFILE, boxline(""));
  output.logTab(0,LOGFILE, boxline("LAWLESS"));
  output.logTab(0,LOGFILE, boxline("(it's AIMLESS for Laue data)"));
  output.logTab(0,LOGFILE, boxline(PROGRAM_VERSION));
  output.logTab(0,LOGFILE, boxline(""));
  output.logTab(0,LOGFILE, boxline("Scaling & analysis of unmerged intensities"));
  output.logTab(0,LOGFILE, boxline("(now with wavelength normalization!)"));
  output.logTab(0,LOGFILE, boxline(""));
  output.logTab(0,LOGFILE, boxline("Phil Evans MRC LMB, Cambridge"));
  output.logTab(0,LOGFILE, boxline("(with tiny contribs from Aaron Finke, ESS, DK)"));
  output.logTab(0,LOGFILE, boxline(""));
  output.logTab(0,LOGFILE, border+"\n");
}
//--------------------------------------------------------------
void PrintFileInfoToXML(const std::string& StreamName,
                        const std::string& FileName,
                        const Scell& cell,
                        const std::string& SpaceGroupName,
                        const std::vector<std::string>& columnlabels,
                        phaser_io::Output& output)
{
  if (output.doXmlout())
    {
      output.logTab(0, LXML,
                     "<ReflectionFile stream=\""+StreamName+
                     "\" name=\""+FileName+"\">\n");
      output.logTab(0, LXML,
                     cell.xml());
      output.logTab(0, LXML,
                    StringUtil::MakeXMLtag("SpacegroupName", SpaceGroupName));
      if (columnlabels.size() > 0) {
        std::string s;
        for (size_t k=0; k<columnlabels.size(); k++) {
          s += " " + columnlabels[k];
        }
        output.logTab(0, LXML,
                      StringUtil::MakeXMLtag("ColumnLabelsUsed", s));
      }

      output.logTab(0, LXML,
                     "</ReflectionFile>");
    }
}
//--------------------------------------------------------------
void PrintOutlierSettings(const all_controls& controls, phaser_io::Output& output)
{
  output.logTab(0,LOGFILE,"\nOutlier rejection parameters:");

  output.logTab(0,LOGFILE,"In scaling:");
  if (controls.outlierScale.GetOutlierPolicy() == OutlierControl::NOREJECT) {
    output.logTab(0,LOGFILE,"No outlier rejection");
  } else {
    output.logTab(0,LOGFILE,controls.outlierScale.Reject(ALL).format());
    output.logTab(0,LOGFILE,controls.outlierScale.EMaxTest().format());
    output.logTab(0,LOGFILE, "Outlier test weight type: "+
	   WeightType::formatWeightType(controls.outlierScale.weightType()));
  }

  output.logTab(0,LOGFILE,"\nIn merging:");
  if (controls.outlierScale.GetOutlierPolicy() == OutlierControl::NOREJECT) {
    output.logTab(0,LOGFILE,"No outlier rejection");
  } else {
    output.logTab(0,LOGFILE,controls.outlierMerge.Reject(ALL).format());
    output.logTab(0,LOGFILE,controls.outlierMerge.EMaxTest().format());
    output.logTab(0,LOGFILE, "Outlier test weight type: "+
	   WeightType::formatWeightType(controls.outlierMerge.weightType()));
  }
  output.logTab(0,LOGFILE,"\n");
}
//--------------------------------------------------------------
//--------------------------------------------------------------
FitBfactorLines::FitBfactorLines(const std::vector<Batch>& batches,
                                 const Batchgroup& batchgroup,
                                 const std::vector<Run>& RunList,
                                 const int& datasetIndex,
                                 const std::vector<float>& bfacbatch,
                                 const std::vector<int>& nbfacrun)
{
  size_t ngroups = batchgroup.numberofgroups();
  ASSERT (bfacbatch.size() == ngroups);

  bfdecaybatch.assign(bfacbatch.size(), 0.0);
  bsloperun.assign(RunList.size(), 0.0);     // for each run
  b0run.assign(RunList.size(), 0.0);         // for each run
  scales.assign(RunList.size(), 0.0);         // for each run
  std::vector<int> batch0run(RunList.size(),-1);

  for (size_t irun=0;irun<RunList.size();++irun) { // loop runs
    if (RunList[irun].DatasetIndex() == datasetIndex) { // is run in this dataset?
      if (nbfacrun[irun] > 1) {
        // Extract Bfactors for batches belonging to this run
        LinearFit linefit;
        float w = 1.0;
        // B slope is determined in A^2/batch, but using the ranges we can convert it to
        // A^2/degree
        Range batchserialrange;  // range of batch serials in this run
        Range phirange;          // range of phi in this run

        for (int jgroup=0;jgroup<batchgroup.numberofgroups();++jgroup) {  // loop groups
          int ib = batchgroup.batchserial(jgroup);  // batch number for group
          if (batches[ib].Accepted() && (batches[ib].RunIndex() == int(irun))) { // in this run
            if (batch0run[irun] < 0) {  // 1st batch in run
              batch0run[irun] = ib;
            }
            float x = ib - batch0run[irun];
            // x = batch serial in run, y = B
            linefit.add(x, bfacbatch[jgroup], w);
            batchserialrange.update(x);
            phirange.update(batches[ib].MidPhi());
          }
        }
        // scaling from batch serial to phi
        scales[irun] = phirange.AbsRange()/batchserialrange.AbsRange();
        float scale = 1.0;

        RPair r = linefit.result();
        bsloperun[irun] = r.first;
        b0run[irun] = r.second;
      }
    }
  } // end loop runs

  for (size_t i=0;i<ngroups;++i) { // loop groups
    int ib = batchgroup.batchserial(i);
    // this dataset & accepted
    if (batches[ib].datasetindex() == datasetIndex && batches[ib].Accepted()) {
      int irun = batches[ib].RunIndex();
      bfdecaybatch[i] = float(ib-batch0run[irun]) * bsloperun[irun] + b0run[irun];  // from straight line
    }
  } // end loop batch groups
}
//--------------------------------------------------------------
void PrintScalesByBatch(const PxdName& dataset_pxd,
                        const std::vector<Batch>& batches,
                        const Batchgroup& batchgroup,
                        const std::vector<Run>& RunList,
                        const int& datasetIndex,
                        const std::vector<float>& scale0batch,
                        const std::vector<float>& bfacbatch,
                        const std::vector<int>& nbfacrun,
                        const std::vector<MeanSD>& scalebatch,
                        phaser_io::Output& output)
{
  output.logTab(0,LOGFILE,
                std::string("\nScale factors analysed by Batch for each dataset\n")+
                              "================================================\n\n"+
                "Note that 0k below is calculated for the centre of each rotation range,\n"+
                "at theta = 0 (for the B-factor)\n"+
                "Mn(k) is average applied scale, including any input scale\n"+
                "0k is the scale calculated excluding any input scale\n");
  output.logTab(0,LOGFILE, batchgroup.format());

  size_t ngroups = batchgroup.numberofgroups();
  ASSERT (bfacbatch.size() == ngroups);

  // no symbol if too many points
  int symbolsize = GetSymbolSizeforNpoints(ngroups);

  // Fit straight line to B factors within each run
  FitBfactorLines fit(batches, batchgroup, RunList, datasetIndex, bfacbatch, nbfacrun);
  std::vector<float> bfdecaybatch = fit.DecayBatch();

  output.logTab(0,LOGFILE,
                "\nBdecay comes from a straight line fit to the B-factors within each run");
  for (size_t irun=0;irun<RunList.size();++irun) { // loop runs
    if (RunList[irun].DatasetIndex() == datasetIndex) { // is run in this dataset?
      if (nbfacrun[irun] > 1) {  // more than 1 Bfactor in run
        output.logTabPrintf(1,LOGFILE,
                            "For run number %4d, slope of B (A^2/degree) %8.3f\n",
                            RunList[irun].RunNumber(), fit.Slope(irun));
      }
    }
  }

  // $TABLE  start
  TableGraph table(" === Scales v rotation range, "+dataset_pxd.dname());
  table.StoreID("Graph-ScalesVsRotationRange");

  TableGraphPlot graph("Mn(k) & 0k (theta=0) v. batch");
  std::string description =
    "Mn(k) (red) is the mean scale over the whole resolution range for each Batch.";
  description += " 0k is the scale at the lowest resolution, ie excluding the relative B-factor";
  graph.SetDescription(description);
  graph.AddLine(TableGraphPlotline(1,5,"red","",symbolsize));  // Mn(k)
  graph.AddLine(TableGraphPlotline(1,6,"blue","",symbolsize));  // 0k
  graph.SetYaxis("", true);  // Y from zero
  // Breaks in X axis
  int xcolbr = 4;  // column for real batch number
  Xbreaks xbreaks(batches, datasetIndex);
  std::vector<Range> xbreaklist = xbreaks.get_breaks();
  graph.SetXbreak(xcolbr, xbreaklist);

  Range xrange(xbreaks.get_batchnumberrange());  // overall batch number range
  graph.SetXaxis("", false, xrange, true);
  table.AddGraph(graph);

  graph.init("Relative Bfactor & Decay v. batch");
  description = "The relative B-factor is largely a radiation damage correction.";
  description += " Negative values below perhaps -10 may indicate sever radiation damage.";
  description += " Bdecay is a straight-line fit to the B-factors";
  graph.SetDescription(description);
  graph.AddLine(TableGraphPlotline(1,8,"red","",symbolsize));  // Bfactor
  graph.AddLine(TableGraphPlotline(1,9,"blue","",symbolsize));  // Bdecay
  graph.SetXbreak(xcolbr, xbreaklist);
  graph.SetXaxis("", false, xrange, true);
  graph.SetYaxis("", false);  // Y not from zero
  table.AddGraph(graph);

  std::vector<std::string> collabels;
  collabels.push_back("N");         // 1
  collabels.push_back("Run");       // 2
  collabels.push_back("Phi");       // 3
  collabels.push_back("Batch");     // 4
  collabels.push_back("Mn(k)");     // 5
  collabels.push_back("0k");        // 6
  collabels.push_back("Number");    // 7
  collabels.push_back("Bfactor");   // 8
  collabels.push_back("Bdecay");    // 9
  int nc = collabels.size();

  std::vector<bool> Zero(nc, false);
  Zero[4] = true;
  Zero[9] = true;
  Zero[10] = true;
  table.StoreColumnFields(collabels, Zero,
                          "%5d%5d%8.2f%8d%11.4f%11.4f%10d%10.4f%10.4f\n");

  int n=1;
  for (size_t i=0;i<ngroups;++i) {  // print even batches that have no reflections
    int ib = batchgroup.batchserial(i);
    //  batch in this dataset
    if (batches[ib].datasetindex() == datasetIndex && batches[ib].Accepted()) {
      // ... but not rejected batches
      if (scalebatch[i].Count() > 0) {
        table.Line(nc, ib+1,
                   RunList[batches[ib].RunIndex()].RunNumber(),
                   batches[ib].MidPhi(), batches[ib].num(),
                   scalebatch[i].Mean(),
                   scale0batch[i], scalebatch[i].Count(),
                   bfacbatch[i], bfdecaybatch[i]);
        n++;
      }
    }
  }
  table.CloseTable();
  output.logTab(0,LOGFILE, "\n"+table.format());
  output.logTab(0,LOGFILE,table.RawLabels());
  output.logTab(0,LXML,table.XMLformat());
}
//--------------------------------------------------------------
void PrintDeviationsByBatch(const PxdName& dataset_pxd,
                            const std::vector<Batch>& batches,
                            const Batchgroup& batchgroup,
                            const int& datasetIndex,
                            const std::vector<MeanSD>& imeanbatch,
                            const std::vector<MeanSD>& rmsDbatch,
                            const std::vector<Rfactor>& rmergebatch,
                            const std::vector<Rfactor>& rmergebatchsmoothed,
                            const std::vector<int>& rejectedbatch,
                            const std::vector<float>& batchcompleteness,
                            const std::vector<float>& batchanomcompleteness,
                            const std::vector<float>& batchmultiplicity,
                            const std::vector<double>& maxresbatch,
                            const std::vector<double>& maxresbatchsmoothed,
                            const std::vector<MeanSD>&  meanChiSqBatch,
                            const std::vector<MeanSD>&  meanChiSqBatch2,
                            const double& MinimumIoverSigma,
                            const int& nbatchsmooth,
                            const ResoRange& ResRange,
                            const RejectFlags& rejflags,
                            phaser_io::Output& output)
{
  output.logTab(0,LOGFILE,
                std::string("\n\nAgreement between batches\n")+
                                "=========================\n\n");
  output.logTab(0,LOGFILE,
                std::string(" Rmerge in this table is the difference from Mn(Imean),\n")+
                "  but in later tables Rmerge is the difference from Mn(I+),Mn(I-)\n");

  std::string corelimit = StringUtil::Strip(clipper::String(rejflags.sdrej, 3, 2) + "sd");
  output.logTab(0,LOGFILE,
                "Mean Chi^2 values are calculated for all data, and also");
  output.logTab(0,LOGFILE,
                " excluding observations with large deviations >"+corelimit+
                ", including for reflections measured only twice (Chi^2c)");

  // Smoothed values
  bool smoothR = false;
  bool smoothMaxRes = false;
  if (nbatchsmooth > 1) {
    smoothR = true;
    ASSERT (maxresbatch.size() == maxresbatchsmoothed.size());
    for (size_t i=0;i<maxresbatch.size();++i) {
      if (std::abs(maxresbatch[i] - maxresbatchsmoothed[i]) > 0.002) {
        smoothMaxRes = true;  // don't use if all same as unsmoothed
      }
    }
    if (smoothMaxRes) {
      output.logTabPrintf(0,LOGFILE,
                          "\n SmRmerge and SmMaxRes in table are smoothed over %3d batch groups\n",
                          nbatchsmooth);
    } else {
      output.logTabPrintf(0,LOGFILE,
                          "\n SmRmerge in table is smoothed over %3d batches\n",
                          nbatchsmooth);
    }
  }

  // Get range of resolutions for maximum resolution plot
  float highres = ResRange.ResHigh();
  float lowres = 0.0;
  for (size_t i=0;i<maxresbatch.size();++i) {
    if (imeanbatch[i].Count() > 0) {
      lowres = Max(lowres, maxresbatch[i]);  // lowest resolution limit
    }
  }
  // choose suitable limits
  float res1=0.0;
  float rinc = 0.5;  // step
  float offset = 0.1;
  while (res1<highres-offset) {res1 += rinc;}
  res1 = Max(0.0, res1-rinc);
  float res2 = res1;
  while (res2<lowres+offset) {res2 += rinc;}
  //  std::cout << "graph reso "<< highres <<" "<<lowres<<" "<<res1<<" "<<res2
  //        <<" " << maxresbatch.size()<<"\n"; //^
  //  int nb = 0;
  //  for (size_t i=0;i<batches.size();++i) {
  //    if (batches[i].Accepted()) { // ... but not rejected batches
  //      nb++;
  //    } // count actual batches
  //  }
  scala::Range yrange(res1, res2);

  // no symbol if too many points
  size_t ngroups = batchgroup.numberofgroups();
  int symbolsize = GetSymbolSizeforNpoints(ngroups);

  TableGraph table
    (" Analysis against all Batches for all runs, "+dataset_pxd.dname());
  table.StoreID("Graph-StatsVsBatch");

  TableGraphPlot graph("Rmerge v Batch for all runs");
  std::string description =
    "Increase of Rmerge towards the end of a run probably indicates radiation damage";

  if (smoothR) { // smoothed, 2 lines, smoothed first
    graph.AddLine(TableGraphPlotline(1,15,"red",
                                     "",symbolsize,false,
                                     "Solid",2)); // smoothed
    graph.AddLine(TableGraphPlotline(1,6,"blue","",symbolsize)); // unsmoothed
    description += ". The red line is smoothed over adjacent batches";
  } else {
    graph.AddLine(TableGraphPlotline(1,6,"blue","",symbolsize)); // unsmoothed
  }
  graph.SetDescription(description);
  graph.SetYaxis("", true);  // Y from zero
  // Breaks in X axis
  int xcolbr = 2;  // column for real batch number
  Xbreaks xbreaks(batches, datasetIndex);
  std::vector<Range> xbreaklist = xbreaks.get_breaks();
  // overall batch number range, for XML plot
  Range xrange(xbreaks.get_batchnumberrange());
  Range xnrange; // dummy for $TABLE range

  graph.SetXbreak(xcolbr, xbreaklist, xrange);
  graph.SetXaxis("", false, xnrange, true);
  table.AddGraph(graph);

  graph.init("Filtered Mean(Chi^2)(<"+corelimit+") v Batch");
  description = std::string("Goodness of fit Mean(chi^2) should be ~1.0 over all ranges");
  graph.SetDescription(description);
  graph.AddLine(TableGraphPlotline(1,14,"red","",symbolsize));  // eg < 5sd
  graph.SetXaxis("", false, xnrange, true);
  graph.SetYaxis("", true);  // Y from zero
  graph.SetXbreak(xcolbr, xbreaklist, xrange);
  table.AddGraph(graph);

  graph.init("Mean(Chi^2), Mean(Chi^2)(<"+corelimit+") v Batch");
  description = std::string("Goodness of fit Mean(chi^2) should be ~1.0 over all ranges");
  graph.SetDescription(description);
  graph.AddLine(TableGraphPlotline(1,13,"blue","",symbolsize)); // all
  graph.AddLine(TableGraphPlotline(1,14,"red","",symbolsize));  // eg < 5sd
  graph.SetXaxis("", false, xnrange, true);
  graph.SetYaxis("", true);  // Y from zero
  graph.SetXbreak(xcolbr, xbreaklist, xrange);
  table.AddGraph(graph);

  graph.init("Cumulative %completeness & Anom%cmpl v Batch");
  description = std::string("In conjunction with radiation damage indicators, ")+
    "cumulative completeness may indicate a suitable point for cutting back poor data."+
    " The blue line is the completeness of anomalous differences";
  graph.SetDescription(description);
  graph.AddLine(TableGraphPlotline(1,10,"blue","",symbolsize));  // anomalous completeness
  graph.AddLine(TableGraphPlotline(1,9,"red","",symbolsize));  // completeness
  graph.SetXaxis("", false, xnrange, true);
  graph.SetYaxis("", true);  // Y from zero
  graph.SetXbreak(xcolbr, xbreaklist, xrange);
  table.AddGraph(graph);

  std::string s = "Maximum resolution limit, I/sigma > "+
    StringUtil::Strip(StringUtil::ftos(MinimumIoverSigma,5,1));
  graph.init(s);
  description =
    "The resolution limit estimate is the point at which I/sig(I) falls below "+
    StringUtil::Strip(StringUtil::ftos(MinimumIoverSigma,5,1))+
    ", red line smoothed over adjacent batches. "+
    "A sharp increase probably indicates radiation damage";
  graph.SetDescription(description);

  if (smoothMaxRes) {
    graph.AddLine(TableGraphPlotline(1,16,"red","",symbolsize,false,"Solid",2));  // smoothed
    graph.AddLine(TableGraphPlotline(1,11,"blue","",symbolsize));  // unsmoothed
  } else {
    graph.AddLine(TableGraphPlotline(1,11,"blue","",symbolsize));  // unsmoothed
  }
  graph.SetXaxis("", false, xnrange, true);
  graph.SetYaxis("", false, yrange);
  graph.SetXbreak(xcolbr, xbreaklist, xrange);
  table.AddGraph(graph);

  graph.init("Imean & RMS Scatter");
  description = "Major deviations may indicate bad images";
  graph.SetDescription(description);

  graph.AddLine(TableGraphPlotline(1,3,"red","",symbolsize));
  graph.AddLine(TableGraphPlotline(1,4,"blue","",symbolsize));
  graph.SetXaxis("", false, xnrange, true);
  graph.SetYaxis("", true);  // Y from zero
  graph.SetXbreak(xcolbr, xbreaklist, xrange);
  table.AddGraph(graph);

  graph.init("Imean/RMS scatter");
  description = "Major deviations may indicate bad images";
  graph.SetDescription(description);
  graph.AddLine(TableGraphPlotline(1,5,"red","",symbolsize));
  graph.SetXaxis("", false, xnrange, true);
  graph.SetYaxis("", true);  // Y from zero
  graph.SetXbreak(xcolbr, xbreaklist, xrange);
  table.AddGraph(graph);

  graph.init("Number of rejects");
  graph.AddLine(TableGraphPlotline(1,8,"red","",symbolsize));
  graph.SetXaxis("", false, xnrange, true);
  graph.SetYaxis("", true);  // Y from zero
  graph.SetXbreak(xcolbr, xbreaklist, xrange);
  table.AddGraph(graph);

  graph.init("Cumulative multiplicity");
  graph.AddLine(TableGraphPlotline(1,12,"red","",symbolsize));
  graph.SetXaxis("", false, xnrange, true);
  graph.SetYaxis("", true);  // Y from zero
  graph.SetXbreak(xcolbr, xbreaklist, xrange);
  table.AddGraph(graph);

  std::vector<std::string> collabels;
  collabels.push_back("N");          // 1
  collabels.push_back("Batch");      // 2
  collabels.push_back("Mn(I)");      // 3
  collabels.push_back("RMSdev");     // 4
  collabels.push_back("I/rms");      // 5
  collabels.push_back("Rmerge");     // 6
  collabels.push_back("Number");     // 7
  collabels.push_back("Nrej");       // 8
  collabels.push_back("Cm%poss");    // 9
  collabels.push_back("AnoCmp");     // 10
  collabels.push_back("MaxRes");     // 11
  collabels.push_back("CMlplc");     // 12
  collabels.push_back("Chi^2");      // 13
  collabels.push_back("Chi^2c");     // 14 without outliers
  if (smoothR) { // if we have smoothed stats as well
    collabels.push_back("SmRmerge");     // 15
  }
  if (smoothMaxRes) {
    collabels.push_back("SmMaxRes");     // 16
  }
  int nc = collabels.size();

  bool z[] =
    {false, false, true, true, true, true, false, false, true, true, true, false, true, true, true, true};
  std::vector<bool> Zero(z, z+nc);
  std::string lineformat = "%5d %7d %8.1f %8.1f %6.2f %7.3f %9d %5d %7.1f %7.1f %6.2f %6.2f %7.2f %7.2f";
  if  (smoothR) {lineformat += " %8.3f";}
  if  (smoothMaxRes) {lineformat += " %8.2f";}
  lineformat += "\n";
  table.StoreColumnFields(collabels, Zero, lineformat);

  int n=1;
  for (size_t i=0;i<ngroups;++i) { // print all groups even if they have no observations
    int ib = batchgroup.batchserial(i);
    //  batch in this dataset
    if (batches[ib].datasetindex() == datasetIndex && batches[ib].Accepted()) {
      // ... but not rejected batches
      if (rmergebatch[i].result().count > 0) {
        float r = 0.0;
        if (rmsDbatch[i].Mean() > 0) {
          r = imeanbatch[i].Mean()/sqrt(rmsDbatch[i].Mean());
        }
        if (nbatchsmooth == 1) { // no smoothed stats
          table.Line(nc, ib+1, batches[ib].num(),
                     imeanbatch[i].Mean(), sqrt(rmsDbatch[i].Mean()),
                     r,
                     rmergebatch[i].R(),
                     rmergebatch[i].result().count,
                     rejectedbatch[i],
                     100.*batchcompleteness[i],
                     100.*batchanomcompleteness[i],
                     maxresbatch[i],
                     batchmultiplicity[i],
                     meanChiSqBatch[i].Mean(),
                     meanChiSqBatch2[i].Mean());
        } else if (smoothMaxRes) {
          table.Line(nc, ib+1, batches[ib].num(),
                     imeanbatch[i].Mean(), sqrt(rmsDbatch[i].Mean()),
                     r,
                     rmergebatch[i].R(),
                     rmergebatch[i].result().count,
                     rejectedbatch[i],
                     100.*batchcompleteness[i],
                     100.*batchanomcompleteness[i],
                     maxresbatch[i],
                     batchmultiplicity[i],
                     meanChiSqBatch[i].Mean(),
                     meanChiSqBatch2[i].Mean(),
                     rmergebatchsmoothed[i].R(),
                     maxresbatchsmoothed[i]);
        } else{
          table.Line(nc, ib+1, batches[ib].num(),
                     imeanbatch[i].Mean(), sqrt(rmsDbatch[i].Mean()),
                     r,
                     rmergebatch[i].R(),
                     rmergebatch[i].result().count,
                     rejectedbatch[i],
                     100.*batchcompleteness[i],
                     100.*batchanomcompleteness[i],
                     maxresbatch[i],
                     batchmultiplicity[i],
                     meanChiSqBatch[i].Mean(),
                     meanChiSqBatch2[i].Mean(),
                     rmergebatchsmoothed[i].R());
        }
        n++;
      }
    }  // rejected batches
  } // batch loop
  table.CloseTable();
  output.logTab(0,LOGFILE, "\n"+table.format());
  output.logTab(0,LOGFILE,table.RawLabels());
  output.logTab(0,LXML,table.XMLformat());
}
//--------------------------------------------------------------
void PrintComparisonToReferenceByBatch(const PxdName& dataset_pxd,
                                       const std::vector<Batch>& batches,
                                       const Batchgroup& batchgroup,
                                       const int& datasetIndex,
                                       const int& nbatchsmooth,
                                       const std::vector<Rfactor>& rreferencebatch,
                                       const std::vector<MeanValue>& ccreferencebatch,
                                       const std::vector<int>& numberinCC,
                                       const std::vector<Rfactor>& rreferencebatchsmoothed,
                                       const std::vector<MeanValue>&
                                           ccreferencebatchsmoothed,
                                       const std::vector<MeanValue>& meanIrefbatch,
                                       const std::vector<MeanValue>& meanIobsbatch,
                                       phaser_io::Output& output)
{
  output.logTab(0,LOGFILE,
                std::string("\n\nAgreement with reference data, analysed by batch\n")+
                                "================================================\n\n");

  output.logTab(0,LOGFILE,
                std::string("Rref   is Sum(Iobs - k.Iref) / Sum(Iobs)\n")+
                "CCref  is mean of CC(Iobs, k.Iref) averaged over resolution ranges\n\n");

  // Smoothed values
  bool smoothR = false;
  if (nbatchsmooth > 1) {
    smoothR = true;
    output.logTabPrintf(0,LOGFILE,
                        "\n SmRref and SmCCref in table are smoothed over %3d batches\n",
                        nbatchsmooth);
  }

  //  int nb = 0;
  //  for (size_t i=0;i<batches.size();++i) {
  //    if (batches[i].Accepted()) { // ... but not rejected batches
  //      nb++;
  //    } // count actual batches
  //  }

  // no symbol if too many points
  size_t ngroups = batchgroup.numberofgroups();
  int symbolsize = GetSymbolSizeforNpoints(ngroups);

  TableGraph table
    (" Comparison to reference data for all Batches for all runs, "+dataset_pxd.dname());
  table.StoreID("Graph-RefStatsVsBatch");

  TableGraphPlot graph("Rref and CCref v Batch for all runs");
  if (smoothR) { // smoothed, 2 lines
    graph.AddLine(TableGraphPlotline(1,7,"red","",symbolsize,false)); // R smoothed
    graph.AddLine(TableGraphPlotline(1,3,"blue","",symbolsize,false)); // R unsmoothed
    TableGraphPlotline ccline(1,8,"black","",symbolsize,false);
    TableGraphPlotline ccline2(1,5,"green","",symbolsize,false); // CC unsmoothed
    ccline.SetRHaxis();
    graph.AddLine(ccline);
    ccline2.SetRHaxis();
    graph.AddLine(ccline2); // CC unsmoothed
  } else {
    graph.AddLine(TableGraphPlotline(1,3,"blue","",symbolsize)); // unsmoothed
    TableGraphPlotline ccline2(1,5,"green","",symbolsize,false); // CC unsmoothed
    ccline2.SetRHaxis();
    graph.AddLine(ccline2); // CC unsmoothed
  }
  graph.SetYaxis("", true);  // Y from zero
  graph.SetRightYaxis("", true,Range(0.0,1.0));
  // Breaks in X axis
  int xcolbr = 2;  // column for real batch number
  Xbreaks xbreaks(batches, datasetIndex);
  std::vector<Range> xbreaklist = xbreaks.get_breaks();
  graph.SetXbreak(xcolbr, xbreaklist);
  Range xrange(xbreaks.get_batchnumberrange());  // overall batch number range
  graph.SetXaxis("", false, xrange, true);
  table.AddGraph(graph);

  graph.init("<Iobs>, <Iref> v Batch for all runs");
  int c1 = 7;
  if (smoothR) {c1 = 9;}
  graph.AddLine(TableGraphPlotline(1,c1,"red","",symbolsize)); // <Iobs>
  graph.AddLine(TableGraphPlotline(1,c1+1,"blue","",symbolsize)); // <Iref>
  graph.SetYaxis("", true);  // Y from zero
  graph.SetXaxis("", false, xrange, true);
  table.AddGraph(graph);

  std::vector<std::string> collabels;
  collabels.push_back("N");          // 1
  collabels.push_back("Batch");      // 2
  collabels.push_back("Rref");       // 3
  collabels.push_back("Number");     // 4
  collabels.push_back("CCref");      // 5
  collabels.push_back("CCnumber");   // 6
  if (smoothR) { // if we have smoothed stats as well
    collabels.push_back("SmRref");   // 7
    collabels.push_back("SmCCref");  // 8
  }
  collabels.push_back("<Iobs>");     // 7, 9
  collabels.push_back("<Iref>");     // 8, 10
  int nc = collabels.size();

  bool z[] =
    {false, false, true, true, true, true, true, true, true, true};
  std::vector<bool> Zero(z, z+nc);
  std::string lineformat = "%8.3f %8d %8.3f %8d"; // excluding 1st 2 columns
  if (smoothR) {lineformat += " %8.3f %8.3f";}
  lineformat += " %9d %9d\n";
  table.StoreColumnFields(collabels, Zero, "%5d %7d "+lineformat);

  // Overall
  Rfactor Rf;
  MeanValue CC;
  int nCC = 0;

  int n=1;
  for (size_t i=0;i<ngroups;++i) { // print all groups even if they have no observations
    int ib = batchgroup.batchserial(i);
    //  batch in this dataset
    if (batches[ib].datasetindex() == datasetIndex && batches[ib].Accepted()) {
      // ... but not rejected batches
      if (smoothR) { // smoothed stats
        table.Line(nc, ib+1, batches[ib].num(),
                   rreferencebatch[i].R(),
                   rreferencebatch[i].result().count,
                   ccreferencebatch[i].Mean(),
                   numberinCC[i],
                   rreferencebatchsmoothed[i].R(),
                   ccreferencebatchsmoothed[i].Mean(),
                   Nint(meanIobsbatch[i].Mean()),
                   Nint(meanIrefbatch[i].Mean()));
      } else {
        table.Line(nc, ib+1, batches[ib].num(),
                   rreferencebatch[i].R(),
                   rreferencebatch[i].result().count,
                   ccreferencebatch[i].Mean(),
                   numberinCC[i],
                   Nint(meanIobsbatch[i].Mean()),
                   Nint(meanIrefbatch[i].Mean()));
      }
      n++;
      Rf += rreferencebatch[i];
      CC += ccreferencebatch[i];
      nCC += numberinCC[i];
      //^
      //      double ratio = meanIobsbatch[i].Mean()/meanIrefbatch[i].Mean();
      //      std::cout << "<Iref>, <Iobs> " << ratio <<" "<<
      //        meanIrefbatch[i].Mean() <<" "<<
      //        meanIobsbatch[i].Mean()<<std::endl;
      //^-
    }  // rejected batches
  } // batch loop
  table.CloseTable();
  output.logTab(0,LOGFILE, "\n"+table.format());
  lineformat = "Overall:      "+lineformat;
  output.logTabPrintf(0,LOGFILE,lineformat.c_str(),
                      Rf.R(), Rf.result().count, CC.Mean(), nCC,
                      Rf.R(), CC.Mean());
  output.logTab(0,LOGFILE,table.RawLabels());
  output.logTab(0,LXML,table.XMLformat());
}
//--------------------------------------------------------------
void PrintComparisonToReferenceByReso(const PxdName& dataset_pxd,
				      const ResoRange& ResRange,
				      const std::vector<Rfactor>& rreferencereso,
				      std::vector<std::vector<correl_coeff> > ccreferencebatch,
				      const std::vector<MeanValue>& meanIrefreso,
				      const std::vector<MeanValue>& meanIobsreso,
				      phaser_io::Output& output)
{
  output.logTab(0,LOGFILE,
                std::string("\n\nAgreement with reference data, analysed by resolution\n")+
                                "=====================================================\n\n");

  int nresbin = ccreferencebatch[0].size();
  std::vector<correl_coeff> ccreferencereso(nresbin);
  size_t nbatch = ccreferencebatch.size();
  for (int mres=0;mres<nresbin;mres++) {
    for (size_t kb=0; kb<nbatch; kb++) { 
      ccreferencereso[mres] += ccreferencebatch[kb][mres];
    }
  }

  output.logTab(0,LOGFILE,
                std::string("Rref   is Sum(Iobs - k.Iref) / Sum(Iobs)\n")+
                "CCref  is  CC(Iobs, k.Iref) \n\n");

  int symbolsize = 1;
  TableGraph table
    (" Comparison to reference data by resolution for all runs, "+dataset_pxd.dname());
  table.StoreID("Graph-RefStatsVsReso");
  TableGraphPlot graph("Rref and CCref v Resolution for all runs");
  graph.AddLine(TableGraphPlotline(1,4,"blue",""));  // R-factor
  TableGraphPlotline ccline(1,6,"green","",symbolsize,false); // CC
  ccline.SetRHaxis();
  graph.AddLine(ccline);
  graph.SetYaxis("", true);  // Y from zero
  graph.SetRightYaxis("", true,Range(0.0,1.0));
  Range xrange = ResRange; // x axis range to full resolution limit
  xrange.first() = 0.0;    // from 0
  graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
  table.AddGraph(graph);

  graph.init("<Iobs>, <Iref> v Resolution for all runs");
  int c1 = 8;
  graph.AddLine(TableGraphPlotline(1,c1,"red","",symbolsize)); // <Iobs>
  graph.AddLine(TableGraphPlotline(1,c1+1,"blue","",symbolsize)); // <Iref>
  graph.SetYaxis("", true);  // Y from zero
  graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
  table.AddGraph(graph);

  std::vector<std::string> collabels;
  collabels.push_back("N");         // 1
  collabels.push_back("1/d^2");     // 2
  collabels.push_back("Dmid");      // 3
  collabels.push_back("Rref");      // 4
  collabels.push_back("Number");    // 5
  collabels.push_back("CCref");     // 6
  collabels.push_back("CCnumber");  // 7
  collabels.push_back("<Iobs>");    // 8
  collabels.push_back("<Iref>");    // 9
  int nc = collabels.size();

  bool z[] =
    {false, false, true, true, true, true, true, true, true};
  std::vector<bool> Zero(z, z+nc);
  std::string lineformat = " %8.4f %8.3f %8d %8.3f %8d %9d %9d\n";
  table.StoreColumnFields(collabels, Zero, "%5d %7.3f "+lineformat);

  // Overall
  Rfactor Rf;
  MeanValue CC;
  int nCC = 0;

  int n=1;
  for (size_t i=0;i<nresbin;++i) {
    table.Line(nc, n++, ResRange.middle(i),
               ResRange.middleA(i),
	       rreferencereso[i].R(),
	       rreferencereso[i].result().count,
	       ccreferencereso[i].CC(),
	       ccreferencereso[i].Number(),
	       Nint(meanIobsreso[i].Mean()),
	       Nint(meanIrefreso[i].Mean()));
      n++;
      Rf += rreferencereso[i];
      CC.Add(ccreferencereso[i].CC());
      nCC += ccreferencereso[i].Number();;
  } // reso loop
  table.CloseTable();
  output.logTab(0,LOGFILE, "\n"+table.format());
  lineformat = "Overall:      "+lineformat;
  output.logTabPrintf(0,LOGFILE,lineformat.c_str(),
                      Rf.R(), Rf.result().count, CC.Mean(), nCC,
                      Rf.R(), CC.Mean());
  output.logTab(0,LOGFILE,table.RawLabels());
  output.logTab(0,LXML,table.XMLformat());
}
//--------------------------------------------------------------
void PrintDeviationsByResolution(const PxdName& dataset_pxd,
                                 const ResoRange& ResRange,const bool& Anom,
                                 const std::vector<Rfactor>& rmergeRes,
                                 const std::vector<Rfactor>& rmergeResFull,
                                 const std::vector<Rfactor>& rmeasRes,
                                 const std::vector<Rfactor>& rpimRes,
                                 const std::vector<MeanSD>&  imeanRes,
                                 const std::vector<MeanSD>&  rmsDRes,
                                 const std::vector<MeanSD>&  avSdRes,
                                 const std::vector<MeanSD>&  mnIsdRes,
                                 const std::vector<MeanSD>&  biasRes,
                                 const std::vector<MeanSD>&  biasIRes,
                                 const std::vector<MeanSD>&  meanChiSqRes,
                                 const std::vector<MeanSD>&  meanChiSqRes2,
                                 const double& MinimumIoverSigma,
                                 SummaryStatistics& summarystatistics,
                                 const RejectFlags& rejflags,
                                 phaser_io::Output& output)
{

  output.logTab(0,LOGFILE,
                std::string
                ("\n\nAnalysis by 4sinTheta/Lambda^2 bins (all statistics use Mn(I+),Mn(I-)etc)\n")+
        "=========================================================================\n");

  std::string corelimit = StringUtil::Strip(clipper::String(rejflags.sdrej, 3, 2) + "sd");

  output.logTab(0,LOGFILE,
                std::string("\n Rmrg    :- conventional Rmerge = Sum(|Ihl - < Ih >|)/Sum(< Ih >)\n")+
                " Rcum    :- Rmrg up to this range\n"+
                " Rfull   :- Rmrg for fully-recorded observations only\n"+
                " Rmeas   :- multiplicity-independent R = Sum(Sqrt(N/(N-1))(|Ihl - < Ih >|))/Sum(< Ih >)\n"+
                " Rpim    :- Precision-indicating R = Sum(Sqrt(1/(N-1))(|Ihl - < Ih >|))/Sum(< Ih >)\n"+
                " Nmeas   :- Number of observations used in statistics\n"+
                " AvI     :- unmerged Ihl averaged in bin < Ihl >\n"+
                " RMSdev  :- rms scatter of observations from mean < Ih >\n"+
                " I/RMS   :- < Ihl > / rms scatter  = Av_I/RMSdev\n"+
                " sd      :- average standard deviation derived from experimental SDs, after\n"+
                "             application of SdFac SdB SdAdd 'correction' terms\n"+
                " Mn(I/sd):- average < merged< Ih >/sd(< Ih >) > ~= signal/noise\n"+
                " Frcbias :- partial bias = Mean( Mn(If) - Ip )/Mean( Mn(I) )\n"+
                "             for mixed sets only (If is a full if present, else the\n"+
                "             partial with the smallest number of parts)\n"+
                " Chi^2   :- mean goodness of fit, all data\n"+
                " Chi^2c  :- mean goodness of fit, excluding large differences\n\n"
                );
  if (Anom) {
    output.logTab(0,LOGFILE,
          "All statistics in this table are within I+ or I- sets (anomalous on)");
  } else {
    output.logTab(0,LOGFILE,
          "All statistics in this table are relative to the overall mean I+/- (anomalous off)");
  }

  output.logTab(0,LOGFILE,
                "Mean Chi^2 values are calculated for all data, and also");
  output.logTab(0,LOGFILE,
                " excluding observations with large deviations >"+corelimit+
                ", including for reflections measured only twice");


  TableGraph table(" Analysis against resolution, "+dataset_pxd.dname());
  table.StoreID("Graph-StatsVsResolution");

  Range xrange = ResRange; // x axis range to full resolution limit
  xrange.first() = 0.0;    // from 0
  std::vector<Range> yranges(5);   // for each graph
  // Get y ranges for each graph (if loggraph would accept just an xrange, wouldn't need to do this)
  // (as qloggraph does)
  bool anyFulls = false;
  for (int i=0;i<ResRange.Nbins();++i) {
    float frcbias = 0.0;
    if (biasIRes[i].Count() > 0) {
      frcbias = biasRes[i].Mean()/biasIRes[i].Mean();
    }
    // Graphs are
    // 1. I/sigma, Mean Mn(I)/sd(Mn(I))  (cols 13,14)
    // 2. Rmerge, Rfull, Rmeas, Rpim v Resolution  (cols 4,5,7,8)
    // 3. Average I, RMSdeviation and Sd (cols 10,11,12)
    // 4. Fractional bias (col 15)
    // 5. Chi^2 and filtered Chi^2 (cols 16,17)
    double sd = sqrt(Max(rmsDRes[i].Mean(), 0.0));
    double Iovsd = 0.0;
    if ( sd > 0.0) {Iovsd = imeanRes[i].Mean()/sd;}

    yranges[0].update(Iovsd);  // I/sigma
    yranges[0].update(mnIsdRes[i].Mean()); //Mn(I/sd)

    yranges[1].update(rmergeRes[i].R());   // Rmerge
    yranges[1].update(rmergeResFull[i].R());   // Rfull
    yranges[1].update(rmeasRes[i].R());    // Rmeas
    yranges[1].update(rpimRes[i].R());     // Rpim

    yranges[2].update(Max(0.0,imeanRes[i].Mean())); // AvI
    yranges[2].update(sqrt(Max(0.0,rmsDRes[i].Mean())));  // RMSdeviation

    yranges[3].update(frcbias);
    yranges[4].update(Max(0.0, sqrt(rmsDRes[i].Mean()))); // RMSdev
    yranges[4].update(Max(0.0,avSdRes[i].Mean()));  // Sd
    if (rmergeRes[i].R() > 0.0) {
      anyFulls = true;
    }
  } // end line loop


  TableGraphPlot graph("I/sigma, Mean Mn(I)/sd(Mn(I))");  // 1st graph
  std::string description =
    std::string("I/sigma = I/RMS = I/(rms scatter before merging). ")+
    "Mean(I/sd) after averaging, ~= signal/noise";
  graph.SetDescription(description);

  graph.AddLine(TableGraphPlotline(2,13)); // column numbers for x,y
  graph.AddLine(TableGraphPlotline(2,14));
  graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
  graph.SetYaxis("", true, yranges[0]);  // y axis from 0 to maximum
  table.AddGraph(graph);

  corelimit = StringUtil::Strip(clipper::String(rejflags.sdrej, 3, 2) + "sd");
  graph.init("Filtered Mean(Chi^2)(<"+corelimit+") v Resolution");
  description = "Goodness of fit Mean(chi^2) should be close to 1.0";
  graph.SetDescription(description);
  graph.AddLine(TableGraphPlotline(2,17));  // eg < 5sd
  graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
  graph.SetYaxis("", true);  // y axis from 0
  table.AddGraph(graph);

  graph.init("Mean(Chi^2), Mean(Chi^2)(<"+corelimit+") v Resolution");
  description = "Goodness of fit Mean(chi^2) should be close to 1.0";
  graph.SetDescription(description);
  graph.AddLine(TableGraphPlotline(2,16));  // all
  graph.AddLine(TableGraphPlotline(2,17));  // eg < 5sd
  graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
  graph.SetYaxis("", true);  // y axis from 0
  table.AddGraph(graph);

  graph.init("Rmerge, Rfull, Rmeas, Rpim v Resolution");  // 2nd graph
  description = "Rmerge, classic R-factor; ";
  if (anyFulls) {
    description +=  "Rfull, for fully-recordeds only; ";
  }
  description += "Rmeas, multiplicity-weighted R; Rpim, precision-indicating R";
  graph.SetDescription(description);
  graph.AddLine(TableGraphPlotline(2,4));
  graph.AddLine(TableGraphPlotline(2,5));
  graph.AddLine(TableGraphPlotline(2,7));
  graph.AddLine(TableGraphPlotline(2,8));
  graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
  graph.SetYaxis("", true, yranges[1]);  // y axis from 0 to maximum
  table.AddGraph(graph);

  graph.init("Average I, RMSdeviation and Sd");  // 3rd graph
  description = "RMSdev is RMS scatter, sd is average corrected sig(I) estimate";
  graph.SetDescription(description);
  graph.AddLine(TableGraphPlotline(2,10)); // AvI
  TableGraphPlotline line1(2,11);
  line1.SetRHaxis();
  graph.AddLine(line1); // RMSdev
  line1.init(2,12);
  line1.SetRHaxis();
  graph.AddLine(line1); // Mn(I/sd)
  graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
  graph.SetYaxis("", true, yranges[2]);  // y axis from 0 to maximum
  graph.SetRightYaxis("", true, yranges[4]);
  table.AddGraph(graph);

  graph.init("Fractional bias");  // 4th graph
  description =
    std::string("FrcBias = Mean( Mn(If) - Ip )/Mean( Mn(I) ), ")+
                "where If is a full if present, else the partials "+
                "with the smallest number of parts";
  graph.SetDescription(description);
  graph.AddLine(TableGraphPlotline(2,15));
  graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
  graph.SetYaxis("", false, yranges[3]);  // y axis from minimum to maximum
  table.AddGraph(graph);

  // Define column labels
  std::vector<std::string> collabels;
  collabels.push_back("N");         // 1
  collabels.push_back("1/d^2");     // 2
  collabels.push_back("Dmid");      // 3
  collabels.push_back("Rmrg");      // 4
  collabels.push_back("Rfull");     // 5
  collabels.push_back("Rcum");      // 6
  collabels.push_back("Rmeas");     // 7
  collabels.push_back("Rpim");      // 8
  collabels.push_back("Nmeas");     // 9
  collabels.push_back("AvI");       // 10
  collabels.push_back("RMSdev");    // 11
  collabels.push_back("sd");        // 12
  collabels.push_back("I/RMS");     // 13
  collabels.push_back("Mn(I/sd)");  // 14
  collabels.push_back("FrcBias");   // 15
  collabels.push_back("Chi^2");     // 16
  collabels.push_back("Chi^2c");    // 17
  int nc = collabels.size();
  // which columns should have "0.0" replaced by "-"
  bool z[] =
    {false, false, false, true, true, true, true, true, false, false, true, true,
    true, true, true, true, true};
  std::vector<bool> Zero(z, z+nc);
  std::string fmt = "%7.3f%7.3f%7.3f%7.3f%7.3f%9d%9d%7d%7d%7.1f%9.1f%9.3f%7.2f%7.2f\n"; // excluding 1st 3 columns
  // store labels, zero flags and format
  table.StoreColumnFields(collabels, Zero, "%3d%8.4f%7.2f"+fmt);

  Rfactor Rcum;  // cumulative R
  Rfactor Rfull, Rmeas, Rpim;
  MeanSD Imean, rmsD, avSd, mnIsd, bias, biasI, chisq, chisq2;

  int n=1;
  double Iovsd = 0.0;
  for (int i=0;i<ResRange.Nbins();++i) {
    Rcum += rmergeRes[i];
    float frcbias = 0.0;
    if (biasIRes[i].Count() > 0) {
      frcbias = biasRes[i].Mean()/biasIRes[i].Mean();
    }
    double sd = sqrt(Max(rmsDRes[i].Mean(), 0.0));
    Iovsd = 0.0;
    if ( sd > 0.0) {Iovsd = imeanRes[i].Mean()/sd;}
    // Store each table line
    table.Line(nc, n++, ResRange.middle(i),
               ResRange.middleA(i),
               rmergeRes[i].R(), rmergeResFull[i].R(), Rcum.R(),
               rmeasRes[i].R(), rpimRes[i].R(),
               rmergeRes[i].result().count,
               Nint(imeanRes[i].Mean()), Nint(sqrt(rmsDRes[i].Mean())),
               Nint(avSdRes[i].Mean()), Iovsd,
               mnIsdRes[i].Mean(), frcbias,
               meanChiSqRes[i].Mean(), meanChiSqRes2[i].Mean());
    // Totals
    Rfull += rmergeResFull[i];
    Rmeas += rmeasRes[i];
    Rpim  += rpimRes[i];
    Imean += imeanRes[i];
    rmsD += rmsDRes[i];
    avSd += avSdRes[i];
    mnIsd += mnIsdRes[i];
    bias += biasRes[i];
    biasI += biasIRes[i];
    chisq += meanChiSqRes[i];
    chisq2 += meanChiSqRes2[i];
  }

  table.CloseTable();
  float frcbias = 0.0;
  if (biasI.Count() > 0) {
    frcbias = bias.Mean()/biasI.Mean();
  }
  // Output table to log file and XML
  output.logTab(0,LOGFILE, "\n"+table.format());
  output.logTab(0,LXML,table.XMLformat());

  fmt = "Overall:          "+fmt;
  Iovsd = 0.0;
  if (rmsD.Mean() > 0.0) {
    Iovsd = Imean.Mean()/sqrt(rmsD.Mean());
  }
  output.logTabPrintf(0,LOGFILE,fmt.c_str(),
                      Rcum.R(), Rfull.R(), Rcum.R(),
                      Rmeas.R(), Rpim.R(), Rcum.result().count,
                      Nint(Imean.Mean()), Nint(sqrt(rmsD.Mean())), Nint(avSd.Mean()),
                      Iovsd, mnIsd.Mean(), frcbias, chisq.Mean(), chisq2.Mean());
  output.logTab(0,LOGFILE,table.RawLabels());
  // Store things in summary object
  summarystatistics.StoreMnIsd(mnIsd.Mean(), mnIsdRes[0].Mean(),
                               mnIsdRes[ResRange.Nbins()-1].Mean());
  // Resolution "limit" from Mn(I/sd)
  summarystatistics.StoreMnIsigresolimit
    (ResolutionLimit(mnIsdRes, ResRange, MinimumIoverSigma, ResolutionLimit::NONE));
  // Resolution "limit" from Mn(I/sd) for I/sig > 2 (special for Frank von Delft)
  summarystatistics.StoreMnIsigresolimit2
    (ResolutionLimit(mnIsdRes, ResRange, 2.0, ResolutionLimit::NONE));
  // Filtered mean Chi^2
  summarystatistics.StoreMnChisqc(chisq2.Mean(), meanChiSqRes2[0].Mean(),
                                  meanChiSqRes2[ResRange.Nbins()-1].Mean());
}
//--------------------------------------------------------------
void PrintDeviationsByRun(const PxdName& dataset_pxd,
                          const ResoRange& ResRange,
                          const std::vector<Run>& runlist,
                          std::vector<std::vector<Rfactor> >& rmergeRun,
                          phaser_io::Output& output)
{
  // Count runs in this dataset: these have non-null R-factors
  ASSERT (runlist.size() == rmergeRun.size());
  int nruns = runlist.size();
  if (nruns == 1) return;  // don't bother for one run
  int nrd = 0;  // non-zero runs
  std::vector<bool> runpresent(nruns, false);
  for (int irun=0;irun<nruns;++irun) {
    ASSERT (int(rmergeRun[irun].size()) == ResRange.Nbins());
    bool nonzero = false;
    for (int i=0;i<ResRange.Nbins();++i) {
      if (rmergeRun[irun][i].result().count > 0) {
        nonzero = true;
        break;
      }
    }
    if (nonzero) {
      nrd++;
      runpresent[irun] =  true;
    }
  } // end loop runs
  if (nrd <= 1) {
    return; // quit if only one run
  }

  const int MAXRUNSTOPRINT = 16;
  if (nrd > MAXRUNSTOPRINT) {
    output.logTab(0,LOGFILE,
          "\nRmeas by resolution for each run suppressed because there are too many runs\n");
    return;
  }

  Range xrange = ResRange; // x axis range to full resolution limit
  xrange.first() = 0.0;    // from 0
  Range yrange;
  // Get y range (if loggraph would accept just an xrange, wouldn't need to do this)
  for (int irun=0;irun<nruns;++irun) {
    for (int i=0;i<ResRange.Nbins();++i) {
      yrange.update(rmergeRun[irun][i].R());   // Rmerge
    }
  }

  output.logTab(0,LOGFILE,
                std::string("\nRmeas by resolution for each run")+
                "\n--------------------------------\n");


  TableGraph table(" Analysis against resolution for each run in dataset, "+
                   dataset_pxd.dname());
  table.StoreID("Graph-StatsByRun");

  TableGraphPlot graph("Rmeas v. resolution for each run");  // 1st graph

  // Graph for each run, column 2 + 4...
  std::vector<int> cln(nrd+1, 2);
  int k = 1;
  for (int irun=0;irun<nruns;++irun) {
    if (runpresent[irun]) {
      cln[k] = 3 + k;
      graph.AddLine(TableGraphPlotline(2,cln[k])); // column numbers for x,y
      k++;
    }
  }
  graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
  graph.SetYaxis("", true, yrange);  // y axis from 0 to maximum
  table.AddGraph(graph);

  std::vector<std::string> collabels;
  collabels.push_back("N");         // 1
  collabels.push_back("1/d^2");     // 2
  collabels.push_back("Dmid");      // 3
  std::vector<bool> Zero(3,false);
  std::string fmt = ""; // excluding 1st 3 columns
  for (int irun=0;irun<nruns;++irun) {
    if (runpresent[irun]) {
      std::string label = StringUtil::Strip("Run"+
                                            StringUtil::itos(runlist[irun].RunNumber(),4));
      collabels.push_back(label);
      Zero.push_back(true);
      fmt += "%7.3f";
    }
  }
  fmt += "\n";
  table.StoreColumnFields(collabels, Zero, "%3d%8.4f%7.2f"+fmt);
  int nc = collabels.size();
  std::vector<Rfactor> RmeasR(nrd); // overall values

  int n=1;
  std::vector<double> Rfacs(nrd);  // Rs for printing

  for (int i=0;i<ResRange.Nbins();++i) {
    int k = 0;
    for (int irun=0;irun<nruns;++irun) {
      if (runpresent[irun]) {
        Rfacs[k] = rmergeRun[irun][i].R();
        RmeasR[k] += rmergeRun[irun][i]; // overall for each run
        k++;
      }
    }
    table.Line(Rfacs, 3,
               n++, ResRange.middle(i), ResRange.middleA(i));
  } // end loop resolution bins
  table.CloseTable();

  // Output table to log file and XML
  output.logTab(0,LOGFILE, "\n"+table.format());
  output.logTab(0,LXML,table.XMLformat());

  // Overall values
  std::string s = "Overall:          ";
  for (int k=0;k<nrd;++k) {
    s += StringUtil::ftos(RmeasR[k].R(),7,3);
  }
  output.logTab(0,LOGFILE, s);
  output.logTab(0,LOGFILE,table.RawLabels());
}
//--------------------------------------------------------------
void PrintDeviationsByResolutionOv(const PxdName& dataset_pxd,
                                   const ResoRange& ResRange,
                                   const std::vector<Rfactor>& rmergeRes,
                                   const std::vector<Rfactor>& rmeasRes,
                                   const std::vector<Rfactor>& rpimRes,
                                   const std::vector<Rfactor>& rmergeResOv,
                                   const std::vector<Rfactor>& rmeasResOv,
                                   const std::vector<Rfactor>& rpimResOv,
                                   SummaryStatistics& summarystatistics,
                                   phaser_io::Output& output)
// Statistics against overall mean I+- (only if ANOMALOUS ON) & separate
{
  output.logTab(0,LOGFILE,
                std::string("\n\nBy 4sinTheta/Lambda^2 bins (statistics with and without anomalous)\n")+
                                "==================================================================\n");
  output.logTab(0,LOGFILE,
                "\nStatistics labelled 'Ov' are relative to the overall mean I+/-, ignoring anomalous");
  output.logTab(0,LOGFILE,
                "Other statistics are with either I+ or I- sets, for acentrics, ie with anomalous\n\n");
  TableGraph table(" Analysis against resolution, with & without anomalous (Ov), "+dataset_pxd.dname());
  table.StoreID("Graph-StatsAllVsResolution");

  Range xrange = ResRange; // x axis range to full resolution limit
  xrange.first() = 0.0;    // from 0
  Range yrange;   // for each graph (just one here)
  // Get y ranges for each graph (if loggraph would accept just an xrange, wouldn't nned to do this)
  for (int i=0;i<ResRange.Nbins();++i) {
    yrange.update(rmergeRes[i].R());
    yrange.update(rmergeResOv[i].R());
    yrange.update(rmeasRes[i].R());
    yrange.update(rmeasResOv[i].R());
    yrange.update(rpimRes[i].R());
    yrange.update(rpimResOv[i].R());
  }
  yrange.first() = 0.0;  // from 0

  TableGraphPlot graph("Rmerge, Rmeas, Rpim v Resolution");
  std::string description = "within I+/I- sets, and over all data (Ov values)";
  graph.SetDescription(description);
  graph.AddLine(TableGraphPlotline(2,4));
  graph.AddLine(TableGraphPlotline(2,5));
  graph.AddLine(TableGraphPlotline(2,8));
  graph.AddLine(TableGraphPlotline(2,9));
  graph.AddLine(TableGraphPlotline(2,10));
  graph.AddLine(TableGraphPlotline(2,11));
  graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
  graph.SetYaxis("", true, yrange);  // y axis from 0 to maximum
  table.AddGraph(graph);

  std::vector<std::string> collabels;
  collabels.push_back("N");         // 1
  collabels.push_back("1/d^2");     // 2
  collabels.push_back("Dmid");      // 3
  collabels.push_back("Rmrg");      // 4
  collabels.push_back("RmrgOv");    // 5
  collabels.push_back("Rcum");      // 6
  collabels.push_back("RcumOv");    // 7
  collabels.push_back("Rmeas");     // 8
  collabels.push_back("RmeasOv");   // 9
  collabels.push_back("Rpim");      // 10
  collabels.push_back("RpimOv");    // 11
  collabels.push_back("Nmeas");     // 12
  bool z[] = {false, false, false, true, true, true, true, true, true, true, true, false};
  std::vector<bool> Zero(z, z+12);
  std::string fmt = "%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%9d\n"; // excluding 1st 3 columns
  table.StoreColumnFields(collabels, Zero, "%3d%8.4f%7.2f"+fmt);
  int nc = collabels.size();
  Rfactor Rcum, RcumOv;  // cumulative R
  Rfactor Rmeas, Rpim, RmeasOv, RpimOv;

  int n=1;
  for (int i=0;i<ResRange.Nbins();++i) {
    Rcum += rmergeRes[i];
    RcumOv += rmergeResOv[i];
    table.Line(nc, n++, ResRange.middle(i), ResRange.middleA(i), // 1,2,3
               rmergeRes[i].R(), rmergeResOv[i].R(), // 4,5
               Rcum.R(), RcumOv.R(), // 6,7
               rmeasRes[i].R(), rmeasResOv[i].R(),  //8,9
               rpimRes[i].R(), rpimResOv[i].R(),    //10,11
               rmergeResOv[i].result().count);     //12
    // Totals
    Rmeas += rmeasRes[i];
    Rpim  += rpimRes[i];
    RmeasOv += rmeasResOv[i];
    RpimOv  += rpimResOv[i];
  }
  table.CloseTable();

  // Output table to log file and XML
  output.logTab(0,LOGFILE, "\n"+table.format());
  output.logTab(0,LXML,table.XMLformat());

  fmt = "Overall:          "+fmt;
  output.logTabPrintf(0,LOGFILE,fmt.c_str(),
                      Rcum.R(), RcumOv.R(),
                      Rcum.R(), RcumOv.R(),
                      Rmeas.R(), RmeasOv.R(),
                      Rpim.R(), RpimOv.R(),
                      RcumOv.result().count);
  output.logTab(0,LOGFILE,table.RawLabels());
  // Store things in summary object
  summarystatistics.StoreRmergeReso(Rcum, rmergeRes[0], rmergeRes[ResRange.Nbins()-1]);
  summarystatistics.StoreRmeasReso(Rmeas, rmeasRes[0], rmeasRes[ResRange.Nbins()-1]);
  summarystatistics.StoreRpimReso(Rpim, rpimRes[0], rpimRes[ResRange.Nbins()-1]);
  summarystatistics.StoreRmergeResoOv(RcumOv, rmergeResOv[0], rmergeResOv[ResRange.Nbins()-1]);
  summarystatistics.StoreRmeasResoOv(RmeasOv, rmeasResOv[0], rmeasResOv[ResRange.Nbins()-1]);
  summarystatistics.StoreRpimResoOv(RpimOv, rpimResOv[0], rpimResOv[ResRange.Nbins()-1]);
}
//--------------------------------------------------------------
void PrintDeviationsByIntensity(const PxdName& dataset_pxd,
                                const IntensityBin& Irange, const bool& Anom,
                                const std::vector<Rfactor>& rmergeInt,
                                const std::vector<Rfactor>& rmeasInt,
                                const std::vector<Rfactor>& rpimInt,
                                const std::vector<MeanSD>&  imeanInt,
                                const std::vector<MeanSD>&  rmsDInt,
                                const std::vector<MeanSD>&  avSdInt,
                                const std::vector<MeanSD>&  mnIsdInt,
                                const std::vector<MeanSD>&  biasInt,
                                const std::vector<MeanSD>&  biasIInt,
                                phaser_io::Output& output)
{
  output.logTab(0,LOGFILE,
                std::string("\n\nBy intensity bins\n")+
                "=================\n");
  if (Anom) {
    output.logTab(0,LOGFILE,
          "All statistics in this table are within I+ or I- sets (anomalous on)");
  } else {
    output.logTab(0,LOGFILE,
          "All statistics in this table are relative to the overall mean I+/- (anomalous off)");
  }

  TableGraph table(" Analysis against intensity, "+dataset_pxd.dname());
  table.StoreID("Graph-StatsVsIntensity");

  TableGraphPlot graph("Rmerge v Intensity");
  std::string description = "The important values are in the top bin: ";
  description += " Rmerge: "+StringUtil::ftos(rmergeInt.back().R(), 7,3);
  description += " Rmeas: "+StringUtil::ftos(rmeasInt.back().R(), 7,3);
  description += " Rpim: "+StringUtil::ftos(rpimInt.back().R(), 7,3);
  graph.SetDescription(description);
  graph.AddLine(TableGraphPlotline(1,2));
  graph.AddLine(TableGraphPlotline(1,4));
  graph.AddLine(TableGraphPlotline(1,5));
  graph.SetYaxis("", true);  // y axis from 0 to maximum
  table.AddGraph(graph);

  std::vector<std::string> collabels;
  collabels.push_back("Imax");      // 1
  collabels.push_back("Rmrg");      // 2
  collabels.push_back("Rcum");      // 3
  collabels.push_back("Rmeas");     // 4
  collabels.push_back("Rpim");      // 5
  collabels.push_back("Nmeas");     // 6
  collabels.push_back("AvI");       // 7
  collabels.push_back("RMSdev");    // 8
  collabels.push_back("sd");        // 9
  collabels.push_back("I/RMS");     // 10
  collabels.push_back("Mn(I/sd)");  // 11
  collabels.push_back("FrcBias");   // 12

  bool z[] = {false, true, true, true, true, false, true, true, true, true, true, true};
  std::vector<bool> Zero(z, z+12);
  std::string fmt = "%7.3f%7.3f%7.3f%7.3f%9d%9d%7d%7d%7.1f%9.1f%9.3f\n"; // excluding 1st column
  table.StoreColumnFields(collabels, Zero, "%10.0f"+fmt);
  int nc = collabels.size();
  Rfactor Rcum;  // cumulative R
  Rfactor Rmeas, Rpim;
  MeanSD Imean, rmsD, avSd, mnIsd, bias, biasI;

  double Iovsd = 0.0;
  for (int i=0;i<Irange.NumberBins();++i) {
    Rcum += rmergeInt[i];
    float frcbias = 0.0;
    if (biasIInt[i].Count() > 0) {
      frcbias = biasInt[i].Mean()/biasIInt[i].Mean();
    }
    Iovsd = 0.0;
    if (rmsDInt[i].Mean() > 0.0) {
      Iovsd = imeanInt[i].Mean()/sqrt(rmsDInt[i].Mean());
    }
    table.Line(nc, Irange.mean(i),
               rmergeInt[i].R(), Rcum.R(),
               rmeasInt[i].R(), rpimInt[i].R(),
               rmergeInt[i].result().count,
               Nint(imeanInt[i].Mean()), Nint(sqrt(rmsDInt[i].Mean())),
               Nint(avSdInt[i].Mean()), Iovsd,
               mnIsdInt[i].Mean(), frcbias);
    // Totals
    Rmeas += rmeasInt[i];
    Rpim  += rpimInt[i];
    Imean += imeanInt[i];
    rmsD += rmsDInt[i];
    avSd += avSdInt[i];
    mnIsd += mnIsdInt[i];
    bias += biasInt[i];
    biasI += biasIInt[i];
  }
  table.CloseTable();
  // Output table to log file and XML
  output.logTab(0,LOGFILE, "\n"+table.format());
  output.logTab(0,LXML,table.XMLformat());

  float frcbias = 0.0;
  if (biasI.Count() > 0) {
    frcbias = bias.Mean()/biasI.Mean();
  }
  fmt = "Overall:  "+fmt+"\n";
  Iovsd = 0.0;
  if (rmsD.Mean() > 0.0) {
    Iovsd = Imean.Mean()/sqrt(rmsD.Mean());
  }
  output.logTabPrintf(0,LOGFILE,fmt.c_str(),
                   Rcum.R(), Rcum.R(),
                      Rmeas.R(), Rpim.R(), Rcum.result().count,
                   Nint(Imean.Mean()), Nint(sqrt(rmsD.Mean())), Nint(avSd.Mean()),
                   Iovsd, mnIsd.Mean(), frcbias);
  output.logTab(0,LOGFILE,table.RawLabels());
}
//--------------------------------------------------------------
void PrintCompletenessMultiplicity(const PxdName& dataset_pxd,
                                   const ResoRange& ResRange,
                                   const hkl_symmetry& symmetry,
                                   const Scell& cell,
                                   std::vector<int>& NumRef,
                                   std::vector<int>& NumObs,
                                   std::vector<int>& NumRefSphere,
                                   std::vector<int>& NumCentric,
                                   std::vector<int>& NumACentric,
                                   std::vector<int>& NumAnom,
                                   std::vector<int>& NumAnomSphere,
                                   std::vector<double>& SNumAnomPairs,
                                   SummaryStatistics& summarystatistics,
                                   phaser_io::Output& output)
//    std::vector<int> NumRef        // Number of unique reflections
//    std::vector<int> NumObs        // Number of observations
//    std::vector<int> NumRefSphere  // Number unique in sphere
//    std::vector<int> NumCentric    // Number unique centric
//    std::vector<int> NumACentric   // Number unique acentric
//    std::vector<int> NumAnom       // number unique anomalous
//    std::vector<int> NumAnomSphere // number unique in sphere
//    std::vector<double> SNumAnomPairs // anomalous pairs
{
  output.logTab(0,LOGFILE,
std::string("\n\nCompleteness and multiplicity, including reflections measured only once\n")+
                "=======================================================================\n\n\n"+
                " %poss is completeness in the shell, C%poss in cumulative to that resolution\n"+
                " The anomalous completeness values (AnomCmpl) are the percentage of possible anomalous "+
                "differences measured\n"+
                " AnomFrc is the % of measured acentric reflections for which an anomalous "+
                "difference has been measured\n"+
                " Anomalous multiplicity AnoMlt is calculated for reflections with both I+ and I- measured,"+
                " and is defined as:\n"+
                " Sum{[Min(n+, n-) + Dn/(Dn+1)]}/NanomMeasured, where n+, n- are the number of measurements "+
                " of I+, I-, Dn = |n+ - n-|\n\n");

  // Get number of reflections in each resolution bin in complete sphere
  int nbins = ResRange.Nbins();
  std::vector<int> Nrefres(nbins,0);  // number
  std::vector<int> Nrefacen(nbins,0);  // number of acentrics
  // Get numbers in each resolution shell
  NumberComplete(ResRange, symmetry, cell, Nrefres, Nrefacen);

  Range xrange = ResRange; // x axis range to full resolution limit
  xrange.first() = 0.0;    // from 0
  std::vector<Range> yranges(2);   // for each graph
  // Get y ranges for each graph (if loggraph would accept just an xrange, wouldn't need to do this)
  for (int i=0;i<ResRange.Nbins();++i) {
    yranges[1].update(FractionN(1.0, NumObs[i], NumRef[i]));
    yranges[1].update(FractionN(1.0, SNumAnomPairs[i], NumACentric[i]));
  }
  yranges[0].first() = 0.0;  // maximum completeness
  yranges[0].last() = 100.0;  // maximum completeness

  TableGraph table(" Completeness & multiplicity v. resolution, "+dataset_pxd.dname());
  table.StoreID("Graph-CompletenessVsResolution");

  TableGraphPlot graph("Completeness v Resolution ");
  std::string description = "%poss, completeness in shell; C%poss, cumulative completeness. ";
  description +=  "\nAnomalous completeness (AnomCmpl) is the percentage of possible anomalous differences measured. ";
  description +=  "\nAnomFrc is the % of measured acentric reflections for which an anomalous difference has been measured";
  graph.SetDescription(description);
  graph.AddLine(TableGraphPlotline(2,7));
  graph.AddLine(TableGraphPlotline(2,8));
  graph.AddLine(TableGraphPlotline(2,10));
  graph.AddLine(TableGraphPlotline(2,11));
  graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
  graph.SetYaxis("", true, yranges[0]);  // y axis from 0 to maximum
  table.AddGraph(graph);

  graph.init("Multiplicity v Resolution");
  description = "Total multiplicity and anomalous multiplicity";
  graph.SetDescription(description);
  graph.AddLine(TableGraphPlotline(2,9));
  graph.AddLine(TableGraphPlotline(2,12));
  graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
  graph.SetYaxis("", true, yranges[1]);  // y axis from 0 to maximum
  table.AddGraph(graph);

  std::vector<std::string> collabels;
  collabels.push_back("N");         // 1
  collabels.push_back("1/d^2");     // 2
  collabels.push_back("Dmid");      // 3
  collabels.push_back("Nmeas");     // 4
  collabels.push_back("Nref");      // 5
  collabels.push_back("Ncent");     // 6
  collabels.push_back("%poss");     // 7
  collabels.push_back("C%poss");    // 8
  collabels.push_back("Mlplct");    // 9
  collabels.push_back("AnoCmp");    // 10
  collabels.push_back("AnoFrc");    // 11
  collabels.push_back("AnoMlt");    // 12
  bool z[] = {false, false, false, false, false, false, true, true, true, true, true, true};
  std::vector<bool> Zero(z, z+12);
  std::string fmt = "%9d%9d%9d%7.1f%7.1f%7.1f%9.1f%7.1f%7.1f\n"; // excluding 1st 3 columns
  table.StoreColumnFields(collabels, Zero, "%3d%8.4f%7.2f"+fmt);
  int nc = collabels.size();

  int Nobs=0;
  int Nref=0;
  int NrefSphere=0;
  int Ncen=0;
  int Nacen=0;
  int Nrefcmpl=0;
  int Nanom=0;
  int NanomSphere=0;
  int Nanomref=0;
  float SNumanompairs = 0.0;
  float poss;
  float cumposs;
  float anomcmpl;
  float anomfrc;
  float anommult;
  float poss0, possN, mult0, multN, anomcmpl0, anomcmplN, anommult0, anommultN;

  int n=1;
  for (int i=0;i<ResRange.Nbins();++i) {
    poss = FractionN(100., NumRefSphere[i], Nrefres[i]);
    if (i==0) poss0 = poss;
    if (i==ResRange.Nbins()-1) possN = poss;
    Nobs +=  NumObs[i];
    NrefSphere +=  NumRefSphere[i];
    Nref +=  NumRef[i];
    Nrefcmpl += Nrefres[i];
    cumposs = FractionN(100., NrefSphere, Nrefcmpl);  // cumulative completeness
    anomcmpl = FractionN(100., NumAnomSphere[i], Nrefacen[i]);
    anomfrc = FractionN(100., NumAnom[i], NumACentric[i]);
    anommult = FractionN(1.0, SNumAnomPairs[i], NumAnom[i]);
    table.Line(nc, n++, ResRange.middle(i), ResRange.middleA(i),
               NumObs[i], NumRef[i], NumCentric[i], poss, cumposs,
               FractionN(1.0, NumObs[i], NumRef[i]),
               anomcmpl, anomfrc, anommult);
    // Totals
    Ncen += NumCentric[i];
    Nanom += NumAnom[i];
    NanomSphere += NumAnomSphere[i];
    Nanomref += Nrefacen[i];
    SNumanompairs += SNumAnomPairs[i];
    Nacen += NumACentric[i];
    // Inner & outer
    if (i==0) {
      poss0 = poss;
      mult0 = FractionN(1.0, NumObs[i], NumRef[i]);
      anomcmpl0 = anomcmpl;
      anommult0 = FractionN(1.0, SNumAnomPairs[i], NumAnom[i]);
    } else if (i==ResRange.Nbins()-1) {
      possN = poss;
      multN = FractionN(1.0, NumObs[i], NumRef[i]);
      anomcmplN = anomcmpl;
      anommultN = FractionN(1.0, SNumAnomPairs[i], NumAnom[i]);
    }

  }
  table.CloseTable();

  // Output table to log file and XML
  output.logTab(0,LOGFILE, "\n"+table.format());
  output.logTab(0,LXML,table.XMLformat());

  anomcmpl = FractionN(100., NanomSphere, Nanomref);
  std::string leader = "Overall:          ";
  fmt = leader+fmt;
  output.logTabPrintf(0,LOGFILE,fmt.c_str(),
                      Nobs, Nref, Ncen, cumposs, cumposs, FractionN(1.0, Nobs, Nref),
                      anomcmpl, FractionN(100., Nanom, Nacen), FractionN(1.0, SNumanompairs, Nanom));
  int lab1 = leader.size(); // 1st character in column to use labels
  leader.assign(lab1,' ');
  output.logTab(0,LOGFILE, leader+table.RawLabels().substr(lab1, table.RawLabels().size()-lab1));
  output.logTabPrintf(0,LOGFILE,"\n");
  // Store things in summary object
  summarystatistics.StoreNumbers(Nobs, NumObs[0], NumObs[ResRange.Nbins()-1],
                                 Nref, NumRef[0], NumRef[ResRange.Nbins()-1]);
  summarystatistics.StoreCmplMult(cumposs, poss0, possN,
                                  FractionN(1.0, Nobs, Nref), mult0, multN);
  summarystatistics.StoreAnomCmplMult(anomcmpl, anomcmpl0, anomcmplN,
                                      FractionN(1.0, SNumanompairs, Nacen), anommult0, anommultN);
}
//--------------------------------------------------------------
void PrintHalfDatasetCorrelations(const PxdName& dataset_pxd,
                                  const ResoRange& ResRange,
                                  const HalfDataset& halfDatasetScores,
                                  SummaryStatistics& summarystatistics,
                                  phaser_io::Output& output)
// Anomstatus    first estimate of whether there is anomalous in any dataset
{
  ASSERT (ResRange.Nbins() == halfDatasetScores.NresBin());

  output.logTab(0,LOGFILE,
    std::string("\n\nCorrelation coefficients for anomalous differences & Imean between random half-datasets (CC1/2)\n")+
                "===============================================================================================");
  output.logTab(0,LOGFILE,std::string("\n")+
   " CC(1/2) values (for Imean and anomalous differences) are calculated by splitting the data randomly in half\n"+
                " CC(1/2)v for Imean is calculated from variances, see Assmann, Brehm & Diederichs(2016),J.Appl.Cryst.49,1021-1028)");
  output.logTab(0,LOGFILE,
 std::string(" The RMS Correlation Ratio (RCR) is calculated from a scatter plot of pairs of DeltaI(anom)\n")+
 " from the two subsets (halves) by comparing the RMS value (excluding extremes) projected on the line \n"+
 " with slope = 1 ('correlation') with the RMS value perpendicular to this ('error').\n"+
 " This ratio will be > 1 if there is a significant anomalous signal\n");

  output.logTab(0,LOGFILE,
    "\n Rsplit = (1/Sqrt(2)) Sum (|I1 - I2|)/0.5*Sum(I1 + I2) where I1,I2 are the half-dataset intensities as for CC(1/2)");
  output.logTab(0,LOGFILE,
                " Note that internal R-factors of any sort are deprecated as metrics for assessment of effective resolution\n");

  ResolutionLimit overallresolimit = halfDatasetScores.OverallResoLimit();
  ResolutionLimit anomresolimit = halfDatasetScores.AnomalousResoLimit();

  double highres = overallresolimit.HighResolution();
  double anomhighres = anomresolimit.HighResolution();

  bool curvefitted = (overallresolimit.Fittype() != ResolutionLimit::NONE);
  // should both be fitted or not fitted
  ASSERT (curvefitted == (anomresolimit.Fittype() != ResolutionLimit::NONE));

  bool validfitCC = overallresolimit.valid() && overallresolimit.sufficientData();
  bool validfitanom = anomresolimit.valid() && anomresolimit.sufficientData();
  int ncurves = 0; // number of fitted curves
  if (validfitCC) {ncurves++;}
  if (validfitanom) {ncurves++;}

  output.logTab(0,LOGFILE,
                std::string("\n Estimates of maximum resolution for intensities and")+
                " anomalous differences,\n  based on the point at which CC(1/2) falls"+
                " below a threshold");
  std::string s = "\n Curve fitting as suggested by Ed Pozharski to a tanh function\n";
  s +="  of the form (1/2)(1 - tanh(z)) where z = (s - d0)/r,\n";
  s += "    s = 1/d^2,";
  s += " d0 is the value of s at the half-falloff value, and r controls the steepness of falloff";
  output.logTab(0,LOGFILE, s);
  output.logTab(0,LOGFILE,
                "\nEstimate of resolution limit for intensities:\n"+overallresolimit.format(false));
  output.logTab(0,LOGFILE,
                "\nEstimate of resolution limit for significant anomalous differences:\n"+
                anomresolimit.format(true));

  TableGraph table(" Correlations CC(1/2) within dataset, "+dataset_pxd.dname());
  table.StoreID("Graph-CChalf");

  Range xrange = ResRange; // x axis range to full resolution limit
  xrange.first() = 0.0;    // from 0
  std::vector<Range> yranges(3);   // for each graph
  // Get y ranges for each graph
  for (int i=0;i<ResRange.Nbins();++i) {
    yranges[0].update(halfDatasetScores.CCanom(i).result().val);
    yranges[0].update(halfDatasetScores.CC_Imean(i).result().val);
    yranges[1].update(halfDatasetScores.RMScorrelRatio(i));
    //    yranges[1].update(halfDatasetScores.RMScorrelRatioCen(i));
    yranges[2].update(halfDatasetScores.rsplit(i).result().val);
  }

  yranges[0].first() = std::min(yranges[0].first(), 0.0);  // CC
  yranges[0].last() = 1.0;  // CC

  std::string title = " CC(1/2) v resolution, max resolution "+
    StringUtil::Strip(StringUtil::ftos(highres,8,2))+
    ", anom "+ StringUtil::Strip(StringUtil::ftos(anomhighres,8,2));

  TableGraphPlot graph(title);
  double cchalflimit = overallresolimit.Limit();
  double ccanomlimit = anomresolimit.Limit();
  std::string description = "Resolution estimate: "+overallresolimit.formatbrief(false)+
    +"\n"+"Anomalous resolution: "+ anomresolimit.formatbrief(true);
  graph.SetDescription(description);
  graph.AddLine(TableGraphPlotline(2,4));
  graph.AddLine(TableGraphPlotline(2,7));
  graph.AddLine(TableGraphPlotline(2,9));
  for (int i=0;i<ncurves;++i) { // ncurves may == 0
    graph.AddLine(TableGraphPlotline(2,11+i,"","",0));
  }
  graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
  graph.SetYaxis("", false, yranges[0]);  // y axis from 0 to maximum
  table.AddGraph(graph);

  graph.init(" RMS anomalous correlation ratio ");
  description = "RMS correlation ratio > 1.0 indicates significant anomalous differences";
  graph.SetDescription(description);
  graph.AddLine(TableGraphPlotline(2,6));
  graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
  graph.SetYaxis("", true, yranges[1]);  // y axis from 0 to maximum
  table.AddGraph(graph);

  graph.init(" Rsplit ");
  description
    = "Rsplit is the R-factor equivalent of CC(1/2), but less useful";
  graph.SetDescription(description);
  graph.AddLine(TableGraphPlotline(2,9));
  graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
  graph.SetYaxis("", true, yranges[2]);  // y axis from 0 to maximum
  table.AddGraph(graph);

  std::vector<std::string> collabels;
  collabels.push_back("N");         // 1
  collabels.push_back("1/d^2");     // 2
  collabels.push_back("Dmid");      // 3
  collabels.push_back("CCanom");    // 4
  collabels.push_back("Nanom");     // 5
  collabels.push_back("RCRanom");   // 6
  collabels.push_back("CC1/2");     // 7
  collabels.push_back("NCC1/2");    // 8
  collabels.push_back("CC1/2v");    // 9
  collabels.push_back("Rsplit");    // 10
  if (validfitCC) {collabels.push_back("CCfit");}      // 11
  if (validfitanom) {collabels.push_back("CCanomfit");}  // 12
  int nc = collabels.size();
  bool z[] = {false, false, false, true, false, true, true,
              false, true, true, true, true};
  std::vector<bool> Zero(z, z+nc);
  std::string fmt1 = "%7.3f%9d   %7.3f %7.3f%9d %7.3f %8.3f"; // excluding 1st 3 columns
  std::string fmt = fmt1;
  for (int i=0;i<ncurves;++i) { // ncurves may == 0
    fmt += "%10.3f";
  }
  fmt += "\n";
  fmt1 += "\n";
  table.StoreColumnFields(collabels, Zero, "%3d%8.4f%7.2f"+fmt);

  int n=1;
  for (int i=0;i<ResRange.Nbins();++i) {
    std::vector<double> vcc;
    if (curvefitted) {
      double s = ResRange.middle(i);
      if (validfitCC) {
        vcc.push_back(overallresolimit.fitvalue(s));
      }
      if (validfitanom) {
        vcc.push_back(anomresolimit.fitvalue(s));
      }
    }
    table.Line(vcc, nc-ncurves, n++, ResRange.middle(i), ResRange.middleA(i),
               halfDatasetScores.CCanom(i).result().val,
               halfDatasetScores.CCanom(i).result().count,
               halfDatasetScores.RMScorrelRatio(i),
               halfDatasetScores.CC_Imean(i).result().val,
               halfDatasetScores.CC_Imean(i).result().count,
               halfDatasetScores.CC_half(i),
               halfDatasetScores.rsplit(i).result().val);
  }
  table.CloseTable();

  // Output table to log file and XML
  output.logTab(0,LOGFILE, "\n"+table.format());
  output.logTab(0,LXML,table.XMLformat());

  // Totals
  std::string leader = "Overall:          ";
  fmt = leader+fmt1;
  output.logTabPrintf(0,LOGFILE,fmt.c_str(),
                      halfDatasetScores.CCanom().result().val,
                      halfDatasetScores.CCanom().result().count,
                      halfDatasetScores.RMScorrelRatio(),
                      halfDatasetScores.CC_Imean().result().val,
                      halfDatasetScores.CC_Imean().result().count,
                      halfDatasetScores.CC_half(),
                      halfDatasetScores.rsplit().result().val,
                      0.0, 0.0);
  int lab1 = leader.size(); // 1st character in column to use labels
  leader.assign(lab1,' ');
  output.logTab(0,LOGFILE, leader+table.RawLabels().substr(lab1, table.RawLabels().size()-lab1));
  output.logTabPrintf(0,LOGFILE,"\n");
  // Store things in summary object
  summarystatistics.StoreImeanCorrel(halfDatasetScores.CC_Imean().result().val,
                                     halfDatasetScores.CC_Imean(0).result().val,
                                     halfDatasetScores.CC_Imean(ResRange.Nbins()-1).result().val);
  summarystatistics.StoreAnomCorrel(halfDatasetScores.CCanom().result().val,
                                    halfDatasetScores.CCanom(0).result().val,
                                    halfDatasetScores.CCanom(ResRange.Nbins()-1).result().val);
  summarystatistics.StoreAnomRCR(halfDatasetScores.RMScorrelRatio(),
                                 halfDatasetScores.RMScorrelRatio(0),
                                 halfDatasetScores.RMScorrelRatio(ResRange.Nbins()-1));
  summarystatistics.StoreHalfdatsetCCresolimit
    (halfDatasetScores.OverallResoLimit());
  summarystatistics.StoreHalfdatsetCCanomresolimit
    (halfDatasetScores.AnomalousResoLimit());
//@//  // all by shell, overall
//@//  summarystatistics.StoreHalfdatasetThings
//@//    (halfDatasetScores.allCC_Imean(), halfDatasetScores.CC_Imean(),
//@//     halfDatasetScores.allCCanom(), halfDatasetScores.CCanom(),
//@//     halfDatasetScores.allrsplit(), halfDatasetScores.rsplit());
}
//--------------------------------------------------------------
void PrintAnisotropyAnalysis(const PxdName& dataset_pxd,
                             const ResoRange& ResRange,
                             const HalfDataset& halfDatasetScores,
                             const std::vector<std::vector<MeanSD> >& mnIsdResCone,
                             const AnisotropicAnalysis& anisoanal,
                             const double& MinimumIoverSigma,
                             SummaryStatistics& summarystatistics,
                             phaser_io::Output& output)
{
  ASSERT (ResRange.Nbins() == halfDatasetScores.NresBin());

  if (anisoanal.IsCubic()) {
    output.logTab(0,LOGFILE,"\nNo anisotropy analysis for cubic system\n");
    return;
  }

  std::string s = std::string("\n\nAnalysis of anisotropy of data\n")+
    "==============================\n\n"+
    "Mn(I/sd) and half-dataset correlation coefficients CC(1/2) are analysed by resolution\n";
  output.logTab(0,LOGFILE,s+anisoanal.formattype());

  std::vector<std::string> axlabels(3);
  std::vector<std::string> axesformat = anisoanal.Axesformat();
  bool isplane = anisoanal.IsPlane();
  if (anisoanal.AreGeneralAxes()) { // General directions
    output.logTab(1,LOGFILE,
                  "Principal axes:");
    axlabels[0] = "d1";
    axlabels[1] = "d2";
    axlabels[2] = "d3";
    for (int i=0;i<3;++i) {
      output.logTab(2,LOGFILE,axlabels[i]+": "+axesformat[i]);
    }
  } else if (isplane) {
    output.logTab(1,LOGFILE,"Directions for analysis:");
    axlabels[0] = "d12";
    axlabels[1] = "";
    axlabels[2] = "d3";
    output.logTab(2,LOGFILE,"Plane "+axlabels[0]+": "+axesformat[0]);
    output.logTab(2,LOGFILE,axlabels[2]+": "+axesformat[2]);
  } else {
    output.logTab(1,LOGFILE,
                  "Principal axes are along a*, b*, c*");
    axlabels[0] = "a*";
    axlabels[1] = "b*";
    axlabels[2] = "c*";
  }

  DVect3  eigenvalues = anisoanal.EigenValuesOrth();
  output.logTabPrintf(0,LOGFILE,
                      "\nEigenvalues of [B](orth) along principal axes : %8.3f %8.3f %8.3f\n",
                      eigenvalues[0], eigenvalues[1], eigenvalues[2]);
  output.logTabPrintf(0,LOGFILE,
      "Difference between maximum and minimum anisotropic B (= 8 pi^2 U) %7.1f\n",
                      anisoanal.BfactorDifference());

  std::vector<ResolutionLimit> resolutionlimits = halfDatasetScores.AnisoResoLimits();
  bool curvefitted = false;
  int nax = resolutionlimits.size();
  int ii = 1;
  if (isplane) {ii = 2;}
  std::vector<bool> validfit(nax, false);  // valid fit for each axis
  int ncurvefits = 0;
  for (int i=0; i<nax; i++) {
    validfit[i] = resolutionlimits.at(i).valid() && resolutionlimits.at(i).sufficientData();
    if (isplane && (i == 1)) {
      // for isplane, skip axis 1 (2nd axis)
      validfit[i] = false;
    } else {
      if (validfit[i]) {ncurvefits++;} // number of plotted curve fits
    }
    if (validfit[i]) {curvefitted = true;}
  }

  std::vector<double> highres(nax);
  for (int i=0; i<nax; i++) {
    highres[i] = resolutionlimits[i].HighResolution();
  }

  std::string smaxres = "Estimated maximum resolution limits, ";
  if (isplane) {
    smaxres += axlabels[0]+":"+StringUtil::ftos(highres[0],6,2)+", "+
      axlabels[2]+":"+StringUtil::ftos(highres[2],6,2);
  } else {
    smaxres += axlabels[0]+":"+StringUtil::ftos(highres[0],6,2)+", "+
      axlabels[1]+":"+StringUtil::ftos(highres[1],6,2)+", "+
      axlabels[2]+":"+StringUtil::ftos(highres[2],6,2);
  }
  output.logTab(0,LOGFILE, "\n "+smaxres);

  if (curvefitted) {
    output.logTab(0,LOGFILE,
                  " Columns 'CCft' are values from curve-fitting as for overall analysis");
  }

  Range xrange = ResRange; // x axis range to full resolution limit
  xrange.first() = 0.0;    // from 0
  std::vector<Range> yranges(2);   // for each graph
  for (int i=0;i<ResRange.Nbins();++i) {
    for (int k=0;k<3;++k) {
      if (!(isplane && k == 1)) {
        yranges[0].update(halfDatasetScores.CCaniso(k,i).result().val);
      }
    }
    yranges[1].update(mnIsdResCone[0][i].Mean());
    yranges[1].update(mnIsdResCone[1][i].Mean());
    yranges[1].update(mnIsdResCone[2][i].Mean());
  }
  yranges[0].first() = std::min(yranges[0].first(), 0.0); // CC
  yranges[0].last() = 1.0; // CC

  TableGraph table(" Anisotropy analysis of CC(1/2) and I/sd, "+dataset_pxd.dname());
  table.StoreID("Graph-Anisotropy");

  std::vector<int> cln;
  int j1;
  if (isplane) { // only two directions if plane
    int c[] = {2,4,5};
    cln.assign(c,c+3);
    j1 = 10;
  } else { // 3 directions
    int c[] = {2,4,5,6};
    cln.assign(c,c+4);
    j1 = 13;
  }
  for (int i=0;i<ncurvefits;++i) {
    cln.push_back(i+j1);
  }

  TableGraphPlot graph(" Anisotropic CC(1/2) v resolution");
  std::string description = "CC(1/2) in cones along principle directions. " + smaxres; // maximum resolution limits
  graph.SetDescription(description);

  for (size_t i=1; i<cln.size(); i++) { // from 1
    if (cln[i] >= j1) {
      graph.AddLine(TableGraphPlotline(2,cln[i],"","",0,false,"Solid",1)); // curve fit, no symbols
    } else {
      graph.AddLine(TableGraphPlotline(2,cln[i]));
    }
  }
  graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
  graph.SetYaxis("", false, yranges[0]);  // y axis from 0 to maximum
  table.AddGraph(graph);

  if (isplane) {
    int c2[] = {2,6,7};
    cln.assign(c2,c2+3);
  } else {
    int c2[] = {2,7,8,9};
    cln.assign(c2,c2+4);
  }
  graph.init(" Anisotropic Mn(I/sd) v resolution");
  description = "Mn(I/sd) in cones along principle directions"; // maximum resolution limits
  graph.SetDescription(description);

  for (size_t i=1; i<cln.size(); i++) { // from 1
    graph.AddLine(TableGraphPlotline(2,cln[i]));
  }
  graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
  graph.SetYaxis("", true, yranges[1]);  // y axis from 0 to maximum
  table.AddGraph(graph);

  if (isplane) {
    int c3[] = {2,8,9};
    cln.assign(c3,c3+3);
  } else {
    int c3[] = {2,10,11,12};
    cln.assign(c3,c3+4);
  }
  graph.init(" Projected CC(1/2) v resolution");
  description = "CC(1/2) projected on principle directions"; // maximum resolution limits
  graph.SetDescription(description);
  for (size_t i=1; i<cln.size(); i++) { // from 1
    graph.AddLine(TableGraphPlotline(2,cln[i]));
  }
  graph.SetXaxis("", true, xrange);  // x axis is 1/d^2
  graph.SetYaxis("", true, yranges[0]);  // y axis from 0 to maximum
  table.AddGraph(graph);

  std::vector<std::string> collabels;
  collabels.push_back("N");         // 1
  collabels.push_back("1/d^2");     // 2
  collabels.push_back("Dmid");      // 3
  collabels.push_back("CC:"+axlabels[0]);     // 4
  if (!isplane) collabels.push_back("CC:"+axlabels[1]);     // 5
  collabels.push_back("CC:"+axlabels[2]);     // 6
  collabels.push_back("I/sd:"+axlabels[0]);     // 7
  if (!isplane) collabels.push_back("I/sd:"+axlabels[1]);     // 8
  collabels.push_back("I/sd:"+axlabels[2]);     // 9
  collabels.push_back("CCp1");     // 10
  if (!isplane) collabels.push_back("CCp2");     // 11
  collabels.push_back("CCp3");     // 12
  if (curvefitted) {
    for (int i=0;i<nax;++i) {
      if (validfit[i]) {
        collabels.push_back("CCft:"+axlabels[i]);     // 13, 14, 15
      }
    }
  }

  int nc = collabels.size();

  bool z[] = {false, false, false, true, true, true, true, true, true, true, true,
              true, true, true, true};
  std::vector<bool> Zero(z, z+nc);
  std::string fmt1 = " %8.3f%8.3f%8.3f%9.2f%9.2f%9.2f%8.3f%8.3f%8.3f"; // excluding 1st 3 columns
  if (isplane) {
    fmt1 = " %8.3f%8.3f%9.2f%9.2f%8.3f%8.3f"; // excluding 1st 3 columns
  }
  std::string fmt = fmt1;
  if (curvefitted) {
    for (int i=0;i<ncurvefits;++i) {
      fmt += "%9.3f";
    }  }
  fmt += "\n";

  table.StoreColumnFields(collabels, Zero, "%3d%8.4f%7.2f"+fmt);

  std::vector<MeanSD> mnIsd(3);
  int n=1;
  int ncol = nc - ncurvefits;
  for (int i=0;i<ResRange.Nbins();++i) {
    double s = ResRange.middle(i);
    std::vector<double> vcc;
    if (isplane) {
      if (curvefitted) {
        // fitted values
        if (validfit[0]) {
          vcc.push_back(resolutionlimits[0].fitvalue(s));
        }
        if (validfit[2]) {
          vcc.push_back(resolutionlimits[2].fitvalue(s));
        }
      }
      if (vcc.size() == 0) {
        table.Line(nc, n++, ResRange.middle(i), ResRange.middleA(i),
                   halfDatasetScores.CCaniso(0,i).result().val,
                   halfDatasetScores.CCaniso(2,i).result().val,
                   mnIsdResCone[0][i].Mean(),
                   mnIsdResCone[2][i].Mean(),
                   halfDatasetScores.CCanisoProjection(0,i).result().val,
                   halfDatasetScores.CCanisoProjection(2,i).result().val);
      } else {
        // note vcc may be null
        table.Line(vcc, ncol, n++, ResRange.middle(i), ResRange.middleA(i),
                   halfDatasetScores.CCaniso(0,i).result().val,
                   halfDatasetScores.CCaniso(2,i).result().val,
                   mnIsdResCone[0][i].Mean(),
                   mnIsdResCone[2][i].Mean(),
                   halfDatasetScores.CCanisoProjection(0,i).result().val,
                   halfDatasetScores.CCanisoProjection(2,i).result().val);
      }
    } else {
      if (curvefitted) {
        // fitted values
        for (int j=0;j<3;++j) {
          if (validfit[j]) {
            vcc.push_back(resolutionlimits[j].fitvalue(s));
          }
        }

      }
      if (vcc.size() == 0) {
        table.Line(nc, n++, ResRange.middle(i), ResRange.middleA(i),
                   halfDatasetScores.CCaniso(0,i).result().val,
                   halfDatasetScores.CCaniso(1,i).result().val,
                   halfDatasetScores.CCaniso(2,i).result().val,
                   mnIsdResCone[0][i].Mean(),
                   mnIsdResCone[1][i].Mean(),
                   mnIsdResCone[2][i].Mean(),
                   halfDatasetScores.CCanisoProjection(0,i).result().val,
                   halfDatasetScores.CCanisoProjection(1,i).result().val,
                   halfDatasetScores.CCanisoProjection(2,i).result().val);
      } else {
        table.Line(vcc, ncol, n++, ResRange.middle(i), ResRange.middleA(i),
                   halfDatasetScores.CCaniso(0,i).result().val,
                   halfDatasetScores.CCaniso(1,i).result().val,
                   halfDatasetScores.CCaniso(2,i).result().val,
                   mnIsdResCone[0][i].Mean(),
                   mnIsdResCone[1][i].Mean(),
                   mnIsdResCone[2][i].Mean(),
                   halfDatasetScores.CCanisoProjection(0,i).result().val,
                   halfDatasetScores.CCanisoProjection(1,i).result().val,
                   halfDatasetScores.CCanisoProjection(2,i).result().val);
      }
    }

    //^
    //    std::cout
    //      << "  " <<      halfDatasetScores.CCaniso(0,i).result().count
    //      << "  " <<      halfDatasetScores.CCaniso(1,i).result().count
    //      << "  " <<      halfDatasetScores.CCaniso(2,i).result().count
    //      << "  " <<      mnIsdResCone[0][i].Count()
    //      << "  " <<      mnIsdResCone[1][i].Count()
    //      << "  " <<      mnIsdResCone[2][i].Count() << "\n";
    //^-
    mnIsd[0] += mnIsdResCone[0][i];
    mnIsd[1] += mnIsdResCone[1][i];
    mnIsd[2] += mnIsdResCone[2][i];
  }
  table.CloseTable();

  // Output table to log file and XML
  output.logTab(0,LOGFILE, "\n"+table.format());
  output.logTab(0,LXML,table.XMLformat());

  // Totals
  std::string leader = "Overall:          ";
  fmt = leader+fmt1;

  if (isplane) {
    s = FormatOutput::logTabPrintf(0,fmt.c_str(),
                                   halfDatasetScores.CCaniso(0).result().val,
                                   halfDatasetScores.CCaniso(2).result().val,
                                   mnIsd[0].Mean(),
                                   mnIsd[2].Mean(),
                                   halfDatasetScores.CCanisoProjection(0).result().val,
                                   halfDatasetScores.CCanisoProjection(2).result().val);
  } else {
    s = FormatOutput::logTabPrintf(0,fmt.c_str(),
                                   halfDatasetScores.CCaniso(0).result().val,
                                   halfDatasetScores.CCaniso(1).result().val,
                                   halfDatasetScores.CCaniso(2).result().val,
                                   mnIsd[0].Mean(),
                                   mnIsd[1].Mean(),
                                   mnIsd[2].Mean(),
                                   halfDatasetScores.CCanisoProjection(0).result().val,
                                   halfDatasetScores.CCanisoProjection(1).result().val,
                                   halfDatasetScores.CCanisoProjection(2).result().val);
  }
  std::string s2;
  if (ncurvefits > 0) {
    for (int i=0;i<ncurvefits;++i) {
      s2 += "    0.0  ";
    }
  }

  output.logTab(0,LOGFILE, s+s2+"\n");

  int lab1 = leader.size(); // 1st character in column to use labels
  leader.assign(lab1,' ');
  output.logTab(0,LOGFILE, leader+
                table.RawLabels().substr(lab1, table.RawLabels().size()-lab1));
  output.logTabPrintf(0,LOGFILE,"\n");
  // Store resolution limit stuff in summary statistics
  summarystatistics.StoreHalfdatsetCCAnisoresolimit
    (halfDatasetScores.AnisoResoLimits());
  std::vector<ResolutionLimit> reslimisig(3);
  for (int i=0;i<3;++i) {
    if (!(isplane && (i==1))) {
      // find resolution "limit" where Mn(I/sd) falls below MinimumIoverSigma
      reslimisig[i].init(mnIsdResCone[i], ResRange, MinimumIoverSigma,
                         ResolutionLimit::NONE);
    }
  }
  summarystatistics.StoreMnIsigAnisoresolimit(reslimisig);
}
//--------------------------------------------------------------
void PrintUnmergedHeaderStuff(const scala::hkl_unmerge_list& hkl_list,
                              phaser_io::Output& output,
                              const int& verbose)
// Optional summary printing
{
  if (hkl_list.MultiLattice()) {
    output.logTab(0,LOGFILE," ");
    output.logTabPrintf(1,LOGFILE,
                        "Multiple lattice data, number of lattices %2d\n",
                        (hkl_list.NumberofMainLattices()));
  }
  if (verbose > 0)  {
    if (verbose == 1 || !hkl_list.IsReady()) {
      output.logTab(0,LOGFILE,
                    "      ResolutionRange    NobsParts  Nbatches  Ndatasets");
      output.logTabPrintf(0,LOGFILE,
                          "     %8.2f %6.2f  %10d%7d%10d\n",
                          hkl_list.ResRange().ResLow(), hkl_list.ResRange().ResHigh(),
                          hkl_list.num_accepted_parts(), hkl_list.num_batches(),
                          hkl_list.num_datasets());

    } else {
      //        output.logTabPrintf(0,LOGFILE,
      //                            "\nSummary of reflection list\n");
      output.logTabPrintf(0,LOGFILE,
                          "\n   Resolution range accepted: %8.2f    %8.2f\n",
                          hkl_list.ResLimRange().ResLow(),
                          hkl_list.ResLimRange().ResHigh());

      output.logTabPrintf(0,LOGFILE,
                          "\n   Number of reflections  =    %10d\n",
                          hkl_list.num_reflections_valid());
      output.logTabPrintf(0,LOGFILE,
                          "   Number of observations =    %10d\n",
                          hkl_list.num_observations());
      output.logTabPrintf(0,LOGFILE,
                          "   Number of parts        =    %10d\n",
                          hkl_list.num_accepted_parts());
      output.logTabPrintf(0,LOGFILE,
                          "   Number of batches      =    %10d\n",
                          hkl_list.num_batches());
      output.logTabPrintf(0,LOGFILE,
                          "   Number of datasets     =    %10d\n",
                          hkl_list.num_datasets());
    }
    if (verbose <= 3) {
      int ndatasets = hkl_list.num_datasets();
      std::vector<Batch> batches = hkl_list.Batches();
      int nbatches = batches.size();
      // 1st, last batch
      std::vector<std::pair<int,int> >  rejectedbatches;
      // Make list of rejected batch ranges
      int b1 = -1;
      for (int ib=0;ib<nbatches;ib++) {
        if (!batches[ib].Accepted()) {
          if (b1 < 0) {b1 = ib;}
        } else {
          if (b1 >= 0) {
            // store range
            if (b1 == ib-1) {
              // just one
              rejectedbatches.push_back
                (std::pair<int,int>(batches[b1].num(),0));
            } else {
              rejectedbatches.push_back
                (std::pair<int,int>(batches[b1].num(),batches[ib-1].num()));
            }
            b1 = -1;
          }
        }
      } // end loop batches
      if (b1 >= 0) {
        int ib = nbatches;
        // store range
        if (b1 == ib-1) {
          // just one
          rejectedbatches.push_back
            (std::pair<int,int>(batches[b1].num(),0));
        } else {
          rejectedbatches.push_back
            (std::pair<int,int>(batches[b1].num(),batches[ib-1].num()));
        }
      }
      int nrb = rejectedbatches.size();
      std::vector<Dataset> datasets = hkl_list.AllDatasets();
      std::vector<Run> runlist = hkl_list.RunList();
      for (int k=0; k<ndatasets; k++) {
        if (!datasets[k].accepted()) {
          output.logTab(0,LOGFILE,
                        "\nRejected dataset "+datasets[k].pxdname().format());
        } else {
          output.logTab(0,LOGFILE,"\n"+datasets[k].format());
          //      output.logTab(0,LOGFILE,"\n"+datasets[k].formatPrint());
          //      output.logTab(3,LOGFILE,"Cell: "+datasets[k].cell().formatPrint());
          //      output.logTabPrintf(3,LOGFILE,
          //                          "Wavelength %8.5f A\n", datasets[k].wavelength());
          if (verbose == 3) {
            for (size_t i=0;i<runlist.size();i++) {
              if (runlist[i].DatasetIndex() == k) {
                output.logTab(0,LOGFILE,runlist[i].formatPrintBrief());
                std::string rejlist;
                for (int k=0;k<nrb;k++) {
                  if (runlist[i].IsInList(rejectedbatches[k].first)) {
                    // Rejected batches in this run
                    if (rejlist.size() > 0) rejlist += ", ";
                    if (rejectedbatches[k].second == 0) {
                      rejlist +=
                        StringUtil::Strip(clipper::String(rejectedbatches[k].first,6));
                    } else {
                      rejlist +=
                        StringUtil::Strip(clipper::String(rejectedbatches[k].first,6)+
                                          "-"+clipper::String(rejectedbatches[k].second,6));
                    }
                  }
                }
                if (rejlist.size() > 0) {
                  output.logTab(3,LOGFILE,"Excluded batches: "+rejlist);}
              }
            }
          }
        }
      }  // end loop datasets
    } else if (verbose > 3) {
      std::vector<Dataset> datasets = hkl_list.AllDatasets();
      for (size_t k=0; k<datasets.size(); k++) {
        if (!datasets[k].accepted()) {
          output.logTab(0,LOGFILE,
                        "\nRejected dataset "+datasets[k].pxdname().format());
        } else {
          output.logTab(0,LOGFILE,
                        datasets[k].formatPrint());
        }
      }
    }
    output.logTabPrintf(0,LOGFILE, "\n   Average unit cell: ");
    output.logTab(0,LOGFILE,hkl_list.Cell().formatPrint());
    output.logTabPrintf(1,LOGFILE,"");

    if (verbose > 3) {
      std::vector<Run> runlist = hkl_list.RunList();
      std::vector<Dataset> datasets = hkl_list.AllDatasets();
      for (size_t i=0;i<runlist.size();i++) {
        output.logTab(0,LOGFILE,runlist[i].formatPrint(datasets));}
    }
  }

  // XML things
  if (verbose > 2) {
    std::string hklstream = "HKLIN";
    // Is there more than one unique file number?
    std::vector<Run> runlist = hkl_list.RunList();
    int fn = runlist[0].FileNumber();
    bool OneFile = true;
    if (runlist.size() > 1) {
      for (size_t i=1;i<runlist.size();i++) {
        if (runlist[i].FileNumber() != fn) {OneFile = false;}
      }
    }
    output.logTab(0,LXML,"<ReflectionData>");
    float resmax = hkl_list.ResLimRange().ResHigh();
    output.logTab(1,LXML,
                  StringUtil::MakeXMLtag("ResolutionHigh",
                                         resmax,8,2));
    output.logTab(1,LXML,
                  StringUtil::MakeXMLtag("NumberReflections",
                  hkl_list.num_reflections_valid(),10));
    output.logTab(1,LXML,
                  StringUtil::MakeXMLtag("NumberObservations",
                                         hkl_list.num_observations(),10));
    output.logTab(1,LXML,
                  StringUtil::MakeXMLtag("NumberParts",
                                         hkl_list.num_accepted_parts(),10));
    int numberoflattices = hkl_list.NumberofMainLattices();
    output.logTab(1,LXML,
                  StringUtil::MakeXMLtag("NumberLattices",
                                         numberoflattices));
    output.logTab(1,LXML,
                  StringUtil::MakeXMLtag("NumberBatches",
                                         hkl_list.num_batches(),10));
    output.logTab(1,LXML,
                  StringUtil::MakeXMLtag("NumberDatasets",
                                         hkl_list.num_datasets(),10));

    int ndatasets = hkl_list.num_datasets();
    std::vector<Dataset> datasets = hkl_list.AllDatasets();
    for (int k=0; k<ndatasets; k++) {
      output.logTabPrintf(1,LXML, "<Dataset  name=\"%s\">\n",
                          datasets[k].formatNames().c_str());
      for (size_t i=0;i<runlist.size();i++) {
        if (runlist[i].DatasetIndex() == k) {
          output.logTab(2,LXML,
              StringUtil::MakeXMLtag("Wavelength", datasets[k].wavelength()));
          output.logTabPrintf(2,LXML,"<Run> <number> %3d </number>\n",
                              runlist[i].RunNumber());
          output.logTab(2,LXML,
                        StringUtil::MakeXMLtag("BatchRange",
                        StringUtil::itos(runlist[i].BatchRange().first,8)+
                        StringUtil::itos(runlist[i].BatchRange().second,8)));

          output.logTab(2,LXML,StringUtil::MakeXMLtag("BatchOffset",
                                                      runlist[i].BatchNumberOffset(),8));
          if (!OneFile) {
            hklstream =
              StringUtil::Strip("HKLIN"+clipper::String(runlist[i].FileNumber()));
          }
          output.logTab(2,LXML,StringUtil::MakeXMLtag("FileStream", hklstream));
          if (numberoflattices > 1) {
            output.logTab(2,LXML,
                          StringUtil::MakeXMLtag("Lattice",runlist[i].LatticeNumber()));
          }
          output.logTabPrintf(2,LXML,"</Run>\n",i+1);
        }
      }
      output.logTabPrintf(1,LXML, "</Dataset>\n");
    }
    output.logTab(0,LXML,"</ReflectionData>");
  }
}
//--------------------------------------------------------------
void PrintPartialCounts(const scala::hkl_unmerge_list& hkl_list,
                        const all_controls& controls,
                        phaser_io::Output& output)
// Print counts of partials rejected etc, also to XML
{
  if (hkl_list.num_observations_part() > 0) {
      output.logTab(0,LOGFILE,controls.partials.format());
      if (hkl_list.num_observations_scaled() > 0) {
        output.logTabPrintf(0,LOGFILE,"Number of scaled partials =  %10d\n",
                            hkl_list.num_observations_scaled());
        output.logTab(0,LXML,StringUtil::MakeXMLtag("NumberScaledPartials",
                                            hkl_list.num_observations_scaled()));
      }
      output.logTabPrintf(0,LOGFILE,"\n");
      // Rejects
      output.logTabPrintf(0,LOGFILE,
                          "%8d  partial sets rejected with total fraction too small\n",
                          hkl_list.num_observations_rejected_FracTooSmall());
      output.logTab(0,LXML,StringUtil::MakeXMLtag("NumberPartialsTooSmall",
                          hkl_list.num_observations_rejected_FracTooSmall()));
      output.logTabPrintf(0,LOGFILE,
                          "%8d  partial sets rejected with total fraction too large\n",
                          hkl_list.num_observations_rejected_FracTooLarge());
      output.logTab(0,LXML,StringUtil::MakeXMLtag("NumberPartialsTooLarge",
                          hkl_list.num_observations_rejected_FracTooLarge()));
      output.logTabPrintf(0,LOGFILE,
                          "%8d  partial sets rejected with gaps\n",
                          hkl_list.num_observations_rejected_Gap());
      output.logTab(0,LXML,StringUtil::MakeXMLtag("NumberPartialsGap",
                          hkl_list.num_observations_rejected_Gap()));
    }
}
//--------------------------------------------------------------
