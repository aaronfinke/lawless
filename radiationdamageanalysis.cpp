//
// radiationdamageanalysis.cpp
//
// Implement decay analyses following Graeme Winter, as he programmed in Chef
//

#include "radiationdamageanalysis.hh"
#include "runthings.hh"
#include "range.hh"
#include "Output.hh"
#include "tablegraph.hh"
#include "string_util.hh"

using phaser_io::LOGFILE;
using phaser_io::LXML;

namespace scala {
  //--------------------------------------------------------------
  RadiationDamageAnalysis::RadiationDamageAnalysis
  (const hkl_unmerge_list& hkl_list,
   const int& jrun,
   const Batchgroup& batchgroup)
  {
    init(hkl_list, jrun, batchgroup);
  }
  //--------------------------------------------------------------
  void RadiationDamageAnalysis::init(const hkl_unmerge_list& hkl_list,
                                     const int& jrun,
                                     const Batchgroup& batchgroup)
  // jrun is run serial number
  // nbatchgroup is number of batches to group together, sometimes 1
  {
    ASSERT (hkl_list.num_runs() == 1); // the present code assumes a single run
    irun = jrun;
    ASSERT ((irun >= 0) && (irun <  hkl_list.num_runs()));

    Run thisrun = hkl_list.RunList()[irun];
    runnum = thisrun.RunNumber();

    // all batch numbers including rejected ones
    std::vector<int> batchnumberlist = thisrun.BatchList(false);
    int nbatches = batchnumberlist.back() - batchnumberlist[0];
    if (nbatches <= 0) return;  // nothing to do

    pbatchgroup = &batchgroup;
    ntimebin = pbatchgroup->numberofgroups();
    if (ntimebin <= 0) {return;}
    int datasetIndex = thisrun.DatasetIndex();

    // resolution ranges
    ResoRange resrange = hkl_list.ResRange(); // store resolution range
    nresbin = resrange.Nbins();
    nresbin = Min(nresbin, 8);
    resrange.SetNbins(nresbin);

    rfactor.assign(nresbin, std::vector<Rfactor>(ntimebin));
    cc.assign(nresbin, std::vector<correl_coeff>(ntimebin));

    std::vector<observation> observations;
    observation this_obs;
    reflection this_refl;
    hkl_list.rewind();  // Just in case
    double w = 1.0;

    while (hkl_list.next_reflection(this_refl) >= 0) {
      // resolution bin
      int rbin = resrange.tbin(this_refl.invresolsq());
      if (rbin >= 0) { // test that reflection is in range
        observations.clear();
        while (this_refl.next_observation(this_obs) >= 0) {
          observations.push_back(this_obs); // list of accepted observations
        }
        if (observations.size() > 1) {
          // we have a list of observations, now do a double loop over
          // different observations
          for (size_t j=0; j<observations.size()-1; j++) {
            int batchj = batchIndex(observations[j].Batch());
            float Ij = std::abs(observations[j].kI());
            for (size_t i=j+1; i<observations.size(); i++) {
              int batchi = batchIndex(observations[i].Batch());
              int bintime = Min(Max(batchi, batchj), ntimebin-1);
              //              ASSERT ((bintime >= 0) && (bintime < ntimebin));
              //              if (!((bintime >= 0) && (bintime < ntimebin))){
              //                std::cout <<batchi<<" "<<batchj<<" "<<bintime<<
              //                  " "<<bintime<<"\n";
              //              } //^
              float Ii = std::abs(observations[i].kI());
              // R-factor and CC
              rfactor[rbin][bintime].add(std::abs(Ii-Ij),
                                         0.5*(Ii+Ij), w);
              cc[rbin][bintime].add(Ii, Ij, w);
            }}
        }
      }
    } // end reflection loop

    // Cumulative R-factors and CC, up to time (batch, dose) point
    for (int rbin=0;rbin<nresbin;++rbin) { // for each resolution bin
      for (int bintime=1;bintime<ntimebin;++bintime) {
        rfactor[rbin][bintime] += rfactor[rbin][bintime-1];
        cc[rbin][bintime] += cc[rbin][bintime-1];
      }}
  }
  //--------------------------------------------------------------
  int RadiationDamageAnalysis::batchIndex(const int& batchnum) const
  // get group index of this batch in the run, including rejected batches
  {
    return pbatchgroup->batchgroup(batchnum);
  }
  //--------------------------------------------------------------
  void RadiationDamageAnalysis::plot
  (const std::vector<float>& batchcompleteness,
   phaser_io::Output& output) const
  // batchcompleteness  for each batch group
  {
    std::string s =
      std::string("\nCumulative radiation damage analysis\n")+
      "====================================\n\n"+
      "At present this analysis is done only if there is a single run\n"+
      "Note that this analysis will not be useful if the multiplicity is low\n"+

"\nRcp is the cumulative pairwise residual devised by Graeme Winter for the program CHEF,\n"+
      "inspired by Diederichs Rd statistic (Acta Cryst D62, 96-101 (2005))\n";
    output.logTab(0,LOGFILE, s);

    s =
      std::string("\n")+
      "         Rcp(k) = Sum(||Ii - Ij||)/Sum(0.5*(Ii + Ij))\n"+
      " where i & j are the batch numbers (proxy for radiation dose) and"+
      " k = Max(i, j)\n"+
      " ie a pairwise R-factor up to batch k\n"+
      "\nCmPoss is cumulative completeness";
    output.logTab(0,LOGFILE, s);

    s =
      std::string("\n")+
      "Batches are binned in groups of "+
      StringUtil::Strip(StringUtil::itos(pbatchgroup->numberInGroup(), 4))+
      ", ~= "+StringUtil::Strip(StringUtil::ftos(pbatchgroup->groupWidth(), 7, 1))+
      " degrees";
    output.logTab(0,LOGFILE, s);

    std::string title = "Radiation damage analysis for run "+
      StringUtil::itos(runnum,3);

    TableGraph table(title);
    table.StoreID("Graph-RadiationDamageAnalysis");

    std::vector<std::string> collabels;
    collabels.push_back("N");           // 1
    collabels.push_back("Batch");       // 2
    collabels.push_back("CmPoss");      // 3
    for (int rbin=0;rbin<nresbin;++rbin) { // for each resolution bin
      std::string label = "R"+StringUtil::Strip(StringUtil::itos(rbin+1,4));
      collabels.push_back(label);
    }
    collabels.push_back("Rcp");
    int nc = collabels.size();
    int nc0 = 3;

    TableGraphPlot graph("Rcp v. batch");
    std::string description =
      "Rcp is cumulative pairwise residual; an increase may indicate radiation damage. ";
    description += "Cumulative completeness may help to choose a suitable cut-off point in case of damage";
    graph.SetDescription(description);

    graph.AddLine(TableGraphPlotline(2,nc,"blue","",0,false,"Solid",1));  // R
    TableGraphPlotline tgpl(2,3,"red","",0,false,"Solid",1);  // %completeness
    tgpl.SetRHaxis();
    graph.AddLine(tgpl);  // %completeness
    graph.SetYaxis("", true);  // Y from zero
    graph.SetRightYaxis("", true, Range(0.0, 1.05));  // Y from zero
    table.AddGraph(graph);

    graph.init("Rcp v. batch, in shells");
    description = "Rcp in resolution shells;";
    graph.SetDescription(description);
    for (int rbin=0;rbin<nresbin;++rbin) { // for each resolution bin
      graph.AddLine(TableGraphPlotline(2,rbin+4,"","",0,false,"Solid",1));   // R
      graph.SetYaxis("", true);  // Y from zero
    }
    table.AddGraph(graph);

    std::vector<bool> Zero(nc, true);
    Zero[0] = false;
    Zero[1] = false;
    Zero[2] = false;

    std::string fmt = "%5d%5d%9.3f";
    for (int i=0;i<nc-nc0;i++) {
      fmt += "%6.3f";
    }
    table.StoreColumnFields(collabels, Zero, fmt);

    std::vector<Batch> batches = pbatchgroup->batches();

    int n=1;
    std::vector<double> vals;
    for (int bintime=0;bintime<ntimebin;++bintime) { // loop time|dose bins == batchgroup
      // batch serial number in run, for completeness
      int jbatch = pbatchgroup->batchserial(bintime);
      vals.clear();
      Rfactor rall;  // all resolution bins
      for (int rbin=0;rbin<nresbin;++rbin) { // for each resolution bin
        vals.push_back(rfactor[rbin][bintime].R());
        rall += rfactor[rbin][bintime];
      }
      vals.push_back(rall.R());

      int batch = batches[jbatch].num();
      table.Line(vals, nc0, n, batch, batchcompleteness.at(bintime));
      n++;
    }
    table.CloseTable();
    output.logTab(0,LOGFILE, "\n"+table.format());
    output.logTab(0,LXML,table.XMLformat());
  }
} // namespace scala
