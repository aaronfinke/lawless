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
   const int& jrun, const int& nbatchgroup)
  {
    init(hkl_list, jrun, nbatchgroup);
  }
  //--------------------------------------------------------------
  void RadiationDamageAnalysis::init(const hkl_unmerge_list& hkl_list,
				     const int& jrun, const int& nbatchgroup)
  // jrun is run serial number
  // nbatchgroup is number of batches to group together, usually 1
  {
    irun = jrun;
    ASSERT ((irun >= 0) && (irun <  hkl_list.num_runs()));

    Run thisrun = hkl_list.RunList()[irun];
    runnum = thisrun.RunNumber();
    int nbatches = thisrun.Nbatches();

    phibinsize = 1.0; // 1 degree bins
    batchgroup = nbatchgroup;
    // reset batchgroup to something sensible (~1 degree) if not set
    if (batchgroup <= 0) {
      double delphi = thisrun.PhiRange().AbsRange()/double(nbatches);
      if (delphi > 0.0001) {
	batchgroup = Max(1, Nint(phibinsize/delphi));
      } else {
	batchgroup = 1;
      }
    }

    ntimebin = (nbatches+batchgroup-1)/batchgroup;
    int datasetIndex = thisrun.DatasetIndex();
    batch0 = thisrun.Batch0();

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
	  observations.push_back(this_obs);
	}
	if (observations.size() > 1) {
	  // we have a list of observations, now do a double loop over
	  // different observations
	  for (size_t j=0; j<observations.size()-1; j++) { 
	    int batchj = observations[j].Batch();
	    float Ij = std::abs(observations[j].kI());
	    for (size_t i=j+1; i<observations.size(); i++) { 
	      int batchi = observations[i].Batch();
	      int bintime = (Max(batchi, batchj) - batch0)/batchgroup;
	      ASSERT ((bintime >= 0) && (bintime < ntimebin));
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
  void RadiationDamageAnalysis::plot
  (const std::vector<float>& batchcompleteness,
   phaser_io::Output& output) const
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
      StringUtil::Strip(StringUtil::itos(batchgroup, 4))+
      ", ~= "+StringUtil::Strip(StringUtil::ftos(phibinsize, 7, 1))+
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
    graph.AddLine(TableGraphPlotline(2,nc,"red"));  // R
    graph.AddLine(TableGraphPlotline(2,3,"blue"));  // %completeness
    graph.SetYaxis("", true);  // Y from zero
    table.AddGraph(graph);

    graph.init("Rcp v. batch, in shells");
    for (int rbin=0;rbin<nresbin;++rbin) { // for each resolution bin
      graph.AddLine(TableGraphPlotline(2,rbin+4));   // R
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

    int n=1;
    std::vector<double> vals;
    for (int bintime=0;bintime<ntimebin;++bintime) { // loop time|dose bins
      // batch serial number in run, for completeness
      int jbatch = bintime * batchgroup;
      vals.clear();
      Rfactor rall;  // all resolution bins
      for (int rbin=0;rbin<nresbin;++rbin) { // for each resolution bin
	vals.push_back(rfactor[rbin][bintime].R());
	rall += rfactor[rbin][bintime];
      }
      vals.push_back(rall.R());

      int batch = bintime*batchgroup + batch0;
      table.Line(vals, nc0, n, batch, batchcompleteness[jbatch]);

      n++;
    }
    table.CloseTable();
    output.logTab(0,LOGFILE, "\n"+table.format());
    output.logTab(0,LXML,table.XMLformat());

  }
} // namespace scala
