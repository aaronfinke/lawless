// comparetoreference.cpp


#include "comparetoreference.hh"
#include "tablegraph.hh"
#include "string_util.hh"

namespace scala {
  // ------------------------------------------------------------
  CompareToReference::CompareToReference(const hkl_unmerge_list& hkl_list,
					 const ReferenceList& hklreflist,
					 const ResoRange& ResRange)
  {
    // Resolution ranges
    int nresbin =  ResRange.Nbins();
    resrange = ResRange;
    // Number of datasets
    ndatasets = hkl_list.num_datasets();
    if (ndatasets <= 1) {return;}
    if (hklreflist.IsEmpty()) {return;}

    //  datasets
    std::vector<Dataset> datasets = hkl_list.AllDatasets();
    pxdnames.resize(ndatasets);
    dnames.resize(ndatasets);
    for (int id=0;id<ndatasets;id++) {
      pxdnames[id] = datasets[id].pxdname();
      dnames[id] = datasets[id].Dname();
    }

    ccref.resize(ndatasets);  // CC
    rfref.resize(ndatasets);  // Rfactor
    for (int i=0;i<ndatasets;++i) {
      ccref[i].resize(nresbin);
      rfref[i].resize(nresbin);
    }

    reflection this_refl;
    SelectedObservations allobs;     // all I+ and I-
    IsigI avIsig;

    hkl_list.rewind();
    
    // * * * * Loop reflections
    while (hkl_list.next_reflection(this_refl) >= 0)  {
      Rtype invresolsq = this_refl.invresolsq();
      // Resolution bin
      int mres = ResRange.bin(invresolsq);

      // Reference data      
      IsigI Isref = hklreflist.Isig(this_refl.hkl());
      if (Isref.sigI() <= 0.0) {continue;}

      for (int idts=0;idts<ndatasets;++idts) {
	if (hkl_list.dataset(idts).accepted()) {
	  allobs.init(this_refl, idts, ALL);  // data for this dataset
	  avIsig = allobs.Average();  // average I, 1/variance weight
	  rfref[idts][mres].add(avIsig.I()-Isref.I(), avIsig.I(), 1.0);
	  ccref[idts][mres].add(avIsig.I(), Isref.I(), 1.0);
	}
      }
    }

  }
  // ------------------------------------------------------------
  void CompareToReference::printTable(phaser_io::Output& output) {
    output.logTab(0,LOGFILE,
		  std::string(
	"\n\nFor all datasets, agreement with reference data, analysed by resolution\n")+
            "==================================================================\n\n");

    int nresbin = resrange.Nbins();
    output.logTab(0,LOGFILE,
                std::string("Rref   is Sum(Iobs - k.Iref) / Sum(Iobs)\n")+
                "CCref  is  CC(Iobs, k.Iref) \n");

    output.logTab(0, LOGFILE,"\nDatasets:\n");
    for (int id=0;id<ndatasets;++id) {
      std::string s("    ");
      output.logTabPrintf(0, LOGFILE,"%4d     %s\n",
                          id+1, pxdnames[id].format().c_str());
    }

    int symbolsize = 1;
    TableGraph table
      ("-=- Comparison to reference data by resolution for all datasets");
    table.StoreID("Graph-DatasetRefStatsVsReso");
    TableGraphPlot graph("Rref and CCref v Resolution for all datasets");

    for (int idts=0;idts<ndatasets;++idts) {
      graph.AddLine(TableGraphPlotline(2,idts*3+4));  // default colour
      graph.AddLine(TableGraphPlotline(2,idts*3+5));  // default colour
    }
    graph.SetXaxis("", true);  // x axis is 1/d^2
    Range yrange(0.0, 1.0);   // Y from 0 to 1
    graph.SetYaxis("", true, yrange);
    table.AddGraph(graph);

    // Column labels, zero-field flags & format
    std::vector<std::string> collabels;
    std::vector<bool> Zero;
    std::string fmt = "%3d%8.4f%6.2f ";    // leader
    std::string fmtn  = "%8.3f%8.3f%8d";   // each part
    collabels.push_back("N");         // 1
    collabels.push_back("1/d^2");     // 2
    collabels.push_back("Dmid");      // 3
    Zero.push_back(false);
    Zero.push_back(false);
    Zero.push_back(false);

    for (int idts=0;idts<ndatasets;++idts) {
      std::string dn = StringUtil::itos(idts+1);
      collabels.push_back("R-"+dn);  // eg R-dn  (Rfactor)
      collabels.push_back("CC-"+dn);  // eg CC-dn (CC)
      collabels.push_back("N-"+dn);   // eg N1-2   (number)
      Zero.push_back(true);
      Zero.push_back(true);
      Zero.push_back(true);
      fmt += fmtn;
    }
    table.StoreColumnFields(collabels, Zero, fmt);

    totalrfref.assign(ndatasets, Rfactor());  // totals
    totalccref.assign(ndatasets, correl_coeff());  // totals
    for (int ir=0;ir<resrange.Nbins();++ir) {
      int nonzero = 0; // number of non-zero item
      for (size_t i=0;i<ndatasets;++i) {
        if (ccref[i][ir].result().count > 0) {nonzero++;}
      }
      if (nonzero > 0) {
        table.StartLine();
        table.AddToLine(ir+1);
        table.AddToLine(resrange.middle(ir));
        table.AddToLine(resrange.middleA(ir));
        for (size_t i=0;i<ndatasets;++i) {
          table.AddToLine(rfref[i][ir].result().val);
          table.AddToLine(ccref[i][ir].result().val);
          table.AddToLine(ccref[i][ir].result().count);
          totalrfref[i] += rfref[i][ir];
          totalccref[i] += ccref[i][ir];
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
    for (size_t i=0;i<ndatasets;++i) {
      snprintf(buf, 256, fmtn.c_str(),
	       totalrfref[i].result().val, totalrfref[i].result().count,
	       totalccref[i].result().val, totalccref[i].result().count);
      line += std::string(buf);
    }
    output.logTab(0,LOGFILE, line);

    std::string dtstag = "CompareToReference";
    output.logTab(0,LXML, "<"+dtstag+">");
    for (size_t i=0;i<ndatasets;++i) {
      output.logTabPrintf(1,LXML, "<Dataset  name=\"%s\">\n", dnames[i].c_str());
      std::string xm = StringUtil::MakeXMLtag("Rreference",
				  StringUtil::ftos(totalrfref[i].result().val));
      output.logTab(2,LXML,xm);
      xm = StringUtil::MakeXMLtag("CCreference",
				  StringUtil::ftos(totalccref[i].result().val));
      output.logTab(2,LXML,xm);
      xm = StringUtil::MakeXMLtag("Number",
				  StringUtil::itos(totalrfref[i].result().count));
      output.logTab(2,LXML,xm);
      output.logTabPrintf(1,LXML, "</Dataset>\n");
    }
    output.logTab(0,LXML, "</"+dtstag+">");
  }
}
