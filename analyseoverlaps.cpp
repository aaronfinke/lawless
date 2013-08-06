// analyseoverlaps.cpp
//
// Check all overlapped observations for self overlaps, ie overlaps with the same
// or symmetry-related spots (but not Friedel-related).
// Now the data are scaled, these overlap set can be combined into a pseudo-singleton
//


#include "analyseoverlaps.hh"
#include "string_util.hh"
#include "tablegraph.hh"

using phaser_io::LOGFILE;
using phaser_io::LXML;

namespace scala {
//--------------------------------------------------------------
  void Analyseoverlaps::init(hkl_unmerge_list& hkl_list, const int& datasetIndex,
			     const bool& verbose, phaser_io::Output& output)
  // Reclassify self-overlaps as singletons
  // Accumulate statistics on overlaps for this dataset (-1 for all)
  {
    datasetindex = datasetIndex;
    resrange = hkl_list.ResRange();
    int nresbin = resrange.Nbins();
    numbermergedres.assign(nresbin, 0);
    numbermerged = 0;
    numbernosingletonsres.assign(nresbin, 0);
    numbernoanomsingletonsres.assign(nresbin, 0);
    numberres.assign(nresbin, 0);
    numberacentricres.assign(nresbin, 0);

    reflection this_refl;
    observation this_obs;
    hkl_list.rewind();
    int isym, jsym;

    while (hkl_list.next_reflection(this_refl) >= 0)  {
      bool reflupdated = false;
      Hkl hkl_reduced = this_refl.hkl();
      bool Centric = hkl_list.symmetry().is_centric(this_refl.hkl());
      int mres = resrange.bin(this_refl.invresolsq());

      numberres[mres]++;
      if (!Centric) {numberacentricres[mres]++;}

      int nsingle = 0;
      int nsingle_minus = 0;
      int nsingle_plus = 0;

      while (this_refl.next_observation(this_obs) >= 0) {
	if (datasetindex < 0 ||
	    this_obs.datasetIndex() == datasetindex) {  // select on dataset
	  isym = this_obs.Isym();
	  int is1 = isym%2;
	  if (!this_obs.IsSingleton()) { // overlap
	    // Conditions for this_obs being truly a singleton are
	    //  (1) same reduced hkl
	    //  (2) same parity of ISYM (ie all I+ or all I-)
	    
	    // for the "principle" hkl
	    std::vector<LatticeIndexInfo> lathkl = this_obs.lathkl();
	    bool same = true;
	    for (size_t l=0; l<lathkl.size(); l++) { 
	      Hkl hkl = hkl_list.symmetry().put_in_asu(lathkl[l].hkl, jsym);
	      if (!(jsym%2 == is1) || !(hkl ==  hkl_reduced)) {
		// these are not the same
		same = false;
	      }
	    } // end loop overlaps
	    
	    if (same) { // all overlapped reflections are equivalent, so merge
	      if (verbose) {
		std::string s = "Self-overlap reset as singleton ";
		s += hkl_reduced.format();
		s += " Batch "+StringUtil::itos(this_obs.Batch(),6)+"\n";
		s += " Overlaps (lattice; hkl): "+
		  StringUtil::itos(this_obs.MainLatticeNumber(),2)+
		  "; "+this_obs.hkl_original().format();
		for (size_t l=0; l<lathkl.size(); l++) { 
		  s += " "+StringUtil::itos(lathkl[l].latnum,2)+
		    "; "+lathkl[l].hkl.format();
		}
		output.logTab(0, LOGFILE, s);
	      }
	      lathkl.clear();
	      this_obs.StoreLathkl(lathkl); // clear overlaps
	      this_refl.replace_observation(this_obs);
	      reflupdated = true;
	      // Count them, overall and by resolution
	      numbermerged++;
	      numbermergedres[mres]++;
	    }
	  } // end not singleton
	  
	  if (this_obs.IsSingleton()) { // singleton, maybe just reclassified
	    nsingle++;
	    if (!Centric) {
	      if (is1 == 0) { // I-
		nsingle_minus++;
	      } else {
		nsingle_plus++;
	      }
	    }
	  }
	} // end dataset
      } // end loop observations
      if (reflupdated) {
	hkl_list.replace_reflection(this_refl);
      }
      if (nsingle == 0) { // no singletons
	numbernosingletonsres[mres]++;
      }
      if (!Centric && ((nsingle_minus == 0) || (nsingle_plus == 0))) {
	numbernoanomsingletonsres[mres]++;  // not singletons for both I+ and I-
      }
    } // end loop reflections
  }
  //--------------------------------------------------------------
  void Analyseoverlaps::PrintOverlapTable(phaser_io::Output& output) const
  {
    output.logTab(0,LOGFILE,
		  "\nAnalysis of distribution of overlapped and singleton observations\n");
    output.logTab(0,LOGFILE,
		  "=================================================================\n\n");
    output.logTab(0,LOGFILE,
		  std::string(
      "NselfOvlp         number of self-overlaps reclassified as singletons\n")+
      "NnoSingle         number of unique reflections with no singletons\n"+
      "NnoSingleAnom     number of acentric reflections without singletons for both I+ and I-\n"+
      "\n%selfOvlp, %noSingle, %noSingleAnom   percentages of measured reflections (acentric for Anom)\n");

    TableGraph table("Overlaps and singletons vs. resolution");
    table.StoreID("Graph-Overlaps");

    TableGraphPlot graph("SelfOverlaps, #noSingletons, #noAnomSingletons");
    graph.AddLine(TableGraphPlotline(2,4)); // column numbers for x,y
    graph.AddLine(TableGraphPlotline(2,5));
    graph.AddLine(TableGraphPlotline(2,6));
    table.AddGraph(graph);

    graph.init("%SelfOverlaps, %noSingletons, %noAnomSingletons");
    graph.AddLine(TableGraphPlotline(2,7));
    graph.AddLine(TableGraphPlotline(2,8));
    graph.AddLine(TableGraphPlotline(2,9));
    table.AddGraph(graph);

    std::vector<std::string> collabels;
    collabels.push_back("N");         // 1
    collabels.push_back("1/d^2");     // 2
    collabels.push_back("Dmid");      // 3
    collabels.push_back("NselfOvlp");       // 4
    collabels.push_back("NnoSingle");       // 5 
    collabels.push_back("NnoSingleAnom");   // 6
    collabels.push_back("%selfOvlp");       // 7
    collabels.push_back("%noSingle");       // 8
    collabels.push_back("%noSingleAnom");   // 9
    int nc = collabels.size();
    // which columns should have "0.0" replaced by "-"
    bool z[] =
      {false, false, false, true, true, true, true, true, true};
    std::vector<bool> Zero(z, z+nc);
    std::string fmt = "%9d%10d%15d%11.1f%10.1f%15.1f\n"; // excluding 1st 3 columns
    // store labels, zero flags and format
    table.StoreColumnFields(collabels, Zero, "%3d%8.4f%7.2f"+fmt); 

    int numbermergedall = 0;
    int numbernosingletons = 0;
    int numbernoanomsingletons = 0;
    int numberall = 0;
    int numberacentric = 0;

    int n=1;
    for (int i=0;i<resrange.Nbins();++i) {
      float fself = 0.0;
      if (numberres[i] > 0) {
	fself = 100.*float(numbermergedres[i])/
	  float(numberres[i]);
      }
      float fnosingle = 0.0;
      if (numberres[i] > 0) {
	fnosingle = 100.*float(numbernosingletonsres[i])/
	  float(numberres[i]);
      }
      float fnosingleanom = 0.0;
      if (numberacentricres[i] > 0) {
	fnosingleanom = 100.*float(numbernoanomsingletonsres[i])/
	  float(numberacentricres[i]);
      }

      table.Line(nc, n++, resrange.middle(i),
		 resrange.middleA(i),
		 numbermergedres[i],
		 numbernosingletonsres[i],
		 numbernoanomsingletonsres[i],
		 fself, fnosingle, fnosingleanom);
      numbermergedall += numbermergedres[i];
      numbernosingletons += numbernosingletonsres[i];
      numbernoanomsingletons += numbernoanomsingletonsres[i];
      numberall += numberres[i];
      numberacentric += numberacentricres[i];
    }
    table.CloseTable();
    // Output table to log file and XML
    output.logTab(0,LOGFILE, "\n"+table.format());
    output.logTab(0,LXML,table.XMLformat());
    
    float fself = 0.0;
    if (numberall > 0) {
      fself = 100.*float(numbermergedall)/
	  float(numberall);
    }
    float fnosingle = 0.0;
    if (numberall > 0) {
      fnosingle = 100.*float(numbernosingletons)/
	float(numberall);
    }
    float fnosingleanom = 0.0;
    if (numberacentric > 0) {
      fnosingleanom = 100.*float(numbernoanomsingletons)/
	float(numberacentric);
    }

    fmt = "Overall:          "+fmt;
    output.logTabPrintf(0,LOGFILE,fmt.c_str(),
			numbermergedall, numbernosingletons,
			numbernoanomsingletons,
			fself, fnosingle, fnosingleanom);
  }
//--------------------------------------------------------------
} // namespace scala
