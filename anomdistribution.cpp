//
// anomdistribution.cpp 
//

#include "anomdistribution.hh"
#include "selectedobservations.hh"
#include "tablegraph.hh"
#include "jiffy.hh"
#include "string_util.hh"

using phaser_io::LOGFILE;
using phaser_io::itos;


namespace scala {
  AllAnomDistributions::AllAnomDistributions(const hkl_unmerge_list& hkl_list,
					     const SDmodel& SDM,
					     const all_controls& controls,
					     const ResoRange& ResRange,
					     const Normalise& NormRes)
  // Analyse distribution of anomalous differences to get estimate
  // of maximum likely values, for all datasets
  //
  {
    reflection this_refl;

    // Resolution ranges
    int nresbin =  ResRange.Nbins();
    resrange = ResRange;
    // Number of datasets
    ndatasets = hkl_list.num_datasets();
    bool correlAnom;

    // rms DelAnom by dataset & resolution
    anomdistributions.resize(ndatasets);
    // Set number of resolution bins
    for (int id=0;id<ndatasets;++id)
      {anomdistributions[id].SetNresbin(nresbin);}
    //  datasets
    std::vector<Xdataset> xdatasets = hkl_list.AllXdatasets();
    pxdnames.resize(ndatasets);
    wavelengths.resize(ndatasets);
    dnames.resize(ndatasets);
    for (int id=0;id<ndatasets;id++) {
      pxdnames[id] = xdatasets[id].pxdname();
      wavelengths[id] = xdatasets[id].wavelength();
      dnames[id] = pxdnames[id].dname();
    }

    basedataset = controls.datasetcontrol.BaseDataset();
    int ncorrel = ndatasets*(ndatasets - 1)/2;    // number of anomalous cc
    int ndispcc = (ndatasets-2)*(ndatasets - 1)/2; // number of dispersive cc
    if (ncorrel > 0) {
      cca.resize(ncorrel);  // CC between anomalous differences
      ccadtsindex.resize(ncorrel);
      for (int i=0;i<ncorrel;++i) {cca[i].resize(nresbin);}
      int k=0;
      for (int j=0;j<ndatasets-1;++j) {
	for (int i=j+1;i<ndatasets;++i) {
	  // for each CC, store pair of dataset indices
	  ccadtsindex[k++] = std::pair<int,int>(j,i);
	}}
      ASSERT (k == ncorrel);
      if (ndispcc > 0) {
	ccd.resize(ndispcc); // CC between dispersive differences
	ccddtsindex.resize(ndispcc);
	for (int i=0;i<ndispcc;++i) {ccd[i].resize(nresbin);}
	// Choose base dataset as the one with the shortest wavelength
	// unless specified
	if (basedataset < 0) {
	  float wvl = 10000.;
	  for (int id=0;id<ndatasets;id++) {
	    if (wavelengths[id] < wvl) {
	      wvl = wavelengths[id];
	      basedataset = id;
	    }
	  }
	  if (wvl < 0.001) basedataset = 0;
	}
	// Store dataset index pairs for each CC
	k = 0;
	for (int j=0;j<ndatasets-1;++j) {
	  if (j != basedataset) {
	    for (int i=j+1;i<ndatasets;++i) {
	      if (i != basedataset) {
		ccddtsindex[k++] = std::pair<int,int>(j,i);
	      }
	    }
	  }
	}
	ASSERT (k == ndispcc);
      }
    }

    float danom;
    std::vector<float> danomdts(ndatasets); // DelAnom for each dataset
    std::vector<float> Imeandts(ndatasets); // <I> for each dataset, for dispersive values
    SelectedObservations Selobs;
    int nacc = 0;

    while (hkl_list.next_reflection(this_refl) >= 0)  {
      Rtype invresolsq = this_refl.invresolsq();
      bool Centric = hkl_list.symmetry().is_centric(this_refl.hkl());
      // Resolution bin
      int mres = ResRange.bin(invresolsq);
      //  Apply current SD correction to reflection (all observations)
      SDM.CorrectReflection(this_refl);
      correlAnom = true;  // flag for at least 2 I+ & 2 I- observations

      // loop datasets
      nacc = 0;
      danomdts.assign(ndatasets,0.0);
      for (int id=0;id<ndatasets;id++) {
	Selobs.init(this_refl, id, ALL);
	if (Selobs.Number() > 0) {
	  Imeandts[id] = Selobs.Average().I();
	  if (!Centric) {
	    Selobs.init(this_refl, id, IPLUS);
	    if (Selobs.Number() > 0) {
	      IsigI Iplus = Selobs.Average();
	      if (Selobs.Number() <= 1) {correlAnom = false;}
	      Selobs.init(this_refl, id, IMINUS);
	      if (Selobs.Number() > 0) {
		IsigI Iminus = Selobs.Average();
		if (Selobs.Number() <= 1) {correlAnom = false;}
		// DelAnom
		float sig = Iplus.sigI()*Iplus.sigI() + Iminus.sigI()*Iminus.sigI();
		if (sig > 0.0) {
		  // Store delAnom, & count reflections used for
		  // half-dataset correlations (ie with n+ & n- > 1)
		  danom = (Iplus.I()-Iminus.I());
		  anomdistributions[id].Add(mres, danom, correlAnom);
		  danomdts[id] = danom;
		  nacc++;
		}
	      }
	    }
	  }
	}
      } // end loop datasets
      if (nacc > 1) {
	// correlations
	AddCorrelations(danomdts, Imeandts, mres);
      }
    } // end loop reflections
  }
  // ------------------------------------------------------------
  void AnomDistribution::SetNresbin(const int& Nresbin)
  // Store number of resolution bins & clear arrays
  {
    nresbin = Nresbin;
    rmsDelAnom.assign(nresbin, MeanSD());
    nDelAnom.assign(nresbin,0);
  }
  // ------------------------------------------------------------
  void AnomDistribution::Add(const int& mres,
	   const float& delAnom, const bool& correlAnom)
  // Store delAnom, & count reflections used for
  // half-dataset correlations (ie with n+ & n- > 1, correlAnom true)
  //  for resolution range mres
  {
    rmsDelAnom[mres].Add(delAnom);
    // count reflections used for half-dataset correlations
    if (correlAnom) nDelAnom[mres]++;
  }
  // ------------------------------------------------------------
  void AllAnomDistributions::AddCorrelations
  (const std::vector<float>& danomdts, const std::vector<float>& Imeandts,
   const int& mres)
  // Add in to correlation sums, resolution bin mres
  {
    ASSERT (int(danomdts.size()) == ndatasets);
    // Anomalous differences
    int k = 0;
    for (int j=0;j<ndatasets-1;++j) {
      for (int i=j+1;i<ndatasets;++i) {
	if (danomdts[i]!=0.0 && danomdts[j]!=0.0) {
	  cca[k][mres].add(danomdts[i], danomdts[j]);
	}
	k++;
      }}
    if (ndatasets > 2) { // no correlation if < 3 datasets
      k = 0;
      for (int j=0;j<ndatasets-1;++j) {
	if (j != basedataset) {
	  float dj = Imeandts[j] - Imeandts[basedataset];
	  for (int i=j+1;i<ndatasets;++i) {
	    if (i != basedataset) {
	      if (dj!=0.0 && (Imeandts[i] - Imeandts[basedataset])!=0.0) {
		ccd[k][mres].add(Imeandts[i] - Imeandts[basedataset], dj);
	      }
	      k++;
	    }
	  }
	}
      }
    }
  }
  // ------------------------------------------------------------
  void AllAnomDistributions::Print(phaser_io::Output& output) const
  // Print correlation tables
  {
    if (ndatasets < 2) return;
    output.logTab(0, LOGFILE,
  "\nCorrelation coefficients for anomalous & dispersive differences between different datasets");
    output.logTab(0, LOGFILE,
		  "==========================================================================================\n");
    output.logTab(0, LOGFILE,"\nDatasets and wavelengths:\n");
    for (int id=0;id<ndatasets;++id) {
      std::string s("    ");
      if (id == basedataset) s = "base";
      output.logTabPrintf(0, LOGFILE,"%4d %4s %8.5f %s\n",
			  id+1, s.c_str(), wavelengths[id],
			  pxdnames[id].format().c_str());
    }
    
    std::string title = ">>> Correlation of Anomalous Differences between datasets";
    std::string graphtitle = "Anom CCs v resln -";
    for (int id=0;id<ndatasets;++id) {graphtitle += " "+dnames[id];}
    std::string cl1 = "1st dataset         ";
    std::string cl2 = "2nd dataset         ";
    std::vector<correl_coeff> allcc;  // totals
    output.logTab(0, LOGFILE,
		  FormatTable(title, graphtitle, cl1, cl2, cca,
			      ccadtsindex, false, allcc));

    // Format cross-correlation table
    title = "\nOverall correlation of Anomalous Differences between datasets\n";
    title += "      (Numbers in brackets)\n\n";
    output.logTab(0, LOGFILE, CrossCorrelation(title, allcc, false));

    if (ndatasets < 3) return;

    // Dispersive differences
    title = "Correlation of Dispersive Differences between datasets";
    graphtitle = "Dispersive CCs v resln -";
    for (int id=0;id<ndatasets;++id) {graphtitle += " "+dnames[id];}
    cl1 = "1st difference      ";
    cl2 = "2nd difference      ";
    output.logTab(0, LOGFILE,
		  FormatTable(title, graphtitle, cl1, cl2, ccd,
			      ccddtsindex, true, allcc));
    // Format cross-correlation table
    title = "\nCorrelation between datasets of Dispersive Differences from base set\n";
    title += "      (Numbers in brackets)\n\n";
    output.logTab(0, LOGFILE, CrossCorrelation(title, allcc, true));
  }
  // ------------------------------------------------------------
  std::string AllAnomDistributions::FormatTable(const std::string& title,
						const std::string& graphtitle,
						const std::string& ccl1,
						const std::string& ccl2,
				const std::vector<std::vector<correl_coeff> >& cc,
				const std::vector<std::pair<int,int> >& ccidx,
				const bool& diff,
  			        std::vector<correl_coeff>& allcc) const
  // diff = true for dispersive differences
  {
    TableGraph table(title);
    std::string outstring = table.formatTitle()+"\n";
    std::vector<int> cln(1,2);
    int ng = cc.size();     // number of graphs

    for (int id=0;id<ng;++id) {
      cln.push_back(id*2+4);
    }
    outstring += table.Graph(graphtitle,"N", cln)+"\n";

    // Column labels, zero-field flags & format
    std::vector<std::string> collabels;
    std::vector<bool> Zero;
    std::string fmt = "%3d%8.4f%6.2f ";  // leader
    std::string fmtn  = "%8.3f%8d";      // each part
    collabels.push_back("N");         // 1
    collabels.push_back("1/d^2");     // 2
    collabels.push_back("Dmid");      // 3
    Zero.push_back(false);
    Zero.push_back(false);
    Zero.push_back(false);
    std::string cl1 = ccl1; // header lines
    std::string cl2 = ccl2; // header lines
    for (size_t i=0;i<cc.size();++i) { // loop anomalous CCs for dataset pairs
      int idx1 = ccidx[i].first;  // dataset indices for this CC
      int idx2 = ccidx[i].second;
      std::string cl;
      if (diff) {
	// dispersive difference
	std::string bd = itos(basedataset+1);
	cl = StringUtil::Strip(itos(idx1+1)+bd+"-"+itos(idx2+1)+bd);
	std::string dl1 = StringUtil::Strip((dnames[idx1]+"-"+dnames[basedataset]));
	std::string dl2 = StringUtil::Strip((dnames[idx2]+"-"+dnames[basedataset]));
	cl1 += StringUtil::CentreString(dl1,16);
	cl2 += StringUtil::CentreString(dl2,16);
      } else {
	// anomalous
	cl = StringUtil::Strip(itos(idx1+1)+"-"+itos(idx2+1));
	cl1 += StringUtil::CentreString(dnames[idx1],16);
	cl2 += StringUtil::CentreString(dnames[idx2],16);
      }
      collabels.push_back("CC"+cl);  // eg CC1-2  (CC)
      collabels.push_back("N"+cl);   // eg N1-2   (number)
      Zero.push_back(true);
      Zero.push_back(true);
      fmt += fmtn;
    }
    fmt += "\n";
    outstring += table.ColumnFields(collabels, Zero, fmt, false)+"\n";
    // extra header lines
    outstring += cl1+"\n"+cl2+"\n$$\n";

    allcc.assign(cc.size(), correl_coeff());  // totals
    for (int ir=0;ir<resrange.Nbins();++ir) {
      int nonzero = 0; // number of non-zero item
      for (size_t i=0;i<cc.size();++i) {
	if (cc[i][ir].result().count > 0) {nonzero++;}
      }
      if (nonzero > 0) {
	table.StartLine();
	table.AddToLine(ir+1);
	table.AddToLine(resrange.middle(ir));
	table.AddToLine(resrange.middleA(ir));
	for (size_t i=0;i<cc.size();++i) {
	  table.AddToLine(cc[i][ir].result().val);
	  table.AddToLine(cc[i][ir].result().count);
	  allcc[i] += cc[i][ir];
	}
	outstring += table.GetLine();
      }
    }  // end loop res bins
    outstring += table.CloseTable()+"\n";
    std::string line = "Overall           ";
    char buf[256];
    for (size_t i=0;i<cc.size();++i) {
      sprintf(buf, fmtn.c_str(),
	      allcc[i].result().val, allcc[i].result().count);
      line += std::string(buf);
    }
    outstring += line+"\n";
    return outstring;
  }
  // ------------------------------------------------------------
  std::string AllAnomDistributions::CrossCorrelation(const std::string& title,
						     const std::vector<correl_coeff>& allcc,
						     const bool& diff) const
  // format CC of anomalous differences between datasets as table
  {
    std::string outstring = title;
    std::string line1 = "            ";
    std::string line2 = line1;
    std::vector<std::string> ldf;
    for (int id=0;id<ndatasets;++id) { // loop datasets
      if (diff) {
	if(id != basedataset) {
	  ldf.push_back(StringUtil::Strip((dnames[id]+"-"+dnames[basedataset])));
	}
      } else {
	ldf.push_back(dnames[id]);
      }
    }
    int ncc = ldf.size();

    for (int id=1;id<ncc;++id) { // loop datasets or differences from 2nd
      line1 += StringUtil::CentreString(ldf.at(id),8);
      line2 += StringUtil::CentreString("*",8);
    }
    outstring += line1+"\n"+line2+"\n";
    int k=0;
    for (int j=0;j<ncc-1;++j) {  // loop lines
      line1 = StringUtil::LeftString(ldf[j], 10)+"*";
      line2 = "           ";
      for (int i=1;i<ncc;++i) {
	if (i < j+1) {
	  line1 += "        ";
	  line2 += "        ";
	} else {
	  line1 += StringUtil::ftos(allcc[k].result().val, 8, 3);
	  line2 += StringUtil::CentreString
	    (StringUtil::Strip("("+itos(allcc[k].result().count)+")"), 8);
	  k++;
	}
      }
      outstring += line1+"\n"+line2+"\n";
    }  // end loop lines
    return outstring+"\n";
  }
  // ------------------------------------------------------------
  // ------------------------------------------------------------
} // namespace scala
