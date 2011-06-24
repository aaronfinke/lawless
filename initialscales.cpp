//
// initialscales.cpp
//
// Initial rough scaling by making average intensities equal
//


#include "initialscales.hh"
#include "foxholmes.hh"
#include <assert.h>
#define ASSERT assert


namespace scala {
// ---------------------------------------------------------
  InitialData::InitialData(const clipper::Array2d<double>& AvI)
  {
    // Store address of data array
  // Data is a 2D array(nrotranges, nresbins)
    avi = &AvI;
    np = avi->rows();  // number of parameters = number of rotation ranges
    next = -1;
  }
// ---------------------------------------------------------
  bool InitialData::ObsArray(std::vector<DPair>& obs) const
  // return one observation, I, sigma pair, length fnp
  // Return false if end of data
  {
    ASSERT (int(obs.size()) >= avi->rows());
    double sd = 1.0;  // dummy sigma
    while (++next < avi->cols()) {
      for (int i=0;i<avi->rows();++i) {
	obs[i] = DPair((*avi)(i,next), sd);
      }
      return true;
    }
    next = -1;
    return false;
  }
// ---------------------------------------------------------
  void InitialScales(hkl_unmerge_list& hkl_list, ScaleModel& AllScales,
		     const all_controls& controls,
		     phaser_io::Output& output)
  // Get initial estimates of primary scales, from making intensity
  // averages equal
  {
    output.logTab(0,LOGFILE,
		  "\n========= Initial scaling =========\n\n");

    std::vector<Run> runlist = hkl_list.RunList();   // runs
    int nruns = runlist.size();
    // true if a run has batch scales
    std::vector<bool> batch_scale_run(nruns, false);
    // Number of ranges, for batch scale = Nbatches, else number of scales-1
    std::vector<int> nranges_run(nruns,0);

    int nrotranges = 0; // total number of ranges
    std::vector<int> idxrun(nruns); // index to 1st rotation range for each run
    // Set up rotation ranges for each run
    for (int irun=0;irun<nruns;++irun) {
      if (AllScales.primary_scale(irun).IsBatchScale()) {
	// Batch scale
	batch_scale_run[irun] = true;
      }
      // Store number of rotation ranges
      //   for batch mode = Nbatches
      nranges_run[irun] = AllScales.primary_scale(irun).Nintervals();
      idxrun[irun] = nrotranges;
      nrotranges += nranges_run[irun];
      runlist[irun].PhiRange().SetNbin(nranges_run[irun]); // set up binning on phi
      //^
      //      std::cout << "Phirange run " << irun << " " << runlist[irun].PhiRange().format();
      //^-
    }
    // Only one range, bail out
    if (nrotranges <= 1) return;

    // Resolution ranges
    ResoRange resrange = hkl_list.ResRange();
    if (controls.analysis.NresoBins() > 0) {
      // set number of ranges
      resrange.SetNbins(controls.analysis.NresoBins());
    }
    reflection this_refl;
    observation this_obs;
    int index;
    
    int nrbins =  resrange.Nbins();
    // sum->mean I (rotation, resolution)
    clipper::Array2d<double> sumI(nrotranges, nrbins);
    // numbers
    clipper::Array2d<int>      nI(nrotranges, nrbins);

    for (int i=0;i<nrotranges;++i) {
      for (int j=0;j<nrbins;++j) {
	sumI(i,j) = 0.0; nI(i,j) = 0;}
    }

    hkl_list.rewind();
    while (hkl_list.next_reflection(this_refl) >= 0)  {  // loop reflections
      int ires = resrange.bin(this_refl.invresolsq());
      while ((index = this_refl.next_observation(this_obs)) >= 0) {
	// loop observations
	int irun = this_obs.run();
	Rtype phi = this_obs.phi();
	int irot = runlist[irun].PhiRange().bin(phi) + idxrun[irun];
	sumI(irot, ires) += this_obs.I(); 
	nI(irot, ires)++;
      }
    }

    // Compute average intensities
    std::vector<int> numobsrotrange(nrotranges, 0); // number of observations for rotrange
    for (int i=0;i<nrotranges;++i) {
      for (int j=0;j<nrbins;++j) {
	if (nI(i,j) > 0) {
	  sumI(i,j) /= double(nI(i,j));
	  numobsrotrange[i] += nI(i,j);
	}
      }
   }

   InitialData data(sumI);  // make data accessible to scale refinement

   FoxHolmes fh(data);
   int Ncyc =  5;  // number of cycles
   phaser::protocolPtr cPtr(new phaser::ProtocolScale(Ncyc));   // default protocols
   //^   std::cout << " InitScales::Minimizer::run::Number of cycles " << cPtr->getNCYC() << "\n";
   phaser::Minimizer Min;

   Min.run(fh, cPtr, output);

   // Print initial scales
   const int nperline = 10;
       
   std::vector<double> gscales(nrotranges);      // inverse scales (g)
   std::vector<double> scales(nrotranges, 0.0);  // scales
   int irr = 0;
   for (int irun=0;irun<nruns;++irun) {
     output.logTabPrintf(0, LOGFILE,
			 "\nInitial scales for run %5d\n", runlist[irun].RunNumber());
     int i1 = idxrun[irun];
     int i2 = nrotranges;
     if (irun < nruns-1) i2 = idxrun[irun+1];  // not last run
     for (int i=i1;i<i2;++i) {
       gscales[irr] = fh.getRefinePars()[irr];
       if (gscales[irr] != 0.0) scales[irr] = 1./gscales[irr];
       if ((i-i1) > 0 && (i-i1)%nperline == 0) output.logTabPrintf(0,LOGFILE,"\n");
       output.logTabPrintf(0,LOGFILE," %9.3f", scales[irr]);
       irr++;
     }
     output.logTabPrintf(0,LOGFILE,"\n");
   }
   // Store initial scales   
   AllScales.SetInitialScales(gscales, numobsrotrange);

  }  // InitialScales

} // namespace scala
