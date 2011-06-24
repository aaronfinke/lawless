//
//
//  AnalyseAnom
//

#include "analyseanom.hh"
#include "selectedobservations.hh"
using phaser_io::LOGFILE;

namespace scala
{
  std::vector<float> AnalyseAnom(const hkl_unmerge_list& hkl_list,
				 const SDmodel& SDM, all_controls& controls, const bool& plot,
				 phaser_io::Output& output)
  // Analyse anomalous differences & adjust anomalous rejection criterion (in controls)
  //
  // For now just use the crude method as in Scala,
  // ie get slope of DelAnom normal probability plot and add a multiplier of that 
  //
  // Returns mid-slopes of DelAnom normal probability plot for each dataset
  {
    int Ndatasets = hkl_list.num_datasets();
    std::vector<float> slopes(Ndatasets,1.0);
    if (!controls.Anomalous) return slopes;

    controls.outlierMerge.SetNdatasets(Ndatasets);

    NormalProbPlot NPPlot;
    if (plot) {
      NPPlot.init("ANOMPLOT", true, "Anomalous differences","");
    }

    std::vector<NormalProbAnal> normalprobanal;
    int nrej = AccumulateDelanomNormProb(SDM, hkl_list, controls.outlierMerge,
					 normalprobanal);
    nrej = nrej;
    // Magic factor to multiply slope by to inflate SD correction
    const float FACTOR = 3.7;

    output.logTab(0,LOGFILE,"\n\nNormal probability analysis of anomalous differences");
    output.logTab(0,LOGFILE,    "====================================================\n");
    output.logTabPrintf(1,LOGFILE,
			"     All data                     Data within expected delta %5.2f\n",
			normalprobanal[0].DltLim());
    if (Ndatasets > 1) {
      output.logTab(1,LOGFILE,
		    "  Slope Intercept   Number            Slope Intercept   Number    Dataset");
    } else {
      output.logTab(1,LOGFILE,
		    "Slope Intercept     Number            Slope Intercept   Number");
    }
    // loop datasets
    for (int id=0;id<Ndatasets;id++) {
      //      float slope = normalprobanal[id].Slope();  // central slope
      output.logTabPrintf(1,LOGFILE,"%7.2f %7.2f %10d          %7.2f %7.2f %10d ",
			  normalprobanal[id].Slope(0.0,true),
			  normalprobanal[id].Intercept(0.0,true),
			  normalprobanal[id].Number(0.0,true),
			  normalprobanal[id].Slope(),
			  normalprobanal[id].Intercept(),
			  normalprobanal[id].Number());
      if (Ndatasets > 1) {
	output.logTab(0,LOGFILE,
		      hkl_list.xdataset(id).pxdname().format());
      } else {
	output.logTab(0,LOGFILE,"  ");
      }
    }
    output.logTabPrintf(0,LOGFILE,"\nOutlier rejection limits for I+ v I-\n");
    for (int id=0;id<Ndatasets;id++) {
      slopes[id] = normalprobanal[id].Slope();
      // Dataset name
      if (Ndatasets > 1) {
	std::string pxdname = hkl_list.xdataset(id).pxdname().format();
	output.logTab(0,LOGFILE, "\n-- For dataset "+pxdname);
      }
      // Update reject limits
      controls.outlierMerge.Reject(BOTH, id).sdrej  *= FACTOR * normalprobanal[id].Slope();
      controls.outlierMerge.Reject(BOTH, id).sdrej2 *= FACTOR * normalprobanal[id].Slope();
      output.logTabPrintf(0,LOGFILE,
			  "    have been adjusted by a factor %7.2f * %7.2f\n\n",
			  FACTOR, normalprobanal[id].Slope());
      output.logTab(0,LOGFILE,
		    controls.outlierMerge.Reject(BOTH, id).format());

      std::string dname = hkl_list.xdataset(id).pxdname().dname();
      if (plot) normalprobanal[id].Plot(NPPlot, dname);
    }  // end loop datasets
    if (plot) NPPlot.ClosePlot();
    return slopes;
  }
  //--------------------------------------------------------------
  int AccumulateDelanomNormProb(const SDmodel& SDM,
				const hkl_unmerge_list& hkl_list,
				const OutlierControl& outliercontrol,
				std::vector<NormalProbAnal>& normalprobanal)
  // Go through reflection list
  //   apply current SD correction model
  //   reject outliers within I+ or I- sets
  //   accumulate deltaAnom into normalprobanal
  // 
  {  
    int Ndatasets = hkl_list.num_datasets();
    normalprobanal.resize(Ndatasets);

    reflection this_refl;
    int Nrej = 0;
    int nref = 0;
    hkl_list.rewind();
  
    while (hkl_list.next_reflection(this_refl) >= 0)  {
      if (!(hkl_list.symmetry().is_centric(this_refl.hkl()))) {
	// Acentric only
	//  Apply current SD correction to reflection (all observations)
	SDM.CorrectReflection(this_refl);
	nref++;
	// loop datasets
	for (int id=0;id<Ndatasets;id++) {
	  // Anomalous, treat I+ & I- separately
	  SelectedObservations Selobs(this_refl, id, IPLUS);
	  // Reject outliers within sets
	  Nrej += Selobs.Outliers(outliercontrol.Reject(IPLUS));
	  if (Selobs.Number() > 0) {
	    IsigI Iplus = Selobs.Average();
	    Selobs = SelectedObservations(this_refl, id, IMINUS);
	    // Reject outliers
	    Nrej += Selobs.Outliers(outliercontrol.Reject(IMINUS));
	    if (Selobs.Number() > 0) {
	      IsigI Iminus = Selobs.Average();
	      // DelAnom
	      float sig = Iplus.sigI()*Iplus.sigI() + Iminus.sigI()*Iminus.sigI();
	      if (sig > 0.0) {
		normalprobanal[id].AddDelta((Iplus.I() - Iminus.I())/sqrt(sig));
	      }
	    }
	  }
	}
      } // end loop datasets
    } // end loop reflections
    //^
    return Nrej;
  }
}
