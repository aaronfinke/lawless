//
//
//  AnalyseAnom
//

#include "analyseanom.hh"
#include "selectedobservations.hh"
using phaser_io::LOGFILE;
using phaser_io::LXML;

namespace scala
{
  AnalyseAnom::AnalyseAnom(const hkl_unmerge_list& hkl_list,
                           const SDmodel& SDM, all_controls& controls,
                           const ResoRange& ResRange,
                           const bool& plot,
                           phaser_io::Output& output)
  // Analyse anomalous differences & adjust anomalous rejection criterion (in controls)
  //
  // For now just use the crude method as in Scala,
  // ie get slope of DelAnom normal probability plot and add a multiplier of that
  //
  // if plot == true, write ANOMPLOT file
  {
    ndatasets = hkl_list.num_datasets();
    nresbin = ResRange.Nbins();
    slopes.assign(ndatasets,1.0);
    rmsdelanom.resize(ndatasets);
    for (int id=0;id<ndatasets;++id) {
      rmsdelanom[id].assign(nresbin,MeanSD());
    }

    controls.outlierMerge.SetNdatasets(ndatasets);

    NormalProbPlot NPPlot;
    if (plot) {
      NPPlot.init("ANOMPLOT", true, "Anomalous differences","");
    }

    std::vector<NormalProbAnal> normalprobanal;
    int nrej = AccumulateDelanomNormProb(SDM, hkl_list, controls.outlierMerge,
                                         ResRange,
                                         normalprobanal);
    nrej = nrej;
    // Magic factor to multiply slope by to inflate SD correction
    const float FACTOR = 3.7;

    output.logTab(0,LOGFILE,"\n\nNormal probability analysis of anomalous differences");
    output.logTab(0,LOGFILE,    "====================================================\n");
    output.logTabPrintf(1,LOGFILE,
                        "     All data                     Data within expected delta %5.2f\n",
                        normalprobanal[0].DltLim());
    if (ndatasets > 1) {
      output.logTab(1,LOGFILE,
                    "  Slope Intercept   Number            Slope Intercept   Number    Dataset");
    } else {
      output.logTab(1,LOGFILE,
                    "Slope Intercept     Number            Slope Intercept   Number");
    }
    // loop datasets
    for (int id=0;id<ndatasets;id++) {
      //      float slope = normalprobanal[id].Slope();  // central slope
      output.logTabPrintf(1,LOGFILE,"%7.2f %7.2f %10d          %7.2f %7.2f %10d ",
                          normalprobanal[id].Slope(0.0,true),
                          normalprobanal[id].Intercept(0.0,true),
                          normalprobanal[id].Number(0.0,true),
                          normalprobanal[id].Slope(),
                          normalprobanal[id].Intercept(),
                          normalprobanal[id].Number());
      if (ndatasets > 1) {
        output.logTab(0,LOGFILE,
                      hkl_list.dataset(id).formatNames());
      } else {
        output.logTab(0,LOGFILE,"  ");
      }
    }
    output.logTabPrintf(0,LOGFILE,"\nOutlier rejection limits for I+ v I-\n");
    for (int id=0;id<ndatasets;id++) {
      slopes[id] = normalprobanal[id].Slope();
      // Dataset name
      if (ndatasets > 1) {
        std::string pxdlabel = hkl_list.dataset(id).formatNames();
        output.logTab(0,LOGFILE, "\n-- For dataset "+pxdlabel);
      }
      // Update reject limits

      RejectFlags rejflags = controls.outlierMerge.Reject(BOTH, id);
      rejflags.sdrej  *= FACTOR * normalprobanal[id].Slope();
      rejflags.sdrej2 *= FACTOR * normalprobanal[id].Slope();
      controls.outlierMerge.SetReject(rejflags, BOTH, id);
      output.logTabPrintf(0,LOGFILE,
                          "    have been adjusted by a factor %7.2f * %7.2f\n\n",
                          FACTOR, normalprobanal[id].Slope());
      output.logTab(0,LOGFILE,
                    controls.outlierMerge.Reject(BOTH, id).format());

      std::string dname = hkl_list.dataset(id).Dname();
      if (plot) normalprobanal[id].Plot(NPPlot, dname);
    }  // end loop datasets
    if (plot) {
      NPPlot.ClosePlot();
      output.logTab(0,LXML,NPPlot.formatXML());
    }
  }
  //--------------------------------------------------------------
  int AnalyseAnom::AccumulateDelanomNormProb(const SDmodel& SDM,
                                             const hkl_unmerge_list& hkl_list,
                                             const OutlierControl& outliercontrol,
                                             const ResoRange& ResRange,
                                             std::vector<NormalProbAnal>& normalprobanal)
  // Go through reflection list
  //   apply current SD correction model
  //   reject outliers within I+ or I- sets
  //   accumulate deltaAnom into normalprobanal
  //
  {
    int ndatasets = hkl_list.num_datasets();
    normalprobanal.resize(ndatasets);

    reflection this_refl;
    int Nrej = 0;
    int nref = 0;
    hkl_list.rewind();

    while (hkl_list.next_reflection(this_refl) >= 0)  {
      if (!(hkl_list.symmetry().is_centric(this_refl.hkl()))) {
        Rtype invresolsq = this_refl.invresolsq();
        // Resolution bin
        int mres = ResRange.bin(invresolsq);
        // Acentric only
        //  Apply current SD correction to reflection (all observations)
        SDM.CorrectReflection(this_refl);
        nref++;
        // loop datasets
        for (int id=0;id<ndatasets;id++) {
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
                float danom = Iplus.I() - Iminus.I();
                normalprobanal[id].AddDelta(danom/sqrt(sig));
                rmsdelanom[id][mres].Add(danom);
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
