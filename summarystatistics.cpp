//
// summarystatistics.cpp
//
// Overall summary data class (== "Table 1")
//

#include "summarystatistics.hh"
#include "string_util.hh"

using phaser_io::LOGFILE;
using phaser_io::RESULT;
using phaser_io::LXML;

namespace scala {
  // ------------------------------------------------------------
  void SummaryStatistics::StoreResRanges(const ResoRange& overall,
					 const ResoRange& inner,
					 const ResoRange& outer)  {
    resRange = Store3val<ResoRange>(overall, inner, outer);
  }
  // ------------------------------------------------------------
  // store Rmerge, overall, inner shell, outer shell
  void SummaryStatistics::StoreRmergeReso(const Rfactor& overall,
					  const Rfactor& inner,
					  const Rfactor& outer)  {
    rmerge = Store3val<Rfactor>(overall, inner, outer);
  }
  // ------------------------------------------------------------
  // store Rmeas, overall, inner shell, outer shell
  void SummaryStatistics::StoreRmeasReso(const Rfactor& overall,
					 const Rfactor& inner,
					 const Rfactor& outer)  {
    rmeas = Store3val<Rfactor>(overall, inner, outer);
  }
  // ------------------------------------------------------------
  // store Rpim, overall, inner shell, outer shell1
  void SummaryStatistics::StoreRpimReso(const Rfactor& overall,
					const Rfactor& inner,
					const Rfactor& outer)  {
    rpim = Store3val<Rfactor>(overall, inner, outer);
  }
  // ------------------------------------------------------------
  // store Rmerge, overall, inner shell, outer shell, v. overall mean I+-
  void SummaryStatistics::StoreRmergeResoOv(const Rfactor& overall,
					    const Rfactor& inner,
					    const Rfactor& outer)  {
    rmergeOv = Store3val<Rfactor>(overall, inner, outer);
  }
  // ------------------------------------------------------------
  // store Rmeas, overall, inner shell, outer shell, v. overall mean I+-
  void SummaryStatistics::StoreRmeasResoOv(const Rfactor& overall,
					   const Rfactor& inner,
					   const Rfactor& outer)  {
    rmeasOv = Store3val<Rfactor>(overall, inner, outer);
  }
  // ------------------------------------------------------------
  // store Rpim, overall, inner shell, outer shell, v. overall mean I+-
  void SummaryStatistics::StoreRpimResoOv(const Rfactor& overall,
					  const Rfactor& inner,
					  const Rfactor& outer)  {
    rpimOv = Store3val<Rfactor>(overall, inner, outer);
  }
  // ------------------------------------------------------------
  // store Rpim, overall, inner shell, outer shell, v. overall mean I+-
  // Numbers
  void SummaryStatistics::StoreNumbers(const int& overallNobs,
				       const int& innerNobs,
				       const int& outerNobs,
				       const int& overallNuniq,
				       const int& innerNuniq,
				       const int& outerNuniq) {
    Nobs = Store3val<int>(overallNobs, innerNobs, outerNobs);
    Nuniq = Store3val<int>(overallNuniq, innerNuniq, outerNuniq);
  }
  // ------------------------------------------------------------  
  // Mean(I/sd)
  void SummaryStatistics::StoreMnIsd(const float& overallMnIsd,
				     const float& innerMnIsd,
				     const float& outerMnIsd) {
    MnIsd = Store3val<float>(overallMnIsd, innerMnIsd, outerMnIsd);
  }
  // ------------------------------------------------------------  
  // <I> half-dataset CC
  void SummaryStatistics::StoreImeanCorrel(const float& overall,
					   const float& inner,
					   const float& outer) {
    Icorrelation = Store3val<float>(overall, inner, outer);
  }
  // ------------------------------------------------------------  
  // Completeness & multiplicity
  void SummaryStatistics::StoreCmplMult(const float& overallComplete,
					const float& innerComplete,
					const float& outerComplete,
					const float& overallMult,
					const float& innerMult,
					const float& outerMult) {
    complete = Store3val<float>(overallComplete, innerComplete, outerComplete);
    multiplicity = Store3val<float>(overallMult, innerMult, outerMult);
  }
  // ------------------------------------------------------------  
  // Anomalous completeness & multiplicity
  void SummaryStatistics::StoreAnomCmplMult(const float& overallComplete,
					    const float& innerComplete,
					    const float& outerComplete,
					    const float& overallMult,
					    const float& innerMult,
					    const float& outerMult) {
    anomcomplete = Store3val<float>(overallComplete, innerComplete, outerComplete);
    anommultiplicity = Store3val<float>(overallMult, innerMult, outerMult);
  }
  // ------------------------------------------------------------  
  // Anomalous correlation of half-datasets
  void SummaryStatistics::StoreAnomCorrel(const float& overall,
					  const float& inner,
					  const float& outer) {
    anomcorrelation = Store3val<float>(overall, inner, outer);
  }
  // ------------------------------------------------------------  
  // Anomalous RMS correlation ratio of half-datasets
  void SummaryStatistics::StoreAnomRCR(const float& overall,
				       const float& inner,
				       const float& outer) {
    anomRCR = Store3val<float>(overall, inner, outer);
  }
  // ------------------------------------------------------------  
  // slope of anomalous normal probability
  void SummaryStatistics::StoreAnomNPslope(const float& anomnpslope) {
    anomNPslope = anomnpslope;
  }
  // ------------------------------------------------------------  
  // Range of SD corrections
  void SummaryStatistics::StoreSDcorrectioRange(const float& minsdcorrfulls, const float& maxsdcorrfulls,
						const float& minsdcorrpartials, const float& maxsdcorrpartials)
  {
    minSDcorrFulls = minsdcorrfulls;
    maxSDcorrFulls = maxsdcorrfulls;
    minSDcorrPartials = minsdcorrpartials;
    maxSDcorrPartials = maxsdcorrpartials;
  }
  // ------------------------------------------------------------  
  // resolution limit estimates
  // overall limit from half-dataset CCs
  void SummaryStatistics::StoreHalfdatsetCCresolimit
  (const ResolutionLimit& OverallResolimitCC)
  {
    overallresolimitCC = OverallResolimitCC;
  }
  // ------------------------------------------------------------  
  // overall limit from Mn(I/sd)
  void SummaryStatistics::StoreMnIsigresolimit
  (const ResolutionLimit& OverallResolimitIsig)
  {
    overallresolimitIsig =OverallResolimitIsig;
  }
  // ------------------------------------------------------------  
  // anisotropic limits from half-dataset CCs
  void SummaryStatistics::StoreHalfdatsetCCAnisoresolimit
  (const std::vector<ResolutionLimit>& AnisoresolimitCC)
  {
    anisoresolimitCC = AnisoresolimitCC;
  }
  // ------------------------------------------------------------  
  // anisotropic limits from Mn(I/sd)
  void SummaryStatistics::StoreMnIsigAnisoresolimit
  (const std::vector<ResolutionLimit>& AnisoresolimitIsig)
  {
    anisoresolimitIsig = AnisoresolimitIsig;
  }
  // ------------------------------------------------------------  
  std::string ResoLimitWarning(const ResolutionLimit& reslimit)
  {
    if (reslimit.Status() > 0) {
      return " == maximum resolution";
    } else if (reslimit.Status() < 0) {
      return "WARNING: weak data, all data below threshold";
    }
    return "";
  }
  // ------------------------------------------------------------  
  void SummaryStatistics::PrintSummaryTable(const bool& Result, phaser_io::Output& output)
  // print the final summary table as RESULT if Result == true
  {
    phaser_io::outStream OUTSTREAM = LOGFILE;
    if (Result) {OUTSTREAM = RESULT;}

    output.logTabPrintf(0,OUTSTREAM,"Summary data for  ");
    output.logTab(0,OUTSTREAM,pxdname.formatPrint());
    output.logTab(0,OUTSTREAM,
     "\n                                           Overall  InnerShell  OuterShell");
    output.logTabPrintf(0,OUTSTREAM,
			"Low resolution limit                  %10.2f%10.2f%10.2f\n",
			resRange[0].ResLow(),resRange[1].ResLow(),resRange[2].ResLow());
    output.logTabPrintf(0,OUTSTREAM,
			"High resolution limit                 %10.2f%10.2f%10.2f\n\n",
			resRange[0].ResHigh(),resRange[1].ResHigh(),resRange[2].ResHigh());
    // Always write out R-factors within I+/I- sets and overall
    output.logTabPrintf(0,OUTSTREAM,
			"Rmerge  (within I+/I-)                %10.3f%10.3f%10.3f\n",
			rmerge[0].R(), rmerge[1].R(), rmerge[2].R());
    output.logTabPrintf(0,OUTSTREAM,
			"Rmerge  (all I+ and I-)               %10.3f%10.3f%10.3f\n",
			rmergeOv[0].R(), rmergeOv[1].R(), rmergeOv[2].R());
    output.logTabPrintf(0,OUTSTREAM,
			"Rmeas (within I+/I-)                  %10.3f%10.3f%10.3f\n",
			rmeas[0].R(), rmeas[1].R(), rmeas[2].R());
    output.logTabPrintf(0,OUTSTREAM,
			"Rmeas (all I+ & I-)                   %10.3f%10.3f%10.3f\n",
			  rmeasOv[0].R(), rmeasOv[1].R(), rmeasOv[2].R());
    output.logTabPrintf(0,OUTSTREAM,
			"Rpim (within I+/I-)                   %10.3f%10.3f%10.3f\n",
			rpim[0].R(), rpim[1].R(), rpim[2].R());
    output.logTabPrintf(0,OUTSTREAM,
			"Rpim (all I+ & I-)                    %10.3f%10.3f%10.3f\n",
			rpimOv[0].R(), rpimOv[1].R(), rpimOv[2].R());
    output.logTabPrintf(0,OUTSTREAM,
			"Rmerge in top intensity bin           %10.3f        -         - \n",
			RmergeTopI.R());
    output.logTabPrintf(0,OUTSTREAM,
			"Total number of observations          %10d%10d%10d\n",
			Nobs[0], Nobs[1], Nobs[2]);
    output.logTabPrintf(0,OUTSTREAM,
			"Total number unique                   %10d%10d%10d\n",
			Nuniq[0], Nuniq[1], Nuniq[2]);
    output.logTabPrintf(0,OUTSTREAM,
			"Mean((I)/sd(I))                       %10.1f%10.1f%10.1f\n",
			MnIsd[0], MnIsd[1], MnIsd[2]);
    output.logTabPrintf(0,OUTSTREAM,
			"<I> correlation between half-sets     %10.3f%10.3f%10.3f\n",
			Icorrelation[0], Icorrelation[1], Icorrelation[2]);

    output.logTabPrintf(0,OUTSTREAM,
			"Completeness                          %10.1f%10.1f%10.1f\n",
			complete[0], complete[1], complete[2]);
    output.logTabPrintf(0,OUTSTREAM,
			"Multiplicity                          %10.1f%10.1f%10.1f\n",
			multiplicity[0], multiplicity[1], multiplicity[2]);
    if (Anom) {
      output.logTabPrintf(0,OUTSTREAM,
			  "\nAnomalous completeness                %10.1f%10.1f%10.1f\n",
			  anomcomplete[0], anomcomplete[1], anomcomplete[2]);
      output.logTabPrintf(0,OUTSTREAM,
			  "Anomalous multiplicity                %10.1f%10.1f%10.1f\n",
			  anommultiplicity[0], anommultiplicity[1], anommultiplicity[2]);
      output.logTabPrintf(0,OUTSTREAM,
			  "DelAnom correlation between half-sets %10.3f%10.3f%10.3f\n",
			  anomcorrelation[0], anomcorrelation[1], anomcorrelation[2]);
      output.logTabPrintf(0,OUTSTREAM,
			  "Mid-Slope of Anom Normal Probability  %10.3f       -         -  \n",
			  anomNPslope);
    }

    // Resolution limit estimates
    output.logTab(0,OUTSTREAM,
		  "\nEstimates of resolution limits: overall");
    output.logTabPrintf(1,OUTSTREAM,
		"from half-dataset correlation coefficient > %5.2f: limit = %5.2fA %s\n",
			overallresolimitCC.Limit(),
			overallresolimitCC.HighResolution(),
			ResoLimitWarning(overallresolimitCC).c_str());

    output.logTabPrintf(1,OUTSTREAM,
		"from Mn(I/sd) > %5.2f:                             limit = %5.2fA %s\n",
			overallresolimitIsig.Limit(),
			overallresolimitIsig.HighResolution(),
			ResoLimitWarning(overallresolimitIsig).c_str());
    // Anisotropy analysis
    output.logTab(0,OUTSTREAM,
		  "\nEstimates of resolution limits along reciprocal lattice axes:");
    std::vector<std::string> axisname(3);
    axisname[0] = "a*";
    axisname[1] = "b*";
    axisname[2] = "c*";
    for (int jax=0;jax<3;++jax) {
      if (anisoresolimitCC[jax].Status() > -2) { // valid direction
	output.logTab(0,OUTSTREAM, "  Along axis "+axisname[jax]);
	output.logTabPrintf(1,OUTSTREAM,
	    "from half-dataset correlation coefficient > %5.2f: limit = %5.2fA %s\n",
			    anisoresolimitCC[jax].Limit(),
			    anisoresolimitCC[jax].HighResolution(),
			    ResoLimitWarning(anisoresolimitCC[jax]).c_str());
	output.logTabPrintf(1,OUTSTREAM,
	    "from Mn(I/sd) > %5.2f:                             limit = %5.2fA %s\n",
			    anisoresolimitIsig[jax].Limit(),
			    anisoresolimitIsig[jax].HighResolution(),
			    ResoLimitWarning(anisoresolimitIsig[jax]).c_str());
      }
    }

    output.logTab(0,OUTSTREAM,
		  "\nAverage unit cell: "+averageCell.format());
    output.logTab(0,OUTSTREAM,
		  "\nSpace group: "+spacegroupname);
    output.logTabPrintf(0,OUTSTREAM,
			"\nAverage mosaicity: %6.2f\n", averageMosaicity);
    output.logTabPrintf(0,OUTSTREAM,
     "\nMinimum and maximum SD correction factors: Fulls %6.2f %6.2f Partials %6.2f %6.2f\n",  
			minSDcorrFulls, maxSDcorrFulls, minSDcorrPartials, maxSDcorrPartials);
  }  // PrintSummaryTable 
  // ------------------------------------------------------------  
  void AllSummaryStatistics::AddSummaryStatistics(const SummaryStatistics& summarystatistics)
  //! store statistics for one dataset
  {
    allsummarystatistics.push_back(summarystatistics);
  }
  // ------------------------------------------------------------  
  //! print the final summary table for one dataset (idts), as RESULT if Result true
  void AllSummaryStatistics::PrintOneSummaryTable(const int& idts,
				     const bool& Result, phaser_io::Output& output)
  {
    if (idts >= 0 && idts < int(allsummarystatistics.size())) {
      allsummarystatistics[idts].PrintSummaryTable(Result, output);
    }
  }
  // ------------------------------------------------------------  
  //! print the final summary table for all datasets as RESULT if Result true
  void AllSummaryStatistics::PrintSummaryTable(const bool& Result,
					       const bool& Anom, phaser_io::Output& output)
  {
    // // allsummarystatistics[idts] 

    phaser_io::outStream OUTSTREAM = LOGFILE;
    if (Result) {OUTSTREAM = RESULT;}
    int ndts = allsummarystatistics.size(); // number of datasets
    if (ndts == 0) {return;}
    if (ndts == 1) {
      allsummarystatistics[0].PrintSummaryTable(Result, output);
      return;
    }

    // Column widths for Overall/Inner/Outer = 9/8/8 for each derivative
    // + 2 between them

    output.logTab(0,OUTSTREAM,"\nSummary data for datasets");
    for (int idts=0;idts<ndts;++idts) {
      output.logTab(0,OUTSTREAM,"   "+allsummarystatistics[idts].pxdname.formatPrint());
    }
    output.logTab(0,OUTSTREAM, "\n                                    ", false);
    for (int idts=0;idts<ndts;++idts) {
      output.logTab(0,OUTSTREAM,
	    StringUtil::CentreString(allsummarystatistics[idts].pxdname.format(),27),
		    false);
    }
    output.logTab(0,OUTSTREAM," ");

    output.logTab(0,OUTSTREAM, "\n                                    ", false);
    for (int idts=0;idts<ndts;++idts) {
      output.logTab(0,OUTSTREAM, "    Overall   Inner   Outer", false);
    }
    output.logTab(0,OUTSTREAM," ");

    output.logTab(0,OUTSTREAM, "Low resolution limit                ", false);
    for (int idts=0;idts<ndts;++idts) {
      output.logTabPrintf(0,OUTSTREAM, "  %9.2f%8.2f%8.2f",
			  allsummarystatistics[idts].resRange[0].ResLow(),
			  allsummarystatistics[idts].resRange[1].ResLow(),
			  allsummarystatistics[idts].resRange[2].ResLow());
    }
    output.logTab(0,OUTSTREAM," ");

    output.logTab(0,OUTSTREAM, "High resolution limit               ", false);
    for (int idts=0;idts<ndts;++idts) {
      output.logTabPrintf(0,OUTSTREAM, "  %9.2f%8.2f%8.2f",
			  allsummarystatistics[idts].resRange[0].ResHigh(),
			  allsummarystatistics[idts].resRange[1].ResHigh(),
			  allsummarystatistics[idts].resRange[2].ResHigh());
    }
    output.logTab(0,OUTSTREAM," ");
    // Always write out R-factors within I+/I- sets and overall
    output.logTab(0,OUTSTREAM, "\nRmerge  (within I+/I-)              ", false);
    for (int idts=0;idts<ndts;++idts) {
      output.logTabPrintf(0,OUTSTREAM, "  %9.3f%8.3f%8.3f",
			allsummarystatistics[idts].rmerge[0].R(),
			allsummarystatistics[idts].rmerge[1].R(),
			allsummarystatistics[idts].rmerge[2].R());
    }
    output.logTab(0,OUTSTREAM," ");
    output.logTab(0,OUTSTREAM, "Rmerge  (all I+ and I-)             ", false);
    for (int idts=0;idts<ndts;++idts) {
      output.logTabPrintf(0,OUTSTREAM, "  %9.3f%8.3f%8.3f",
			  allsummarystatistics[idts].rmergeOv[0].R(),
			  allsummarystatistics[idts].rmergeOv[1].R(),
			  allsummarystatistics[idts].rmergeOv[2].R());
    }
    output.logTab(0,OUTSTREAM," ");
    output.logTab(0,OUTSTREAM, "Rmeas (within I+/I-)                ", false);
    for (int idts=0;idts<ndts;++idts) {
      output.logTabPrintf(0,OUTSTREAM, "  %9.3f%8.3f%8.3f",
			  allsummarystatistics[idts].rmeas[0].R(),
			  allsummarystatistics[idts].rmeas[1].R(),
			  allsummarystatistics[idts].rmeas[2].R());
    }
    output.logTab(0,OUTSTREAM," ");
    output.logTab(0,OUTSTREAM, "Rmeas (all I+ & I-)                 ", false);
    for (int idts=0;idts<ndts;++idts) {
      output.logTabPrintf(0,OUTSTREAM, "  %9.3f%8.3f%8.3f",
			  allsummarystatistics[idts].rmeasOv[0].R(),
			  allsummarystatistics[idts].rmeasOv[1].R(),
			  allsummarystatistics[idts].rmeasOv[2].R());
    }
    output.logTab(0,OUTSTREAM," ");
    output.logTab(0,OUTSTREAM, "Rpim (within I+/I-)                 ", false);
    for (int idts=0;idts<ndts;++idts) {
      output.logTabPrintf(0,OUTSTREAM, "  %9.3f%8.3f%8.3f",
			  allsummarystatistics[idts].rpim[0].R(),
			  allsummarystatistics[idts].rpim[1].R(),
			  allsummarystatistics[idts].rpim[2].R());
    }
    output.logTab(0,OUTSTREAM," ");
    output.logTab(0,OUTSTREAM, "Rpim (all I+ & I-)                  ", false);
    for (int idts=0;idts<ndts;++idts) {
      output.logTabPrintf(0,OUTSTREAM, "  %9.3f%8.3f%8.3f",
			  allsummarystatistics[idts].rpimOv[0].R(),
			  allsummarystatistics[idts].rpimOv[1].R(),
			  allsummarystatistics[idts].rpimOv[2].R());
    }
    output.logTab(0,OUTSTREAM," ");
    output.logTab(0,OUTSTREAM, "Rmerge in top intensity bin         ", false);
    for (int idts=0;idts<ndts;++idts) {
      output.logTabPrintf(0,OUTSTREAM, "  %9.3f      -       - ",
			allsummarystatistics[idts].RmergeTopI.R());
    }
    output.logTab(0,OUTSTREAM," ");
    output.logTab(0,OUTSTREAM, "Total number of observations        ", false);
    for (int idts=0;idts<ndts;++idts) {
      output.logTabPrintf(0,OUTSTREAM, "  %9d%8d%8d",
			allsummarystatistics[idts].Nobs[0],
			allsummarystatistics[idts].Nobs[1],
			allsummarystatistics[idts].Nobs[2]);
    }
    output.logTab(0,OUTSTREAM," ");
    output.logTab(0,OUTSTREAM, "Total number unique                 ", false);
    for (int idts=0;idts<ndts;++idts) {
      output.logTabPrintf(0,OUTSTREAM, "  %9d%8d%8d",
			allsummarystatistics[idts].Nuniq[0],
			allsummarystatistics[idts].Nuniq[1],
			allsummarystatistics[idts].Nuniq[2]);
    }
    output.logTab(0,OUTSTREAM," ");
    output.logTab(0,OUTSTREAM, "Mean((I)/sd(I))                     ", false);
    for (int idts=0;idts<ndts;++idts) {
      output.logTabPrintf(0,OUTSTREAM, "  %9.1f%8.1f%8.1f",
			allsummarystatistics[idts].MnIsd[0],
			allsummarystatistics[idts].MnIsd[1],
			allsummarystatistics[idts].MnIsd[2]);
    }
    output.logTab(0,OUTSTREAM," ");
    output.logTab(0,OUTSTREAM, "Mn(I) correlation between half-sets ", false);
    for (int idts=0;idts<ndts;++idts) {
      output.logTabPrintf(0,OUTSTREAM,"  %9.3f%8.3f%8.3f",
			  allsummarystatistics[idts].Icorrelation[0],
			  allsummarystatistics[idts].Icorrelation[1],
			  allsummarystatistics[idts].Icorrelation[2]);
    }
    output.logTab(0,OUTSTREAM," ");
    output.logTab(0,OUTSTREAM, "Completeness                        ", false);
    for (int idts=0;idts<ndts;++idts) {
      output.logTabPrintf(0,OUTSTREAM, "  %9.1f%8.1f%8.1f",
			  allsummarystatistics[idts].complete[0],
			  allsummarystatistics[idts].complete[1],
			  allsummarystatistics[idts].complete[2]);
    }
    output.logTab(0,OUTSTREAM," ");
    output.logTab(0,OUTSTREAM, "Multiplicity                        ", false);
    for (int idts=0;idts<ndts;++idts) {
      output.logTabPrintf(0,OUTSTREAM, "  %9.1f%8.1f%8.1f",
			allsummarystatistics[idts].multiplicity[0],
			allsummarystatistics[idts].multiplicity[1],
			  allsummarystatistics[idts].multiplicity[2]);
    }
    output.logTab(0,OUTSTREAM," ");
    if (Anom) {
      output.logTab(0,OUTSTREAM,
		    "\nAnomalous completeness              ", false);
      for (int idts=0;idts<ndts;++idts) {
	output.logTabPrintf(0,OUTSTREAM, "  %9.1f%8.1f%8.1f",
			    allsummarystatistics[idts].anomcomplete[0],
			    allsummarystatistics[idts].anomcomplete[1],
			    allsummarystatistics[idts].anomcomplete[2]);
      }
      output.logTab(0,OUTSTREAM," ");
      output.logTabPrintf(0,OUTSTREAM,
			  "Anomalous multiplicity              ", false);
      for (int idts=0;idts<ndts;++idts) {
	output.logTabPrintf(0,OUTSTREAM, "  %9.1f%8.1f%8.1f",
			  allsummarystatistics[idts].anommultiplicity[0],
			  allsummarystatistics[idts].anommultiplicity[1],
			    allsummarystatistics[idts].anommultiplicity[2]);
      }
      output.logTab(0,OUTSTREAM," ");
      output.logTab(0,OUTSTREAM, "DelAnom correlation between half-sets ",false);
      for (int idts=0;idts<ndts;++idts) {
	output.logTabPrintf(0,OUTSTREAM,"  %9.3f%8.3f%8.3f",
			    allsummarystatistics[idts].anomcorrelation[0],
			    allsummarystatistics[idts].anomcorrelation[1],
			    allsummarystatistics[idts].anomcorrelation[2]);
      }
      output.logTab(0,OUTSTREAM," ");
      output.logTab(0,OUTSTREAM, "Mid-Slope of Anom Normal Probability  ",false);
      for (int idts=0;idts<ndts;++idts) {
	output.logTabPrintf(0,OUTSTREAM,"  %9.3f     -       -  ",
			    allsummarystatistics[idts].anomNPslope);
      }
      output.logTab(0,OUTSTREAM," ");
    }
    // Resolution limit estimates
    output.logTab(0,OUTSTREAM,
		  "\nEstimates of resolution limits: overall");
    for (int idts=0;idts<ndts;++idts) {
      output.logTab(1,OUTSTREAM,
		    "Dataset: "+allsummarystatistics[idts].pxdname.format());
      output.logTabPrintf(2,OUTSTREAM,
		  "from half-dataset correlation coefficient > %5.2f: limit = %5.2fA %s\n",
			  allsummarystatistics[idts].overallresolimitCC.Limit(),
			  allsummarystatistics[idts].overallresolimitCC.HighResolution(),
			  ResoLimitWarning(allsummarystatistics[idts].overallresolimitCC).c_str());

      output.logTabPrintf(2,OUTSTREAM,
		"from Mn(I/sd) > %5.2f:                             limit = %5.2fA %s\n",
			allsummarystatistics[idts].overallresolimitIsig.Limit(),
			allsummarystatistics[idts].overallresolimitIsig.HighResolution(),
			ResoLimitWarning(allsummarystatistics[idts].overallresolimitIsig).c_str());
    }
    output.logTab(0,OUTSTREAM," ");

    // Anisotropy analysis
    output.logTab(0,OUTSTREAM,
		  "\nEstimates of resolution limits along reciprocal lattice axes:");

    std::vector<std::string> axisname(3);
    axisname[0] = "a*";
    axisname[1] = "b*";
    axisname[2] = "c*";
    for (int idts=0;idts<ndts;++idts) {
      output.logTab(1,OUTSTREAM,
		    "Dataset: "+allsummarystatistics[idts].pxdname.format());
      for (int jax=0;jax<3;++jax) {
	if (allsummarystatistics[idts].anisoresolimitCC[jax].Status() > -2) { // valid direction
	  output.logTab(2,OUTSTREAM, "  Along axis "+axisname[jax]);
	  output.logTabPrintf(2,OUTSTREAM,
			      "from half-dataset correlation coefficient > %5.2f: limit = %5.2fA %s\n",
			      allsummarystatistics[idts].anisoresolimitCC[jax].Limit(),
			      allsummarystatistics[idts].anisoresolimitCC[jax].HighResolution(),
			      ResoLimitWarning(allsummarystatistics[idts].anisoresolimitCC[jax]).c_str());
	  output.logTabPrintf(2,OUTSTREAM,
			      "from Mn(I/sd) > %5.2f:                             limit = %5.2fA %s\n",
			      allsummarystatistics[idts].anisoresolimitIsig[jax].Limit(),
			      allsummarystatistics[idts].anisoresolimitIsig[jax].HighResolution(),
			      ResoLimitWarning(allsummarystatistics[idts].anisoresolimitIsig[jax]).c_str());
	}
      }
    }
    output.logTab(0,OUTSTREAM," ");

    for (int idts=0;idts<ndts;++idts) {
      output.logTab(0,OUTSTREAM,
		    "\nDataset: "+allsummarystatistics[idts].pxdname.format()+"\n");
      output.logTab(1,OUTSTREAM,
		    "Average unit cell: "+allsummarystatistics[idts].averageCell.format());
      output.logTab(1,OUTSTREAM,
		  "Space group: "+allsummarystatistics[idts].spacegroupname);
      output.logTabPrintf(1,OUTSTREAM,
			  "Average mosaicity: %6.2f\n",
			  allsummarystatistics[idts].averageMosaicity);
      output.logTabPrintf(1,OUTSTREAM,
     "Minimum and maximum SD correction factors: Fulls %6.2f %6.2f Partials %6.2f %6.2f\n",  
			  allsummarystatistics[idts].minSDcorrFulls,
			  allsummarystatistics[idts].maxSDcorrFulls,
			  allsummarystatistics[idts].minSDcorrPartials,
			  allsummarystatistics[idts].maxSDcorrPartials);
    }
  }
  // ------------------------------------------------------------  
} // namespace scala
