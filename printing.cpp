// printing.cpp

#include "printing.hh"
#include "version.hh"
#include "tablegraph.hh"
#include "numbercomplete.hh"
#include "halfdataset.hh"
#include "string_util.hh"

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
  // Make version line
  int sp = 52; // number of spaces in line
  int n1 = (sp - PROGRAM_VERSION.size())/2;
  int n2 = sp - PROGRAM_VERSION.size() - n1;
  std::string pad1(n1,' ');
  std::string pad2(n2,' ');
  std::string line = "        *"+pad1+PROGRAM_VERSION+pad2+"*\n";

  output.logTabPrintf(0,LOGFILE,
          "\n        ******************************************************\n");
  output.logTabPrintf(0,LOGFILE,
          "        *                                                    *\n");
  output.logTabPrintf(0,LOGFILE,
          "        *                      AIMLESS                       *\n");
  output.logTab(0,LOGFILE, line);
  output.logTabPrintf(0,LOGFILE,
          "        *                                                    *\n");
  output.logTabPrintf(0,LOGFILE,
          "        *     Scaling & analysis of unmerged intensities     *\n");
  output.logTabPrintf(0,LOGFILE,
          "        *     Phil Evans MRC LMB, Cambridge                  *\n");
  output.logTabPrintf(0,LOGFILE,
          "        *                                                    *\n");
  output.logTabPrintf(0,LOGFILE,
          "        ******************************************************\n\n");
}
//--------------------------------------------------------------
void PrintFileInfoToXML(const std::string& StreamName,
			const std::string& FileName,
			const Scell& cell,
			const std::string& SpaceGroupName,
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
		     "<SpacegroupName> "+SpaceGroupName+"</SpacegroupName>");
      output.logTab(0, LXML,
		     "</ReflectionFile>");
    }
}
//--------------------------------------------------------------
void PrintOutlierSettings(const all_controls& controls, phaser_io::Output& output)
{
  output.logTab(0,LOGFILE,"\nOutlier rejection parameters:");
  output.logTab(0,LOGFILE,"In scaling:");
  output.logTab(0,LOGFILE,controls.outlierScale.Reject(ALL).format());
  output.logTab(0,LOGFILE,controls.outlierScale.EMaxTest().format());
  output.logTab(0,LOGFILE,"\nIn merging:");
  output.logTab(0,LOGFILE,controls.outlierMerge.Reject(ALL).format());
  output.logTab(0,LOGFILE,controls.outlierMerge.EMaxTest().format());
  output.logTab(0,LOGFILE,"\n");
}
//--------------------------------------------------------------
// Print scale factors
void PrintScales(const ScaleModel& AllScales, phaser_io::Output& output)
{
}
//--------------------------------------------------------------
void PrintScalesByBatch(const PxdName& dataset_pxd,
			const std::vector<Batch>& batches, const std::vector<Run>& RunList,
			const int& datasetIndex,
			const std::vector<float>& scale0batch,
			const std::vector<float>& bfacbatch,
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
 
  
  // $TABLE  start 
  TableGraph table(" >>> Scales v rotation range, "+dataset_pxd.dname());
  output.logTab(0,LOGFILE,table.formatTitle());
  int c[] = {1,6,7};
  std::vector<int> cln(c,c+3);
  output.logTab(0,LOGFILE,
		table.Graph("Mn(k) & 0k (theta=0) v. batch","N",cln));
  int c2[] = {1,5};
  cln.assign(c2, c2+2);
  output.logTab(0,LOGFILE,
		table.Graph("Relative Bfactor v. batch","A",cln));
  std::vector<std::string> collabels;
  collabels.push_back("N");         // 1
  collabels.push_back("Run");       // 2
  collabels.push_back("Phi");       // 3
  collabels.push_back("Batch");     // 4
  collabels.push_back("Bfactor");   // 5
  collabels.push_back("Mn(k)");     // 6
  collabels.push_back("0k");        // 7
  collabels.push_back("Number");    // 8
  std::vector<bool> Zero(8, false);
  output.logTab(0,LOGFILE,
  		table.ColumnFields(collabels, Zero,
   				   "%5d%5d%8.2f%8d%10.2f%10.4f%10.4f%10d\n")); 

  int nc = collabels.size();

  int n=1;
  for (size_t i=0;i<batches.size();++i) {
    if (scalebatch[i].Count() > 0) {
      //      std::cout << n<< " " << batches[i].MidPhi()<< " " << batches[i].num()<< " " <<
      //	bfacbatch[i]<< " " << scalebatch[i].Mean()<< " " << scale0batch[i]<< " " << scalebatch[i].Count() << "\n";;

      output.logTab(0,LOGFILE,
		    table.Line(nc, n,
			       RunList[batches[i].RunIndex()].RunNumber(),
			       batches[i].MidPhi(), batches[i].num(),
			       bfacbatch[i], scalebatch[i].Mean(),
			       scale0batch[i], scalebatch[i].Count()));
      n++;
    }
  }
  output.logTab(0,LOGFILE,
		table.CloseTable());
  output.logTab(0,LOGFILE,table.RawLabels());
}
//--------------------------------------------------------------
void PrintDeviationsByBatch(const PxdName& dataset_pxd,
			    const std::vector<Batch>& batches,
			    const int& datasetIndex,
			    const std::vector<MeanSD>& imeanbatch,
			    const std::vector<MeanSD>& rmsDbatch,
			    const std::vector<Rfactor>& rmergebatch,
			    const std::vector<Rfactor>& rmergebatchsmoothed,
			    const std::vector<int>& rejectedbatch,
			    const std::vector<float>& batchcompleteness,
			    const std::vector<float>& batchanomcompleteness,
			    const std::vector<double>& maxresbatch,
			    const std::vector<double>& maxresbatchsmoothed,
			    const double& MinimumIoverSigma,
			    const int& nbatchsmooth,
			    const ResoRange& ResRange,
			    phaser_io::Output& output)
{
  output.logTab(0,LOGFILE,
		std::string("\n\nAgreement between batches\n")+
  			        "=========================\n\n");
  output.logTab(0,LOGFILE,
		std::string(" Rmerge in this table is the difference from Mn(Imean),\n")+
		"  but in later tables Rmerge is the difference from Mn(I+),Mn(I-)\n");
  
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
			  "\n SmRmerge and SmMaxRes in table are smoothed over %3d batches\n",
			  nbatchsmooth);
    } else {
      output.logTabPrintf(0,LOGFILE,
			  "\n SmRmerge in table is smoothed over %3d batches\n",
			  nbatchsmooth);
    }
  }

  // Get range of resolutions for maximum resolution plot
  ASSERT (batches.size() == maxresbatch.size());
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
  //	    <<" " << maxresbatch.size()<<"\n"; //^
  int nb = 0;
  for (size_t i=0;i<batches.size();++i) {
    if (imeanbatch[i].Count() > 0) {nb++;} // count actual batches
  }
  scala::Range xrange(1,nb);
  scala::Range yrange(res1, res2);

  TableGraph table
    (" Analysis against all Batches for all runs, "+dataset_pxd.dname());
  output.logTab(0,LOGFILE,table.formatTitle());
  int ncl = (smoothR) ? 3 : 2;
  std::vector<int> cln(ncl);
  cln[0] = 1;
  cln[1] = 6; // just unsmoothed
  if (smoothR) {
    cln[1] = 12; // smoothed one first so it comes out on top
    cln[2] = 6;
  }
  output.logTab(0,LOGFILE,
		table.Graph("Rmerge v Batch for all runs","N",cln));

  int c1[] = {1,9,10};
  cln.assign(c1,c1+3);
  output.logTab(0,LOGFILE,
		table.Graph("Cumulative %completeness & Anom%cmpl v Batch","N",cln));

  ncl = (smoothMaxRes) ? 3 : 2;
  cln.resize(ncl);
  cln[0] = 1;
  cln[1] = 11;
  if (smoothMaxRes) {
    cln[1] = 13;
    cln[2] = 11;
  }
  std::string s = "Maximum resolution limit, I/sigma > "+StringUtil::ftos(MinimumIoverSigma,5,1);
  output.logTab(0,LOGFILE,
		table.Graph(s,GraphAxesType(xrange,yrange,false),cln));
  int c3[] = {1,3,4};
  cln.assign(c3, c3+3);
  output.logTab(0,LOGFILE,
		table.Graph("Imean & RMS Scatter","N",cln));
  int c4[] = {1,5};
  cln.assign(c4, c4+2);
  output.logTab(0,LOGFILE,
		table.Graph("Imean/RMS scatter","N",cln));
  int c5[] = {1,8};
  cln.assign(c5, c5+2);
  output.logTab(0,LOGFILE,
		table.Graph("Number of rejects","N",cln));
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
  if (smoothR) { // if we have smoothed stats as well
    collabels.push_back("SmRmerge");     // 12
  }
  if (smoothMaxRes) {
    collabels.push_back("SmMaxRes");     // 13
  }  
  int nc = collabels.size();

  bool z[] = {false, false, true, true, true, true, false, false, true, true, false, true, true};
  std::vector<bool> Zero(z, z+nc);
  std::string lineformat = "%5d %7d %8.1f %8.1f %6.2f %7.3f %9d %5d %7.1f %6.1f %6.2f";
  if  (smoothR) {lineformat += " %8.3f";}
  if  (smoothMaxRes) {lineformat += " %8.2f";}
  lineformat += "\n";
  output.logTab(0,LOGFILE,
		table.ColumnFields(collabels, Zero, lineformat));

  int n=1;
  for (size_t i=0;i<batches.size();++i) {
    if (imeanbatch[i].Count() > 0) {
      float r = 0.0;
      if (rmsDbatch[i].Count() > 0) {
	r = imeanbatch[i].Mean()/sqrt(rmsDbatch[i].Mean());
      }

      if (nbatchsmooth == 1) { // no smoothed stats
	output.logTab(0,LOGFILE,
		      table.Line(nc, n, batches[i].num(),
				 imeanbatch[i].Mean(), sqrt(rmsDbatch[i].Mean()),
				 r,
				 rmergebatch[i].R(),
				 rmergebatch[i].result().count,
				 rejectedbatch[i],
				 100.*batchcompleteness[i],
				 100.*batchanomcompleteness[i],
				 maxresbatch[i]));
      } else if (smoothMaxRes) {
	output.logTab(0,LOGFILE,
		      table.Line(nc, n, batches[i].num(),
				 imeanbatch[i].Mean(), sqrt(rmsDbatch[i].Mean()),
				 r,
				 rmergebatch[i].R(), 
				 rmergebatch[i].result().count,
				 rejectedbatch[i],
				 100.*batchcompleteness[i],
				 100.*batchanomcompleteness[i],
				 maxresbatch[i],
				 rmergebatchsmoothed[i].R(),
				 maxresbatchsmoothed[i]));
      } else{
	output.logTab(0,LOGFILE,
		      table.Line(nc, n, batches[i].num(),
				 imeanbatch[i].Mean(), sqrt(rmsDbatch[i].Mean()),
				 r,
				 rmergebatch[i].R(), 
				 rmergebatch[i].result().count,
				 rejectedbatch[i],
				 100.*batchcompleteness[i],
				 100.*batchanomcompleteness[i],
				 maxresbatch[i],
				 rmergebatchsmoothed[i].R()));
      }
      n++;
    }
  }
  output.logTab(0,LOGFILE,
		table.CloseTable());
  output.logTab(0,LOGFILE,table.RawLabels());
} 
//--------------------------------------------------------------
void PrintDeviationsByResolution(const PxdName& dataset_pxd,
				 const ResoRange& ResRange,
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
				 const double& MinimumIoverSigma,
				 SummaryStatistics& summarystatistics,
				 phaser_io::Output& output)
{
  output.logTab(0,LOGFILE,
		std::string("\n Rmrg    :- conventional Rmerge = Sum(|Ihl - <Ih>|)/Sum(<Ih>)\n")+
		" Rcum    :- Rmrg up to this range\n"+
		" Rfull   :- Rmrg for fully-recorded observations only\n"+
		" Rmeas   :- multiplicity-independent R = Sum(Sqrt(N/(N-1))(|Ihl - <Ih>|))/Sum(<Ih>)\n"+
		" Rpim    :- Precision-indicating R = Sum(Sqrt(1/(N-1))(|Ihl - <Ih>|))/Sum(<Ih>)\n"+
		" Nmeas   :- Number of observations used in statistics\n"+
		" Av_I    :- unmerged Ihl averaged in bin <Ihl>\n"+
		" RMSdev  :- rms scatter of observations from mean <Ih>\n"+
		" I/RMS   :- <Ihl> / rms scatter  = Av_I/RMSdev\n"+
		" sd      :- average standard deviation derived from experimental SDs, after\n"+
		"             application of SdFac SdB SdAdd 'correction' terms\n"+
		" Mn(I/sd):- average < merged<Ih>/sd(<Ih>) > ~= signal/noise\n"+
		" Frcbias :- partial bias = Mean( Mn(If) - Ip )/Mean( Mn(I) )\n"+
		"             for mixed sets only (If is a full if present, else the\n"+
		"             partial with the smallest number of parts)\n\n");

  output.logTab(0,LOGFILE,
		std::string("\n\nBy 4sinTheta/Lambda^2 bins (all statistics use Mn(I+),Mn(I-)etc)\n")+
                	    "----------------------------------------------------------------\n");
  TableGraph table(" Analysis against resolution, "+dataset_pxd.dname());
  Range xrange = ResRange; // x axis range to full resolution limit
  xrange.first() = 0.0;    // from 0
  std::vector<Range> yranges(4);   // for each graph
  // Get y ranges for each graph (if loggraph would accept just an xrange, wouldn't need to do this)
  for (int i=0;i<ResRange.Nbins();++i) {
    float frcbias = 0.0;
    if (biasIRes[i].Count() > 0) {
      frcbias = biasRes[i].Mean()/biasIRes[i].Mean();
    }
    // Graphs are
    // 1. I/sigma, Mean Mn(I)/sd(Mn(I))  (cols 13,14)
    // 2. Rmerge, Rfull, Rmeas, Rpim v Resolution  (cols 4,5,6,7)
    // 3. Average I, RMSdeviation and Sd (cols 10,11,12)
    // 4. Fractional bias (col 15)
    yranges[0].update(imeanRes[i].Mean()/sqrt(rmsDRes[i].Mean()));  // I/sigma
    yranges[0].update(mnIsdRes[i].Mean()); //Mn(I/sd)
    yranges[1].update(rmergeRes[i].R());   // Rmerge
    yranges[1].update(rmergeResFull[i].R());   // Rfull
    yranges[1].update(rmeasRes[i].R());    // Rmeas
    yranges[1].update(rpimRes[i].R());     // Rpim
    yranges[2].update(imeanRes[i].Mean()); // AvI
    yranges[2].update(sqrt(rmsDRes[i].Mean()));  // RMSdeviation
    yranges[2].update(avSdRes[i].Mean());  // Sd
    yranges[3].update(frcbias);
  } // end line loop

  output.logTab(0,LOGFILE,table.formatTitle());
  int c[] = {2,13,14};
  std::vector<int> cln(c,c+3);
  output.logTab(0,LOGFILE,
		table.Graph("I/sigma, Mean Mn(I)/sd(Mn(I))",GraphAxesType(xrange,yranges[0],true),cln));
  int c2[] = {2,4,5,6,7};
  cln.assign(c2, c2+5);
  output.logTab(0,LOGFILE,
		table.Graph("Rmerge, Rfull, Rmeas, Rpim v Resolution",GraphAxesType(xrange,yranges[1],true),cln));
  int c3[] = {2,10,11,12};
  cln.assign(c3, c3+4);
  output.logTab(0,LOGFILE,
		table.Graph("Average I, RMSdeviation and Sd",GraphAxesType(xrange,yranges[2],true),cln));
  int c4[] = {2,15};
  cln.assign(c4, c4+2);
  output.logTab(0,LOGFILE,
		table.Graph("Fractional bias",GraphAxesType(xrange,yranges[3],false),cln));
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
  bool z[] =
    {false, false, false, true, true, true, true, true, false, true, true, true, true, true, true};
  std::vector<bool> Zero(z, z+15);
  std::string fmt = "%7.3f%7.3f%7.3f%7.3f%7.3f%9d%9d%7d%7d%7.1f%9.1f%9.3f\n"; // excluding 1st 3 columns
  output.logTab(0,LOGFILE,
		table.ColumnFields(collabels, Zero, "%3d%8.4f%7.2f"+fmt)); 
  int nc = collabels.size();
  Rfactor Rcum;  // cumulative R
  Rfactor Rfull, Rmeas, Rpim;
  MeanSD Imean, rmsD, avSd, mnIsd, bias, biasI;

  int n=1;
  for (int i=0;i<ResRange.Nbins();++i) {
    Rcum += rmergeRes[i];
    float frcbias = 0.0;
    if (biasIRes[i].Count() > 0) {
      frcbias = biasRes[i].Mean()/biasIRes[i].Mean();
    }
    //    output.logTab(0,LOGFILE,
    //		  table.Line(nc, n++, ResRange.boundsS(i).second,
    //			     ResRange.boundsA(i).second,
    output.logTab(0,LOGFILE,
		  table.Line(nc, n++, ResRange.middle(i),
			     ResRange.middleA(i),
			     rmergeRes[i].R(), rmergeResFull[i].R(), Rcum.R(),
			     rmeasRes[i].R(), rpimRes[i].R(),
			     rmergeRes[i].result().count,
			     Nint(imeanRes[i].Mean()), Nint(sqrt(rmsDRes[i].Mean())),
			     Nint(avSdRes[i].Mean()), imeanRes[i].Mean()/sqrt(rmsDRes[i].Mean()),
			     mnIsdRes[i].Mean(), frcbias));
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
  }
  output.logTab(0,LOGFILE,
		table.CloseTable());
  float frcbias = 0.0;
  if (biasI.Count() > 0) {
    frcbias = bias.Mean()/biasI.Mean();
  }
  fmt = "Overall:          "+fmt;
  // //  fmt = "Overall:          "+fmt+"\n";
  output.logTabPrintf(0,LOGFILE,fmt.c_str(),
		      Rcum.R(), Rcum.R(), Rfull.R(),
		      Rmeas.R(), Rpim.R(), Rcum.result().count,
		   Nint(Imean.Mean()), Nint(sqrt(rmsD.Mean())), Nint(avSd.Mean()),
		   Imean.Mean()/sqrt(rmsD.Mean()), mnIsd.Mean(), frcbias);
  output.logTab(0,LOGFILE,table.RawLabels());
  // Store things in summary object
  summarystatistics.StoreRmergeReso(Rcum, rmergeRes[0], rmergeRes[ResRange.Nbins()-1]);
  summarystatistics.StoreRmeasReso(Rmeas, rmeasRes[0], rmeasRes[ResRange.Nbins()-1]);
  summarystatistics.StoreRpimReso(Rpim, rpimRes[0], rpimRes[ResRange.Nbins()-1]);
  summarystatistics.StoreMnIsd(mnIsd.Mean(), mnIsdRes[0].Mean(),
			       mnIsdRes[ResRange.Nbins()-1].Mean());
  // Resolution "limit" from Mn(I/sd)
  summarystatistics.StoreMnIsigresolimit
    (ResolutionLimit(mnIsdRes, ResRange, MinimumIoverSigma));
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
                	        "------------------------------------------------------------------\n");
  output.logTab(0,LOGFILE,
		"\nStatistics labelled 'Ov' are relative to the overall mean I+/-, ignoring anomalous");
  output.logTab(0,LOGFILE,
		"Other statistics are with either I+ or I- sets, for acentrics, ie with anomalous\n\n");
  TableGraph table(" Analysis against resolution, with & without anomalous (Ov), "+dataset_pxd.dname());

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

  output.logTab(0,LOGFILE,table.formatTitle());
  int c[] = {2,4,5,8,9,10,11};
  std::vector<int> cln(c,c+7);
  output.logTab(0,LOGFILE,
		table.Graph("Rmerge, Rmeas, Rpim v Resolution",GraphAxesType(xrange,yrange,true),cln));
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
  output.logTab(0,LOGFILE,
		table.ColumnFields(collabels, Zero, "%3d%8.4f%7.2f"+fmt)); 
  int nc = collabels.size();
  Rfactor Rcum, RcumOv;  // cumulative R
  Rfactor Rmeas, Rpim, RmeasOv, RpimOv;

  int n=1;
  for (int i=0;i<ResRange.Nbins();++i) {
    Rcum += rmergeRes[i];
    RcumOv += rmergeResOv[i];
    //  output.logTab(0,LOGFILE,
    //		table.Line(nc, n++, ResRange.boundsS(i).second, ResRange.boundsA(i).second,
  output.logTab(0,LOGFILE,
		table.Line(nc, n++, ResRange.middle(i), ResRange.middleA(i), // 1,2,3
			   rmergeRes[i].R(), rmergeResOv[i].R(), // 4,5
			   Rcum.R(), RcumOv.R(), // 6,7
			   rmeasRes[i].R(), rmeasResOv[i].R(),  //8,9
			   rpimRes[i].R(), rpimResOv[i].R(),    //10,11
			   rmergeResOv[i].result().count));     //12
    // Totals
    Rmeas += rmeasRes[i];
    Rpim  += rpimRes[i];
    RmeasOv += rmeasResOv[i];
    RpimOv  += rpimResOv[i];
  }
  output.logTab(0,LOGFILE,
		table.CloseTable());
  // //   fmt = "Overall:          "+fmt+"\n";
  fmt = "Overall:          "+fmt;
  output.logTabPrintf(0,LOGFILE,fmt.c_str(),
		      Rcum.R(), RcumOv.R(),
		      Rcum.R(), RcumOv.R(),
		      Rmeas.R(), RmeasOv.R(),
		      Rpim.R(), RpimOv.R(),
		      Rcum.result().count);
  output.logTab(0,LOGFILE,table.RawLabels());
  // Store things in summary object
  summarystatistics.StoreRmergeResoOv(RcumOv, rmergeResOv[0], rmergeResOv[ResRange.Nbins()-1]);
  summarystatistics.StoreRmeasResoOv(RmeasOv, rmeasResOv[0], rmeasResOv[ResRange.Nbins()-1]);
  summarystatistics.StoreRpimResoOv(RpimOv, rpimResOv[0], rpimResOv[ResRange.Nbins()-1]);
}
//--------------------------------------------------------------
void PrintDeviationsByIntensity(const PxdName& dataset_pxd,
				const IntensityBin& Irange,
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
		"-----------------\n");
  TableGraph table(" Analysis against intensity, "+dataset_pxd.dname());
  output.logTab(0,LOGFILE,table.formatTitle());
  int c[] = {1,2,4,5};
  std::vector<int> cln(c,c+4);
  output.logTab(0,LOGFILE,
		table.Graph("Rmerge v Intensity","N",cln));
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
  output.logTab(0,LOGFILE,
		table.ColumnFields(collabels, Zero, "%10.0f"+fmt)); 
  int nc = collabels.size();
  Rfactor Rcum;  // cumulative R
  Rfactor Rmeas, Rpim;
  MeanSD Imean, rmsD, avSd, mnIsd, bias, biasI;

  for (int i=0;i<Irange.NumberBins();++i) {
    Rcum += rmergeInt[i];
    float frcbias = 0.0;
    if (biasIInt[i].Count() > 0) {
      frcbias = biasInt[i].Mean()/biasIInt[i].Mean();
    }
  output.logTab(0,LOGFILE,
		table.Line(nc, Irange.bounds(i).second,
			   rmergeInt[i].R(), Rcum.R(),
			   rmeasInt[i].R(), rpimInt[i].R(),
			   rmergeInt[i].result().count,
			   Nint(imeanInt[i].Mean()), Nint(sqrt(rmsDInt[i].Mean())),
			   Nint(avSdInt[i].Mean()), imeanInt[i].Mean()/sqrt(rmsDInt[i].Mean()),
			   mnIsdInt[i].Mean(), frcbias));
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
  output.logTab(0,LOGFILE,
		table.CloseTable());
  float frcbias = 0.0;
  if (biasI.Count() > 0) {
    frcbias = bias.Mean()/biasI.Mean();
  }
  fmt = "Overall:  "+fmt+"\n";
  output.logTabPrintf(0,LOGFILE,fmt.c_str(),
		   Rcum.R(), Rcum.R(),
		      Rmeas.R(), Rpim.R(), Rcum.result().count,
		   Nint(Imean.Mean()), Nint(sqrt(rmsD.Mean())), Nint(avSd.Mean()),
		   Imean.Mean()/sqrt(rmsD.Mean()), mnIsd.Mean(), frcbias);
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
		"  %poss is completeness in the shell, C%poss in cumulative to that resolution\n"+
		" The anomalous completeness values (AnomCmpl) are the percentage of possible anomalous "+
		"differences measured\n"+
		" AnomFrc is the % of measured acentric reflections for which an anomalous "+
		"difference has been measured\n\n");

  // Get number of reflections in each resolution bin in complete sphere
  int nbins = ResRange.Nbins();
  std::vector<int> Nrefres(nbins,0);  // number
  std::vector<int> Nrefacen(nbins,0);  // number of acentrics
  // Get numbers in each resolution shell
  NumberComplete(ResRange, symmetry, cell, Nrefres, Nrefacen);

  Range xrange = ResRange; // x axis range to full resolution limit
  xrange.first() = 0.0;    // from 0
  std::vector<Range> yranges(2);   // for each graph
  // Get y ranges for each graph (if loggraph would accept just an xrange, wouldn't nned to do this)
  for (int i=0;i<ResRange.Nbins();++i) {
    yranges[1].update(FractionN(1.0, NumObs[i], NumRef[i]));
    yranges[1].update(FractionN(1.0, SNumAnomPairs[i], NumACentric[i]));
  }
  yranges[0].first() = 0.0;  // maximum completeness 
  yranges[0].last() = 100.0;  // maximum completeness 
  TableGraph table(" Completeness & multiplicity v. resolution, "+dataset_pxd.dname());
  output.logTab(0,LOGFILE,table.formatTitle());
  int c[] = {2,7,8,10,11};
  std::vector<int> cln(c,c+5);
  output.logTab(0,LOGFILE,
		table.Graph("Completeness v Resolution ",GraphAxesType(xrange,yranges[0],true),cln));
  int c2[] = {2,9,12};
  cln.assign(c2, c2+3);
  output.logTab(0,LOGFILE,
		table.Graph("Multiplicity v Resolution",GraphAxesType(xrange,yranges[1],true),cln));
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
  output.logTab(0,LOGFILE,
		table.ColumnFields(collabels, Zero, "%3d%8.4f%7.2f"+fmt)); 
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
    //    output.logTab(0,LOGFILE,
    //		table.Line(nc, n++, ResRange.boundsS(i).second, ResRange.boundsA(i).second,
    output.logTab(0,LOGFILE,
		table.Line(nc, n++, ResRange.middle(i), ResRange.middleA(i),
			   NumObs[i], NumRef[i], NumCentric[i], poss, cumposs,
			   FractionN(1.0, NumObs[i], NumRef[i]),
			   anomcmpl, anomfrc, FractionN(1.0, SNumAnomPairs[i], NumACentric[i])));
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
      anommult0 = FractionN(1.0, SNumAnomPairs[i], NumACentric[i]);
    } else if (i==ResRange.Nbins()-1) {
      possN = poss;
      multN = FractionN(1.0, NumObs[i], NumRef[i]);
      anomcmplN = anomcmpl;
      anommultN = FractionN(1.0, SNumAnomPairs[i], NumACentric[i]);
    }

  }
  output.logTab(0,LOGFILE,
		table.CloseTable());

    anomcmpl = FractionN(100., NanomSphere, Nanomref);
  std::string leader = "Overall:          ";
  fmt = leader+fmt;
  output.logTabPrintf(0,LOGFILE,fmt.c_str(),
		      Nobs, Nref, Ncen, cumposs, cumposs, FractionN(1.0, Nobs, Nref),
		      anomcmpl, FractionN(100., Nanom, Nacen), FractionN(1.0, SNumanompairs, Nacen));
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
{
  ASSERT (ResRange.Nbins() == halfDatasetScores.NresBin());

  output.logTab(0,LOGFILE,
    std::string("\n\nCorrelation coefficients for anomalous differences & Imean between random half-datasets\n")+
                    "=======================================================================================\n\n"+
" The RMS Correlation Ratio (RCR) is calculated from a scatter plot of pairs of DeltaI(anom)\n"+
" from the two subsets by comparing the RMS value (excluding extremes) projected on the line \n"+
" with slope = 1 ('correlation') with the RMS value perpendicular to this ('error').\n"+
" This ratio will be > 1 if there is an anomalous signal\n");

  TableGraph table(" Correlations within dataset, "+dataset_pxd.dname());
  output.logTab(0,LOGFILE,table.formatTitle());

  Range xrange = ResRange; // x axis range to full resolution limit
  xrange.first() = 0.0;    // from 0
  std::vector<Range> yranges(2);   // for each graph
  // Get y ranges for each graph (if loggraph would accept just an xrange, wouldn't nned to do this)
  for (int i=0;i<ResRange.Nbins();++i) {
    yranges[1].update(halfDatasetScores.RMScorrelRatio(i));
    yranges[1].update(halfDatasetScores.RMScorrelRatioCen(i));
  }
  yranges[0].first() = 0.0;  // CC
  yranges[0].last() = 1.0;  // CC
  int c[] = {2,4,6,10};
  std::vector<int> cln(c,c+4);
  output.logTab(0,LOGFILE,
		table.Graph(" Anom & Imean CCs v resolution - ",
		GraphAxesType(xrange,yranges[0],true),cln));
  int c2[] = {2,8,9};
  cln.assign(c2,c2+3);
  output.logTab(0,LOGFILE,
		table.Graph(" RMS correlation ratio ",GraphAxesType(xrange,yranges[1],true),cln));
  std::vector<std::string> collabels;
  collabels.push_back("N");         // 1
  collabels.push_back("1/d^2");     // 2
  collabels.push_back("Dmid");      // 3
  collabels.push_back("CCanom");     // 4
  collabels.push_back("Nanom");      // 5
  collabels.push_back("CCcen");     // 6
  collabels.push_back("Ncen");     // 7
  collabels.push_back("RCRanom");    // 8
  collabels.push_back("RCRcen");    // 9
  collabels.push_back("CCImean");    // 10
  collabels.push_back("NImean");    // 11
  bool z[] = {false, false, false, true, false, true, false, true, true, true, false};
  std::vector<bool> Zero(z, z+11);
  std::string fmt = "%7.3f%9d%7.3f%9d   %7.3f%7.3f %7.3f%9d\n"; // excluding 1st 3 columns
  output.logTab(0,LOGFILE,
		table.ColumnFields(collabels, Zero, "%3d%8.4f%7.2f"+fmt)); 
  int nc = collabels.size();

  int n=1;
  for (int i=0;i<ResRange.Nbins();++i) {
    //  output.logTab(0,LOGFILE,
    //		table.Line(nc, n++, ResRange.boundsS(i).second, ResRange.boundsA(i).second,
  output.logTab(0,LOGFILE,
		table.Line(nc, n++, ResRange.middle(i), ResRange.middleA(i),
			   halfDatasetScores.CCanom(i).result().val,
			   halfDatasetScores.CCanom(i).result().count,
			   halfDatasetScores.CCanomCen(i).result().val,
			   halfDatasetScores.CCanomCen(i).result().count,
			   halfDatasetScores.RMScorrelRatio(i),
			   halfDatasetScores.RMScorrelRatioCen(i),
			   halfDatasetScores.CC_Imean(i).result().val,
			   halfDatasetScores.CC_Imean(i).result().count));
  }
  output.logTab(0,LOGFILE,
		table.CloseTable());
  // Totals
  std::string leader = "Overall:          ";
  fmt = leader+fmt;
  output.logTabPrintf(0,LOGFILE,fmt.c_str(),
	       halfDatasetScores.CCanom().result().val,
	       halfDatasetScores.CCanom().result().count,
	       halfDatasetScores.CCanomCen().result().val,
	       halfDatasetScores.CCanomCen().result().count,
	       halfDatasetScores.RMScorrelRatio(),
	       halfDatasetScores.RMScorrelRatioCen(),
	       halfDatasetScores.CC_Imean().result().val,
	       halfDatasetScores.CC_Imean().result().count);
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
}
//--------------------------------------------------------------
void PrintAnisotropyAnalysis(const PxdName& dataset_pxd,
			     const ResoRange& ResRange,
			     const HalfDataset& halfDatasetScores,
			     const std::vector<std::vector<MeanSD> >& mnIsdResCone,
			     const double& coneangledegrees,
			     const double& MinimumIoverSigma,
			     SummaryStatistics& summarystatistics,
			     phaser_io::Output& output)
{
  ASSERT (ResRange.Nbins() == halfDatasetScores.NresBin());

  std::string s = std::string("\n\nAnalysis of anisotropy of data\n")+
    "==============================\n\n"+
    "Mn(I/sd) and half-dataset correlation coefficients are analysed by resolution\n"+
    "within cones around the three reciprocal lattice axes,"+
    " cone half-angle = %5.1f degrees\n";
  output.logTabPrintf(0,LOGFILE, s.c_str(), coneangledegrees);

  Range xrange = ResRange; // x axis range to full resolution limit
  xrange.first() = 0.0;    // from 0
  std::vector<Range> yranges(2);   // for each graph
  for (int i=0;i<ResRange.Nbins();++i) {
    yranges[1].update(mnIsdResCone[0][i].Mean());
    yranges[1].update(mnIsdResCone[1][i].Mean());
    yranges[1].update(mnIsdResCone[2][i].Mean());
  }
  yranges[0].first() = 0.0; // CC
  yranges[0].last() = 1.0; // CC
  TableGraph table(" Anisotropy analysis, "+dataset_pxd.dname());
  output.logTab(0,LOGFILE,table.formatTitle());
  int c[] = {2,4,5,6};
  std::vector<int> cln(c,c+4);
  output.logTab(0,LOGFILE,
		table.Graph(" Imean CCs v resolution - ",GraphAxesType(xrange,yranges[0],true),cln));
  int c2[] = {2,7,8,9};
  cln.assign(c2,c2+4);
  output.logTab(0,LOGFILE,
		table.Graph(" Mn(I/sd) v resolution - ",GraphAxesType(xrange,yranges[1],true),cln));
  std::vector<std::string> collabels;
  collabels.push_back("N");         // 1
  collabels.push_back("1/d^2");     // 2
  collabels.push_back("Dmid");      // 3
  collabels.push_back("CC_a*");     // 4
  collabels.push_back("CC_b*");     // 5
  collabels.push_back("CC_c*");     // 6
  collabels.push_back("Mn(I/sd)a*");     // 7
  collabels.push_back("Mn(I/sd)b*");     // 8
  collabels.push_back("Mn(I/sd)c*");     // 9
  bool z[] = {false, false, false, true, true, true, true, true, true};
  std::vector<bool> Zero(z, z+9);
  std::string fmt = " %8.3f%8.3f%8.3f %11.2f%11.2f%11.2f\n"; // excluding 1st 3 columns
  output.logTab(0,LOGFILE,
		table.ColumnFields(collabels, Zero, "%3d%8.4f%7.2f"+fmt)); 
  int nc = collabels.size();

  std::vector<MeanSD> mnIsd(3);
  int n=1;
  for (int i=0;i<ResRange.Nbins();++i) {
    //    output.logTab(0,LOGFILE,
    //		  table.Line(nc, n++, ResRange.boundsS(i).second, ResRange.boundsA(i).second,
    output.logTab(0,LOGFILE,
		  table.Line(nc, n++, ResRange.middle(i), ResRange.middleA(i),
			     halfDatasetScores.CCaniso(0,i).result().val,
			     halfDatasetScores.CCaniso(1,i).result().val,
			     halfDatasetScores.CCaniso(2,i).result().val,
			     mnIsdResCone[0][i].Mean(),
			     mnIsdResCone[1][i].Mean(),
			     mnIsdResCone[2][i].Mean()));
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
  output.logTab(0,LOGFILE,
		table.CloseTable());
  // Totals
  std::string leader = "Overall:          ";
  fmt = leader+fmt;
  output.logTabPrintf(0,LOGFILE,fmt.c_str(),
		      halfDatasetScores.CCaniso(0).result().val,
		      halfDatasetScores.CCaniso(1).result().val,
		      halfDatasetScores.CCaniso(2).result().val,
		      mnIsd[0].Mean(),
		      mnIsd[1].Mean(),
		      mnIsd[2].Mean());

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
    // find resolution "limit" where Mn(I/sd) falls below MinimumIoverSigma
    reslimisig[i].init(mnIsdResCone[i], ResRange, MinimumIoverSigma);
  }
  summarystatistics.StoreMnIsigAnisoresolimit(reslimisig);
}
//--------------------------------------------------------------
void PrintUnmergedHeaderStuff(const scala::hkl_unmerge_list& hkl_list,
			      phaser_io::Output& output,
			      const int& verbose)
// Optional summary printing
{
  if (verbose > 0)  {
    if (verbose == 1 || !hkl_list.IsReady()) {
      output.logTab(0,LOGFILE,
		    "      ResolutionRange    NobsParts  Nbatches  Ndatasets");
      output.logTabPrintf(0,LOGFILE,
			  "     %8.2f %6.2f  %10d%7d%10d\n",
			  hkl_list.ResRange().ResLow(), hkl_list.ResRange().ResHigh(),
			  hkl_list.num_parts(), hkl_list.num_batches(),
			  hkl_list.num_datasets());
      
    } else {
      //        output.logTabPrintf(0,LOGFILE,
      //			    "\nSummary of reflection list\n");
      output.logTabPrintf(0,LOGFILE,
			  "\n   Resolution range accepted: %8.2f    %8.2f\n",
			  hkl_list.ResRange().ResLow(),
			  hkl_list.ResRange().ResHigh());
      
      output.logTabPrintf(0,LOGFILE,
			  "\n   Number of reflections  =    %10d\n",
			  hkl_list.num_reflections_valid());
      output.logTabPrintf(0,LOGFILE,
			  "   Number of observations =    %10d\n",
			  hkl_list.num_observations());
      output.logTabPrintf(0,LOGFILE,
			  "   Number of parts        =    %10d\n",
			  hkl_list.num_parts());
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
      }
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
      std::vector<Xdataset> datasets = hkl_list.AllXdatasets();
      std::vector<Run> runlist = hkl_list.RunList();
      for (int k=0; k<ndatasets; k++) {
	output.logTab(0,LOGFILE,"\n"+datasets[k].pxdname().formatPrint());
	output.logTab(3,LOGFILE,"Cell: "+datasets[k].cell().formatPrint());
	output.logTabPrintf(3,LOGFILE,
			    "Wavelength %8.5f A\n", datasets[k].wavelength());
	if (verbose == 3) {
	  for (size_t i=0;i<runlist.size();i++) {
	    if (runlist[i].DatasetIndex() == k) {
	      output.logTab(0,LOGFILE,runlist[i].formatPrintBrief(datasets));
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
    } else if (verbose > 3) {
      std::vector<Xdataset> datasets = hkl_list.AllXdatasets();
      for (size_t k=0; k<datasets.size(); k++) {
	output.logTab(0,LOGFILE,
		      datasets[k].formatPrint());
      }
    }
    output.logTabPrintf(1,LOGFILE, "Average unit cell: ");
    output.logTab(0,LOGFILE,hkl_list.Cell().formatPrint());
    output.logTabPrintf(1,LOGFILE,"");
      
    if (verbose > 3) {
      std::vector<Run> runlist = hkl_list.RunList();
      std::vector<Xdataset> datasets = hkl_list.AllXdatasets();
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
    output.logTabPrintf(1,LXML,
			"<NumberReflections>  %10d </NumberReflections>\n",
			hkl_list.num_reflections_valid());
    output.logTabPrintf(1,LXML,
			"<NumberObservations> %10d </NumberObservations>\n",
			hkl_list.num_observations());
    output.logTabPrintf(1,LXML,
			"<NumberParts>        %10d </NumberParts>\n",
			hkl_list.num_parts());
    output.logTabPrintf(1,LXML,
			"<NumberBatches>      %10d </NumberBatches>\n", 
			hkl_list.num_batches());
    output.logTabPrintf(1,LXML,
			"<NumberDatasets>     %10d </NumberDatasets>\n",
			hkl_list.num_datasets());

    int ndatasets = hkl_list.num_datasets();
    std::vector<Xdataset> datasets = hkl_list.AllXdatasets();
    for (int k=0; k<ndatasets; k++) {
      output.logTabPrintf(1,LXML, "<Dataset  name=\"%s\">\n",
			  datasets[k].pxdname().format().c_str());
      for (size_t i=0;i<runlist.size();i++) {
	if (runlist[i].DatasetIndex() == k) {
	  output.logTabPrintf(2,LXML,"<Run> <number> %3d </number>\n",i+1);
	  output.logTabPrintf(2,LXML,
			      "<BatchRange> %8d %8d </BatchRange>\n",
			      runlist[i].BatchRange().first, runlist[i].BatchRange().second);
	  output.logTabPrintf(2,LXML,"<BatchOffset> %8d </BatchOffset>\n",
			      runlist[i].BatchNumberOffset());
	  if (!OneFile) {
	    hklstream =
	      StringUtil::Strip("HKLIN"+clipper::String(runlist[i].FileNumber()));
	  }
	  output.logTab(2,LXML,"<FileStream> "+hklstream+" </FileStream>");
	  output.logTabPrintf(2,LXML,"</Run>\n",i+1);
	}
      }
      output.logTabPrintf(1,LXML, "</Dataset>\n");
    }
    output.logTab(0,LXML,"</ReflectionData>");
  }
}
//--------------------------------------------------------------
