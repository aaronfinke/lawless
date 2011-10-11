// writeoutputfiles.cpp
//
// Write output reflection files as specified by outputcontrols
//

#include "writeoutputfiles.hh"
#include "mergedlist.hh"
#include "scala_util.hh"
#include "writeunmerged.hh"

namespace scala {
// ---------------------------------------------------------
  void WriteMergedOutputFiles(const std::string& title,
			const MergedList& mergedlist, 
			const OutputControls& outputcontrols,
			phaser_io::Output& output)
  // Merged file output
  {
    // MTZ output
    if (outputcontrols.MTZoutputMerged()) {
      WriteMergedOutput(+1, title, mergedlist, outputcontrols, output);
    }
    // SCA output
    if (outputcontrols.SCAoutputMerged()) {
      WriteMergedOutput(+2, title, mergedlist, outputcontrols, output);
    }
  }
// ---------------------------------------------------------
  void WriteUnmergedOutputFiles(const std::string& title,
				const hkl_unmerge_list& hkl_list,
				const SDmodel& SDM,
				const float& Imax,
				const OutputControls& outputcontrols,
				phaser_io::Output& output)
  // Unmerged file output
  {
    // MTZ output
    if (outputcontrols.MTZoutputUnmerged()) {
      bool summedpartials = true;  // observations not parts (ie summed partials)
      int datasetindex = -1;  // all datasets
      std::string filename = outputcontrols.Mtzunmergedfilename("");
      if (datasetindex < 0) {
	output.logTab(0, LOGFILE,
		      "\n==== Writing unmerged data for all datasets to file "+
		      filename);
      } else {
	PxdName pxdname = hkl_list.xdataset(datasetindex).pxdname();
	output.logTab(0, LOGFILE,
		      "\n==== Writing unmerged data for dataset "+pxdname.format()+
		    " to file "+filename+"\n");
      }
      int nref =
	MtzIO::WriteUnmergedMTZ(hkl_list, SDM, summedpartials, datasetindex,
				filename, title);
      output.logTabPrintf(0, LOGFILE,
			  "\nNumber of observations written = %8d\n", nref);
    }

    // SCA output
    if (outputcontrols.SCAoutputUnmerged()) {
      for (int idts=0;idts<hkl_list.num_datasets();++idts) {  // loop datasets
	PxdName pxdname = hkl_list.xdataset(idts).pxdname();
	std::string filedname = pxdname.dname(); // append to filename if > 1 dataset
	if (hkl_list.num_datasets() <= 1) {
	  filedname = "";
	}
	std::string filename = outputcontrols.Scaunmergedfilename(filedname);
	output.logTab(0, LOGFILE,
		      "\n==== Writing unmerged data for dataset "+pxdname.format()+
		      " to file "+filename+"\n");
	int nref =
	  MtzIO::WriteUnmergedSCA(hkl_list, SDM, idts, filename, Imax);
	output.logTabPrintf(0, LOGFILE,
			    "\nNumber of observations written = %8d\n", nref);
      } // end loop datasets
    } // SCA
  }
// ---------------------------------------------------------
  void WriteMergedOutput(const int& outputformat,
			 const std::string& title,
			 const MergedList & mergedlist,
			 const OutputControls& outputcontrols,
			 phaser_io::Output& output)
  // Create file containing columns
  //  h,k,l
  //  IMEAN, SIGIMEAN
  //  I(+), SIGI(+), I(-), SIGI(-)
  //
  //  outputformat  = +1  MTZ, +2 Scalepack
  {
    ASSERT (outputformat >= +1 && outputformat <= +2);
    int nfiles = mergedlist.NumberDatasets();  // one file/dataset
    // combine datasets together: only for MTZ:  NOT WRITTEN YET!
    if (!outputcontrols.Split() && outputformat == +1) nfiles = 1;

    // File name(s)
    std::vector<Xdataset> xdatasets = mergedlist.Datasets();

    for (int idts=0;idts<mergedlist.NumberDatasets();++idts) {
      PxdName pxdname = xdatasets[idts].pxdname();
      std::string filedname = "";
      if (nfiles > 1) {
	filedname = pxdname.dname();  // multiple files identified by dataset name
      }
      std::string filename;

      if (outputformat == +2) {  // Scalepack format
	filename = outputcontrols.Scamergedfilename(filedname);
      } else if (outputformat == +1) { // MTZ format
	filename = outputcontrols.Mtzmergedfilename(filedname);
      }
      output.logTab(0, LOGFILE,
		    "\n==== Writing merged data for dataset "+pxdname.format()+
		    " to file "+filename);

      int nref;
      if (outputformat == +1) {
	nref = mergedlist.WriteDatasetToMTZ(filename, idts);
      } else if (outputformat == +2) {
	nref = mergedlist.WriteDatasetToSCA(filename, idts);
      }

      float rmax = mergedlist.InvResMax(idts);
      if (rmax > 0.0) rmax = 1.0/sqrt(rmax);
      output.logTabPrintf(0, LOGFILE,
	  "\nNumber of reflections written %8d maximum resolution %8.3f\n",
			  nref, rmax);
    }
  }
} // namespace scala
