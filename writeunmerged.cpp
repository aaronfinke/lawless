// writeunmerged.cpp
//


// Clipper
#include "clipper/clipper.h"
#include "clipper/clipper-ccp4.h"
using clipper::Message;
using clipper::Message_fatal;

#define ASSERT assert
#include <assert.h>


// CCP4
#include "ccp4/csymlib.h"    // CCP4 symmetry stuff
#include "ccp4/ccp4_general.h"
#include "ccp4/cmtzlib.h"    // CCP4 MTZlib headers (namespace CMtz)

#include "hkl_datatypes.hh"
#include "hkl_controls.hh"
#include "scala_util.hh"
#include "mtz_utils.hh"
#include "hkl_symmetry.hh"
#include "writeunmerged.hh"
#include "openinputfile.hh"
#include "observationflags.hh"
#include "cellgroup.hh"
#include "file_util.hh"
#include "string_util.hh"

using namespace CMtz;

namespace MtzIO
{
  WriteUnmerged::WriteUnmerged()
  {

  }
  //--------------------------------------------------------------
  void WriteUnmerged::OptAddCol(const bool& coln,
				MTZCOL* col[], int& ic, MTZ* mtzout, MTZSET* baseset,
				const char* label, const char* type)
    // Conditional column addition, only if coln > 0    
  {
    if (coln)
	col[ic++] = MtzAddColumn(mtzout, baseset, label, type);
  }
  //--------------------------------------------------------------
  int WriteUnmerged::writeUnmergedMTZ(const scala::hkl_unmerge_list& hkl_list,
				      const SDmodel& SDM,
				      const bool& summedpartials,
				      const int& datasetIndex,
				      const std::string& filename_out,
				      const std::string& title)
  // Write unmerged MTZ file from hkl_unmerge_list object
  // returns number written
  //
  //  summedpartials   true to output scaled observations with summed partials
  //                   false to output original data
  //  datasetIndex     dataset index to output, -1 all data
  //  filename_out     output filename
  //  title
  {
    nref = 0;
    nmultiple = 0;
    // Initialise MTZ data structure and output file
    MTZ* mtzout = MtzMalloc(0,0);
    mtzout->refs_in_memory = 0;   // not in memory
    mtzout->fileout = MtzOpenForWrite(filename_out.c_str());
    if (mtzout->fileout == NULL)
      {Message::message(Message_fatal("Can't open file "+filename_out));}
    int nlattices = hkl_list.NumberofLattices();

    // Write title
    ccp4_lwtitl(mtzout, title.c_str(), 0);

    // Symmetry
    hkl_symmetry NewSymm = hkl_list.symmetry();
    // Check for rhombohedral lattice (H or R)
    char HorR = 'H'; // default H setting
    if (NewSymm.lattice_type() == 'H' || NewSymm.lattice_type() == 'R') {
      if (RhombohedralAxes(hkl_list.Cell().UnitCell())) { // true if not H
	HorR = 'R';
      }}
    char spg_status = hkl_list.MtzSym().spg_confidence;
    mtzout->mtzsymm = spg_to_mtz(NewSymm.GetSpaceGroup(), HorR, spg_status);
    // for cell constraints
    CCtbxSym::CellGroup CG(NewSymm.GetSpaceGroup());

    // Crystals & datasets
    int Nxd = hkl_list.num_datasets();
    std::vector<std::string> xtl;
    float ucell[6];
    std::string lastxname = "";
    MTZXTAL* xtal;
    MTZSET* set;

    // Put in HKL_base dataset explicitly
    scala::Scell HKLcell = CG.constrain(hkl_list.cell());  // constrained cell
    for (int i=0;i<6;i++)
      {ucell[i] = HKLcell[i];}
    xtal = MtzAddXtal(mtzout, "HKL_base", "HKL_base",
		      ucell);
    MTZSET* baseset = MtzAddDataset(mtzout, xtal, "HKL_base", 0.0);

    int j1=0;
    int j2=Nxd;
    if (datasetIndex >= 0) {
      j1=datasetIndex;
      j2=datasetIndex+1;
    }

    std::map<std::string, MTZXTAL*> xtals;     // list of crystals indexed by Xname

    for (int jxd=j1; jxd<j2;jxd++) {
      // loop datasets
      //   Note that the MTZ library insists that a crystal contains datasets,
      //   rather than the other way round, and that the Cell is a property of the crystal
      Dataset dataset = hkl_list.dataset(jxd);
      std::vector<PxdName> pxdnames = dataset.pxdnames(); // names for this dataset (maybe only one)
      for (size_t ixt=0;ixt<pxdnames.size();++ixt) { // loop crystals
	std::string xname = pxdnames[ixt].xname();
	if (xtals[xname] == 0) { // new crystal
	  HKLcell = CG.constrain(hkl_list.dataset(jxd).cell());
	  for (int i=0;i<6;i++) {
	    ucell[i] = HKLcell[i];
	  }
	  xtal = MtzAddXtal(mtzout, xname.c_str(),
			    hkl_list.dataset(jxd).Pname().c_str(),
			    ucell);
	  xtals[xname] = xtal; // store pointer
	} else { // this dataset belongs to a crystal we have already
	  xtal = xtals[xname];
	}
	// Add this Xdataset
	set = MtzAddDataset(mtzout, xtal,
			    hkl_list.dataset(jxd).Dname().c_str(),
			    hkl_list.dataset(jxd).wavelength(xname));
      } // end loop Xdatasets (ie crystals)
    } // end loop datasets

    // Columns, all in base dataset
    data_flags  col_sel = hkl_list.DataFlags();
    //    col_sel.print(); //^ Debug
    int maxhkloverlap = hkl_list.MaxHKLoverlap();
    ASSERT (maxhkloverlap == col_sel.n_latinfo);

    // How many columns?
    MTZCOL* col[MAXNCOLUMNS];  // MAXNCOLUMNS defined in openinputfile.hh

    int ic=0;
    col[ic++] = MtzAddColumn(mtzout, baseset, "H", "H");
    col[ic++] = MtzAddColumn(mtzout, baseset, "K", "H");
    col[ic++] = MtzAddColumn(mtzout, baseset, "L", "H");
    col[ic++] = MtzAddColumn(mtzout, baseset, "M/ISYM", "Y");
    col[ic++] = MtzAddColumn(mtzout, baseset, "BATCH", "B");
    col[ic++] = MtzAddColumn(mtzout, baseset, "I", "J");
    col[ic++] = MtzAddColumn(mtzout, baseset, "SIGI", "Q");

    if (summedpartials) {
      OptAddCol(true, col, ic, mtzout, baseset, "SCALEUSED", "R");
      OptAddCol(true, col, ic, mtzout, baseset, "SIGSCALEUSED", "R");
      // output of observations, write NPART as number of images used for observation
      OptAddCol(true, col, ic, mtzout, baseset, "NPART", "I");
    }
    // Optional columns, add in only if present (ie read from input file)
    // also increment ic
    if (!summedpartials) {  // only for raw output
      OptAddCol(col_sel.is_Ipr, col, ic, mtzout, baseset, "IPR", "J");
      OptAddCol(col_sel.is_sigIpr, col, ic, mtzout, baseset, "SIGIPR", "Q");
    }
    OptAddCol(col_sel.is_fractioncalc, col, ic, mtzout, baseset, "FRACTIONCALC", "R");
    OptAddCol(col_sel.is_Xdet, col, ic, mtzout, baseset, "XDET", "R");
    OptAddCol(col_sel.is_Ydet, col, ic, mtzout, baseset, "YDET", "R");
    OptAddCol(col_sel.is_Rot, col, ic, mtzout, baseset, "ROT", "R");
    OptAddCol(col_sel.is_Width, col, ic, mtzout, baseset, "WIDTH", "R");
    OptAddCol(col_sel.is_LP, col, ic, mtzout, baseset, "LP", "R");
    if (!summedpartials) {  // only for raw output
      // MPART only if in input file
      OptAddCol(col_sel.is_Mpart, col, ic, mtzout, baseset, "MPART", "I");
      OptAddCol(col_sel.is_ObsFlag, col, ic, mtzout, baseset, "FLAG", "I");
      OptAddCol(col_sel.is_BgPkRatio, col, ic, mtzout, baseset, "BGPKRATIOS", "R");
    }
    OptAddCol(col_sel.is_time, col, ic, mtzout, baseset, "TIME", "R");
    if (!summedpartials) {  // only for raw output
      OptAddCol(col_sel.is_scale, col, ic, mtzout, baseset, "SCALE", "R");
      OptAddCol(col_sel.is_sigscale, col, ic, mtzout, baseset, "SIGSCALE", "R");
    }
    if (nlattices > 1) { // multilattice
      OptAddCol(true, col, ic, mtzout, baseset, "LATTNUM", "I");
      for (int l=0;l<maxhkloverlap;++l) {
	std::string ns = StringUtil::Strip(StringUtil::itos(l+1,3));
	OptAddCol(true, col, ic, mtzout, baseset, ("LATTNUM"+ns).c_str(), "I");
	OptAddCol(true, col, ic, mtzout, baseset, ("H"+ns).c_str(), "H");
	OptAddCol(true, col, ic, mtzout, baseset, ("K"+ns).c_str(), "H");
	OptAddCol(true, col, ic, mtzout, baseset, ("L"+ns).c_str(), "H");
	OptAddCol(true, col, ic, mtzout, baseset, ("SCALE"+ns).c_str(), "R");
      }
    }
    int NumCol = ic;
    
    // List is sorted on the first 5 columns
    MtzSetSortOrder(mtzout, col);
    
    // History    
    char history[MTZRECORDLENGTH];
    char date[11];
    char time[9];
    CCP4::ccp4_utils_date(date);
    CCP4::ccp4_utils_time(time);
    std::string text = "AIMLESS, "+std::string(date)+" "+
      std::string(time);
    strcpy(history, text.c_str()); 
    int Nhist = MtzAddHistory(mtzout, &history, 1);
    Nhist = Nhist;

    std::vector<int> nobsbatch;
    hkl_list.rewind();

    if (summedpartials) {
      // write out summed observations, for selected dataset(s)
      nobsbatch =  writeObservations(hkl_list, SDM, NumCol, datasetIndex,
				     mtzout, col);
    } else {
      // write out unsummed parts, for all datasets
      nobsbatch =  writeParts(hkl_list, NumCol, mtzout, col);
    }
    for (size_t jbat=0;jbat<nobsbatch.size();jbat++)  {
      nref += nobsbatch[jbat];
    }
    // Batches: construct linked list in mtzout object
    MTZBAT* batch;
    MTZBAT* previous_batch;
    int nbat = 0;

    for (int jbat=0;jbat<hkl_list.num_batches();jbat++)  {
      // Only output accepted batches
      if (hkl_list.batch(jbat).Accepted() && nobsbatch[jbat] > 0) {
	batch = MtzMallocBatch(); // make space for batch data
	
	if (nbat == 0) {
	  mtzout->batch = batch; // pointer to first batch
	} else {
	  // Link previous batch to this one
	  previous_batch->next = batch;
	}
	nbat++;
	previous_batch = batch;

	*batch = hkl_list.batch(jbat).batchdata(); // copy data
	// reset time limits if no time data
	if (!col_sel.is_time) {
	  batch->time1 = 0.0;
	  batch->time2 = 0.0;
	}
	// Fix up cell, constrain to symmetry
	scala::Scell Bcell = scala::Scell(batch->cell);
	// Is this a valid cell?
	bool valid = true;
	for (int i=0;i<6;i++) {
	  if (Bcell[i] <= 0.001) {valid = false;}
	}
	if (valid) {
	  scala::Scell Bcell = CG.constrain(scala::Scell(batch->cell));
	  for (int i=0;i<6;i++) {batch->cell[i] = Bcell[i];}
	}
	// Update NBsetid if required, ie index in dataset list
	// look it up in new mtzout structure

	int nbsetid = hkl_list.batch(jbat).DatasetID();
	PxdName pxdname =  hkl_list.dataset(hkl_list.batch(jbat).datasetindex()).pxdname(nbsetid);
	std::string path = "/"+pxdname.xname()+"/"+pxdname.dname();
	batch->nbsetid = MtzSetLookup(mtzout, path.c_str())->setid;  // setid
	batch->next = NULL;  // for last one
      }
    }

    ccp4_lhprt(mtzout,1);
    if (!MtzPut(mtzout, " ")) {
      Message::message(Message_fatal("Can't write file "+filename_out));
    }
    MtzFree(mtzout);

    return nref;
  }
  //--------------------------------------------------------------
  int WriteUnmerged::writeUnmergedSCA(const scala::hkl_unmerge_list& hkl_list,
				      const SDmodel& SDM,
				      const int& datasetIndex,
				      const std::string& filename_out,
				      const float& maxintensity)
  // Write unmerged scalepack file from hkl_unmerge_list object
  // returns number written
  // Skip multiples
  //
  //  datasetIndex     dataset index to output, -1 all data
  //  filename_out     output filename
  //  Imax             maximum intensity
  {
    nref = 0;
    nmultiple = 0;
    hkl_list.rewind();
    data_flags  col_sel = hkl_list.DataFlags();
    reflection this_refl;
    observation this_obs;
    int index;

    FILE* scafile = OpenFile(filename_out, true);
    if (scafile == NULL) {
      Message::message(Message_fatal("Can't write file "+filename_out));
    }

    SpaceGroup SG = hkl_list.symmetry().GetSpaceGroup();
    int nsym = SG.num_symops();

    fprintf(scafile, "%5d %s\n", nsym,
	    StringUtil::Strip(SG.Symbol_hm()).c_str());
    for (int k=0;k<nsym;++k) {
      clipper::Symop symop = SG.Symop(k);
      for (int i=0;i<3;++i) for(int j=0;j<3;++j) {
	fprintf(scafile, "%3d", Nint(symop.rot()(i,j)));
      }
      fprintf(scafile, "\n");
      for (int i=0;i<3;++i) {
	fprintf(scafile, "%3d", Nint(symop.trn()[i]));
      }
      fprintf(scafile, "\n");
    }

    float scale = 999900.0;  // keep scaled intensity in format %8.1f
    if (maxintensity > scale) {
      scale = scale/maxintensity;
    } else {
      scale = 1.0;
    }

    while (hkl_list.next_reflection(this_refl) >= 0)  { // loop reflections
      SDM.CorrectReflection(this_refl);
      scala::Hkl hkl = this_refl.hkl();
      bool Centric = hkl_list.symmetry().is_centric(hkl);

      // loop observations
      while ((index = this_refl.next_observation(this_obs)) >= 0) {
	if (this_obs.IsSingleton()) {
	  if (datasetIndex < 0 || this_obs.datasetIndex() == datasetIndex) {
	    scala::Hkl hkl_orig = this_obs.hkl_original();
	    int isym = this_obs.Isym();
	    int iasym = ((isym-1)/2+1);
	    int batch = this_obs.Batch();
	    int icn = 0; // centric
	    int ispndle=0;  // dummy here
	    if (!Centric) {
	      if (isym%2 == 0) icn = 2;  // I-
	      else icn = 1;             // I+
	    }
	    fprintf(scafile, "%4d%4d%4d%4d%4d%4d%6d%2d%2d%3d%8.1f%8.1f\n",
		    hkl_orig.h(), hkl_orig.k(), hkl_orig.l(),
		    hkl.h(), hkl.k(), hkl.l(),
		    batch, icn, ispndle, iasym,
		    scale*this_obs.kI(), scale*this_obs.ksigI());
	    nref++;
	  }
	} else {
	  nmultiple++; // omitted
	}
      } // end loop observations
    }  // end loop reflections
    return nref;
  }
  //--------------------------------------------------------------
  std::vector<int>  WriteUnmerged::writeParts(const scala::hkl_unmerge_list& hkl_list,
					      const int& NumCol,
					      MTZ* mtzout,
					      MTZCOL* col[])
  // write out all parts, all datasets
  {
    // Count observation parts in each batch
    std::vector<int> nobsbatch(hkl_list.num_batches(),0);
    // Write all data
    float data[MAXNCOLUMNS];
    data_flags  col_sel = hkl_list.DataFlags();

    for (int i=0;i<hkl_list.num_parts();i++)  {
      scala::observation_part part = hkl_list.find_part(i);
      scala::Hkl hkl = part.hkl();
      data[0] = hkl.h();
      data[1] = hkl.k();
      data[2] = hkl.l();
      // Packed M/ISYM, M = 1 for partial
      int isym = part.isym();
      int M_Isym = isym;
      if (part.Npart() != 1) M_Isym = 256 + isym;
      data[3] = M_Isym;
     
      data[4] = part.batch();
      data[5] = part.Ic();
      data[6] = part.sigIc();
      
      // Optional columns
      int ic = 7;
      if (col_sel.is_Ipr) data[ic++] = part.Ipr();
      if (col_sel.is_sigIpr) data[ic++] = part.sigIpr();
      if (col_sel.is_fractioncalc) data[ic++] = part.fraction_calc();
      if (col_sel.is_Xdet) data[ic++] = part.Xdet();
      if (col_sel.is_Ydet) data[ic++] = part.Ydet();
      if (col_sel.is_Rot) data[ic++] = part.phi();
      if (col_sel.is_Width) data[ic++] = part.width();
      if (col_sel.is_LP) data[ic++] = part.LP();
      if (col_sel.is_Mpart) {
	int Mpart = 0;
	if (part.Npart() > 1) Mpart = 100*part.Npart() + part.Ipart();
	data[ic++] = Mpart;
      }
      if (col_sel.is_ObsFlag) data[ic++] = part.ObsFlag().Flags();
      if (col_sel.is_BgPkRatio) data[ic++] = part.ObsFlag().BgPk();
      if (col_sel.is_time) data[ic++] = part.time();
      // Optional scale: dummy for now
      if (col_sel.is_scale) data[ic++] = 1.0;
      if (col_sel.is_sigscale) data[ic++] = 0.0;
      ASSERT (ic == NumCol);

      // count observation parts by batch serial
      nobsbatch.at(hkl_list.batch_serial(Nint(data[4])))++;
      
      ccp4_lwrefl(mtzout, data, col, NumCol, i+1);
    }
    return nobsbatch;
  }
  //--------------------------------------------------------------
  std::vector<int>  WriteUnmerged::writeObservations(const scala::hkl_unmerge_list& hkl_list,
						     const SDmodel& SDM,
						     const int& NumCol,
						     const int& datasetIndex,
						     MTZ* mtzout,
						     MTZCOL* col[])
  // write out summed observations, for selected dataset(s)
  // omitting rejections
  {
    // Count observations in each batch
    std::vector<int> nobsbatch(hkl_list.num_batches(),0);
    // Write all data
    float data[MAXNCOLUMNS];
    data_flags  col_sel = hkl_list.DataFlags();
    int maxhkloverlap = hkl_list.MaxHKLoverlap();

    reflection this_refl;
    observation this_obs;
    int index;
    int ic;
    int i = 0;

    while (hkl_list.next_reflection(this_refl) >= 0)  { // loop reflections
      SDM.CorrectReflection(this_refl);
      scala::Hkl hkl = this_refl.hkl();
      data[0] = hkl.h();
      data[1] = hkl.k();
      data[2] = hkl.l();

      // loop observations
      while ((index = this_refl.next_observation(this_obs)) >= 0) {
	if (datasetIndex < 0 || this_obs.datasetIndex() == datasetIndex) {
	  // Packed M/ISYM
	  // always M = 0 for "full" since partials have been summed (but see NPART)
	  int isym = this_obs.Isym();
	  ic = 3;
	  data[ic++] = isym;
	  data[ic++] = this_obs.Batch();
	  data[ic++] = this_obs.kI();
	  data[ic++] = this_obs.ksigI();
	  // Applied scale
	  float g = this_obs.Gscale();
	  if (g != 0.0) g = 1.0f/g;
	  data[ic++] = g;
	  data[ic++] = 0.0;
	  data[ic] = this_obs.num_parts(); // NPART = number of parts
	  // negate for scaled partial
	  if (this_obs.PartFlag() == SCALE) {data[ic] = -data[ic];}
	  ic++;

	  // Optional columns
	  if (col_sel.is_fractioncalc) data[ic++] = this_obs.TotalFraction();
	  std::pair<float,float> xydet = this_obs.XYdet();
	  if (col_sel.is_Xdet) data[ic++] = xydet.first;
	  if (col_sel.is_Ydet) data[ic++] = xydet.second;
	  if (col_sel.is_Rot) data[ic++]  = this_obs.phi();
	  if (col_sel.is_Width) data[ic++] = this_obs.width();
	  if (col_sel.is_LP) data[ic++] = this_obs.LP();
	  if (col_sel.is_time) data[ic++] = this_obs.time();

	  if (col_sel.is_latnum) {
	    data[ic++] = this_obs.MainLatticeNumber();
	    for (int j=0;j<maxhkloverlap*5;++j) {
	      data[ic+j] = 0.0; // clear multilattice columns
	    }
	    std::vector<LatticeIndexInfo> lathkl = this_obs.lathkl();
	    if (lathkl.size() > 0) {
	      nmultiple++;
	      for (size_t l=0; l<lathkl.size(); l++) { 
		data[ic++] = lathkl[l].latnum;
		data[ic++] = lathkl[l].hkl[0];
		data[ic++] = lathkl[l].hkl[1];
		data[ic++] = lathkl[l].hkl[2];
		data[ic++] = lathkl[l].gscale;
	      }
	    }
	  }
	  
	  // count observation parts by batch serial
	  nobsbatch.at(hkl_list.batch_serial(Nint(data[4])))++;
	  
	  ccp4_lwrefl(mtzout, data, col, NumCol, i+1);
	  i++;
	}  // end dataset selection
      } // end loop observations
    }  // end loop reflections
    return nobsbatch;
  }
  //--------------------------------------------------------------
}
