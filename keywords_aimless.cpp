// keywords_aimless.cpp

#include "keywords_aimless.hh"
#include "jiffy.hh"
#include "scala_util.hh"
#include "hkl_datatypes.hh"
#include "scaletypes.hh"

// Clipper
#include <clipper/clipper.h>
#include "clipper/clipper-ccp4.h"
using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;

namespace phaser_io {

//--------------------------------------------------------------
ANOMALOUS::ANOMALOUS() : CCP4base(), InputBase()
{
  Add_Key("ANOM");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);

  anomalous = false;
}
//--------------------------------------------------------------
Token_value ANOMALOUS::parse(std::istringstream& input_stream)
//  ANOMALOUS [ON] [OFF]
{
  anomalous = true;
  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (keyIs("ON"))
	{anomalous = true;}
      else if (keyIs("OFF"))
	{anomalous = false;}
      else
	{throw SyntaxError
	    (keywords, "key not ON or OFF");}
    }
  }
  return ENDLINE;
}
//--------------------------------------------------------------
SCALES::SCALES() : CCP4base(), InputBase()
{
  Add_Key("SCAL");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);

  nspecs = 1;
  specs.clear();
  // Set up default SCALES parameters
  scala::ScaleSpecification spec_default;
  spec_default.run = -1;        // all runs
  spec_default.batch = false;   // smooth
  spec_default.nscales = -1;
  spec_default.spacing = 5.0;   // rotation spacing 5
  spec_default.nbfac = -1;
  spec_default.bspacing = 20.;  // brotation spacing 20
  ////  spec_default.sec_abs = scala::NONE;  // for now
  spec_default.sec_abs = scala::SecondaryScale::SECONDARY;
  spec_default.lmax = 4;        // lmax for spherical harmonics
  spec_default.pole = -1;       // unspecified pole
  spec_default.ntilex = -1;
  spec_default.detectorscaletype = scala::DetectorScale::NONE;
  specs.push_back(spec_default);
}
//--------------------------------------------------------------
Token_value SCALES::parse(std::istringstream& input_stream)
// SCALES [RUN <irun>]
// [BATCH || ROTATION [<nscales> || SPACING <spacing>]
// BFACTOR ON || OFF BROTATION [<nscales> || SPACING <spacing>]
// [SECONDARY  [<Lmax> [<LmaxOdd>]]]
// [ABSORPTION [<Lmax> [<LmaxOdd>]] [POLE [h|k|l]]]
// [CONSTANT]
// [TILE [<Ntilex> [<Ntiley>]] [CCD | FLAT | PIXEL]]
{
  // Read one SCALES specification
  int irun = -1;
  scala::ScaleSpecification spec;
  //**  if (nspecs > 0) spec = specs[0];
  int expectingNumber = 0; // = 0 not expecting number, +1 expecting number
			   // = -1 maybe expecting number
  int batch_spec = 0;
  int bfac_spec = 0;
  int secabs = 0;         // SECONDARY or ABSORPTION given
			  // = +1 looking for Lmax, = +2 Lmax read, = +3 LmaxOdd read
  int tile = -1;

  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (expectingNumber > 0) {throw SyntaxError
	  (keywords, "SCALES: expecting number not "+string_value);}
      secabs = 0;
      if (keyIs("RUN")) {
	// RUN <Irun>, next token must be number
	irun = 0;
	expectingNumber = +1;
      } else if (keyIs("CONSTANT")) {
	spec.SetConstant(irun);
      } else if (keyIs("BATCH")) {
	if (batch_spec != 0) {throw SyntaxError
	    (keywords, "SCALES: can't have BATCH & ROTATION)");}
	spec.batch = true;
	batch_spec = -1;
	expectingNumber = 0;
      } else if (keyIs("ROTATION")) {
	if (batch_spec < 0) {throw SyntaxError
	    (keywords, "SCALES: can't have BATCH & ROTATION)");}
	spec.batch = false;
	batch_spec = +1;
	expectingNumber = -1;
      } else if (keyIs("SPACING") &&  bfac_spec == 0) {
	if (batch_spec < 0) {throw SyntaxError
	    (keywords, "SCALES: can't have BATCH & ROTATION)");}
	if (batch_spec == +2) {throw SyntaxError
	    (keywords, "SCALES: can't have ROTATION number & SPACING)");}
	spec.batch = false;
	batch_spec = +3;
	expectingNumber = +1;
      } else if (keyIs("BFACTOR")) {
	// BFACTOR default to ON
	expectingNumber = 0;
      } else if (keyIs("ON")) {
      } else if (keyIs("OFF")) {
	spec.nbfac = 0;  // switch off B-factors
      } else if (keyIs("BROTATION")) {
	if (batch_spec < 0) {throw SyntaxError
	    (keywords, "SCALES: can't have BATCH & BROTATION)");}
	spec.batch = false;
	bfac_spec = +1;
	expectingNumber = -1;
      } else if (keyIs("SPACING")  &&  bfac_spec > 0) {
	if (bfac_spec == +2) {throw SyntaxError
	    (keywords, "SCALES: can't have BROTATION number & SPACING)");}
	spec.batch = false;
	bfac_spec = +3;
	expectingNumber = +1;
      } else if (keyIs("SECONDARY")) {
	spec.sec_abs = scala::SecondaryScale::SECONDARY;
	expectingNumber = -1;
	secabs = +1;
      } else if (keyIs("ABSORPTION")) {
	spec.sec_abs = scala::SecondaryScale::ABSORPTION;
	expectingNumber = -1;
	secabs = +2;
      } else if (keyIs("POLE")) {
	if (secabs == +1) {throw SyntaxError
	    (keywords, "SCALES: can't have SECONDARY & POLE)");}
	secabs = -1;
	expectingNumber = 0;
      } else if (secabs < 0) {
	// POLE read, look for h|k|l
	if (keyIs("H")) {
	  spec.pole = +1;
	} else if (keyIs("K")) {
	  spec.pole = +2;
	} else if (keyIs("L")) {
	  spec.pole = +3;
	} else {
	  throw SyntaxError
	    (keywords, "SCALES: POLE must be h, k or l)");
	}
      } else if (keyIs("TILE")) {
	tile = 0;  // expecting Ntilex next
	spec.detectorscaletype = scala::DetectorScale::AUTOMATIC;
	expectingNumber = -1;
      } else if (tile >= 0) {
      	if (keyIs("FLAT")) {
	  spec.detectorscaletype = scala::DetectorScale::FLAT;
	} else if (keyIs("CCD")) {
	  spec.detectorscaletype = scala::DetectorScale::CCD;
	} else if (keyIs("PIXEL")) {
	  spec.detectorscaletype = scala::DetectorScale::PIXEL;
	}
      }
    } else if (tokenIs(1,NUMBER)) {
      if (expectingNumber == 0) {throw SyntaxError
	  (keywords, "SCALES: not expecting number"+
	   clipper::String(number_value));}
      expectingNumber = 0;
      if (irun == 0) {
	// Run number
	spec.run = Nint(number_value);
	irun = spec.run;
      } else if (batch_spec == +1) {
	// ROTATION <Nscales>
	spec.nscales = Nint(number_value);
	batch_spec = +2;
      } else if (batch_spec == +3) {
	// ROTATION SPACING
	spec.nscales = -1;
	spec.spacing = number_value;
	batch_spec = +4;
      } else if (bfac_spec == +1) {
	// BROTATION <Nbfacs>
	spec.nbfac = Nint(number_value);
	bfac_spec = +2;
      } else if (bfac_spec == +3) {
	// BROTATION SPACING
	spec.nbfac = -1;
	spec.bspacing = number_value;
	bfac_spec = +4;
      } else if (secabs == 1) {
	// Read Lmax
	spec.lmax = Nint(number_value);
	//  force even
	spec.lmax = (spec.lmax/2)*2;
	secabs = +2;
      } else if (secabs == 2) {
	// Read LmaxOdd
	spec.lmaxodd = Nint(number_value);
	//  force odd & < lmax
	spec.lmaxodd = ((Min(spec.lmaxodd, spec.lmax)+1)/2)*2-1;
	secabs = +3;
      } else if (tile == 0) {
	spec.ntilex = Nint(number_value);
	spec.ntiley = spec.ntilex;
	tile = +1;
	expectingNumber = -1;
      } else if (tile == +1) {
	spec.ntiley = Nint(number_value);
	tile = +2;
      }
    } else {
      throw SyntaxError
	(keywords, "SCALES: invalid syntax");
    }
  }
  if (secabs == 2) {
    // set lmaxodd if not defined & lmax is
    spec.lmaxodd = ((spec.lmax+1)/2)*2-1;
  }
  if (secabs != 0 && spec.lmax == 0) {
    spec.sec_abs = scala::SecondaryScale::NONE;
  }

  if (irun >= 0) {
    // Run specified, increment specification count
    specs.push_back(spec);
  } else {
    // If no run specified, overwrite 1st slot, ie the default
    specs[0] = spec;
  }
  nspecs = specs.size();  
  return ENDLINE;
}
//--------------------------------------------------------------
RUNSET::RUNSET() : CCP4base(), InputBase()
{
  Add_Key("RUN");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value RUNSET::parse(std::istringstream& input_stream)
// Syntax:
// NOT DONE RUN <runnumber> DATASET <datasetname> | <crystalname>/<datasetname>
//     assign a dataset to a run
//   RUN <runnumber> [FILE|SERIES  <Jfile>]  BATCH <b1> TO <b2>
//     assign batch range to a run: if FILE key present, b1 etc
//     refer to original (file) batch numbers (SERIES is the equivalent for a
//     series of files defined with a wild-card)
//   RUN <runnumber> ALL
{
  int Select = -1; // 0 nothing; +1 batch selection; -1 dataset selection
  int fileSeries = 0; // >0 if FILE|SERIES keywords given
  bool brange = false;
  bool all = false;
  std::vector<int> bnum;
  int runnum = -1;

  while (get_token(input_stream) != ENDLINE)  {
    if (Select < 0) {
      if (tokenIs(1,NUMBER)) {
	runnum = Nint(number_value);
      } else {
	  throw SyntaxError(keywords,"RUN must be followed by a run number");
      }
      Select = 0;
    } else if (tokenIs(1,NAME) && (keyIs("FILE") || keyIs("SERIES"))) {
      // Keyword FILE or SERIES (synonymous)
      fileSeries = Nint(get1num(input_stream));
      if (fileSeries < 0) {
	throw SyntaxError(keywords,"FILE | SERIES value must be > 0");
      }
    } else if (Select == 0) {
      if (tokenIs(1,NAME)) {
	if (keyIs("BATCH")) {
	  Select = +1;
	  continue;
	} else if (keyIs("ALL")) {
	  all = true;
	  Select = +1;
	  continue;
	} else if (keyIs("DATASET")) {
	  Select = -1;
	  continue;
	} else {
	  throw SyntaxError(keywords,"unrecognised keyword");
	}
      }
    } else if (Select == +1) {
      // Batch selection
      if (tokenIs(1,NAME)) {
	if (keyIs("TO")) {
	  // Range
	  if (brange) {
	    throw SyntaxError
	      (keywords,"only one batch range allowed per RUN command");
	  }
	  brange = true;
	} else {
	  throw SyntaxError(keywords,"unrecognised keyword");
	}
      } else {
	// gather numbers
	bnum.push_back(Nint(number_value));
      }
    } else if (Select == -1) {
      throw SyntaxError(keywords,"RUN DATASET option not yet implemented");
    }
  }
  // Line finished, store results
  if (Select == 0) {
    throw SyntaxError(keywords,"subkeyword BATCH must be given");
    //    throw SyntaxError(keywords,"subkeyword BATCH or DATASET must be given");
  }
  if (all) {
    batchranges.AddRange(0, 999999, fileSeries, runnum);
  } else if (brange) {
    if (bnum.size() != 2) {
      throw SyntaxError(keywords,"batch range must be given as 'n1 TO n2'");
    }
    batchranges.AddRange(bnum[0], bnum[1], fileSeries, runnum);
  } else {
    // List, fail must have range
    throw SyntaxError(keywords,"batch range must be given as 'n1 TO n2'");
  }  
  return skip_line(input_stream);
}

//--------------------------------------------------------------
RESO::RESO() : CCP4base(), InputBase()
{
  Add_Key("RESO");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
  Valid = false;
}
//--------------------------------------------------------------
Token_value RESO::parse(std::istringstream& input_stream)
{
  // Syntax:
  //  RESOlution [RUN <Irun>]  <high> [<low>]  either order
  //  RESOlution [RUN <Irun>] HIGH <high>
  //  RESOlution [RUN <Irun>] LOW  <low>
  //
  // Resolution by run may override earlier run commands 
  double high = input_resolution_range.min();
  double  low = input_resolution_range.max();
  int found = 0;
  int irun = -1;

  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (keyIs("HIGH")) {
	high = get1num(input_stream);
	found = +1;
      } else if (keyIs("LOW")) {
	low = get1num(input_stream);
	found = 0;
      } else if (keyIs("RUN")) {
	found = -1; // flag for run number
      }
    } else {
      if (found < 0) {
	irun = Nint(number_value);
	found = 0;
      } else if (found == 0) {
	high = number_value;
	found = +1;
      } else {
	low = high;
	high = number_value;
	found = 0;
      }
    }
  }
  if (irun < 0) {
    setRESO(ResoRange(low,high));
  } else {
    run_resolution_ranges.push_back(std::pair<int,ResoRange>(irun, ResoRange(low,high)));
  }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
void RESO::setRESO(const ResoRange& reso_range)
{
  input_resolution_range = reso_range;
  Valid = true;
}
//--------------------------------------------------------------
ResoRange RESO::getRESO() const {return input_resolution_range;}
//--------------------------------------------------------------
//--------------------------------------------------------------
PARTIALS::PARTIALS()
{
  Add_Key("PARTIALS");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
  // Default values
  minfraclim = 0.95;
  maxfraclim = 1.05;
  minscalefrac = 0.0;
  check = true;
  maxgap = 0;
}
//--------------------------------------------------------------
Token_value PARTIALS::parse(std::istringstream& input_stream)
// Syntax: PARTials [TEST <lower_limit> <upper_limit>] [CORRECT <minimum_fraction>]
//  [[NO]CHECK] [[NO]GAP <gap limit>]
{
  bool limits = false;
  int nlim = 0;
  bool scalelim = false;
  int nsclim = 0;
  bool gapval = false;

  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      limits = false;
      scalelim = false;
      gapval = false;
      if (keyIs("TEST")) {
	limits = true;
      } else if (keyIs("CORRECT")) {
	scalelim = true;
      } else if (keyIs("CHECK")) {
	check = true;
      } else if (keyIs("NOCHECK")) {
	check = false;
      } else if (keyIs("GAP")) {
	gapval = true;
	maxgap = 1;
      } else if (keyIs("NOGAP")) {
	maxgap = 0;
      } else
	{throw SyntaxError
	    (keywords, "unrecognised keyword");}
    } else if (tokenIs(1,NUMBER)) {
      if (limits) {
	if (nlim == 0)
	  {minfraclim = number_value;}
	else if (nlim == 1)
	  {maxfraclim = number_value;}
	else
	  {throw SyntaxError
	      (keywords, "ERROR in TEST two numbers must be given");}
	nlim++;
      } else if (scalelim) {
	if (nsclim == 0)
	  {minscalefrac = number_value;}
	else
	  {throw SyntaxError
	      (keywords,
	       "ERROR in CORRECT one numbers must be given");}
	nsclim++;
      } else if (gapval) {
	maxgap = Nint(number_value);
	if (maxgap < 0 || maxgap > 5) {
	  throw SyntaxError(keywords,"unreasonable MaxGap ");
	}
      } else
	{throw SyntaxError
	    (keywords, "unexpected number");}
    }
  }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
EXCLUDE::EXCLUDE() : CCP4base(), InputBase()
{
  Add_Key("EXCLUDE");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value EXCLUDE::parse(std::istringstream& input_stream)
// Syntax:
//   EXCLUDE DATASET <datasetname> | <crystalname>/<datasetname>
//     exclude dataset
//   EXCLUDE BATCH [FILE|SERIES  <Jfile>]  <b1> <b2> <b3> ... | <b1> TO <b2>
//     exclude batch list or range: if FILE key present, b1 etc
//     refer to original (file) batch numbers (SERIES is the equivalent for a
//     series of files defined with a wild-card)
{
  int Select = 0; // 0 nothing; +1 batch selection; -1 dataset selection
  int fileSeries = 0; // >0 if FILE|SERIES keywords given
  bool brange = false;
  std::vector<int> bnum;

  while (get_token(input_stream) != ENDLINE)  {
    if (tokenIs(1,NAME) && (keyIs("FILE") || keyIs("SERIES"))) {
	// Keyword FILE or SERIES (synonymous)
	fileSeries = Nint(get1num(input_stream));
	if (fileSeries < 0) {
	  throw SyntaxError(keywords,"FILE | SERIES value must be > 0");
	}
    } else if (Select == 0) {
      if (tokenIs(1,NAME)) {
	if (keyIs("BATCH")) {
	  Select = +1;
	  continue;
	} else if (keyIs("DATASET")) {
	  Select = -1;
	  continue;
	} else {
	  throw SyntaxError(keywords,"unrecognised keyword");
	}
      }
    } else if (Select == +1) {
      // Batch selection
      if (tokenIs(1,NAME)) {
	if (keyIs("TO")) {
	  // Range
	  if (brange) {
	    throw SyntaxError
	      (keywords,"only one batch range allowed per EXCLUDE command");
	  }
	  brange = true;
	} else {
	  throw SyntaxError(keywords,"unrecognised keyword");
	}
      } else {
	// gather numbers
	bnum.push_back(Nint(number_value));
      }
    } else if (Select == -1) {
      throw SyntaxError(keywords,"EXCLUDE DATASET option not yet implemented");
    }
  }
  // Line finished, store results
  if (Select == 0) {
    throw SyntaxError(keywords,"subkeyword BATCH must be given");
    //    throw SyntaxError(keywords,"subkeyword BATCH or DATASET must be given");
  }
  if (brange) {
    if (bnum.size() != 2) {
      throw SyntaxError(keywords,"batch range must be given as 'n1 TO n2'");
    }
    batchexclude.AddRange(bnum[0], bnum[1], fileSeries);
  } else {
    // List
    for (size_t i=0;i<bnum.size();i++) {
      batchexclude.AddBatch(bnum[i], fileSeries);
    }
  }  
  return skip_line(input_stream);
}
//--------------------------------------------------------------
BINS::BINS()  : nrbins(-1), nibins(10)
{
  Add_Key("BINS");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value BINS::parse(std::istringstream& input_stream)
{
  // Syntax:
  // BINS RESOLUTION <NresoBins> INTENSITY <NintensityBins>

  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (keyIs("RESOLUTION")) nrbins = Nint(get1num(input_stream));
      else if (keyIs("INTENSITY")) nibins = Nint(get1num(input_stream));
    }
  }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
REJECT::REJECT()
{
  Add_Key("REJECT");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value REJECT::parse(std::istringstream& input_stream)
{
  //  Outlier rejection criteria
  //
  //  [SCALE|MERGE]  use these values for scaling|merging steps
  //            if not specified, use for both
  //  [COMBINE] compare observations across all datasets [default]
  //  [SEPARATE]  outlier checks only within datasets
  //   sdrej    sd multiplier for maximum deviation from scale-weighted mean I
  //  [sdrej2]  special value for reflections measured twice
  //  [REJECT|KEEP|SMALLER|LARGER]] these flags control what to do in the 
  //            merging stage if there are two observations which disagree
  //            In scaling, both are always rejected.
  //            REJECT      reject both
  //            KEEP        keep both [default]
  //            LARGER      reject larger
  //            SMALLER     reject smaller
  //  [ALL sdreja [sdrej2a]]
  //            if ALL is present and the ANOMALOUS option is on, then
  //            test all data (ie all I+ & I-) values against these limits
  //            after testing within the I+ & I- set
  //            NOALL (== ALL 0) switches off this test
  //            Only applies to merging step (scaling step checks all anyway)
  //  EMAX <Emax> maximum normalised F accepted
  //            
  float sdrej = 0.0;
  float sdrej2 = 0.0;
  float sdreja = 0.0;
  float sdrej2a = 0.0;
  float emax = -1.0;
  bool  emaxgiven = false;
  bool combine = true;

  int merge = 0;  // expecting values for MERGE && SCALE, = +1 for MERGE, = -1 for SCALE
  // enum Reject2Policy {REJECT, KEEP, REJECTLARGER, REJECTSMALLER};
  scala::RejectFlags::Reject2Policy rej2policy = scala::RejectFlags::KEEP;
  bool anom = false;
  bool first = true;  // first of pair

  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (keyIs("SCALE")) {
	if (merge == 0) merge = -1;
	else merge = 0;   // if both MERGE & SCALE are given
      } else if (keyIs("MERGE")) {
	if (merge == 0) merge = +1;
	else merge = 0;
      } else if (keyIs("COMBINE")) {
	combine = true;
      } else if (keyIs("SEPARATE")) {
	combine = true;
      } else if (keyIs("REJECT")) {
	rej2policy = scala::RejectFlags::REJECT;
      } else if (keyIs("KEEP")) {
	rej2policy = scala::RejectFlags::KEEP;
      } else if (keyIs("REJECTLARGER")) {
	rej2policy = scala::RejectFlags::REJECTLARGER;
      } else if (keyIs("REJECTSMALLER")) {
	rej2policy = scala::RejectFlags::REJECTSMALLER;
      } else if (keyIs("ALL")) {
	anom = true;
	first = true;
      } else if (keyIs("EMAX")) {
	emaxgiven = true;
      } else {
	throw SyntaxError
	  (keywords, "REJECT: unrecognised keyword");
      }
      } else if (tokenIs(1,NUMBER)) {
      if (anom) {
	if (first) {
	  sdreja = number_value;
	  sdrej2a = sdreja;
	  first = false;
	} else {
	  sdrej2a = number_value;
	}
      } else if (emaxgiven) {
	emax = number_value;
	emaxgiven = false;
      } else {
	if (first) {
	  sdrej = number_value;
	  sdrej2 = sdrej;
	  first = false;
	} else {
	  sdrej2 = number_value;
	}
      }
    } else {
      throw SyntaxError
	(keywords, "REJECT: invalid syntax");
    }
  }
  
  // Store them
  outliercontrolsscale.Combine() = combine;
  outliercontrolsmerge.Combine() = combine;
  //  SCALE or both
  if (merge <= 0) {
    outliercontrolsscale.Reject(scala::ALL) = scala::RejectFlags(sdrej, sdrej2, rej2policy);
    outliercontrolsscale.Reject(scala::BOTH) = scala::RejectFlags(sdreja, sdrej2a, rej2policy);
    if (emax > 0.0) outliercontrolsscale.SetEmax(emax);
  }
  //  MERGE or both
  if (merge >= 0) {
    outliercontrolsmerge.Reject(scala::ALL) = scala::RejectFlags(sdrej, sdrej2, rej2policy);
    outliercontrolsmerge.Reject(scala::BOTH) = scala::RejectFlags(sdreja, sdrej2a, rej2policy);
    if (emax > 0.0) outliercontrolsmerge.SetEmax(emax);
  }
  
  return skip_line(input_stream);
}
//--------------------------------------------------------------
void REJECT::analyse()
{
  // Check that rejects between anomalous related reflections is not
  // more stringent that within
  if (outliercontrolsscale.Reject(scala::BOTH).sdrej != 0.0) {
    if (outliercontrolsscale.Reject(scala::BOTH).sdrej < outliercontrolsscale.Reject(scala::ALL).sdrej) {
      Message::message(Message_warn
		       (std::string("\nWARNING: on REJECT command\n")+
			" Scaling outlier rejection limit for all data is "+
			"tighter than that within I+|I- sets\n"+
			" ==== This is not sensible"));
    }
  }
  if (outliercontrolsmerge.Reject(scala::BOTH).sdrej != 0.0) {
    if (outliercontrolsmerge.Reject(scala::BOTH).sdrej < outliercontrolsmerge.Reject(scala::ALL).sdrej) {
      Message::message(Message_warn
		       (std::string("\nWARNING: on REJECT command\n")+
			" Merging outlier rejection limit for all data is "+
			"tighter than that within I+|I- sets\n"+
			" ==== This is not sensible"));
    }
  }
}
//--------------------------------------------------------------
TIE::TIE()
{
  Add_Key("TIE");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
  // Defaults
  tiesd_surface = 0.001; 
  tiesd_rotation = -1.0;
  tiesd_bfactor = -1.0; 
  tiesd_zerob = -1.0; 
  // SDs for CCD tile, parameters r, w, A, x0|y0
  // SDs for r, w, xy0 are relative to maximum radius
  float sd_tile[] = {0.1, 0.1, 0.05, 0.05};
  tiesd_tile.assign(sd_tile, sd_tile+4);
}
//--------------------------------------------------------------
Token_value TIE::parse(std::istringstream& input_stream)
{
  // TIE <parameter> <sd> [<sd2> etc]
  //  <parameter> = SURFACE    for SECONDARY or ABSORPTION
  //              = ROTATION   for primary scale parameters (eg BATCH)
  //              = BFACTOR    for B-factors
  //              = ZEROB      for B-factors tied to B = 0
  //              = TILE       for tile correction parameters (4 sds)
  bool invalid = false;
  while (get_token(input_stream) != ENDLINE)  {
    if (tokenIs(1,NAME)) {
      if (keyIs("SURFACE") || keyIs("SECONDARY") || keyIs("ABSORPTION")) {
	tiesd_surface = get1num(input_stream);
	if (tiesd_surface < 0.0) invalid = true;
      } else if (keyIs("ROTATION") || keyIs("SCALE")) {
	tiesd_rotation = get1num(input_stream);
	if (tiesd_rotation < 0.0) invalid = true;
      } else if (keyIs("BFACTOR")) {
	tiesd_bfactor = get1num(input_stream);
	if (tiesd_bfactor < 0.0) invalid = true;
      } else if (keyIs("ZEROB")) {
	tiesd_zerob = get1num(input_stream);
	if (tiesd_zerob < 0.0) invalid = true;
      } else if (keyIs("TILE")) {
	// Expect 4 numbers
	tiesd_tile.clear();
	tiesd_tile.push_back(get1num(input_stream));
	tiesd_tile.push_back(get1num(input_stream));
	tiesd_tile.push_back(get1num(input_stream));
	tiesd_tile.push_back(get1num(input_stream));
	ASSERT (tiesd_tile.size() == 4);
	for (size_t i=0;i<tiesd_tile.size();++i) {
	  if (tiesd_tile[i] < 0.0) invalid = true;
	}
      }
    }
  }
  if (invalid) {
    throw SyntaxError
      (keywords, "TIE: sd must be > 0");
  }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
NAME::NAME() : CCP4base(), InputBase()
{
  Add_Key("NAME");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);  
}
//--------------------------------------------------------------
Token_value NAME::parse(std::istringstream& input_stream)
{
  // Syntax: NAME PROJECT <Pname> CRYSTAL <Xname> DATASET <Dname>
  std::string pname;
  std::string xname;
  std::string dname;
  while (get_token(input_stream) != ENDLINE)
    {
      if (keyIs("PROJECT")) {
	if (get_token(input_stream) != ENDLINE)
	  {pname = string_value;}
      }
      else if (keyIs("CRYSTAL")) {
	if (get_token(input_stream) != ENDLINE)
	  {xname = string_value;}
      }
      else if (keyIs("DATASET")) {
	if (get_token(input_stream) != ENDLINE)
	  {dname = string_value;}
      }
    }
  pxdname = scala::PxdName(pname, xname, dname);
  return ENDLINE;
}
//--------------------------------------------------------------
scala::PxdName NAME::pxdname;
//--------------------------------------------------------------
void NAME::setNAME(const scala::PxdName& PXDname) 
// Only reset blank fields
  {
    if (pxdname.pname() == "") {
      pxdname.pname() = PXDname.pname();
    }
    if (pxdname.xname() == "") {
      pxdname.xname() = PXDname.xname();
    }
    if (pxdname.dname() == "") {
      pxdname.dname() = PXDname.dname();
    }
  }
//--------------------------------------------------------------
REFINE::REFINE() : CCP4base(), InputBase()
{
  Add_Key("REFINE");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);  
  // default values set in refinecontrol
}
//--------------------------------------------------------------
Token_value REFINE::parse(std::istringstream& input_stream)
// Syntax: REFINE BFGS|FH CYCLE [Ncyc1>] <Ncycles> CONVERGE <convergeLimit>
//            SELECT <IovSDmin> <E2min> [<E2max>]
//    CYCLE
//         Ncyc1   number of cycles in 1st stage [default 2]
//         Ncycles number of cycles in main scaling [default 10]
//    BFGS use BFGS minimiser
//    FH   use Fox_Holmes least-squares
//    CONVERGE set convergence limit, multiplier of SD (for FH only)
// Selection criteria for scaling reflections:
//    IovSDmin   <I>/sd'(<I>) limit for 1st pass scaling
//    E2min      |E^2| limit for 2nd pass scaling
{
  int nn = -1; // counter for CYCLES
  int ns = -1; // counter for SELECT
  int nc = -1; // counter for CONVERGE
  int expectingNumber = -1; // = 0 not expecting number, +1 expecting number
			   // = -1 maybe expecting number

  while (get_token(input_stream) != ENDLINE)  {
    if (tokenIs(1,NAME)) {
      if (expectingNumber > 0) {
	throw SyntaxError
	  (keywords, "syntax error expecting a number");
      }
      if (keyIs("BFGS")) {
	refinecontrol.BFGS() = true;
	expectingNumber = 0;
      } else if (keyIs("FH")) {
	refinecontrol.BFGS() = false;
	expectingNumber = 0;
      } else if (keyIs("CYCLES")) {
	if (nn < 0) {
	  expectingNumber = +1;
	  nn = 0;
	} else {
	  expectingNumber = -1;
	  nn = 0;
	}
      } else if (keyIs("CONVERGE")) {
	  expectingNumber = +1;
	  nc = 0;
      } else if (keyIs("SELECT")) {
	if (ns < 0) {
	  expectingNumber = +1;
	  ns = 0;
	}
      }
    } else if (tokenIs(1,NUMBER)) {
      if (expectingNumber == 0) {
	throw SyntaxError
	  (keywords, "syntax error not expecting a number");
      }
      if (nn == 0) {
	refinecontrol.Ncycles() = Nint(number_value);
	nn++;
	expectingNumber = -1;
      } else if (nn == 1) {
	refinecontrol.Ncyc1() = refinecontrol.Ncycles();
	refinecontrol.Ncycles() = Nint(number_value);
	nn = -1;
	expectingNumber = 0;
      } else if (nc >= 0) {
	refinecontrol.Converge() = number_value;
	nc = -1;
	expectingNumber = 0;
      } else if (ns == 0) {
	refinecontrol.IovSDmin() = number_value;
	ns++;
      } else if (ns == 1) {
	refinecontrol.E2min() = number_value;
	ns++;
      } else if (ns == 2) {
	refinecontrol.E2max() = number_value;
	ns = -1;
	expectingNumber = 0;
      }
    }
  }
  return ENDLINE;
}
//--------------------------------------------------------------
TITLE::TITLE() : CCP4base(), InputBase()
{
  Add_Key("TITLE");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);  
}
//--------------------------------------------------------------
Token_value TITLE::parse(std::istringstream& input_stream)
// Syntax: TITLE <run title>
{
  title = getLine(input_stream);
  return ENDLINE;
}
//--------------------------------------------------------------
ONLYMERGE::ONLYMERGE() : CCP4base(), InputBase()
{
  Add_Key("ONLYMERGE");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
  onlymerge = false;
}
//--------------------------------------------------------------
Token_value ONLYMERGE::parse(std::istringstream& input_stream)
// Syntax: ONLYMERGE
{
  onlymerge = true;
  return ENDLINE;
}
//--------------------------------------------------------------
BLANK::BLANK() : CCP4base(), InputBase()
{
  Add_Key("BLANK");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);

  nullResolutionfraction = -1.0;
  nullNegativeReject = -1.0;
}
//--------------------------------------------------------------
Token_value BLANK::parse(std::istringstream& input_stream)
{
  // Read parameters for detecting blank images
  // Syntax:
  //  BLANK [RESO[FRACTION] <fraction>] [NEGATIVE <negativefraction>]
  //    <fraction>   blank image test will be done out to this fraction
  //                 of the maximum resolution, default 0.5
  //    <negativefraction> fraction of reflections with negative intensity
  //                 above which the image is considered to be blank
  //                 Default 0.3
  while (get_token(input_stream) != ENDLINE)  {
    if (tokenIs(1,NAME))      {
      if (keyIs("RESO")) nullResolutionfraction = get1num(input_stream);
      if (keyIs("NEGATIVE")) nullNegativeReject = get1num(input_stream);
    }
  }
  if (nullResolutionfraction > 0.0) {
    if (nullResolutionfraction < 0.1 || nullResolutionfraction > 1.0001) {
      throw SyntaxError("ERROR in BLANK command, Resolutionfraction must be between 0.1 and 1.0",
			"");
    }}
  if (nullNegativeReject > 0.0) {
    if (nullNegativeReject < 0.01 || nullNegativeReject > 0.501) {
      throw SyntaxError("ERROR in BLANK command, negative fraction must be between 0.01 and 0.5",
			"");
    }}
  return skip_line(input_stream);
}
//--------------------------------------------------------------
SDCORRECTION::SDCORRECTION() : CCP4base(), InputBase()
{
  Add_Key("SDCORRECTION");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);

  refine = true;
  refine_set = false;
  allsame = false;
  fixsdb = false;
  sdinput.clear();
  runnumbers.clear();
  damp = -1.0;
  tietype = -1;  // default tietype
  targets.assign(3,0.0);
  sdtargets.assign(3,0.0);
}
//--------------------------------------------------------------
Token_value SDCORRECTION::parse(std::istringstream& input_stream)
// Read parameters for SD correction
// Syntax:
//   SDCORRECTION [[NO]REFINE] [INDIVIDUAL|SAME] [FIXSDB]
//     [RUN <RunNumber>] [FULL | PARTIAL] <SdFac> [<SdB>] <SdAdd>
//     DAMP <dampfactor>
//     TIE [<parameter> <value> <sd>] | NOTIE
//  <parameter> is "SdFac" "SdB" or "SdAdd" (case insensitive)
{
  int expectingNumber = -1; // = 0 not expecting number, +1 expecting number
			   // = -1 maybe expecting number
  int irun = -1;
  int k = 0;
  bool dampset = false;
  int tieset = -1;
  double SDfac;
  double SDb;
  double SDadd;
  scala::SDcorrection sdcfull;
  scala::SDcorrection sdcpartial;
  bool found = false;

  int fullpart = 0;  // +1 full, -1 partial, 0 both

  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (expectingNumber > 0) {throw SyntaxError
	  (keywords, "SDCORRECTION: expecting number not "+string_value);}
      expectingNumber = -1;
      // Do we have some values ready to be stored?
      if (k > 0) {
	if (fullpart >= 0) {
	  sdcfull.Set(SDfac, SDb, SDadd);
	}
	if (fullpart <= 0) {
	  sdcpartial.Set(SDfac, SDb, SDadd);
	}
	k = 0;
	found = true;
      }
      if (keyIs("REFINE")) {
	refine = true;
	refine_set = true;
      } else if (keyIs("NOREFINE")) {
	refine = false;
	refine_set = true;
      } else if (keyIs("INDIVIDUAL")) {
	allsame = false;
      } else if (keyIs("SAME")) {
	allsame = true;
      } else if (keyIs("FIXSDB")) {
	fixsdb = true;
      } else if (keyIs("RUN")) {
	expectingNumber = +1;
	irun = 0;
      } else if (keyIs("FULL")) {
	expectingNumber = +1; // should be followed by 2 or 3 numbers
	fullpart = +1;
	k = 0;
      } else if (keyIs("PARTIAL")) {
	expectingNumber = +1;
	fullpart = -1;
	k = 0;
      } else if (keyIs("DAMP")) {
	expectingNumber = +1;
	dampset = true;
      } else if (keyIs("TIE")) {
	expectingNumber = -1;
	tieset = 0;
	tietype = -1;  // default tietype
      } else if (keyIs("NOTIE")) {
	expectingNumber = 0;
	tieset = -1;
	tietype = 0;  // no ties tietype
      } else if (keyIs("SDFAC")) {
	if (tieset != 0) {
	  throw SyntaxError
	    (keywords, "SDFAC not expected except after TIE");
	}
	tieset = +1;
	expectingNumber = +1;
      } else if (keyIs("SDB")) {
	if (tieset != 0) {
	  throw SyntaxError
	    (keywords, "SDB not expected except after TIE");
	}
	tieset = +2;
	expectingNumber = +1;
      } else if (keyIs("SDADD")) {
	if (tieset != 0) {
	  throw SyntaxError
	    (keywords, "SDADD not expected except after TIE");
	}
	tieset = +3;
	expectingNumber = +1;
      } else {
	throw SyntaxError
	  (keywords, "unrecognised keyword");
      }
    } else if (tokenIs(1,NUMBER)) {
      if (expectingNumber == 0) {throw SyntaxError
	  (keywords, "SDCORRECTION: not expecting number"+
	   clipper::String(number_value));}
      if (dampset) {
	damp = number_value;
	dampset = false;
      } else if (tieset > 0) {
	if (tieset < 10) {  // target value
	  targets.at(tieset-1) = number_value;
	  tieset += 10; // reset for SD
	} else { // SD
	  sdtargets.at(tieset-11) = number_value;
	  tieset = -1; // reset
	  expectingNumber = 0;
	}
	tietype = +1;
      } else if (irun == 0) {
	// Run number
	irun = Nint(number_value);
      } else {
	// Interpret number as k'th SD parameter (0,1,2)
	if (k == 0) {
	  SDfac = number_value;
	  SDadd = 0.0;
	  expectingNumber = -1;
	} else if (k == 1) {
	  SDb = 0.0;
	  SDadd = number_value; // into Sdadd in case there is no SdB
	  expectingNumber = -1;
	} else if (k == 2) {
	  SDb = SDadd; // swap into SDb
	  SDadd = number_value;
	  expectingNumber = -1;
	} else {
	  throw SyntaxError
	    (keywords, "SDCORRECTION:: more than 3 values");
	}
	k++;
      }
    }
  }
  // Do we have some values ready to be stored?
  if (k > 0) {
    if (fullpart >= 0) {
      sdcfull.Set(SDfac, SDb, SDadd);
    }
    if (fullpart <= 0) {
      sdcpartial.Set(SDfac, SDb, SDadd);
    }
    k = 0;
    found = true;
  }

  if (tieset >= 10) {
    throw SyntaxError
      (keywords, "SDCORRECTION:: TIE needs two values for target and weight");
  }
  if (found) {
    if (irun == 0) irun = -1;
    if (sdinput.size() == 0) { // none stored yet
      // Store if values read
      sdinput.push_back
	(std::pair<scala::SDcorrection,scala::SDcorrection>(sdcfull, sdcpartial));
      runnumbers.push_back(irun);
    } else {
      // we already have some, check new one against what we have
      found = false;
      for (size_t i=0;i<sdinput.size();++i) {
	if (irun == runnumbers[i]) { // found
	  if (fullpart < 0) { // new values for partial only
	    sdcfull = sdinput[i].first;
	  } else if (fullpart > 0) { // new values for full only
	    sdcpartial = sdinput[i].second;
	  }
	  // replace values
	  sdinput[i] =
	    std::pair<scala::SDcorrection,scala::SDcorrection>(sdcfull, sdcpartial);
	  found = true;
	  break;
	}
      }
      if (!found) {
	// Store new values in new slot
	sdinput.push_back
	  (std::pair<scala::SDcorrection,scala::SDcorrection>(sdcfull, sdcpartial));
	runnumbers.push_back(irun);
      }
    }
  }
  ASSERT (sdinput.size() == runnumbers.size());
  return skip_line(input_stream);
}
//--------------------------------------------------------------
//! return number of input corrections, = 0 none, = -1 overall
int SDCORRECTION::SDC_NumberInput() const 
{
  if (runnumbers.size() == 1) {
    if (runnumbers[0] < 0) return -1;
  }
  return sdinput.size();
}
//--------------------------------------------------------------
//! return ties, target & SDs, = 0 no tie, = -1 defaults, = +1 set here
int SDCORRECTION::SDCties(std::vector<double>& Targets,
			  std::vector<double>& SDtarget) const
{ 
  Targets = targets;
  SDtarget = sdtargets;
  return tietype;
}
//--------------------------------------------------------------
void SDCORRECTION::analyse()
// if explicit values have been given, default refine flag to false
// unless explicitly set
{
  if (SDC_NumberInput() > 1 && 	allsame == true) {
    throw SyntaxError
      (keywords, "SDCORRECTION:: multiple values givem with SAME flag");
  }
  if (refine_set) return;  // explicit refine flag set
  if (SDC_NumberInput() != 0) {
    refine = false;
  }
}
//--------------------------------------------------------------
INTENSITIES::INTENSITIES() : CCP4base(), InputBase()
{
  Add_Key("INTENSITIES");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);

  selecticolflag = +1; // default COMBINE
  ipowercomb = 3;
  imid = -1.0;
}
//--------------------------------------------------------------
Token_value INTENSITIES::parse(std::istringstream& input_stream)
// Read parameters for SD correction
// Syntax:
// INTENSITIES [SUMMATION | PROFILE | COMBINE [<Imid>] [POWER <Ipower>] ]
{
  int expectingNumber = 0; // = 0 not expecting number, +1 expecting number
			   // = -1 maybe expecting number
  int type = -1;  // -1 none, = 0 COMBINE, = +1 POWER
  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (expectingNumber > 0) {throw SyntaxError
	  (keywords, "INTENSITIES: expecting number not "+string_value);}
      if (keyIs("SUMMATION") || keyIs("INTEGRATED")) {
	selecticolflag = 0;
      } else if (keyIs("PROFILE")) {
	selecticolflag = -1;
      } else if (keyIs("COMBINE")) {
	selecticolflag = +1;
	expectingNumber = -1;
	type = 0;
      } else if (keyIs("POWER")) {
	expectingNumber = +1;
	type = +1;	
      } else {
	throw SyntaxError
	  (keywords, "unrecognised keyword");
      }
    } else if (tokenIs(1,NUMBER)) {
      if (expectingNumber == 0) {throw SyntaxError
	  (keywords, "INTENSITIES: not expecting number"+
	   clipper::String(number_value));}
      if (type == 0) {
	// COMBINE set Imid
	imid = number_value;
      } else if (type == +1) {
	// POWER set power
	ipowercomb = Nint(number_value);
	if (ipowercomb < 1 || ipowercomb > 5) {
	  throw SyntaxError
	    (keywords, "INTENSITIES:: unreasonable POWER value");
	}
      } else {
	  throw SyntaxError
	    (keywords, "INTENSITIES: no number expected");
      }
      type = -1;
    }
  }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
KEEP::KEEP() : CCP4base(), InputBase()
{
  Add_Key("KEEP");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value KEEP::parse(std::istringstream& input_stream)
// Read parameters for SD correction
// Syntax:
// KEEP [OVERLOADS|BGRATIO <bgratio_max>|PKRATIO <pkratio_max>|
//       GRADIENT <bg_gradient_max>|EDGE]
{
  int expectingNumber = 0; // = 0 not expecting number, +1 expecting number
			   // = -1 maybe expecting number
  int type = -1;  // -1 none, +1 BGRATIO, +2 PKRATIO, +3 GRADIENT
  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (expectingNumber > 0) {throw SyntaxError
	  (keywords, "KEEP: expecting number not "+string_value);}
      if (keyIs("OVERLOADS")) {
	observationflagcontrol.SetAcceptOverload();
	expectingNumber = 0;
      } else if (keyIs("BGRATIO")) {
	type = +1;
	expectingNumber = +1;
      } else if (keyIs("PKRATIO")) {
	type = +2;
	expectingNumber = +1;
      } else if (keyIs("GRADIENT")) {
	type = +3;
	expectingNumber = +1;
      } else if (keyIs("EDGE")) {
	observationflagcontrol.SetAcceptEdge();
	expectingNumber = 0;
      } 
    } else if (tokenIs(1,NUMBER)) {
      if (expectingNumber == 0) {throw SyntaxError
	  (keywords, "INTENSITIES: not expecting number"+
	   clipper::String(number_value));}
      if (type == +1) {
	observationflagcontrol.SetBGRlimit(number_value);
	expectingNumber = 0;
      } else if (type == +2) {
	observationflagcontrol.SetPKRlimit(number_value);
	expectingNumber = 0;
      } else if (type == +3) {
	observationflagcontrol.SetGradlimit(number_value);
	expectingNumber = 0;
      } else {
	  throw SyntaxError
	    (keywords, "KEEP: no number expected");
      }
      type = -1;
    }
  }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
OUTPUT::OUTPUT() : CCP4base(), InputBase()
{
  Add_Key("OUTPUT");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value OUTPUT::parse(std::istringstream& input_stream)
// Read parameters for SD correction
// Syntax:
//  OUTPUT [MTZ] [NO]MERGED | UNMERGED [SPLIT | TOGETHER]
//        [POLISH MERGED | UNMERGED]
{
  bool mtz = true;
  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (keyIs("MERGED") || keyIs("AVERAGE")) {
	if (mtz) {
	  outputcontrols.SetMTZoutputType(scala::OutputControls::MERGED);
	} else {
	  outputcontrols.SetSCAoutputType(scala::OutputControls::MERGED);
	}
      } else if (keyIs("NOMERGED")) {
	if (mtz) {
	  outputcontrols.SetMTZoutputType(scala::OutputControls::NONE);
	} else {
	  outputcontrols.SetSCAoutputType(scala::OutputControls::NONE);
	}
      } else if (keyIs("UNMERGED")) {
	if (mtz) {
	  outputcontrols.SetMTZoutputType(scala::OutputControls::UNMERGED);
	} else {
	  outputcontrols.SetSCAoutputType(scala::OutputControls::UNMERGED);
	}
      } else if (keyIs("NONE")) {
	if (mtz) {
	  outputcontrols.SetMTZoutputType(scala::OutputControls::NONE);
	} else {
	  outputcontrols.SetSCAoutputType(scala::OutputControls::NONE);
	}
      } else if (keyIs("MTZ")) {
	mtz = true;
      } else if (keyIs("POLISH") || keyIs("SCALEPACK")) {
	mtz = false;
	outputcontrols.SetSCAoutputType(scala::OutputControls::MERGED);
      } else if (keyIs("SPLIT")) {
	outputcontrols.Split() = true;
      } else if (keyIs("TOGETHER")) {
	outputcontrols.Split() = false;
      } else {
	  throw SyntaxError
	    (keywords, "OUTPUT: unrecognised keyword");
      }
    }
  }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
DUMP::DUMP() : CCP4base(), InputBase()
{
  Add_Key("DUMP");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
  dumpfilename = "SCALES";
}
//--------------------------------------------------------------
Token_value DUMP::parse(std::istringstream& input_stream)
//  Dump scale model from file [default SCALES]
// Syntax:
//  DUMP <filename>
{
  while (get_token(input_stream) != ENDLINE) {
    dumpfilename = string_value;
  }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
RESTORE::RESTORE() : CCP4base(), InputBase()
{
  Add_Key("RESTORE");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
  restorefilename = "SCALES";
  restore = false;
}
//--------------------------------------------------------------
Token_value RESTORE::parse(std::istringstream& input_stream)
// Read parameters for SD correction
// Syntax:
//  RESTORE <filename>
{
  restore = true;
  while (get_token(input_stream) != ENDLINE) {
    restorefilename = string_value;
  }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
ANALYSIS::ANALYSIS()  : coneangledegrees(20.0),
			minimumhalfdatasetcc(0.5),
			minimumioversigma(2.0),
			minimumbatchioversigma(1.0),
			smoothstatisticsrange(-1.0),
			detector(false)
{
  Add_Key("ANALYSIS");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value ANALYSIS::parse(std::istringstream& input_stream)
{
  // Syntax:
  //  ANALYSIS  CONE <angle>
  //    CCMINIMUM <MinimumHalfdatasetCC>
  //    ISIGMINIMUM <MinimumIoverSigma>
  //    BATCHISIGMINIMUM <MinimumBatchIoverSigma>
  //    SMOOTHSTATISTICS <SmoothStatisticsRange>
  //
  // Cone angle is the half-angle (degrees) for cones around each reciprocal axis
  // MinimumHalfdatasetCC  minimum CC for resolution warning
  // MinimumIoverSigma          minimum <<I>/sd(<I>)> for resolution warning
  // MinimumBatchIoverSigma     minimum <I/sd(I)> for resolution warning by batch, from unmerged I
  // SmoothStatisticsRange angle in degrees over which (roghly) to smooth
  //            batch statistics, <0 to default to automatic setting

  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (keyIs("CONE")) coneangledegrees = get1num(input_stream);
      else if (keyIs("CCMINIMUM"))   minimumhalfdatasetcc = get1num(input_stream);
      else if (keyIs("ISIGMINIMUM")) minimumioversigma = get1num(input_stream);
      else if (keyIs("BATCHISIGMINIMUM")) minimumbatchioversigma = get1num(input_stream);
      else if (keyIs("DETECTOR")) detector = true;
      else if (keyIs("NODETECTOR")) detector = false;
    }
  }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
} // phaser_io
