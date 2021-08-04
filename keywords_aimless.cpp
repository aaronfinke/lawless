// keywords_aimless.cpp

#include "keywords_aimless.hh"
#include "jiffy.hh"
#include "scala_util.hh"
#include "hkl_datatypes.hh"
#include "scaletypes.hh"
#include "string_util.hh"
#include "report_errors.hh"

// Clipper
#include <clipper/clipper.h>
#include "clipper/clipper-ccp4.h"
using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;

namespace phaser_io {

//--------------------------------------------------------------
  void ReportSyntaxError(const std::string& keywords, const std::string& message)
  {
    ReportErrors::printFatalError
      ("Syntax error, keywords:\n"+keywords+"\n"
                               +"Message: "+message+"\n");
    throw SyntaxError(keywords, message);
  }
//--------------------------------------------------------------
ANOMALOUS::ANOMALOUS() : CCP4base(), InputBase()
{
  Add_Key("ANOM");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);

  anomalous = false;
  given = false;
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
        {ReportSyntaxError
            (keywords, "key not ON or OFF");}
      given = true;
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
  // Set up default SCALES parameters, see scaletypes.hh
  scala::ScaleSpecification spec_default;
  specs.push_back(spec_default);  // and store it
  notile = false;
}
//--------------------------------------------------------------
Token_value SCALES::parse(std::istringstream& input_stream)
// SCALES [RUN <irun>]
// [BATCH || ROTATION [<nscales> || SPACING <spacing>]
// BFACTOR ON || OFF BROTATION [<nscales> || SPACING <spacing>]
// [SECONDARY  [<Lmax> [<LmaxOdd>]]]
// [ABSORPTION [<Lmax> [<LmaxOdd>]] [POLE [h|k|l]]]
// [CONSTANT]
// [[NO]TILE [<Ntilex> [<Ntiley>]] [CCD[n] | FLAT | PIXEL]]
//    CCD1, CCD2, CCD3 are 3 different tile models, default CCD2 if just CCD
{
  // Read one SCALES specification
  int irun = -1;
  scala::ScaleSpecification spec;
  int expectingNumber = 0; // = 0 not expecting number, +1 expecting number
                           // = -1 maybe expecting number
  int batch_spec = 0;
  int bfac_spec = 0;
  int secabs = 0;       // SECONDARY or ABSORPTION given
                        // = +1 looking for Lmax, = +2 Lmax read, = +3 LmaxOdd read
  int seclmaxset = 0;    // = 0 not set, = +1 lmax set, = -2 lmax & lmaxodd set
  int sectype = 0;       // = 0 None, = +1 SECONDARY, -1 ABSORPTION
  int tile = -1;

  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (expectingNumber > 0) {ReportSyntaxError
          (keywords, "SCALES: expecting number not "+string_value);}
      secabs = 0;
      if (keyIs("RUN")) {
        // RUN <Irun>, next token must be number
        irun = 0;
        expectingNumber = +1;
      } else if (keyIs("CONSTANT")) {
        spec.SetConstant(irun);
      } else if (keyIs("BATCH")) {
        if (batch_spec != 0) {ReportSyntaxError
            (keywords, "SCALES: can't have BATCH & ROTATION)");}
        spec.batch = true;
        batch_spec = -1;
        expectingNumber = 0;
      } else if (keyIs("ROTATION")) {
        if (batch_spec < 0) {ReportSyntaxError
            (keywords, "SCALES: can't have BATCH & ROTATION)");}
        spec.batch = false;
        batch_spec = +1;
        expectingNumber = -1;
      } else if (keyIs("SPACING") &&  bfac_spec == 0) {
        if (batch_spec < 0) {ReportSyntaxError
            (keywords, "SCALES: can't have BATCH & ROTATION)");}
        if (batch_spec == +2) {ReportSyntaxError
            (keywords, "SCALES: can't have ROTATION number & SPACING)");}
        spec.batch = false;
        batch_spec = +3;
        expectingNumber = +1;
      } else if (keyIs("BFACTOR")) {
        // BFACTOR default to ON (-1)
        expectingNumber = 0;
      } else if (keyIs("ON")) {
        if (spec.nbfac == -1) {
          spec.nbfac = -2;  // switch on B-factors
        }
      } else if (keyIs("OFF")) {
        spec.nbfac = 0;  // switch off B-factors
      } else if (keyIs("BROTATION")) {
        if (batch_spec < 0) {ReportSyntaxError
            (keywords, "SCALES: can't have BATCH & BROTATION)");}
        spec.batch = false;
        bfac_spec = +1;
        expectingNumber = -1;
      } else if (keyIs("SPACING")  &&  bfac_spec > 0) {
        if (bfac_spec == +2) {ReportSyntaxError
            (keywords, "SCALES: can't have BROTATION number & SPACING)");}
        spec.batch = false;
        bfac_spec = +3;
        expectingNumber = +1;
      } else if (keyIs("SECONDARY")) {
        spec.sec_abs = scala::SecondaryScale::SECONDARY;
        expectingNumber = -1;
        secabs = +1;
        sectype = +1;
      } else if (keyIs("ABSORPTION")) {
        spec.sec_abs = scala::SecondaryScale::ABSORPTION;
        expectingNumber = -1;
        secabs = +1;
        sectype = -1;
      } else if (keyIs("POLE")) {
        if (sectype == +1) {ReportSyntaxError
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
          ReportSyntaxError
            (keywords, "SCALES: POLE must be h, k or l)");
        }
      } else if (keyIs("NOTILE")) {
        spec.detectorscaletype = scala::DetectorScale::NONE;
        notile = true;
      } else if (keyIs("TILE")) {
        tile = 0;  // expecting Ntilex next
        spec.detectorscaletype = scala::DetectorScale::AUTOMATIC;
        expectingNumber = -1;
        notile = false;
      } else if (tile >= 0) {
        if (keyIs("FLAT")) {
          spec.detectorscaletype = scala::DetectorScale::FLAT;
        } else if (keyIs("CCD1")) {
          spec.detectorscaletype = scala::DetectorScale::CCD1;
        } else if (keyIs("CCD3")) {
          spec.detectorscaletype = scala::DetectorScale::CCD3;
        } else if (keyIs("CCD2") || keyIs("CCD")) {
          spec.detectorscaletype = scala::DetectorScale::CCD2;
        } else if (keyIs("PIXEL")) {
          spec.detectorscaletype = scala::DetectorScale::PIXEL;
        }
      }
    } else if (tokenIs(1,NUMBER)) {
      if (expectingNumber == 0) {ReportSyntaxError
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
      } else if (secabs == +1) {
        // Read Lmax
        spec.lmax = Nint(number_value);
        //  force even
        spec.lmax = (spec.lmax/2)*2;
        secabs = +2;
        seclmaxset = +1;
        expectingNumber = -1;
      } else if (secabs == +2) {
        // Read LmaxOdd
        spec.lmaxodd = Nint(number_value);
        //  force odd & < lmax
        spec.lmaxodd = ((Min(spec.lmaxodd, spec.lmax)+1)/2)*2-1;
        secabs = +3;
        seclmaxset = -1;
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
      ReportSyntaxError
        (keywords, "SCALES: invalid syntax");
    }
  }
  if (seclmaxset > 0) {
    // set lmaxodd if not defined & lmax is
    spec.lmaxodd = ((spec.lmax+1)/2)*2-1;
  }
  if (secabs != 0 && spec.lmax == 0) {
    spec.sec_abs = scala::SecondaryScale::NONE;
  }

  spec.isdefault = false;  // explicit spec, not default
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
  runsetselection.runsettype = scala::RunSelection::AUTO;
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
          ReportSyntaxError(keywords,"RUN must be followed by a run number");
      }
      Select = 0;
    } else if (tokenIs(1,NAME) && (keyIs("FILE") || keyIs("SERIES"))) {
      // Keyword FILE or SERIES (synonymous)
      fileSeries = Nint(get1num(input_stream));
      if (fileSeries < 0) {
        ReportSyntaxError(keywords,"FILE | SERIES value must be > 0");
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
          ReportSyntaxError(keywords,"unrecognised keyword");
        }
      } else { // number, unexpected
        ReportSyntaxError(keywords,"unexpected number");
      }
    } else if (Select == +1) {
      // Batch selection
      if (tokenIs(1,NAME)) {
        if (keyIs("TO")) {
          // Range
          if (brange) {
            ReportSyntaxError
              (keywords,"only one batch range allowed per RUN command");
          }
          brange = true;
        } else {
          ReportSyntaxError(keywords,"unrecognised keyword");
        }
      } else {
        // gather numbers
        bnum.push_back(Nint(number_value));
      }
    } else if (Select == -1) {
      ReportSyntaxError(keywords,"RUN DATASET option not yet implemented");
    }
  }
  // Line finished, store results
  if (Select == 0) {
    ReportSyntaxError(keywords,"subkeyword BATCH must be given");
    //    ReportSyntaxError(keywords,"subkeyword BATCH or DATASET must be given");
  }
  if (all) {
    runsetselection.batchranges.AddRange(0, 999999, fileSeries, runnum);
    runsetselection.runsettype = scala::RunSelection::EXPLICIT;
  } else if (brange) {
    if (bnum.size() != 2) {
      ReportSyntaxError(keywords,"batch range must be given as 'n1 TO n2'");
    }
    runsetselection.batchranges.AddRange(bnum[0], bnum[1], fileSeries, runnum);
    runsetselection.runsettype = scala::RunSelection::EXPLICIT;
  } else {
    // List, fail must have range
    ReportSyntaxError(keywords,"batch range must be given as 'n1 TO n2'");
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
  double high = input_resolution_range.ResHigh();
  double  low = input_resolution_range.ResLow();
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
        {ReportSyntaxError
            (keywords, "unrecognised keyword");}
    } else if (tokenIs(1,NUMBER)) {
      if (limits) {
        if (nlim == 0)
          {minfraclim = number_value;}
        else if (nlim == 1)
          {maxfraclim = number_value;}
        else
          {ReportSyntaxError
              (keywords, "ERROR in TEST two numbers must be given");}
        nlim++;
      } else if (scalelim) {
        if (nsclim == 0)
          {minscalefrac = number_value;}
        else
          {ReportSyntaxError
              (keywords,
               "ERROR in CORRECT one numbers must be given");}
        nsclim++;
      } else if (gapval) {
        maxgap = Nint(number_value);
        if (maxgap < 0 || maxgap > 5) {
          ReportSyntaxError(keywords,"unreasonable MaxGap ");
        }
      } else
        {ReportSyntaxError
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
          ReportSyntaxError(keywords,"FILE | SERIES value must be > 0");
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
          ReportSyntaxError(keywords,"unrecognised keyword");
        }
      }
    } else if (Select == +1) {
      // Batch selection
      if (tokenIs(1,NAME)) {
        if (keyIs("TO")) {
          // Range
          if (brange) {
            ReportSyntaxError
              (keywords,"only one batch range allowed per EXCLUDE command");
          }
          brange = true;
        } else {
          ReportSyntaxError(keywords,"unrecognised keyword");
        }
      } else {
        // gather numbers
        bnum.push_back(Nint(number_value));
      }
    } else if (Select == -1) {
      ReportSyntaxError(keywords,"EXCLUDE DATASET option not yet implemented");
    }
  }
  // Line finished, store results
  if (Select == 0) {
    ReportSyntaxError(keywords,"subkeyword BATCH must be given");
    //    ReportSyntaxError(keywords,"subkeyword BATCH or DATASET must be given");
  }
  if (brange) {
    if (bnum.size() != 2) {
      ReportSyntaxError(keywords,"batch range must be given as 'n1 TO n2'");
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
  // BINS [RESOLUTION] <NresoBins> INTENSITY <NintensityBins>

  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (keyIs("RESOLUTION")) nrbins = Nint(get1num(input_stream));
      else if (keyIs("INTENSITY")) nibins = Nint(get1num(input_stream));
    } else if (tokenIs(1,NUMBER)) {
      nrbins = Nint(number_value);
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
  set = false;
}
//--------------------------------------------------------------
Token_value REJECT::parse(std::istringstream& input_stream)
{
  //  Outlier rejection criteria
  //
  //  [SCALE|MERGE]  use these values for scaling|merging steps
  //            if not specified, use for both
  //  [COMBINE] compare observations across all datasets
  //  [SEPARATE]  outlier checks only within datasets [default]
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
  //  NONE no outlier rejection

  // default to current values
  float sdrej = outliercontrolsmerge.Reject(scala::ALL).sdrej;
  float sdrej2 = outliercontrolsmerge.Reject(scala::ALL).sdrej2;
  float sdreja = outliercontrolsmerge.Reject(scala::BOTH).sdrej;
  float sdrej2a = outliercontrolsmerge.Reject(scala::BOTH).sdrej2;
  float emax = -1.0;
  int  emaxgiven = -1;
  bool combine = false;
  float batchrejectfactor = -1.0;  // no batch rejection

  int merge = 0;  // expecting values for MERGE && SCALE, = +1 for MERGE, = -1 for SCALE
  // enum Reject2Policy {REJECT, KEEP, REJECTLARGER, REJECTSMALLER};
  scala::RejectFlags::Reject2Policy rej2policy = scala::RejectFlags::KEEP;
  bool anom = false;
  bool first = true;  // first of pair
  bool batchreject = false; // REJECT BATCH
  bool none = false;  // do some rejection

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
      } else if (keyIs("LARGER")) {
        rej2policy = scala::RejectFlags::REJECTLARGER;
      } else if (keyIs("SMALLER")) {
        rej2policy = scala::RejectFlags::REJECTSMALLER;
      } else if (keyIs("ALL")) {
        anom = true;
        first = true;
      } else if (keyIs("EMAX")) {
        emaxgiven = 0;
      } else if (keyIs("BATCH")) {
        batchreject = true;
      } else if (keyIs("NONE")) {
        none = true;
      } else {
        ReportSyntaxError
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
      } else if (emaxgiven >= 0) {
        emax = number_value;
        emaxgiven = +1;
      } else if (batchreject) {
        batchrejectfactor = number_value;
        batchreject = false;
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
      ReportSyntaxError
        (keywords, "REJECT: invalid syntax");
    }
  }

  // Store them
  outliercontrolsscale.Combine() = combine;
  outliercontrolsmerge.Combine() = combine;
  //  SCALE or both
  if (merge <= 0) {
    outliercontrolsscale.SetReject(scala::RejectFlags(sdrej, sdrej2, rej2policy,
                                                      batchrejectfactor), scala::ALL);
    outliercontrolsscale.SetReject(scala::RejectFlags(sdreja, sdrej2a, rej2policy,
                                                      batchrejectfactor), scala::BOTH);
    if (emaxgiven >= 0) outliercontrolsscale.SetEmax(emax);
    if (none) {
      outliercontrolsscale.SetOutlierPolicy(scala::OutlierControl::NOREJECT);
    }
  }
  //  MERGE or both
  if (merge >= 0) {
    outliercontrolsmerge.SetReject(scala::RejectFlags(sdrej, sdrej2, rej2policy,
                                                      batchrejectfactor), scala::ALL);
    outliercontrolsmerge.SetReject(scala::RejectFlags(sdreja, sdrej2a, rej2policy,
                                                      batchrejectfactor), scala::BOTH);
    if (emaxgiven >= 0) outliercontrolsmerge.SetEmax(emax);
    if (none) {
      outliercontrolsmerge.SetOutlierPolicy(scala::OutlierControl::NOREJECT);
    }
  }
  set = true;

  return skip_line(input_stream);
}
//--------------------------------------------------------------
void REJECT::analyse()
{
  // Check that rejects between anomalous related reflections is not
  // more stringent that within
  if (outliercontrolsscale.Reject(scala::BOTH).sdrej != 0.0) {
    float sdrejALL = outliercontrolsscale.Reject(scala::ALL).sdrej;    // within I+ or I-
    float sdrejBOTH = std::abs(outliercontrolsscale.Reject(scala::BOTH).sdrej);  // between I+ and I-
    if (sdrejALL > sdrejBOTH) {
      Message::message(Message_warn
                       (std::string("\nWARNING: on REJECT command\n")+
                        " Scaling outlier rejection limit for all data ("+
                        StringUtil::Strip(StringUtil::ftos(sdrejBOTH, 6, 2)) +
                        ") is tighter than that within I+|I- sets ("+
                        StringUtil::Strip(StringUtil::ftos(sdrejALL, 6, 2))+")\n"+
                        " ==== This is not sensible"));
    }
  }
  if (outliercontrolsmerge.Reject(scala::BOTH).sdrej != 0.0) {
    float sdrejALL = outliercontrolsmerge.Reject(scala::ALL).sdrej;    // within I+ or I-
    float sdrejBOTH = std::abs(outliercontrolsmerge.Reject(scala::BOTH).sdrej);  // between I+ and I-
    if (sdrejALL > sdrejBOTH) {
      Message::message(Message_warn
                       (std::string("\nWARNING: on REJECT command\n")+
                        " Merging outlier rejection limit for all data ("+
                        StringUtil::Strip(StringUtil::ftos(sdrejBOTH, 6, 2)) +
                        ") is tighter than that within I+|I- sets ("+
                        StringUtil::Strip(StringUtil::ftos(sdrejALL, 6, 2))+")\n"+
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
  tiesd_surface = 0.005;
  tiesd_rotation = -1.0;
  tiesd_bfactor = -1.0;
  tiesd_zerob = -1.0;
  // SDs for CCD tile, parameters r, w, A, x0|y0, Fourier components
  // SDs for r, w, xy0 are relative to maximum radius (edge)
  double tie_tile[] = {0.01, 0.01, 0.001, 0.01, 0.002, 0.70, 0.40};
  ties_tile.assign(tie_tile, tie_tile+7);
}
//--------------------------------------------------------------
Token_value TIE::parse(std::istringstream& input_stream)
{
  // TIE <parameter> <sd> [<sd2> etc]
  //  <parameter> = SURFACE    for SECONDARY or ABSORPTION
  //              = ROTATION   for primary scale parameters (eg BATCH)
  //              = BFACTOR    for B-factors
  //              = ZEROB      for B-factors tied to B = 0
  //              = TILE       for tile correction parameters (5 sds)
  //              = TARGETTILE targets for tile correction r & w
  std::vector<double> tsd;      // for tile, up to 5 numbers
  std::vector<double> targets;  // for tile, up to 2 numbers
  int tilesd = -1;
  while (get_token(input_stream) != ENDLINE)  {
    if (tokenIs(1,NAME)) {
      if (keyIs("SURFACE") || keyIs("SECONDARY") || keyIs("ABSORPTION")) {
        tiesd_surface = get1num(input_stream);
      } else if (keyIs("ROTATION") || keyIs("SCALE")) {
        tiesd_rotation = get1num(input_stream);
      } else if (keyIs("BFACTOR")) {
        tiesd_bfactor = get1num(input_stream);
      } else if (keyIs("ZEROB")) {
        tiesd_zerob = get1num(input_stream);
      } else if (keyIs("TILE")) {
        tilesd = 0;
      } else if (keyIs("TARGETTILE")) {
        tilesd = +1;
      }
    } else if (tokenIs(1,NUMBER)) {
      // only for TILE
      if (tilesd < 0) {
        ReportSyntaxError
          (keywords, "TIE: unexpected number when not TILE");
      }
      if (tilesd == 0) {
        tsd.push_back(number_value);
      } else if (tilesd == +1) {
        targets.push_back(number_value);
      }
    }
  }
  ASSERT (ties_tile.size() == 7);
  for (size_t i=0;i<tsd.size();++i) {
    ties_tile[i] = tsd[i];  // override defaults
  }
  if (targets.size() > 0) {
    ASSERT (targets.size() == 2);
    ties_tile[5] = targets[0];
    ties_tile[6] = targets[1];
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
// Syntax: REFINE BFGS|REFERENCE|FH CYCLE [Ncyc1>] <Ncycles> CONVERGE <convergeLimit>
//            SELECT <IovSDmin> <E2min> [<E2max>]
//            PARALLEL [AUTO] | <nproc> | <fproc>
//    CYCLE
//         Ncyc1   number of cycles in 1st stage [default 2]
//         Ncycles number of cycles in main scaling [default 10]
//    BFGS use BFGS minimiser
//    REFERENCE scale to reference dataset
//    FH   use Fox_Holmes least-squares
//    CONVERGE set convergence limit, multiplier of SD (for FH only)
// Selection criteria for scaling reflections:
//    IovSDmin   <I>/sd'(<I>) limit for 1st pass scaling
//    E2min      |E^2| limit for 2nd pass scaling
// If OpenMP is enabled:
//    PARALLEL  number of processors to use in scaling, or
//              fraction of available processors to use, or
//              AUTO determine a "best" number of processors to use
//    If PARALLEL is specified without an argument, then AUTO is assumed
{
  int nn = -1; // counter for CYCLES
  int ns = -1; // counter for SELECT
  int nc = -1; // counter for CONVERGE

  // for NPROC: 0 default, -1 AUTO, +1 waiting for number, +2 number read
  int nproc = 0;
  float fproc = +1.0; // value read for NPROC

  int expectingNumber = -1; // = 0 not expecting number, +1 expecting number
                           // = -1 maybe expecting number

  while (get_token(input_stream) != ENDLINE)  {
    if (tokenIs(1,NAME)) {
      if (expectingNumber > 0) {
        ReportSyntaxError
          (keywords, "syntax error expecting a number");
      }
      if (keyIs("BFGS")) {
        refinecontrol.setMethodBFGS();
        expectingNumber = 0;
      } else if (keyIs("FH")) {
        refinecontrol.setMethodFH();
        expectingNumber = 0;
      } else if (keyIs("REFERENCE")) {
        refinecontrol.setMethodReference();
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
      } else if (keyIs("PARALLEL")) {
        expectingNumber = -1;
        nproc = +1;
      } else if (keyIs("AUTO")) {
        expectingNumber = 0;
        nproc = -1;
        fproc = -1.0;
      }
    } else if (tokenIs(1,NUMBER)) {
      if (expectingNumber == 0) {
        ReportSyntaxError
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
      } else if (nproc == +1) {
        fproc = number_value;
        nproc = +2;
      }
    }
  }
  if (nproc != 0) {
    // NPROC set: 0 default, -1 AUTO, +1 NPROC but no number, +2 number read
    // NPROC set
    if (nproc < 0) {
      fproc = -1.0;  // AUTO, determine later
    } else if (nproc == +1) {
      fproc = -1.0;  // no argument given, assume AUTO
    }
    refinecontrol.SetNprocs(fproc); // set number of processors
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
      ReportSyntaxError("ERROR in BLANK command, Resolutionfraction must be between 0.1 and 1.0",
                        "");
    }}
  if (nullNegativeReject > 0.0) {
    if (nullNegativeReject < 0.01 || nullNegativeReject > 0.501) {
      ReportSyntaxError("ERROR in BLANK command, negative fraction must be between 0.01 and 0.5",
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
  allsame = true;
  fixsdb = false;
  sdinput.clear();
  runnumbers.clear();
  damp = -1.0;
  tietype = -1;  // default tietype
  targets.assign(3,0.0);
  sdtargets.assign(3,0.0);
  weighttype = scala::WeightType::VARIANCE;
  //weighttype = scala::WeightType::SCALE;
  sampleSD = false;
  set = false;
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
//     SIMILAR <sd1> <sd2> <sd3>   for SDfac, [SDb,] SDadd
//     WEIGHT VARIANCE | UNIT | SQRTSCALE | SCALE | SAMPLESD  set weighting scheme
//        for averaging Ih in calculating deviations
//        VARIANCE  w = 1/var(I)  [default]
//        UNIT      w = 1
//        SCALE     w = g = 1/scale
//        SQRTSCALE w = 1/sqrt(g) = sqrt(scale)
//        SAMPLESD use sample SD in final averaging
//
// NB LINEAR and GAUSSIAN options don't work - do not use
{
  int expectingNumber = -1; // = 0 not expecting number, +1 expecting number
                           // = -1 maybe expecting number
  int irun = -1;
  int k = 0;
  bool dampset = false;
  int tieset = -1;
  int similarset = -1;
  double SDfac;
  double SDb;
  double SDadd;
  scala::SDcorrection sdcfull;
  scala::SDcorrection sdcpartial;
  bool found = false;
  bool weight = false;  // WEIGHT keyword found

  int fullpart = 0;  // +1 full, -1 partial, 0 both

  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (expectingNumber > 0) {ReportSyntaxError
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
        refine = +1;  // default refine
        refine_set = true;
      } else if (keyIs("NOREFINE")) {
        refine = 0;
        refine_set = true;
      } else if (keyIs("LINEAR")) {
        refine = -1;
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
        k = 0;
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
        similarset = -1; // TIE & SIMILAR are exclusive
      } else if (keyIs("SIMILAR")) {
        similarset = 0;
        tieset = -1;   // TIE & SIMILAR are exclusive
        expectingNumber = -1;
        tietype = +2;
      } else if (keyIs("NOTIE")) {
        expectingNumber = 0;
        tieset = -1;
        tietype = 0;  // no ties tietype
      } else if (keyIs("SDFAC")) {
        if (tieset != 0) {
          ReportSyntaxError
            (keywords, "SDFAC not expected except after TIE");
        }
        tieset = +1;
        expectingNumber = +1;
      } else if (keyIs("SDB")) {
        if (tieset != 0) {
          ReportSyntaxError
            (keywords, "SDB not expected except after TIE");
        }
        tieset = +2;
        expectingNumber = +1;
      } else if (keyIs("SDADD")) {
        if (tieset != 0) {
          ReportSyntaxError
            (keywords, "SDADD not expected except after TIE");
        }
        tieset = +3;
        expectingNumber = +1;
      } else if (keyIs("WEIGHT")) {
        weight = true;
      } else if (keyIs("VARIANCE")) {
        weighttype = scala::WeightType::VARIANCE;
      } else if (keyIs("UNIT")) {
        weighttype = scala::WeightType::UNIT;
      } else if (keyIs("SQRTSCALE")) {
        weighttype = scala::WeightType::SQRTSCALE;
      } else if (keyIs("SCALE")) {
        weighttype = scala::WeightType::SCALE;
      } else if (keyIs("SAMPLESD")) {
        sampleSD = true;
      } else {
        ReportSyntaxError
          (keywords, "unrecognised keyword");
      }
    } else if (tokenIs(1,NUMBER)) {
      if (expectingNumber == 0) {ReportSyntaxError
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
        expectingNumber = -1;
      } else if (similarset >= 0) {
        // target sd values
        sdtargets[similarset] = number_value;
        similarset++;
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
          ReportSyntaxError
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
    ReportSyntaxError
      (keywords, "SDCORRECTION:: TIE needs two values for target and weight");
  }

  if (similarset >=0) {
    // should be either = 0, use defaults, or = 3
    if (similarset == 0) {
      // use defaults
      sdtargets[0] = 0.2;
      sdtargets[1] = 3.0;
      sdtargets[2] = 0.04;
    } else if (!(similarset == 3)) {
      ReportSyntaxError
        (keywords, "SDCORRECTION:: SIMILAR needs 0 or 3 numbers for target SDs");
    }
    // Set some starting target defaults
    targets[0] = 1.5;
    targets[1] = 0.0;
    targets[2] = 0.04;
    tietype = +2;
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
  set = true;
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
//! return ties, target & SDs, = 0 no tie, = -1 defaults, = +1 set here, +2 similarity
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
  if (SDC_NumberInput() > 1 &&  allsame == true) {
    ReportSyntaxError
      (keywords, "SDCORRECTION:: multiple values given with SAME flag");
  }
  if (refine_set) return;  // explicit refine flag set
  if (SDC_NumberInput() != 0) {
    refine = 0;
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
      if (expectingNumber > 0) {ReportSyntaxError
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
        ReportSyntaxError
          (keywords, "unrecognised keyword");
      }
    } else if (tokenIs(1,NUMBER)) {
      if (expectingNumber == 0) {ReportSyntaxError
          (keywords, "INTENSITIES: not expecting number"+
           clipper::String(number_value));}
      if (type == 0) {
        // COMBINE set Imid
        imid = number_value;
      } else if (type == +1) {
        // POWER set power
        ipowercomb = Nint(number_value);
        if (ipowercomb < 1 || ipowercomb > 5) {
          ReportSyntaxError
            (keywords, "INTENSITIES:: unreasonable POWER value");
        }
      } else {
          ReportSyntaxError
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
      if (expectingNumber > 0) {ReportSyntaxError
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
      } else if (keyIs("MISFIT")) {
        observationflagcontrol.SetAcceptMisfit();
        expectingNumber = 0;
      }
    } else if (tokenIs(1,NUMBER)) {
      if (expectingNumber == 0) {ReportSyntaxError
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
          ReportSyntaxError
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
//        [SCALEPACK (aka POLISH) MERGED | UNMERGED]
//        [ORIGINAL | REDUCED]
//
// Scalepack output is always split
{
  bool mtz = true;  // reading keywords for MTZ
  bool merged = true;  // reading keywords for MERGED
  int isplitmtzmerged = 0;   // no default here, +1 split, -1 together
  int isplitmtzunmerged = 0; // no default here

  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (keyIs("MERGED") || keyIs("AVERAGE")) {
        if (mtz) {
          outputcontrols.SetMTZoutputType(scala::OutputControls::MERGED);
        } else {
          outputcontrols.SetSCAoutputType(scala::OutputControls::MERGED);
        }
        merged = true;
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
        merged = false;
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
        if (merged) { // MTZ MERGED SPLIT
          isplitmtzmerged = +1;
        } else { // MTZ UNMERGED SPLIT
          isplitmtzunmerged = +1;
        }
      } else if (keyIs("TOGETHER")) {
        if (merged) { // MTZ MERGED SPLIT
          isplitmtzmerged = -1;
        } else { // MTZ UNMERGED SPLIT
          isplitmtzunmerged = -1;
        }
      } else if (keyIs("ORIGINAL")) {
        outputcontrols.originalHKL() = true;
      } else if (keyIs("REDUCED")) {
        outputcontrols.originalHKL() = false;
      } else {
          ReportSyntaxError
            (keywords, "OUTPUT: unrecognised keyword");
      }
    }
  }
  if (isplitmtzmerged != 0) {
    outputcontrols.SplitMerged() = (isplitmtzmerged > 0);  // true if +1 SPLIT
    if (isplitmtzmerged < 0) {
      ReportSyntaxError
            (keywords, "OUTPUT MERGED TOGETHER option not yet available");
    }
  }
  if (isplitmtzunmerged != 0) {
    outputcontrols.SplitUnmerged() = (isplitmtzunmerged > 0);  // true if +1 SPLIT
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
                        minimumhalfdatasetcc(0.3),
                        minimumhalfdatasetanomcc(0.15),
                        minimumioversigma(1.5),
                        minimumbatchioversigma(1.0),
                        smoothstatisticsrange(-1.0),
                        batchgrouprange(1.0),
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
  //    CCANOMMINIMUM <MinimumHalfdatasetAnomCC>
  //    ISIGMINIMUM <MinimumIoverSigma>
  //    BATCHISIGMINIMUM <MinimumBatchIoverSigma>
  //    SMOOTHSTATISTICS <SmoothStatisticsRange>
  //    GROUPBATCH  <BatchGroupRange>
  //
  // Cone angle is the half-angle (degrees) for cones around each reciprocal axis
  // MinimumHalfdatasetCC  minimum CC for resolution warning
  // MinimumIoverSigma          minimum <<I>/sd(<I>)> for resolution warning
  // MinimumBatchIoverSigma     minimum <I/sd(I)> for resolution warning by batch, from unmerged I
  // SmoothStatisticsRange angle in degrees over which (roughly) to smooth
  //            batch statistics, <0 to default to automatic setting
  // BatchGroupRange        phi range for grouping batches in analysis, < 0 individual

  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (keyIs("CONE")) {coneangledegrees = get1num(input_stream);}
      else if (keyIs("CCMINIMUM")) {minimumhalfdatasetcc = get1num(input_stream);}
      else if (keyIs("CCANOMMINIMUM")) {minimumhalfdatasetanomcc = get1num(input_stream);}
      else if (keyIs("ISIGMINIMUM")) {minimumioversigma = get1num(input_stream);}
      else if (keyIs("BATCHISIGMINIMUM")) {minimumbatchioversigma = get1num(input_stream);}
      else if (keyIs("SMOOTHSTATISTICS")) {smoothstatisticsrange = get1num(input_stream);}
      else if (keyIs("GROUPBATCH")) {batchgrouprange = get1num(input_stream);}
      else if (keyIs("DETECTOR")) {detector = true;}
      else if (keyIs("NODETECTOR")) {detector = false;}
      else {
          ReportSyntaxError(keywords,"unrecognised keyword");
      }
    }
  }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
INITIAL::INITIAL()  : unity(false), minimum_overlap(-0.05), maximum_gap(2)
{
  Add_Key("INITIAL");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value INITIAL::parse(std::istringstream& input_stream)
{
  // Syntax:
  //  INITIAL UNITY  set all initial scales to unity
  //  INITIAL MEAN   set all initial scales from mean intensities [default]
  //  INITIAL MINIMUM_OVERLAP  read minimum overlap fraction
  //  INITIAL MAXIMUM_GAP  read maximum contiguous "gaps" in rotation ranges

  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (keyIs("UNITY")) {
        unity = true;
      } else if (keyIs("MEAN")) {
        unity = false;
      } else if (keyIs("MINIMUM_OVERLAP")) {
        minimum_overlap = get1num(input_stream);
      } else if (keyIs("MAXIMUM_GAP")) {
        maximum_gap = Nint(get1num(input_stream));
      }
    }
  }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
XMLOUT::XMLOUT() : CCP4base(), InputBase()
{
  Add_Key("XMLOUT");
  name = "";
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value XMLOUT::parse(std::istringstream& input_stream)
{
  name = StringUtil::Unquote(getLine(input_stream));
  return ENDLINE;
}
//--------------------------------------------------------------
ROGUES::ROGUES() : CCP4base(), InputBase()
{
  Add_Key("ROGUES");
  name = "";
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value ROGUES::parse(std::istringstream& input_stream)
{
  name = StringUtil::Unquote(getLine(input_stream));
  return ENDLINE;
}
//--------------------------------------------------------------
HKLIN::HKLIN() : CCP4base(), InputBase()
{
  Add_Key("HKLIN");
  name = "";
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value HKLIN::parse(std::istringstream& input_stream)
{
  name = StringUtil::Unquote(getLine(input_stream));
  return ENDLINE;
}
//--------------------------------------------------------------
HKLOUT::HKLOUT() : CCP4base(), InputBase()
{
  Add_Key("HKLOUT");
  name = "";
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value HKLOUT::parse(std::istringstream& input_stream)
{
  name = StringUtil::Unquote(getLine(input_stream));
  return ENDLINE;
}
//--------------------------------------------------------------
UNMERGEDOUT::UNMERGEDOUT() : CCP4base(), InputBase()
{
  Add_Key("UNMERGEDOUT");
  name = "";
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value UNMERGEDOUT::parse(std::istringstream& input_stream)
{
  name = StringUtil::Unquote(getLine(input_stream));
  return ENDLINE;
}
//--------------------------------------------------------------
HKLREF::HKLREF() : CCP4base(), InputBase()
{
  Add_Key("HKLREF");
  name = "";
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value HKLREF::parse(std::istringstream& input_stream)
{
  name = StringUtil::Unquote(getLine(input_stream));
  return ENDLINE;
}
//--------------------------------------------------------------
bool SetLabel(const std::string& lkey,
              const std::string& label,
              int& nlab,
              std::string& Ilabel,
              std::string& sigIlabel)
// return false if error
{
  if ((lkey == "F" || lkey == "I") ||
      (lkey == "" && nlab == 0)) {
    // store I|F label
    Ilabel = label;
  } else if ((lkey == "SIGF" || lkey == "SIGI") ||
             (lkey == "" && nlab == 1)) {
    // store sigI|F label
    sigIlabel = label;
  } else {
    return false;
  }
  nlab += 1;
  return true;
}
//--------------------------------------------------------------
LABREF::LABREF() : CCP4base(), InputBase()
{
  Add_Key("LABR");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value LABREF::parse(std::istringstream& input_stream)
{
  // Syntax: LABREF [F|I = ] <F|Ilabel> [[SIGF|I = ] <sigF|Ilabel>]
  int nlab = 0;
  int proglab = -1; // Initially no "program label" read
  std::string label;
  std::string lkey;

  while (get_token(input_stream) != ENDLINE) {
    if (proglab >= 0) {
      // Have potential program label, is next field "="?
      if (curr_tok == ASSIGN) {
        lkey = label;
        label = "";
        proglab = +1;
      } else {
        // No, store as actual label
        if (proglab > 0) label = string_value;
        if (! SetLabel(lkey, label, nlab, FIlabel, sigFIlabel))
          throw SyntaxError(keywords,"Too many column labels");
        if (proglab > 0) {
          // "=" was read, just had label
          proglab = -1;
        } else {
          // current string is next key or label
          proglab = 0;
          label = string_value;
          lkey = "";
        }
      }
    } else {
      label = string_value;
      lkey = "";
      proglab = 0;
    }
  }
  if (proglab == 0)
    if (! SetLabel(lkey, label, nlab, FIlabel, sigFIlabel))
      throw SyntaxError(keywords,"Too many column labels");
  if (sigFIlabel != "" && FIlabel == "") {
      throw SyntaxError(keywords,"Can't specify just SIG");
  }

  return skip_line(input_stream);
}
//--------------------------------------------------------------
XYZIN::XYZIN() : CCP4base(), InputBase()
{
  Add_Key("XYZIN");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value XYZIN::parse(std::istringstream& input_stream)
{
  name = StringUtil::Unquote(getLine(input_stream));
  return ENDLINE;
}
//--------------------------------------------------------------
USESDPARAMETER::USESDPARAMETER() : CCP4base(), InputBase()
{
  Add_Key("USESDPARAMETER");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
  parametersdusage = scala::ScaleSpecification::DIAGONAL;
}
//--------------------------------------------------------------
Token_value USESDPARAMETER::parse(std::istringstream& input_stream)
// Syntax: USESDPARAMETER [NO | DIAGONAL | COVARIANCE]
//    (default DIAGONAL if not explicit)
{
  bool OK = true;
  bool nokey = true;
  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (keyIs("NO")) {
        parametersdusage = scala::ScaleSpecification::NONE;
        nokey = false;
      } else if (keyIs("DIAGONAL")) {
        parametersdusage = scala::ScaleSpecification::DIAGONAL;
        nokey = false;
      } else if (keyIs("COVARIANCE")) {
        parametersdusage = scala::ScaleSpecification::COVARIANCE;
        nokey = false;
      } else {
        OK = false;
      }
    } else {
      OK = false;
    }
  }
  if (!OK) {
    throw SyntaxError(keywords,
      "Unrecognised keyword, should be NO | DIAGONAL | COVARIANCE");
  }
  if (nokey) { // default
    parametersdusage = scala::ScaleSpecification::DIAGONAL;
  }

  return ENDLINE;
}
//--------------------------------------------------------------
LINK::LINK() : CCP4base(), InputBase()
{
  Add_Key("LINK");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
  Add_Key("UNLINK");
  //Add to CCP4base;
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value LINK::parse(std::istringstream& input_stream)
// Syntax: [UN]LINK [SURFACE] [ALL] | <run2> TO <run1>
{
  std::string command = stoup(string_value);
  bool link;
  if (command == "LINK") {
    link = true;
  } else if (command == "UNLINK") {
    link = false;
  }
  bool surface = true;  // for now always surface
  bool all = false;

  int expectingNumber = -1; // = 0 not expecting number, +1 expecting number
                           // = -1 maybe expecting number

  std::vector<int> runs;

  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (keyIs("SURFACE")) {
        surface = true;
      } else if (keyIs("ALL")) {
        all = true;
        expectingNumber = 0;
      } else if (keyIs("TO")) {
        all = false;
      }
    } else if (tokenIs(1,NUMBER)) {
      runs.push_back(Nint(number_value));
    }
  }

  // if no runs specified assume ALL
  if (runs.size() == 0) {
    all = true;
  } else {
    if (runs.size() != 2) {
      ReportSyntaxError(keywords, "should be LINK <run2> TO <run1>");
    }
  }

  if (all) {
    // link or unlink all
    if (link) {
      linkspecs.setLinkAll();
    } else {
      linkspecs.setUnlinkAll();
    }
  } else {
    std::pair<int, int> runs2(runs.at(0), runs.at(1));
    if (link) {
      if (linkspecs.linkAll()) {
        ReportSyntaxError(keywords, "can't add explicit LINKs to 'LINK ALL'");
      }
      linkspecs.addLink(runs2);
    } else {
      if (linkspecs.unlinkAll()) {
        ReportSyntaxError(keywords, "can't add explicit UNLINKs to 'UNLINK ALL'");
      }
      // may have LINK ALL (linkall true)
      linkspecs.addUnlink(runs2);
    }
  }

  return ENDLINE;
}
//--------------------------------------------------------------
PLOT::PLOT() : CCP4base(), InputBase()
{
  Add_Key("PLOT");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
  xmgroutput = true;
}
//--------------------------------------------------------------
Token_value PLOT::parse(std::istringstream& input_stream)
// Option to suppress XMgrace file output
//    NORMPLOT, ANOMPLOT, ROGUEPLOT, CORRELPLOT
//
// Syntax: PLOT [NOXMGR]
{
  bool OK = true;
  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (keyIs("NOXMGR")) {
        xmgroutput = false;
      } else if (keyIs("XMGR")) {
        xmgroutput = true;
      } else {
        OK = false;
      }
    } else {
      OK = false;
    }
  }
  if (!OK) {
    throw SyntaxError(keywords,
                      "Unrecognised keyword, should be NOXMGR | XMGR");
  }

  return ENDLINE;
}
//--------------------------------------------------------------
CELL::CELL() : CCP4base(), InputBase()
{
  Add_Key("CELL");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value CELL::parse(std::istringstream& input_stream)
{
  int i = 0;
  std::vector<Dtype> vcell(6);
  while (get_token(input_stream) != ENDLINE)  {
    if (tokenIs(1,NUMBER)) {
      vcell[i++] = number_value;
    }
  }
  if (i != 6) {
    ReportSyntaxError(keywords,
                      "ERROR in command CELL: six numbers must be given");
  }
  cell = scala::Scell(vcell);
  return skip_line(input_stream);
}
//--------------------------------------------------------------
} // phaser_io
