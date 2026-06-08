//
// scalerefine.cpp
//
// Main scaling routine
//


#include "scalerefine.hh"
#include "refinescale.hh"
#include "hessian.hh"
#include "report_errors.hh"
#include "string_util.hh"

#include <assert.h>
#define ASSERT assert



namespace scala {
// ---------------------------------------------------------
  void ScaleRefine(hkl_unmerge_list& hkl_list,
                   const ReferenceList& hklreflist,
                   ScaleModel& AllScales,
                   const SDmodel& SDM,
                   const all_controls& controls, const int& Ncycles,
                   const bool& print, phaser_io::Output& output)
  //  Main scaling
  // print = true to print scales
  //
  {
    // default protocols
    phaser::protocolPtr cPtr(new phaser::ProtocolScale(Ncycles));
    phaser::Minimizer Min;
    // Set up up refinement object:
    //  store addresses of reflection & scale objects
    if (controls.refinecontrol.BFGS()) {
      RefineScale refscl(hkl_list, AllScales, SDM,
                         controls.refinecontrol.Nprocs());
      Min.run(refscl, cPtr, output);    // run minimiser
      AllScales.NormaliseParameters();  // Normalise result
      // Calculate variance/covariance matrix and store in AllScales
      calculateParameterVariances(refscl, AllScales);
    } else if (controls.refinecontrol.Reference()) {
      RefineScaleRef refscl(hkl_list, hklreflist,
                            AllScales, SDM,
                            controls.refinecontrol.Nprocs());
      Min.run(refscl, cPtr, output);    // run minimiser
      //AllScales.NormaliseParameters();  // Normalise result
      // Calculate variance/covariance matrix and store in AllScales
      calculateParameterVariances(refscl, AllScales);
    }

    // NB this clears the count
    int negscalecount = AllScales.NegativeSecScaleOccurred();
    if (negscalecount > 0) {
      std::string message = StringUtil::itos(negscalecount)+
        " observations had a negative secondary scale, so parameters were scaled down";
      ReportErrors::printWarning(message, "NegativeSecScaleOccurred", true);
    }

    if (print) AllScales.PrintScales(output);

  }  // ScaleRefine
  // ---------------------------------------------------------
  // Calculate variance/covariance matrix and store in AllScales
  void calculateParameterVariances(RefineScale& refscl, ScaleModel& AllScales)
  {
    int minfiltered = AllScales.Nfilter();
    if (minfiltered <= 0) {return;} // no parameters

    TNT::Fortran_Matrix<floatType> HM;
    // Sum(wDel^2) (target) and Hessian
    bool is_diagonal;
    floatType wd2 = refscl.hessianFn(HM, is_diagonal);
    Hessian H(HM);
    //^^
    //    std::cout << "calculateParameterVariances H\n" << H.format() <<"\n"; //^-
    // Invert with scaling and filtering, replace HM
    HM = H.FilteredInverse(minfiltered);

    int nobs = refscl.NobservationsAll();  // number of observations
    int nparams = AllScales.Nparameters(); // number of parameters

    int nminusm = nobs - nparams; // n - m
    if (nminusm < 1) {return;}  // almost no data!

    // Store in scale model
    AllScales.storeParameterVariances(HM, wd2, nminusm);
  }
  // ---------------------------------------------------------
  // Calculate variance/covariance matrix and store in AllScales
  void calculateParameterVariances(RefineScaleRef& refscl,
                                   ScaleModel& AllScales)
  {
    int minfiltered = AllScales.Nfilter();
    if (minfiltered <= 0) {return;} // no parameters

    TNT::Fortran_Matrix<floatType> HM;
    // Sum(wDel^2) (target) and Hessian
    bool is_diagonal;
    floatType wd2 = refscl.hessianFn(HM, is_diagonal);
    Hessian H(HM);
    //^^
    //    std::cout << "calculateParameterVariances H\n" << H.format() <<"\n"; //^-
    // Invert with scaling and filtering, replace HM
    HM = H.FilteredInverse(minfiltered);

    int nobs = refscl.NobservationsAll();  // number of observations
    int nparams = AllScales.Nparameters(); // number of parameters

    int nminusm = nobs - nparams; // n - m
    if (nminusm < 1) {return;}  // almost no data!

    // Store in scale model
    AllScales.storeParameterVariances(HM, wd2, nminusm);
  }
} // namespace scala
