//
// scalerefinefh.cpp
//
// Main scaling routine, traditional filtered Newton method (Fox-Holmes)
//


#include "scalerefinefh.hh"
#include "refinescale.hh"
#include "Minimizer.h"
#include "hessian.hh"

#include <assert.h>
#define ASSERT assert


// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;


namespace scala {
// ---------------------------------------------------------
  void ScaleRefineFH(hkl_unmerge_list& hkl_list, ScaleModel& AllScales,
                     const all_controls& controls, const int& Ncycles,
                     phaser_io::Output& output)
  //  Main scaling
  //
  {
    // Set up up refinement object:
    //  store addresses of reflection & scale objects
    RefineScale refscl(hkl_list, AllScales, SDmodel(), 1);

    int minFiltered = 2;
    if (AllScales.NBfactors() == 0) minFiltered = 1;  // case of no Bfactors

    floatType target;   // residual
    TNT::Vector<floatType> gradient;
    TNT::Fortran_Matrix<floatType> hessian;
    bool is_diagonal;
    // scale vector (diagonal matrix) to condition matrix, ~= sqrt(diagonal H)
    TNT::Vector<floatType> ubs;

    for (int icyc=0;icyc<Ncycles;++icyc) {  // loop cycles
      // Calculate Hessian & gradient
      target = refscl.hessianFn(hessian, is_diagonal); // store gradient
      refscl.gradientFn(gradient);
      Hessian H(hessian);
      //^
      std::cout << "Cycle " << icyc << "  residual " << target << "\n"; //^

      // Get scales to condition Hessian
      // Get filtered inverse Hessian
      hessian = H.FilteredInverse(minFiltered);
      TNT::Vector<floatType> shift = hessian * gradient; // shift vector
      // Apply shift
      int Np = AllScales.Nparameters();
      std::vector<double> params = AllScales.GetParameters();
      TNT::Vector<floatType> newparams(Np);
      for (int i=0;i<Np;++i) {
        params[i] -= shift[i];
        newparams[i] = params[i];
      }
      refscl.applyShift(newparams);
      // no needed?      AllScales.SetParameters(params, refscl.Nobservations());
      AllScales.NormaliseParameters();  // Normalise result
      //      AllScales.PrintScales(output);
    }
    //    AllScales.NormaliseParameters();  // Normalise result
    AllScales.PrintScales(output);

  }  // ScaleRefine

} // namespace scala
