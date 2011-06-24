//
// scalerefinefh.cpp
//
// Main scaling routine, traditional filtered Newton method (Fox-Holmes)
//


#include "scalerefinefh.hh"
#include "refinescale.hh"
#include <assert.h>
#define ASSERT assert

#include "Minimizer.h"

#include "symm_eigen.h"

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;


namespace scala {
// ---------------------------------------------------------
  Hessian::Hessian(const TNT::Fortran_Matrix<floatType>& H)
  {
    hessian = H;
    ASSERT (H.dim(1) == H.dim(2));
    Np = H.dim(1);  // dimension = number of parameters
  }
// --------------------------------------------------------
  TNT::Vector<floatType> Hessian::Ubs()
  // Return scale vector to condition Hessian matrix
  {
    if (ubs.size() == 0) {
      // not already calculated
      ubs.newsize(Np);
      // For now, just sqrt(diagonal)
      for (int i=0;i<Np;++i) {
	if (hessian(i+1,i+1) < 0.0) {
	  Message::message(Message_fatal("Negative Hessian diagonal"));
	}
	ubs[i] = hessian(i+1,i+1);
	if (ubs[i] > 0.0) {ubs[i] = sqrt(ubs[i]);}
      }
    //^
      //    std::cout << "UBS  ";
      //    for (int i=0;i<Np;++i) {
      //      std::cout << " " << ubs[i];
      //    }
      //    std::cout << "\n";
    //^-
    }

    return ubs;
  }
// ---------------------------------------------------------
  TNT::Fortran_Matrix<floatType> Hessian::FilteredInverse(const int& MinFiltered)
  // Filtered inverse rescaled Hessian
  //   MinFiltered  minimum number of eigenvalues filtered
  {
    int filtered = 0;
    if (MinFiltered > 0) {filtered = 1;}
    //^    std::cout << "Unscaled Hessian:\n" << format() << "\n";

    // Scale Hessian
    for (int i=0;i<Np;++i) {
      if (ubs[i] > 0.0) {
	for (int j=0;j<Np;++j) {
	  if (ubs[j] > 0.0) {
	    hessian(i+1,j+1) /= (ubs[i]*ubs[j]);
	  }
	}
      }
    }
    //    std::cout << "Scaled Hessian:\n" << format() << "\n";
    hessian = SymmetricPseudoinverse<floatType>
      (hessian,filtered,false,MinFiltered).getInv();
    //^
    //^    std::cout << "Number of eigenvalues filtered " << filtered << "\n";
    //^    std::cout << "Inverse Hessian:\n" << format() << "\n";

    // rescale inverse Hessian
    for (int i=0;i<Np;++i) {
      if (ubs[i] > 0.0) {
	for (int j=0;j<Np;++j) {
	  if (ubs[j] > 0.0) {
	    hessian(i+1,j+1) /= (ubs[i]*ubs[j]);
	  }
	}
      }
    }
    //^    std::cout << "Rescaled inverse Hessian:\n" << format() << "\n";

    return hessian;
  }

// ---------------------------------------------------------
  std::string Hessian::format() const
  {
    std::string s("");
    for (int i=0;i<Np;++i) {
      for (int j=0;j<Np;++j) {
	double h = hessian(i+1,j+1);
	s += " "+clipper::String(h,8,3);
      }
      s += "\n";
    }
    return s;
  }
// ---------------------------------------------------------
// ---------------------------------------------------------
// ---------------------------------------------------------

// ---------------------------------------------------------
  void ScaleRefineFH(hkl_unmerge_list& hkl_list, ScaleModel& AllScales,
		     const all_controls& controls, const int& Ncycles,
		     phaser_io::Output& output)
  //  Main scaling
  // 
  {
    // Set up up refinement object:
    //  store addresses of reflection & scale objects 
    RefineScale refscl(hkl_list, AllScales);
    
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
      ubs = H.Ubs();
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
