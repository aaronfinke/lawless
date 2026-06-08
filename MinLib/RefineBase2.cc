//(c) 2000-2007 Cambridge University Technical Services Ltd
//All rights reserved
#include <limits>
#include "RefineBase2.h"

namespace phaser {

floatType RefineBase2::finiteDiffGradient(TNT::Vector<floatType>& Gradient, floatType fracLarge)
// Fraction of large shift to use for each type of function should be
// investigated by numerical tests, varying by, say, factors of two.
{
  floatType f,fplus,sz;
  TNT::Vector<floatType> x(npars_ref),old_x(npars_ref);
  TNT::Vector<floatType> largeShifts = getLargeShifts();
  reparLargeShifts(largeShifts);
  f = targetFn();
  TNT::Vector<floatType> unrepar_oldx = getRefinePars();
  old_x = reparRefinePars(unrepar_oldx);
  x = old_x;
  for (int i = 0; i < npars_ref; i++)
  {
    sz = fracLarge*largeShifts[i];
    x[i] = old_x[i] + sz;
    TNT::Vector<floatType> unrepar_x = reparRefineParsInv(x);
    applyShift(unrepar_x);
    fplus = targetFn();
    x[i] = old_x[i];
    Gradient[i] = (fplus - f)/sz;
  }
  unrepar_oldx = reparRefineParsInv(old_x);
  applyShift(unrepar_oldx);
  return f;
}

floatType RefineBase2::finiteGDiffHessian(TNT::Fortran_Matrix<floatType>& Hessian, floatType fracLarge)
// Fraction of large shift to use for each type of function should be
// investigated by numerical tests, varying by, say, factors of two.
{
  floatType f,sz;
  Hessian.newsize(npars_ref,npars_ref);
  TNT::Vector<floatType> g(npars_ref),gplus(npars_ref),unrepar_g(npars_ref);
  TNT::Vector<floatType> x(npars_ref),old_x(npars_ref),unrepar_oldx(npars_ref);
  TNT::Vector<floatType> largeShifts = getLargeShifts();
  reparLargeShifts(largeShifts);
  for (int i = 1; i <= npars_ref; i++)
  for (int j = 1; j <= npars_ref; j++)
    Hessian(i,j) = 0.;
  unrepar_oldx = getRefinePars();
  old_x = reparRefinePars(unrepar_oldx);
  x = old_x;
  f = gradientFn(unrepar_g);
  g = reparGradient(unrepar_oldx,unrepar_g);
  for (int i = 0; i < npars_ref; i++)
  {
    sz = fracLarge*largeShifts[i];
    x[i] = old_x[i] + sz;
    TNT::Vector<floatType> unrepar_x = reparRefineParsInv(x);
    applyShift(unrepar_x);
    gradientFn(unrepar_g);
    gplus = reparGradient(unrepar_x,unrepar_g);
    x[i] = old_x[i];
    for (int j = 0; j < npars_ref; j++)
    {
      if (i == j)
      {
        Hessian(i+1,i+1) = (gplus[i] - g[i])/sz;
      }
      else
      {
        Hessian(i+1,j+1) += (gplus[j] - g[j])/(2*sz);
        Hessian(j+1,i+1) += (gplus[j] - g[j])/(2*sz);
      }
    }
  }
  unrepar_oldx = reparRefineParsInv(old_x);
  applyShift(unrepar_oldx);
  return f;
}

floatType RefineBase2::finiteGDiffDiagHessian(TNT::Fortran_Matrix<floatType>& Hessian, floatType fracLarge)
// Fraction of large shift to use for each type of function should be
// investigated by numerical tests, varying by, say, factors of two.
// This routine is inefficient, costing as much as finiteGDiffHessian
// to compute full Hessian, but may be useful for testing purposes.
{
  Hessian.newsize(npars_ref,npars_ref);
  TNT::Vector<floatType> g(npars_ref),gplus(npars_ref),unrepar_g(npars_ref);
  TNT::Vector<floatType> largeShifts = getLargeShifts();
  reparLargeShifts(largeShifts);
  for (int i = 1; i <= npars_ref; i++)
    for (int j = 1; j <= npars_ref; j++)
      Hessian(i,j) = 0.;
  TNT::Vector<floatType> unrepar_oldx = getRefinePars();
  TNT::Vector<floatType> old_x = reparRefinePars(unrepar_oldx);
  TNT::Vector<floatType> x = old_x;
  floatType f = gradientFn(unrepar_g);
  g = reparGradient(unrepar_oldx,unrepar_g);
  for (int i = 0; i < npars_ref; i++)
  {
    floatType sz = fracLarge*largeShifts[i];
    x[i] = old_x[i] + sz;
    TNT::Vector<floatType> unrepar_x = reparRefineParsInv(x);
    applyShift(unrepar_x);
    gradientFn(unrepar_g);
    gplus = reparGradient(unrepar_x,unrepar_g);
    gradientFn(gplus);
    x[i] = old_x[i];
    Hessian(i+1,i+1) = (gplus[i] - g[i])/sz;
  }
  unrepar_oldx = reparRefineParsInv(old_x);
  applyShift(unrepar_oldx);
  return f;
}

floatType RefineBase2::finiteFDiffHessian(TNT::Fortran_Matrix<floatType>& Hessian, floatType fracLarge)
// Fraction of large shift to use for each type of function should be
// investigated by numerical tests, varying by, say, factors of two.
{
  floatType f,fplusi,fminusi,fplusij,fplusj,szi,szj;
  Hessian.newsize(npars_ref,npars_ref);
  TNT::Vector<floatType> largeShifts = getLargeShifts();
  reparLargeShifts(largeShifts);
  for (int i = 1; i <= npars_ref; i++)
  for (int j = 1; j <= npars_ref; j++)
    Hessian(i,j) = 0.;
  f = targetFn();
  TNT::Vector<floatType> unrepar_oldx = getRefinePars();
  TNT::Vector<floatType> old_x = reparRefinePars(unrepar_oldx);
  TNT::Vector<floatType> x = old_x;
  for (int i = 0; i < npars_ref; i++)
  {
    szi = fracLarge*largeShifts[i];
    x[i] = old_x[i] - szi;
    TNT::Vector<floatType> unrepar_x = reparRefineParsInv(x);
    applyShift(unrepar_x);
    fminusi = targetFn();
    x[i] = old_x[i] + szi;
    unrepar_x = reparRefineParsInv(x);
    applyShift(unrepar_x);
    fplusi = targetFn();
    x[i] = old_x[i];
    Hessian(i+1,i+1) = (fplusi - 2*f + fminusi)/(szi*szi);
    for (int j = i+1; j < npars_ref; j++)
    {
      x[i] = old_x[i] + szi;
      szj = fracLarge*largeShifts[j];
      x[j] = old_x[j] + szj;
      unrepar_x = reparRefineParsInv(x);
      applyShift(unrepar_x);
      fplusij = targetFn();
      x[i] = old_x[i];
      unrepar_x = reparRefineParsInv(x);
      applyShift(unrepar_x);
      fplusj = targetFn();
      x[j] = old_x[j];
      Hessian(i+1,j+1) = Hessian(j+1,i+1) = (fplusij - fplusi - fplusj + f)/(szi*szj);
    }
  }
  unrepar_oldx = reparRefineParsInv(old_x);
  applyShift(unrepar_oldx);
  return f;
}

floatType RefineBase2::finiteFDiffDiagHessian(TNT::Fortran_Matrix<floatType>& Hessian, floatType fracLarge)
// Fraction of large shift to use for each type of function should be
// investigated by numerical tests, varying by, say, factors of two.
{
  floatType f,fplus,fminus,sz;
  Hessian.newsize(npars_ref,npars_ref);
  TNT::Vector<floatType> x(npars_ref),old_x(npars_ref);
  TNT::Vector<floatType> largeShifts = getLargeShifts();
  reparLargeShifts(largeShifts);
  for (int i = 1; i <= npars_ref; i++)
    for (int j = 1; j <= npars_ref; j++)
      Hessian(i,j) = 0.;
  f = targetFn();
  TNT::Vector<floatType> unrepar_oldx = getRefinePars();
  old_x = reparRefinePars(unrepar_oldx);
  x = old_x;
  for (int i = 0; i < npars_ref; i++)
  {
    sz = fracLarge*largeShifts[i];
    x[i] = old_x[i] + sz;
    TNT::Vector<floatType> unrepar_x = reparRefineParsInv(x);
    applyShift(unrepar_x);
    fplus = targetFn();
    x[i] = old_x[i] - sz;
    unrepar_x = reparRefineParsInv(x);
    applyShift(unrepar_x);
    fminus = targetFn();
    x[i] = old_x[i];
    Hessian(i+1,i+1) = (fplus - 2*f + fminus)/(sz*sz);
  }
  unrepar_oldx = reparRefineParsInv(old_x);
  applyShift(unrepar_oldx);
  return f;
}

} //phaser
