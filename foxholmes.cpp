//  FoxHolmes.cpp
//

// Test routine for simple-minded Fox-Holmes scaling

#define ASSERT assert
#include <assert.h>


#include "foxholmes.hh"

namespace scala {
// ---------------------------------------------------------
FoxHolmes::FoxHolmes(const scala::InitialData& Data)
{
  data = &Data;   // copy pointer to data object
  npar = data->Npar();
  scales = std::vector<double>(npar, 1.0);
  gradient.newsize(npar);
  gradientOK = false;
}
// ---------------------------------------------------------
int FoxHolmes::MeanI(const std::vector<DPair>& y,
		       double& mnI, double& sumwg2) const
// Mean I <I> 
// returns number which have non-zero sd
//
// On exit:
//  mnI   <I>
//  sumwg2  Sum(w g^2)
{
  // <I> = Sum(w g I) / Sum (w g^2)
  double sd, w;
  double sumwgI = 0.0;
  sumwg2 = 0.0;
  int n = 0;
  for (int i=0;i<npar;i++)
    {
      sd = y[i].second;
      if (sd > 0.00001)
	{
	  w = 1./(sd*sd);
	  sumwgI += w * scales[i] * y[i].first;
	  sumwg2 += w * scales[i] * scales[i];
	  n++;
	}
    }
  mnI = 0.0;
  if (sumwg2 > 0.0) 
    {mnI =  sumwgI/sumwg2;}
  return n;
}
// ---------------------------------------------------------
floatType FoxHolmes::targetFn()
{
  // R = 0.5 * Sum( w (I - g<I>)^2
  if (!gradientOK)
    {
      TNT::Fortran_Matrix<floatType> H;
      TargetGradientHessian(false, false, H);
    }
  return target;
}
// ---------------------------------------------------------
floatType    FoxHolmes::gradientFn(TNT::Vector<floatType>& grad)
// dR/dp = Sum [ - w (I - g <I>) d(g<I>)/dp
// d(g<I>)/dp = g d<I>/dp  + <I> dg/dp
// d<I>/dpi = (Ii- 2g<I>)/g^2    in this case (one observation/paraemeter)
{
  if (!gradientOK)
    {
      TNT::Fortran_Matrix<floatType> H;
      TargetGradientHessian(true, false, H);
    }
  grad = gradient;
  //  std::cout << "%% Gradient: \n";
  //  for (int i=0;i<npar;i++) {
  //      std::cout << " " << gradient[i];
  //  }
  //  std::cout << "\n";

  return target;
}
// ---------------------------------------------------------
floatType FoxHolmes::hessianFn(TNT::Fortran_Matrix<floatType>& H,
				 bool& is_diagonal)
{
  is_diagonal = false;
  TargetGradientHessian(true, true, H);
  return target;
}
// ---------------------------------------------------------
void FoxHolmes::TargetGradientHessian(bool DoGradient,
					bool DoHessian,
					TNT::Fortran_Matrix<floatType>& H)
// Private function to calculate target function, gradient & Hessian
//  gradient is stored locally, Hessian is returned
// DoHessian implies DoGradient
{
  target = 0.0;
  if (DoHessian)
    {
      DoGradient = true;
      H.newsize(npar,npar);
      for (int i=0;i<npar;i++) { // loop parameters
	for (int j=0;j<npar;j++) // loop parameters
	  {H(i+1,j+1) = 0.0;}
      }
    }
  if (DoGradient) {
    for (int i=0;i<npar;i++) // loop parameters
      {gradient[i] = 0.0;}
  }
  
  double sd, di, w, mnI, sumwg2;
  std::vector<double> dmnIdp(npar);

  std::vector<DPair> y(npar);
  std::vector<double> dmnIgldgi(npar);

  while (data->ObsArray(y)) {  // Loop observations
    // y(npar) is array of I,sigma pairs
    int n = MeanI(y, mnI, sumwg2);  // Mean I with current scales
    if (n > 0) {
      for (int i=0;i<npar;i++) {  // Loop parameters
	sd = y[i].second;
	if (sd > 0.00001) {
	  w = 1./(sd*sd);
	  di = y[i].first - scales[i] * mnI;
	  
	  // target function = 0.5 * Sum(w(I-g<I>)^2)
	  target += 0.5 * w * di * di;
	  if (DoGradient) {	   // d<I>/dgi = w (Ii - 2 gi <I>)/gi^2
	    dmnIdp[i] = w * (y[i].first - 2.*scales[i]*mnI)/sumwg2;
	  }
	} else {
	  if (DoGradient) {dmnIdp[i] = 0.0;}
	}
      }
      
      if (DoGradient)	{
	for (int l=0;l<npar;l++) { // loop l observations
	  sd = y[l].second;
	  if (sd > 0.00001) {
	    for (int i=0;i<npar;i++) { // loop parameters
	      // observation l, parameter i
	      //  d(gl<I>)/dgi = gl d<I>/dgi [+ <I> if l=i]
	      dmnIgldgi[i] = scales[l] * dmnIdp[i];
	      if (i == l) dmnIgldgi[i] += mnI;
	      if (dmnIgldgi[i] != 0.0) {
		w = 1./(sd*sd);
		gradient[i] += - w * (y[l].first - scales[l] * mnI) * dmnIgldgi[i];
			    
		if (DoHessian) {
		  /// Approximation, diagonal matrix
		  ///		  H(i+1,i+1) += w * dmnIgldgi[i] * dmnIgldgi[i];
		  for (int j=0;j<=i;j++) { // loop parameters again
		    H(i+1,j+1) += w * dmnIgldgi[i] * dmnIgldgi[j];
		  }
		}
	      }
	    } // end loop parameters
	  }
	}  // end loop observations in set
      }
    }
  }  // end loop sets

  if (DoHessian) {
    for (int i=0;i<npar-1;i++) {
      for (int j=i+1;j<npar;j++) 
	{H(i+1,j+1) = H(j+1,i+1);}
    }
  }
  //^
  //  std::cout << "** Parameters: ";
  //  for (int i=0;i<npar;i++) {std::cout << scales[i] << " ";}
  //  std::cout <<"\n";
  //  std::cout << "** Target function : " << target << "\n";
  //-!
  //  if (DoGradient) {
  //    std::cout << "** Gradient: \n";
  //    for (int i=0;i<npar;i++){ 
  //      std::cout << " " << gradient[i];
  //    }
  //    std::cout << "\n";
  //  }
  //  if (DoHessian) {
  //    std::cout << "** Hessian: \n";
  //    for (int i=0;i<npar;i++) {
  //      for (int j=0;j<npar;j++) 
  //	{std::cout << " " << H(i+1,j+1);}
  //      std::cout << "\n";
  //    }
  //  }
  if (DoGradient) gradientOK = true;
}
// ---------------------------------------------------------
void FoxHolmes::applyShift(TNT::Vector<floatType>& newg)
{
  //^
  //  std::cout << "<< Apply shift:\n";
  for (int i=0;i<npar;i++)  {
    scales[i] = newg[i];
  //    std::cout << " " << scales[i];
  }
  //  std::cout << "\n";
  gradientOK = false;
}
// ---------------------------------------------------------
void FoxHolmes::logCurrent(outStream where, Output& output)
{
  /*
  output.logTab(1,where,"\nScales:");
  for (int i=0;i<npar;i++) 
    {output.logTabPrintf(2,where," %7.3f", scales[i]);}
  output.logTab(1,where,"\n");
  */
}
// ---------------------------------------------------------
std::string  FoxHolmes::whatAmI(int& i)
{
  return "Scale "+itos(i);
}
// ---------------------------------------------------------
std::vector<bounds>     FoxHolmes::getLowerBounds()
{
  std::vector<bounds > Lower(npar);
  for (int i=0;i<npar;i++) {
    Lower[i].on(0.0);
  }
  return Lower;
}
// ---------------------------------------------------------
std::vector<bounds>     FoxHolmes::getUpperBounds()
{
  std::vector<bounds > Upper(npar);
  for (int i=0;i<npar;i++) {
    Upper[i].off();
  }
  return Upper;
}
// ---------------------------------------------------------
TNT::Vector<floatType> FoxHolmes::getRefinePars()
{
  TNT::Vector<floatType> pars(npar);
  for (int i=0;i<npar;i++) {
    pars[i] = scales[i];
  }
  return pars;
}
// ---------------------------------------------------------
TNT::Vector<floatType> FoxHolmes::getLargeShifts()
{
  TNT::Vector<floatType> large(npar);
  for (int i=0;i<npar;i++) {
    large[i] = 1.0;
  }
  return large;
}
// ---------------------------------------------------------
bool1D FoxHolmes::getRefineMask(protocolPtr protocol)
{
  return bool1D(npar, true);
}
// ---------------------------------------------------------
}
