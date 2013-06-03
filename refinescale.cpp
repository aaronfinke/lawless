//  refinescale.cpp
//

// 
//  Main Fox-Holmes scaling with all parameters, using BFGS minimiser

#include <assert.h>
#define ASSERT assert

#include "refinescale.hh"
#include "tie.hh"
#include "timer.hh"

#if _OPENMP
#include <omp.h>
#endif

namespace scala {
  // ---------------------------------------------------------
  RefineScale::RefineScale(const hkl_unmerge_list& Hkl_list,
			   ScaleModel& Scalemodel,
			   const int& Nprocs)
  {
    hkl_list = &Hkl_list;     // pointer to data list
    scalemodel = &Scalemodel; // pointer to scales object
    npar = scalemodel->Nparameters();
    params = scalemodel->GetParameters();  // initial parameters
    gradient.newsize(npar);
    gradientOK = false;
    nprocs=Nprocs;
  }
  // ---------------------------------------------------------
  // ---------------------------------------------------------
  // ---------------------------------------------------------
  floatType RefineScale::targetFn()
  {
    // R = 0.5 * Sum( w (I - g<I>)^2
    if (!gradientOK)  {
      //^      std::cout << "targetFn NPROCS " << nprocs <<"\n";

      target = 0.0;
  
      double w, mnI;
      std::vector<double> ghl;

      reflection this_refl;
      observation this_obs;
      std::vector<observation> obs_used; // accepted observations for this reflection
      int index;
      double sd;
      int nobs = 0;
      hkl_list->rewind();
      int Nref = hkl_list->num_reflections(); // total in file

      ////////////////////////
      int jref, jr, l, myid;
      int nused;
      double sumwgI, sumwg2, g, di;
      Rtype invresolsq;
      ////////////////////////

      ////////////////////////////////////////
      int jcount=0;
      std::vector<int> jref_array(Nref, 0);
      std::vector<int> nrefparTEMP(npar*nprocs, 0);
 

      jref=-1; // reflection index counter
      while (++jref!=Nref) { // ---- loop reflections
	// This call to get a reflection should be thread-safe
	jr = hkl_list->get_accepted_reflection(jref, this_refl); // jr is accepted reflection index
	if (jr == -1) break; // end of list
	jref=jr;
	jref_array[jcount]=(jr);
	jcount++;
      }
      ////////////////////////////////////////
 
#pragma omp parallel for default(shared) firstprivate(jcount)		\
  private(jref, jr, myid,						\
	  invresolsq, index, this_refl,					\
	  obs_used, ghl,						\
	  this_obs, sd, nused, sumwgI, sumwg2, w, g, mnI, l, di)

      for (jref=0;jref<jcount;++jref) { // ---- loop reflections
	this_refl = hkl_list->get_reflection(jref_array[jref]); // jr is accepted reflection index
	if (jref_array[jref] >= 0) { 

	  myid = 0;
#if _OPENMP
	  if (nprocs==1) {
	    myid=0;
	  } else {
	    myid=omp_get_thread_num();
	  }
#endif

	  invresolsq = this_refl.invresolsq();
     
	  // Get Mean I & sum wg^2 with current scale,
	  // store intensities & derivatives for each observation
	  obs_used.clear();
	  ghl.clear();
	  nused = 0;
	  sumwgI = 0.0;
	  sumwg2 = 0.0;
    
	  while ((index = this_refl.next_observation(this_obs)) >= 0) {
	    // Loop accepted observations
	    //  Get scale for observation & derivative vector
	    sd = this_obs.sigI();
	    if (sd > 0.00001) {
	      // get scale
	      g = scalemodel->ScaleObs(this_obs, invresolsq);
	      w = 1./(sd*sd);
	      sumwgI += w * g * this_obs.I();
	      sumwg2 += w * g * g;
	      nused++;
	      // Accepted observation
	      obs_used.push_back(this_obs);
	      ghl.push_back(g);
	    }
	  }
	  if (nused > 1) {
	    mnI = 0.0;
	    if (sumwg2 > 0.0)  {mnI =  sumwgI/sumwg2;}  // mnI = <Ih>

	    for (l=0;l!=nused;++l) {  // loop observations
	      sd = obs_used[l].sigI();
	      w = 1./(sd*sd);
	      di = (obs_used[l].I() - ghl[l] * mnI);   // deviation
	      // target function = 0.5 * Sum(w(I-g<I>)^2)
#pragma omp atomic
	      target += 0.5 * w * di * di;
#pragma omp atomic
	      nobs++;

	    }  // end loop observations

	    // got target & d<Ih>/dp
	  }  // end, at least 2 observations
	}
      } // end loop reflections

      // Add ties (restraints) into target etc. Return target

      std::vector<floatType> dRdpi;
      std::vector<TieHessian> Htie;
      floatType Rtie = scalemodel->TieValues(false, false, params,
					     dRdpi, Htie);
    
      target = target + Rtie;
    }
    return target;
  }
  // ---------------------------------------------------------
  floatType    RefineScale::gradientFn(TNT::Vector<floatType>& grad)
  // dR/dp = Sum [ - w (I - g <I>) d(g<I>)/dp
  // d(g<I>)/dp = g d<I>/dp  + <I> dg/dp
  // d<I>/dpi = Sum[w (Ii- 2g<I>)]/Sum w g^2 
  {
    if (!gradientOK) {
      TNT::Fortran_Matrix<floatType> H;
      TargetGradientHessian(true, false, H);
    }
    grad = gradient;
    return target;
  }
  // ---------------------------------------------------------
  floatType RefineScale::hessianFn(TNT::Fortran_Matrix<floatType>& H,
				   bool& is_diagonal)
  {
    is_diagonal = false;
    TargetGradientHessian(true, true, H);
    return target;
  }
  // ---------------------------------------------------------
  void RefineScale::TargetGradientHessian(bool DoGradient,
					  bool DoHessian,
					  TNT::Fortran_Matrix<floatType>& H)
  // Private function to calculate target function, gradient & Hessian
  //  gradient is stored locally, Hessian is returned
  // DoHessian implies DoGradient
  // floatType == double
  //
  //  OpenMP threading added by Ronan Keegan Oct 2011
  // 
  {
    bool DEBUG = false;
    //    bool DEBUG = true;
    //^    std::cout << "NPROCS " << nprocs <<"\n";
    //^
    //^    std::cout << "Npar " << npar <<"\n";

    // it is slightly faster to accumulate Hessian in one-dimensional array Hv and
    // then copy it into the 2D array H
    // Accumulating the Hessian is the rate-limiting step
    std::vector<floatType> Hv(npar*npar*nprocs,0.0);
  
    //  Timer timer;
    target = 0.0;
    if (DoHessian) {
      // Set up & clear Hessian
      DoGradient = true;
      H.newsize(npar,npar);
      for (int i=0;i<npar;i++) { // loop parameters
	for (int j=0;j<npar;j++) // loop parameters
	  {H(i+1,j+1) = 0.0;}
      }
    }
    ASSERT (DoGradient);
    // Clear gradient
    for (int i=0;i<npar;i++) // loop parameters
      {gradient[i] = 0.0;}
  
    double w, mnI;
    std::vector<double> dmnIdp(npar);  // d(<Ih>)/dp  vector
    std::vector<double> dmnIgldp(npar); // d(ghl <Ih>)/dp  should be double
    std::vector<double> dghldp;         // d(ghl)/dp   vector
    std::vector<std::vector<double> > dghldp_obs;         // d(ghl)/dp   vector
    std::vector<double> ghl;
    nrefpar.assign(npar,0);
  
    reflection this_refl;
    observation this_obs;
    std::vector<observation> obs_used; // accepted observations for this reflection
    int index;
    double sd;
    int nobs = 0;
    hkl_list->rewind();
    int Nref = hkl_list->num_reflections(); // total in file
    double wdmnIgldp;
    int ip1;
    int jp;
    // //  while (hkl_list->next_reflection(this_refl) >= 0)  {  // loop reflections
  
    ////////////////////////
    int jref, jr, l, ip, myid;
    int nused;
    double sumwgI, sumwg2, g, di, d;
    Rtype invresolsq;
    ////////////////////////

    ////////////////////////////////////////
    int jcount=0;
    std::vector<int> jref_array(Nref, 0);
    std::vector<int> nrefparTEMP(npar*nprocs, 0);
    std::vector<double> gradientTEMP(npar*nprocs,0.0); 
  
  
    jref=-1; // reflection index counter
    while (++jref!=Nref) { // ---- loop reflections
      // This call to get a reflection should be thread-safe
      jr = hkl_list->get_accepted_reflection(jref, this_refl); // jr is accepted reflection index
      if (jr == -1) break; // end of list
      jref=jr;
      jref_array[jcount]=(jr);
      jcount++;
    }
    ////////////////////////////////////////
 
    //int jref=-1; // reflection index counter
    //while (++jref<Nref) { // ---- loop reflections
    // This call to get a reflection should be thread-safe

#pragma omp parallel for default(shared) firstprivate(jcount, dmnIdp, dmnIgldp) \
  private(jref, jr, myid,						\
	  invresolsq, index, this_refl,					\
	  obs_used, ghl, dghldp_obs, dghldp, wdmnIgldp,			\
	  this_obs, sd, nused, sumwgI, sumwg2, w, g, mnI, l, di, d, ip,	\
	  jp, ip1)

    for (jref=0;jref<jcount;++jref) { // ---- loop reflections
      this_refl = hkl_list->get_reflection(jref_array[jref]); // jr is accepted reflection index
      if (jref_array[jref] >= 0) { 

	myid = 0;
#if _OPENMP
	if (nprocs==1) {
	  myid=0;
	} else {
	  myid=omp_get_thread_num();
	}
#endif

	invresolsq = this_refl.invresolsq();
     
	// Get Mean I & sum wg^2 with current scale,
	// store intensities & derivatives for each observation
	obs_used.clear();
	ghl.clear();
	dghldp_obs.clear();  // derivatives
	nused = 0;
	sumwgI = 0.0;
	sumwg2 = 0.0;

	while ((index = this_refl.next_observation(this_obs)) >= 0) {
	  // Loop accepted observations
	  //  Get scale for observation & derivative vector
	  sd = this_obs.sigI();
	  if (sd > 0.00001) {
	    // get scale and partial derivative vector
	    g = scalemodel->ScaleObs(this_obs, invresolsq, dghldp);
	    ///	ASSERT (int(dghldp.size()) == npar);
	      w = 1./(sd*sd);
	      sumwgI += w * g * this_obs.I();
	      sumwg2 += w * g * g;
	      nused++;
	      //std::cout << "\nhkl, invresolsq " << this_refl.hkl().format() << " " << invresolsq << "\n";
	      //printf("nused=%d sd=%.3lf w=%.3lf g=%.3lf sumwgI=%.3lf this_obs.I()=%.3lf\n"\
	      //            , nused, sd, w, g, sumwgI, this_obs.I());
	      //for (int i=0;i<npar;++i) {std::cout << "  " << dghldp[i];}
	      //   	    std::cout << "\n";
	      // Accepted observation
	      obs_used.push_back(this_obs);
	      ghl.push_back(g);
	      dghldp_obs.push_back(dghldp);
	      //	  if (DEBUG) {
	      //	    std::cout << "\nhkl, Obs, ghl, dghldp " << this_refl.hkl().format()
	      //		      << " " << this_obs.I() << " " << g << "\n";
	      //	    for (int i=0;i<npar;++i) {std::cout << "  " << dghldp[i];}
	      //	    std::cout << "\n";
	      //	  }
	  }
	}
	if (nused > 1) {
	  mnI = 0.0;
	  if (sumwg2 > 0.0)  {mnI =  sumwgI/sumwg2;}  // mnI = <Ih>
	  dmnIdp.assign(npar,0.0);

	  for (l=0;l!=nused;++l) {  // loop observations
	    sd = obs_used[l].sigI();
	    w = 1./(sd*sd);
	    di = (obs_used[l].I() - ghl[l] * mnI);   // deviation
	    // target function = 0.5 * Sum(w(I-g<I>)^2)
#pragma omp atomic
	    target += 0.5 * w * di * di;
#pragma omp atomic
	    nobs++;

	    // d<Ih>/dp = Sum [ whl (d(ghl)/dp) (Ihl - 2 ghl <Ih>)] Sum(whl ghl^2)
	    d = w * (obs_used[l].I() - 2.*ghl[l]*mnI)/sumwg2;
	    for (ip=0;ip!=npar;++ip) {
	      if (dghldp_obs[l][ip] != 0.0) {
		// d<Ih>/dp = Sum [ whl (d(ghl)/dp) (Ihl - 2 ghl <Ih>)] Sum(whl ghl^2)
		dmnIdp[ip] += d * dghldp_obs[l][ip];
	      }
	    }
	  }  // end loop observations
      
	  // got target & d<Ih>/dp
	  for (l=0;l!=nused;++l) { // loop observations in set
	    sd = obs_used[l].sigI();
	    w = 1./(sd*sd);
	    di = (obs_used[l].I() - ghl[l] * mnI);   // deviation
	    //	    if (DEBUG) {std::cout << "\nObs " << l <<"|";}

	    for (ip=0;ip!=npar;++ip) {   // Loop parameters
	      // d(ghl<Ih>)/dp = ghl (dIh/dp)   +   Ih (dghl/dp)
	      dmnIgldp[ip] = ghl[l] * dmnIdp[ip];
	      if (dghldp_obs[l][ip] != 0.0) {
		dmnIgldp[ip] +=  mnI * dghldp_obs[l][ip];
		if (nprocs==1) {
		  nrefpar[ip]++;  // count contributions to this parameter
		} else {
		  nrefparTEMP[ip+(myid*npar)]++;  // count contributions to this parameter
		}
	      }
	      // Gradient vector d = - Sum [ whl delI d(ghl<Ih>)/dp]
	      if (dmnIgldp[ip] != 0.0) {
		wdmnIgldp =  w * dmnIgldp[ip];
		if (nprocs==1) {
		  gradient[ip] += - di * wdmnIgldp;
		} else {
		  gradientTEMP[ip+(myid*npar)] += - di * wdmnIgldp;
		}
		//		if (DEBUG) {
		//		  std::cout << " " <<
		//		    (- w * di * dmnIgldp[ip]) << " " << gradient[ip];
		//		}
	    
		if (DoHessian) {
		  ip1 = ip*npar;
		  for (jp=0;jp<=ip;jp++) { // loop parameters again (to ip)
		    //-//	H(ip+1,jp+1) += w * dmnIgldp[ip] * dmnIgldp[jp];  // half matrix
		    //Hv[ip1+jp] += wdmnIgldp * dmnIgldp[jp];  // half matrix
		    Hv[(ip1+jp)+(myid*npar*npar)] += wdmnIgldp * dmnIgldp[jp];  // half matrix
		  }
		}
	      }
	    } // end loop parameters
	    //	    if (DEBUG) {std::cout << "\n";}
	  }  // end loop observations in set
	}  // end, at least 2 observations
      }
    } // end loop reflections

    // Add ties (restraints) into target etc. Return target

    //    if (DEBUG) {
    //      printf("\n3 >>>>>>>>>>>> target=%.3lf nobs=%d\n",target, nobs);
    //      if (DoGradient) {
    //	for (int i=0;i<npar;i++){
    //	  printf("4 >>>>>>>>>>>> nrefpar[%d]=%d gradient[%d]=%.3lf\n",i, nrefpar[i], i, gradient[i]);
    //	  int iwg = i*(npar+1);
    //	  printf("4 >>>>>>>>>>>> Hv[%d]=%.3lf\n",iwg, Hv[iwg]);
    //	} 
    //	if (DoHessian) {
    //	  printf("Hessian:\n");
    //	  for (int j=0;j<npar;j++){
    //	    for (int i=0;i<=j;i++){
    //	      int iwg = j*npar+i;
    //	      printf(" %.3lf", Hv[iwg]);
    //	    }
    //	    printf("\n");
    //	  }
    //	}
    //      }
    //    }

    std::vector<floatType> dRdpi;
    std::vector<TieHessian> Htie;
    floatType Rtie = scalemodel->TieValues(DoGradient, DoHessian, params,
					   dRdpi, Htie);

    /////////////////////////////////////////////////
      if (DoGradient && nprocs!=1) {
	// Sum the values for gradient and nrefpar from all processes
	for (int i=0;i!=npar*nprocs;++i) { // loop no. parameters * no. of processes
	  gradient[i%npar] += gradientTEMP[i];
	  nrefpar[i%npar] += nrefparTEMP[i];
	}
      }
      /////////////////////////////////////////////////

	if (DoGradient) {
	  ///    ASSERT (int(dRdpi.size()) == npar);
	  for (int ip=0;ip<npar;++ip) {
	    gradient[ip] += dRdpi[ip];
	  }
	  if(DoHessian) {
	    // Symmetrise Hessian
	    //-//      for (int i=0;i<npar-1;i++) {
	    //-//      	for (int j=i+1;j<npar;j++) {
	    //-//      	  H(i+1,j+1) = H(j+1,i+1);} // other half 
	    for (int i=0;i<npar;i++) {
	      for (int j=0;j<=i;j++) {
		for (int proc=0; proc<nprocs; proc++) {
		  H(j+1,i+1) += Hv[(i*npar+j)+(proc*npar*npar)];
		}
		H(i+1,j+1) = H(j+1,i+1); // other half
	      }
	    }
	    if (DEBUG) {
	      printf("Full Hessian:\n");
	      for (int j=0;j<npar;j++){
		for (int i=0;i<npar;i++){
		  printf(" %.3lf", H(i+1,j+1));
		}
		printf("\n");
	      }
	      if (npar == 2) {
		double det = H(1,1) * H(2,2) - H(1,2) * H(2,1);
		printf("Determinant %.4f\n", det);
	      }
	    }

	    for (size_t l=0;l<Htie.size();++l) {
	      int i = Htie[l].index1 + 1; // +1 for fortran-like indexing
	      int j = Htie[l].index2 + 1; // +1 for fortran-like indexing
	      H(i,j) += Htie[l].Hij;
	      if (i != j) H(j,i) += Htie[l].Hij;
	    }
	  }
	}
	//^
	//  std::cout << "\n** Parameters: ";
	//  for (int i=0;i<npar;i++) {std::cout << params[i] << " ";}
	//  std::cout <<"\n";
	//  std::cout << "** Residual: " << target
	//	    << " Restraint residual " << Rtie
	//	    << " Total " << target+Rtie
	//	    << "  Nobservations " << nobs<< "\n";
	//-!

	if (DoGradient) gradientOK = true;
	target = target + Rtie;
	//  double dtime = timer.Stop();
	//^  if (DoGradient) {
	//    if (DoHessian) {
	//      std::cout << "Time for gradient + Hessian: " << dtime << " secs\n";
	//    } else {
	//      std::cout << "Time for gradient:           " << dtime << " secs\n";
	//    }
	//  } else {
	//    std::cout << "Time for target function:    " << dtime << " secs\n";
	//  }
  }
  // ---------------------------------------------------------
  void RefineScale::applyShift(TNT::Vector<floatType>& newpar)
  {
    //    std::cout << "applyShift: pars"; //^
    for (int i=0;i<npar;i++)  {
      params[i] = newpar[i];
      //      std::cout << " " << params[i]; //^
    }
    //    std::cout <<"\n"; //^
    scalemodel->SetParameters(params, nrefpar);
    gradientOK = false;
  }
  // ---------------------------------------------------------
  void RefineScale::logCurrent(outStream where, Output& output)
  {
    /*
      output.logTab(1,where,"\nScales:");
      for (int i=0;i<npar;i++) 
      {output.logTabPrintf(2,where," %7.3f", params[i]);}
      output.logTab(1,where,"\n");
    */
  }
  // ---------------------------------------------------------
  std::string  RefineScale::whatAmI(int& i)
  {
    return "Scale "+itos(i);
  }
  // ---------------------------------------------------------
  std::vector<bounds>     RefineScale::getLowerBounds()
  {
    std::vector<bounds> Lower(npar);
    double bound;
    for (int i=0;i<npar;i++) {
      if (scalemodel->GetLowerBound(i, bound)) {
	Lower[i].on(bound);
      } else {
	Lower[i].off();
      }
    }
    return Lower;
  }
  // ---------------------------------------------------------
  std::vector<bounds>     RefineScale::getUpperBounds()
  {
    std::vector<bounds > Upper(npar);
    double bound;
    for (int i=0;i<npar;i++) {
      if (scalemodel->GetUpperBound(i, bound)) {
	Upper[i].on(bound);
      } else {
	Upper[i].off();
      }
    }
    return Upper;
  }
  // ---------------------------------------------------------
  TNT::Vector<floatType> RefineScale::getRefinePars()
  {
    TNT::Vector<floatType> pars(npar);
    for (int i=0;i<npar;i++) {
      pars[i] = params[i];
    }
    return pars;
  }
  // ---------------------------------------------------------
  TNT::Vector<floatType> RefineScale::getLargeShifts()
  {
    TNT::Vector<floatType> large(npar);
    for (int i=0;i<npar;i++) {
      ///      large[i] = 1.0;
      large[i] = scalemodel->GetLargeShift(i);
    }
    return large;
  }
  // ---------------------------------------------------------
  bool1D RefineScale::getRefineMask(protocolPtr protocol)
  {
    return bool1D(npar, true);
  }
  // ---------------------------------------------------------
}
