// datatypes.cpp

#include <math.h>
#include <cmath>

#include "sdctypes.hh"
#include "hkl_unmerge.hh"

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;


namespace scala {
  // minimum fraction of input variance for correction
  const double SDcorrection::MINVARINFRAC = 0.1;
  const double SDcorrection::MINSDFAC = 0.2;     // minimum SDfac
  const int SDcorrection::NPARALL;               // number of parameters = 3

  //--------------------------------------------------------------
  SDcorrection::SDcorrection() : sdfac(1.0), sdb(0.0), sdadd(0.0) {ResetRange();}
  // Reset minimum & maximum
  //--------------------------------------------------------------
  SDcorrection::SDcorrection(const double& SDfac,
			     const double& SDb, const double& SDadd,
			     const bool& fixSDb)
    : sdfac(SDfac), sdb(SDb), sdadd(SDadd), fixsdb(fixSDb),
      weights(std::vector<double>(3,0.0))
  {
    sdadd2 = sdadd * sdadd;
    SetRefineParameters();
    ClearRestraints();
    ResetRange();
  }
  //--------------------------------------------------------------
  void SDcorrection::Set(const double& SDfac, const double& SDb, const double& SDadd)
  {
    sdfac = SDfac; sdb = SDb;
    sdadd = SDadd;
    sdadd2 = sdadd*sdadd;
    if (sdadd < 0.0) sdadd2 = -sdadd2; // store SdAdd^2, negated if necessary
  }
  //--------------------------------------------------------------
  double SDcorrection::SigmaPrime(const double& sigma, const double& gscale,
				  const double& Iav) const
  // Correct sigma on Ihl scale
  // sigma  uncorrected sigma(Ihl)
  // gscale inverse scale for Ihl
  // Iav    average <Ih>
  // 
  // Ih = gscale * Iav
  // sigma' = SdFac * sqrt ( sigma^2 + SdB * Ih + (SdAdd * Iav)^2)
  //        = SdFac * sqrt ( sigma^2 + SdB * Ih + SdAdd2 * Iav^2)
  {
    double var = sigma*sigma;
    double gi = gscale * Iav;
    double si2 = sdadd2 * gi * gi; // sdadd2 may be negative
    double ssc = Max(MINVARINFRAC*var, (var + sdb*gscale*gi + si2));
    return Max(sdfac, MINSDFAC) * sqrt(ssc);
  }
  //--------------------------------------------------------------
  float SDcorrection::Correct(observation& Observation, const float& Iav) const
  // in-place correction, using Iav as intensity
  // returns original uncorrected but scaled sd(I)
  {
    float sig0 = Observation.ksigI();  // scaled sigI
    // sd' on Ihl scale 
    double sd = SigmaPrime(Observation.sigI(), Observation.Gscale(), Iav);
    double corr = sd/Observation.sigI();
    mincorr = Min(corr, mincorr);
    maxcorr = Max(corr, maxcorr);
    Rtype fsd = sd;
    //^
    //    if (fsd <= 0.0) {
    //      std::cout << "sdfac " << sdfac << " " << fsd <<"\n";;
    //    }  //^
    Observation.set_sigI(fsd);
    return sig0;
  }
  //--------------------------------------------------------------
  IsigI SDcorrection::Correct(const IsigI& Is, const double& gscale,
				const float& Iav) const
  {
    double sd = SigmaPrime(Is.sigI(), gscale, Iav);
    double corr = sd/Is.sigI();
    mincorr = Min(corr, mincorr);
    maxcorr = Max(corr, maxcorr);
    return IsigI(Is.I(), sd);
  }
  //--------------------------------------------------------------
  std::string SDcorrection::format() const
  {
    return
      "SdFac = "+clipper::String(Max(sdfac, MINSDFAC),6,3)+
      ",  SdB = "+clipper::String(sdb,6,3)+
      ",  SdAdd = "+clipper::String(sdadd,6,3);
  }
  //--------------------------------------------------------------
  std::vector<double> SDcorrection::GetParameters() const
  {
    std::vector<double> params;
    double sdf2 = sdfac*sdfac;
    params.push_back(sdf2);
    if (!fixsdb) params.push_back(sdf2*sdb);
    params.push_back(sdf2*sdadd2);
    p = params; // save
    return params;
  }
  //--------------------------------------------------------------
  std::vector<double> SDcorrection::GetRealParameters() const
  //  return SdFac, [SDb,] SDadd
  {
    std::vector<double> realparams;
    realparams.push_back(sdfac);
    if (!fixsdb) realparams.push_back(sdb);
    realparams.push_back(sdadd);
    return realparams;
  }
  //--------------------------------------------------------------
  void SDcorrection::SetRefineParameters()
  // Set internal vector p (2 or 3) from Sdfac etc
  //  vector p is (p1,p2,p3), p2 may be missing 
  {
    p.clear();
    double sdf2 = sdfac*sdfac;
    p.push_back(sdf2);
    if (!fixsdb) p.push_back(sdf2*sdb);
    double s2 = sdadd*sdadd;
    if (sdadd < 0.0) s2 = -s2;
    p.push_back(sdf2*s2);
    return;
  }
  //--------------------------------------------------------------
  void SDcorrection::SetParameters(const std::vector<double>& params)
  // Set all parameters from vector (2 or 3)
  //  vector params is (p1,p2,p3), p2 may be missing 
  {
    p = params;  // save
    if (params[0] > 0.0) {
      sdfac = sqrt(params[0]);
    } else {
      sdfac = MINSDFAC;
      p[0] = sdfac*sdfac;
    }
    sdfac = Max(sdfac, MINSDFAC);
    int k=1;
    if (!fixsdb) sdb = params[k++]/(sdfac*sdfac);	 
    sdadd2 = params[k]/(sdfac*sdfac);
    sdadd = sqrt(std::abs(sdadd2));
    if (sdadd2 < 0.0) sdadd = -sdadd;
  }
  //--------------------------------------------------------------
  std::vector<double> SDcorrection::GetShifts(const double& scale) const
  // Initial shifts for each parameter type, scaled by "scale"
  {
    const double SDFAC_SHIFT = 0.3;
    const double SDB_SHIFT =  2;
    const double SDADD_SHIFT = 0.01;
    std::vector<double> shifts;
    shifts.push_back(scale * SDFAC_SHIFT);
    if (!fixsdb) shifts.push_back(scale * SDB_SHIFT);
    shifts.push_back(scale * SDADD_SHIFT);
    return shifts;
  }
  //--------------------------------------------------------------
  std::vector<double> SDcorrection::GetDerivatives(const float& sigma,
						   const float& Iav) const
  // On entry:
  //  sigma    uncorrected sigma, scaled to Iav
  //  Iav      average I <Ih>
  //
  // Return d(sigma')/dp vector
  //  parameters are p1,p2,p3
  //   s'^2 = p1 s^2 + p2 Ih + p3 Ih^2
  //   ds'^2/dp1 = s^2
  //   ds'^2/dp2 = Ih
  //   ds'^2/dp3 = Ih^2
  //   ds'/dp = ds'/d(s'^2) * ds(s'^2)/dp
  //   ds'/d(s'^2) = 1/2s'
  //
  // Uses vector p = (p1,[p2,]p3) not sdfac etc
  {
    std::vector<double> dp;
    double var = sigma*sigma;
    double Ih = Iav;
    // s'^2 = p1 s^2 + p2 Ih + p3 Ih^2
    double sd2 = p[0] * var;
    int k=1;
    if (!fixsdb) {
      sd2 += p[k++] * Ih;
    }
    sd2 += p[k] * Ih * Ih;
    //    if (sd2 < MINVARINFRAC*var) {
    //      // corrected sd^2 too small
    //      sd2 = MINVARINFRAC*var;
    //    }
    double sd = sqrt(std::abs(sd2));
    if (sd2 < 0.0) sd = -sd;
    //  ds'/d(s'^2) = 1/2s'
    double dsdsp = 1.0/(2.*sd);
    // ds'^2/dp1 = s^2
    dp.push_back(dsdsp * var);
    if (!fixsdb) {
      //  ds'^2/dp2 = Iav
      dp.push_back(dsdsp * Ih);
    }
    //  ds'^2/dp3 =Iav^2
    dp.push_back(dsdsp * Ih * Ih);
    return dp;
  }
  //--------------------------------------------------------------
  std::vector<double> SDcorrection::LowerBounds() const
  //! Get vector of lower bounds (0.0 means no bound)
  {
    std::vector<double> bounds;
    double b1 = 0.1; // Sdfac
    bounds.push_back(b1*b1); // p1
    if (!fixsdb) bounds.push_back(0.0);  // no lowerbound on SdB
    double b3 = 0.1;  // SdAdd
    bounds.push_back(-b1*b3*b3);  // p3
    return bounds;
  }
  //--------------------------------------------------------------
  std::vector<double> SDcorrection::UpperBounds() const
  //! Get vector of upper bounds (0.0 means no bound)
  {
    std::vector<double> bounds;
    double b1 = 10.0; // Sdfac
    bounds.push_back(b1*b1); // p1
    if (!fixsdb) bounds.push_back(0.0); // no upper bound on SdB
    double b3 = 0.1;  // SdAdd
    bounds.push_back(b1*b3*b3);  // p3
    return bounds;
  }
  //--------------------------------------------------------------
  std::vector<double> SDcorrection::LargeShifts() const
  //! Get vector of "large shifts"
  {
    std::vector<double> large;
    double l1 = 0.5; // Sdfac
    large.push_back(l1*l1); // p1
    if (!fixsdb) {
      double l2 = 1.0; // SdB
      large.push_back(l2);
    }
    double l3 = 0.02;
    large.push_back(l3*l3);  // for SDadd^2
    return large;
  }
  //--------------------------------------------------------------
    //! clear all restraints
  void SDcorrection::ClearRestraints()
  {
    sdtargets.assign(NPARALL,0.0);  // always 3
    weights.assign(NPARALL,0.0);  // always 3
  }
  //--------------------------------------------------------------
  //! set targets & SD (=0 for no restraint), for 3 parameters always
  //! weight = 1/SD^2
  void SDcorrection::SetRestraints(const std::vector<double>& Targets,
				   const std::vector<double>& SDtarget)
  // If FixSdB (ie 2 parameters), this can still specify 3 target/weight
  // pairs & the 2nd one will be ignored
  // Target for SdAdd is given here for SdAdd, but stored as target for SdAdd^2 = sdadd2
  {
    bool OK = true;
    ASSERT (Targets.size() == SDtarget.size());
    ASSERT (int(Targets.size()) == NPARALL);
    rawtargets = Targets;
    sdtargets = SDtarget;
    targets = rawtargets;
    weights = sdtargets;  // in case they are = 0.0
    
    int k = 0;
    // SdFac
    // Check for sensible values
    if (sdtargets[k] > 0.0 && 
	(targets[k] < 0.01 || targets[k] > 100.0)) {
      OK = false;
    }
    k++;
    if (!fixsdb) { // SdB present
      // SdB
      if (sdtargets[k] > 0.0 && 
	  (targets[k] < -20 || targets[k] > 20.0)) {
	OK = false;
      }
    }
    // SdAdd
    k++;
    if (sdtargets[k] > 0.0 && 
	(targets[k] < -0.2 || targets[k] > 0.2)) {
      OK = false;
    }
    for (k=0;k<int(sdtargets.size());++k) {
      if (sdtargets[k] < 0.0) {OK = false;}
    }
    if (!OK) {
      std::string s = "SDcorrection::SetRestraints unreasonable values:";
      for (k=0;k<int(sdtargets.size());++k) {
	s += " "+clipper::String(targets[k])+" "+clipper::String(sdtargets[k]);
      }
      Message::message(Message_fatal(s));
    }
    // Store weights = 1/SD^2
    for (k=0;k<int(sdtargets.size());++k) {
      if (sdtargets[k] > 0.0) {
	weights[k] = sdtargets[k]*sdtargets[k];
	if (k==int(sdtargets.size())-1) { // SdAdd
	  // for SdAdd use target to SdAdd^2 = sdadd2
	  targets[k] = targets[k]*targets[k];
	  // and square SD as well (not necessarily correct!)
	  weights[k] = sdtargets[k]*sdtargets[k];
	}
	weights[k] = 1.0/(weights[k]*weights[k]); // weight = 1/SD^2
      }
    }
   }
  //--------------------------------------------------------------
  //! set default values for restraints (on SdB only)
  void SDcorrection::SetDefaultRestraints()
  {
    // Always 3 values, even if SdB [1] will be ignored
    double dtarget[] = {0.0,0.0,0.0};
    double dweight[] = {0.0,10.0,0.0};  // SdB only
    SetRestraints(std::vector<double>(dtarget, dtarget+3),
		  std::vector<double>(dweight, dweight+3));
  }
  //--------------------------------------------------------------
  //! return restraint values
  void SDcorrection::GetRestraints(std::vector<double>& Targets,
				   std::vector<double>& SDtarget) const
  {
    Targets = rawtargets;
    SDtarget = sdtargets;
  }
  //--------------------------------------------------------------
  //! return contribution to restraint residual R2, summed over parameters
  //  R2 = 1/2 Sum(parameters) (weight (value - target)^2)
  double SDcorrection::RestraintR() const
  {
    double R2 = 0.0;
    int k=0;
    if (weights[k] > 0.0) { // SdFac
      R2 += weights[k] * (sdfac - targets[k]) * (sdfac - targets[k]);
    }
    k++;
    if (!fixsdb) { // SdB
      if (weights[k] > 0.0) { // SdFac
	R2 += weights[k] * (sdb - targets[k]) * (sdb - targets[k]);
      }}
    k++;
    if (weights[k] > 0.0) { // SdAdd
	R2 += weights[k] * (sdadd2 - targets[k]) * (sdadd2 - targets[k]);
    }
    return 0.5*R2;
  }
  //--------------------------------------------------------------
  //! return derivative vector dR2/dp and Hessian, 2 or 3 parameters
   void SDcorrection::RestraintDerivatives
   (std::vector<double>& dp, clipper::Array2d<double>& H) const
  {
    // this is not elegant!
    int npar = Nparams();

    // dq(l)/dp(j)  partial derivatives for l, j = 0,npar-1
    //  only for parameters used (ie maybe not SdB)
    clipper::Array2d<double> dqdp(npar,npar,0.0);
    // NB l = 0,npar-1, k = 0,NPARALL-1, ie different if fixsdb

    int l = 0;
    int k = 0; // q[k]    
    int j = 0; // p[j]
    // --- p1 components, wk * Delk * dqk/dpj
    //  q1 (SdFac),p1  dq1/dp1 = 1/(2 q1)
    if (weights[k] > 0.0) {
      dqdp(l,j) = 1.0 / (2.0*sdfac);
    }
    k++;l++;
    if (!fixsdb) {
      if (weights[k] > 0.0) {
	// q2 (SdB), p1 , dq2/dp1 = -p2/p1
	dqdp(l,j) =  -p[1] / p[0];
      }
      l++;
    }
    k++;
    // q3 (SdAdd^2), p1, dq3^2/dp1 = -p3/p1^2
    if (weights[k] > 0.0) {
      dqdp(l,j) =  -p[2]/(p[0]*p[0]);
    }  // 
    // --- p2 components
    l = 1;  // for p2, only q2, dq1/dp2 = 0 = dq3^2/dp2
    k = 1;
    int l3 = 1;  // index for SdAdd (q3)
    if (!fixsdb) {
      j++;
      //  q2 (SdB),p2  dq1/dp2 = 1/p1
      if (weights[k] > 0.0) {
	dqdp(l,j) = 1.0 / p[0];
      }
      l++;      
      l3 = 2;  // if SdB is used, SdAdd is 3rd parameter
    }
    k++;
    // --- p3 components
    j++;
    // k = 1 or 2, only q3, dq1/dp3 = 0 = dq2/dp3
    if (weights[k] > 0.0) {
      // q3 (SdAdd), p3 , dq3^2/dp3 = 1/p1
      dqdp(l,j) = 1.0 / p[0];
    }
    ASSERT (j+1 == Nparams());

    dp.assign(npar, 0.0);
    H.resize(npar,npar,0.0);  // Hessian ~= d2R/dp2
    double wD;  // weight * delta

    for (j=0;j<npar;++j) { // loop parameter j
      l = 0;
      for (k=0;k<NPARALL;++k) { // loop parameter k
	if (k != 1 || !fixsdb) {
	  if (k == 0) { 
	    wD = weights[k] * (sdfac - targets[k]);
	  } else if (k == 1) {
	    wD = weights[k] * (sdb - targets[k]);
	  } else if (k == 2) {
	    wD = weights[k] * (sdadd2 - targets[k]);
	  }
	  // Gradient = - w dR/dp  (NB negative!)
	  dp[j] -= wD * dqdp(l,j);
	  // Hessian
	  for (int i=0;i<npar;++i) { // loop parameter i
	    H(i,j) += weights[k] * dqdp(l,i) * dqdp(l,j);
	  }
	  l++;
	} // not SdB
      } // end loop k
    } // end loop j 
    //^
    //    std::cout << "Gradient: ";
    //    for (j=0;j<npar;++j) {std::cout << " " << dp[j];}
    //    std::cout <<"\ndq/dp: ";
    //    for (j=0;j<npar;++j) {
    //      for (k=0;k<npar;++k) {std::cout << " " << dqdp(k,j);}
    //    }
    //    std::cout <<"\nHessian:\n";
    //    for (j=0;j<npar;++j) {
    //      for (int i=0;i<npar;++i) {std::cout << " " << H(i,j);}
    //      std::cout <<"\n";
    //    }
    //^-
    // All done
  }
  //--------------------------------------------------------------
  std::string SDcorrection::FormatSave() const
  {
    std::string ds = "SDcorrection V1 { ";
    if (fixsdb) {
      ds += "FIXSDB  ";
    } else {
      ds += "FREESDB ";
    }
    ds += clipper::String(sdfac)+" "+clipper::String(sdb)+" "+clipper::String(sdadd);
    ds += " }";
    return ds;
  }
  //--------------------------------------------------------------
  void SDcorrection::Restore(Fileread& FR)
  {
    FR.ReadTag("SDcorrection"); // fails if tag does not match
    if (FR.GetTag() != "V1") {  // version check
      clipper::Message::message(Message_fatal
	("SDcorrection::Restore incompatible version"));
    }
    FR.Skip();
    std::string fsdb = FR.GetTag();
    if (fsdb == "FREESDB") {
      fixsdb = false;
    } else if (fsdb == "FIXSDB") {
      fixsdb = true;
    }
    sdfac = FR.Double();
    sdb = FR.Double();
    sdadd = FR.Double();
    if (!FR.CheckEnd()) {
      clipper::Message::message(Message_warn
	("SDcorrection Restore unexpected tag "+FR.Tag()));
    }

    sdadd2 = sdadd*sdadd;
    if (sdadd < 0.0) sdadd2 = -sdadd2; // store SdAdd^2, negated if necessary
    ResetRange();
  }
  //--------------------------------------------------------------
}
