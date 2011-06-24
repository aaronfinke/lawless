// score_datatypes.cpp

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;

#include "csymlib.h"    // CCP4 symmetry stuff

#include "rotation.hh"
#include "score_datatypes.hh"
#include "scala_util.hh"
#include "spacegroupreindex.hh"

using phaser_io::LOGFILE;


namespace scala
{
  //--------------------------------------------------------------
  void LinearFit::add(const float& x, const float& y, const float& w)
  {
    sumw   += w;
    sumwx  += w*x;
    sumwy  += w*y;
    sumwxx += w*x*x;
    sumwxy   += w*x*y;
    np += 1;
    maxy = Max(y ,maxy);
    miny = Min(y ,miny);
    //  std::cout << "LF::add "<< x << " " << y << " " << sumwy << "\n";
  }
  //--------------------------------------------------------------
  RPair LinearFit::result() const
  {
    if (np > 0)
      {
	double d = sumw*sumwxx - sumwx*sumwx;
	double slope = (sumw*sumwxy - sumwx*sumwy)/d;
	double intercept = (sumwxx*sumwy - sumwx*sumwxy)/d;
	return RPair(slope, intercept);
      }
    else
      {return RPair(0.0,0.0);}
  }
  //--------------------------------------------------------------
  double LinearFit::slope(const float& b) const
    // Return slope only, fixed intercept b (eg = 0)
  {
    if (np > 0)
      {
	return (sumwxy - b * sumwx)/sumwxx;
      }
    else
      {return 0.0;}
  }
  //--------------------------------------------------------------
  bool operator < (const IKode& a,const IKode& b)
  // for sort on resolutiom
  {return a.sSqr < b.sSqr;}
  //--------------------------------------------------------------
  void correl_coeff::zero()
  {
    sum_wx  = 0.0;
    sum_wx2 = 0.0;
    sum_wy  = 0.0;
    sum_wy2 = 0.0;
    sum_wxy = 0.0;
    sum_w   = 0.0;

    sum_wwx = 0.0;
    sum_wwy = 0.0;
    sum_w2x = 0.0;
    sum_w2y = 0.0;
    n = 0;
    nw = 0;
  }
  //--------------------------------------------------------------
  void correl_coeff::add(const Rtype& x, const Rtype& y, const Rtype& w)
  // Add in contributions to correlation coefficient
  {
    sum_wx  += w * x;
    //^
    //^  printf("CCadd: x,sum_wx %8.4f %8.3f\n", x, sum_wx);
    sum_wx2 += w * x * x;
    sum_wy  += w * y;
    sum_wy2 += w * y * y;
    sum_wxy += w * x * y;
    sum_w   += w;
    n += 1;
  }
  //--------------------------------------------------------------
  void correl_coeff::add(const IsigI& ix, const IsigI& iy, const Rtype& w)
  // Add in contributions to correlation coefficient
  {
    Rtype x = ix.I();
    Rtype y = iy.I();
    sum_wx  += w * x;
    //^
    //^  printf("CCadd: x,sum_wx %8.4f %8.3f\n", x, sum_wx);
    sum_wx2 += w * x * x;
    sum_wy  += w * y;
    sum_wy2 += w * y * y;
    sum_wxy += w * x * y;
    sum_w   += w;
    n += 1;
  }
  //--------------------------------------------------------------
  void correl_coeff::add(const Rtype& x, const Rtype& y,
			 const Rtype& w, const Rtype& sx,const Rtype& sy)
  // Add in contributions to correlation coefficient with x,y weights
  // w is weight, sx, sy are sigma(x,y)
  {
    double wx = w/(sx*sx);
    double wy = w/(sy*sy);
    double wxy = w/(sx*sy);

    sum_wx  += wx * x ;
    sum_wx2 += wx * x * x;
    sum_wy  += wy * y;
    sum_wy2 += wy * y * y;
    sum_wxy += wxy * x * y;
    sum_w   += wxy;

    sum_wwx += wxy * x;
    sum_wwy += wxy * y;
    sum_w2x += wx;
    sum_w2y += wy;

    nw += 1;
  }

  //--------------------------------------------------------------
  void correl_coeff::dump() const
  {
    std::cout << "CCdump:\n";
    std::cout  << "sum_wx sum_wx2 sum_wy sum_wy2 sum_wxy sum_w\n"
	       << sum_wx << " " <<  sum_wx2 << " " <<  sum_wy << " "
	       <<  sum_wy2 << " " <<  sum_wxy << " " <<  sum_w << "\n";

    std::cout << "sum_wwx sum_wwy sum_w2x sum_w2y\n"
	      << sum_wwx << " " <<  sum_wwy << " " <<  sum_w2x << " "
	      <<  sum_w2y << "\n";
    std::cout << "n, nw\n"
	      << n << " " <<  nw << "\n";
  }
  //--------------------------------------------------------------
  // Number of contributions
  int correl_coeff::Number() const
  {
    //  One of n & nw should be == 0
    if (n*nw != 0)
      Message::message(Message_fatal("correl_coef.result: mixed weighting schemes") );
    return n+nw;  // either
  }
  //--------------------------------------------------------------
  void correl_coeff::SetCC()
  // set CC = 1.0 as dummy
  {
    sum_wx  = 1.0;
    sum_wx2 = 1.0;
    sum_wy  = 1.0;
    sum_wy2 = 1.0;
    sum_wxy = 1.0;
    sum_w   = 4.0;

    sum_wwx = 0.0;
    sum_wwy = 0.0;
    sum_w2x = 0.0;
    sum_w2y = 0.0;
    n = 4;
    nw = 0;    
  }
  //--------------------------------------------------------------
  ValCount correl_coeff::result() const
  // Calculate correlation coefficient from sums
  {
    // Can't mix up two sorts of "add", ie CC weighted by w &
    // CC weighted by wx, wy
    //  One of n & nw should be == 0
    if (n*nw != 0)
      Message::message(Message_fatal("correl_coef.result: mixed weighting schemes") );
    double denom = 0.0;
    if (n > 4) {
      // Same or unit weighting for x & y
      if (sum_w > 0.0) 
	denom = (sum_wx2 - (sum_wx * sum_wx)/sum_w) *
	  (sum_wy2 - (sum_wy * sum_wy)/sum_w);
      if (denom > 0.0)
	return ValCount((sum_wxy - sum_wx * sum_wy / sum_w)/sqrt(denom), n);
    }
    if (nw > 4) {
      if (sum_w2x > 0.0 && sum_w2y > 0.0) 
	denom = (sum_wx2 - (sum_wx * sum_wx)/sum_w2x) *
	  (sum_wy2 - (sum_wy * sum_wy)/sum_w2y);
      if (denom > 0.0)
	return ValCount((sum_wxy - sum_wx*sum_wwy/sum_w2x
			 - sum_wy*sum_wwx/sum_w2y
			 + sum_wx*sum_wy*sum_w/(sum_w2x*sum_w2y))
			/sqrt(denom), nw);
    }
    return ValCount(0.0,n+nw);
  }
  //--------------------------------------------------------------
  correl_coeff& correl_coeff::operator +=(const correl_coeff& other)
  {
    sum_wx  += other.sum_wx;
    sum_wx2 += other.sum_wx2;
    sum_wy  += other.sum_wy;
    sum_wy2 += other.sum_wy2;
    sum_wxy += other.sum_wxy;
    sum_w   += other.sum_w;

    sum_wwx += other.sum_wwx;
    sum_wwy += other.sum_wwy;
    sum_w2x += other.sum_w2x;
    sum_w2y += other.sum_w2y;
    nw += other.nw;
    n += other.n;

    return *this;
  }
  //--------------------------------------------------------------
  correl_coeff& operator +
  (const correl_coeff& a, const correl_coeff& b)
  {
    correl_coeff c = a;
    return c += b;
  }
  //--------------------------------------------------------------
  void MSdiff::zero()
  {
    sum_ndf2 = 0.0;
    sum_w = 0.0;
    n_f = 0;
  }
  //--------------------------------------------------------------
  // add in contribution 
  void MSdiff::add(const IsigI& Is1, const IsigI& Is2, const double& w,
		   const double& VarK)
  {
    double I1 = Is1.I();
    double sig1 = Is1.sigI();
    double I2 = Is2.I();
    double sig2 = Is2.sigI();
    sum_ndf2 += w*(I1-I2)*(I1-I2)/(sig1*sig1 + sig2*sig2 + VarK);
    n_f++;
    sum_w += w;
  }
  //--------------------------------------------------------------
  MSdiff& MSdiff::operator +=(const MSdiff& other)
  {
    sum_ndf2 += other.sum_ndf2;
    sum_w += other.sum_w;
    n_f += other.n_f;
  
    return *this;
  }
  //--------------------------------------------------------------
  MSdiff& operator +
  (const MSdiff& a, const MSdiff& b)
  {
    MSdiff c = a;
    return c += b;
  }
  //--------------------------------------------------------------
  Rfactor& Rfactor::operator +=(const Rfactor& other)
  {
    sum_df += other.sum_df;
    sum_f  += other.sum_f;
    n_f += other.n_f;
  
    return *this;
  }
  //--------------------------------------------------------------
  Rfactor& operator +
  (const Rfactor& a, const Rfactor&  b)
  {
    Rfactor c = a;
    return c += b;
  }
  //--------------------------------------------------------------
  PairSet::PairSet(const int& I1, const int& I2, const int& SymElement)
    : weight(1.0), SymElmt(SymElement)
  {
    pairs.push_back(IndexPair(I1,I2));
    fac12.push_back(RPair(1.0,1.0));
    obsindices.push_back(I1);
    obsindices.push_back(I2);
  }
  //--------------------------------------------------------------
  PairSet::PairSet(const int& I1, const int& I2, const int& SymElement,
		   const double& wt, const float& fac1, const float& fac2)
    : weight(wt), SymElmt(SymElement)
  {
    pairs.push_back(IndexPair(I1,I2));
    fac12.push_back(RPair(fac1, fac2));
    obsindices.push_back(I1);
    obsindices.push_back(I2);
  }

  //--------------------------------------------------------------
  bool PairSet::AddPair(const int& I1, const int& I2, const int& SymElement,
			const float& fac1, const float& fac2)
  {
    if (SymElement != SymElmt) return false;  // Symmetry element does not match
    int Npair = pairs.size();
    if (Npair > 0)       {
      // Check if it matches
      bool found = false;
      for (int l=0;l<Npair;l++)	{
	// Check for common observation
	if (I1 == pairs[l].first || I2 == pairs[l].first ||
	    I1 == pairs[l].second || I2 == pairs[l].second) {
	  found = true;
	  break;
	}
      }
      if (!found) return false;
    }
    // Matches or first one, add to list
    pairs.push_back(IndexPair(I1,I2));
    fac12.push_back(RPair(fac1, fac2));

    int Nobs = obsindices.size();
    if (Nobs == 0)  {
      obsindices.push_back(I1);
      obsindices.push_back(I2);
    } else  {
      bool l1 = true;
      bool l2 = true;
      for (int i=0;i<Nobs;i++)
	{
	  if (I1 == obsindices[i]) l1 = false;
	  if (I2 == obsindices[i]) l2 = false;
	}
      if (l2) obsindices.push_back(I2);
      if (l1) obsindices.push_back(I1);
    }
    return true;
  }
  //--------------------------------------------------------------
  // check all subgroups, reindex H lattices to R
  // Only use if original lattice type was R
  void ResetRhombohedralHexLatticestoR
  (std::vector<scala::PossibleSpaceGroup>& groups)
  {
    for (size_t i=0;i<groups.size();i++) {
      groups[i].ReindexRhombohedral(true);  // change H settings to R
    }
  }
  //--------------------------------------------------------------
  PossibleSpaceGroup::PossibleSpaceGroup
  (const std::string& Sname, const int& number,
   const std::string& Rname, const std::string& Scond, const std::string& ScondLG,
   const ReindexOp& Sreindex, const Chirality& Schiral,
   const double& probin, const std::vector<int>& ZonesinGroup)
    : name(Sname), refname(Rname), sgnumber(number),
      condition(Scond), conditionLG(ScondLG), reindex(Sreindex),
      chiral(Schiral), prob(probin), accepted(true)
  {
    sysabsprob = probin;
    centroprob = 1.0;
    lauegroupprob = 1.0;
    zonesingroup = ZonesinGroup;
  }
  //--------------------------------------------------------------
  void PossibleSpaceGroup::StoreAcenProb(const double& probacen)
  {
    if (chiral == CENTROSYMMETRIC)
      centroprob = 1.0-probacen;
    else
      centroprob = probacen;
    prob = lauegroupprob*sysabsprob*centroprob;
  }
  //--------------------------------------------------------------
  ReindexOp PossibleSpaceGroup::Reindex(const bool& Orig) const
  // If true, original reindex
  // else if false, total reindex operator from original to reference setting
  {
    if (Orig) {return reindex_orig;}
    // total reindex operator H =
    //     H(original->Lauegroup) * H(LaueGroup->SpaceGroup)
    return reindex_orig * reindex;
  }
  //--------------------------------------------------------------
  void PossibleSpaceGroup::StoreLaueGroupProb(const double& problg,const double& confidence)
  {
    lauegroupprob = problg;
    prob = lauegroupprob*sysabsprob*centroprob;
    lauegroupconfidence =  confidence;
  }
  //--------------------------------------------------------------
  bool PossibleSpaceGroup::IsPIForthorhombic() const
  // Return true if P,I or F orthorhombic (to allow for non-"reference" setting option
  {
    // Not clever!
    if (sgnumber >= 16 && sgnumber <= 74) {
      if ((lauegroupname[0] == 'P') || (lauegroupname[0] == 'I') || (lauegroupname[0] == 'F')) 
	{return true;}
    }
    return false;
  }
  //--------------------------------------------------------------
  // Return true if group is I2 etc (mI class)
  bool PossibleSpaceGroup::IsI2() const
  {
    return (lauegroupname ==  "I 1 2/m 1");
  }
  //--------------------------------------------------------------
  std::string PossibleSpaceGroup::Name(const bool& nameSG) const
  // basic spacegroup name if SGname true, else reference space group name 
  {
    return (nameSG ? name : refname);
  }
  //--------------------------------------------------------------
  // List of reflection conditions
  //  eg "0kl: k+l = 2n, h00: h = 2n"
  // for basic spacegroup if SGname true, else reference space group setting
  std::string PossibleSpaceGroup::Condition(const bool& nameSG) const
  {
    return (nameSG ? conditionLG : condition);
  }
  //--------------------------------------------------------------
  bool PossibleSpaceGroup::ReindexRhombohedral(const bool& OrigIsR)
  // Reindex rhombohedral H lattice to back to originalrhombohedral R setting
  //  OrigIsR  true if original lattice was R
  // Do nothing if not H lattice
  // Return true if done
  {
    if (name[0] != 'H') return false; // not H lattice
    ReindexOp htor;
    if (OrigIsR) {
      // H -> R reindex operator
      htor = reindex_orig.inverse();  // back to original
    } else {
      htor = SpacegroupReindexOp("H 3", "R 3");  // H to R general
    }

    //^
    //    std::cout << "\n\n=== ReindexRhombohedral\nH to R reindex:\n" << htor.as_hkl() << "\n";
    //    std::cout << "\nCurrent reindex solution (Laue group) to reference setting [H]:\n"
    //    	      << reindex.as_hkl() << "\n";
    //    std::cout << "\nInverse:\n"
    //    	      << reindex.inverse().as_hkl() << "\n";
    //    std::cout << "\nCurrent reindex operator input to Laue group [H]:\n"
    //    	      << reindex_orig.as_hkl() << "\n";
    //    std::cout << "\nInverse:\n"
    //    	      << reindex_orig.inverse().as_hkl() << "\n";
    //    std::cout << "[HtoR] * [H]\n" << (htor * reindex_orig).as_hkl() << "\n";
    //^-
    // Change H to R in most names
    name = SGnameHtoR(name, 'R');         // Spacegroup name
    //refname;      // Reference spacegroup name stays as "H"
    lauegroupname  = SGnameHtoR(lauegroupname, 'R');  // Lauegroup name
    pointgroupname = SGnameHtoR(pointgroupname, 'R');  // Pointgroup name
    // condition & conditionLG should be ""
    // Rhombohedral space group
    SpaceGroup RSG(pointgroupname);
    // reindex operator from R to reference H setting
    ///    reindex = htor * reindex;
      ///    reindex.FindSimplest(RSG);
    // reindex operator original setting to Laue group
    reindex_orig = reindex_orig * htor;
    reindex_orig.FindSimplest(RSG);
    //^
    //    std::cout <<
    //      "ReindexRhombohedral: reindex operator original setting to Laue group "
    //    	      << reindex_orig.as_hkl() << "\n===\n\n"; 
    //^-
    return true;
  }
  //--------------------------------------------------------------
} // end namespace scala

// end scala_datatypes.
