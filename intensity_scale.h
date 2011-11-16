/* sfscale.cpp: structure factor anisotropic scaling implementation */
//C Copyright (C) 2000-2006 Kevin Cowtan and University of York
//L
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA

#ifndef CLIPPER_ISCALE
#define CLIPPER_ISCALE


/*#include "../core/hkl_datatypes.h"
#include "../core/xmap.h"
#include "../core/nxmap_operator.h"
#include "clipper/contrib/sfscale.h"
#include "../core/resol_targetfn.h"*/

#include "clipper/core/hkl_datatypes.h"
#include "clipper/core/xmap.h"
#include "clipper/core/nxmap_operator.h"
#include "clipper/contrib/sfscale.h"
#include "clipper/core/resol_targetfn.h"

#include "intensity_target.h"


namespace clipper {

  // Base class for intensity scaling methods (NDS)
  template<class T> class Iscale_base {
  public:
    virtual bool operator() ( HKL_data<datatypes::I_sigI<T> >& io ) = 0;
    virtual ~Iscale_base() {}  //!< destructor
  };

  // Intensity anisotropic scaling (NDS)

  template<class T> class Iscale_aniso : public Iscale_base<T> {
  public:
    // constructor: takes rejection criterion for I/sigI
    Iscale_aniso( double nsig = 0.0 ) : nsig_(nsig) {}
    // Scale Io to isotropic (approximate)
    bool operator() ( HKL_data<datatypes::I_sigI<T> >& io );
    const U_aniso_orth& u_aniso_orth() const { return u; }
  private:
    U_aniso_orth u;
    double nsig_;
  };


template<class T> bool Iscale_aniso<T>::operator() ( HKL_data<datatypes::I_sigI<T> >& io )  //(NDS)
{
  typedef HKL_info::HKL_reference_index HRI;
  // expand to P1 in order to preserve symmetry
  const HKL_info& hkls = io.hkl_info();
  Spacegroup spgrp1( Spacegroup::P1 );
  HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
  HKL_data<datatypes::I_sigI<T> > io1( hkl1 ), is1( hkl1 ), ic1( hkl1 );
  for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
    datatypes::I_sigI<T> i = io[ih.hkl()];
    if ( i.I() >= nsig_ * i.sigI() ) is1[ih] = io1[ih] = i;
  }

  // perform aniso scaling 3 times to allow aniso scale from previous
  // cycle to correct for missing data on the next cycle:
  // start with unscaled data and iterate
  BasisFn_log_aniso_gaussian bfn;
  std::vector<ftype> param( 7, 0.0 ), params( 12, 1.0 );
  for ( int c = 0; c < 3; c++ ) {
    // create artificial I's from mean I with resolution
    TargetFn_meanInth<datatypes::I_sigI<T> > tfns( is1, 1.0 );  //check powers (NDS)  was 2.0
    BasisFn_spline bfns( is1, 12 );
    ResolutionFn rfns( hkl1, bfns, tfns, params );
    for ( HRI ih = hkl1.first(); !ih.last(); ih.next() )
      //ic1[ih] = datatypes::I_sigI<T>( sqrt(rfns.f(ih)), 1.0 );  
	  ic1[ih] = datatypes::I_sigI<T>( rfns.f(ih), 1.0 );  


    // do the aniso scaling
    TargetFn_scaleLogI1I2<datatypes::I_sigI<T>,datatypes::I_sigI<T> >
      tfn( io1, ic1 );     // check this exists (NDS)
    ResolutionFn rfn( hkl1, bfn, tfn, param );
    param = rfn.params();

    // set trace to zero (i.e. no isotropic correction)
    ftype dp = (param[1]+param[2]+param[3])/3.0;
    param[1] -= dp;
    param[2] -= dp;
    param[3] -= dp;
    // create new scaled list correct with current estimate
    is1 = io1;
    for ( HRI ih = hkl1.first(); !ih.last(); ih.next() )
      if ( !is1[ih].missing() )
	is1[ih].scale( exp( 0.5*bfn.f(ih.hkl(),hkl1.cell(),param) ) );  //was 0.5 NDS
    u = bfn.u_aniso_orth( param );
    //std::cout << c << " | " << param[1] << " " << param[2] << " " << param[3] << " " << param[4] << " " << param[5] << " " << param[6] << "\n";
  }

  // store the results
  for ( HRI ih = hkls.first(); !ih.last(); ih.next() )
    if ( !io[ih].missing() )
      io[ih].scale( exp( 0.5*bfn.f(ih.hkl(),hkls.cell(),param) ) ); //was 0.5 NDS
  return true;
}


} // namespace clipper

#endif
