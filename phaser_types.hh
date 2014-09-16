
#ifndef PHASER_TYPES
#define PHASER_TYPES

#include <vector>

#include <boost/smart_ptr.hpp>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/misc_functions.h>


namespace phaser { namespace af = scitbx::af; }
namespace phaser { namespace fn = scitbx::fn; }
typedef double floatType;

namespace phaser {


    typedef std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<floatType> > > > > > float6D;
    typedef std::vector<std::vector<std::vector<std::vector<std::vector<floatType> > > > > float5D;
    typedef std::vector<std::vector<std::vector<std::vector<floatType> > > > float4D;
    typedef std::vector<std::vector<std::vector<floatType> > > float3D;
    typedef std::vector<std::vector<floatType> > float2D;
    typedef std::vector<floatType> float1D;

    typedef std::complex<floatType> cmplxType;
    typedef std::vector<std::vector<std::vector<std::vector<cmplxType> > > > cmplx4D;
    typedef std::vector<std::vector<std::vector<cmplxType> > > cmplx3D;
    typedef std::vector<std::vector<cmplxType> > cmplx2D;
    typedef std::vector<cmplxType> cmplx1D;

    typedef std::vector<std::vector<std::vector<std::vector<bool> > > > bool4D;
    typedef std::vector<std::vector<std::vector<bool> > > bool3D;
    typedef std::vector<std::vector<bool> > bool2D;
    typedef std::vector<bool> bool1D;

    typedef std::vector<std::vector<std::vector<std::vector<unsigned> > > > unsigned4D;
    typedef std::vector<std::vector<std::vector<unsigned> > > unsigned3D;
    typedef std::vector<std::vector<unsigned> > unsigned2D;
    typedef std::vector<unsigned> unsigned1D;

    typedef std::vector<std::vector<std::vector<std::vector<std::string> > > > string4D;
    typedef std::vector<std::vector<std::vector<std::string> > > string3D;
    typedef std::vector<std::vector<std::string> > string2D;
    typedef std::vector<std::string> string1D;

    typedef std::vector<std::vector<std::vector<std::vector<int> > > > int4D;
    typedef std::vector<std::vector<std::vector<int> > > int3D;
    typedef std::vector<std::vector<int> > int2D;
    typedef std::vector<int> int1D;

    typedef scitbx::vec3<floatType> dvect3;
    typedef scitbx::mat3<floatType> dmat33;
    typedef scitbx::vec3<int> ivect3;
    typedef scitbx::mat3<int> imat33;

    typedef double dFloatType;
    typedef std::vector<std::vector<std::vector<std::vector<dFloatType> > > > dFloat4D;
    typedef std::vector<std::vector<std::vector<dFloatType> > > dFloat3D;
    typedef std::vector<std::vector<dFloatType> > dFloat2D;
    typedef std::vector<dFloatType> dFloat1D;

    typedef std::complex<double> dCmplxType;
    typedef std::vector<std::vector<std::vector<std::vector<dCmplxType> > > > dCmplx4D;
    typedef std::vector<std::vector<std::vector<dCmplxType> > > dCmplx3D;
    typedef std::vector<std::vector<dCmplxType> > dCmplx2D;
    typedef std::vector<dCmplxType> dCmplx1D;


}



#endif
