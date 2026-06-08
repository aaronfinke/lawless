#ifndef __PHASER_MATH__
#define __PHASER_MATH__
//#include <phaser/include/FloatType.h>
//#include <phaser/include/Phaser.h>

#include "phaser_types.h"
#include <cmath>


namespace phaser {

  floatType Fn_alogI0(floatType);
  floatType eBesselI0(floatType);
  floatType eBesselI1(floatType);
  floatType besrat(floatType&);
  floatType sim(floatType&);
  floatType pRiceF(floatType,floatType,floatType);
  floatType pRiceI(floatType,floatType,floatType);
  floatType pWoolfsonF(floatType,floatType,floatType);
  floatType pWoolfsonI(floatType,floatType,floatType);



// The getalogch function is now encapsulated in the auxiliary class, Alogch, instantiated
// by DataMR for the sake of thread safety
class Alogch
{
   int tbllen;
   int index;
   floatType tableArg;
   floatType argcut1,argcut2;
   floatType argfac1,argfac2;
   floatType log2;
   std::vector<std::pair<floatType,floatType> > lnchTable1;
   std::vector<std::pair<floatType,floatType> > lnchTable2;
   floatType absArg;
   std::vector<std::pair<floatType,floatType> > getlnchTable(int tbllen, floatType argmax);

public:
   Alogch();
   floatType getalogch(floatType&);
};


// The getalogI0 function is now encapsulated in the auxiliary class, AlogchI0, instantiated
// by DataMR for the sake of thread safety
class AlogchI0
{
  // Finely-sampled linear interpolation trades a small amount
  // of memory for significant speed improvement, compared with
  // quadratic interpolation
   int tbllen1, tbllen2, tbllen3, tbllen4;
   floatType ONE;
   floatType argcut1, argcut2, argcut3;
   floatType argfac1, argfac2;
   floatType argfac3;
   floatType invArgFac;
   std::vector<std::pair<floatType,floatType> > lnI0Table1;
   std::vector<std::pair<floatType,floatType> > lnI0Table2;
   std::vector<std::pair<floatType,floatType> > lnI0Table3;
   std::vector<std::pair<floatType,floatType> > lnInvI0Table;

   std::vector<std::pair<floatType,floatType> > getlnInvI0Table(int tbllen, floatType argcut);
   std::vector<std::pair<floatType,floatType> > getlnI0Table(int tbllen, floatType argmax);

public:
   AlogchI0();
   floatType getalogI0(floatType& arg);
};



}//end namespace phaser

#endif
