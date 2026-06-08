//(c) 2000-2007 University of Cambridge
//All rights reserved
#ifndef __RefineClass__
#define __RefineClass__
#include "phaser_types.h"
#include "ProtocolBase.h"
#include "RefineBase.h"
#include "Output.h"

namespace phaser {

class Minimizer
{
  public:
    Minimizer() {}

  public:
    void run(RefineBase&,protocolPtr,Output&);
    void run(RefineBase&,af::shared<protocolPtr>,Output&,bool=true,bool=false);
    //python wrapper
    //    void run(RefineSAD&,ProtocolSAD,Output&);

  private:
    floatType  descent(RefineBase&,Output&,floatType,bool&,bool&);
    floatType  bfgs(RefineBase&,Output&,floatType,bool&,bool&,bool&,bool&,TNT::Vector<floatType>&,
                    TNT::Vector<floatType>&,TNT::Fortran_Matrix<floatType>&);
    floatType  newton(RefineBase&,Output&,floatType,bool&,bool&);
    floatType  LineSearch(RefineBase&,Output&,TNT::Vector<floatType>,floatType,TNT::Vector<floatType>,
                          TNT::Vector<floatType>&,floatType,floatType,floatType,bool&,bool&);

    floatType  ShiftScore(floatType, RefineBase&, TNT::Vector<floatType>&, TNT::Vector<floatType>&,
                           TNT::Vector<floatType>&, Output* poutput = NULL, bool bcount = true);
    void ShiftX(floatType, RefineBase&, TNT::Vector<floatType>&, TNT::Vector<floatType>&,
                        TNT::Vector<floatType>&);

    bool GetQuadraticCoefficients(floatType, floatType, floatType, floatType,
                                  floatType, floatType, floatType&, floatType&, floatType&);
    bool GetQuadraticMinimum(floatType&, floatType, floatType, floatType, floatType, floatType,
                             floatType, floatType, floatType);

    int mincount;
};


} //phaser
#endif
