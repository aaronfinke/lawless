//(c) 2000-2007 Cambridge University Technical Services Ltd
//All rights reserved
#ifndef __RefineBaseClass__
#define __RefineBaseClass__
#include "ProtocolBase.h"
#include "Output.h"
#include "Errors.h"
#include "fmat.h"
#include "jiffy.h"
using namespace phaser_io;

namespace phaser {

class bounds
{
  public:
    bounds() {bounded = false; limit = 0;}
    bounds(bool b,double f) {bounded = b; limit = f;}
    void off() {bounded = false; limit = 0;}
    void on(double f) {bounded = true; limit = f;}
  public:
    bool bounded;
    double limit;
};
 
class reparams
{
  public:
    reparams() {reparamed = false; offset = 0;}
    reparams(bool b,double f) {reparamed = b; offset = f;}
    void off() {reparamed = false; offset = 0;}
    void on(double f) {reparamed = true; offset = f;}
  public:
    bool reparamed;
    double offset;
};

typedef std::vector<reparams> repars1D;
typedef std::vector<bounds> bounds1D;

//abstract - reminder! not objects of an abstract
//base class can be instantiated
class RefineBase  //abstract 
{
  public:
    RefineBase() { npars_ref = npars_all = 0; refinePar.clear(); }
    virtual ~RefineBase() {}
  
    //pure virtual functions
    virtual double      targetFn() = 0;
    virtual double      gradientFn(TNT::Vector<double>&) = 0;
    virtual double      hessianFn(TNT::Fortran_Matrix<double>&,bool&) = 0;
    virtual void        applyShift(TNT::Vector<double>&) = 0;
    virtual bool1D      getRefineMask(protocolPtr) = 0;
    virtual std::string whatAmI(int&) = 0;
    virtual void        logCurrent(outStream,Output&) = 0;
    virtual TNT::Vector<double>  getRefinePars() = 0;
    virtual TNT::Vector<double>  getLargeShifts() = 0;

    //virtual functions
    virtual void        cleanUp(outStream,Output&);
    virtual void        rejectOutliers(outStream,Output&);
    virtual repars1D    getRepar();
    virtual double      getMaxDistSpecial(TNT::Vector<double>&, TNT::Vector<double>&, TNT::Vector<double>&);
    virtual void        logProtocolPars(outStream,Output&);
    virtual bounds1D    getLowerBounds();
    virtual bounds1D    getUpperBounds();
    virtual void        setProtocol(protocolPtr,Output&);

    //functions
    int    numRefinePars();
    double getMaxDist(TNT::Vector<double>&, TNT::Vector<double>&, TNT::Vector<double>&);
    double getMaxStep(TNT::Vector<double>&, TNT::Vector<double>&);
    int    filterGradient(TNT::Vector<double>&,TNT::Vector<double>&,double);
    bool   fixHessian(TNT::Fortran_Matrix<double>&);
    void   studyParams(Output&);
    void   logVector(outStream,std::string,Output&,TNT::Vector<double>&);
    void   logHessian(outStream,Output&,TNT::Fortran_Matrix<double>&);
    void   reparBounds(bounds1D&);
    TNT::Vector<double> reparRefinePars(TNT::Vector<double>);
    TNT::Vector<double> reparRefineParsInv(TNT::Vector<double>);
    TNT::Vector<double> reparGradient(TNT::Vector<double>&,TNT::Vector<double>);
    void   reparHessian(TNT::Vector<double>&,TNT::Vector<double>&,TNT::Fortran_Matrix<double>&);
    void   reparLargeShifts(TNT::Vector<double>&);
 
  protected:
    int    npars_ref;
    int    npars_all;
    bool1D refinePar;
};

} //phaser

#endif

