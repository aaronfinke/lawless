# Dependencies:
#  CCTBX
#  TNT

DIR = ${PWD}

CC = cc
CXX=c++
CFLAGS= -O2
LFLAGS= -O2
#CFLAGS= -g
#LFLAGS= -g

# Point to source-build CCP4
CCP4=/usr/local/CCP4versions/Current
CLIB=${CCP4}/lib
CLIBS=${CLIB}/libccp4

CCFLAGS = -I${CLIBS}

XCPPFLAGS= 

# CCTBX locations
CCTBX=${CLIB}/cctbx
CCTBX_sources=${CCTBX}/cctbx_sources

# Includes

#*# version before 6.1.24
#*#ITBX=\
#*#-I${CCTBX_sources}/boost \
#*#-I${CCP4}/include \
#*#-I${CCTBX_sources}/cctbx/include
#*#-I${CCTBX_sources}/scitbx/include

# 6.4.0
ITBX=\
-I${CCP4}/include/cctbx -I${CCP4}/include

# 6.3.0
#ITBX=\
#-I${CCTBX_sources}/boost \
#-I${CCP4}/include \
#-I${CCTBX_sources}/cctbx_project

ILIB = 

# TNT
###TNT = /Users/pre/Projects/Xtal/src/TNT
TNT = ${CCTBX_sources}/tntbx/include
ITNT = -I${TNT}

INCLUDE=-I.  ${ILIB} ${ITBX} ${ITNT} 

# Libraries
LTBX= -L${CLIB} -lcctbx 

LDLIBS=${LTBX} 

CPPFLAGS=${XCPPFLAGS}  ${INCLUDE}

MINLIB=MinLib

# Objects
LIBOBJ= RefineBase.o RefineBase2.o Minimizer.o PhaserError.o Errors.o

SOURCE= RefineBase.cc RefineBase2.cc Minimizer.cc PhaserError.cc Errors.cc

HDR= Errors.h        PhaserError.h   RefineBase2.h   maths.h         symm_eigen.h \
Minimizer.h     ProtocolBase.h  fmat.h          phaser_types.h  vec.h \
Output.h        RefineBase.h    jiffy.h         subscript.h

lib: libbfgs
libbfgs.a: libbfgs

libbfgs: ${LIBOBJ}
	ar -r libbfgs.a  ${LIBOBJ}

.SUFFIXES: .cpp .o

.cpp.o: ${HDR}
	${CXX} -c $(CFLAGS) $(CPPFLAGS) $<

.cc.o: ${HDR}
	${CXX} -c $(CFLAGS) $(CPPFLAGS) $<

.c.o:
	${CC} -c $(CFLAGS) $(CCFLAGS) $<

clean: 
	rm *.o ${EXE} libbfgs.a

####################################################################
# Export etc
EXPORTDIR = ../../Export

export:
	/bin/rm -rf ${EXPORTDIR}/${MINLIB}
	mkdir ${EXPORTDIR}/${MINLIB}
	cp ${SOURCE} ${EXPORTDIR}/${MINLIB}/
	cp ${HDR} ${EXPORTDIR}/${MINLIB}/
	cp Makefile ${EXPORTDIR}/${MINLIB}/
