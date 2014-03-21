
DIR = ${PWD}

CC = cc
CXX=c++
LXX=c++

OPT = -g -Wreturn-type
#OPT = -g -Wall
#OPT = -O2 -Wall
#OPT = -O2 -Wreturn-type
#STATIC =
STATIC = -static

CFLAGS= ${OPT}
LFLAGS= ${OPT} ${STATIC}


#  CCP4 locations
#CCP4=${CCP4}  # default
CCP4SRC=/public/xtal/CCP4/ccp4-64/CurrentSource
#CCP4SRC=/public/xtal/CCP4/ccp4-64/ccp4-6.4.0/ccp4-src-6.4.0/ccp4-6.4.0
#CCP4SRC=/public/xtal/CCP4/ccp4-64/ccp4-6.4.0/ccp4-6.4.0-gfortran
#CCP4SRC=/public/xtal/CCP4/ccp4-64/ccp4-6.4.0/ccp4-6.4.0-jhbuild/ccp4-src-6.4.0/ccp4-6.4.0
CCP4=${CCP4SRC}
CLIB=${CCP4SRC}/lib
CLIBS=${CLIB}/libccp4
Clipper=${CLIB}/clipper

CCFLAGS = -I${CLIBS}

XCPPFLAGS= -falign-loops=16

#XCPPFLAGS= -ftemplate-depth-120 -fcoalesce-templates -Wno-long-double

# CCTBX locations
CCTBX=${CLIB}/cctbx
CCTBX_sources=${CCTBX}/cctbx_sources

# Clipper location
CLPR=${CLIB}/clipper
#CLPR=${Clipper}

# Includes
ICCP4=-I${CLIBS}

ICLPR= \
-I${CCP4SRC}/include

#*# version before 6.1.24
#*#ITBX=\
#*#-I${CCTBX_sources}/boost \
#*#-I${CCP4}/include \
#*#-I${CCTBX_sources}/cctbx/include
#*#-I${CCTBX_sources}/scitbx/include

# 6.4.0
ITBX=\
-I${CCP4SRC}/include/cctbx

# 6.3.0
#ITBX=\
#-I${CCTBX_sources}/boost \
#-I${CCP4}/include \
#-I${CCTBX_sources}/cctbx_project


#MINLIB = ./Minimise
MINLIB = ./MinLib

ILIB = -I${MINLIB}


# TNT
##TNT = /Users/pre/Projects/Xtal/src/TNT
TNT = ${CCTBX_sources}/tntbx/include
ITNT = -I${TNT}

INCLUDE=-I. ${ILIB} ${ITBX} ${ICLPR} ${ICCP4} ${ITNT}

# Libraries
LMIN=${MINLIB}/libbfgs.a

LCCP4=-L${CLIB} -lrfftw -lfftw -lccp4f -lccp4c -lmmdb -lm

LCLPR=-L${CLIB} -lclipper-ccp4 -lclipper-contrib -lclipper-minimol -lclipper-mmdb -lclipper-core 

LTBX= -L${CLIB} -lcctbx 
#LTBX= -L${CCTBX_build}/libtbx -L${CCTBX_sources}/libtbx -lcctbx  -lm

SLIB = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
LDLIBS=${LMIN} ${LTBX} ${LCLPR} ${LCCP4} ${SLIB}

CPPFLAGS=${XCPPFLAGS}  ${INCLUDE}

EXE = aimless

OBJ =   \
	InputAll_aimless.o analyseanom.o analysesd.o \
	anomdistribution.o applyscales.o  \
	cumulativecompleteness.o  file_util.o foxholmes.o \
	fracdevanal.o globalcontrols_aimless.o halfdataset.o \
	initialscales.o intensitybin.o keywords_aimless.o  \
	mergedlist.o normalprobanal.o normprobfunc.o optimisesdcorr.o \
	plotfiles.o printing.o refinescale.o refinesdcorrection.o \
	reject.o resolutionlimit.o samplegaussian.o scalemodel.o \
	scalerefine.o scalerefinefh.o scaletypes.o sdanalysis.o \
	sdctypes.o sdmodel.o selectedobservations.o \
	selectscalingreflections.o selectsdcorrreflections.o \
	sphericalharmonic.o statistics.o summarystatistics.o \
	tie.o writeoutputfiles.o writeunmerged.o \
	simplex-lib.o optimisecombine.o tile.o imagearray.o \
	restore.o anisotropy.o intensity_scale.o \
	dataset.o observationstatuscontrol.o analyseoverlaps.o \
	histogram.o radiationdamageanalysis.o \
	hkl_merge.o anisotropicmodel.o \
	referencelist.o referencescalemodel.o refinereferencescale.o \
	refinetargets.o

COBJ = \
	CCP4base.o Errors.o Output.o Preprocessor.o \
	cctbx_utils.o cellgroup.o checkcompatiblesymmetry.o \
	columnlabels.o controls.o eprob.o fileread.o \
	getsubgroups.o hash.o hkl_controls.o hkl_datatypes.o \
	hkl_symmetry.o hkl_unmerge.o icering.o \
	interpretcommandline.o jiffy.o latsym.o crystaltype.o \
	linearlsq.o matvec_utils.o \
	mtz_merge_io.o mtz_unmerge_io.o mtz_utils.o normalise.o \
	numbercomplete.o observationflags.o \
	openinputfile.o pointgroup.o probfunctions.o range.o \
	rotation.o runthings.o scala_util.o score_datatypes.o \
	spacegroupreindex.o spline.o string_util.o timer.o \
	zone.o tablegraph.o

SOURCE = ${EXE}.cpp \
	InputAll_aimless.cpp analyseanom.cpp analysesd.cpp \
	anomdistribution.cpp applyscales.cpp  \
	cumulativecompleteness.cpp  file_util.cpp foxholmes.cpp \
	fracdevanal.cpp globalcontrols_aimless.cpp halfdataset.cpp \
	initialscales.cpp intensitybin.cpp keywords_aimless.cpp  \
	mergedlist.cpp normalprobanal.cpp normprobfunc.cpp optimisesdcorr.cpp \
	plotfiles.cpp printing.cpp refinescale.cpp refinesdcorrection.cpp \
	reject.cpp resolutionlimit.cpp samplegaussian.cpp scalemodel.cpp \
	scalerefine.cpp scalerefinefh.cpp scaletypes.cpp sdanalysis.cpp \
	sdctypes.cpp sdmodel.cpp selectedobservations.cpp \
	selectscalingreflections.cpp selectsdcorrreflections.cpp \
	sphericalharmonic.cpp statistics.cpp summarystatistics.cpp \
	tie.cpp writeoutputfiles.cpp writeunmerged.cpp \
	simplex-lib.cpp optimisecombine.cpp tile.cpp imagearray.cpp \
	restore.cpp anisotropy.cpp intensity_scale.cpp \
	dataset.cpp observationstatuscontrol.cpp analyseoverlaps.cpp \
	histogram.cpp radiationdamageanalysis.cpp \
	hkl_merge.cpp anisotropicmodel.cpp \
	referencelist.cpp referencescalemodel.cpp refinereferencescale.cpp \
	refinetargets.cpp


# Common with Pointless
CSOURCE = \
	CCP4base.cpp Errors.cpp Output.cpp Preprocessor.cpp \
	cctbx_utils.cpp cellgroup.cpp checkcompatiblesymmetry.cpp \
	columnlabels.cpp controls.cpp eprob.cpp fileread.cpp \
	getsubgroups.cpp hash.cpp hkl_controls.cpp hkl_datatypes.cpp \
	hkl_symmetry.cpp hkl_unmerge.cpp icering.cpp \
	interpretcommandline.cpp jiffy.cpp latsym.cpp crystaltype.cpp \
	linearlsq.cpp matvec_utils.cpp \
	mtz_merge_io.cpp mtz_unmerge_io.cpp mtz_utils.cpp normalise.cpp \
	numbercomplete.cpp observationflags.cpp \
	openinputfile.cpp pointgroup.cpp probfunctions.cpp range.cpp \
	rotation.cpp runthings.cpp scala_util.cpp score_datatypes.cpp \
	spacegroupreindex.cpp spline.cpp string_util.cpp timer.cpp \
	zone.cpp tablegraph.cpp

# Headers
# Just for Aimless
HDR = ${EXE}.hh \
	InputAll_aimless.hh analyseanom.hh analysesd.hh \
	anomdistribution.hh applyscales.hh  \
	cumulativecompleteness.hh  file_util.hh foxholmes.hh \
	fracdevanal.hh globalcontrols_aimless.hh halfdataset.hh \
	initialscales.hh intensitybin.hh keywords_aimless.hh  \
	mergedlist.hh normalprobanal.hh normprobfunc.hh optimisesdcorr.hh \
	plotfiles.hh printing.hh refinescale.hh refinesdcorrection.hh \
	reject.hh resolutionlimit.hh samplegaussian.hh scalemodel.hh \
	scalerefine.hh scalerefinefh.hh scaletypes.hh sdanalysis.hh \
	sdctypes.hh sdmodel.hh selectedobservations.hh \
	selectscalingreflections.hh selectsdcorrreflections.hh \
	sphericalharmonic.hh statistics.hh summarystatistics.hh \
	tie.hh writeoutputfiles.hh writeunmerged.hh \
	optimisecombine.hh tile.hh imagearray.hh \
	restore.hh anisotropy.hh  \
	dataset.hh observationstatuscontrol.hh analyseoverlaps.hh \
	simplex-lib.h intensity_scale.h intensity_target.h  \
	weighttype.hh histogram.hh radiationdamageanalysis.hh \
	hkl_merge.hh anisotropicmodel.hh \
	referencelist.hh referencescalemodel.hh refinereferencescale.hh \
	refinetargets.hh


# common with Pointless
CHDR = \
	CCP4base.hh Errors.hh Output.hh Preprocessor.hh \
	cctbx_utils.hh cellgroup.hh checkcompatiblesymmetry.hh \
	columnlabels.hh controls.hh eprob.hh fileread.hh \
	getsubgroups.hh hash.hh hkl_controls.hh hkl_datatypes.hh \
	hkl_symmetry.hh hkl_unmerge.hh icering.hh \
	interpretcommandline.hh jiffy.hh latsym.hh crystaltype.hh \
	linearlsq.hh matvec_utils.hh \
	mtz_merge_io.hh mtz_unmerge_io.hh mtz_utils.hh normalise.hh \
	numbercomplete.hh observationflags.hh \
	openinputfile.hh pointgroup.hh probfunctions.hh range.hh \
	rotation.hh runthings.hh scala_util.hh score_datatypes.hh \
	spacegroupreindex.hh spline.hh string_util.hh timer.hh \
	zone.hh tablegraph.hh \
	InputBase.hh InputAll.hh util.hh globalcontrols.hh \
	phaser_types.hh alt_hkl_datatypes.h version.hh \
	ProtocolScale.hh

OBJECTS = ${EXE}.o ${OBJ} ${COBJ}

MINLIBS = ${MINLIB}

${EXE}: ${OBJECTS} ${LMIN}
	${LXX} ${LFLAGS} -o ${EXE} ${OBJECTS}  ${LDLIBS}

${LMIN}:
	cd ${MINLIB}; make 


TTOBJ= t.o  eprob.o string_util.o range.o fileread.o
TTEXE= t

t: ${TTOBJ}
	${CXX} $(LFLAGS) -o ${TTEXE} ${TTOBJ}   ${LDLIBS}

TESTOBJ= test.o hkl_datatypes.o matvec_utils.o rotation.o scala_util.o \
 string_util.o range.o cone.o resolutionlimit.o


TESTEXE= test

test: ${TESTOBJ}
	${CXX} $(LFLAGS) -o ${TESTEXE} ${TESTOBJ}  ${LDLIBS}

BOBJ= simulateB.o hkl_datatypes.o  hkl_unmerge.o \
	scala_util.o Output.o matvec_utils.o jiffy.o range.o \
	hkl_symmetry.o rotation.o controls.o hash.o lattice.o runthings.o \
	observationflags.o icering.o writeunmergedmtz.o cellgroup.o mtz_utils.o

BEXE= simulateB

${BEXE}: ${BOBJ}
	${CXX} $(LFLAGS) -o ${BEXE} ${BOBJ}  ${LDLIBS}

UOBJ= simulate.o hkl_datatypes.o  hkl_unmerge.o \
	scala_util.o Output.o matvec_utils.o jiffy.o range.o \
	hkl_symmetry.o rotation.o hash.o lattice.o runthings.o \
	observationflags.o icering.o writeunmerged.o cellgroup.o \
	mtz_utils.o string_util.o  score_datatypes.o spacegroupreindex.o \
	probfunctions.o  cctbx_utils.o pointgroup.o zone.o latsym.o \
	getsubgroups.o samplegaussian.o eprob.o file_util.o \
	sdmodel.o sdctypes.o selectedobservations.o normalise.o \
	keywords_aimless.o CCP4base.o Errors.o  Preprocessor.o controls.o \
	linearlsq.o hkl_merged_list.o spline.o mtz_merge_io.o \
	globalcontrols_aimless.o \
	columnlabels.o 


UEXE= simulate

${UEXE}: ${UOBJ}
	${CXX} $(LFLAGS) -o ${UEXE} ${UOBJ}  ${LDLIBS}

SOBJ= simulate2.o hkl_datatypes.o score_datatypes.o  normprobfunc.o rotation.o \
	scala_util.o Output.o matvec_utils.o jiffy.o \
	normalprobanal.o plotfiles.o file_util.o range.o \
	fracdevanal.o intensitybin.o samplegaussian.o 

####
VEXE= simdelta 

VOBJ= simdelta.o samplegaussian.o scala_util.o \
  hkl_datatypes.o string_util.o  matvec_utils.o range.o rotation.o \
  score_datatypes.o spacegroupreindex.o hkl_symmetry.o latsym.o lattice.o \
  pointgroup.o zone.o getsubgroups.o probfunctions.o Output.o jiffy.o \
  cctbx_utils.o sdctypes.o


${VEXE}: ${VOBJ}
	${CXX} $(LFLAGS) -o ${VEXE} ${VOBJ}  ${LDLIBS}
####
V1EXE= simdelta1 

V1OBJ= simdelta1.o samplegaussian.o scala_util.o \
  hkl_datatypes.o string_util.o  matvec_utils.o range.o rotation.o \
  score_datatypes.o spacegroupreindex.o hkl_symmetry.o latsym.o lattice.o \
  pointgroup.o zone.o getsubgroups.o probfunctions.o Output.o jiffy.o \
  cctbx_utils.o sdctypes.o


${V1EXE}: ${V1OBJ}
	${CXX} $(LFLAGS) -o ${V1EXE} ${V1OBJ}  ${LDLIBS}
####

QEXE= studyparams

QOBJ= studyparams.o tablegraph.o string_util.o range.o

${QEXE}: ${QOBJ}
	${CXX} $(LFLAGS) -o ${QEXE} ${QOBJ}  ${LDLIBS}

SEXE= simulate2

${SEXE}: ${SOBJ}
	${CXX} $(LFLAGS) -o ${SEXE} ${SOBJ}  ${LDLIBS}



.SUFFIXES: .cpp .o

.cpp.o: ${HDR}
	${CXX} -c $(CFLAGS) $(CPPFLAGS) $<

.c.o:
	${CC} -c $(CFLAGS) $(CCFLAGS) $<

clean: 
	-rm *.o ${EXE} *.gch

allclean: clean
	cd ${MINLIB}; make clean

####################################################################
# Export etc

export:
	/bin/rm -rf ../Export/
	mkdir ../Export
	cp ${SOURCE} ../Export/
	cp ${CSOURCE} ../Export/
	cp ${HDR} ../Export/
	cp ${CHDR} ../Export/
	cp Makefile ../Export/Makefile.make
	cp Makefile.osx ../Export/Makefile.osx.make
	cp aimless.html ../Export/
	cd ${MINLIB}; make export

tar: export
	cd ../Export; tar cvf ../aimless.tar .
	cd ..; gzip -f aimless.tar

