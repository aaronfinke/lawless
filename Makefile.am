SUBDIRS = MinLib

aimless_SOURCES = \
        aimless.cpp \
        InputAll_aimless.cpp analyseanom.cpp analysesd.cpp \
        anomdistribution.cpp applyscales.cpp \
        cumulativecompletenes.cpp  file_util.cpp foxholmes.cpp \
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
        simplex-lib.cpp \
        CCP4base.cpp Errors.cpp Output.cpp Preprocessor.cpp \
        cctbx_utils.cpp cellgroup.cpp checkcompatiblesymmetry.cpp \
        columnlabels.cpp controls.cpp eprob.cpp fileread.cpp \
        getsubgroups.cpp hash.cpp hkl_controls.cpp hkl_datatypes.cpp \
        hkl_merged_list.cpp hkl_symmetry.cpp hkl_unmerge.cpp icering.cpp \
        interpretcommandline.cpp jiffy.cpp latsym.cpp lattice.cpp \
        linearlsq.cpp matvec_utils.cpp \
        mtz_merge_io.cpp mtz_unmerge_io.cpp mtz_utils.cpp normalise.cpp \
        numbercomplete.cpp observationflags.cpp \
        openinputfile.cpp pointgroup.cpp probfunctions.cpp range.cpp \
        rotation.cpp runthings.cpp scala_util.cpp score_datatypes.cpp \
        spacegroupreindex.cpp spline.cpp string_util.cpp timer.cpp \
        zone.cpp tablegraph.cpp optimisecombine.cpp tile.cpp imagearray.cpp \
        restore.cpp anisotropy.cpp intensity_scale.cpp \
        aimless.hh \
        InputAll.hh InputAll_aimless.hh ProtocolSDref.hh ProtocolScale.hh \
        analyseanom.hh analysesd.hh anomdistribution.hh \
        applyscales.hh cumulativecompleteness.hh file_util.hh \
        foxholmes.hh fracdevanal.hh globalcontrols.hh \
        globalcontrols_aimless.hh halfdataset.hh initialscales.hh \
        intensitybin.hh keywords_aimless.hh mergedlist.hh normalprobanal.hh \
        normprobfunc.hh optimisesdcorr.hh plotfiles.hh printing.hh \
        refinescale.hh refinesdcorrection-bfgs.hh refinesdcorrection.hh \
        reject.hh resolutionlimit.hh samplegaussian.hh scalemodel.hh \
        scalerefine.hh scalerefinefh.hh scaletypes.hh sdanalysis.hh \
        sdctypes.hh sdmodel.hh selectedobservations.hh \
        selectscalingreflections.hh selectsdcorrreflections.hh \
        sphericalharmonic.hh spline.hh statistics.hh summarystatistics.hh \
        tie.hh version.hh writeoutputfiles.hh writeunmerged.hh \
        simplex-lib.h alt_hkl_datatypes.h \
        CCP4base.hh Errors.hh InputBase.hh Output.hh Preprocessor.hh \
        cctbx_utils.hh cellgroup.hh checkcompatiblesymmetry.hh columnlabels.hh \
        controls.hh eprob.hh fileread.hh getsubgroups.hh hash.hh \
        hkl_controls.hh hkl_datatypes.hh hkl_merged_list.hh hkl_symmetry.hh \
        hkl_unmerge.hh icering.hh interpretcommandline.hh jiffy.hh latsym.hh \
        lattice.hh linearlsq.hh matvec_utils.hh mtz_merge_io.hh \
        mtz_unmerge_io.hh mtz_utils.hh normalise.hh numbercomplete.hh \
        observationflags.hh openinputfile.hh phaser_types.hh pointgroup.hh \
        probfunctions.hh range.hh rotation.hh runthings.hh scala.hh \
        scala_util.hh score_datatypes.hh spacegroupreindex.hh spline.hh \
        string_util.hh timer.hh util.hh zone.hh tablegraph.hh \
	optimisecombine.hh tile.hh imagearray.hh restore.hh \
        anisotropy.hh intensity_scale.h


bin_PROGRAMS = aimless

INCLUDES = -I@srcdir@/MinLib $(PTHREAD_CFLAGS) $(CLIPPER_CXXFLAGS) $(CCTBX_CXXFLAGS) $(CCP4_CXXFLAGS) $(FFTW_CXXFLAGS)

aimless_LDADD = $(PTHREAD_CFLAGS) $(CLIPPER_LIBS) $(CCTBX_LIBS) $(CCP4_LIBS) $(FFTW_LIBS) $(PTHREAD_LIBS) $(FCLIBS) ./MinLib/libbfgs.a

EXTRA_DIST = macros/ccp4.m4 macros/cctbx.m4 macros/clipper.m4 aimless.html
