cmake_minimum_required(VERSION 2.8)
project(aimless C CXX)
find_package(CCP4 REQUIRED mmdb ccp4c clipper-core clipper-mmdb clipper-minimol
                           clipper-contrib clipper-ccp4 cctbx)

include_directories(${CCP4_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/MinLib)

set(MINLIB_SOURCES
    MinLib/RefineBase.cc MinLib/RefineBase2.cc MinLib/Minimizer.cc
    MinLib/PhaserError.cc MinLib/Errors.cc)

# Unique to Aimless
set(AIMLESS_SOURCES
    InputAll_aimless.cpp aimless.cpp analyseanom.cpp analysesd.cpp
    anomdistribution.cpp applyscales.cpp
    cumulativecompleteness.cpp  file_util.cpp foxholmes.cpp
    fracdevanal.cpp globalcontrols_aimless.cpp halfdataset.cpp
    initialscales.cpp intensitybin.cpp keywords_aimless.cpp
    mergedlist.cpp normalprobanal.cpp normprobfunc.cpp optimisesdcorr.cpp
    plotfiles.cpp printing.cpp refinescale.cpp refinesdcorrection.cpp
    reject.cpp resolutionlimit.cpp samplegaussian.cpp scalemodel.cpp
    scalerefine.cpp scalerefinefh.cpp scaletypes.cpp sdanalysis.cpp
    sdctypes.cpp sdmodel.cpp selectedobservations.cpp
    selectscalingreflections.cpp selectsdcorrreflections.cpp
    sphericalharmonic.cpp statistics.cpp summarystatistics.cpp
    tie.cpp writeoutputfiles.cpp writeunmerged.cpp
    simplex-lib.cpp optimisecombine.cpp tile.cpp imagearray.cpp
    restore.cpp anisotropy.cpp intensity_scale.cpp
)

# Shared (duplicated) with Pointless
set(SHARED_SOURCES
    CCP4base.cpp Errors.cpp Output.cpp Preprocessor.cpp
    cctbx_utils.cpp cellgroup.cpp checkcompatiblesymmetry.cpp
    columnlabels.cpp controls.cpp eprob.cpp fileread.cpp
    getsubgroups.cpp hash.cpp hkl_controls.cpp hkl_datatypes.cpp
    hkl_merged_list.cpp hkl_symmetry.cpp hkl_unmerge.cpp icering.cpp
    interpretcommandline.cpp jiffy.cpp latsym.cpp lattice.cpp
    linearlsq.cpp matvec_utils.cpp
    mtz_merge_io.cpp mtz_unmerge_io.cpp mtz_utils.cpp normalise.cpp
    numbercomplete.cpp observationflags.cpp
    openinputfile.cpp pointgroup.cpp probfunctions.cpp range.cpp
    rotation.cpp runthings.cpp scala_util.cpp score_datatypes.cpp
    spacegroupreindex.cpp spline.cpp string_util.cpp timer.cpp
    zone.cpp tablegraph.cpp
)

add_executable(aimless ${AIMLESS_SOURCES} ${SHARED_SOURCES} ${MINLIB_SOURCES})

target_link_libraries(aimless ${CCP4_LIBRARIES})

install(TARGETS aimless DESTINATION bin)

