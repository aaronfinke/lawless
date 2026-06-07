cmake_minimum_required(VERSION 2.8)
project(aimless C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)

set(SRC   /Users/aaronfinke/ccp4-8.0-src/checkout)
set(CCTBX /Users/aaronfinke/miniforge3/envs/cctbx)

# Libraries: clipper/ccp4c/mmdb2 from CCP4-9; cctbx from conda cctbx env (arm64 dylib)
set(CCP4_LIBRARIES
    /Applications/ccp4-9/lib/libclipper-ccp4.dylib
    /Applications/ccp4-9/lib/libclipper-contrib.dylib
    /Applications/ccp4-9/lib/libclipper-minimol.dylib
    /Applications/ccp4-9/lib/libclipper-mmdb.dylib
    /Applications/ccp4-9/lib/libclipper-core.dylib
    /Applications/ccp4-9/lib/libccp4c.dylib
    /Applications/ccp4-9/lib/libmmdb2.dylib
    ${CCTBX}/lib/libcctbx.dylib
)

#OpenMP
FIND_PACKAGE(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Headers from ccp4-8.0-src source tree
include_directories(${SRC}/mmdb2)
include_directories(${SRC}/clipper)
include_directories(${SRC}/libccp4)
include_directories(${SRC}/boost_1_74_0)
include_directories(${SRC}/cctbx-packages/tntbx/include)
# cctbx/scitbx headers from conda cctbx env (complete, arm64-compiled set)
include_directories(${CCTBX}/include)
include_directories(${CCTBX}/lib/python3.12/site-packages)
include_directories(${CMAKE_SOURCE_DIR}/MinLib)

set(MINLIB_SOURCES
    MinLib/RefineBase.cc MinLib/RefineBase2.cc MinLib/Minimizer.cc
    MinLib/PhaserError.cc MinLib/Errors.cc)

# Unique to Aimless
set(AIMLESS_SOURCES
    InputAll_aimless.cpp aimless.cpp analyseanom.cpp analysesd.cpp
    analyseoverlaps.cpp anomdistribution.cpp applyscales.cpp
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
    radiationdamageanalysis.cpp
    anisotropicmodel.cpp histogram.cpp hkl_merge.cpp referencelist.cpp
    referencescalemodel.cpp refinereferencescale.cpp refinetargets.cpp
    rejectbatches.cpp fitresolutiondata.cpp hessian.cpp radialfunction.cpp
    simpleminimise.cpp comparesds.cpp
    batchgroup.cpp linearsdcorrection.cpp weighttype.cpp
    refinereferencescale.cpp refinescaleref.cpp
    secondaryscalestats.cpp comparetoreference.cpp
)

# Shared (duplicated) with Pointless
set(SHARED_SOURCES
    CCP4base.cpp Errors.cpp Output.cpp Preprocessor.cpp
    cctbx_utils.cpp cellgroup.cpp checkcompatiblesymmetry.cpp
    columnlabels.cpp controls.cpp crystaltype.cpp dataset.cpp eprob.cpp
    fileread.cpp getsubgroups.cpp hash.cpp hkl_controls.cpp hkl_datatypes.cpp
    hkl_symmetry.cpp hkl_unmerge.cpp icering.cpp
    interpretcommandline.cpp jiffy.cpp latsym.cpp
    linearlsq.cpp matvec_utils.cpp
    mtz_merge_io.cpp mtz_unmerge_io.cpp mtz_utils.cpp normalise.cpp
    numbercomplete.cpp observationflags.cpp observationstatuscontrol.cpp
    openinputfile.cpp pointgroup.cpp probfunctions.cpp range.cpp
    rotation.cpp runthings.cpp scala_util.cpp score_datatypes.cpp
    spacegroupreindex.cpp spline.cpp string_util.cpp timer.cpp
    zone.cpp tablegraph.cpp runcorrelations.cpp report_errors.cpp
)

add_executable(aimless ${AIMLESS_SOURCES} ${SHARED_SOURCES} ${MINLIB_SOURCES})

target_link_libraries(aimless ${CCP4_LIBRARIES})

# The CCP4-9 dylibs have Jenkins build-machine paths as their install names.
# Remap them to the actual locations at post-build so the binary runs directly from build/.
set(JENKINS_LIB /Users/Shared/ccb63/Jenkins/workspace/CCP4/release-9-osx-clang-m1/devtools/install/lib)
set(CCP4_LIB    /Applications/ccp4-9/lib)
foreach(_lib
    libclipper-ccp4.2    libclipper-contrib.2  libclipper-minimol.2
    libclipper-mmdb.2    libclipper-core.2     libccp4c.8.0          libmmdb2.0)
  add_custom_command(TARGET aimless POST_BUILD
    COMMAND install_name_tool -change
      ${JENKINS_LIB}/${_lib}.dylib
      ${CCP4_LIB}/${_lib}.dylib
      $<TARGET_FILE:aimless>
    COMMENT "Fixing rpath: ${_lib}")
endforeach()

install(TARGETS aimless DESTINATION bin)
