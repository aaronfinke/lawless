cmake_minimum_required(VERSION 3.16)
project(aimless C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)

# ---------------------------------------------------------------------------
# Dependencies
#
# Linux (Ubuntu/Debian): everything comes from distribution packages on the
#   system paths, and none of the variables below is needed.
# macOS: the clipper, ccp4c and mmdb2 libraries, Eigen and the tntbx headers
#   come from a CCP4 installation (CCP4_DIR); cctbx from a conda env with
#   cctbx-base (CCTBX_DIR); boost headers from Homebrew. CCP4 does not install
#   the clipper, mmdb2 or libccp4 headers, so CCP4_SRC must point at a
#   directory holding them: a CCP4 source checkout, or the unpacked libccp4,
#   mmdb2 and clipper tarballs from https://ftp.ccp4.ac.uk/opensource/.
#
# Each variable may also be given as an environment variable of the same
# name. CCP4_DIR falls back to $CCP4 and then to the newest /Applications/ccp4-*;
# CCTBX_DIR falls back to $CONDA_PREFIX.
# ---------------------------------------------------------------------------
if(NOT DEFINED CCP4_DIR)
  if(DEFINED ENV{CCP4_DIR})
    set(_ccp4dir "$ENV{CCP4_DIR}")
  elseif(DEFINED ENV{CCP4})
    set(_ccp4dir "$ENV{CCP4}")
  elseif(APPLE)
    file(GLOB _ccp4dirs /Applications/ccp4-*)
    if(_ccp4dirs)
      list(SORT _ccp4dirs)
      list(GET _ccp4dirs -1 _ccp4dir)
    endif()
  endif()
endif()
if(NOT DEFINED CCTBX_DIR)
  if(DEFINED ENV{CCTBX_DIR})
    set(_cctbxdir "$ENV{CCTBX_DIR}")
  else()
    set(_cctbxdir "$ENV{CONDA_PREFIX}")
  endif()
endif()
set(CCP4_DIR  "${_ccp4dir}"     CACHE PATH "CCP4 installation (libraries, Eigen, tntbx)")
set(CCP4_SRC  "$ENV{CCP4_SRC}"  CACHE PATH "Directory holding the clipper, mmdb2 and libccp4 headers")
set(CCTBX_DIR "${_cctbxdir}"    CACHE PATH "Prefix of a conda env containing cctbx-base")

# Header search locations, most specific first. In CCP4_SRC the packages sit
# one level down (clipper/, mmdb2-2.0.22/, ...), so each subdirectory is a hint.
set(_hints)
if(CCP4_SRC)
  foreach(_d clipper libccp4 mmdb2 boost_1_74_0 eigen-3.4.0
             cctbx-packages/tntbx/include)
    list(APPEND _hints ${CCP4_SRC}/${_d})
  endforeach()
  file(GLOB _tarballs LIST_DIRECTORIES true
       ${CCP4_SRC}/clipper-* ${CCP4_SRC}/libccp4-* ${CCP4_SRC}/mmdb2-*)
  list(APPEND _hints ${_tarballs})
endif()
set(_ccp4hints)
set(_cctbxhints)
if(CCP4_DIR)
  # tntbx is only installed inside CCP4's Python (a framework on macOS)
  file(GLOB _ccp4hints LIST_DIRECTORIES true
       ${CCP4_DIR}/Frameworks/Python.framework/Versions/*/lib/python*/site-packages/tntbx/include
       ${CCP4_DIR}/lib/python*/site-packages/tntbx/include)
  list(APPEND _ccp4hints ${CCP4_DIR}/include)
endif()
if(CCTBX_DIR)
  # conda installs the cctbx and scitbx headers alongside the Python modules
  file(GLOB _cctbxhints LIST_DIRECTORIES true ${CCTBX_DIR}/lib/python*/site-packages)
  list(APPEND _cctbxhints ${CCTBX_DIR}/include)
endif()

find_path(CLIPPER_INCLUDE_DIR clipper/clipper.h       HINTS ${_hints})
find_path(CCP4C_INCLUDE_DIR   ccp4/cmtzlib.h          HINTS ${_hints})
find_path(MMDB2_INCLUDE_DIR   mmdb2/mmdb_manager.h    HINTS ${_hints})
find_path(BOOST_INCLUDE_DIR   boost/version.hpp       HINTS ${_hints}
          PATHS /opt/homebrew/opt/boost/include /usr/local/opt/boost/include)
find_path(EIGEN3_INCLUDE_DIR  Eigen/Core              HINTS ${_hints} ${_ccp4hints}
          PATH_SUFFIXES eigen3)
find_path(TNT_INCLUDE_DIR     tnt_array1d.h           HINTS ${_hints} ${_ccp4hints}
          PATH_SUFFIXES tnt)
find_path(JAMA_INCLUDE_DIR    jama_eig.h              HINTS ${_hints} ${_ccp4hints}
          PATH_SUFFIXES jama)
find_path(CCTBX_INCLUDE_DIR   cctbx/sgtbx/basic.h     HINTS ${_cctbxhints})
find_path(SCITBX_INCLUDE_DIR  scitbx/math/eigensystem.h HINTS ${_cctbxhints})

set(_libvars)
foreach(_lib clipper-ccp4 clipper-contrib clipper-minimol clipper-mmdb
             clipper-core ccp4c mmdb2)
  string(MAKE_C_IDENTIFIER ${_lib} _id)
  string(TOUPPER ${_id} _id)
  find_library(${_id}_LIBRARY ${_lib} HINTS ${CCP4_DIR}/lib)
  list(APPEND _libvars ${_id}_LIBRARY)
endforeach()
find_library(CCTBX_LIBRARY cctbx HINTS ${CCTBX_DIR}/lib)
list(APPEND _libvars CCTBX_LIBRARY)

set(_missing)
foreach(_v CLIPPER_INCLUDE_DIR CCP4C_INCLUDE_DIR MMDB2_INCLUDE_DIR
           BOOST_INCLUDE_DIR EIGEN3_INCLUDE_DIR TNT_INCLUDE_DIR JAMA_INCLUDE_DIR
           CCTBX_INCLUDE_DIR SCITBX_INCLUDE_DIR ${_libvars})
  if(${_v})
    message(STATUS "${_v}: ${${_v}}")
  else()
    list(APPEND _missing ${_v})
  endif()
endforeach()
if(_missing)
  string(REPLACE ";" " " _missing "${_missing}")
  message(FATAL_ERROR "Not found: ${_missing}\n"
    "CCP4_DIR='${CCP4_DIR}' CCP4_SRC='${CCP4_SRC}' CCTBX_DIR='${CCTBX_DIR}'\n"
    "See \"Building\" in README.md.")
endif()

set(CCP4_LIBRARIES)
foreach(_v ${_libvars})
  list(APPEND CCP4_LIBRARIES ${${_v}})
endforeach()

# boost/timer.hpp is an #error in recent boost unless this is defined
add_compile_definitions(BOOST_TIMER_ENABLE_DEPRECATED)

#OpenMP
FIND_PACKAGE(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Boost comes before cctbx: a conda include/ can hold a partial boost tree
include_directories(${MMDB2_INCLUDE_DIR} ${CLIPPER_INCLUDE_DIR} ${CCP4C_INCLUDE_DIR}
                    ${BOOST_INCLUDE_DIR} ${TNT_INCLUDE_DIR} ${JAMA_INCLUDE_DIR})
# Eigen3 (header-only) for GPR wavelength normalisation
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${CCTBX_INCLUDE_DIR} ${SCITBX_INCLUDE_DIR})
if(CCTBX_DIR AND EXISTS ${CCTBX_DIR}/include)
  include_directories(${CCTBX_DIR}/include)
endif()
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
    plotfiles.cpp printing.cpp probe.cpp refinescale.cpp refinesdcorrection.cpp
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

# Keep the library directories in the installed binary's rpath, so that an
# installed lawless finds libcctbx (an @rpath install name) as build/ does.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

add_executable(lawless ${AIMLESS_SOURCES} ${SHARED_SOURCES} ${MINLIB_SOURCES})

target_link_libraries(lawless ${CCP4_LIBRARIES})

# The CCP4 macOS dylibs carry the build machine's paths as their install
# names. Where such a path does not exist here, point the binary at the
# library's actual location instead.
if(APPLE)
  foreach(_lib ${CCP4_LIBRARIES})
    get_filename_component(_real ${_lib} REALPATH)
    execute_process(COMMAND otool -D ${_real}
                    OUTPUT_VARIABLE _out OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REGEX REPLACE ".*\n" "" _id "${_out}")
    if(_id MATCHES "^/" AND NOT EXISTS ${_id})
      get_filename_component(_dir ${_real} DIRECTORY)
      get_filename_component(_name ${_id} NAME)
      add_custom_command(TARGET lawless POST_BUILD
        COMMAND install_name_tool -change ${_id} ${_dir}/${_name}
                $<TARGET_FILE:lawless>
        COMMENT "Fixing install name: ${_name}")
    endif()
  endforeach()
endif()

# build/aimless -> lawless, so scripts that still call the old name keep working.
add_custom_command(TARGET lawless POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E create_symlink lawless
    $<TARGET_FILE_DIR:lawless>/aimless
  COMMENT "Linking aimless -> lawless")

install(TARGETS lawless DESTINATION bin)
