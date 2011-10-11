
# Process this file with autoconf to produce a configure script.

dnl package-name version-number bugs-to:
AC_INIT(aimless, 0.1.2, pre@lmb-mrc.ac.uk) 

AC_CANONICAL_TARGET

dnl package-name version-number, the gzip'd tar file is named
dnl package-version.tar.gz.  
AM_INIT_AUTOMAKE(aimless, 0.1.2)

AC_CONFIG_SRCDIR([aimless.cpp])

AM_MAINTAINER_MODE

AC_FC_LIBRARY_LDFLAGS

AC_LANG(C++)

# Checks for programs.
AC_PROG_CXX(xlc++ CC cxx c++ icc g++)
AC_HEADER_STDC
AC_PROG_INSTALL
AC_PROG_RANLIB

AX_PTHREAD(,
   [ AC_MSG_ERROR([Failed to find threading library])
   ] )

AX_OPENMP([ AC_MSG_WARN([OpenMP enabled])
   AC_DEFINE([_OPENMP],[1],[Define for OPENMP compilation])]
   CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
   ,[AC_MSG_WARN([OpenMP disabled])])
   
AC_AIMLESS_OPTIONS

AC_ARG_ENABLE(debug,
  AC_HELP_STRING( [--enable-debug], [ enable debugging options ] ),
  [ac_cv_enable_debug=yes ],
  [ac_cv_enable_debug=no ] )

AM_PATH_CCP4(,
  [ if test "x$ccp4_prefix" != x ; then
     if test "x$ac_CCP4_CXXFLAGS" != x ; then 
       case "$ac_CCP4_CXXFLAGS" in
        *include/ccp4 )
        ac_CCP4_LDOPTS=`echo "$ac_CCP4_CXXFLAGS" | sed s#include/ccp4#lib#g`
        ;;
        *)
        ac_CCP4_LDOPTS="$ac_CCP4_CXXFLAGS"
        esac
        ac_CCP4_LDOPTS=`echo "$ac_CCP4_LDOPTS" | sed s#-I#-L#`
        CCP4_CXXFLAGS="$ac_CCP4_CXXFLAGS" 
        CCP4_LIBS="$ac_CCP4_LDOPTS -lccp4f -lccp4c $ac_CCP4_LDOPTS/../ccif -lccif -lmmdb"
        AC_MSG_WARN([assuming library $CCP4_LIBS headers $CCP4_CXXFLAGS])
        have_ccp4=yes
      else
        CCP4_CXXFLAGS="-I$ccp4_prefix/lib/src -I$ccp4_prefix/lib/libccp4/ccp4"
        CCP4_LIBS="-L$ccp4_prefix/lib/src -L$ccp4_prefix/lib -lccp4f -lccp4c -lccif -lmmdb"
        AC_MSG_WARN([assuming library $CCP4_LIBS headers $CCP4_CXXFLAGS])
        have_ccp4=yes
      fi
    else
      AC_MSG_ERROR([Failed to find ccp4 libs])
    fi
  ]
)

AM_PATH_CCTBX(, 
 [ if test "x$cctbx_prefix" != x ; then
     if test "x$ac_CCTBX_CXXFLAGS" != x ; then 
       if test "x$cctbx_build_prefix" != "x$cctbx_prefix" ; then
         CCTBX_CXXFLAGS="$ac_CCTBX_CXXFLAGS -I$cctbx_build_prefix/cctbx_project/include -I$cctbx_build_prefix/include" 
         CCTBX_LIBS="$CCTBX_LIBS -L$cctbx_build_prefix/lib -L$cctbx_build_prefix/cctbx_project/lib -lcctbx"
         AC_MSG_WARN([assuming library $CCTBX_LIBS headers $CCTBX_CXXFLAGS])
         have_cctbx=yes
       else
         AC_MSG_ERROR([Unable to locate cctbx build directory])
       fi
      else
        AC_MSG_ERROR([Unable to locate cctbx source files])
      fi
    else
      AC_MSG_ERROR([Failed to find cctbx libs])
    fi
  ]
) 

dnl Check for FFTW.
AM_PATH_FFTW(, AC_MSG_ERROR([locating necessary fftw]) )

AM_PATH_CLIPPER(,
     [ if test "x$clipper_prefix" != x ; then
     if test "x$ac_CLIPPER_CXXFLAGS" != x ; then
       case "$ac_CLIPPER_CXXFLAGS" in
        *include/clipper )
        ac_CLIPPER_LDOPTS=`echo "$ac_CLIPPER_CXXFLAGS" | sed s#include/clipper#lib#g`
        ;;
        *)
        ac_CLIPPER_LDOPTS="$ac_CLIPPER_CXXFLAGS"
        esac
        ac_CLIPPER_LDOPTS=`echo "$ac_CLIPPER_LDOPTS" | sed s#-I#-L#`
        CLIPPER_CXXFLAGS="$ac_CLIPPER_CXXFLAGS"
        CLIPPER_LIBS="$ac_CLIPPER_LDOPTS -lclipper-ccp4 -lclipper-contrib -lclipper-minimol -lclipper-mmdb -lclipper-core"
        AC_MSG_WARN([assuming library $CLIPPER_LIBS headers $CLIPPER_CXXFLAGS])
        have_clipper=yes
      else
        CLIPPER_CXXFLAGS="-I$clipper_prefix/lib/clipper"
        CLIPPER_LIBS="-L$clipper_prefix/lib -L$clipper_prefix/lib/clipper -lclipper-ccp4 -lclipper-contrib -lclipper-minimol -lclipper-mmdb -lclipper-core"
        AC_MSG_WARN([assuming library $CLIPPER_LIBS headers $CLIPPER_CXXFLAGS])
        have_clipper=yes
      fi
    else
     AC_MSG_ERROR([Failed to find clipper ])
    fi
   ]    
) 

# Checks for typedefs, structures, and compiler characteristics.
#
# Does clipper use these?  I'm guessing not, commented.
# AC_C_CONST
# AC_C_INLINE
# AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor pow rint sqrt])

# test for function in m
AC_CHECK_LIB([m], [sqrt]) 

AC_SUBST(AR_FLAGS)
AM_CONDITIONAL(DEBUG, test "$ac_cv_enable_debug" = yes)
AM_CONDITIONAL(AMDEP, test "x$enable_dependency_tracking" != xno)

AC_SUBST(aimless_LDFLAGS)

AC_CONFIG_FILES([MinLib/Makefile
                Makefile])

AC_OUTPUT
