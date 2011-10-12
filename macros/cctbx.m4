m4_define([_AM_PATH_CCTBX_BUILD],
[
AC_ARG_WITH(cctbx-build,
  AC_HELP_STRING( [--with-cctbx-build=PFX], [location of cctbx build, if different from source] ),
  [
    test "$withval" = no || with_cctbx_build=yes
    test "$withval" = yes || cctbx_build_prefix="$withval" ],
  [ with_cctbx_build=yes ] ) #dnl required

saved_CXXFLAGS="$CXXFLAGS"
CCTBX_LIBS=""
CCTBX_CXXFLAGS=""

test "x$cctbx_build_prefix" = "x" && cctbx_build_prefix="$cctbx_prefix"

if test "x$cctbx_build_prefix" != x; then
ac_cctbx_dirs='
.
lib
include
cctbx_project/lib
cctbx_project/include'
ac_cctbx_os='
.
intel-linux
sun-sunos
alpha-tru64
itanium-linux
mac-linux
alpha-linux
mac-osx
build/intel-linux
build/sun-sunos
build/alpha-tru64
build/itanium-linux
build/mac-linux
build/alpha-linux
build/mac-osx'
for os_dir in $ac_cctbx_os; do
for ac_dir in $ac_cctbx_dirs; do
  if test -r "$cctbx_build_prefix/$os_dir/$ac_dir/scitbx/array_family/operator_traits_builtin.h"; then
    ac_CCTBX_CXXFLAGS="$ac_CCTBX_CXXFLAGS -I$cctbx_build_prefix/$os_dir/$ac_dir"
    break 2
    fi
  done
  done
for os_dir in $ac_cctbx_os; do
for ac_dir in $ac_cctbx_dirs; do
  for ac_extension in a so sl dylib; do
  if test -r "$cctbx_build_prefix/$os_dir/$ac_dir/libcctbx.$ac_extension"; then
    ac_CCTBX_LDOPTS="-L$cctbx_build_prefix/$os_dir/$ac_dir -lcctbx"
    break 3
    fi
  done
  done
  done
fi
])

m4_define([_AM_PATH_CCTBX_EXTRA],
[

AC_ARG_ENABLE(boost,
  AC_HELP_STRING( [--enable-boost], [enable linking to boost python interface (default no)] ),
  [
    case $enableval in
      yes) enable_boost=yes ;;
      no) enable_boost=no ;;
      *) enable_boost=yes
    esac ],
  [  enable_boost=no ] #dnl default is yes
)

if test $enable_boost = yes; then
AC_ARG_VAR(BOOST, [boost top dir -optional])

if test "x$ac_cv_env_BOOST_set" != xset; then
  if test "x$cctbx_prefix" != x; then
    BOOST="$cctbx_prefix/../boost"
  fi
fi

ac_CCTBX_CXXFLAGS="$ac_CCTBX_CXXFLAGS -I$BOOST"
#extend for systems that need it
case "$host_os" in
  *osf* | *64* | *irix* )
    ac_CCTBX_CXXFLAGS="$ac_CCTBX_CXXFLAGS -I$BOOST/boost/compatibility/ccp_c_headers"
esac

ac_CCTBX_LDOPTS="$ac_CCTBX_LDOPTS -lboost_python"
i
fi

CXXFLAGS="$CXXFLAGS $ac_CCTBX_CXXFLAGS"
AC_ARG_VAR(CCTBX_VERSION, [year of release of cctbx version -optional])
if test "x$ac_cv_env_CCTBX_VERSION_set" = xset; then
  ac_CCTBX_CXXFLAGS="$ac_CCTBX_CXXFLAGS -DCCTBX_VERSION=$CCTBX_VERSION"
else
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
#include "cctbx/sgtbx/symbols.h";
#include <iostream>;],
[cctbx::sgtbx::space_group_symbols SGS(19); std::cout << SGS.extended_hermann_mauguin() << std::endl;])],
[ac_cv_HAVE_EXTENDED_HERMANN_MANGUIN=yes
CCTBX_VERSION=2004
ac_CCTBX_CXXFLAGS="$ac_CCTBX_CXXFLAGS -DCCTBX_VERSION=2004"
],
[ac_cv_HAVE_EXTENDED_HERMANN_MANGUIN=no
CCTBX_VERSION=2007
ac_CCTBX_CXXFLAGS="$ac_CCTBX_CXXFLAGS -DCCTBX_VERSION=2007"
])
fi
])
 
# AM_PATH_CCTBX([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
AC_DEFUN([AM_PATH_CCTBX],
[
AC_PROVIDE([AM_PATH_CCTBX])

AC_ARG_WITH(cctbx,
  AC_HELP_STRING( [--with-cctbx=PFX], [location of cctbx libraries] ),
  [
    test "$withval" = no || with_cctbx=yes 
    test "$withval" = yes || cctbx_prefix="$withval" ],
  [ with_cctbx=yes ] ) #dnl required

if test "x${with_cctbx}" = xyes ; then  
AS_IF([test "x$CCTBX_LIBS" != x && test "x$CCTBX_CXXFLAGS" != x ],
[
  have_cctbx=yes
],
[
saved_LIBS="$LIBS"
saved_CXXFLAGS="$CXXFLAGS"
CCTBX_LIBS=""
CCTBX_CXXFLAGS=""

if test "x$cctbx_prefix" != x; then
ac_cctbx_dirs='
.
lib
include
build/cctbx/lib
cctbx/include
scitbx/include
tntbx/include
cctbx_project
boost'
for ac_dir in $ac_cctbx_dirs; do
  if test -r "$cctbx_prefix/$ac_dir/cctbx/miller.h"; then
    ac_CCTBX_CXXFLAGS="-I$cctbx_prefix/$ac_dir"
    break
    fi
  done
for ac_dir in $ac_cctbx_dirs; do
  if test -r "$cctbx_prefix/$ac_dir/scitbx/vec3.h"; then
    ac_CCTBX_CXXFLAGS="$ac_CCTBX_CXXFLAGS -I$cctbx_prefix/$ac_dir"
    break
    fi
  done
for ac_dir in $ac_cctbx_dirs; do
  if test -r "$cctbx_prefix/$ac_dir/boost/type_traits.hpp"; then
    ac_CCTBX_CXXFLAGS="$ac_CCTBX_CXXFLAGS -I$cctbx_prefix/$ac_dir"
    break
    fi
  done
for ac_dir in $ac_cctbx_dirs; do
  if test -r "$cctbx_prefix/$ac_dir/tnt_array1d.h"; then
    ac_CCTBX_CXXFLAGS="$ac_CCTBX_CXXFLAGS -I$cctbx_prefix/$ac_dir"
    break
    fi
  done
else
 ac_CCTBX_CXXFLAGS=""
 ac_CCTBX_LDOPTS="-lcctbx"
fi

_AM_PATH_CCTBX_EXTRA
_AM_PATH_CCTBX_BUILD

AC_MSG_CHECKING([for CCTBX and BOOST])

LIBS="$ac_CCTBX_LDOPTS $saved_LIBS"
CXXFLAGS="$ac_CCTBX_CXXFLAGS $saved_CXXFLAGS"
#
# AC_TRY_LINK uses the c compiler (set by AC_LANG), so we will
# temporarily reassign $CC to the c++ compiler.
#
AC_LANG_PUSH(C++)
AC_TRY_LINK([#include "cctbx/miller.h"] ,[  cctbx::miller::index<int> a;  ], have_cctbx=yes, have_cctbx=no)
AC_LANG_POP(C++)  # the language we have just quit
AC_MSG_RESULT($have_cctbx)

 LIBS="$saved_LIBS"
 CXXFLAGS="$saved_CXXFLAGS"
]) # user override

AS_IF([test x$have_cctbx = xyes],
 [
   test "x$CCTBX_CXXFLAGS" = x && CCTBX_CXXFLAGS="$ac_CCTBX_CXXFLAGS"
   test "x$CCTBX_LIBS" = x && CCTBX_LIBS="$ac_CCTBX_LDOPTS"
   ifelse([$1], , :, [$1]) ],
 [
   ifelse([$2], , :, [$2]) ]
)
fi #dnl --with-cctbx

AC_SUBST(CCTBX_CXXFLAGS)
AC_SUBST(CCTBX_LIBS)
])
