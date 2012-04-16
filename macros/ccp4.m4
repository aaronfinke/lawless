
# AM_PATH_CCP4([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
AC_DEFUN([AM_PATH_CCP4],
[
AC_PROVIDE([AM_PATH_CCP4])

AC_ARG_WITH(ccp4,
  AC_HELP_STRING( [--with-ccp4=PFX], [location of ccp4c] ),
  [
    test "$withval" = no || with_ccp4=yes 
    test "$withval" = yes || ccp4_prefix="$withval" ],
  [ with_ccp4=yes ] ) 

if test x$with_ccp4 = xyes ; then
#user override
AS_IF([test "x$CCP4_LIBS" != x && test "x$CCP4_CXXFLAGS" != x ],
[
  have_ccp4=yes
],
[
AC_MSG_CHECKING([for ccp4_errno in CCP4])

saved_LIBS="$LIBS"
saved_CXXFLAGS="$CXXFLAGS"
CCP4_LIBS=""
CCP4_CXXFLAGS=""

if test "x$ccp4_prefix" != x; then
ac_ccp4_dirs='
.
include
lib'
for ac_dir in $ac_ccp4_dirs; do
  if test -r "$ccp4_prefix/$ac_dir/ccp4/ccp4_errno.h"; then
    ac_CCP4_CXXFLAGS="-I$ccp4_prefix/$ac_dir"
    break
    fi
  done
for ac_dir in $ac_ccp4_dirs; do
  for ac_extension in a so sl dylib; do
  if test -r "$ccp4_prefix/$ac_dir/libccp4c.$ac_extension"; then
    ac_CCP4_LDOPTS="-L$ccp4_prefix/$ac_dir -lccp4f -lccp4c"
    break 2
    fi
  done
  done
for ac_dir in $ac_ccp4_dirs; do
  for ac_extension in a so sl dylib; do
  if test -r "$ccp4_prefix/$ac_dir/libccif.$ac_extension"; then
    ac_CCP4_LDOPTS="$ac_CCP4_LDOPTS -L$ccp4_prefix/$ac_dir -lccif"
    break 2
    fi
  done
  done
for ac_dir in $ac_ccp4_dirs; do
  for ac_extension in a so sl dylib; do
  if test -r "$ccp4_prefix/$ac_dir/libmmdb.$ac_extension"; then
    ac_CCP4_LDOPTS="$ac_CCP4_LDOPTS -L$ccp4_prefix/$ac_dir -lmmdb"
    break 2
    fi
  done
  done
else
 ac_CCP4_CXXFLAGS=""
 ac_CCP4_LDOPTS="-lccp4f -lccp4c -lccif -lmmdb"
fi


LIBS="$ac_CCP4_LDOPTS $saved_LIBS $FCLIBS"
CXXFLAGS="$ac_CCP4_CXXFLAGS $saved_CXXFLAGS"
#
# AC_TRY_LINK uses the c compiler (set by AC_LANG), so we will
# temporarily reassign $CC to the c++ compiler.
#
AC_LANG_PUSH(C++)
AC_TRY_LINK([#include "ccp4/ccp4_errno.h"],
  [int a = ccp4_errno;  CCP4::ccp4_error("conftest"); ], have_ccp4=yes, have_ccp4=no)
AC_LANG_POP(C++)  # the language we have just quit
AC_MSG_RESULT($have_ccp4)

 LIBS="$saved_LIBS"
 CXXFLAGS="$saved_CXXFLAGS"

]) # user override

AS_IF([test x$have_ccp4 = xyes],
 [
   test "x$CCP4_CXXFLAGS" = x && CCP4_CXXFLAGS="$ac_CCP4_CXXFLAGS"
   test "x$CCP4_LIBS" = x && CCP4_LIBS="$ac_CCP4_LDOPTS"
   ifelse([$1], , :, [$1]) ],
 [
   ifelse([$2], , :, [$2]) ]
)

fi #dnl --with-ccp4

AC_SUBST(CCP4_CXXFLAGS)
AC_SUBST(CCP4_LIBS)
])

