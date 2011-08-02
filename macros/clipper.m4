
# AM_PATH_CLIPPER([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
AC_DEFUN([AM_PATH_CLIPPER],
[
AC_PROVIDE([AM_PATH_CLIPPER])

AC_ARG_WITH(clipper,
  AC_HELP_STRING( [--with-clipper=PFX], [location of clipper] ),
  [
    test "$withval" = no || with_clipper=yes 
    test "$withval" = yes || clipper_prefix="$withval" ],
  [ with_clipper=yes ] ) 

if test x$with_clipper = xyes ; then
#user override
AS_IF([test "x$CLIPPER_LIBS" != x && test "x$CLIPPER_CXXFLAGS" != x ],
[
  have_clipper=yes
],
[
AC_MSG_CHECKING([for clipper_errno in CLIPPER])

saved_LIBS="$LIBS"
saved_CXXFLAGS="$CXXFLAGS"
CLIPPER_LIBS=""
CLIPPER_CXXFLAGS=""

if test "x$clipper_prefix" != x; then
ac_clipper_dirs='
.
include
include/clipper
lib/clipper/clipper/
lib'
for ac_dir in $ac_clipper_dirs; do
  if test -r "$clipper_prefix/$ac_dir/clipper/clipper.h"; then
    ac_CLIPPER_CXXFLAGS="-I$clipper_prefix/$ac_dir"
    break
    fi
  done
for ac_dir in $ac_clipper_dirs; do
  for ac_extension in a so sl dylib; do
  if test -r "$clipper_prefix/$ac_dir/libclipper-core.$ac_extension"; then
    ac_CLIPPER_LDOPTS="-L$clipper_prefix/$ac_dir -lclipper-ccp4 -lclipper-contrib -lclipper-minimol -lclipper-mmdb -lclipper-core"
    break 2
    fi
  done
  done
else
 ac_CLIPPER_CXXFLAGS=""
 ac_CLIPPER_LDOPTS="-lclipper-ccp4 -lclipper-contrib -lclipper-minimol -lclipper-mmdb -lclipper-core"
fi


LIBS="$ac_CLIPPER_LDOPTS $saved_LIBS"
CXXFLAGS="$ac_CLIPPER_CXXFLAGS $saved_CXXFLAGS"
#
# AC_TRY_LINK uses the c compiler (set by AC_LANG), so we will
# temporarily reassign $CC to the c++ compiler.
#
AC_LANG_PUSH(C++)
AC_TRY_LINK([#include "clipper/clipper.h"],
  [clipper::String a("test");  ], have_clipper=yes, have_clipper=no)
AC_LANG_POP(C++)  # the language we have just quit
AC_MSG_RESULT($have_clipper)

 LIBS="$saved_LIBS"
 CXXFLAGS="$saved_CXXFLAGS"

]) # user override

AS_IF([test x$have_clipper = xyes],
 [
   test "x$CLIPPER_CXXFLAGS" = x && CLIPPER_CXXFLAGS="$ac_CLIPPER_CXXFLAGS"
   test "x$CLIPPER_LIBS" = x && CLIPPER_LIBS="$ac_CLIPPER_LDOPTS"
   ifelse([$1], , :, [$1]) ],
 [
   ifelse([$2], , :, [$2]) ]
)

fi #dnl --with-clipper

AC_SUBST(CLIPPER_CXXFLAGS)
AC_SUBST(CLIPPER_LIBS)
])

