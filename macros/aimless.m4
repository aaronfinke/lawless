# AC_AIMLESS_OPTIONS      
# ------------------
# target compilation optons
#
AC_DEFUN([AC_AIMLESS_OPTIONS],
[
#specifics for various machines
test "${target_os:+set}" = set || target_os="$host_os"
case "$target_os" in
  *osf* | *64* )
    if test "`basename $CXX`" = cxx; then
      case "$CXXFLAGS" in
        *strict_ansi* ) ;;
        * )
          CXXFLAGS="$CXXFLAGS -ieee -std strict_ansi -alternative_tokens -timplicit_local -no_implicit_include"
        esac
    fi
  ;;
  *linux* )
    if test "`basename $CXX`" = icpc ; then
      case `$CXX -dumpversion 2>&1` in
        11.0 )
          case "$CXXFLAGS" in
            *-fno-inline*) ;;
            *) CXXFLAGS="$CXXFLAGS -fno-inline"
          esac
        ;;
        *)
      esac
    fi
   ;;
  *irix* )
    if test "`basename $CXX`" = CC; then
      case "$CXXFLAGS" in
        *LANG:std* ) ;;
        * )
          CXXFLAGS="$CXXFLAGS -LANG:std"
        esac
    fi
  ;;
  *darwin* )
    if test "`basename $CXX`" = gcc || test "`basename $CXX`" = g++; then
      case `$CXX -v 2>&1` in 
       *3.1*)
# problem with PIC relocation tables for 3.1
        case "$CXXFLAGS" in
         *-O* | *-O1* | *-O2* | *-O3* )
          CXXFLAGS=`echo $CXXFLAGS | sed s%-O[[\ 123]]%-O0%g`
          ;;
         *-O0* ) ;;
         * )
          CXXFLAGS="$CXXFLAGS -O0"
        esac
        ;;
      *4.0*)
        case "$CXXFLAGS" in
          *-fno-inline* ) ;;
          *) CXXFLAGS="$CXXFLAGS -fno-inline"
          esac
        ;;
      esac
    fi
  ;;
  * )
esac

] )

