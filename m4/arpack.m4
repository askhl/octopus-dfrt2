dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/acx_arpack.html
dnl
AC_DEFUN([ACX_ARPACK], [
AC_REQUIRE([ACX_BLAS])
acx_arpack_ok=no

dnl We cannot use ARPACK if BLAS is not found
if test "x$acx_blas_ok" != xyes; then
  acx_arpack_ok=noblas
fi

dnl Get fortran linker name of ARPACK function to check for.
dnl if not compiling with fortran, convert the names
m4_if(_AC_LANG, Fortran, [dnaupd=dnaupd], [AC_F77_FUNC(dnaupd)])

dnl Check if the library was given in the command line
if test $acx_arpack_ok = no; then
  AC_ARG_WITH(arpack, [AS_HELP_STRING([--with-arpack=<lib>], [use ARPACK library <lib> http://forge.scilab.org/index.php/p/arpack-ng/])])
  case $with_arpack in
    yes | "") ;;
    no) acx_arpack_ok=disable ;;
    -* | */* | *.a | *.so | *.so.* | *.o) LIBS_ARPACK="$with_arpack" ;;
    *) LIBS_ARPACK="-l$with_arpack" ;;
  esac
fi

dnl Backup LIBS 
acx_arpack_save_LIBS="$LIBS"
LIBS="$LIBS_ARPACK $LIBS_BLAS $LIBS $FLIBS"

dnl First, check LIBS_ARPACK environment variable
if test $acx_arpack_ok = no; then
  AC_MSG_CHECKING([for $dnaupd in $LIBS_ARPACK])
  AC_LINK_IFELSE($dnaupd, [acx_arpack_ok=yes], [])
  if test $acx_arpack_ok = no; then
    AC_MSG_RESULT([$acx_arpack_ok])
  else
    AC_MSG_RESULT([$acx_arpack_ok ($LIBS_ARPACK)])
  fi
fi

dnl Generic ARPACK library?
for arpack in arpack parpack; do
  if test $acx_arpack_ok = no; then
    AC_CHECK_LIB($arpack, $dnaupd,
      [acx_arpack_ok=yes; LIBS_ARPACK="$LIBS_ARPACK -l$arpack"], [], [$FLIBS])
  fi
done

AC_SUBST(LIBS_ARPACK)
LIBS="$acx_arpack_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_arpack_ok" = xyes; then
  AC_DEFINE(HAVE_ARPACK,1,[Defined if you have ARPACK library.])
  $1
else
    AC_MSG_WARN([Could not find ARPACK library. 
               *** Will compile without ARPACK support])
  $2
fi
])dnl ACX_ARPACK
