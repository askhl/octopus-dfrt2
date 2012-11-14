# Check whether zdotc routine in BLAS works properly.
# The test program here may either give nonzero results or crash.

AC_DEFUN([ACX_ZDOTC], [

dnl may be needed for cross-compiling, when compilation and compute nodes are different in architecture (yet report the same for host and build)
  acx_enable_zdotc_test=yes
  AC_ARG_ENABLE(zdotc-test, [AS_HELP_STRING([--disable-zdotc-test], [Assume zdotc works and do not perform a test.])], [acx_enable_zdotc_test=${enableval}])

  if test "x$acx_enable_zdotc_test" = "xyes"; then

  AC_MSG_CHECKING(whether zdotc works)
  AC_LANG_ASSERT(Fortran)

  acx_blas_save_LIBS="$LIBS"
  LIBS="$LIBS $LIBS_BLAS"

  acx_zdotc_ok=yes

  AC_RUN_IFELSE([AC_LANG_PROGRAM([],[
complex*16, allocatable :: f1(:), f2(:)
complex*16 :: result1, result2
complex*16, external :: zdotc
integer :: nn, ii

nn=100

allocate(f1(nn), f2(nn))

result1 = cmplx(0.0d0,0.0d0)
do ii=1, nn
  f1(ii) = cmplx(exp(dble(ii)/dble(nn)), 0.0d0)
  f2(ii) = cmplx(cos(dble(ii)/dble(nn)), 0.0d0)
  result1 = result1 + exp(dble(ii)/dble(nn))*cos(dble(ii)/dble(nn))
end do

result2 = cmplx(0.0d0,0.0d0)
result2 = zdotc(nn,f1,1,f2,1)

open(1, file='conftest.out')
if(abs(result1-result2) .lt. 1d-6) then
  write(1, '(a)') 'success'
endif
])], [  
  if test "x$acx_zdotc_ok" = "xyes"; then
    if test "x`cat conftest.out`" != "xsuccess"; then
      acx_zdotc_ok=no
      # program didn't crash, but gave wrong answer
    fi
  fi
], [acx_zdotc_ok=no], [acx_zdotc_ok=yes;echo -n "cross-compiling; assuming... "])


  AC_MSG_RESULT([$acx_zdotc_ok])
  LIBS="$acx_blas_save_LIBS"

  if test "x$acx_zdotc_ok" = "xno"; then
    AC_MSG_ERROR([

    Blas was found, but the zdotc subroutine does not work
    properly. This probably means that your Blas implementation is
    incompatible with your Fortran compiler.

    ])
  fi

  fi
])
