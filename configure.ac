dnl Process this file with autoconf to produce a configure script.
AC_INIT(nlopt, 1.2, stevenj@alum.mit.edu)
AC_CONFIG_SRCDIR(api/nlopt.h)
SHARED_VERSION_INFO="2:0:2" # CURRENT:REVISION:AGE

AM_INIT_AUTOMAKE(1.7)
AM_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([m4])
AM_MAINTAINER_MODE
AC_SUBST(SHARED_VERSION_INFO)
AC_DISABLE_SHARED dnl to hell with shared libraries

dnl Checks for programs.
AC_PROG_CC
AC_PROG_CC_STDC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_LIBTOOL_WIN32_DLL
AC_PROG_LIBTOOL

NLOPT_SUFFIX=""
AC_ARG_WITH(cxx,
	[AC_HELP_STRING([--with-cxx], [include C++-based routines])],
        with_cxx=$withval,with_cxx=no)
AM_CONDITIONAL(WITH_CXX, test "x$with_cxx" = xyes)
AC_PROG_CXX
if test "x$with_cxx" = xyes; then
   AC_DEFINE([WITH_CXX], 1, [Define if compiled including C++-based routines])
   NLOPT_SUFFIX="_cxx"
fi
AC_SUBST(NLOPT_SUFFIX)

AC_ARG_WITH(windows-f77-mangling, [AC_HELP_STRING([--with-windows-f77-mangling],[use common Win32 Fortran interface styles])], with_windows_f77_mangling=$withval, with_windows_f77_mangling=no)
if test "$with_windows_f77_mangling" = "yes"; then
	AC_DEFINE(WINDOWS_F77_MANGLING,1,[Use common Windows Fortran mangling styles for the Fortran interfaces.])
fi

dnl Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDC
AC_HEADER_TIME
AC_CHECK_HEADERS([unistd.h getopt.h stdint.h])
AC_C_CONST
AC_C_INLINE

dnl find 32-bit unsigned integer type for random-number generator
AC_CHECK_SIZEOF(unsigned int)
AC_CHECK_SIZEOF(unsigned long)
AC_CHECK_TYPES(uint32_t, [], [], [$ac_includes_default
#ifdef HAVE_STDINT_H
#  include <stdint.h>
#endif])

dnl Checks for libraries and functions
AC_CHECK_LIB(m, sin)
AC_CHECK_FUNCS([BSDgettimeofday gettimeofday time qsort_r])

AC_MSG_CHECKING([for isnan])
AC_TRY_LINK([#include <math.h>
], if (!isnan(3.14159)) isnan(2.7183);, ok=yes, ok=no)
if test "$ok" = "yes"; then
	AC_DEFINE(HAVE_ISNAN,1,[Define if the isnan() function/macro is available.])
fi
AC_MSG_RESULT(${ok})

AC_MSG_CHECKING([for isinf])
AC_TRY_LINK([#include <math.h>
], if (!isinf(3.14159)) isinf(2.7183);, ok=yes, ok=no)
if test "$ok" = "yes"; then
	AC_DEFINE(HAVE_ISINF,1,[Define if the isinf() function/macro is available.])
fi
AC_MSG_RESULT(${ok})

dnl -----------------------------------------------------------------------

if test "x$with_cxx" = xyes; then
  test -r $srcdir/lbfgs/ap.cpp && test -r $srcdir/lbfgs/ap.h && test -r $srcdir/lbfgs/l-bfgs-b.cpp && test -r $srcdir/lbfgs/l-bfgs-b.h && have_lbfgs=yes
fi
AM_CONDITIONAL(WITH_NOCEDAL, test -n "$have_lbfgs")
if test -n "$have_lbfgs"; then
   AC_DEFINE(WITH_NOCEDAL, [1], [Define if we have the non-free Nocedal LBFGS code])
fi

dnl -----------------------------------------------------------------------
dnl Compiling Octave plug-in

AC_ARG_VAR(OCT_INSTALL_DIR, [where to install GNU Octave .oct plug-ins])
AC_ARG_VAR(M_INSTALL_DIR, [where to install GNU Octave .m plug-ins])
AC_ARG_VAR(MKOCTFILE, [name of mkoctfile program to compile Octave plug-ins])

AC_ARG_WITH(octave,
	[AC_HELP_STRING([--without-octave], [don't compile Octave plugin])],
        with_octave=$withval,with_octave=yes)

AC_CHECK_PROGS(MKOCTFILE, mkoctfile, echo)

if test x"$with_octave" = xno; then
	OCT_INSTALL_DIR=""
elif test "$MKOCTFILE" = "echo"; then
	AC_MSG_WARN([can't find mkoctfile: won't be able to compile GNU Octave plugin])
	OCT_INSTALL_DIR=""
elif test x"$OCT_INSTALL_DIR" = "x"; then
	# try to find installation directory
	AC_CHECK_PROGS(OCTAVE, octave, echo)
	AC_CHECK_PROGS(OCTAVE_CONFIG, octave-config, echo)

	AC_MSG_CHECKING(where octave plugins go)
	OCT_INSTALL_DIR=`octave-config --oct-site-dir 2> /dev/null | grep '/'`
	if test -z "$OCT_INSTALL_DIR"; then
		OCT_INSTALL_DIR=`octave-config --print OCTFILEDIR 2> /dev/null | grep '/'`
	fi
	if test -z "$OCT_INSTALL_DIR"; then
		OCT_INSTALL_DIR=`echo "path" | $OCTAVE -q 2> /dev/null | grep "/oct/" | head -1`
	fi
	if test -z "$OCT_INSTALL_DIR"; then
		OCT_INSTALL_DIR=`echo "DEFAULT_LOADPATH" | $OCTAVE -q 2> /dev/null | tr ':' '\n' | grep "site/oct" | head -1`
        fi
	if test -n "$OCT_INSTALL_DIR"; then
		AC_MSG_RESULT($OCT_INSTALL_DIR)
	else
		AC_MSG_RESULT(unknown)
		AC_MSG_WARN([can't find where to install octave plugins: won't be able to compile octave plugin])
	fi

	AC_MSG_CHECKING(where octave scripts go)
	M_INSTALL_DIR=`octave-config --m-site-dir 2> /dev/null | grep '/'`
	if test -z "$M_INSTALL_DIR"; then
		M_INSTALL_DIR=`octave-config --print FCNFILEDIR 2> /dev/null | grep '/'`
	fi
	if test -z "$M_INSTALL_DIR"; then
		M_INSTALL_DIR=`echo "path" | $OCTAVE -q 2> /dev/null | grep "/m" | head -1`
	fi
	if test -z "$M_INSTALL_DIR"; then
		M_INSTALL_DIR=`echo "DEFAULT_LOADPATH" | $OCTAVE -q 2> /dev/null | tr ':' '\n' | grep "site/m" | head -1`
        fi
	if test -n "$M_INSTALL_DIR"; then
		AC_MSG_RESULT($M_INSTALL_DIR)
	else
		AC_MSG_RESULT(unknown)
		AC_MSG_WARN([can't find where to install octave scripts: won't be able to install octave plugin])
		OCT_INSTALL_DIR=""
	fi
fi

if test x"$OCT_INSTALL_DIR" != "x"; then
if test "$enable_shared" = no; then
     AC_MSG_CHECKING([whether mkoctfile can link with static objects])
     rm -f conftest*
     echo 'void foo(void) {};' > conftestfoo.c
     $CC -c conftestfoo.c >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD
     cat > conftest.cc <<EOF
#include <octave/oct.h>
extern "C" void foo(void);
DEFUN_DLD(foobar, args, nargout, "foo bar") { foo(); }
EOF
     if $MKOCTFILE conftest.cc conftestfoo.o >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD; then
	AC_MSG_RESULT(yes)
     else
	AC_MSG_RESULT(no)
        AC_MSG_WARN([mkoctfile requires --enable-shared; won't compile Octave plugin])
	OCT_INSTALL_DIR=""
     fi
fi
fi

AM_CONDITIONAL(WITH_OCTAVE, test x"$OCT_INSTALL_DIR" != "x")
AC_SUBST(OCT_INSTALL_DIR)
AC_SUBST(M_INSTALL_DIR)
AC_SUBST(MKOCTFILE)

dnl -----------------------------------------------------------------------
dnl Compiling Matlab plug-in

AC_ARG_WITH(matlab,
	[AC_HELP_STRING([--without-matlab], [don't compile Matlab plugin])],
        with_matlab=$withval,with_matlab=yes)

AC_ARG_VAR(MEX_INSTALL_DIR, [where to install Matlab .mex plug-ins])
AC_ARG_VAR(MEX, [name of mex program to compile Matlab plug-ins])
AC_CHECK_PROGS(MEX, mex, echo)
if test x"$with_matlab" = xno; then
     MEX_INSTALL_DIR=""
elif test "$MEX" = "echo"; then
     AC_MSG_WARN([can't find mex: won't be able to compile Matlab plugin])
     MEX_INSTALL_DIR=""
else
     AC_MSG_CHECKING([for extension of compiled mex files])
     rm -f conftest*
     cat > conftest.c <<EOF
#include <mex.h>
void mexFunction(int nlhs, mxArray *plhs[[]],
                 int nrhs, const mxArray *prhs[[]]) { }
EOF
     if $MEX conftest.c >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD; then
        MEXSUFF=`ls conftest.m* | head -1 | cut -d'.' -f2`
	AC_MSG_RESULT($MEXSUFF)
	AC_CHECK_PROGS(MATLAB, matlab, echo)
     else
        AC_MSG_WARN([$MEX failed to compile a simple file; won't compile Matlab plugin])
	MEX_INSTALL_DIR=""
	MATLAB=echo
     fi

     if test x"$MATLAB" != xecho; then
     if test "$enable_shared" = no; then
     AC_MSG_CHECKING([whether mex can link with static objects])
     rm -f conftest*
     echo 'void foo(void) {};' > conftestfoo.c
     $CC -c conftestfoo.c >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD
     cat > conftest.c <<EOF
#include <mex.h>
extern void foo(void);
void mexFunction(int nlhs, mxArray *plhs[[]],
                 int nrhs, const mxArray *prhs[[]]) { foo(); }
EOF
     if $MEX conftest.c conftestfoo.o >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD; then
	AC_MSG_RESULT(yes)
     else
	AC_MSG_RESULT(no)
        AC_MSG_WARN([mex requires --enable-shared; won't compile Matlab plugin])
	MEX_INSTALL_DIR=""
	MATLAB=echo
     fi
     fi
     fi

     if test x"$MATLAB" != xecho; then
       # try to find installation directory
       if test x"$MEX_INSTALL_DIR" = "x"; then
          AC_MSG_CHECKING(for MATLAB mex installation dir)
	  MEX_INSTALL_DIR=`matlab -nodisplay -nodesktop -nojvm -r 'path;quit' | grep toolbox/local |sed 's,^[[^/]]*,,g' |sort |head -1`
	  AC_MSG_RESULT($MEX_INSTALL_DIR)
	  if test x`basename "$MEX_INSTALL_DIR"` != xlocal; then
	     MEX_INSTALL_DIR=""
	  fi
          if test x"$MEX_INSTALL_DIR" = "x"; then
            AC_MSG_WARN([can't find reasonable Matlab installation directory; Matlab plugins will not be compiled unless you manually specify MEX_INSTALL_DIR])
          fi
       fi
     else
	MEX_INSTALL_DIR=""
     fi
fi
AM_CONDITIONAL(WITH_MATLAB, test x"$MEX_INSTALL_DIR" != "x")
AC_SUBST(MEX_INSTALL_DIR)
AC_SUBST(MEX)
AC_SUBST(MEXSUFF)

dnl -----------------------------------------------------------------------
dnl Check for broken Solaris HUGE_VAL macro under gcc 3.4.x and similar

AC_MSG_CHECKING([for working HUGE_VAL])
AC_TRY_COMPILE([#include <math.h>], [double x = -HUGE_VAL;], 
[AC_MSG_RESULT([ok])],
[AC_TRY_COMPILE([#include <math.h>
#ifdef __GNUC__
#undef HUGE_VAL
#define HUGE_VAL __builtin_huge_val()
#endif], [double x = -HUGE_VAL;], 
[AC_MSG_RESULT([__builtin_huge_val()])
AC_DEFINE(REPLACEMENT_HUGE_VAL,[__builtin_huge_val()],
          [replacement for broken HUGE_VAL macro, if needed])],
[AC_MSG_RESULT([unknown])
AC_MSG_ERROR([broken HUGE_VAL macro with this compiler, unknown workaround])])])

dnl -----------------------------------------------------------------------
dnl Debugging

AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug],[compile with extra runtime checks for debugging])], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
	AC_DEFINE(DEBUG,1,[Define to enable extra debugging code.])
fi

dnl override CFLAGS selection when debugging
if test "${enable_debug}" = "yes"; then
        CFLAGS="-g"
        CXXFLAGS="-g"
        FFLAGS="-g"
fi

dnl add gcc warnings, in debug/maintainer mode only
if test "$enable_debug" = yes || test "$USE_MAINTAINER_MODE" = yes; then
if test $ac_cv_prog_gcc = yes; then
   if test "$ac_test_CFLAGS" != "set"; then
      CFLAGS="$CFLAGS -Wall -W -Wcast-qual -Wpointer-arith -Wcast-align -pedantic -Wshadow -Wbad-function-cast -Wwrite-strings -Wstrict-prototypes -Wredundant-decls -Wnested-externs -Wundef -Wconversion -Wmissing-prototypes -Wmissing-declarations"
   fi
   if test "$ac_test_CXXFLAGS" != "set"; then
      CXXFLAGS="$CXXFLAGS -Wall -W -Wcast-qual -Wpointer-arith -Wcast-align -pedantic -Wshadow -Wwrite-strings -Wredundant-decls -Wundef -Wconversion"
   fi
fi
fi

vers=`echo ${VERSION}.0.0 | cut -d. -f1`
AC_DEFINE_UNQUOTED(MAJOR_VERSION, $vers, [Major version number.])
vers=`echo ${VERSION}.0.0 | cut -d. -f2`
AC_DEFINE_UNQUOTED(MINOR_VERSION, $vers, [Minor version number.])
vers=`echo ${VERSION}.0.0 | cut -d. -f3`
AC_DEFINE_UNQUOTED(BUGFIX_VERSION, $vers, [Bugfix version number.])

dnl -----------------------------------------------------------------------

AC_CONFIG_FILES([
   Makefile
   nlopt.pc
   api/Makefile
   util/Makefile
   octave/Makefile
   direct/Makefile
   cdirect/Makefile
   stogo/Makefile
   praxis/Makefile
   lbfgs/Makefile
   luksan/Makefile
   crs/Makefile
   mlsl/Makefile
   mma/Makefile
   cobyla/Makefile
   newuoa/Makefile
   neldermead/Makefile
   auglag/Makefile
   bobyqa/Makefile
   isres/Makefile
   test/Makefile
])

AC_OUTPUT
