dnl Process this file with autoconf to produce a configure script.
AC_INIT(nlopt, 0.1, stevenj@alum.mit.edu)
AC_CONFIG_SRCDIR(api/nlopt.h)
SHARED_VERSION_INFO="0:0:0" # CURRENT:REVISION:AGE

AM_INIT_AUTOMAKE(1.7)
AM_CONFIG_HEADER(config.h)
AM_MAINTAINER_MODE
AC_SUBST(SHARED_VERSION_INFO)
AC_DISABLE_SHARED dnl to hell with shared libraries

dnl Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_CC_STDC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_LIBTOOL_WIN32_DLL
AC_PROG_LIBTOOL

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
AC_CHECK_FUNCS([BSDgettimeofday gettimeofday time])

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
dnl Compiling Octave plug-in

AC_ARG_VAR(OCT_INSTALL_DIR, [where to install GNU Octave .oct plug-ins])
AC_ARG_VAR(M_INSTALL_DIR, [where to install GNU Octave .m plug-ins])
AC_ARG_VAR(MKOCTFILE, [name of mkoctfile program to compile Octave plug-ins])

AC_CHECK_PROGS(MKOCTFILE, mkoctfile, echo)
if test "$MKOCTFILE" = "echo"; then
	AC_MSG_WARN([can't find mkoctfile: won't be able to compile GNU Octave plugin])
elif test x"$OCT_INSTALL_DIR" = "x"; then
	# try to find installation directory
	AC_CHECK_PROGS(OCTAVE, octave, echo)
	AC_MSG_CHECKING(for Octave loadpath)
	OCTAVE_LOADPATH=`echo "DEFAULT_LOADPATH" | $OCTAVE -q | cut -d'=' -f2`
	AC_MSG_RESULT($OCTAVE_LOADPATH)
	AC_MSG_CHECKING(where Octave plugins go)
	OCT_INSTALL_DIR=`echo "$OCTAVE_LOADPATH" | tr ':' '\n' | grep "site/oct" | head -1`
	if test -n "$OCT_INSTALL_DIR"; then
		AC_MSG_RESULT($OCT_INSTALL_DIR)
	else
		AC_MSG_RESULT(unknown)
		AC_MSG_WARN([can't find where to install GNU Octave plugins])
	fi
	AC_MSG_CHECKING(where Octave scripts go)
	M_INSTALL_DIR=`echo "$OCTAVE_LOADPATH" | tr ':' '\n' | grep "site/m" | head -1`
	if test -n "$M_INSTALL_DIR"; then
		AC_MSG_RESULT($M_INSTALL_DIR)
	else
		AC_MSG_RESULT(unknown)
		AC_MSG_WARN([can't find where to install GNU Octave scripts])
	fi
fi
AM_CONDITIONAL(WITH_OCTAVE, test x"$OCT_INSTALL_DIR" != "x")
AC_SUBST(OCT_INSTALL_DIR)
AC_SUBST(M_INSTALL_DIR)
AC_SUBST(MKOCTFILE)

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
   lbfgs/Makefile
   subplex/Makefile
   test/Makefile
])

AC_OUTPUT
