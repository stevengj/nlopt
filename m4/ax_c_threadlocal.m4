dnl @synopsis AX_C_THREADLOCAL([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl @summary determine C keyword for threadlocal storage
dnl
dnl This macro tries to discover a C keyword to declare variables
dnl as having thread-local storage.  Most commonly, this is either
dnl __thread [gcc] or __declspec(thread) [Windows].
dnl
dnl On success, it #defines the THREADLOCAL preprocessor symbol to
dnl the appropriate keyword.  You would then use it in C code as, e.g.:
dnl     THREADLOCAL int myvariable;
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if an thread-local
dnl keyword is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if one is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action does nothing.
dnl
dnl @version 2010-05-28
dnl @license GPLWithACException
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
AC_DEFUN([AX_C_THREADLOCAL],
[AC_ARG_WITH(threadlocal,
  [AC_HELP_STRING([--without-threadlocal], [no thread-local storage keyword])],
  with_ax_c_threadlocal=$withval, with_ax_c_threadlocal=yes)
 AC_CACHE_CHECK([for C thread-local keyword], ax_cv_c_threadlocal,
[if test "x$with_ax_c_threadlocal" = xno; then
   ax_cv_c_threadlocal=disabled
 else
   ax_cv_c_threadlocal=unsupported
   AC_LANG_SAVE
   AC_LANG_C
   for ax_kw in __thread "__declspec(thread)"; do
     AC_TRY_COMPILE([], [static $ax_kw int x = 0;], 
                        [ax_cv_c_threadlocal=$ax_kw; break])
   done
   AC_LANG_RESTORE
 fi
])
 ax_kw="$ax_cv_c_threadlocal"
 if test "x$ax_kw" = xunsupported; then ax_kw=""; fi
 if test "x$ax_kw" = xdisabled; then ax_kw=""; fi
 AC_DEFINE_UNQUOTED(THREADLOCAL, $ax_kw, [Define to C thread-local keyword, or to nothing if this is not supported in your compiler.])
 if test "$ax_cv_c_threadlocal" = unsupported; then
   m4_default([$2],:)
 else
   m4_default([$1],:)
 fi
])
