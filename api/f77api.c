/* Copyright (c) 2007-2008 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

#include <stdlib.h>

#include "nlopt.h"
#include "nlopt-util.h"

/*-----------------------------------------------------------------------*/
/* wrappers around f77 procedures */

typedef void (*nlopt_f77_func)(double *val, const int *n, const double *x,
			       double *gradient, const int *need_gradient,
			       void *func_data);

typedef struct {
     nlopt_f77_func f;
     void *f_data;
} f77_func_data;

static double f77_func_wrap(int n, const double *x, double *grad, void *data)
{
     f77_func_data *d = (f77_func_data *) data;
     double val;
     int need_gradient = grad != 0;
     d->f(&val, &n, x, grad, &need_gradient, d->f_data);
     return val;
}

/*-----------------------------------------------------------------------*/
/* rather than trying to detect the Fortran name-mangling scheme with
   autoconf, we just include wrappers with all common name-mangling
   schemes ... this avoids problems and also allows us to work with
   multiple Fortran compilers on the same machine . 

   Note that our Fortran function names do not contain underscores;
   otherwise, we would need to deal with the additional headache that
   g77 appends two underscores in that case. */

#ifndef WINDOWS_F77_MANGLING

/* name + underscore is by far the most common (gfortran, g77, Intel, ...) */
#  define F77(a, A) a ## _
#  include "f77funcs.h"

/* AIX and HPUX use just the lower-case name */
#  undef F77
#  define F77(a, A) a
#  include "f77funcs.h"

/* old Cray UNICOS used just the upper-case name */
#  undef F77
#  define F77(a, A) A
#  include "f77funcs.h"

#else /* WINDOWS_F77_MANGLING */

/* Various mangling conventions common (?) under Windows. */

/* name + underscore for gfortran, g77, ...? */
#  define F77(a, A) a ## _
#  include "f77funcs.h"

/* Digital/Compaq/HP Visual Fortran, Intel Fortran.  stdcall attribute
   is apparently required to adjust for calling conventions (callee
   pops stack in stdcall).  See also:
       http://msdn.microsoft.com/library/en-us/vccore98/html/_core_mixed.2d.language_programming.3a_.overview.asp
*/
#  undef F77
#  if defined(__GNUC__)
#    define F77(a, A) __attribute__((stdcall)) A
#  elif defined(_MSC_VER) || defined(_ICC) || defined(_STDCALL_SUPPORTED)
#    define F77(a, A) __stdcall A
#  else
#    define F77(a, A) A /* oh well */
#  endif
#  include "f77funcs.h"

#endif /* WINDOWS_F77_MANGLING */
