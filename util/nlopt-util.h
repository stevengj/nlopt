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

#ifndef NLOPT_UTIL_H
#define NLOPT_UTIL_H

#include <stdlib.h>
#include <math.h>
#include "config.h"

/* workaround for Solaris + gcc 3.4.x bug (see configure.ac) */
#if defined(__GNUC__) && defined(REPLACEMENT_HUGE_VAL)
#  undef HUGE_VAL
#  define HUGE_VAL REPLACEMENT_HUGE_VAL
#endif

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

int nlopt_isinf(double x);

/* re-entrant qsort */
extern void nlopt_qsort_r(void *base_, size_t nmemb, size_t size, void *thunk,
			  int (*compar)(void *, const void *, const void *));

/* seconds timer */
extern double nlopt_seconds(void);
extern unsigned long nlopt_time_seed(void);

/* pseudorandom number generation by Mersenne twister algorithm */
extern void nlopt_init_genrand(unsigned long s);
extern double nlopt_urand(double a, double b);
extern int nlopt_iurand(int n);

/* Sobol' low-discrepancy-sequence generation */
typedef struct nlopt_soboldata_s *nlopt_sobol;
extern nlopt_sobol nlopt_sobol_create(unsigned sdim);
extern void nlopt_sobol_destroy(nlopt_sobol s);
extern void nlopt_sobol_next01(nlopt_sobol s, double *x);
extern void nlopt_sobol_next(nlopt_sobol s, double *x,
			    const double *lb, const double *ub);
extern void nlopt_sobol_skip(nlopt_sobol s, unsigned n, double *x);

/* stopping criteria */
typedef struct {
     int n;
     double minf_max;
     double ftol_rel;
     double ftol_abs;
     double xtol_rel;
     const double *xtol_abs;
     int nevals, maxeval;
     double maxtime, start;
} nlopt_stopping;
extern int nlopt_stop_f(const nlopt_stopping *stop, double f, double oldf);
extern int nlopt_stop_ftol(const nlopt_stopping *stop, double f, double oldf);
extern int nlopt_stop_x(const nlopt_stopping *stop, 
			const double *x, const double *oldx);
extern int nlopt_stop_xs(const nlopt_stopping *stop, 
			 const double *xs, const double *oldxs,
			 const double *scale_min, const double *scale_max);
extern int nlopt_stop_evals(const nlopt_stopping *stop);
extern int nlopt_stop_time(const nlopt_stopping *stop);
extern int nlopt_stop_evalstime(const nlopt_stopping *stop);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
