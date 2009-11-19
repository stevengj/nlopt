/* Copyright (c) 2007-2009 Massachusetts Institute of Technology
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

#ifndef NLOPT_H
#define NLOPT_H

#include <stddef.h> /* for ptrdiff_t */

/* for Windows compilers, you should add a line
           #define NLOPT_DLL
   when using NLopt from a DLL, in order to do the proper
   Windows importing nonsense. */
#if defined(NLOPT_DLL) && (defined(_WIN32) || defined(__WIN32__))
/* annoying Windows syntax for calling functions in a DLL */
#  define NLOPT_EXTERN extern __declspec(dllimport)
#else
#  define NLOPT_EXTERN extern
#endif

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

typedef double (*nlopt_func)(int n, const double *x,
			     double *gradient, /* NULL if not needed */
			     void *func_data);

typedef enum {
     /* Naming conventions:

        NLOPT_{G/L}{D/N}_* 
	    = global/local derivative/no-derivative optimization, 
              respectively 
 
	*_RAND algorithms involve some randomization.

	*_NOSCAL algorithms are *not* scaled to a unit hypercube
	         (i.e. they are sensitive to the units of x)
	*/

     NLOPT_GN_DIRECT = 0,
     NLOPT_GN_DIRECT_L,
     NLOPT_GN_DIRECT_L_RAND,
     NLOPT_GN_DIRECT_NOSCAL,
     NLOPT_GN_DIRECT_L_NOSCAL,
     NLOPT_GN_DIRECT_L_RAND_NOSCAL,

     NLOPT_GN_ORIG_DIRECT,
     NLOPT_GN_ORIG_DIRECT_L,

     NLOPT_GD_STOGO,
     NLOPT_GD_STOGO_RAND,

     NLOPT_LD_LBFGS_NOCEDAL,

     NLOPT_LD_LBFGS,

     NLOPT_LN_PRAXIS,

     NLOPT_LD_VAR1,
     NLOPT_LD_VAR2,

     NLOPT_LD_TNEWTON,
     NLOPT_LD_TNEWTON_RESTART,
     NLOPT_LD_TNEWTON_PRECOND,
     NLOPT_LD_TNEWTON_PRECOND_RESTART,

     NLOPT_GN_CRS2_LM,

     NLOPT_GN_MLSL,
     NLOPT_GD_MLSL,
     NLOPT_GN_MLSL_LDS,
     NLOPT_GD_MLSL_LDS,

     NLOPT_LD_MMA,

     NLOPT_LN_COBYLA,

     NLOPT_LN_NEWUOA,
     NLOPT_LN_NEWUOA_BOUND,

     NLOPT_LN_NELDERMEAD,
     NLOPT_LN_SBPLX,

     NLOPT_LN_AUGLAG,
     NLOPT_LD_AUGLAG,
     NLOPT_LN_AUGLAG_EQ,
     NLOPT_LD_AUGLAG_EQ,

     NLOPT_LN_BOBYQA,

     NLOPT_GN_ISRES,

     NLOPT_NUM_ALGORITHMS /* not an algorithm, just the number of them */
} nlopt_algorithm;

extern const char *nlopt_algorithm_name(nlopt_algorithm a);

typedef enum {
     NLOPT_FAILURE = -1, /* generic failure code */
     NLOPT_INVALID_ARGS = -2,
     NLOPT_OUT_OF_MEMORY = -3,
     NLOPT_ROUNDOFF_LIMITED = -4,

     NLOPT_SUCCESS = 1, /* generic success code */
     NLOPT_MINF_MAX_REACHED = 2,
     NLOPT_FTOL_REACHED = 3,
     NLOPT_XTOL_REACHED = 4,
     NLOPT_MAXEVAL_REACHED = 5,
     NLOPT_MAXTIME_REACHED = 6
} nlopt_result;

NLOPT_EXTERN nlopt_result nlopt_minimize(
     nlopt_algorithm algorithm,
     int n, nlopt_func f, void *f_data,
     const double *lb, const double *ub, /* bounds */
     double *x, /* in: initial guess, out: minimizer */
     double *minf, /* out: minimum */
     double minf_max, double ftol_rel, double ftol_abs,
     double xtol_rel, const double *xtol_abs,
     int maxeval, double maxtime);

NLOPT_EXTERN nlopt_result nlopt_minimize_constrained(
     nlopt_algorithm algorithm,
     int n, nlopt_func f, void *f_data,
     int m, nlopt_func fc, void *fc_data, ptrdiff_t fc_datum_size,
     const double *lb, const double *ub, /* bounds */
     double *x, /* in: initial guess, out: minimizer */
     double *minf, /* out: minimum */
     double minf_max, double ftol_rel, double ftol_abs,
     double xtol_rel, const double *xtol_abs,
     int maxeval, double maxtime);

NLOPT_EXTERN nlopt_result nlopt_minimize_econstrained(
     nlopt_algorithm algorithm,
     int n, nlopt_func f, void *f_data,
     int m, nlopt_func fc, void *fc_data, ptrdiff_t fc_datum_size,
     int p, nlopt_func h, void *h_data, ptrdiff_t h_datum_size,
     const double *lb, const double *ub, /* bounds */
     double *x, /* in: initial guess, out: minimizer */
     double *minf, /* out: minimum */
     double minf_max, double ftol_rel, double ftol_abs,
     double xtol_rel, const double *xtol_abs,
     double htol_rel, double htol_abs,
     int maxeval, double maxtime);

NLOPT_EXTERN void nlopt_srand(unsigned long seed);
NLOPT_EXTERN void nlopt_srand_time(void);

NLOPT_EXTERN void nlopt_version(int *major, int *minor, int *bugfix);

NLOPT_EXTERN void nlopt_get_local_search_algorithm(nlopt_algorithm *deriv,
					     nlopt_algorithm *nonderiv,
					     int *maxeval);
NLOPT_EXTERN void nlopt_set_local_search_algorithm(nlopt_algorithm deriv,
					     nlopt_algorithm nonderiv,
					     int maxeval);

NLOPT_EXTERN int nlopt_get_stochastic_population(void);
NLOPT_EXTERN void nlopt_set_stochastic_population(int pop);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
