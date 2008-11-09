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
#include <math.h>
#include <string.h>
#include <float.h>

#include "nlopt.h"
#include "nlopt-util.h"
#include "config.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))

/*************************************************************************/

#ifdef INFINITY
#  define MY_INF INFINITY
#else
#  define MY_INF HUGE_VAL
#endif

int nlopt_isinf(double x) {
     return fabs(x) >= HUGE_VAL * 0.99
#ifdef HAVE_ISINF
	  || isinf(x)
#endif
	  ;
}

#ifndef HAVE_ISNAN
static int my_isnan(double x) { return x != x; }
#  define isnan my_isnan
#endif

/*************************************************************************/

void nlopt_version(int *major, int *minor, int *bugfix)
{
     *major = MAJOR_VERSION;
     *minor = MINOR_VERSION;
     *bugfix = BUGFIX_VERSION;
}

/*************************************************************************/

static const char nlopt_algorithm_names[NLOPT_NUM_ALGORITHMS][256] = {
     "DIRECT (global, no-derivative)",
     "DIRECT-L (global, no-derivative)",
     "Randomized DIRECT-L (global, no-derivative)",
     "Unscaled DIRECT (global, no-derivative)",
     "Unscaled DIRECT-L (global, no-derivative)",
     "Unscaled Randomized DIRECT-L (global, no-derivative)",
     "Original DIRECT version (global, no-derivative)",
     "Original DIRECT-L version (global, no-derivative)",
     "Subplex (local, no-derivative)",
#ifdef WITH_CXX
     "StoGO (global, derivative-based)",
     "StoGO with randomized search (global, derivative-based)",
#else
     "StoGO (NOT COMPILED)",
     "StoGO randomized (NOT COMPILED)",
#endif
#ifdef WITH_NOCEDAL_LBFGS
     "original NON-FREE L-BFGS code by Nocedal et al. (local, deriv.-based)",
#else
     "original NON-FREE L-BFGS code by Nocedal et al. (NOT COMPILED)",
#endif
     "Limited-memory BFGS (L-BFGS) (local, derivative-based)",
     "Principal-axis, praxis (local, no-derivative)",
     "Limited-memory variable-metric, rank 1 (local, derivative-based)",
     "Limited-memory variable-metric, rank 2 (local, derivative-based)",
     "Truncated Newton (local, derivative-based)",
     "Truncated Newton with restarting (local, derivative-based)",
     "Preconditioned truncated Newton (local, derivative-based)",
     "Preconditioned truncated Newton with restarting (local, derivative-based)",
     "Controlled random search (CRS2) with local mutation (global, no-derivative)",
     "Multi-level single-linkage (MLSL), random (global, no-derivative)",
     "Multi-level single-linkage (MLSL), random (global, derivative)",
     "Multi-level single-linkage (MLSL), quasi-random (global, no-derivative)",
     "Multi-level single-linkage (MLSL), quasi-random (global, derivative)",
     "Method of Moving Asymptotes (MMA) (local, derivative)",
     "COBYLA (Constrained Optimization BY Linear Approximations) (local, no-derivative)",
     "NEWUOA unconstrained optimization via quadratic models (local, no-derivative)",
     "Bound-constrained optimization via NEWUOA-based quadratic models (local, no-derivative)",
     "Nelder-Mead simplex algorithm"
};

const char *nlopt_algorithm_name(nlopt_algorithm a)
{
     if (a < 0 || a >= NLOPT_NUM_ALGORITHMS) return "UNKNOWN";
     return nlopt_algorithm_names[a];
}

/*************************************************************************/

static int nlopt_srand_called = 0;
void nlopt_srand(unsigned long seed) {
     nlopt_srand_called = 1;
     nlopt_init_genrand(seed);
}

void nlopt_srand_time() {
     nlopt_srand(nlopt_time_seed());
}

/*************************************************************************/

/* crude heuristics for initial step size of nonderivative algorithms */
static double initial_step(int n, 
			   const double *lb, const double *ub, const double *x)
{
     int i;
     double step = HUGE_VAL;
     for (i = 0; i < n; ++i) {
	  if (!nlopt_isinf(ub[i]) && !nlopt_isinf(lb[i])
	      && (ub[i] - lb[i]) * 0.25 < step && ub[i] > lb[i])
	       step = (ub[i] - lb[i]) * 0.25;
	  if (!nlopt_isinf(ub[i]) 
	      && ub[i] - x[i] < step && ub[i] > x[i])
	       step = ub[i] - x[i];
	  if (!nlopt_isinf(lb[i]) 
	      && x[i] - lb[i] < step && x[i] > lb[i])
	       step = x[i] - lb[i];
     }
     if (nlopt_isinf(step))
	  for (i = 0; i < n; ++i) {
	       if (!nlopt_isinf(ub[i]) 
		   && ub[i] - x[i] < step && ub[i] > x[i] + 1e-4)
		    step = ub[i] - x[i];
	       if (!nlopt_isinf(lb[i]) 
		   && x[i] - lb[i] < step && x[i] > lb[i] + 1e-4)
		    step = x[i] - lb[i];
	  }
     if (nlopt_isinf(step)) {
	  step = 0;
	  for (i = 0; i < n; ++i)
	       if (fabs(x[i]) * 0.25 > step)
		    step = fabs(x[i]) * 0.25;
	  if (step == 0)
	       step = 1;
     }
     return step;
}
	      
/*************************************************************************/

typedef struct {
     nlopt_func f;
     void *f_data;
     const double *lb, *ub;
} nlopt_data;

#include "subplex.h"
#include "praxis.h"

static double f_subplex(int n, const double *x, void *data_)
{
     int i;
     nlopt_data *data = (nlopt_data *) data_;
     double f;

     /* subplex does not support bound constraints, but it supports
	discontinuous objectives so we can just return Inf for invalid x */
     for (i = 0; i < n; ++i)
	  if (x[i] < data->lb[i] || x[i] > data->ub[i])
	       return MY_INF;

     f = data->f(n, x, NULL, data->f_data);
     return (isnan(f) || nlopt_isinf(f) ? MY_INF : f);
}

static double f_noderiv(int n, const double *x, void *data_)
{
     nlopt_data *data = (nlopt_data *) data_;
     return data->f(n, x, NULL, data->f_data);
}

#include "direct.h"

static double f_direct(int n, const double *x, int *undefined, void *data_)
{
     nlopt_data *data = (nlopt_data *) data_;
     double f;
     f = data->f(n, x, NULL, data->f_data);
     *undefined = isnan(f) || nlopt_isinf(f);
     return f;
}

#ifdef WITH_CXX
#  include "stogo.h"
#endif

#include "cdirect.h"

#ifdef WITH_NOCEDAL
#  include "l-bfgs-b.h"
#endif

#include "luksan.h"

#include "crs.h"

#include "mlsl.h"
#include "mma.h"
#include "cobyla.h"
#include "newuoa.h"
#include "neldermead.h"

/*************************************************************************/

/* for "hybrid" algorithms that combine global search with some
   local search algorithm, most of the time we anticipate that people
   will stick with the default local search algorithm, so we
   don't add this as a parameter to nlopt_minimize.  Instead, the user
   can call nlopt_{set/get}_hybrid_local_algorithm to get/set the defaults. */

/* default local-search algorithms */
static nlopt_algorithm local_search_alg_deriv = NLOPT_LD_MMA;
static nlopt_algorithm local_search_alg_nonderiv = NLOPT_LN_COBYLA;

static int local_search_maxeval = -1; /* no maximum by default */

void nlopt_get_local_search_algorithm(nlopt_algorithm *deriv,
				      nlopt_algorithm *nonderiv,
				      int *maxeval)
{
     *deriv = local_search_alg_deriv;
     *nonderiv = local_search_alg_nonderiv;
     *maxeval = local_search_maxeval;
}

void nlopt_set_local_search_algorithm(nlopt_algorithm deriv,
				      nlopt_algorithm nonderiv,
				      int maxeval)
{
     local_search_alg_deriv = deriv;
     local_search_alg_nonderiv = nonderiv;
     local_search_maxeval = maxeval;
}

/*************************************************************************/

/* same as nlopt_minimize, but xtol_abs is required to be non-NULL */
static nlopt_result nlopt_minimize_(
     nlopt_algorithm algorithm,
     int n, nlopt_func f, void *f_data,
     int m, nlopt_func fc, void *fc_data, ptrdiff_t fc_datum_size,
     const double *lb, const double *ub, /* bounds */
     double *x, /* in: initial guess, out: minimizer */
     double *minf, /* out: minimum */
     double minf_max, double ftol_rel, double ftol_abs,
     double xtol_rel, const double *xtol_abs,
     int maxeval, double maxtime)
{
     int i;
     nlopt_data d;
     nlopt_stopping stop;

     /* some basic argument checks */
     if (!minf || !f || n < 0 || m < 0
	  || (m > 0 && !fc)) return NLOPT_INVALID_ARGS;
     if (n == 0) { /* trivial case: no degrees of freedom */
	  *minf = f(n, x, NULL, f_data);
	  return NLOPT_SUCCESS;
     }
     else if (n < 0 || !lb || !ub || !x)
	  return NLOPT_INVALID_ARGS;

     /* nonlinear constraints are only supported with MMA or COBYLA */
     if (m != 0 && algorithm != NLOPT_LD_MMA && algorithm != NLOPT_LN_COBYLA) 
	  return NLOPT_INVALID_ARGS;

     d.f = f;
     d.f_data = f_data;
     d.lb = lb;
     d.ub = ub;

     /* make sure rand generator is inited */
     if (!nlopt_srand_called)
	  nlopt_srand_time(); /* default is non-deterministic */

     /* check bound constraints */
     for (i = 0; i < n; ++i)
	  if (lb[i] > ub[i] || x[i] < lb[i] || x[i] > ub[i])
	       return NLOPT_INVALID_ARGS;

     stop.n = n;
     stop.minf_max = (isnan(minf_max) 
		      || (nlopt_isinf(minf_max) && minf_max < 0))
	  ? -MY_INF : minf_max;
     stop.ftol_rel = ftol_rel;
     stop.ftol_abs = ftol_abs;
     stop.xtol_rel = xtol_rel;
     stop.xtol_abs = xtol_abs;
     stop.nevals = 0;
     stop.maxeval = maxeval;
     stop.maxtime = maxtime;
     stop.start = nlopt_seconds();

     switch (algorithm) {
	 case NLOPT_GN_DIRECT:
	 case NLOPT_GN_DIRECT_L: 
	 case NLOPT_GN_DIRECT_L_RAND: 
	      return cdirect(n, f, f_data, lb, ub, x, minf, &stop, 0.0, 
			     (algorithm != NLOPT_GN_DIRECT)
			     + 3 * (algorithm == NLOPT_GN_DIRECT_L_RAND ? 2 : (algorithm != NLOPT_GN_DIRECT))
			     + 9 * (algorithm == NLOPT_GN_DIRECT_L_RAND ? 1 : (algorithm != NLOPT_GN_DIRECT)));
	      
	 case NLOPT_GN_DIRECT_NOSCAL:
	 case NLOPT_GN_DIRECT_L_NOSCAL: 
	 case NLOPT_GN_DIRECT_L_RAND_NOSCAL: 
	      return cdirect_unscaled(n, f, f_data, lb, ub, x, minf, 
				      &stop, 0.0, 
				      (algorithm != NLOPT_GN_DIRECT)
				      + 3 * (algorithm == NLOPT_GN_DIRECT_L_RAND ? 2 : (algorithm != NLOPT_GN_DIRECT))
				      + 9 * (algorithm == NLOPT_GN_DIRECT_L_RAND ? 1 : (algorithm != NLOPT_GN_DIRECT)));
	      
	 case NLOPT_GN_ORIG_DIRECT:
	 case NLOPT_GN_ORIG_DIRECT_L: 
	      switch (direct_optimize(f_direct, &d, n, lb, ub, x, minf,
				      maxeval, -1, 0.0, 0.0,
				      pow(xtol_rel, (double) n), -1.0,
				      stop.minf_max, 0.0,
				      NULL, 
				      algorithm == NLOPT_GN_ORIG_DIRECT
				      ? DIRECT_ORIGINAL
				      : DIRECT_GABLONSKY)) {
		  case DIRECT_INVALID_BOUNDS:
		  case DIRECT_MAXFEVAL_TOOBIG:
		  case DIRECT_INVALID_ARGS:
		       return NLOPT_INVALID_ARGS;
		  case DIRECT_INIT_FAILED:
		  case DIRECT_SAMPLEPOINTS_FAILED:
		  case DIRECT_SAMPLE_FAILED:
		       return NLOPT_FAILURE;
		  case DIRECT_MAXFEVAL_EXCEEDED:
		  case DIRECT_MAXITER_EXCEEDED:
		       return NLOPT_MAXEVAL_REACHED;
		  case DIRECT_GLOBAL_FOUND:
		       return NLOPT_MINF_MAX_REACHED;
		  case DIRECT_VOLTOL:
		  case DIRECT_SIGMATOL:
		       return NLOPT_XTOL_REACHED;
		  case DIRECT_OUT_OF_MEMORY:
		       return NLOPT_OUT_OF_MEMORY;
	      break;
	 }

	 case NLOPT_GD_STOGO:
	 case NLOPT_GD_STOGO_RAND:
#ifdef WITH_CXX
	      if (!stogo_minimize(n, f, f_data, x, minf, lb, ub, &stop,
				  algorithm == NLOPT_GD_STOGO
				  ? 0 : 2*n))
		   return NLOPT_FAILURE;
	      break;
#else
	      return NLOPT_FAILURE;
#endif

	 case NLOPT_LN_SUBPLEX: {
	      int iret;
	      double *scale = (double *) malloc(sizeof(double) * n);
	      if (!scale) return NLOPT_OUT_OF_MEMORY;
	      for (i = 0; i < n; ++i)
		   scale[i] = initial_step(1, lb+i, ub+i, x+i);
	      iret = nlopt_subplex(f_subplex, minf, x, n, &d, &stop, scale);
	      free(scale);
	      switch (iret) {
		  case -2: return NLOPT_INVALID_ARGS;
		  case -10: return NLOPT_MAXTIME_REACHED;
		  case -1: return NLOPT_MAXEVAL_REACHED;
		  case 0: return NLOPT_XTOL_REACHED;
		  case 1: return NLOPT_SUCCESS;
		  case 2: return NLOPT_MINF_MAX_REACHED;
		  case 20: return NLOPT_FTOL_REACHED;
		  case -200: return NLOPT_OUT_OF_MEMORY;
		  default: return NLOPT_FAILURE; /* unknown return code */
	      }
	      break;
	 }

	 case NLOPT_LN_PRAXIS:
	      return praxis_(0.0, DBL_EPSILON, 
			     initial_step(n, lb, ub, x), n, x, f_subplex, &d,
			     &stop, minf);

#ifdef WITH_NOCEDAL
	 case NLOPT_LD_LBFGS_NOCEDAL: {
	      int iret, *nbd = (int *) malloc(sizeof(int) * n);
	      if (!nbd) return NLOPT_OUT_OF_MEMORY;
	      for (i = 0; i < n; ++i) {
		   int linf = nlopt_isinf(lb[i]) && lb[i] < 0;
		   int uinf = nlopt_isinf(ub[i]) && ub[i] > 0;
		   nbd[i] = linf && uinf ? 0 : (uinf ? 1 : (linf ? 3 : 2));
	      }
	      iret = lbfgsb_minimize(n, f, f_data, x, nbd, lb, ub,
				     MIN(n, 5), 0.0, ftol_rel, 
				     xtol_abs ? *xtol_abs : xtol_rel,
				     maxeval);
	      free(nbd);
	      if (iret <= 0) {
		   switch (iret) {
		       case -1: return NLOPT_INVALID_ARGS;
		       case -2: default: return NLOPT_FAILURE;
		   }
	      }
	      else {
		   *minf = f(n, x, NULL, f_data);
		   switch (iret) {
		       case 5: return NLOPT_MAXEVAL_REACHED;
		       case 2: return NLOPT_XTOL_REACHED;
		       case 1: return NLOPT_FTOL_REACHED;
		       default: return NLOPT_SUCCESS;
		   }
	      }
	      break;
	 }
#endif

	 case NLOPT_LD_LBFGS: 
	      return luksan_plis(n, f, f_data, lb, ub, x, minf, &stop);

	 case NLOPT_LD_VAR1: 
	 case NLOPT_LD_VAR2: 
	      return luksan_plip(n, f, f_data, lb, ub, x, minf, &stop,
		   algorithm == NLOPT_LD_VAR1 ? 1 : 2);

	 case NLOPT_LD_TNEWTON: 
	 case NLOPT_LD_TNEWTON_RESTART: 
	 case NLOPT_LD_TNEWTON_PRECOND: 
	 case NLOPT_LD_TNEWTON_PRECOND_RESTART: 
	      return luksan_pnet(n, f, f_data, lb, ub, x, minf, &stop,
				 1 + (algorithm - NLOPT_LD_TNEWTON) % 2,
				 1 + (algorithm - NLOPT_LD_TNEWTON) / 2);

	 case NLOPT_GN_CRS2_LM:
	      return crs_minimize(n, f, f_data, lb, ub, x, minf, &stop, 0);

	 case NLOPT_GN_MLSL:
	 case NLOPT_GD_MLSL:
	 case NLOPT_GN_MLSL_LDS:
	 case NLOPT_GD_MLSL_LDS:
	      return mlsl_minimize(n, f, f_data, lb, ub, x, minf, &stop,
				   (algorithm == NLOPT_GN_MLSL ||
				    algorithm == NLOPT_GN_MLSL_LDS)
				   ? local_search_alg_nonderiv
				   : local_search_alg_deriv,
				   local_search_maxeval,
				   algorithm >= NLOPT_GN_MLSL_LDS);

	 case NLOPT_LD_MMA:
	      return mma_minimize(n, f, f_data, m, fc, fc_data, fc_datum_size,
				  lb, ub, x, minf, &stop,
				  local_search_alg_deriv, 1e-12, 100000);

	 case NLOPT_LN_COBYLA:
	      return cobyla_minimize(n, f, f_data, 
				     m, fc, fc_data, fc_datum_size,
				     lb, ub, x, minf, &stop,
				     initial_step(n, lb, ub, x));
				     
	 case NLOPT_LN_NEWUOA:
	      return newuoa(n, 2*n+1, x, 0, 0, initial_step(n, lb, ub, x),
			    &stop, minf, f_noderiv, &d);
				     
	 case NLOPT_LN_NEWUOA_BOUND:
	      return newuoa(n, 2*n+1, x, lb, ub, initial_step(n, lb, ub, x),
			    &stop, minf, f_noderiv, &d);

	 case NLOPT_LN_NELDERMEAD: {
	      nlopt_result ret;
              double *scale = (double *) malloc(sizeof(double) * n);
              if (!scale) return NLOPT_OUT_OF_MEMORY;
              for (i = 0; i < n; ++i)
		   scale[i] = initial_step(1, lb+i, ub+i, x+i);
              ret = nldrmd_minimize(n, f,f_data, lb,ub,x, minf,scale,&stop);
	      free(scale);
	      return ret;
	 }

	 default:
	      return NLOPT_INVALID_ARGS;
     }

     return NLOPT_SUCCESS;
}

nlopt_result nlopt_minimize_constrained(
     nlopt_algorithm algorithm,
     int n, nlopt_func f, void *f_data,
     int m, nlopt_func fc, void *fc_data, ptrdiff_t fc_datum_size,
     const double *lb, const double *ub, /* bounds */
     double *x, /* in: initial guess, out: minimizer */
     double *minf, /* out: minimum */
     double minf_max, double ftol_rel, double ftol_abs,
     double xtol_rel, const double *xtol_abs,
     int maxeval, double maxtime)
{
     nlopt_result ret;
     if (xtol_abs)
	  ret = nlopt_minimize_(algorithm, n, f, f_data,
				m, fc, fc_data, fc_datum_size, lb, ub,
				x, minf, minf_max, ftol_rel, ftol_abs,
				xtol_rel, xtol_abs, maxeval, maxtime);
     else {
	  int i;
	  double *xtol = (double *) malloc(sizeof(double) * n);
	  if (!xtol) return NLOPT_OUT_OF_MEMORY;
	  for (i = 0; i < n; ++i) xtol[i] = -1;
	  ret = nlopt_minimize_(algorithm, n, f, f_data, 
				m, fc, fc_data, fc_datum_size, lb, ub,
				x, minf, minf_max, ftol_rel, ftol_abs,
				xtol_rel, xtol, maxeval, maxtime);
	  free(xtol);
     }
     return ret;
}

nlopt_result nlopt_minimize(
     nlopt_algorithm algorithm,
     int n, nlopt_func f, void *f_data,
     const double *lb, const double *ub, /* bounds */
     double *x, /* in: initial guess, out: minimizer */
     double *minf, /* out: minimum */
     double minf_max, double ftol_rel, double ftol_abs,
     double xtol_rel, const double *xtol_abs,
     int maxeval, double maxtime)
{
     return nlopt_minimize_constrained(
	  algorithm, n, f, f_data, 0, NULL, NULL, 0,
	  lb, ub, x, minf, minf_max, ftol_rel, ftol_abs,
	  xtol_rel, xtol_abs, maxeval, maxtime);
}
