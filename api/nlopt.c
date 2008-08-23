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
     "Method of Moving Asymptotes (MMA) (local, derivative)"
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
     return (isnan(f) ? MY_INF : f);
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

/*************************************************************************/

/* for "hybrid" algorithms that combine global search with some
   local search algorithm, most of the time we anticipate that people
   will stick with the default local search algorithm, so we
   don't add this as a parameter to nlopt_minimize.  Instead, the user
   can call nlopt_{set/get}_hybrid_local_algorithm to get/set the defaults. */

/* default local-search algorithms */
static nlopt_algorithm local_search_alg_deriv = NLOPT_LD_LBFGS;
static nlopt_algorithm local_search_alg_nonderiv = NLOPT_LN_SUBPLEX;

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

     /* nonlinear constraints are only supported with MMA */
     if (m != 0 && algorithm != NLOPT_LD_MMA) return NLOPT_INVALID_ARGS;

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
	      for (i = 0; i < n; ++i) {
		   if (!nlopt_isinf(ub[i]) && !nlopt_isinf(lb[i]))
			scale[i] = (ub[i] - lb[i]) * 0.01;
		   else if (!nlopt_isinf(lb[i]) && x[i] > lb[i])
			scale[i] = (x[i] - lb[i]) * 0.01;
		   else if (!nlopt_isinf(ub[i]) && x[i] < ub[i])
			scale[i] = (ub[i] - x[i]) * 0.01;
		   else
			scale[i] = 0.01 * x[i] + 0.0001;
	      }
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

	 case NLOPT_LN_PRAXIS: {
	      double h0 = HUGE_VAL;
	      for (i = 0; i < n; ++i) {
		   if (!nlopt_isinf(ub[i]) && !nlopt_isinf(lb[i]))
			h0 = MIN(h0, (ub[i] - lb[i]) * 0.01);
		   else if (!nlopt_isinf(lb[i]) && x[i] > lb[i])
			h0 = MIN(h0, (x[i] - lb[i]) * 0.01);
		   else if (!nlopt_isinf(ub[i]) && x[i] < ub[i])
			h0 = MIN(h0, (ub[i] - x[i]) * 0.01);
		   else
			h0 = MIN(h0, 0.01 * x[i] + 0.0001);
	      }
	      return praxis_(0.0, DBL_EPSILON, h0, n, x, f_subplex, &d,
			     &stop, minf);
	 }

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
