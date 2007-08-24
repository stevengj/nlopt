#include <stdlib.h>
#include <math.h>

#include "nlopt.h"
#include "config.h"

static const char nlopt_algorithm_names[NLOPT_NUM_ALGORITHMS][128] = {
     "DIRECT (global)",
     "DIRECT-L (global)",
     "Subplex (local)",
     "StoGO (global)",
     "Low-storage BFGS (LBFGS) (local)"
};

const char *nlopt_algorithm_name(nlopt_algorithm a)
{
     if (a < 0 || a >= NLOPT_NUM_ALGORITHMS) return "UNKNOWN";
     return nlopt_algorithm_names[a];
}

static int my_isinf(double x) {
     return x == HUGE_VAL
#ifdef HAVE_ISINF
	  || isinf(x)
#endif
	  ;
}

#ifndef HAVE_ISNAN
static int my_isnan(double x) { return x != x; }
#  define isnan my_isnan
#endif

typedef struct {
     nlopt_func f;
     void *f_data;
} nlopt_data;

#include "subplex.h"

static double f_subplex(int n, const double *x, void *data_)
{
     nlopt_data *data = (nlopt_data *) data_;
     return data->f(n, x, NULL, data->f_data);
}

#include "direct.h"

static double f_direct(int n, const double *x, int *undefined, void *data_)
{
     nlopt_data *data = (nlopt_data *) data_;
     double f = data->f(n, x, NULL, data->f_data);
     *undefined = isnan(f) || my_isinf(f);
     return f;
}

#include "stogo.h"
#include "l-bfgs-b.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))

nlopt_result nlopt_minimize(
     nlopt_algorithm algorithm,
     int n, nlopt_func f, void *f_data,
     const double *lb, const double *ub, /* bounds */
     double *x, /* in: initial guess, out: minimizer */
     double *fmin, /* out: minimum */
     double fmin_max, double ftol_rel, double ftol_abs,
     double xtol_rel, const double *xtol_abs,
     int maxeval, double maxtime)
{
     nlopt_data d;
     d.f = f;
     d.f_data = f_data;

     switch (algorithm) {
	 case NLOPT_GLOBAL_DIRECT:
	 case NLOPT_GLOBAL_DIRECT_L:
	      switch (direct_optimize(f_direct, &d, n, lb, ub, x, fmin,
				      maxeval, 500, ftol_rel, ftol_abs,
				      xtol_rel, xtol_rel,
				      DIRECT_UNKNOWN_FGLOBAL, -1.0,
				      NULL, 
				      algorithm == NLOPT_GLOBAL_DIRECT
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
		       return NLOPT_SUCCESS;
		  case DIRECT_VOLTOL:
		  case DIRECT_SIGMATOL:
		       return NLOPT_XTOL_REACHED;
		  case DIRECT_OUT_OF_MEMORY:
		       return NLOPT_OUT_OF_MEMORY;
	      }
	      break;

	 case NLOPT_GLOBAL_STOGO:
	      if (!stogo_minimize(n, f, f_data, x, fmin, lb, ub,
				  maxeval, maxtime))
		   return NLOPT_FAILURE;
	      break;

	 case NLOPT_LOCAL_SUBPLEX: {
	      int iret, i;
	      double *scale = (double *) malloc(sizeof(double) * n);
	      if (!scale) return NLOPT_OUT_OF_MEMORY;
	      for (i = 0; i < n; ++i)
		   scale[i] = fabs(ub[i] - lb[i]);
	      iret = subplex(f_subplex, fmin, x, n, &d, xtol_rel, maxeval,
			     fmin_max, !my_isinf(fmin_max), scale);
	      free(scale);
	      switch (iret) {
		  case -2: return NLOPT_INVALID_ARGS;
		  case -1: return NLOPT_MAXEVAL_REACHED;
		  case 0: return NLOPT_XTOL_REACHED;
		  case 1: return NLOPT_SUCCESS;
		  case 2: return NLOPT_FMIN_MAX_REACHED;
	      }
	      break;
	 }

	 case NLOPT_LOCAL_LBFGS: {
	      int iret, i, *nbd = (int *) malloc(sizeof(int) * n);
	      if (!nbd) return NLOPT_OUT_OF_MEMORY;
	      for (i = 0; i < n; ++i) {
		   int linf = my_isinf(lb[i]) && lb[i] < 0;
		   int uinf = my_isinf(ub[i]) && ub[i] > 0;
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
		   *fmin = f(n, x, NULL, f_data);
		   switch (iret) {
		       case 5: return NLOPT_MAXEVAL_REACHED;
		       case 2: return NLOPT_XTOL_REACHED;
		       case 1: return NLOPT_FTOL_REACHED;
		       default: return NLOPT_SUCCESS;
		   }
	      }
	      break;
	 }

	 default:
	      return NLOPT_INVALID_ARGS;
     }

     return NLOPT_SUCCESS;
}
