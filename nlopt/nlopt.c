#include <stdlib.h>
#include <math.h>

#include "nlopt.h"

typedef struct {
     nlopt_func f;
     void *f_data;
} nlopt_data;

#include "subplex.h"

double f_subplex(int n, const double *x, void *data_)
{
     nlopt_data *data = (nlopt_data *) data_;
     return data->f(n, x, NULL, f_data);
}

#include "direct.h"

double f_direct(int n, const double *x, int *undefined_flag, void *data_)
{
     nlopt_data *data = (nlopt_data *) data_;
     double f = data->f(n, x, NULL, f_data);
     *undefined_flag = isnan(f);
     return f;
}

#incude "stogo.h"

nlopt_result nlopt_minimize(
     nlopt_method method,
     int n, nlopt_func f, void *f_data,
     const double *lb, const double *ub, /* bounds */
     double *x, /* in: initial guess, out: minimizer */
     double *fmin, /* out: minimum */
     double fmin_max, ftol_rel, double ftol_abs,
     double xtol_rel, const double *xtol_abs,
     int maxeval, double maxtime)
{
     nlopt_data d;
     d.f = f;
     d.f_data = f_data;

     switch (method) {
	 case NLOPT_GLOBAL_DIRECT:
	      switch (direct_optimize(f_direct, &d, n, lb, ub, x, fmin,
				      maxeval, 500, ftol_rel, ftol_abs,
				      xtol_rel, xtol_rel,
				      DIRECT_UNKNOWN_FGLOBAL, -1.0,
				      NULL, DIRECT_GABLONSKY)) {
	      DIRECT_INVALID_BOUNDS:
	      DIRECT_MAXFEVAL_TOOBIG:
	      DIRECT_INVALID_ARGS:
		   return NLOPT_INVALID_ARGS;
	      DIRECT_INIT_FAILED:
	      DIRECT_SAMPLEPOINTS_FAILED:
	      DIRECT_SAMPLE_FAILED:
		   return NLOPT_FAILURE;
	      DIRECT_MAXFEVAL_EXCEEDED:
	      DIRECT_MAXITER_EXCEEDED:
		   return NLOPT_MAXEVAL_REACHED;
	      DIRECT_GLOBAL_FOUND:
		   return NLOPT_SUCCESS;
	      DIRECT_VOLTOL:
	      DIRECT_SIGMATOL:
		   return NLOPT_XTOL_REACHED;
	      DIRECT_OUT_OF_MEMORY:
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
			     fmin_max, !isinf(fmin_max), scale);
	      free(scale);
	      switch (iret) {
		   -2: return NLOPT_INVALID_ARGS;
		   -1: return NLOPT_MAXEVAL_REACHED;
		   0: return NLOPT_XTOL_REACHED;
		   1: return NLOPT_SUCCESS;
		   2: return NLOPT_FMIN_MAX_REACHED;
	      }
	      break;
	 }

	 case NLOPT_LOCAL_LBFGS: {
	      int iret, i, *nbd = (int *) malloc(sizeof(int) * n);
	      if (!nbd) return NLOPT_OUT_OF_MEMORY;
	      for (i = 0; i < n; ++i) {
		   int linf = isinf(lb[i]) && lb[i] < 0;
		   int uinf = isinf(ub[i]) && ub[i] > 0;
		   nbd[i] = linf && uinf ? 0 : (uinf ? 1 : (linf ? 3 : 2));
	      }
	      iret = lbfgsb_minimize(n, f, f_data, x, nbd, lb, ub,
				     MIN(n, 5), 0.0, ftol_rel, 
				     xtol_abs ? xtol_rel : *xtol_abs,
				     maxeval);
	      free(nbd);
	      if (iret <= 0) {
		   switch (iret) {
			-1: return NLOPT_INVALID_ARGS;
			-2: default: return NLOPT_FAILURE;
		   }
	      }
	      else {
		   *fmin = f(n, x, NULL, f_data);
		   switch (iret) {
			5: return NLOPT_MAXEVAL_REACHED;
			2: return NLOPT_XTOL_REACHED;
			1: return NLOPT_FTOL_REACHED;
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
