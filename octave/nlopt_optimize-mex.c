/* Copyright (c) 2007-2014 Massachusetts Institute of Technology
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

/* MATLAB MEX interface to NLopt, and in particular to nlopt_optimize */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mex.h"
#include "nlopt.h"

#define ERRID "nlopt:error"
#define CHECK1(cond, msg) if (!(cond)) { nlopt_destroy(opt); mexWarnMsgTxt(msg); return NULL; };
#define CHECK(cond, msg) if (!(cond)) { mxFree(dh); mxFree(dfc); nlopt_destroy(opt); mexErrMsgIdAndTxt(ERRID, msg); }

#define NAMELENGTHMAX 64 /* TMW_NAME_LENGTH_MAX: max length of varname */
typedef struct user_function_data_s {
    char f[NAMELENGTHMAX];
    mxArray *plhs[2];
    mxArray *prhs[3];
    int xrhs, nrhs;
    int verbose, neval;
    struct user_function_data_s *dpre;
    nlopt_opt opt;
} user_function_data;

static const char *output_fields[] = {
    "algorithm", "funcCount", "iterations", "message"};

static bool mx_isscalar(const mxArray *arr)
{
    return (mxIsDouble(arr) && !mxIsComplex(arr) && !mxIsSparse(arr) &&
        mxIsScalar(arr));
}

static bool mx_isvector(const mxArray *arr)
{
    return (mxIsDouble(arr) && !mxIsComplex(arr) && !mxIsSparse(arr) &&
        mxGetNumberOfDimensions(arr) == 2 &&
        (mxGetM(arr) == 1 || mxGetN(arr) == 1));
}

static bool mx_isvector_len(const mxArray *arr, unsigned n)
{
    return (mx_isvector(arr) && mxGetNumberOfElements(arr) == n);
}

static bool mx_isfunction(const mxArray *arr)
{
    return (mxIsFunctionHandle(arr) || (mxIsChar(arr) &&
        mxGetNumberOfDimensions(arr) == 2 && mxGetM(arr) == 1));
}

static double struct_val(const mxArray *s, const char *name, double dflt)
{
     mxArray *val = mxGetField(s, 0, name);
     if (val) {
	  if (!mx_isscalar(val))
		mexErrMsgIdAndTxt(ERRID, "opt.%s must be a real scalar", name);
	  return mxGetScalar(val);
     }
     return dflt;
}

static double *struct_arrval(const mxArray *s, const char *name, unsigned n,
			     double *dflt)
{
     mxArray *val = mxGetField(s, 0, name);
     if (val) {
	  if (!mx_isvector_len(val, n))
		mexErrMsgIdAndTxt(ERRID,
			"opt.%s must be a real vector of length %u", name, n);
	  return mxGetPr(val);
     }
     return dflt;
}

static mxArray *struct_funcval(const mxArray *s, const char *name)
{
     mxArray *val = mxGetField(s, 0, name);
     if (val) {
	  if (!mx_isfunction(val))
		 mexErrMsgIdAndTxt(ERRID,
		 	"opt.%s must be a function handle or name", name);
	  return val;
     }
     return NULL;
}

static mxArray *cell_funcval(const mxArray *c, unsigned i)
{
     mxArray *val = mxGetCell(c, i);
     if (val) {
	  if (!mx_isfunction(val))
		 mexErrMsgIdAndTxt(ERRID,
		 	"opt constraint {%u} must be a function handle or name", i);
	  return val;
     }
     return NULL;
}

static double user_function(unsigned n, const double *x,
			    double *grad, void *data)
{
  user_function_data *d = (user_function_data *) data;
  double f;

  /* x */
  memcpy(mxGetPr(d->prhs[d->xrhs]), x, n * sizeof(double));

  /* [f, g] = feval(objFunc, x) */
  d->plhs[0] = d->plhs[1] = NULL;
  if (mexCallMATLAB(grad ? 2 : 1, d->plhs, d->nrhs, d->prhs, d->f) != 0)
	mexErrMsgIdAndTxt(ERRID, "error calling objective function");

  /* f */
  if (!mx_isscalar(d->plhs[0]))
	mexErrMsgIdAndTxt(ERRID, "objective function must return a real scalar");
  f = mxGetScalar(d->plhs[0]);
  mxDestroyArray(d->plhs[0]);

  /* g */
  if (grad) {
     if (!mx_isvector_len(d->plhs[1], n))
	   mexErrMsgIdAndTxt(ERRID,
	     "objective function must return a gradient vector of length %u", n);
     memcpy(grad, mxGetPr(d->plhs[1]), n * sizeof(double));
     mxDestroyArray(d->plhs[1]);
  }

  d->neval++;
  if (d->verbose) mexPrintf("nlopt_optimize eval #%d: %g\n", d->neval, f);
  if (mxIsNaN(f)) nlopt_force_stop(d->opt);
  return f;
}

static void user_pre(unsigned n, const double *x, const double *v,
		       double *vpre, void *data)
{
  user_function_data *d = ((user_function_data *) data)->dpre;

  /* x and v */
  memcpy(mxGetPr(d->prhs[d->xrhs]), x, n * sizeof(double));
  memcpy(mxGetPr(d->prhs[d->xrhs + 1]), v, n * sizeof(double));

  /* vpre = feval(precondFunc, x, v) */
  d->plhs[0] = NULL;
  if (mexCallMATLAB(1, d->plhs, d->nrhs, d->prhs, d->f) != 0)
	 mexErrMsgIdAndTxt(ERRID, "error calling preconditioned function");

  /* vpre */
  if (!mx_isvector_len(d->plhs[0], n))
	 mexErrMsgIdAndTxt(ERRID,
	 	"preconditioned function must return a vpre vector of length %u", n);
  memcpy(vpre, mxGetPr(d->plhs[0]), n * sizeof(double));
  mxDestroyArray(d->plhs[0]);

  d->neval++;
  if (d->verbose) mexPrintf("nlopt_optimize precond eval #%d\n", d->neval);
}

static nlopt_opt make_opt(const mxArray *s, unsigned n)
{
     nlopt_opt opt = NULL;
     nlopt_algorithm alg;
     double *tmp;
     const mxArray *s_local;

     /* algorithm */
     alg = (nlopt_algorithm) struct_val(s, "algorithm", NLOPT_NUM_ALGORITHMS);
     if (alg < 0 || alg >= NLOPT_NUM_ALGORITHMS)
	    mexErrMsgIdAndTxt(ERRID, "opt.algorithm is invalid: %d", alg);

     /* create object */
     opt = nlopt_create(alg, n);
     if (!opt)
     	mexErrMsgIdAndTxt(ERRID, "error creating nlopt object");

     /* bound constraints */
     tmp = struct_arrval(s, "lower_bounds", n, NULL);
     if (tmp)
         nlopt_set_lower_bounds(opt, tmp);
     else
         nlopt_set_lower_bounds1(opt, -HUGE_VAL);
     tmp = struct_arrval(s, "upper_bounds", n, NULL);
     if (tmp)
         nlopt_set_upper_bounds(opt, tmp);
     else
         nlopt_set_upper_bounds1(opt, +HUGE_VAL);

     /* stopping criteria */
     nlopt_set_stopval(opt, struct_val(s, "stopval", -HUGE_VAL));
     nlopt_set_ftol_rel(opt, struct_val(s, "ftol_rel", 0.0));
     nlopt_set_ftol_abs(opt, struct_val(s, "ftol_abs", 0.0));
     nlopt_set_xtol_rel(opt, struct_val(s, "xtol_rel", 0.0));
     tmp = struct_arrval(s, "xtol_abs", n, NULL);
     if (tmp)
         nlopt_set_xtol_abs(opt, tmp);
     else
         nlopt_set_xtol_abs1(opt, 0.0);
     nlopt_set_maxeval(opt, (int) struct_val(s, "maxeval", 0));
     nlopt_set_maxtime(opt, struct_val(s, "maxtime", 0.0));

     /* stochastic population */
     nlopt_set_population(opt, (unsigned) struct_val(s, "population", 0));

     /* vector storage */
     nlopt_set_vector_storage(opt,
     	(unsigned) struct_val(s, "vector_storage", 0));

     /* initial step size */
	 nlopt_set_initial_step(opt, struct_arrval(s, "initial_step", n, NULL));

     /* local optimization algorithm */
	 s_local = mxGetField(s, 0, "local_optimizer");
     if (s_local) {
	  nlopt_opt opt_local = NULL;
	  CHECK1(mxIsStruct(s_local) && mxIsScalar(s_local),
		 "opt.local_optimizer must be a scalar structure");
	  opt_local = make_opt(s_local, n);
	  CHECK1(opt_local, "error initializing local optimizer options");
	  nlopt_set_local_optimizer(opt, opt_local);
	  nlopt_destroy(opt_local);
     }

     return opt;
}

static void debug_opt(const nlopt_opt opt, const double *x)
{
    unsigned i, n;
    int version[3];
    double *tmp;

    n = nlopt_get_dimension(opt);
    tmp = (double *) mxCalloc(n, sizeof(double));

    nlopt_version(version+0, version+1, version+2);
    mexPrintf("NLOpt v%d.%d.%d\n", version[0], version[1], version[2]);

    mexPrintf(" dimension = %u\n", n);

    mexPrintf(" algorithm = %s\n",
        nlopt_algorithm_name(nlopt_get_algorithm(opt)));

    nlopt_get_lower_bounds(opt, tmp);
    mexPrintf(" lower_bounds = [");
    for (i = 0; i < n; ++i) mexPrintf("%g, ", tmp[i]);
    mexPrintf("]\n");

    nlopt_get_upper_bounds(opt, tmp);
    mexPrintf(" upper_bounds = [");
    for (i = 0; i < n; ++i) mexPrintf("%g, ", tmp[i]);
    mexPrintf("]\n");

    mexPrintf(" stopval = %g\n", nlopt_get_stopval(opt));
    mexPrintf(" ftol_rel = %g\n", nlopt_get_ftol_rel(opt));
    mexPrintf(" ftol_abs = %g\n", nlopt_get_ftol_abs(opt));
    mexPrintf(" xtol_rel = %g\n", nlopt_get_xtol_rel(opt));

    nlopt_get_xtol_abs(opt, tmp);
    mexPrintf(" xtol_abs = [");
    for (i = 0; i < n; ++i) mexPrintf("%g, ", tmp[i]);
    mexPrintf("]\n");

    mexPrintf(" maxeval = %d\n", nlopt_get_maxeval(opt));
    mexPrintf(" maxtime = %g\n", nlopt_get_maxtime(opt));

    nlopt_get_initial_step(opt, x, tmp);
    mexPrintf(" initial_step = [");
    for (i = 0; i < n; ++i) mexPrintf("%g, ", tmp[i]);
    mexPrintf("]\n");

    mexPrintf(" population = %u\n", nlopt_get_population(opt));
    mexPrintf(" vector_storage = %u\n", nlopt_get_vector_storage(opt));

    mxFree(tmp);
}

static const char *translate_result(nlopt_result ret)
{
    switch (ret) {
        case NLOPT_FAILURE: return "NLOPT_FAILURE";
        case NLOPT_INVALID_ARGS: return "NLOPT_INVALID_ARGS";
        case NLOPT_OUT_OF_MEMORY: return "NLOPT_OUT_OF_MEMORY";
        case NLOPT_ROUNDOFF_LIMITED: return "NLOPT_ROUNDOFF_LIMITED";
        case NLOPT_FORCED_STOP: return "NLOPT_FORCED_STOP";
        case NLOPT_SUCCESS: return "NLOPT_SUCCESS";
        case NLOPT_STOPVAL_REACHED: return "NLOPT_STOPVAL_REACHED";
        case NLOPT_FTOL_REACHED: return "NLOPT_FTOL_REACHED";
        case NLOPT_XTOL_REACHED: return "NLOPT_XTOL_REACHED";
        case NLOPT_MAXEVAL_REACHED: return "NLOPT_MAXEVAL_REACHED";
        case NLOPT_MAXTIME_REACHED: return "NLOPT_MAXTIME_REACHED";
        default: return "";
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
     unsigned n, j, m;
     double *x, opt_f, *tol;
     nlopt_result ret;
     mxArray *mx;
     user_function_data d, dpre, *dfc = NULL, *dh = NULL;
     nlopt_opt opt = NULL;

     /* [x, fval, exitflag, output] = nlopt_optimize(opt, x0) */
     if (nrhs == 0 && nlhs == 0) {
     	for (j = 0; j < NLOPT_NUM_ALGORITHMS; ++j)
     		mexPrintf("%2d = %s\n", j,
     			nlopt_algorithm_name((nlopt_algorithm) j));
     	return;
     }
     else if (nrhs != 2 || nlhs > 4)
 		mexErrMsgIdAndTxt(ERRID, "wrong number of arguments");

     /* x0 */
     if (!mx_isvector(prhs[1]))
     	mexErrMsgIdAndTxt(ERRID, "x0 must be a real row or column vector");
     n = mxGetNumberOfElements(prhs[1]);

     /* opt */
     if (!mxIsStruct(prhs[0]) || !mxIsScalar(prhs[0]))
     	mexErrMsgIdAndTxt(ERRID, "opt must be a scalar struct");
     opt = make_opt(prhs[0], n);
     if (!opt)
     	mexErrMsgIdAndTxt(ERRID, "error initializing nlopt options");

     /* random seed */
     if (struct_val(prhs[0], "seed", -1) >= 0)
	   nlopt_srand((unsigned long) struct_val(prhs[0], "seed", 0));

     /* objective function */
     mx = struct_funcval(prhs[0], "min_objective");
     if (!mx) mx = struct_funcval(prhs[0], "max_objective");
     CHECK(mx, "opt.(min|max)_objective must be set");
     if (mxIsChar(mx)) {
	  CHECK(mxGetString(mx, d.f, NAMELENGTHMAX) == 0,
		"error reading function name string from opt.(min|max)_objective");
	  d.nrhs = 1;
	  d.xrhs = 0;
     }
     else {
	  d.prhs[0] = mx;
	  strcpy(d.f, "feval");
	  d.nrhs = 2;
	  d.xrhs = 1;
     }
     d.verbose = (int) struct_val(prhs[0], "verbose", 0);
     d.neval = 0;
     d.opt = opt;
     d.prhs[d.xrhs] = mxCreateDoubleMatrix(1, n, mxREAL);

     /* preconditioned objective function */
     mx = struct_funcval(prhs[0], "pre");
     if (!mx) {
	  dpre.nrhs = 0;  /* indicate no memory allocation */
	  if (struct_funcval(prhs[0], "min_objective"))
	       nlopt_set_min_objective(opt, user_function, &d);
	  else
	       nlopt_set_max_objective(opt, user_function, &d);
     }
     else {
	  if (mxIsChar(mx)) {
	       CHECK(mxGetString(mx, dpre.f, NAMELENGTHMAX) == 0,
                     "error reading function name string from opt.pre");
	       dpre.nrhs = 2;
	       dpre.xrhs = 0;
	  }
	  else {
	       dpre.prhs[0] = mx;
	       strcpy(dpre.f, "feval");
	       dpre.nrhs = 3;
	       dpre.xrhs = 1;
	  }
	  dpre.verbose = d.verbose > 2;
	  dpre.neval = 0;
	  dpre.opt = opt;
	  dpre.prhs[dpre.xrhs] = d.prhs[d.xrhs];
	  dpre.prhs[dpre.xrhs+1] = mxCreateDoubleMatrix(1, n, mxREAL);
	  d.dpre = &dpre;
	  if (struct_funcval(prhs[0], "min_objective"))
	       nlopt_set_precond_min_objective(opt, user_function, user_pre, &d);
	  else
	       nlopt_set_precond_max_objective(opt, user_function, user_pre, &d);
     }

     /* nonlinear inequality constraints */
     mx = mxGetField(prhs[0], 0, "fc");
     if (mx) {
	  CHECK(mxIsCell(mx), "opt.fc must be a cell array");
	  m = mxGetNumberOfElements(mx);
	  dfc = (user_function_data *) mxCalloc(m, sizeof(user_function_data));
	  tol = struct_arrval(prhs[0], "fc_tol", m, NULL);
	  for (j = 0; j < m; ++j) {
	       mxArray *fc = cell_funcval(mx, j);
	       if (mxIsChar(fc)) {
		    CHECK(mxGetString(fc, dfc[j].f, NAMELENGTHMAX) == 0,
		     "error reading function name string from opt.fc");
		    dfc[j].nrhs = 1;
		    dfc[j].xrhs = 0;
	       }
	       else {
		    dfc[j].prhs[0] = fc;
		    strcpy(dfc[j].f, "feval");
		    dfc[j].nrhs = 2;
		    dfc[j].xrhs = 1;
	       }
	       dfc[j].verbose = d.verbose > 1;
	       dfc[j].neval = 0;
	       dfc[j].opt = opt;
	       dfc[j].prhs[dfc[j].xrhs] = d.prhs[d.xrhs];
	       ret = nlopt_add_inequality_constraint(opt, user_function, dfc + j,
		     tol ? tol[j] : 0.0);
	       CHECK(ret == NLOPT_SUCCESS,
		     "error adding inequality constraint opt.fc");
	  }
     }

     /* nonlinear equality constraints */
     mx = mxGetField(prhs[0], 0, "h");
     if (mx) {
	  CHECK(mxIsCell(mx), "opt.h must be a cell array");
	  m = mxGetNumberOfElements(mx);
	  dh = (user_function_data *) mxCalloc(m, sizeof(user_function_data));
	  tol = struct_arrval(prhs[0], "h_tol", m, NULL);
	  for (j = 0; j < m; ++j) {
	       mxArray *h = cell_funcval(mx, j);
	       if (mxIsChar(h)) {
		    CHECK(mxGetString(h, dh[j].f, NAMELENGTHMAX) == 0,
		     "error reading function name string from opt.h");
		    dh[j].nrhs = 1;
		    dh[j].xrhs = 0;
	       }
	       else {
		    dh[j].prhs[0] = h;
		    strcpy(dh[j].f, "feval");
		    dh[j].nrhs = 2;
		    dh[j].xrhs = 1;
	       }
	       dh[j].verbose = d.verbose > 1;
	       dh[j].neval = 0;
	       dh[j].opt = opt;
	       dh[j].prhs[dh[j].xrhs] = d.prhs[d.xrhs];
	       ret = nlopt_add_equality_constraint(opt, user_function, dh + j,
		     tol ? tol[j] : 0.0);
	       CHECK(ret == NLOPT_SUCCESS,
		     "error adding equality constraint opt.h");
	  }
     }

     /* x = x0 */
     plhs[0] = mxDuplicateArray(prhs[1]);
     x = mxGetPr(plhs[0]);

     /* optimize */
     if (d.verbose > 3) debug_opt(opt, x);
     ret = nlopt_optimize(opt, x, &opt_f);

     /* assign ouput arguments */
     if (nlhs > 1) plhs[1] = mxCreateDoubleScalar(opt_f);
     if (nlhs > 2) plhs[2] = mxCreateDoubleScalar((int) ret);
     if (nlhs > 3) {
         plhs[3] = mxCreateStructMatrix(1, 1, 4, output_fields);
         mxSetField(plhs[3], 0, "algorithm",
             mxCreateString(nlopt_algorithm_name(nlopt_get_algorithm(opt))));
         /* TODO: also separately return neval from: pre, fc, h, c, ceq */
         mxSetField(plhs[3], 0, "funcCount", mxCreateDoubleScalar(d.neval));
         /* TODO: is there a way to get number of iterations? */
         /*mxSetField(plhs[3], 0, "iterations", mxCreateDoubleScalar(0));*/
         mxSetField(plhs[3], 0, "message",
             mxCreateString(translate_result(ret)));
     }

     /* cleanup */
     mxFree(dh);
     mxFree(dfc);
     mxDestroyArray(d.prhs[d.xrhs]);
     if (dpre.nrhs > 0) mxDestroyArray(dpre.prhs[dpre.xrhs+1]);
     nlopt_destroy(opt);
}
