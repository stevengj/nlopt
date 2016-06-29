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

#define CHECK0(cond, msg) if (!(cond)) mexErrMsgTxt(msg);
#define CHECK1(cond, msg) if (!(cond)) { nlopt_destroy(opt); mexWarnMsgTxt(msg); return NULL; };
#define CHECK(cond, msg) if (!(cond)) { mxFree(dh); mxFree(dfc); nlopt_destroy(opt); mexErrMsgTxt(msg); }

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

static double struct_val(const mxArray *s, const char *name, double dflt)
{
     mxArray *val = mxGetField(s, 0, name);
     if (val) {
	  CHECK0(mxIsDouble(val) && !mxIsComplex(val) && mxIsScalar(val),
		"opt fields, other than xtol_abs, must be real scalars");
	  return mxGetScalar(val);
     }
     return dflt;
}

static double *struct_arrval(const mxArray *s, const char *name, unsigned n,
			     double *dflt)
{
     mxArray *val = mxGetField(s, 0, name);
     if (val) {
	  CHECK0(mxIsDouble(val) && !mxIsComplex(val)
		&& mxGetNumberOfElements(val) == n,
		"opt vector field is not of length n");
	  return mxGetPr(val);
     }
     return dflt;
}

static mxArray *struct_funcval(const mxArray *s, const char *name)
{
     mxArray *val = mxGetField(s, 0, name);
     if (val) {
	  CHECK0(mxIsChar(val) || mxIsFunctionHandle(val),
		 "opt function field is not a function handle/name");
	  return val;
     }
     return NULL;
}

static mxArray *cell_funcval(const mxArray *c, unsigned i)
{
     mxArray *val = mxGetCell(c, i);
     if (val) {
	  CHECK0(mxIsChar(val) || mxIsFunctionHandle(val),
		 "opt constraint cell is not a function handle/name");
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
  CHECK0(0 == mexCallMATLAB(grad ? 2 : 1, d->plhs, d->nrhs, d->prhs, d->f),
	"error calling user function");

  /* f */
  CHECK0(mxIsDouble(d->plhs[0]) && !mxIsComplex(d->plhs[0])
	&& mxIsScalar(d->plhs[0]),
	"user function must return real scalar");
  f = mxGetScalar(d->plhs[0]);
  mxDestroyArray(d->plhs[0]);

  /* g */
  if (grad) {
     CHECK0(mxIsDouble(d->plhs[1]) && !mxIsComplex(d->plhs[1])
	   && (mxGetM(d->plhs[1]) == 1 || mxGetN(d->plhs[1]) == 1)
	   && mxGetNumberOfElements(d->plhs[1]) == n,
	   "gradient vector from user function is the wrong size");
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
  CHECK0(0 == mexCallMATLAB(1, d->plhs, d->nrhs, d->prhs, d->f),
	 "error calling user function");

  /* vpre */
  CHECK0(mxIsDouble(d->plhs[0]) && !mxIsComplex(d->plhs[0])
	 && (mxGetM(d->plhs[0]) == 1 || mxGetN(d->plhs[0]) == 1)
	 && mxGetNumberOfElements(d->plhs[0]) == n,
	 "vpre vector from user function is the wrong size");
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
     CHECK1(alg >= 0 && alg < NLOPT_NUM_ALGORITHMS, "invalid opt.algorithm");

     /* create object */
     opt = nlopt_create(alg, n);
     CHECK1(opt, "nlopt: out of memory");

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
		 "opt.local_optimizer must be a structure");
	  opt_local = make_opt(s_local, n);
	  CHECK1(opt_local, "error initializing local optimizer");
	  nlopt_set_local_optimizer(opt, opt_local);
	  nlopt_destroy(opt_local);
     }

     return opt;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
     unsigned n, j, m;
     double *x, opt_f, *tol;
     nlopt_result ret;
     mxArray *mx;
     user_function_data d, dpre, *dfc = NULL, *dh = NULL;
     nlopt_opt opt = NULL;

     /* [x, fval, exitflag] = nlopt_optimize(opt, x0) */
     CHECK(nrhs == 2 && nlhs <= 3, "wrong number of arguments");

     /* x0 */
     CHECK(mxIsDouble(prhs[1]) && !mxIsComplex(prhs[1])
	   && (mxGetNumberOfDimensions(prhs[1]) <= 2)
	   && (mxGetM(prhs[1]) == 1 || mxGetN(prhs[1]) == 1),
	   "x must be real row or column vector");
     n = mxGetNumberOfElements(prhs[1]);

     /* opt */
     CHECK(mxIsStruct(prhs[0]) && mxIsScalar(prhs[0]), "opt must be a struct");
     opt = make_opt(prhs[0], n);
     CHECK(opt, "error initializing nlopt options");

     /* objective function */
     mx = struct_funcval(prhs[0], "min_objective");
     if (!mx) mx = struct_funcval(prhs[0], "max_objective");
     CHECK(mx, "either opt.min_objective or opt.max_objective must exist");
     if (mxIsChar(mx)) {
	  CHECK(mxGetString(mx, d.f, NAMELENGTHMAX) == 0,
		"error reading function name string (too long?)");
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
                     "error reading function name string (too long?)");
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
	  CHECK(mxIsCell(mx), "fc must be a Cell array");
	  m = mxGetNumberOfElements(mx);
	  dfc = (user_function_data *) mxCalloc(m, sizeof(user_function_data));
	  tol = struct_arrval(prhs[0], "fc_tol", m, NULL);
	  for (j = 0; j < m; ++j) {
	       mxArray *fc = cell_funcval(mx, j);
	       if (mxIsChar(fc)) {
		    CHECK(mxGetString(fc, dfc[j].f, NAMELENGTHMAX) == 0,
		     "error reading function name string (too long?)");
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
	       CHECK(nlopt_add_inequality_constraint(opt, user_function, dfc + j,
		     tol ? tol[j] : 0.0) == NLOPT_SUCCESS,
		     "nlopt error adding inequality constraint");
	  }
     }

     /* nonlinear equality constraints */
     mx = mxGetField(prhs[0], 0, "h");
     if (mx) {
	  CHECK(mxIsCell(mx), "h must be a Cell array");
	  m = mxGetNumberOfElements(mx);
	  dh = (user_function_data *) mxCalloc(m, sizeof(user_function_data));
	  tol = struct_arrval(prhs[0], "h_tol", m, NULL);
	  for (j = 0; j < m; ++j) {
	       mxArray *h = cell_funcval(mx, j);
	       if (mxIsChar(h)) {
		    CHECK(mxGetString(h, dh[j].f, NAMELENGTHMAX) == 0,
		     "error reading function name string (too long?)");
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
	       CHECK(nlopt_add_equality_constraint(opt, user_function, dh + j,
		     tol ? tol[j] : 0.0) == NLOPT_SUCCESS,
		     "nlopt error adding equality constraint");
	  }
     }

     /* x = x0 */
     plhs[0] = mxDuplicateArray(prhs[1]);
     x = mxGetPr(plhs[0]);

     /* optimize */
     ret = nlopt_optimize(opt, x, &opt_f);

     /* assign ouput arguments */
     if (nlhs > 1) plhs[1] = mxCreateDoubleScalar(opt_f);
     if (nlhs > 2) plhs[2] = mxCreateDoubleScalar((int) ret);

     /* cleanup */
     mxFree(dh);
     mxFree(dfc);
     mxDestroyArray(d.prhs[d.xrhs]);
     if (dpre.nrhs > 0) mxDestroyArray(dpre.prhs[dpre.xrhs+1]);
     nlopt_destroy(opt);
}
