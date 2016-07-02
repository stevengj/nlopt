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

#define CHECK1(cond, msg)               \
    if (!(cond)) {                      \
        nlopt_destroy(opt);             \
        mexWarnMsgIdAndTxt(ERRID, msg); \
        return NULL;                    \
    }

#define CHECK2(cond, msg)              \
    if (!(cond)) {                     \
        mxFree(dh);                    \
        mxFree(dfc);                   \
        nlopt_destroy(opt);            \
        mexErrMsgIdAndTxt(ERRID, msg); \
    }

#define CHECK3(ret, msg)                          \
    if ((ret) != NLOPT_SUCCESS) {                 \
        const char *buf = ret_msg(ret, opt, msg); \
        mxFree(dh);                               \
        mxFree(dfc);                              \
        nlopt_destroy(opt);                       \
        mexErrMsgIdAndTxt(ERRID, buf);            \
    }

#define BUFSIZE 256      /* length of buffer used for error messages */
#define NAMELENGTHMAX 64 /* TMW_NAME_LENGTH_MAX: max length of varname */

typedef struct user_function_data_s {
    char func[NAMELENGTHMAX];
    char type[NAMELENGTHMAX];
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
    return (arr && mxIsDouble(arr) && !mxIsComplex(arr) && !mxIsSparse(arr) &&
        mxIsScalar(arr));
}

static bool mx_isvector(const mxArray *arr)
{
    return (arr && mxIsDouble(arr) && !mxIsComplex(arr) && !mxIsSparse(arr) &&
        mxGetNumberOfDimensions(arr) == 2 &&
        (mxGetM(arr) == 1 || mxGetN(arr) == 1));
}

static bool mx_isvector_len(const mxArray *arr, unsigned n)
{
    return (mx_isvector(arr) && mxGetNumberOfElements(arr) == n);
}

static bool mx_ismatrix_len(const mxArray *arr, unsigned r, unsigned c)
{
    return (arr && mxIsDouble(arr) && !mxIsComplex(arr) && !mxIsSparse(arr) &&
        mxGetNumberOfDimensions(arr) == 2 &&
        mxGetNumberOfElements(arr) == r*c &&
        mxGetM(arr) == r && mxGetN(arr) == c);
}

static bool mx_isfunction(const mxArray *arr)
{
    return arr && (mxIsFunctionHandle(arr) || (mxIsChar(arr) &&
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

static double arr_norm_inf(const double *x, unsigned n)
{
	unsigned i;
	double mx = 0.0;
	for (i = 0; i < n; ++i) {
		if (fabs(x[i]) > mx)
			mx = fabs(x[i]);
	}
	return mx;
}

static double arr_max(const double *x, unsigned n)
{
	unsigned i;
	double mx = -HUGE_VAL;
	for (i = 0; i < n; ++i) {
		if (x[i] > mx)
			mx = x[i];
	}
	return mx;
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
  if (mexCallMATLAB(grad ? 2 : 1, d->plhs, d->nrhs, d->prhs, d->func) != 0)
	mexErrMsgIdAndTxt(ERRID, "error calling %s function", d->type);

  /* f */
  if (!mx_isscalar(d->plhs[0]))
	mexErrMsgIdAndTxt(ERRID, "%s function must return a real scalar", d->type);
  f = mxGetScalar(d->plhs[0]);
  mxDestroyArray(d->plhs[0]);

  /* g */
  if (grad) {
     if (!mx_isvector_len(d->plhs[1], n))
	   mexErrMsgIdAndTxt(ERRID,
	     "%s function must return a gradient vector of length %u", d->type, n);
     memcpy(grad, mxGetPr(d->plhs[1]), n * sizeof(double));
     mxDestroyArray(d->plhs[1]);
  }

  d->neval++;
  if (d->verbose) {
  	mexPrintf("%13s eval #%-3d: f(x) = %12g", d->type, d->neval, f);
  	if (grad) mexPrintf(", norm(g(x)) = %12g", arr_norm_inf(grad,n));
  	mexPrintf("\n");
  }
  if (mxIsNaN(f)) nlopt_force_stop(d->opt);
  return f;
}

static void user_mfunction(unsigned m, double *result,
		       unsigned n, const double *x, double *grad, void *data)
{
	user_function_data *d = (user_function_data *) data;

	/* x */
	memcpy(mxGetPr(d->prhs[d->xrhs]), x, n * sizeof(double));

	/* [result, grad] = feval(conFunc, x) */
	d->plhs[0] = d->plhs[1] = NULL;
	if (mexCallMATLAB(grad ? 2 : 1, d->plhs, d->nrhs, d->prhs, d->func) != 0)
		mexErrMsgIdAndTxt(ERRID,
			"error calling %s constraint function", d->type);

	/* result */
	if (!mx_isvector_len(d->plhs[0], m))
		mexErrMsgIdAndTxt(ERRID,
			"%s constraint function must return a vector of length %u",
			d->type, m);
	memcpy(result, mxGetPr(d->plhs[0]), m * sizeof(double));
	mxDestroyArray(d->plhs[0]);

	/* grad */
	/* NOTE: MATLAB arrays are column-major while C is row-major, so we ask
		the user function to return a transposed jacobian matrix n-by-m,
		so that memcpy() works directly. */
	if (grad) {
		if (!mx_ismatrix_len(d->plhs[1], n, m))
			mexErrMsgIdAndTxt(ERRID,
				"%s constraint function must return a jacobian matrix of size (%u,%u)",
				d->type, n, m);
		memcpy(grad, mxGetPr(d->plhs[1]), n * m * sizeof(double));
		mxDestroyArray(d->plhs[1]);
	}

	d->neval++;
	if (d->verbose) {
		mexPrintf("%13s eval #%-3d: max(c(x)) = %12g\n",
			d->type, d->neval, arr_max(result,m));
	}
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
  if (mexCallMATLAB(1, d->plhs, d->nrhs, d->prhs, d->func) != 0)
	 mexErrMsgIdAndTxt(ERRID, "error calling %s function", d->type);

  /* vpre = H(x)*v */
  if (!mx_isvector_len(d->plhs[0], n))
	 mexErrMsgIdAndTxt(ERRID,
	 	"%s function must return a vpre vector of length %u", d->type, n);
  memcpy(vpre, mxGetPr(d->plhs[0]), n * sizeof(double));
  mxDestroyArray(d->plhs[0]);

  d->neval++;
  if (d->verbose) mexPrintf("%13s eval #%-3d\n", d->type, d->neval);
}

static nlopt_opt make_opt(const mxArray *s, unsigned n)
{
     nlopt_opt opt = NULL;
     nlopt_algorithm alg;
     nlopt_result ret;
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
     ret = (tmp) ? nlopt_set_lower_bounds(opt, tmp) :
         nlopt_set_lower_bounds1(opt, -HUGE_VAL);
     CHECK1(ret, "error setting opt.lower_bounds");
     tmp = struct_arrval(s, "upper_bounds", n, NULL);
     ret = (tmp) ? nlopt_set_upper_bounds(opt, tmp) :
         nlopt_set_upper_bounds1(opt, +HUGE_VAL);
     CHECK1(ret, "error setting opt.upper_bounds");

     /* stopping criteria */
     ret = nlopt_set_stopval(opt, struct_val(s, "stopval", -HUGE_VAL));
     CHECK1(ret, "error setting opt.stopval");
     ret = nlopt_set_ftol_rel(opt, struct_val(s, "ftol_rel", 0.0));
     CHECK1(ret, "error setting opt.ftol_rel");
     ret = nlopt_set_ftol_abs(opt, struct_val(s, "ftol_abs", 0.0));
     CHECK1(ret, "error setting opt.ftol_abs");
     ret = nlopt_set_xtol_rel(opt, struct_val(s, "xtol_rel", 0.0));
     CHECK1(ret, "error setting opt.xtol_rel");
     tmp = struct_arrval(s, "xtol_abs", n, NULL);
     ret = (tmp) ? nlopt_set_xtol_abs(opt, tmp) :
         nlopt_set_xtol_abs1(opt, 0.0);
     CHECK1(ret, "error setting opt.xtol_abs");
     ret = nlopt_set_maxeval(opt, (int) struct_val(s, "maxeval", 0));
     CHECK1(ret, "error setting opt.maxeval");
     ret = nlopt_set_maxtime(opt, struct_val(s, "maxtime", 0.0));
     CHECK1(ret, "error setting opt.maxtime");

     /* stochastic population */
     ret = nlopt_set_population(opt,
    	(unsigned) struct_val(s, "population", 0));
     CHECK1(ret, "error setting opt.population");

     /* vector storage */
     ret = nlopt_set_vector_storage(opt,
         (unsigned) struct_val(s, "vector_storage", 0));
     CHECK1(ret, "error setting opt.vector_storage");

     /* initial step size */
     ret = nlopt_set_initial_step(opt,
         struct_arrval(s, "initial_step", n, NULL));
     CHECK1(ret, "error setting opt.initial_step");

     /* local optimization algorithm */
	 s_local = mxGetField(s, 0, "local_optimizer");
     if (s_local) {
	  nlopt_opt opt_local = NULL;
	  CHECK1(mxIsStruct(s_local) && mxIsScalar(s_local),
		 "opt.local_optimizer must be a scalar structure");
	  opt_local = make_opt(s_local, n);
	  CHECK1(opt_local, "error initializing local optimizer options");
	  ret = nlopt_set_local_optimizer(opt, opt_local);
	  nlopt_destroy(opt_local);
	  CHECK1(ret, "error setting opt.local_optimizer");
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

static const char *ret_msg(nlopt_result ret, const nlopt_opt opt,
    const char *msg)
{
    static char buffer[BUFSIZE];
    memset(buffer, 0, BUFSIZE);
    sprintf(buffer, "%s: %s\n%s", translate_result(ret),
        (nlopt_get_errmsg(opt) ? nlopt_get_errmsg(opt) : "\b\b"),
        (msg ? msg : "\b"));
    return buffer;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
     unsigned n, j, m;
     double *x, opt_f, *tol;
     nlopt_result ret;
     mxArray *mx;
     user_function_data d = {0}, dpre = {0}, dc = {0}, dceq = {0},
         *dfc = NULL, *dh = NULL;
     nlopt_opt opt = NULL;
     bool ismin;

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
     ismin = (mx != NULL);
     if (!mx) mx = struct_funcval(prhs[0], "max_objective");
     CHECK2(mx, "opt.(min|max)_objective must be set");
     if (mxIsChar(mx)) {
	  CHECK2(mxGetString(mx, d.func, NAMELENGTHMAX) == 0,
		"error reading function name string from opt.(min|max)_objective");
	  d.nrhs = 1;
	  d.xrhs = 0;
     }
     else {
	  d.prhs[0] = mx;
	  strcpy(d.func, "feval");
	  d.nrhs = 2;
	  d.xrhs = 1;
     }
     d.verbose = (int) struct_val(prhs[0], "verbose", 0);
     d.neval = 0;
     d.opt = opt;
     d.prhs[d.xrhs] = mxCreateDoubleMatrix(1, n, mxREAL);
     strcpy(d.type, ismin ? "min_objective": "max_objective");

     /* preconditioned objective function */
     mx = struct_funcval(prhs[0], "pre");
     if (!mx) {
	  dpre.nrhs = 0;  /* indicate no memory allocation */
	  ret = (ismin) ?
	       nlopt_set_min_objective(opt, user_function, &d) :
	       nlopt_set_max_objective(opt, user_function, &d);
     }
     else {
	  if (mxIsChar(mx)) {
	       CHECK2(mxGetString(mx, dpre.func, NAMELENGTHMAX) == 0,
                     "error reading function name string from opt.pre");
	       dpre.nrhs = 2;
	       dpre.xrhs = 0;
	  }
	  else {
	       dpre.prhs[0] = mx;
	       strcpy(dpre.func, "feval");
	       dpre.nrhs = 3;
	       dpre.xrhs = 1;
	  }
	  dpre.verbose = d.verbose > 2;
	  dpre.neval = 0;
	  dpre.opt = opt;
	  dpre.prhs[dpre.xrhs] = d.prhs[d.xrhs];
	  dpre.prhs[dpre.xrhs+1] = mxCreateDoubleMatrix(1, n, mxREAL);
	  strcpy(dpre.type, "precond");
	  d.dpre = &dpre;
	  ret = (ismin) ?
	       nlopt_set_precond_min_objective(opt, user_function, user_pre, &d) :
	       nlopt_set_precond_max_objective(opt, user_function, user_pre, &d);
     }
     CHECK3(ret, "error setting objective function");

     /* nonlinear inequality constraints */
     mx = mxGetField(prhs[0], 0, "fc");
     if (mx) {
	  CHECK2(mxIsCell(mx), "opt.fc must be a cell array");
	  m = mxGetNumberOfElements(mx);
	  dfc = (user_function_data *) mxCalloc(m, sizeof(user_function_data));
	  tol = struct_arrval(prhs[0], "fc_tol", m, NULL);
	  for (j = 0; j < m; ++j) {
	       mxArray *fc = cell_funcval(mx, j);
	       if (mxIsChar(fc)) {
		    CHECK2(mxGetString(fc, dfc[j].func, NAMELENGTHMAX) == 0,
		     "error reading function name string from opt.fc");
		    dfc[j].nrhs = 1;
		    dfc[j].xrhs = 0;
	       }
	       else {
		    dfc[j].prhs[0] = fc;
		    strcpy(dfc[j].func, "feval");
		    dfc[j].nrhs = 2;
		    dfc[j].xrhs = 1;
	       }
	       dfc[j].verbose = d.verbose > 1;
	       dfc[j].neval = 0;
	       dfc[j].opt = opt;
	       dfc[j].prhs[dfc[j].xrhs] = d.prhs[d.xrhs];
	       sprintf(dfc[j].type, "fc{%u}", j);
	       ret = nlopt_add_inequality_constraint(opt, user_function, dfc + j,
		     tol ? tol[j] : 0.0);
	       CHECK3(ret, "error adding inequality constraint opt.fc");
	  }
     }

     /* nonlinear equality constraints */
     mx = mxGetField(prhs[0], 0, "h");
     if (mx) {
	  CHECK2(mxIsCell(mx), "opt.h must be a cell array");
	  m = mxGetNumberOfElements(mx);
	  dh = (user_function_data *) mxCalloc(m, sizeof(user_function_data));
	  tol = struct_arrval(prhs[0], "h_tol", m, NULL);
	  for (j = 0; j < m; ++j) {
	       mxArray *h = cell_funcval(mx, j);
	       if (mxIsChar(h)) {
		    CHECK2(mxGetString(h, dh[j].func, NAMELENGTHMAX) == 0,
		     "error reading function name string from opt.h");
		    dh[j].nrhs = 1;
		    dh[j].xrhs = 0;
	       }
	       else {
		    dh[j].prhs[0] = h;
		    strcpy(dh[j].func, "feval");
		    dh[j].nrhs = 2;
		    dh[j].xrhs = 1;
	       }
	       dh[j].verbose = d.verbose > 1;
	       dh[j].neval = 0;
	       dh[j].opt = opt;
	       dh[j].prhs[dh[j].xrhs] = d.prhs[d.xrhs];
	       sprintf(dh[j].type, "h{%u}", j);
	       ret = nlopt_add_equality_constraint(opt, user_function, dh + j,
		     tol ? tol[j] : 0.0);
	       CHECK3(ret, "error adding equality constraint opt.h");
	  }
     }

     /* vector-valued nonlinear inequality constraints */
     mx = struct_funcval(prhs[0], "c");
     if (mx) {
        const mxArray *mx_tol = mxGetField(prhs[0], 0, "c_tol");
        CHECK2(mx_tol && mx_isvector(mx_tol),
            "opt.c_tol must be a real vector of length m (dim of opt.c)");
        m = (unsigned) mxGetNumberOfElements(mx_tol);
        tol = mxGetPr(mx_tol);
        if (mxIsChar(mx)) {
            CHECK2(mxGetString(mx, dc.func, NAMELENGTHMAX) == 0,
                "error reading function name string from opt.c");
            dc.nrhs = 1;
            dc.xrhs = 0;
        }
        else {
            dc.prhs[0] = mx;
            strcpy(dc.func, "feval");
            dc.nrhs = 2;
            dc.xrhs = 1;
        }
        dc.verbose = d.verbose > 1;
        dc.neval = 0;
        dc.opt = opt;
        dc.prhs[dc.xrhs] = d.prhs[d.xrhs];
        strcpy(dc.type, "c");
        ret = nlopt_add_inequality_mconstraint(
            opt, m, user_mfunction, &dc, tol);
        CHECK3(ret, "error adding vector-valued inequality constraint opt.c");
     }

     /* vector-valued nonlinear equality constraints */
     mx = struct_funcval(prhs[0], "ceq");
     if (mx) {
        const mxArray *mx_tol = mxGetField(prhs[0], 0, "ceq_tol");
        CHECK2(mx_tol && mx_isvector(mx_tol),
            "opt.ceq_tol must be a real vector of length m (dim of opt.c)");
        m = (unsigned) mxGetNumberOfElements(mx_tol);
        tol = mxGetPr(mx_tol);
        if (mxIsChar(mx)) {
            CHECK2(mxGetString(mx, dceq.func, NAMELENGTHMAX) == 0,
                "error reading function name string from opt.ceq");
            dceq.nrhs = 1;
            dceq.xrhs = 0;
        }
        else {
            dceq.prhs[0] = mx;
            strcpy(dceq.func, "feval");
            dceq.nrhs = 2;
            dceq.xrhs = 1;
        }
        dceq.verbose = d.verbose > 1;
        dceq.neval = 0;
        dceq.opt = opt;
        dceq.prhs[dceq.xrhs] = d.prhs[d.xrhs];
        strcpy(dceq.type, "ceq");
        ret = nlopt_add_equality_mconstraint(
            opt, m, user_mfunction, &dceq, tol);
        CHECK3(ret, "error adding vector-valued equality constraint opt.ceq");
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
            mxCreateString(ret_msg(ret, opt, NULL)));
     }

     /* cleanup */
     mxFree(dh);
     mxFree(dfc);
     mxDestroyArray(d.prhs[d.xrhs]);
     if (dpre.nrhs > 0) mxDestroyArray(dpre.prhs[dpre.xrhs+1]);
     nlopt_destroy(opt);
}
