/* Matlab MEX interface to NLopt, and in particular to nlopt_minimize */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>

#include "nlopt.h"

#define xSTRIZE(x) #x
#define STRIZE(x) xSTRIZE(x)
#define CHECK(cond, msg) if (!(cond)) mexErrMsgTxt(msg);

static double struct_val_default(const mxArray *s, const char *name, double dflt)
{
     mxArray *val = mxGetField(s, 0, name);
     if (val) {
	  CHECK(mxIsNumeric(val) && !mxIsComplex(val) 
		&& mxGetM(val) * mxGetN(val) == 1,
		"stop fields, other than xtol_abs, must be real scalars");
	  return mxGetScalar(val);
     }
     return dflt;
}

#define FLEN 1024 /* max length of user function name */
#define MAXRHS 1024 /* max nrhs for user function */
typedef struct {
     char f[FLEN];
     mxArray *plhs[2];
     mxArray *prhs[MAXRHS];
     int xrhs, nrhs;
} user_function_data;

static double user_function(int n, const double *x,
			    double *gradient, /* NULL if not needed */
			    void *d_)
{
  user_function_data *d = (user_function_data *) d_;
  double f;

  d->plhs[0] = d->plhs[1] = NULL;
  memcpy(mxGetPr(d->prhs[d->xrhs]), x, n * sizeof(double));

  CHECK(0 == mexCallMATLAB(gradient ? 2 : 1, d->plhs, 
			   d->nrhs, d->prhs, d->f),
	"error calling user function");

  CHECK(mxIsNumeric(d->plhs[0]) && !mxIsComplex(d->plhs[0]) 
	&& mxGetM(d->plhs[0]) * mxGetN(d->plhs[0]) == 1,
	"user function must return real scalar");
  f = mxGetScalar(d->plhs[0]);
  mxDestroyArray(d->plhs[0]);
  if (gradient) {
     CHECK(mxIsDouble(d->plhs[1]) && !mxIsComplex(d->plhs[1])
	   && (mxGetM(d->plhs[1]) == 1 || mxGetN(d->plhs[1]) == 1)
	   && mxGetM(d->plhs[1]) * mxGetN(d->plhs[1]) == n,
	   "gradient vector from user function is the wrong size");
     memcpy(gradient, mxGetPr(d->plhs[1]), n * sizeof(double));
     mxDestroyArray(d->plhs[1]);
  }
  return f;
}				 

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
     nlopt_algorithm algorithm;
     int n, i;
     double *lb, *ub, *x, *x0;
     double minf_max, ftol_rel, ftol_abs, xtol_rel, *xtol_abs, maxtime;
     int maxeval;
     nlopt_result ret;
     mxArray *x_mx;
     double minf = HUGE_VAL;
     user_function_data d;

     CHECK(nrhs == 7 && nlhs <= 3, "wrong number of arguments");

     /* algorithm = prhs[0] */
     CHECK(mxIsNumeric(prhs[0]) && !mxIsComplex(prhs[0]) 
	   && mxGetM(prhs[0]) * mxGetN(prhs[0]) == 1,
	   "algorithm must be real (integer) scalar");
     algorithm = (nlopt_algorithm) (mxGetScalar(prhs[0]) + 0.5);
     CHECK(algorithm >= 0 && algorithm < NLOPT_NUM_ALGORITHMS,
	   "unknown algorithm");

     /* function f = prhs[1] */
     CHECK(mxIsChar(prhs[1]) || mxIsFunctionHandle(prhs[1]), 
	   "f must be a function handle or function name");
     if (mxIsChar(prhs[1])) {
	  CHECK(mxGetString(prhs[1], d.f, FLEN) == 0,
		"error reading function name string (too long?)");
	  d.nrhs = 1;
	  d.xrhs = 0;
     }
     else {
	  d.prhs[0] = prhs[1];
	  strcpy(d.f, "feval");
	  d.nrhs = 2;
	  d.xrhs = 1;
     }
     
     /* Cell f_data = prhs[2] */
     CHECK(mxIsCell(prhs[2]), "f_data must be a Cell array");
     CHECK(mxGetM(prhs[2]) * mxGetN(prhs[2]) + 1 <= MAXRHS,
	   "user function cannot have more than " STRIZE(MAXRHS) " arguments");
     d.nrhs += mxGetM(prhs[2]) * mxGetN(prhs[2]);
     for (i = 0; i < d.nrhs - (1+d.xrhs); ++i)
	  d.prhs[(1+d.xrhs)+i] = mxGetCell(prhs[2], i);

     /* lb = prhs[3] */
     CHECK(mxIsDouble(prhs[3]) && !mxIsComplex(prhs[3])
	   && (mxGetM(prhs[3]) == 1 || mxGetN(prhs[3]) == 1),
	   "lb must be real row or column vector");
     lb = mxGetPr(prhs[3]);
     n = mxGetM(prhs[3]) * mxGetN(prhs[3]);

     /* ub = prhs[4] */
     CHECK(mxIsDouble(prhs[4]) && !mxIsComplex(prhs[4])
	   && (mxGetM(prhs[4]) == 1 || mxGetN(prhs[4]) == 1)
	   && mxGetM(prhs[4]) * mxGetN(prhs[4]) == n,
	   "ub must be real row or column vector of same length as lb");
     ub = mxGetPr(prhs[4]);

     /* x0 = prhs[5] */
     CHECK(mxIsDouble(prhs[5]) && !mxIsComplex(prhs[5])
	   && (mxGetM(prhs[5]) == 1 || mxGetN(prhs[5]) == 1)
	   && mxGetM(prhs[5]) * mxGetN(prhs[5]) == n,
	   "x must be real row or column vector of same length as lb");
     x0 = mxGetPr(prhs[5]);

     /* stopping criteria = prhs[6] */
     CHECK(mxIsStruct(prhs[6]), "stopping criteria must be a struct");
     minf_max = struct_val_default(prhs[6], "minf_max", -HUGE_VAL);
     ftol_rel = struct_val_default(prhs[6], "ftol_rel", 0);
     ftol_abs = struct_val_default(prhs[6], "ftol_abs", 0);
     xtol_rel = struct_val_default(prhs[6], "xtol_rel", 0);
     maxeval = (int) (struct_val_default(prhs[6], "maxeval", -1) + 0.5);
     maxtime = struct_val_default(prhs[6], "maxtime", -1);
     {
	  mxArray *val = mxGetField(prhs[6], 0, "xtol_abs");
	  if (val) {
	       CHECK(mxIsNumeric(val) && !mxIsComplex(val) 
		     && (mxGetM(val) == 1 || mxGetN(val) == 1)
		     && mxGetM(val) * mxGetN(val) == n,
		     "stop.xtol_abs must be real row/col vector of length n");
	       xtol_abs = mxGetPr(val);
	  }
	  else
	       xtol_abs = NULL;
     }


     x_mx = mxCreateDoubleMatrix(1, n, mxREAL);
     x = mxGetPr(x_mx);
     memcpy(x, x0, sizeof(double) * n);

     d.prhs[d.xrhs] = mxCreateDoubleMatrix(1, n, mxREAL);
     
     ret = nlopt_minimize(algorithm,
			  n,
			  user_function, &d,
			  lb, ub, x, &minf,
			  minf_max, ftol_rel, ftol_abs, xtol_rel, xtol_abs,
			  maxeval, maxtime);

     mxDestroyArray(d.prhs[d.xrhs]);

     plhs[0] = x_mx;
     if (nlhs > 1) {
	  plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	  *(mxGetPr(plhs[1])) = minf;
     }
     else
	  mxDestroyArray(d.plhs[0]);
     if (nlhs > 2) {
	  plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
	  *(mxGetPr(plhs[2])) = (int) ret;
     }
}
