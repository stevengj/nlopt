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

/* Matlab MEX interface to NLopt, and in particular to nlopt_minimize_constrained */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>

#include "nlopt.h"
#include "nlopt-util.h"

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

#define FLEN 128 /* max length of user function name */
#define MAXRHS 128 /* max nrhs for user function */
typedef struct {
     char f[FLEN];
     mxArray *plhs[2];
     mxArray *prhs[MAXRHS];
     int xrhs, nrhs;
     int verbose, neval;
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
  d->neval++;
  if (d->verbose) mexPrintf("nlopt_minimize_constrained eval #%d: %g\n", d->neval, f);
  return f;
}				 

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
     nlopt_algorithm algorithm;
     int n, m, i, j;
     double *lb, *ub, *x, *x0;
     double minf_max, ftol_rel, ftol_abs, xtol_rel, *xtol_abs, maxtime;
     int maxeval;
     nlopt_result ret;
     mxArray *x_mx;
     double minf = HUGE_VAL;
     user_function_data d, *dc;

     CHECK(nrhs == 9 && nlhs <= 3, "wrong number of arguments");

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
	  d.prhs[0] = (mxArray *) prhs[1];
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

     /* m = length(fc = prhs[3]) = length(fc_data = prhs[4])  */
     CHECK(mxIsCell(prhs[3]), "fc must be a Cell array");
     CHECK(mxIsCell(prhs[4]), "fc_data must be a Cell array");
     m = mxGetM(prhs[3]) * mxGetN(prhs[3]);
     CHECK(m == mxGetM(prhs[4]) * mxGetN(prhs[4]), "fc and fc_data must have the same length");
     dc = (user_function_data *) malloc(sizeof(user_function_data) * m);

     for (j = 0; j < m; ++j) {
	  mxArray *fc, *fc_data;

	  /* function fc = phrs[3] */
	  fc = mxGetCell(prhs[3], j);
	  CHECK(mxIsChar(fc) || mxIsFunctionHandle(fc),
		"fc must be Cell array of function handles or function names");
	  if (mxIsChar(fc)) {
	       CHECK(mxGetString(fc, dc[j].f, FLEN) == 0,
		     "error reading function name string (too long?)");
	       dc[j].nrhs = 1;
	       dc[j].xrhs = 0;
	  }
	  else {
	       dc[j].prhs[0] = fc;
	       strcpy(dc[j].f, "feval");
	       dc[j].nrhs = 2;
	       dc[j].xrhs = 1;
	  }
	  
	  /* Cell fc_data = prhs[4] */
	  fc_data = mxGetCell(prhs[4], j);
	  CHECK(mxIsCell(fc_data), "fc_data must be a Cell array of Cell arrays");
	  CHECK(mxGetM(fc_data) * mxGetN(fc_data) + 1 <= MAXRHS,
		"user function cannot have more than " STRIZE(MAXRHS) " arguments");
	  dc[j].nrhs += mxGetM(fc_data) * mxGetN(fc_data);
	  for (i = 0; i < dc[j].nrhs - (1+dc[j].xrhs); ++i)
	       dc[j].prhs[(1+dc[j].xrhs)+i] = mxGetCell(fc_data, i);
     }

     /* lb = prhs[5] */
     CHECK(mxIsDouble(prhs[5]) && !mxIsComplex(prhs[5])
	   && (mxGetM(prhs[5]) == 1 || mxGetN(prhs[5]) == 1),
	   "lb must be real row or column vector");
     lb = mxGetPr(prhs[5]);
     n = mxGetM(prhs[5]) * mxGetN(prhs[5]);

     /* ub = prhs[6] */
     CHECK(mxIsDouble(prhs[6]) && !mxIsComplex(prhs[6])
	   && (mxGetM(prhs[6]) == 1 || mxGetN(prhs[6]) == 1)
	   && mxGetM(prhs[6]) * mxGetN(prhs[6]) == n,
	   "ub must be real row or column vector of same length as lb");
     ub = mxGetPr(prhs[6]);

     /* x0 = prhs[7] */
     CHECK(mxIsDouble(prhs[7]) && !mxIsComplex(prhs[7])
	   && (mxGetM(prhs[7]) == 1 || mxGetN(prhs[7]) == 1)
	   && mxGetM(prhs[7]) * mxGetN(prhs[7]) == n,
	   "x must be real row or column vector of same length as lb");
     x0 = mxGetPr(prhs[7]);

     /* stopping criteria = prhs[8] */
     CHECK(mxIsStruct(prhs[8]), "stopping criteria must be a struct");
     minf_max = struct_val_default(prhs[8], "minf_max", -HUGE_VAL);
     ftol_rel = struct_val_default(prhs[8], "ftol_rel", 0);
     ftol_abs = struct_val_default(prhs[8], "ftol_abs", 0);
     xtol_rel = struct_val_default(prhs[8], "xtol_rel", 0);
     maxeval = (int) (struct_val_default(prhs[8], "maxeval", -1) + 0.5);
     maxtime = struct_val_default(prhs[8], "maxtime", -1);
     d.verbose = (int) struct_val_default(prhs[8], "verbose", 0);
     d.neval = 0;
     for (i = 0; i < m; ++i) {
	  dc[i].verbose = d.verbose > 1;
	  dc[i].neval = 0;
     }
     {
	  mxArray *val = mxGetField(prhs[8], 0, "xtol_abs");
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
     for (i = 0; i < m;++i)
	  dc[i].prhs[dc[i].xrhs] = d.prhs[d.xrhs];
     
     ret = nlopt_minimize_constrained(algorithm,
				      n,
				      user_function, &d,
				      m, user_function, dc,
				      sizeof(user_function_data),
				      lb, ub, x, &minf, minf_max, 
				      ftol_rel, ftol_abs, xtol_rel, xtol_abs,
				      maxeval, maxtime);

     mxDestroyArray(d.prhs[d.xrhs]);
     free(dc);

     plhs[0] = x_mx;
     if (nlhs > 1) {
	  plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	  *(mxGetPr(plhs[1])) = minf;
     }
     if (nlhs > 2) {
	  plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
	  *(mxGetPr(plhs[2])) = (int) ret;
     }
}
