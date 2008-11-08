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

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "neldermead.h"
#include "redblack.h"

/* Nelder-Mead simplex algorithm, used as a subroutine for the Rowan's
   subplex algorithm.  Modified to handle bound constraints ala
   Richardson and Kuester (1973), as mentioned below. */

/* heuristic "strategy" constants: */
static const double alpha = 1, beta = 0.5, gamm = 2, delta = 0.5;

/* sort order in red-black tree: keys [f(x), x] are sorted by f(x) */
static int simplex_compare(double *k1, double *k2)
{
     if (*k1 < *k2) return -1;
     if (*k1 > *k2) return +1;
     return k1 - k2; /* tie-breaker */
}

/* pin x to lie within the given lower and upper bounds; this is
   probably a suboptimal way to handle bound constraints, but I don't
   know a better way and it is the method suggested by J. A. Richardson
   and J. L. Kuester, "The complex method for constrained optimization,"
   Commun. ACM 16(8), 487-489 (1973). */
static void pin(int n, double *x, const double *lb, const double *ub) {
     int i;
     for (i = 0; i < n; ++i) {
	  if (x[i] < lb[i]) x[i] = lb[i];
	  else if (x[i] > ub[i]) x[i] = ub[i];
     }
}


#define CHECK_EVAL(xc,fc) 						  \
 if ((fc) <= *minf) {							  \
   *minf = (fc); memcpy(x, (xc), n * sizeof(double));			  \
   if (*minf < stop->minf_max) { ret=NLOPT_MINF_MAX_REACHED; goto done; } \
 }									  \
 stop->nevals++;							  \
 if (nlopt_stop_evals(stop)) { ret=NLOPT_MAXEVAL_REACHED; goto done; }	  \
 if (nlopt_stop_time(stop)) { ret=NLOPT_MAXTIME_REACHED; goto done; }

nlopt_result nldrmd_minimize(int n, nlopt_func f, void *f_data,
			     const double *lb, const double *ub, /* bounds */
			     double *x, /* in: initial guess, out: minimizer */
			     double *minf,
			     const double *xstep, /* initial step sizes */
			     nlopt_stopping *stop)
{
     double *pts; /* (n+1) x (n+1) array of n+1 points plus function val [0] */
     double *c; /* centroid * n */
     double *xcur; /* current point */
     rb_tree t; /* red-black tree of simplex, sorted by f(x) */
     int i, j;
     double ninv = 1.0 / n;
     nlopt_result ret = NLOPT_SUCCESS;

     pts = (double*) malloc(sizeof(double) * (n+1) * (n+1));
     if (!pts) return NLOPT_OUT_OF_MEMORY;
     c = (double*) malloc(sizeof(double) * n * 2);
     if (!c) { free(pts); return NLOPT_OUT_OF_MEMORY; }
     xcur = c + n;

     /* initialize the simplex based on the starting xstep */
     memcpy(pts+1, x, sizeof(double)*n);
     *minf = pts[0] = f(n, pts+1, NULL, f_data);
     CHECK_EVAL(pts+1, pts[0]);
     for (i = 0; i < n; ++i) {
	  double *pt = pts + (i+1)*(n+1);
	  memcpy(pt+1, x, sizeof(double)*n);
	  if (x[i] + xstep[i] <= ub[i])
	       pt[1+i] += xstep[i];
	  else if (x[i] - xstep[i] >= lb[i])
	       pt[1+i] -= xstep[i];
	  else if (ub[i] - x[i] > x[i] - lb[i])
	       pt[1+i] += 0.5 * (ub[i] - x[i]);
	  else
	       pt[1+i] -= 0.5 * (x[i] - lb[i]);
	  pt[0] = f(n, pt+1, NULL, f_data);
	  CHECK_EVAL(pt+1, pt[0]);
     }

     rb_tree_init(&t, simplex_compare);
 restart:
     for (i = 0; i < n + 1; ++i)
	  if (!rb_tree_insert(&t, pts + i*(n+1))) {
	       ret = NLOPT_OUT_OF_MEMORY;
	       goto done;
	  }

     while (1) {
	  rb_node *low = rb_tree_min(&t);
	  rb_node *high = rb_tree_max(&t);
	  double fl = low->k[0], *xl = low->k + 1;
	  double fh = high->k[0], *xh = high->k + 1;
	  double fr;

	  if (nlopt_stop_f(stop, fl, fh)) ret = NLOPT_FTOL_REACHED;

	  /* compute centroid ... if we cared about the perfomance of this,
	     we could do it iteratively by updating the centroid on
	     each step, but then we would have to be more careful about
	     accumulation of rounding errors... anyway n is unlikely to
	     be very large for Nelder-Mead in practical cases */
	  memset(c, 0, sizeof(double)*n);
	  for (i = 0; i < n + 1; ++i) {
	       double *xi = pts + i*(n+1) + 1;
	       if (xi != xh)
		    for (j = 0; j < n; ++j)
			 c[j] += xi[j];
	  }
	  for (i = 0; i < n; ++i) c[i] *= ninv;

	  /* x convergence check: find xcur = max radius from centroid */
	  memset(xcur, 0, sizeof(double)*n);
	  for (i = 0; i < n + 1; ++i) {
               double *xi = pts + i*(n+1) + 1;
	       for (j = 0; j < n; ++j) {
		    double dx = fabs(xi[j] - c[j]);
		    if (dx > xcur[j]) xcur[j] = dx;
	       }
	  }
	  for (i = 0; i < n; ++i) xcur[i] += c[i];
	  if (nlopt_stop_x(stop, c, xcur)) ret = NLOPT_XTOL_REACHED;

	  /* reflection */
	  for (i = 0; i < n; ++i) xcur[i] = c[i] + alpha * (c[i] - xh[i]);
	  pin(n, xcur, lb, ub);
	  fr = f(n, xcur, NULL, f_data);
	  CHECK_EVAL(xcur, fr);

	  if (fr < fl) { /* new best point, expand simplex */
	       for (i = 0; i < n; ++i) xh[i] = c[i] + gamm * (c[i] - xh[i]);
	       pin(n, xh, lb, ub);
	       fh = f(n, xh, NULL, f_data);
	       CHECK_EVAL(xh, fh);
	       if (fh >= fr) { /* expanding didn't improve */
		    fh = fr;
		    memcpy(xh, xcur, sizeof(double)*n);
	       }
	  }
	  else if (fr < rb_tree_pred(high)->k[0]) { /* accept new point */
	       memcpy(xh, xcur, sizeof(double)*n);
	       fh = fr;
	  }
	  else { /* new worst point, contract */
	       double fc;
	       if (fh <= fr)
		    for (i = 0; i < n; ++i) xcur[i] = c[i] - beta*(c[i]-xh[i]);
	       else
		    for (i = 0; i < n; ++i) xcur[i] = c[i] + beta*(c[i]-xh[i]);
	       pin(n, xcur, lb, ub);
	       fc = f(n, xcur, NULL, f_data);
	       CHECK_EVAL(xcur, fc);
	       if (fc < fr && fc < fh) { /* successful contraction */
		    memcpy(xh, xcur, sizeof(double)*n);
		    fh = fc;
	       }
	       else { /* failed contraction, shrink simplex */
		    rb_tree_destroy(&t);
		    rb_tree_init(&t, simplex_compare);
		    for (i = 0; i < n+1; ++i) {
			 double *pt = pts + i * (n+1);
			 if (pt+1 != xl) {
			      for (j = 0; j < n; ++j)
				   pt[1+j] = xl[j] + delta * (pt[1+j] - xl[j]);
			      pt[0] = f(n, pt+1, NULL, f_data);
			      CHECK_EVAL(pt+1, pt[0]);
			 }
		    }
		    goto restart;
	       }
	  }

	  high->k[0] = fh;
	  rb_tree_resort(&t, high);
     }
     
done:
     rb_tree_destroy(&t);
     free(c);
     free(pts);
     return ret;
}
