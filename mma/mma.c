#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "mma.h"

int mma_verbose = 0; /* > 0 for verbose output */

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* magic minimum value for rho in MMA ... the 2002 paper says it should
   be a "fixed, strictly positive `small' number, e.g. 1e-5"
   ... grrr, I hate these magic numbers, which seem like they
   should depend on the objective function in some way ... in particular,
   note that rho is dimensionful (= dimensions of objective function) */
#define MMA_RHOMIN 1e-5

/***********************************************************************/
/* function for MMA's dual solution of the approximate problem */

typedef struct {
     int count; /* evaluation count, incremented each call */
     int n; /* must be set on input to dimension of x */
     const double *x, *lb, *ub, *sigma, *dfdx; /* arrays of length n */
     const double *dfcdx; /* m-by-n array of fc gradients */
     double fval, rho; /* must be set on input */
     const double *fcval, *rhoc; /* arrays of length m */
     double *xcur; /* array of length n, output each time */
     double gval, wval, *gcval; /* output each time (array length m) */
} dual_data;

static double sqr(double x) { return x * x; }

static double dual_func(int m, const double *y, double *grad, void *d_)
{
     dual_data *d = (dual_data *) d_;
     int n = d->n;
     const double *x = d->x, *lb = d->lb, *ub = d->ub, *sigma = d->sigma, 
	  *dfdx = d->dfdx;
     const double *dfcdx = d->dfcdx;
     double rho = d->rho, fval = d->fval;
     const double *rhoc = d->rhoc, *fcval = d->fcval;
     double *xcur = d->xcur;
     double *gcval = d->gcval;
     int i, j;
     double val;

     d->count++;

     val = d->gval = fval;
     d->wval = 0;
     for (i = 0; i < m; ++i) val += y[i] * (gcval[i] = fcval[i]);

     for (j = 0; j < n; ++j) {
	  double u, v, dx, denominv, c, sigma2, dx2;

	  /* first, compute xcur[j] for y.  Because this objective is
	     separable, we can minimize over x analytically, and the minimum
	     dx is given by the solution of a quadratic equation:
	             u dx^2 + 2 v sigma^2 dx + u sigma^2 = 0
	     where u and v are defined by the sums below.  Because of
	     the definitions, it is guaranteed that |u/v| <= sigma,
	     and it follows that the only dx solution with |dx| <= sigma
	     is given by:
	             (v/u) sigma^2 (-1 + sqrt(1 - (u / v sigma)^2))
             (which goes to zero as u -> 0). */

	  u = dfdx[j];
	  v = fabs(dfdx[j]) * sigma[j] + 0.5 * rho;
	  for (i = 0; i < m; ++i) {
	       u += dfcdx[i*n + j] * y[i];
	       v += (fabs(dfcdx[i*n + j]) * sigma[j] + 0.5 * rhoc[i]) * y[i];
	  }
	  u *= (sigma2 = sqr(sigma[j]));
	  if (fabs(u) < 1e-3 * (v*sigma[j])) { /* Taylor exp. for small u */
	       double a = u / (v*sigma[j]);
	       dx = -sigma[j] * (0.5 * a + 0.125 * a*a*a);
	  }
	  else
	       dx = (v/u)*sigma2 * (-1 + sqrt(1 - sqr(u/(v*sigma[j]))));
	  xcur[j] = x[j] + dx;
	  if (xcur[j] > ub[j]) xcur[j] = ub[j];
	  else if (xcur[j] < lb[j]) xcur[j] = lb[j];
	  if (xcur[j] > x[j]+0.9*sigma[j]) xcur[j] = x[j]+0.9*sigma[j];
	  else if (xcur[j] < x[j]-0.9*sigma[j]) xcur[j] = x[j]-0.9*sigma[j];
	  dx = xcur[j] - x[j];
	  
	  /* function value: */
	  dx2 = dx * dx;
	  denominv = 1.0 / (sigma2 - dx2);
	  val += (u * dx + v * dx2) * denominv;

	  /* update gval, wval, gcval (approximant functions) */
	  c = sigma2 * dx;
	  d->gval += (dfdx[j] * c + (fabs(dfdx[j])*sigma[j] + 0.5*rho) * dx2)
	       * denominv;
	  d->wval += 0.5 * dx2 * denominv;
	  for (i = 0; i < m; ++i)
	       gcval[i] += (dfcdx[i*n+j] * c + (fabs(dfcdx[i*n+j])*sigma[j] 
						+ 0.5*rhoc[i]) * dx2)
		    * denominv;
     }

     /* gradient is easy to compute: since we are at a minimum x (dval/dx=0),
	we only need the partial derivative with respect to y, and
	we negate because we are maximizing: */
     if (grad) for (i = 0; i < m; ++i) grad[i] = -gcval[i];
     return -val;
}

/***********************************************************************/

nlopt_result mma_minimize(int n, nlopt_func f, void *f_data,
			  int m, nlopt_func fc,
			  void *fc_data_, ptrdiff_t fc_datum_size,
			  const double *lb, const double *ub, /* bounds */
			  double *x, /* in: initial guess, out: minimizer */
			  double *minf,
			  nlopt_stopping *stop,
			  nlopt_algorithm dual_alg, 
			  double dual_tolrel, int dual_maxeval)
{
     nlopt_result ret = NLOPT_SUCCESS;
     double *xcur, rho, *sigma, *dfdx, *dfdx_cur, *xprev, *xprevprev, fcur;
     double *dfcdx, *dfcdx_cur;
     double *fcval, *fcval_cur, *rhoc, *gcval, *y, *dual_lb, *dual_ub;
     int i, j, k = 0;
     char *fc_data = (char *) fc_data_;
     dual_data dd;
     int feasible;
     
     sigma = (double *) malloc(sizeof(double) * (6*n + 2*m*n + m*7));
     if (!sigma) return NLOPT_OUT_OF_MEMORY;
     dfdx = sigma + n;
     dfdx_cur = dfdx + n;
     xcur = dfdx_cur + n;
     xprev = xcur + n;
     xprevprev = xprev + n;
     fcval = xprevprev + n;
     fcval_cur = fcval + m;
     rhoc = fcval_cur + m;
     gcval = rhoc + m;
     dual_lb = gcval + m;
     dual_ub = dual_lb + m;
     y = dual_ub + m;
     dfcdx = y + m;
     dfcdx_cur = dfcdx + m*n;

     dd.n = n;
     dd.x = x;
     dd.lb = lb;
     dd.ub = ub;
     dd.sigma = sigma;
     dd.dfdx = dfdx;
     dd.dfcdx = dfcdx;
     dd.fcval = fcval;
     dd.rhoc = rhoc;
     dd.xcur = xcur;
     dd.gcval = gcval;

     for (j = 0; j < n; ++j) {
	  if (nlopt_isinf(ub[j]) || nlopt_isinf(lb[j]))
	       sigma[j] = 1.0; /* arbitrary default */
	  else
	       sigma[j] = 0.5 * (ub[j] - lb[j]);
     }
     rho = 1.0;
     for (i = 0; i < m; ++i) {
	  rhoc[i] = 1.0;
	  dual_lb[i] = y[i] = 0.0;
	  dual_ub[i] = HUGE_VAL;
     }

     dd.fval = fcur = *minf = f(n, x, dfdx, f_data);
     stop->nevals++;
     memcpy(xcur, x, sizeof(double) * n);

     feasible = 1;
     for (i = 0; i < m; ++i) {
	  fcval[i] = fc(n, x, dfcdx + i*n, fc_data + fc_datum_size * i);
	  feasible = feasible && (fcval[i] <= 0);
     }
     if (!feasible) { ret = NLOPT_FAILURE; goto done; } /* TODO: handle this */

     while (1) { /* outer iterations */
	  double fprev = fcur;
	  if (nlopt_stop_evals(stop)) ret = NLOPT_MAXEVAL_REACHED;
	  else if (nlopt_stop_time(stop)) ret = NLOPT_MAXTIME_REACHED;
	  else if (*minf < stop->minf_max) ret = NLOPT_MINF_MAX_REACHED;
	  if (ret != NLOPT_SUCCESS) goto done;
	  if (++k > 1) memcpy(xprevprev, xprev, sizeof(double) * n);
	  memcpy(xprev, xcur, sizeof(double) * n);

	  while (1) { /* inner iterations */
	       double min_dual;
	       int feasible_cur, inner_done;
	       nlopt_result reti;

	       /* solve dual problem */
	       dd.rho = rho; dd.count = 0;
	       reti = nlopt_minimize(
		    dual_alg, m, dual_func, &dd,
		    dual_lb, dual_ub, y, &min_dual,
		    -HUGE_VAL, dual_tolrel,0., 0.,NULL, dual_maxeval,
		    stop->maxtime - (nlopt_seconds() - stop->start));
	       if (reti < 0 || reti == NLOPT_MAXTIME_REACHED) {
		    ret = reti;
		    goto done;
	       }

	       dual_func(m, y, NULL, &dd); /* evaluate final xcur etc. */
	       if (mma_verbose) {
		    printf("MMA dual converged in %d iterations to g=%g:\n",
			   dd.count, dd.gval);
		    for (i = 0; i < MIN(mma_verbose, m); ++i)
			 printf("    MMA y[%d]=%g, gc[%d]=%g\n",
				i, y[i], i, dd.gcval[i]);
	       }

	       fcur = f(n, xcur, dfdx_cur, f_data);
	       stop->nevals++;
	       feasible_cur = 1;
	       inner_done = dd.gval >= fcur;
	       for (i = 0; i < m; ++i) {
		    fcval_cur[i] = fc(n, xcur, dfcdx_cur + i*n, 
				      fc_data + fc_datum_size * i);
		    feasible_cur = feasible_cur && (fcval[i] <= 0);
		    inner_done = inner_done && (dd.gcval[i] >= fcval_cur[i]);
	       }

	       /* once we have reached a feasible solution, the
		  algorithm should never make the solution infeasible
		  again, although the constraints may be violated
		  slightly by rounding errors etc. so we must be a
		  little careful about checking feasibility */
	       if (feasible_cur) feasible = 1;

	       if (fcur < *minf) {
		    dd.fval = *minf = fcur;
		    memcpy(fcval, fcval_cur, sizeof(double)*m);
		    memcpy(x, xcur, sizeof(double)*n);
		    memcpy(dfdx, dfdx_cur, sizeof(double)*n);
		    memcpy(dfcdx, dfcdx_cur, sizeof(double)*n*m);
	       }
	       if (nlopt_stop_evals(stop)) ret = NLOPT_MAXEVAL_REACHED;
	       else if (nlopt_stop_time(stop)) ret = NLOPT_MAXTIME_REACHED;
	       else if (*minf < stop->minf_max) ret = NLOPT_MINF_MAX_REACHED;
	       if (ret != NLOPT_SUCCESS) goto done;

	       if (inner_done) break;

	       if (fcur > dd.gval)
		    rho = MIN(10*rho, 1.1 * (rho + (fcur-dd.gval) / dd.wval));
	       for (i = 0; i < m; ++i)
		    if (fcval_cur[i] > dd.gcval[i])
			 rhoc[i] = 
			      MIN(10*rhoc[i], 
				  1.1 * (rhoc[i] + (fcval_cur[i]-dd.gcval[i]) 
					 / dd.wval));
	       
	       if (mma_verbose)
		    printf("MMA inner iteration: rho -> %g\n", rho);
	       for (i = 0; i < MIN(mma_verbose, m); ++i)
		    printf("                 MMA rhoc[%d] -> %g\n", i,rhoc[i]);
	  }

	  if (nlopt_stop_ftol(stop, fcur, fprev))
	       ret = NLOPT_FTOL_REACHED;
	  if (nlopt_stop_x(stop, xcur, xprev))
	       ret = NLOPT_XTOL_REACHED;
	  if (ret != NLOPT_SUCCESS) goto done;
	       
	  /* update rho and sigma for iteration k+1 */
	  rho = MAX(0.1 * rho, MMA_RHOMIN);
	  if (mma_verbose)
	       printf("MMA outer iteration: rho -> %g\n", rho);
	  for (i = 0; i < m; ++i)
	       rhoc[i] = MAX(0.1 * rhoc[i], MMA_RHOMIN);
	  for (i = 0; i < MIN(mma_verbose, m); ++i)
	       printf("                 MMA rhoc[%d] -> %g\n", i, rhoc[i]);
	  if (k > 1) {
	       for (j = 0; j < n; ++j) {
		    double dx2 = (xcur[j]-xprev[j]) * (xprev[j]-xprevprev[j]);
		    double gam = dx2 < 0 ? 0.7 : (dx2 > 0 ? 1.2 : 1);
		    sigma[j] *= gam;
		    if (!nlopt_isinf(ub[j]) && !nlopt_isinf(lb[j])) {
			 sigma[j] = MIN(sigma[j], 10*(ub[j]-lb[j]));
			 sigma[j] = MAX(sigma[j], 0.01*(ub[j]-lb[j]));
		    }
	       }
	       for (j = 0; j < MIN(mma_verbose, n); ++j)
		    printf("                 MMA sigma[%d] -> %g\n", 
			   j, sigma[j]);
	  }
     }

 done:
     free(sigma);
     return ret;
}
