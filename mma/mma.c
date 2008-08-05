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
   should depend on the objective function in some way */
#define MMA_RHOMIN 1e-5

nlopt_result mma_minimize(int n, nlopt_func f, void *f_data,
			  const double *lb, const double *ub, /* bounds */
			  double *x, /* in: initial guess, out: minimizer */
			  double *minf,
			  nlopt_stopping *stop)
{
     nlopt_result ret = NLOPT_SUCCESS;
     double *xcur, rho, *sigma, *dfdx, *dfdx_cur, *xprev, *xprevprev, fcur;
     int j, k = 0;
     
     sigma = (double *) malloc(sizeof(double) * 6*n);
     if (!sigma) return NLOPT_OUT_OF_MEMORY;
     dfdx = sigma + n;
     dfdx_cur = dfdx + n;
     xcur = dfdx_cur + n;
     xprev = xcur + n;
     xprevprev = xprev + n;
     for (j = 0; j < n; ++j)
	  sigma[j] = 0.5 * (ub[j] - lb[j]);
     rho = 1.0;

     fcur = *minf = f(n, x, dfdx, f_data);
     memcpy(xcur, x, sizeof(double) * n);
     stop->nevals++;
     while (1) { /* outer iterations */
	  double fprev = fcur;
	  if (nlopt_stop_evals(stop)) ret = NLOPT_MAXEVAL_REACHED;
	  else if (nlopt_stop_time(stop)) ret = NLOPT_MAXTIME_REACHED;
	  else if (*minf < stop->minf_max) ret = NLOPT_MINF_MAX_REACHED;
	  if (ret != NLOPT_SUCCESS) goto done;
	  if (++k > 1) memcpy(xprevprev, xprev, sizeof(double) * n);
	  memcpy(xprev, xcur, sizeof(double) * n);

	  while (1) { /* inner iterations */
	       double gcur = *minf, w = 0;
	       for (j = 0; j < n; ++j) {
		    /* because we have no constraint functions, minimizing
		       the MMA approximate function can be done analytically */
		    double dx = -sigma[j]*sigma[j]*dfdx[j]
			 / (2*sigma[j]*fabs(dfdx[j]) + rho);
		    xcur[j] = x[j] + dx;
		    if (xcur[j] > x[j] + 0.9*sigma[j])
			 xcur[j] = x[j] + 0.9*sigma[j];
		    else if (xcur[j] < x[j] - 0.9*sigma[j])
			 xcur[j] = x[j] - 0.9*sigma[j];
		    if (xcur[j] > ub[j]) xcur[j] = ub[j];
		    else if (xcur[j] < lb[j]) xcur[j] = lb[j];
		    dx = xcur[j] - x[j];
		    gcur += (sigma[j]*sigma[j]*dfdx[j]*dx
			     + sigma[j]*fabs(dfdx[j])*dx*dx
			     + 0.5*rho*dx*dx) / (sigma[j]*sigma[j]-dx*dx);
		    w += 0.5*rho*dx*dx / (sigma[j]*sigma[j]-dx*dx);
	       }

	       fcur = f(n, xcur, dfdx_cur, f_data);
	       stop->nevals++;
	       if (fcur < *minf) {
		    *minf = fcur;
		    memcpy(x, xcur, sizeof(double)*n);
		    memcpy(dfdx, dfdx_cur, sizeof(double)*n);
	       }
	       if (nlopt_stop_evals(stop)) ret = NLOPT_MAXEVAL_REACHED;
	       else if (nlopt_stop_time(stop)) ret = NLOPT_MAXTIME_REACHED;
	       else if (*minf < stop->minf_max) ret = NLOPT_MINF_MAX_REACHED;
	       if (ret != NLOPT_SUCCESS) goto done;

	       if (gcur >= fcur) break;
	       rho = MIN(10*rho, 1.1 * (rho + (fcur - gcur) / w));
	       if (mma_verbose)
		    printf("MMA inner iteration: rho -> %g\n", rho);
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
	  if (k > 1)
	       for (j = 0; j < n; ++j) {
		    double dx2 = (xcur[j]-xprev[j]) * (xprev[j]-xprevprev[j]);
		    double gam = dx2 < 0 ? 0.7 : (dx2 > 0 ? 1.2 : 1);
		    sigma[j] *= gam;
		    sigma[j] = MIN(sigma[j], 10*(ub[j]-lb[j]));
		    sigma[j] = MAX(sigma[j], 0.01*(ub[j]-lb[j]));
	       }
     }

 done:
     free(sigma);
     return ret;
}
