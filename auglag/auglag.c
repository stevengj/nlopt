#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "auglag.h"

int auglag_verbose = 1;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/***************************************************************************/

typedef struct {
     nlopt_func f; void *f_data;
     int m; nlopt_constraint *fc;
     int p; nlopt_constraint *h;
     double rho, *lambda, *mu;
     double *gradtmp;
     nlopt_stopping *stop;
} auglag_data;

/* the augmented lagrangian objective function */
static double auglag(int n, const double *x, double *grad, void *data)
{
     auglag_data *d = (auglag_data *) data;
     double *gradtmp = grad ? d->gradtmp : NULL;
     double rho = d->rho;
     const double *lambda = d->lambda, *mu = d->mu;
     double L;
     int i, j;

     L = d->f(n, x, grad, d->f_data);

     for (i = 0; i < d->p; ++i) {
	  double h;
	  h = d->h[i].f(n, x, gradtmp, d->h[i].f_data) + lambda[i] / rho;
	  L += 0.5 * rho * h*h;
	  if (grad) for (j = 0; j < n; ++j) grad[j] += (rho * h) * gradtmp[j];
     }

     for (i = 0; i < d->m; ++i) {
	  double fc;
	  fc = d->fc[i].f(n, x, gradtmp, d->fc[i].f_data) + mu[i] / rho;
	  if (fc > 0) {
	       L += 0.5 * rho * fc*fc;
	       if (grad) for (j = 0; j < n; ++j) 
		    grad[j] += (rho * fc) * gradtmp[j];
	  }
     }
     
     d->stop->nevals++;

     return L;
}

/***************************************************************************/

nlopt_result auglag_minimize(int n, nlopt_func f, void *f_data,
			     int m, nlopt_constraint *fc,
			     int p, nlopt_constraint *h,
			     const double *lb, const double *ub, /* bounds */
			     double *x, /* in: initial guess, out: minimizer */
			     double *minf,
			     nlopt_stopping *stop,
			     nlopt_opt sub_opt, int sub_has_fc)
{
     auglag_data d;
     nlopt_result ret = NLOPT_SUCCESS;
     double ICM = HUGE_VAL;
     double *xcur = NULL, fcur;
     int i, feasible;

     /* magic parameters from Birgin & Martinez */
     const double tau = 0.5, gam = 10;
     const double lam_min = -1e20, lam_max = 1e20, mu_max = 1e20;

     d.f = f; d.f_data = f_data;
     d.m = m; d.fc = fc;
     d.p = p; d.h = h;
     d.stop = stop;

     /* whether we handle inequality constraints via the augmented
	Lagrangian penalty function, or directly in the sub-algorithm */
     if (sub_has_fc)
	  d.m = 0;
     else
	  m = 0;

     ret = nlopt_set_min_objective(sub_opt, auglag, &d); if (ret<0) return ret;
     ret = nlopt_set_lower_bounds(sub_opt, lb); if (ret<0) return ret;
     ret = nlopt_set_upper_bounds(sub_opt, ub); if (ret<0) return ret;
     ret = nlopt_remove_inequality_constraints(sub_opt); if (ret<0) return ret;
     ret = nlopt_remove_equality_constraints(sub_opt); if (ret<0) return ret;
     for (i = 0; i < m; ++i) {
	  ret = nlopt_add_inequality_constraint(sub_opt, fc[i].f, fc[i].f_data,
						fc[i].tol);
	  if (ret < 0) return ret;
     }

     xcur = (double *) malloc(sizeof(double) * (2*n + d.p + d.m));
     if (!xcur) return NLOPT_OUT_OF_MEMORY;
     memcpy(xcur, x, sizeof(double) * n);

     d.gradtmp = xcur + n;
     memset(d.gradtmp, 0, sizeof(double) * (n + d.p + d.m));
     d.lambda = d.gradtmp + n;
     d.mu = d.lambda + d.p;

     /* starting rho suggested by B & M */
     {
	  double con2 = 0;
	  d.stop->nevals++;
	  fcur = f(n, xcur, NULL, f_data);
	  feasible = 1;
	  for (i = 0; i < d.p; ++i) {
               double hi = h[i].f(n, xcur, NULL, d.h[i].f_data);
	       feasible = feasible && fabs(hi) <= h[i].tol;
	       con2 += hi * hi;
	  }
	  for (i = 0; i < d.m; ++i) {
               double fci = fc[i].f(n, xcur, NULL, d.fc[i].f_data);
	       feasible = feasible && fci <= fc[i].tol;
	       if (fci > 0) con2 += fci * fci;
	  }
	  d.rho = MAX(1e-6, MIN(10, 2 * fabs(*minf) / con2));
     }

     if (feasible)
	  *minf = fcur;
     else
	  *minf = HUGE_VAL;

     do {
	  double prev_ICM = ICM;
	  
	  ret = nlopt_optimize_limited(sub_opt, xcur, minf,
				       stop->maxeval - stop->nevals,
				       stop->maxtime - (nlopt_seconds() 
							- stop->start));
	  if (ret < 0) break;
	  
	  d.stop->nevals++;
	  fcur = f(n, xcur, NULL, f_data);
	  
	  ICM = 0;
	  feasible = 1;
	  for (i = 0; i < d.p; ++i) {
	       double hi = h[i].f(n, xcur, NULL, d.h[i].f_data);
	       double newlam = d.lambda[i] + d.rho * hi;
	       feasible = feasible && fabs(hi) <= h[i].tol;
	       ICM = MAX(ICM, fabs(hi));
	       d.lambda[i] = MIN(MAX(lam_min, newlam), lam_max);
	  }
	  for (i = 0; i < d.m; ++i) {
	       double fci = fc[i].f(n, xcur, NULL, d.fc[i].f_data);
	       double newmu = d.mu[i] + d.rho * fci;
	       feasible = feasible && fci <= fc[i].tol;
	       ICM = MAX(ICM, fabs(MAX(fci, -d.mu[i] / d.rho)));
	       d.mu[i] = MIN(MAX(0.0, newmu), mu_max);
	  }
	  if (ICM > tau * prev_ICM)
	       d.rho *= gam;

	  if (auglag_verbose) {
	       printf("auglag: ICM=%g, rho=%g\nauglag lambda=", ICM, d.rho);
	       for (i = 0; i < d.p; ++i) printf(" %g", d.lambda[i]);
	       printf("\nauglag mu = ");
	       for (i = 0; i < d.m; ++i) printf(" %g", d.mu[i]);
	       printf("\n");
	  }

	  /* only check f & x convergence for feasible points...
	     for this to be effective on active constraints, the user
	     must set some nonzero tolerance for each constraint */
	  if (feasible && fcur < *minf) {
	       ret = NLOPT_SUCCESS;
	       if (fcur < stop->minf_max) ret = NLOPT_MINF_MAX_REACHED;
	       if (nlopt_stop_ftol(stop, fcur,*minf)) ret = NLOPT_FTOL_REACHED;
	       if (nlopt_stop_x(stop, xcur, x)) ret = NLOPT_XTOL_REACHED;
	       *minf = fcur;
	       memcpy(x, xcur, sizeof(double) * n);
	       if (ret != NLOPT_SUCCESS) break;
	  }

	  if (nlopt_stop_evals(stop)) {ret = NLOPT_MAXEVAL_REACHED; break;}
          if (nlopt_stop_time(stop)) {ret = NLOPT_MAXTIME_REACHED; break;}

	  /* TODO: use some stopping criterion on ICM? */
     } while (1);

     free(xcur);
     return ret;
}
