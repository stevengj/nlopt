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
     int m; nlopt_func fc; char *fc_data; ptrdiff_t fc_datum_size;
     int p; nlopt_func h; char *h_data; ptrdiff_t h_datum_size;
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
	  h = d->h(n, x, gradtmp, d->h_data + d->h_datum_size * i)
	       + lambda[i] / rho;
	  L += 0.5 * rho * h*h;
	  if (grad) for (j = 0; j < n; ++j) grad[j] += (rho * h) * gradtmp[j];
     }

     for (i = 0; i < d->m; ++i) {
	  double fc;
	  fc = d->fc(n, x, gradtmp, d->fc_data + d->fc_datum_size * i)
	       + mu[i] / rho;
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
			     int m, nlopt_func fc,
			     void *fc_data, ptrdiff_t fc_datum_size,
			     int p, nlopt_func h,
			     void *h_data, ptrdiff_t h_datum_size,
			     const double *lb, const double *ub, /* bounds */
			     double *x, /* in: initial guess, out: minimizer */
			     double *minf,
			     nlopt_stopping *stop,
			     nlopt_algorithm sub_alg, int sub_has_fc)
{
     auglag_data d;
     nlopt_result ret = NLOPT_SUCCESS;
     double ICM = HUGE_VAL;
     int i;

     /* magic parameters from Birgin & Martinez */
     const double tau = 0.5, gam = 10;
     const double lam_min = -1e20, lam_max = 1e20, mu_max = 1e20;

     d.f = f; d.f_data = f_data;
     d.m = m; d.fc = fc; d.fc_data = (char *) fc_data;
     d.fc_datum_size = fc_datum_size;
     d.p = p; d.h = h; d.h_data = (char *) h_data;
     d.h_datum_size = h_datum_size;
     d.stop = stop;

     /* whether we handle inequality constraints via the augmented
	Lagrangian penalty function, or directly in the sub-algorithm */
     if (sub_has_fc)
	  d.m = 0;
     else
	  m = 0;

     d.gradtmp = (double *) malloc(sizeof(double) * (n + d.p + d.m));
     if (!d.gradtmp) return NLOPT_OUT_OF_MEMORY;
     memset(d.gradtmp, 0, sizeof(double) * (n + d.p + d.m));
     d.lambda = d.gradtmp + n;
     d.mu = d.lambda + d.p;

     /* starting rho suggested by B & M */
     {
	  double con2 = 0;
	  *minf = f(n, x, NULL, f_data);
	  for (i = 0; i < d.p; ++i) {
               double hi = h(n, x, NULL, d.h_data + i*h_datum_size);
	       con2 += hi * hi;
	  }
	  for (i = 0; i < d.m; ++i) {
               double fci = fc(n, x, NULL, d.fc_data + i*fc_datum_size);
	       if (fci > 0) con2 += fci * fci;
	  }
	  d.rho = MAX(1e-6, MIN(10, 2 * fabs(*minf) / con2));
     }

     *minf = HUGE_VAL;

     do {
	  double prev_ICM = ICM;

	  ret = nlopt_minimize_constrained(sub_alg, n, auglag, &d,
					   m, fc, fc_data, fc_datum_size,
					   lb, ub, x, minf,
					   -HUGE_VAL, 
					   stop->ftol_rel, stop->ftol_abs,
					   stop->xtol_rel, stop->xtol_abs,
					   stop->maxeval - stop->nevals,
					   stop->maxtime - (nlopt_seconds() 
							    - stop->start));
	  if (ret < 0) break;
	  if (nlopt_stop_evals(stop)) {ret = NLOPT_MAXEVAL_REACHED; break;}
          if (nlopt_stop_time(stop)) {ret = NLOPT_MAXTIME_REACHED; break;}
	  
	  *minf = f(n, x, NULL, f_data);
	  
	  ICM = 0;
	  for (i = 0; i < d.p; ++i) {
	       double hi = h(n, x, NULL, d.h_data + i*h_datum_size);
	       double newlam = d.lambda[i] + d.rho * hi;
	       ICM = MAX(ICM, fabs(hi));
	       d.lambda[i] = MIN(MAX(lam_min, newlam), lam_max);
	  }
	  for (i = 0; i < d.m; ++i) {
	       double fci = fc(n, x, NULL, d.fc_data + i*fc_datum_size);
	       double newmu = d.mu[i] + d.rho * fci;
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

	  if (ICM <= stop->htol_abs || ICM <= stop->htol_rel * fabs(*minf)) {
	       ret = NLOPT_FTOL_REACHED;
	       break;
	  }
     } while (1);

     free(d.gradtmp);
     return ret;
}
