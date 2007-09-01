#ifndef CDIRECT_H
#define CDIRECT_H

#include "nlopt-util.h"
#include "nlopt.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

extern nlopt_result cdirect_unscaled(int n, nlopt_func f, void *f_data,
				     const double *lb, const double *ub,
				     double *x,
				     double *fmin,
				     nlopt_stopping *stop,
				     double magic_eps, int which_alg);

extern nlopt_result cdirect(int n, nlopt_func f, void *f_data,
			    const double *lb, const double *ub,
			    double *x,
			    double *fmin,
			    nlopt_stopping *stop,
			    double magic_eps, int which_alg);

extern nlopt_result cdirect_hybrid(int n, nlopt_func f, void *f_data,
			    const double *lb, const double *ub,
			    double *x,
			    double *fmin,
			    nlopt_stopping *stop,
			    nlopt_algorithm local_alg,
			    int local_maxeval,
			    int randomized_div);

extern nlopt_result cdirect_hybrid_unscaled(int n, nlopt_func f, void *f_data,
			    const double *lb, const double *ub,
			    double *x,
			    double *fmin,
			    nlopt_stopping *stop,
			    nlopt_algorithm local_alg,
			    int local_maxeval,
			    int randomized_div);

/* internal routines and data structures: */
extern int cdirect_hyperrect_compare(double *a, double *b);
typedef struct {
     nlopt_func f;
     void *f_data;
     double *x;
     const double *lb, *ub;
} cdirect_uf_data;
extern double cdirect_uf(int n, const double *xu, double *grad, void *d_);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* DIRECT_H */
