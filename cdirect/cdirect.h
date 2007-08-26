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

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* DIRECT_H */
