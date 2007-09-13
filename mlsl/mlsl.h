#ifndef MLSL_H
#define MLSL_H

#include "nlopt.h"
#include "nlopt-util.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

nlopt_result mlsl_minimize(int n, nlopt_func f, void *f_data,
			   const double *lb, const double *ub, /* bounds */
			   double *x, /* in: initial guess, out: minimizer */
			   double *minf,
			   nlopt_stopping *stop,
                           nlopt_algorithm local_alg,
                           int local_maxeval,
                           int lds);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif

