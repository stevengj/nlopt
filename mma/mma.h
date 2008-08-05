#ifndef MMA_H
#define MMA_H

#include "nlopt.h"
#include "nlopt-util.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

extern int mma_verbose;

nlopt_result mma_minimize(int n, nlopt_func f, void *f_data,
			  const double *lb, const double *ub, /* bounds */
			  double *x, /* in: initial guess, out: minimizer */
			  double *minf,
			  nlopt_stopping *stop);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif

