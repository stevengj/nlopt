#ifndef CRS_H
#define CRS_H

#include "nlopt.h"
#include "nlopt-util.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

nlopt_result crs_minimize(int n, nlopt_func f, void *f_data,
			  const double *lb, const double *ub, /* bounds */
			  double *x, /* in: initial guess, out: minimizer */
			  double *minf,
			  nlopt_stopping *stop,
			  int random); /* random or low-discrepancy seq. */

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif

