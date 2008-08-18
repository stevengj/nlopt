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
			  int m, nlopt_func fc,
			  void *fc_data_, ptrdiff_t fc_datum_size,
			  const double *lb, const double *ub, /* bounds */
			  double *x, /* in: initial guess, out: minimizer */
			  double *minf,
			  nlopt_stopping *stop,
			  nlopt_algorithm dual_alg, 
			  double dual_tolrel, int dual_maxeval);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif

