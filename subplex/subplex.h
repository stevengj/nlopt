#ifndef SUBPLEX_H
#define SUBPLEX_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

typedef double (*subplex_func)(int n, const double *x, void *func_data);

extern int subplex(subplex_func f, double *fmin, double *x, int n, void *fdata,
		   double tol, int maxnfe,
		   double fmin_max, int use_fmin_max,
		   const double *scale);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
