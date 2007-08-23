#ifndef NLOPT_H
#define NLOPT_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

typedef double (*nlopt_func)(int n, const double *x,
			     double *gradient, /* NULL if not needed */
			     void *func_data);

typedef enum {
     /* non-gradient methods */
     NLOPT_GLOBAL_DIRECT,
     NLOPT_LOCAL_SUBPLEX

     /* gradient-based methods */
     NLOPT_GLOBAL_STOGO,
     NLOPT_LOCAL_LBFGS,
} nlopt_method;

typedef enum {
     NLOPT_FAILURE = -1, /* generic failure code */
     NLOPT_INVALID_ARGS = -2,
     NLOPT_OUT_OF_MEMORY = -3,

     NLOPT_SUCCESS = 1, /* generic success code */
     NLOPT_FMIN_MAX_REACHED = 2,
     NLOPT_FTOL_REACHED = 3,
     NLOPT_XTOL_REACHED = 4,
     NLOPT_MAXEVAL_REACHED = 5,
     NLOPT_MAXTIME_REACHED = 6
} nlopt_result;

extern nlopt_result nlopt_minimize(
     nlopt_method method,
     int n, nlopt_func f, void *f_data,
     const double *lb, const double *ub, /* bounds */
     double *x, /* in: initial guess, out: minimizer */
     double *fmin, /* out: minimum */
     double fmin_max, ftol_rel, double ftol_abs,
     double xtol_rel, const double *xtol_abs,
     int maxeval, double maxtime);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
