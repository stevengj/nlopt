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
     /* Naming conventions:

        NLOPT_{G/L}{D/N}_* 
	    = global/local derivative/no-derivative optimization, 
              respectively 
 
	*_RAND algorithms involve some randomization.

	*_NOSCAL algorithms are *not* scaled to a unit hypercube
	         (i.e. they are sensitive to the units of x)
	*/

     NLOPT_GN_DIRECT = 0,
     NLOPT_GN_DIRECT_L,
     NLOPT_GN_DIRECT_L_RAND,
     NLOPT_GN_DIRECT_NOSCAL,
     NLOPT_GN_DIRECT_L_NOSCAL,
     NLOPT_GN_DIRECT_L_RAND_NOSCAL,

     NLOPT_GN_ORIG_DIRECT,
     NLOPT_GN_ORIG_DIRECT_L,

     NLOPT_LN_SUBPLEX,

     NLOPT_GD_STOGO,
     NLOPT_GD_STOGO_RAND,

     NLOPT_LD_LBFGS,

     NLOPT_LN_PRAXIS,

     NLOPT_LD_VAR1,
     NLOPT_LD_VAR2,

     NLOPT_LD_TNEWTON,
     NLOPT_LD_TNEWTON_RESTART,
     NLOPT_LD_TNEWTON_PRECOND,
     NLOPT_LD_TNEWTON_PRECOND_RESTART,

     NLOPT_GN_CRS2_LM,

     NLOPT_NUM_ALGORITHMS /* not an algorithm, just the number of them */
} nlopt_algorithm;

extern const char *nlopt_algorithm_name(nlopt_algorithm a);

typedef enum {
     NLOPT_FAILURE = -1, /* generic failure code */
     NLOPT_INVALID_ARGS = -2,
     NLOPT_OUT_OF_MEMORY = -3,

     NLOPT_SUCCESS = 1, /* generic success code */
     NLOPT_MINF_MAX_REACHED = 2,
     NLOPT_FTOL_REACHED = 3,
     NLOPT_XTOL_REACHED = 4,
     NLOPT_MAXEVAL_REACHED = 5,
     NLOPT_MAXTIME_REACHED = 6
} nlopt_result;

extern nlopt_result nlopt_minimize(
     nlopt_algorithm algorithm,
     int n, nlopt_func f, void *f_data,
     const double *lb, const double *ub, /* bounds */
     double *x, /* in: initial guess, out: minimizer */
     double *minf, /* out: minimum */
     double minf_max, double ftol_rel, double ftol_abs,
     double xtol_rel, const double *xtol_abs,
     int maxeval, double maxtime);

extern void nlopt_srand(unsigned long seed);
extern void nlopt_srand_time(void);

extern void nlopt_version(int *major, int *minor, int *bugfix);

extern void nlopt_get_local_search_algorithm(nlopt_algorithm *deriv,
					     nlopt_algorithm *nonderiv,
					     int *maxeval);
extern void nlopt_set_local_search_algorithm(nlopt_algorithm deriv,
					     nlopt_algorithm nonderiv,
					     int maxeval);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
