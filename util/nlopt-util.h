#ifndef NLOPT_UTIL_H
#define NLOPT_UTIL_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/* seconds timer */
extern double nlopt_seconds(void);
extern unsigned long nlopt_time_seed(void);

/* pseudorandom number generation by Mersenne twister algorithm */
extern void nlopt_init_genrand(unsigned long s);
extern double nlopt_urand(double a, double b);
extern int nlopt_iurand(int n);

/* Sobol' low-discrepancy-sequence generation */
typedef struct nlopt_soboldata_s *nlopt_sobol;
extern nlopt_sobol nlopt_sobol_create(unsigned sdim);
extern void nlopt_sobol_destroy(nlopt_sobol s);
extern void nlopt_sobol_next01(nlopt_sobol s, double *x);
extern void nlopt_sobol_next(nlopt_sobol s, double *x,
			    const double *lb, const double *ub);

/* stopping criteria */
typedef struct {
     int n;
     double minf_max;
     double ftol_rel;
     double ftol_abs;
     double xtol_rel;
     const double *xtol_abs;
     int nevals, maxeval;
     double maxtime, start;
} nlopt_stopping;
extern int nlopt_stop_f(const nlopt_stopping *stop, double f, double oldf);
extern int nlopt_stop_ftol(const nlopt_stopping *stop, double f, double oldf);
extern int nlopt_stop_x(const nlopt_stopping *stop, 
			const double *x, const double *oldx);
extern int nlopt_stop_xs(const nlopt_stopping *stop, 
			 const double *xs, const double *oldxs,
			 const double *scale_min, const double *scale_max);
extern int nlopt_stop_evals(const nlopt_stopping *stop);
extern int nlopt_stop_time(const nlopt_stopping *stop);
extern int nlopt_stop_evalstime(const nlopt_stopping *stop);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
