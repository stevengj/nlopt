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

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
