#ifndef PRAXIS_H
#define PRAXIS_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

typedef double (*praxis_func)(int n, const double *x, void *f_data);

double praxis_(double *t0, double *machep, double *h0,
               int *n, double *x, praxis_func f, void *f_data);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* PRAXIS_H */
