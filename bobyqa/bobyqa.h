#ifndef BOBYQA_H
#define BOBYQA_H 1

#include "nlopt-util.h"
#include "nlopt.h"

typedef double (*bobyqa_func)(int n, const double *x, void *func_data);

extern nlopt_result bobyqa(int n, int npt, double *x, 
			   const double *lb, const double *ub,
			   double rhobeg, nlopt_stopping *stop, double *minf,
			   bobyqa_func calfun, void *calfun_data);

#endif /* BOBYQA_H */
