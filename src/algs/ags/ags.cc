// A C-callable front-end to the AGS global-optimization library.
//  -- Vladislav Sovrasov

#include "ags.h"
#include "solver.hpp"

/*
int ags_minimize(int n,
		   objective_func fgrad, void *data,
		   double *x, double *minf,
		   const double *l, const double *u,
#ifdef NLOPT_UTIL_H
		   nlopt_stopping *stop,
#else
		   long int maxeval, double maxtime,
#endif
		   int nrandom)
       */
int ags_minimize()
{
  return 1;
}
