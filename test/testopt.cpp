#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "config.h"

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif
#ifdef HAVE_GETOPT_H
#  include <getopt.h>
#endif
#if TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# if HAVE_SYS_TIME_H
#  include <sys/time.h>
# else
#  include <time.h>
# endif
#endif

#include "nlopt.h"
#include "testfuncs.h"

static double urand(double a, double b)
{
  return a + (rand() * (b - a) / RAND_MAX);
}

static nlopt_algorithm algorithm = NLOPT_GLOBAL_DIRECT;
static double ftol_rel = 0, ftol_abs = 0, xtol_rel = 0;
static int maxeval = 1000;
static double maxtime = 0.0;

static int test_function(int ifunc)
{
  testfunc func;
  int i;
  double *x, fmin, f0;
  nlopt_result ret;
  
  if (ifunc < 0 || ifunc >= NTESTFUNCS) {
    fprintf(stderr, "testopt: invalid function %d\n", ifunc);
    return 0;
  }
  func = testfuncs[ifunc];
  if (!func.has_gradient && algorithm >= NLOPT_GLOBAL_STOGO) {
    fprintf(stderr, 
	    "testopt: A function with gradients is required for %s\n",
	    nlopt_algorithm_name(algorithm));
    return 0;
  }
  x = (double *) malloc(sizeof(double) * func.n * 2);
  if (!x) { fprintf(stderr, "testopt: Out of memory!\n"); return 0; }

  
  printf("-----------------------------------------------------------\n");
  printf("Optimizing %s (%d dims) using %s algorithm\n",
	 func.name, func.n, nlopt_algorithm_name(algorithm));
  printf("Starting guess x = [");
  for (i = 0; i < func.n; ++i)
    printf(" %g", x[i] = urand(func.lb[i], func.ub[i]));
  printf("]\n");
  f0 = func.f(func.n, x, x + func.n, func.f_data);
  printf("Starting function value = %g\n", f0);

  if (testfuncs_verbose && func.has_gradient) {
    printf("checking gradient:\n");
    for (i = 0; i < func.n; ++i) {
      double f;
      x[i] *= 1 + 1e-6;
      f = func.f(func.n, x, NULL, func.f_data);
      x[i] /= 1 + 1e-6;
      printf("  grad[%d] = %g vs. numerical derivative %g\n",
	     i, x[i + func.n], (f - f0) / (x[i] * 1e-6));
    }
  }

  testfuncs_counter = 0;
  ret = nlopt_minimize(algorithm,
		       func.n, func.f, func.f_data,
		       func.lb, func.ub,
		       x, &fmin,
		       HUGE_VAL, ftol_rel, ftol_abs, xtol_rel, NULL,
		       maxeval, maxtime);
  printf("return code %d from nlopt_minimize\n", ret);
  if (ret < 0) {
    fprintf(stderr, "testopt: error in nlopt_minimize\n");
    return 0;
  }
  printf("Found minimum f = %g after %d evaluations.\n", 
	 fmin, testfuncs_counter);
  printf("Minimum at x = [");
  for (i = 0; i < func.n; ++i) printf(" %g", x[i]);
  printf("]\n");
  printf("vs. global minimum f = %g at x = [", func.fmin);
  for (i = 0; i < func.n; ++i) printf(" %g", func.xmin[i]);
  printf("]\n");
  printf("|f - fmin| = %g, |f - fmin| / |fmin| = %e\n",
	 fabs(fmin - func.fmin), fabs(fmin - func.fmin) / fabs(func.fmin));
  
  free(x);
  return 1;
}

static void usage(FILE *f)
{
  fprintf(f, "Usage: testopt [OPTIONS]\n"
	  "Options:\n"
	  "     -h : print this help\n"
	  "     -v : verbose mode\n"
	  " -r <s> : use random seed <s> for starting guesses\n"
	  " -a <n> : use optimization algorithm <n>\n"
	  " -o <n> : use objective function <n>\n"
	  " -e <n> : use at most <n> evals (default: %d)\n",
	  maxeval);
}

int main(int argc, char **argv)
{
  int c;
  
  srand((unsigned) time(NULL));
  testfuncs_verbose = 0;
  
  if (argc <= 1)
    usage(stdout);
  
  while ((c = getopt(argc, argv, "hvra:o:e:")) != -1)
    switch (c) {
    case 'h':
      usage(stdout);
      return EXIT_SUCCESS;
    case 'v':
      testfuncs_verbose = 1;
      break;
    case 'r':
      srand((unsigned) atoi(optarg));
      break;
    case 'a':
      c = atoi(optarg);
      if (c < 0 || c >= NLOPT_NUM_ALGORITHMS) {
	fprintf(stderr, "testopt: invalid algorithm %d\n", c);
	return EXIT_FAILURE;
      }
      algorithm = (nlopt_algorithm) c;
      break;
    case 'o':
      if (!test_function(atoi(optarg)))
	return EXIT_FAILURE;
      break;
    case 'e':
      maxeval = atoi(optarg);
      break;
    default:
      fprintf(stderr, "harminv: invalid argument -%c\n", c);
      usage(stderr);
      return EXIT_FAILURE;
    }
  
  return EXIT_SUCCESS;
}
