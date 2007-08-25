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

#include "nlopt.h"
#include "nlopt-util.h"
#include "testfuncs.h"

static nlopt_algorithm algorithm = NLOPT_GLOBAL_DIRECT;
static double ftol_rel = 0, ftol_abs = 0, xtol_rel = 0, xtol_abs = 0, fmin_max = -HUGE_VAL;
static int maxeval = 1000, iterations = 1;
static double maxtime = 0.0;
static double xinit_tol = -1;

static void listalgs(FILE *f)
{
  int i;
  fprintf(f, "Available algorithms:\n");
  for (i = 0; i < NLOPT_NUM_ALGORITHMS; ++i)
    fprintf(f, "  %2d: %s\n", i, nlopt_algorithm_name((nlopt_algorithm) i));
}

static void listfuncs(FILE *f)
{
  int i;
  fprintf(f, "Available objective functions:\n");
  for (i = 0; i < NTESTFUNCS; ++i)
    fprintf(f, "  %2d: %s (%d dims)\n", i, testfuncs[i].name, testfuncs[i].n);
}

static int test_function(int ifunc)
{
  testfunc func;
  int i, iter;
  double *x, fmin, f0, *xtabs;
  nlopt_result ret;
  double start = nlopt_seconds();
  
  if (ifunc < 0 || ifunc >= NTESTFUNCS) {
    fprintf(stderr, "testopt: invalid function %d\n", ifunc);
    listfuncs(stderr);
    return 0;
  }
  func = testfuncs[ifunc];
  if (!func.has_gradient && algorithm >= NLOPT_GLOBAL_STOGO) {
    fprintf(stderr, 
	    "testopt: A function with gradients is required for %s\n",
	    nlopt_algorithm_name(algorithm));
    return 0;
  }
  x = (double *) malloc(sizeof(double) * func.n * 3);
  if (!x) { fprintf(stderr, "testopt: Out of memory!\n"); return 0; }

  xtabs = x + func.n * 2;
  for (i = 0; i < func.n; ++i) xtabs[i] = xtol_abs;

  
  printf("-----------------------------------------------------------\n");
  printf("Optimizing %s (%d dims) using %s algorithm\n",
	 func.name, func.n, nlopt_algorithm_name(algorithm));
  printf("lower bounds at lb = [");
  for (i = 0; i < func.n; ++i) printf(" %g", func.lb[i]);
  printf("]\n");
  printf("upper bounds at ub = [");
  for (i = 0; i < func.n; ++i) printf(" %g", func.ub[i]);
  printf("]\n");


  printf("Starting guess x = [");
  for (i = 0; i < func.n; ++i) {
    if (xinit_tol < 0) { /* random starting point near center of box */
      double dx = (func.ub[i] - func.lb[i]) * 0.25;
      double xm = 0.5 * (func.ub[i] + func.lb[i]);
      x[i] = nlopt_urand(xm - dx, xm + dx);
    }
    else 
      x[i] = nlopt_urand(-xinit_tol, xinit_tol)
	+ (1 + nlopt_urand(-xinit_tol, xinit_tol)) * func.xmin[i];
    printf(" %g", x[i]);
  }
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

  for (iter = 0; iter < iterations; ++iter) {
    testfuncs_counter = 0;
    ret = nlopt_minimize(algorithm,
			 func.n, func.f, func.f_data,
			 func.lb, func.ub,
			 x, &fmin,
			 fmin_max, ftol_rel, ftol_abs, xtol_rel, xtabs,
			 maxeval, maxtime);
    printf("finished after %g seconds.\n", nlopt_seconds() - start);
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
    printf("|f - fmin| = %g, |f - fmin| / |fmin| = %e\n",
	   fabs(fmin - func.fmin), fabs(fmin - func.fmin) / fabs(func.fmin));
  }
  printf("vs. global minimum f = %g at x = [", func.fmin);
  for (i = 0; i < func.n; ++i) printf(" %g", func.xmin[i]);
  printf("]\n");
  
  free(x);
  return 1;
}

static void usage(FILE *f)
{
  fprintf(f, "Usage: testopt [OPTIONS]\n"
	  "Options:\n"
	  "     -h : print this help\n"
	  "     -L : list available algorithms and objective functions\n"
	  "     -v : verbose mode\n"
	  " -a <n> : use optimization algorithm <n>\n"
	  " -o <n> : use objective function <n>\n"
	  " -0 <x> : starting guess within <x> + (1+<x>) * optimum\n"
	  " -e <n> : use at most <n> evals (default: %d, 0 to disable)\n"
	  " -t <t> : use at most <t> seconds (default: disabled)\n"
	  " -x <t> : relative tolerance <t> on x (default: disabled)\n"
	  " -X <t> : absolute tolerance <t> on x (default: disabled)\n"
	  " -f <t> : relative tolerance <t> on f (default: disabled)\n"
	  " -F <t> : absolute tolerance <t> on f (default: disabled)\n"
	  " -m <m> : minimize f until <m> is reached (default: disabled)\n"
	  " -i <n> : iterate optimization <n> times (default: 1)\n"
	  " -r <s> : use random seed <s> for starting guesses\n"
	  , maxeval);
}

int main(int argc, char **argv)
{
  int c;
  
  nlopt_srand_time();
  testfuncs_verbose = 0;
  
  if (argc <= 1)
    usage(stdout);
  
  while ((c = getopt(argc, argv, "hLv0:r:a:o:i:e:t:x:X:f:F:m:")) != -1)
    switch (c) {
    case 'h':
      usage(stdout);
      return EXIT_SUCCESS;
    case 'L':
      listalgs(stdout);
      listfuncs(stdout);
      return EXIT_SUCCESS;
    case 'v':
      testfuncs_verbose = 1;
      break;
    case 'r':
      nlopt_srand((unsigned long) atoi(optarg));
      break;
    case 'a':
      c = atoi(optarg);
      if (c < 0 || c >= NLOPT_NUM_ALGORITHMS) {
	fprintf(stderr, "testopt: invalid algorithm %d\n", c);
	listalgs(stderr);
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
    case 'i':
      iterations = atoi(optarg);
      break;
    case 't':
      maxtime = atof(optarg);
      break;
    case 'x':
      xtol_rel = atof(optarg);
      break;
    case 'X':
      xtol_abs = atof(optarg);
      break;
    case 'f':
      ftol_rel = atof(optarg);
      break;
    case 'F':
      ftol_abs = atof(optarg);
      break;
    case 'm':
      fmin_max = atof(optarg);
      break;
    case '0':
      xinit_tol = atof(optarg);
      break;
    default:
      fprintf(stderr, "harminv: invalid argument -%c\n", c);
      usage(stderr);
      return EXIT_FAILURE;
    }
  
  return EXIT_SUCCESS;
}
