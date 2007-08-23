#include <stdio.h>
#include <math.h>

#include "testfuncs.h"
#include "config.h"

#define UNUSED(x) (void) x

static double sqr(double x) { return x * x; }

int testfuncs_verbose = 0;
int testfuncs_counter = 0;

static double testfuncs_status(int n, const double *x, double f)
{
     ++testfuncs_counter;
     if (testfuncs_verbose) {
	  int i;
	  printf("f_%d (%g", testfuncs_counter, x[0]);
	  for (i = 1; i < n; ++i) printf(", %g", x[i]);
	  printf(") = %g\n", f);
     }
     return f;
}

#define RETURN(f) return testfuncs_status(n, x, f);

/****************************************************************************/
static double rosenbrock_f(int n, const double *x, double *grad, void *data)
{
     double a = x[1] - x[0] * x[0], b = 1 - x[0];
     UNUSED(data);
     if (grad) {
	  grad[0] = -400 * a * x[0] - 2*b;
	  grad[1] = -200 * a;
     }
     RETURN(100 * sqr(a) + sqr(b));
}

static const double rosenbrock_lb[2] = {-2, -2};
static const double rosenbrock_ub[2] = {2, 2};
static const double rosenbrock_xmin[2] = {1, 1};

/****************************************************************************/
static double mccormic_f(int n, const double *x, double *grad, void *data)
{
     double a = x[0] + x[1], b = x[0] - x[1];
     UNUSED(data);
     if (grad) {
	  grad[0] = cos(a) + 2*b - 1.5;
	  grad[1] = cos(a) - 2*b + 2.5;
     }
     RETURN(sin(a) + sqr(b) - 1.5*x[0] + 2.5*x[1] + 1);
}

static const double mccormic_lb[2] = {-1.5, -3};
static const double mccormic_ub[2] = {4, 4};
static const double mccormic_xmin[2] = {-0.54719, 1.54719};

/****************************************************************************/
static double boxbetts_f(int n, const double *x, double *grad, void *data)
{
     int i;
     double f = 0;
     UNUSED(data);
     if (grad)
	  grad[0] = grad[1] = grad[2] = 0;
     for (i = 1; i <= 10; ++i) {
	  double e0 = exp(-0.1*i*x[0]);
	  double e1 = exp(-0.1*i*x[1]);
	  double e2 = exp(-0.1*i) - exp((double) -i);
	  double g = e0 - e1 - e2 * x[2];
	  f += sqr(g);
	  if (grad) {
	       grad[0] += (2 * g) * (-0.1*i*e0);
	       grad[1] += (2 * g) * (0.1*i*e1);
	       grad[2] += -(2 * g) * e2;
	  }
     }
     RETURN(f);
}

static const double boxbetts_lb[3] = {0.9,9,0.9};
static const double boxbetts_ub[3] = {1.2,11.2,1.2};
static const double boxbetts_xmin[3] = {1,10,1};

/****************************************************************************/
static double paviani_f(int n, const double *x, double *grad, void *data)
{
     int i;
     double f = 0, prod = 1;
     UNUSED(data);
     if (grad) for (i = 0; i < 10; ++i) grad[i] = 0;
     for (i = 0; i < 10; ++i) {
	  double ln1 = log(x[i] - 2);
	  double ln2 = log(10 - x[i]);
	  f += sqr(ln1) + sqr(ln2);
	  if (grad)
	       grad[i] += 2 * ln1 / (x[i] - 2) - 2 * ln2 / (10 - x[i]);
	  prod *= x[i];
     }
     f -= (prod = pow(prod, 0.2));
     if (grad)
	  for (i = 0; i < 10; ++i)
	       grad[i] -= 0.2 * prod / x[i];
     RETURN(f);
}

static const double paviani_lb[10] = {2.001,2.001,2.001,2.001,2.001,2.001,2.001,2.001,2.001,2.001};
static const double paviani_ub[10] = {9.999,9.999,9.999,9.999,9.999,9.999,9.999,9.999,9.999,9.999};
static const double paviani_xmin[10] = {9.350266,9.350266,9.350266,9.350266,9.350266,9.350266,9.350266,9.350266,9.350266,9.350266};

/****************************************************************************/
static double grosenbrock_f(int n, const double *x, double *grad, void *data)
{
     int i;
     double f = 0;
     UNUSED(data);
     if (grad) grad[0] = 0;
     for (i = 0; i < 29; ++i) {
	  double a = x[i+1] - x[i] * x[i], b = 1 - x[i];
	  if (grad) {
	       grad[i] += -400 * a * x[i] - 2*b;
	       grad[i+1] = -200 * a;
	  }
	  f += 100 * sqr(a) + sqr(b);
     }
     RETURN(f);
}

static const double grosenbrock_lb[30] = {-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30};
static const double grosenbrock_ub[30] = {30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30};
static const double grosenbrock_xmin[30] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

/****************************************************************************/
static double goldsteinprice_f(int n, const double *x, double *grad, void *data)
{
     double x0, x1, a1, a12, a2, b1, b12, b2;
     UNUSED(data);
     x0 = x[0]; x1 = x[1];
     a1 = x0+x1+1; a12 = sqr(a1);
     a2 = 19 - 14*x0 + 3*x0*x0 - 14*x1 + 6*x0*x1 + 3*x1*x1;
     b1 = 2*x0-3*x1; b12 = sqr(b1);
     b2 = 18 - 32*x0 + 12*x0*x0 + 48*x1 - 36*x0*x1 + 27*x1*x1;
     if (grad) {
	  grad[0] = (1 + a12 * a2) * (2 * b1 * 2 * b2 
				      + b12 * (-32 + 24*x0 - 36*x1))
	       + (2 * a1 * a2 + a12 * (-14 + 6*x0 + 6*x1)) * (30 + b12 * b2);
	  grad[1] = (1 + a12 * a2) * (2 * b1 * (-3) * b2
				      + b12 * (48 - 36*x0 + 54 * x1))
	       + (2 * a1 * a2 + a12 * (-14 + 6*x0 + 6*x1)) * (30 + b12 * b2);
     }
     RETURN((1 + a12 * a2) * (30 + b12 * b2));
}

static const double goldsteinprice_lb[2] = {-2, -2};
static const double goldsteinprice_ub[2] = {2, 2};
static const double goldsteinprice_xmin[2] = {0, -1};

/****************************************************************************/
/****************************************************************************/

const testfunc testfuncs[NTESTFUNCS] = {
     { rosenbrock_f, NULL, 1, 2,
       rosenbrock_lb, rosenbrock_ub, rosenbrock_xmin,
       0.0, "Rosenbrock function" },
     { mccormic_f, NULL, 1, 2,
       mccormic_lb, mccormic_ub, mccormic_xmin,
       -1.9133, "McCormic function" },
     { boxbetts_f, NULL, 1, 3,
       boxbetts_lb, boxbetts_ub, boxbetts_xmin,
       0.0, "Box and Betts exponential quadratic sum" },
     { paviani_f, NULL, 1, 10,
       paviani_lb, paviani_ub, paviani_xmin,
       -45.778470, "Paviani function" },
     { grosenbrock_f, NULL, 1, 30,
       grosenbrock_lb, grosenbrock_ub, grosenbrock_xmin,
       0.0, "Generalized Rosenbrock function" },
     { goldsteinprice_f, NULL, 1, 2,
       goldsteinprice_lb, goldsteinprice_ub, goldsteinprice_xmin,
       3.0, "Goldstein and Price function" }
};
