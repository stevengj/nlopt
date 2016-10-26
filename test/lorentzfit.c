#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nlopt.h"
#include "nlopt-util.h"

typedef struct {
    int N;
    double *x, *y; /* array of N points (x,y) */
} lorentzdata;

static double sqr(double x) { return x * x; }

static int count = 0;

static double lorentzerr(int n, const double *p, double *grad, void *data)
{
    lorentzdata *d = (lorentzdata *) data;
    double A = p[0], w = p[1], G = p[2];
    double val = 0;
    int i;

    /* loop over the N data points (x,y) and accumulate errors */
    for (i = 0; i < d->N; ++i) {
        /* i-th point (x,y) */
        double x = d->x[i], y = d->y[i];
        /* compute f(x) using current solution p */
        double lor = A / (sqr(x - w) + G*G);
        /* fit function: sum of squared errors */
        val += sqr(y - lor);

        /* gradient of fit function wrt each param */
        if (grad) {
            double deninv =  1.0 / (sqr(x - w) + G*G);
            grad[0] += -2 * (y - lor) * deninv;
            grad[1] += 4*A * (w - x) * (y - lor) * sqr(deninv);
            grad[2] += 4*A * G * (y - lor) * sqr(deninv);
        }
    }

    ++count;
    /*printf("#%d: err(%g,%g,%g) = %g\n", count, A, w, G, val);*/

    return val;
}

int main(void)
{
    lorentzdata d;
    int i;
    /* true parameters of Lorentz function: A, w, G */
    double A = 1, w = 0, G = 1;
    /* control added noise */
    double noise = 0.01;
    /* lower/upper bounds for params */
    double lb[3] = {-HUGE_VAL,-HUGE_VAL,       0};
    double ub[3] = { HUGE_VAL, HUGE_VAL,HUGE_VAL};
    /* initial guess for A,w,G params */
    double p[3] = {0,1,2}, minf = 0;
    nlopt_result res;

    /* seed random number generator */
    nlopt_srand_time();

    /* generate N random points from Lorentz function (with params A,w,G) */
    /* https://en.wikipedia.org/wiki/Cauchy_distribution */
    d.N = 200;
    d.x = (double *) malloc(sizeof(double) * d.N * 2);
    d.y = d.x + d.N;
    for (i = 0; i < d.N; ++i) {
        d.x[i] = nlopt_urand(-0.5, 0.5) * 8*G + w; /* x */
        d.y[i] = A / (sqr(d.x[i] - w) + G*G);      /* y = Lorentz(x) */
        d.y[i] += 2*noise * nlopt_urand(-0.5,0.5); /* add some noise */
    }

    /* solve least-squares curve-fitting problem using different algorithms */
    /* (i.e find best A,w,G params to minimize sum of squared errors) */
    /* min_p sum_i((f(x(i),p) - y(i))^2) */
    printf("params: A=%g, w=%g, G=%g\n", A, w, G);

    count = 0;
    p[0] = 0; p[1] = 1; p[2] = 2;
    res = nlopt_minimize(NLOPT_LN_NEWUOA_BOUND, 3, lorentzerr, &d,
        lb, ub, p, &minf, -HUGE_VAL, 0,0, 1e-6,NULL, 0,0);
    printf("%4d minf=%g at A=%g, w=%g, G=%g\n", count, minf, p[0],p[1],p[2]);

    count = 0;
    p[0] = 0; p[1] = 1; p[2] = 2;
    res = nlopt_minimize(NLOPT_LN_COBYLA, 3, lorentzerr, &d,
        lb, ub, p, &minf, -HUGE_VAL, 0,0, 1e-6,NULL, 0,0);
    printf("%4d minf=%g at A=%g, w=%g, G=%g\n", count, minf, p[0],p[1],p[2]);

    count = 0;
    p[0] = 0; p[1] = 1; p[2] = 2;
    res = nlopt_minimize(NLOPT_LN_NELDERMEAD, 3, lorentzerr, &d,
        lb, ub, p, &minf, -HUGE_VAL, 0,0, 1e-6,NULL, 0,0);
    printf("%4d minf=%g at A=%g, w=%g, G=%g\n", count, minf, p[0],p[1],p[2]);

    count = 0;
    p[0] = 0; p[1] = 1; p[2] = 2;
    res = nlopt_minimize(NLOPT_LN_SBPLX, 3, lorentzerr, &d,
        lb, ub, p, &minf, -HUGE_VAL, 0,0, 1e-6,NULL, 0,0);
    printf("%4d minf=%g at A=%g, w=%g, G=%g\n", count, minf, p[0],p[1],p[2]);

    free(d.x);
    return 0;
}
