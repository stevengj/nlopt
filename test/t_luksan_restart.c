#include <nlopt.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* The restart interval 999*n exceeds INT_MAX just above this dimension.
   A single history vector keeps the large-dimension cases bounded in memory. */
static double quadratic(unsigned n, const double *x, double *grad, void *data)
{
    double value = 0;
    unsigned i;
    ++*(unsigned *)data;
    for (i = 0; i < n; ++i) {
        double d = x[i] - 0.25;
        value += d * d;
        if (grad) grad[i] = 2 * d;
    }
    return value;
}

int main(void)
{
    const nlopt_algorithm algorithms[] = {NLOPT_LD_LBFGS, NLOPT_LD_VAR1, NLOPT_LD_VAR2,
        NLOPT_LD_TNEWTON, NLOPT_LD_TNEWTON_RESTART,
        NLOPT_LD_TNEWTON_PRECOND, NLOPT_LD_TNEWTON_PRECOND_RESTART};
    const unsigned dimensions[] = {4, INT_MAX / 999, INT_MAX / 999 + 1};
    unsigned a, d, i, failures = 0;
    for (a = 0; a < sizeof(algorithms) / sizeof(algorithms[0]); ++a)
        for (d = 0; d < sizeof(dimensions) / sizeof(dimensions[0]); ++d) {
            unsigned n = dimensions[d], calls = 0;
            double *x = (double *)malloc(sizeof(double) * n), value, error = 0;
            nlopt_opt opt = nlopt_create(algorithms[a], n);
            nlopt_result result;
            if (!x || !opt) return 2;
            for (i = 0; i < n; ++i) x[i] = 1;
            if (nlopt_set_min_objective(opt, quadratic, &calls) < 0
                || nlopt_set_vector_storage(opt, 1) < 0
                || nlopt_set_maxeval(opt, 20) < 0
                || nlopt_set_maxtime(opt, 10) < 0
                || nlopt_set_xtol_rel(opt, 1e-10) < 0) return 3;
            result = nlopt_optimize(opt, x, &value);
            for (i = 0; i < n; ++i) if (fabs(x[i] - 0.25) > error) error = fabs(x[i] - 0.25);
            printf("algorithm=%s n=%u result=%d calls=%u f=%.17g max_error=%.17g\n",
                   nlopt_algorithm_to_string(algorithms[a]), n, result, calls, value, error);
            fflush(stdout);
            failures += result <= 0 || !(error <= 1e-6) || !(value <= 1e-6);
            nlopt_destroy(opt);
            free(x);
        }
    printf("failures=%u cases=21\n", failures);
    return failures ? 1 : 0;
}
