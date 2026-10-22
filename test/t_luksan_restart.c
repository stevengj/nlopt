#include <nlopt.h>
#include "luksan.h"
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Exercise large restart intervals without allocating dimension-sized vectors. */
static unsigned check_intervals(void)
{
    const struct { int n, ires1, ires2; long long expected; } cases[] = {
        {0, 999, 0, 0},
        {1, 999, 0, 999},
        {INT_MAX / 999, 999, 0, 2147483367LL},
        {INT_MAX / 999 + 1, 999, 0, 2147484366LL},
        {INT_MAX, 999, 0, 2145336163353LL},
        {INT_MAX, INT_MAX, INT_MAX, 4611686016279904256LL}
    };
    unsigned i, failures = 0;
    for (i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
        long long actual = luksan_restart_interval(cases[i].n, cases[i].ires1, cases[i].ires2);
        if (actual != cases[i].expected) {
            fprintf(stderr, "incorrect restart interval for case %u\n", i);
            ++failures;
        }
    }
    return failures;
}

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

int main(int argc, char **argv)
{
    const nlopt_algorithm algorithms[] = {NLOPT_LD_LBFGS, NLOPT_LD_VAR1, NLOPT_LD_VAR2,
        NLOPT_LD_TNEWTON, NLOPT_LD_TNEWTON_RESTART,
        NLOPT_LD_TNEWTON_PRECOND, NLOPT_LD_TNEWTON_PRECOND_RESTART};
    const unsigned dimensions[] = {4, INT_MAX / 999, INT_MAX / 999 + 1};
    /* Large public-API regressions remain available as an opt-in stress test. */
    unsigned ndims = argc == 2 && !strcmp(argv[1], "--large") ? 3 : 1;
    unsigned a, d, i, failures = check_intervals();
    for (a = 0; a < sizeof(algorithms) / sizeof(algorithms[0]); ++a)
        for (d = 0; d < ndims; ++d) {
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
    printf("failures=%u cases=%u\n", failures, 6 + 7 * ndims);
    return failures ? 1 : 0;
}
