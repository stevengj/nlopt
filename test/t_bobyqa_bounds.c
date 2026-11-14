#include <nlopt.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    double lower[3], upper[3];
    unsigned evaluations, outside;
} objective_data;

static double objective_value(const double *x)
{
    return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
}

static double objective(unsigned n, const double *x, double *grad, void *ptr)
{
    objective_data *data = (objective_data *)ptr;
    unsigned i;
    (void)grad;
    ++data->evaluations;
    for (i = 0; i < n; ++i)
        if (x[i] < data->lower[i] || x[i] > data->upper[i])
            ++data->outside;
    return objective_value(x);
}

static int inside(const double *x, const objective_data *data)
{
    unsigned i;
    for (i = 0; i < 3; ++i)
        if (x[i] < data->lower[i] || x[i] > data->upper[i])
            return 0;
    return 1;
}

int main(int argc, char **argv)
{
    int test_case = argc > 1 ? atoi(argv[1]) : 0;
    objective_data data = {{-60, 1e-4, 1e-4}, {60, 60, 60}, 0, 0};
    double x[3] = {0, 1, 1}, minf = 0, steps[3] = {1, 1, 1};
    nlopt_opt opt;
    nlopt_result first, second;
    int passed;
    unsigned i, first_outside;
    if (test_case < 0 || test_case > 3)
        return EXIT_FAILURE;
    opt = nlopt_create(NLOPT_LN_BOBYQA, 3);
    if (!opt)
        return EXIT_FAILURE;
    if (test_case == 1) {
        for (i = 1; i < 3; ++i) {
            data.lower[i] = -60;
            data.upper[i] = -1e-4;
            x[i] = -1;
        }
    }
    if (test_case == 3) {
        steps[0] = 30;
        steps[1] = steps[2] = -0.749925;
    }
    if (nlopt_set_lower_bounds(opt, data.lower) < 0 ||
        nlopt_set_upper_bounds(opt, data.upper) < 0 ||
        nlopt_set_min_objective(opt, objective, &data) < 0 ||
        nlopt_set_xtol_abs1(opt, 1e-2) < 0 ||
        nlopt_set_xtol_rel(opt, 1e-3) < 0 ||
        nlopt_set_maxeval(opt, 200) < 0 ||
        ((test_case == 2 || test_case == 3) && nlopt_set_initial_step(opt, steps) < 0)) {
        nlopt_destroy(opt);
        return EXIT_FAILURE;
    }
    first = nlopt_optimize(opt, x, &minf);
    passed = first > 0 && data.evaluations > 0 && !data.outside && inside(x, &data) &&
        fabs(objective_value(x) - minf) <= 16 * DBL_EPSILON;
    first_outside = data.outside;
    data.evaluations = data.outside = 0;
    second = nlopt_optimize(opt, x, &minf);
    passed = passed && second > 0 && data.evaluations > 0 && !data.outside &&
        inside(x, &data) && fabs(objective_value(x) - minf) <= 16 * DBL_EPSILON;
    if (!passed)
        fprintf(stderr, "BOBYQA case %d: first=%d second=%d, x=[%.17g, %.17g, %.17g], outside=%u/%u\n",
                test_case, first, second, x[0], x[1], x[2], first_outside, data.outside);
    nlopt_destroy(opt);
    return passed ? EXIT_SUCCESS : EXIT_FAILURE;
}
