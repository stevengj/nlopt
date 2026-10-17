#include <nlopt.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    double lower[3], upper[3];
    unsigned evaluations, constraints, outside;
} problem_data;

static int inside(unsigned n, const double *x, const problem_data *data)
{
    unsigned i;
    for (i = 0; i < n; ++i)
        if (!(data->lower[i] <= x[i] && x[i] <= data->upper[i]))
            return 0;
    return 1;
}

static double value(const double *x)
{
    return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
}

static double objective(unsigned n, const double *x, double *grad, void *ptr)
{
    problem_data *data = (problem_data *)ptr;
    (void)grad;
    ++data->evaluations;
    data->outside += !inside(n, x, data);
    return value(x);
}

static double constraint(unsigned n, const double *x, double *grad, void *ptr)
{
    problem_data *data = (problem_data *)ptr;
    (void)grad;
    ++data->constraints;
    data->outside += !inside(n, x, data);
    return x[1] - x[2];
}

int main(int argc, char **argv)
{
    int test_case = argc > 1 ? atoi(argv[1]) : 0;
    int constraint_kind = argc > 2 ? atoi(argv[2]) : 0;
    problem_data data = {{-60, 1e-4, 1e-4}, {60, 60, 60}, 0, 0, 0};
    double x[3] = {0, 1e-4, 1e-4}, minf, steps[3] = {1, 1, 1};
    unsigned i, repeat;
    nlopt_opt opt;
    if (test_case < 0 || test_case > 7 || constraint_kind < 0 || constraint_kind > 2)
        return EXIT_FAILURE;
    if (test_case >= 4) {
        for (i = 1; i < 3; ++i) {
            data.lower[i] = -60;
            data.upper[i] = -1e-4;
            x[i] = -1e-4;
        }
    }
    if (test_case % 4 >= 2) {
        steps[0] = 30;
        steps[1] = steps[2] = test_case % 4 == 2 ? 0.749925 : -0.749925;
    }
    opt = nlopt_create(NLOPT_LN_COBYLA, 3);
    if (!opt)
        return EXIT_FAILURE;
    if (nlopt_set_lower_bounds(opt, data.lower) < 0 ||
        nlopt_set_upper_bounds(opt, data.upper) < 0 ||
        nlopt_set_min_objective(opt, objective, &data) < 0 ||
        nlopt_set_xtol_abs1(opt, 1e-2) < 0 ||
        nlopt_set_xtol_rel(opt, 1e-3) < 0 ||
        nlopt_set_maxeval(opt, 200) < 0 ||
        (test_case % 4 && nlopt_set_initial_step(opt, steps) < 0) ||
        (constraint_kind == 1 && nlopt_add_inequality_constraint(opt, constraint, &data, 1e-12) < 0) ||
        (constraint_kind == 2 && nlopt_add_equality_constraint(opt, constraint, &data, 1e-12) < 0)) {
        nlopt_destroy(opt);
        return EXIT_FAILURE;
    }
    for (repeat = 0; repeat < 2; ++repeat) {
        nlopt_result result;
        data.evaluations = data.constraints = data.outside = 0;
        result = nlopt_optimize(opt, x, &minf);
        if (result <= 0 || !data.evaluations || (constraint_kind && !data.constraints) ||
            data.outside || !inside(3, x, &data) || fabs(value(x) - minf) > 16 * DBL_EPSILON) {
            fprintf(stderr, "COBYLA case %d constraint %d repeat %u: result=%d, outside=%u, x=[%.17g, %.17g, %.17g]\n",
                    test_case, constraint_kind, repeat, result, data.outside, x[0], x[1], x[2]);
            nlopt_destroy(opt);
            return EXIT_FAILURE;
        }
    }
    nlopt_destroy(opt);
    return EXIT_SUCCESS;
}
