#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <nlopt.h>

/* COBYLA evaluates these functions without gradients. */
static double objective(unsigned n, const double *x, double *grad, void *data)
{
    return x[0] * x[0];
}

static void constraint(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data)
{
    result[0] = *(const double *) data - x[0];
}

static int check(int equality)
{
    const unsigned count = (unsigned) (SIZE_MAX / sizeof(double)) + 1;
    double first = 1, second = 2, x = 3, f;
    nlopt_result (NLOPT_STDCALL *add)(nlopt_opt, unsigned, nlopt_mfunc, void *, const double *) =
        equality ? nlopt_add_equality_mconstraint : nlopt_add_inequality_mconstraint;
    nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, 1);
    int ok = opt && nlopt_set_min_objective(opt, objective, NULL) > 0
        && nlopt_set_xtol_rel(opt, 1e-9) > 0 && nlopt_set_maxeval(opt, 100) > 0
        && nlopt_add_inequality_mconstraint(opt, 1, constraint, &first, NULL) > 0;
    if (ok) {
        /* Reject the oversized addition and retain the original x >= 1. */
        ok = add(opt, count, constraint, &second, NULL) == NLOPT_OUT_OF_MEMORY
            && nlopt_optimize(opt, &x, &f) > 0 && fabs(x - first) < 1e-6;
    }
    if (ok) {
        x = 3;
        ok = add(opt, 1, constraint, &second, NULL) > 0
            && nlopt_optimize(opt, &x, &f) > 0 && fabs(x - second) < 1e-6;
    }
    nlopt_destroy(opt);
    if (!ok)
        fprintf(stderr, "%s constraint tolerance size check failed\n", equality ? "equality" : "inequality");
    return ok;
}

int main(void)
{
    return check(0) && check(1) ? 0 : 1;
}
