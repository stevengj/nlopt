#include <cstdlib>
#include <iostream>
#include <vector>

#include <nlopt.hpp>

static double linear(unsigned, const double *x, double *, void *)
{
    return x[0] + 1.25 * x[1];
}

int main()
{
    nlopt::opt opt(nlopt::LN_COBYLA, 2);
    opt.set_lower_bounds(0.0);
    opt.set_upper_bounds(1.0);
    opt.set_min_objective(linear, NULL);
    opt.set_stopval(1e-8);
    opt.set_maxeval(200);

    // The minimum is at (0, 0). Reflecting a trial step away from a lower
    // bound can leave the solver in the interior until it runs out of evaluations.
    std::vector<double> x(2);
    x[0] = 0.44;
    x[1] = 0.5;
    double value;
    opt.optimize(x, value);
    if (value <= opt.get_stopval())
        return EXIT_SUCCESS;

    std::cerr << "COBYLA stopped at (" << x[0] << ", " << x[1]
              << "), objective " << value << std::endl;
    return EXIT_FAILURE;
}
