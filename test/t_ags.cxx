#include <nlopt.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

struct objective_data {
  nlopt_opt opt;
  int test_case;
  unsigned evaluations;
  double best;
  bool received_gradient;
};

static double objective_value(const double *x, int test_case)
{
  if (test_case == 4)
    return 3.25;
  return 100 * (x[1] - x[0]*x[0]) * (x[1] - x[0]*x[0])
    + (1 - x[0]) * (1 - x[0]);
}

static double objective(unsigned n, const double *x, double *grad, void *data)
{
  (void)n;
  objective_data &d = *static_cast<objective_data *>(data);
  d.received_gradient = d.received_gradient || grad != NULL;
  const double value = objective_value(x, d.test_case);
  d.best = std::min(d.best, value);
  ++d.evaluations;
  if (d.test_case == 1 && d.evaluations == 7)
    nlopt_force_stop(d.opt);
  return value;
}

static double constraint(unsigned n, const double *x, double *grad, void *data)
{
  (void)n;
  (void)grad;
  (void)data;
  return x[0] + x[1] - 0.25;
}

int main(int argc, char *argv[])
{
  const int test_case = argc > 1 ? std::atoi(argv[1]) : 0;
  if (test_case < 0 || test_case > 4)
    return EXIT_FAILURE;

  nlopt_opt opt = nlopt_create(NLOPT_GN_AGS, 2);
  if (!opt)
    return EXIT_FAILURE;
  objective_data data = {opt, test_case, 0,
                        std::numeric_limits<double>::infinity(), false};
  const int maxeval = test_case == 1 ? 100 : test_case == 2 ? 6 : test_case == 3 ? 8 : 7;
  if (nlopt_set_lower_bounds1(opt, -10) < 0 ||
      nlopt_set_upper_bounds1(opt, 5) < 0 ||
      nlopt_set_min_objective(opt, objective, &data) < 0 ||
      nlopt_set_maxeval(opt, maxeval) < 0 ||
      (test_case == 2 && nlopt_add_inequality_constraint(opt, constraint, NULL, 0) < 0)) {
    nlopt_destroy(opt);
    return EXIT_FAILURE;
  }

  double x[] = {0, 1};
  double minf = 0;
  const nlopt_result result = nlopt_optimize(opt, x, &minf);
  const nlopt_result expected = test_case == 1 ? NLOPT_FORCED_STOP : NLOPT_MAXEVAL_REACHED;
  const double tolerance = 1e-12 * std::max(1.0, std::abs(data.best));
  const bool passed = result == expected && data.evaluations > 0 &&
    nlopt_get_numevals(opt) == static_cast<int>(data.evaluations) &&
    !data.received_gradient && std::isfinite(minf) &&
    std::abs(minf - data.best) <= tolerance &&
    std::abs(objective_value(x, test_case) - minf) <= tolerance &&
    x[0] >= -10 && x[0] <= 5 && x[1] >= -10 && x[1] <= 5 &&
    (test_case != 2 || constraint(2, x, NULL, NULL) <= 0);

  if (!passed)
    std::cerr << "AGS case " << test_case << ": status=" << result
              << ", evaluations=" << data.evaluations << ", returned=" << minf
              << ", best evaluated=" << data.best << std::endl;
  nlopt_destroy(opt);
  return passed ? EXIT_SUCCESS : EXIT_FAILURE;
}
