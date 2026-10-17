// COBYLA must return when its trust-region radius becomes so small that the
// step computed in trstlp is nan (issue #370). Before, it looped forever
// without evaluating the objective, ignoring maxeval and maxtime.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <nlopt.hpp>

static int count = 0;

double quadratic(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
  (void)grad;
  (void)data;
  ++count;
  return (x[0] - 1) * (x[0] - 1) + 10 * (x[1] - 2) * (x[1] - 2);
}

// the example from issue #370
double issue370(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
  (void)grad;
  (void)data;
  ++count;
  double t = 2 - std::cos(x[0]) + x[1] * x[1];
  return t * t;
}

static bool run(nlopt::vfunc f, double fopt, int maxeval)
{
  nlopt::opt opt(nlopt::LN_COBYLA, 2);
  opt.set_min_objective(f, NULL);
  opt.set_maxeval(maxeval);
  std::vector<double> x(2, 0.0);
  double minf = HUGE_VAL;
  count = 0;
  try {
    opt.optimize(x, minf);
  } catch (nlopt::roundoff_limited &) {
    minf = opt.last_optimum_value();
  }
  std::cout << "f(" << x[0] << "," << x[1] << ") = " << minf << " after " << count << " evaluations" << std::endl;
  return std::fabs(minf - fopt) < 1e-8 && count <= maxeval;
}

int main() {
  bool ok1 = run(quadratic, 0.0, 2000);
  bool ok2 = run(issue370, 1.0, 668);
  bool ok = ok1 && ok2;
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
