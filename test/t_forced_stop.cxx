// An exception thrown by the objective stops the optimization (issue #229).
// The point where it was thrown must not become the reported optimum, also
// when maximizing: the reported value must be finite and equal to the
// objective at the returned x.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <nlopt.hpp>

static int count = 0;
static double sign = 1; // -1 when maximizing

static double value(const std::vector<double> &x)
{
  return sign * ((x[0] - 1) * (x[0] - 1) + (x[1] - 2) * (x[1] - 2));
}

double objective(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
  (void)data;
  if (++count == 20)
    throw nlopt::forced_stop();
  if (!grad.empty()) {
    grad[0] = sign * 2 * (x[0] - 1);
    grad[1] = sign * 2 * (x[1] - 2);
  }
  return value(x);
}

int main() {
  const nlopt::algorithm algorithms[] = {nlopt::LN_COBYLA, nlopt::LN_BOBYQA, nlopt::LN_NELDERMEAD,
                                         nlopt::LN_SBPLX, nlopt::LN_PRAXIS, nlopt::GN_DIRECT_L,
                                         nlopt::GN_ISRES, nlopt::LD_MMA};
  int failures = 0;
  for (int maximize = 0; maximize < 2; ++maximize) {
    sign = maximize ? -1 : 1;
    for (nlopt::algorithm algorithm : algorithms) {
      nlopt::opt opt(algorithm, 2);
      opt.set_lower_bounds(std::vector<double>(2, -5.0));
      opt.set_upper_bounds(std::vector<double>(2, 5.0));
      if (maximize)
        opt.set_max_objective(objective, NULL);
      else
        opt.set_min_objective(objective, NULL);
      opt.set_maxeval(1000);
      std::vector<double> x(2, 0.0);
      double opt_f;
      count = 0;
      bool stopped = false;
      try {
        opt.optimize(x, opt_f);
      } catch (nlopt::forced_stop &) {
        stopped = true;
      }
      double reported = opt.last_optimum_value();
      bool ok = stopped && std::isfinite(reported) && reported == value(x);
      std::cout << (maximize ? "max " : "min ") << nlopt::algorithm_name(algorithm) << ": "
                << (stopped ? "forced_stop" : "no exception") << ", reported " << reported
                << ", objective at returned x " << value(x) << (ok ? "" : "  FAILED") << std::endl;
      if (!ok) ++failures;
    }
  }
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
