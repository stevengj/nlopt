#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>
#include <nlopt.hpp>

static double constant(const std::vector<double> &, std::vector<double> &, void *)
{
  return 0.0;
}

static std::pair<nlopt::result, int> run(double tolerance)
{
  nlopt::opt opt(nlopt::LN_COBYLA, 2);
  opt.set_min_objective(constant, NULL);
  // A small initial radius reaches the roundoff limit in a short run.
  opt.set_initial_step(std::vector<double>(2, std::ldexp(1.0, -1000)));
  opt.set_xtol_rel(tolerance);
  opt.set_xtol_abs(std::vector<double>(2, tolerance));
  opt.set_maxeval(250);
  std::vector<double> x(2, 0.0);
  double minf;
  nlopt::srand(17);
  try {
    opt.optimize(x, minf);
  } catch (const nlopt::roundoff_limited &) {
    // Roundoff is an expected reason to stop when both tolerances are disabled.
  }
  return std::make_pair(opt.last_optimize_result(), opt.get_numevals());
}

int main()
{
  const std::pair<nlopt::result, int> zero = run(0.0);
  const std::pair<nlopt::result, int> negative = run(-1.0);
  if (negative != zero) {
    std::cerr << "Zero tolerance returned " << zero.first << " after " << zero.second
              << " evaluations; negative tolerance returned " << negative.first
              << " after " << negative.second << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
