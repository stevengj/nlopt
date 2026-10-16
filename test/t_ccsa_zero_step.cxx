// MMA and CCSA on f(x) = x'x with a relative tolerance on f (issue #635).
// Since min f = 0, ftol_rel only stops the run once f underflows to 0. Before
// that the steps become so small that their squares underflow, and the
// update of rho must not divide by zero.

#include <cfenv>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <nlopt.hpp>

double sumsq(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
  (void)data;
  if (!grad.empty()) {
    grad[0] = 2 * x[0];
    grad[1] = 2 * x[1];
  }
  return x[0] * x[0] + x[1] * x[1];
}

int main(int argc, char *argv[]) {
  nlopt::opt opt(argc < 2 ? "LD_MMA" : argv[1], 2);
  opt.set_lower_bounds(std::vector<double>(2, -1.0));
  opt.set_upper_bounds(std::vector<double>(2, 1.0));
  opt.set_min_objective(sumsq, NULL);
  opt.set_ftol_rel(1e-2);
  std::vector<double> x(2, 0.5);
  double minf;

  std::feclearexcept(FE_ALL_EXCEPT);
  nlopt::result result = opt.optimize(x, minf);
  bool divbyzero = std::fetestexcept(FE_DIVBYZERO) != 0;

  std::cout << opt.get_algorithm_name() << ": result " << result
            << ", f(" << x[0] << "," << x[1] << ") = " << minf
            << ", division by zero: " << (divbyzero ? "yes" : "no") << std::endl;
  return (result > 0 && minf < 1e-10 && !divbyzero) ? EXIT_SUCCESS : EXIT_FAILURE;
}
