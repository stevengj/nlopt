#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <nlopt.hpp>

/**
 * Test exception handling. Contributed by Kevin Kofler <kofler@dagopt.com>.
 * This tests a test case known to fail (t_bounded.cxx on LD_SLSQP) with
 * set_exceptions_enabled set to true if the argument is 1, false if the
 * argument is 0, and verifies that the exception is thrown resp. not thrown,
 * as expected.
 */

double myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
  (void)data;
  if (!grad.empty()) {
    grad[0] = 2.0 * x[0];
    grad[1] = 2.0 * x[1];
  }
  return x[0]*x[0]+x[1]*x[1];
}

int main(int argc, char *argv[]) {
  nlopt::opt opt(nlopt::LD_SLSQP, 2);
  if (argc >= 2) {
    bool enable = (bool) atoi(argv[1]);
    opt.set_exceptions_enabled(enable);
    if (opt.get_exceptions_enabled() != enable) {
      std::cerr << "set_exceptions_enabled(" << (enable ? "true" : "false")
                << ") failed" << std::endl;
      return EXIT_FAILURE;
    }
  }
  opt.set_lower_bounds(0.0);
  opt.set_upper_bounds(1.0);
  opt.set_max_objective(myvfunc, nullptr);
  opt.set_maxeval(100000);
  std::vector<double> x = {0.5, 0.5};
  double minf = 0.0;
  std::cout << "exceptions enabled: "
            << (opt.get_exceptions_enabled() ? "true" : "false")
            << std::endl;
  try{
    opt.optimize(x, minf);
    std::cout << "no exception was thrown" << std::endl;
    return opt.get_exceptions_enabled() ? EXIT_FAILURE : EXIT_SUCCESS;
  }
  catch(std::exception &e) {
    std::cerr << "exception: " << e.what() << std::endl;
    return opt.get_exceptions_enabled() ? EXIT_SUCCESS : EXIT_FAILURE;
  }
}

