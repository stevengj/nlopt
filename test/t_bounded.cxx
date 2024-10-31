#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <nlopt.hpp>


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
  nlopt::opt opt(argc < 2 ? nlopt::GD_STOGO : (nlopt::algorithm)atoi(argv[1]), 2);
  opt.set_lower_bounds(0.0);
  opt.set_upper_bounds(1.0);
  opt.set_max_objective(myvfunc, nullptr);
  opt.set_maxeval(100000);
  std::vector<double> x = {0.5, 0.5};
  double minf = 0.0;
  std::cout << "algo: " << opt.get_algorithm_name() << std::endl;
  try{
    opt.optimize(x, minf);
    std::cout << "found minimum at f(" << x[0] << "," << x[1] << ") = "
              << std::setprecision(10) << minf <<std::endl;
    return std::fabs(minf - 2.0) < 2e-2 ? EXIT_SUCCESS : EXIT_FAILURE;
  }
  catch(std::exception &e) {
    std::cerr << "nlopt failed: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
}
