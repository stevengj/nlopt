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
    grad[0] = 2*x[0];
    grad[1] = 2*x[1];
  }
  if (std::isnan(x[0]))
    throw std::invalid_argument("nan");
  return -x[0]*x[0]-x[1]*x[1];
}


int main() {
  nlopt::srand(0);
  for (int repeat = 0; repeat < 100; ++repeat)
  {
    std::cout << "repeat="<<repeat<<std::endl;
    nlopt::opt opt(nlopt::LN_PRAXIS, 2);
    std::vector<double> lb = {-1.0, -10.0};
    std::vector<double> ub = {10.0, 1.0};
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    opt.set_min_objective(myvfunc, NULL);
    opt.set_ftol_rel(1e-6);
    std::vector<double> x = {0.5, 0.5};
    double minf;
    try{
      opt.optimize(x, minf);
      std::cout << "found optimum at f(" << x[0] << "," << x[1] << ") = "
                << std::setprecision(10) << minf <<std::endl;
    }
    catch(std::exception &e) {
      std::cerr << "nlopt failed: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}