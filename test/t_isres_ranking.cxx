// ISRES on problem g06 of Runarsson and Yao (2 variables, 2 inequality
// constraints, minimum -6961.81387558015). The differential variation must
// use the best survivor and the ranked survivors (issue #590). Over 300
// seeds, the runs need on average about 21000 evaluations to reach
// f* + 0.1 with the ranked survivors and about 32600 with storage order.
// The test requires fewer than 250000 evaluations for seeds 1 to 10.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <nlopt.hpp>

static double cube(double a) { return a * a * a; }
static double sqr(double a) { return a * a; }

static double objective(unsigned n, const double *x, double *grad, void *data)
{
  (void)n;
  (void)grad;
  (void)data;
  return cube(x[0] - 10) + cube(x[1] - 20);
}

static void constraints(unsigned m, double *g, unsigned n, const double *x, double *grad, void *data)
{
  (void)m;
  (void)n;
  (void)grad;
  (void)data;
  g[0] = -sqr(x[0] - 5) - sqr(x[1] - 5) + 100;
  g[1] = sqr(x[0] - 6) + sqr(x[1] - 5) - 82.81;
}

int main()
{
  const double fopt = -6961.81387558015;
  const double lb[2] = {13, 0}, ub[2] = {100, 100};
  int total = 0, solved = 0;
  for (int seed = 1; seed <= 10; ++seed) {
    nlopt::srand(seed);
    nlopt::opt opt(nlopt::GN_ISRES, 2);
    opt.set_lower_bounds(std::vector<double>(lb, lb + 2));
    opt.set_upper_bounds(std::vector<double>(ub, ub + 2));
    opt.set_min_objective(objective, NULL);
    opt.add_inequality_mconstraint(constraints, NULL, std::vector<double>(2, 1e-8));
    opt.set_stopval(fopt + 0.1);
    opt.set_maxeval(100000);
    std::vector<double> x(2);
    x[0] = 56.5;
    x[1] = 50;
    double minf;
    nlopt::result result = opt.optimize(x, minf);
    std::cout << "seed " << seed << ": result " << result << ", " << opt.get_numevals()
              << " evaluations, f = " << minf << std::endl;
    total += opt.get_numevals();
    if (result == nlopt::STOPVAL_REACHED)
      ++solved;
  }
  std::cout << solved << " of 10 runs reached f* + 0.1, " << total << " evaluations in total" << std::endl;
  return solved == 10 && total < 250000 ? EXIT_SUCCESS : EXIT_FAILURE;
}
