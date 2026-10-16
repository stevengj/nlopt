// The unscaled DIRECT variants (issue #276). On the unit square the scaling of
// GN_DIRECT, GN_DIRECT_L and GN_DIRECT_L_RAND is the identity, so each of them
// must evaluate the same points as its _NOSCAL variant.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <nlopt.hpp>

static std::vector<std::vector<double> > trace; // all evaluated points

// Branin function, mapped from [-5,10] x [0,15] to the unit square
double branin(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
  (void)grad;
  (void)data;
  trace.push_back(x);
  const double pi = 3.14159265358979323846;
  double u = 15 * x[0] - 5, v = 15 * x[1];
  double t = v - 5.1 / (4 * pi * pi) * u * u + 5 / pi * u - 6;
  return t * t + 10 * (1 - 1 / (8 * pi)) * std::cos(u) + 10;
}

// symmetric in x[0] and x[1], so that rectangles tie and the randomized
// variant makes random choices
double symmetric(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
  (void)grad;
  (void)data;
  trace.push_back(x);
  return (x[0] - 0.3) * (x[0] - 0.3) + (x[1] - 0.3) * (x[1] - 0.3);
}

static std::vector<std::vector<double> > run(nlopt::algorithm algorithm, nlopt::vfunc f)
{
  nlopt::srand(1);
  nlopt::opt opt(algorithm, 2);
  opt.set_lower_bounds(std::vector<double>(2, 0.0));
  opt.set_upper_bounds(std::vector<double>(2, 1.0));
  opt.set_min_objective(f, NULL);
  opt.set_maxeval(300);
  std::vector<double> x(2, 0.5);
  double minf;
  trace.clear();
  opt.optimize(x, minf);
  return trace;
}

int main() {
  const nlopt::algorithm scaled[] = {nlopt::GN_DIRECT, nlopt::GN_DIRECT_L, nlopt::GN_DIRECT_L_RAND};
  const nlopt::algorithm unscaled[] = {nlopt::GN_DIRECT_NOSCAL, nlopt::GN_DIRECT_L_NOSCAL,
                                       nlopt::GN_DIRECT_L_RAND_NOSCAL};
  const nlopt::vfunc funcs[] = {branin, symmetric};
  const char *names[] = {"branin", "symmetric"};
  int failures = 0;
  for (int i = 0; i < 3; ++i) {
    for (int k = 0; k < 2; ++k) {
      std::vector<std::vector<double> > t1 = run(scaled[i], funcs[k]);
      std::vector<std::vector<double> > t2 = run(unscaled[i], funcs[k]);
      bool same = t1 == t2;
      std::cout << nlopt::algorithm_name(unscaled[i]) << ", " << names[k] << ": " << t2.size()
                << " evaluations, " << (same ? "same points as" : "different points from")
                << " the scaled variant" << std::endl;
      if (!same) ++failures;
    }
  }
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
