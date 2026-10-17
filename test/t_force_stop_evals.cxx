// After the objective calls force_stop() at evaluation K, the algorithm must
// return NLOPT_FORCED_STOP without evaluating the objective again, for several
// K (the first evaluation, the initial population, later iterations). A
// finite reported optimum must equal the objective at the returned x.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <nlopt.hpp>

struct probe {
  nlopt::opt *opt;
  int count, stop_at;
};

static double value(const std::vector<double> &x)
{
  return (x[0] - 1) * (x[0] - 1) + (x[1] - 0.5) * (x[1] - 0.5) + 0.1 * std::sin(5 * x[0]);
}

static double objective(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
  probe *p = static_cast<probe *>(data);
  if (++p->count == p->stop_at)
    p->opt->force_stop();
  if (!grad.empty()) {
    grad[0] = 2 * (x[0] - 1) + 0.5 * std::cos(5 * x[0]);
    grad[1] = 2 * (x[1] - 0.5);
  }
  return value(x);
}

int main(int argc, char *argv[])
{
  if (argc < 2) {
    std::cerr << "usage: t_force_stop_evals <algorithm number>" << std::endl;
    return EXIT_FAILURE;
  }
  const nlopt::algorithm algorithm = (nlopt::algorithm)atoi(argv[1]);
  const int stops[] = {1, 2, 5, 15, 40};
  int failures = 0;
  for (int stop_at : stops) {
    nlopt::opt opt(algorithm, 2);
    probe p = {&opt, 0, stop_at};
    opt.set_lower_bounds(std::vector<double>(2, -2.0));
    opt.set_upper_bounds(std::vector<double>(2, 2.0));
    opt.set_min_objective(objective, &p);
    opt.set_maxeval(1000);
    nlopt::srand(7);
    std::vector<double> x(2, 0.5);
    double opt_f;
    bool stopped = false;
    try {
      opt.optimize(x, opt_f);
    } catch (nlopt::forced_stop &) {
      stopped = true;
    }
    double reported = opt.last_optimum_value();
    bool ok;
    if (p.count < stop_at)  // the run ended before evaluation K
      ok = !stopped;
    else
      ok = stopped && p.count == stop_at && (!std::isfinite(reported) || reported == value(x));
    std::cout << nlopt::algorithm_name(algorithm) << ", stop at evaluation " << stop_at << ": "
              << p.count << " evaluations, " << (stopped ? "forced_stop" : "no forced_stop")
              << ", reported " << reported << (ok ? "" : "  FAILED") << std::endl;
    if (!ok)
      ++failures;
  }
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
