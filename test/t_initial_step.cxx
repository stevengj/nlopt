// Initial-step getters and setters of nlopt::opt (issue #545).

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <nlopt.hpp>

static int failures = 0;

static void check(bool ok, const char *what)
{
  if (!ok) {
    std::cout << "FAILED: " << what << std::endl;
    ++failures;
  }
}

static bool equal(const std::vector<double> &a, const std::vector<double> &b)
{
  return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin());
}

int main() {
  nlopt::opt opt(nlopt::LN_NELDERMEAD, 2);
  opt.set_lower_bounds(std::vector<double>{-1.0, 0.0});
  opt.set_upper_bounds(std::vector<double>{1.0, 4.0});
  const std::vector<double> x{0.5, 1.0};

  // no step set: the default step is computed from x and the bounds
  std::vector<double> dx(2);
  opt.get_initial_step(x, dx);
  check(dx[0] > 0 && dx[1] > 0, "default step is positive");
  check(equal(opt.get_initial_step_(x), dx), "get_initial_step_(x) matches get_initial_step(x, dx)");

  // no step set and no x: the getter cannot compute the default step
  std::string message;
  try {
    opt.get_initial_step();
  } catch (std::invalid_argument &e) {
    message = e.what();
  }
  check(message.find("depends on x") != std::string::npos,
        "get_initial_step() without a set step explains that the default depends on x");

  // a step that was set can be read without x
  opt.set_initial_step(0.25);
  check(equal(opt.get_initial_step(), std::vector<double>{0.25, 0.25}), "get_initial_step() after set_initial_step(double)");
  check(equal(opt.get_initial_step_(x), std::vector<double>{0.25, 0.25}), "get_initial_step_(x) after set_initial_step(double)");

  opt.set_initial_step(std::vector<double>{0.1, 0.2});
  std::vector<double> v(2);
  opt.get_initial_step(v);
  check(equal(v, std::vector<double>{0.1, 0.2}), "get_initial_step(v) after set_initial_step(vector)");

  bool thrown = false;
  try {
    opt.set_initial_step(std::vector<double>{0.1, 0.2, 0.3});
  } catch (std::invalid_argument &) {
    thrown = true;
  }
  check(thrown, "set_initial_step with the wrong dimension throws");

  std::cout << (failures ? "initial step tests failed" : "initial step tests passed") << std::endl;
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
