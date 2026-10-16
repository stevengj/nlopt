// Initial-step getters and setters of nlopt::opt (issue #545).

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
#include <nlopt.hpp>

// The default initial step depends on x, so there is no getter without x.
template <typename T, typename = void>
struct has_get_initial_step_without_x : std::false_type {};
template <typename T>
struct has_get_initial_step_without_x<T, decltype((void)std::declval<const T &>().get_initial_step())>
  : std::true_type {};
static_assert(!has_get_initial_step_without_x<nlopt::opt>::value,
              "get_initial_step() needs x to compute the default step");

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

  // default step, computed from x and the bounds
  std::vector<double> dx(2);
  opt.get_initial_step(x, dx);
  check(dx[0] > 0 && dx[1] > 0, "default step is positive");
  check(equal(opt.get_initial_step_(x), dx), "get_initial_step_(x) matches get_initial_step(x, dx)");

  opt.set_initial_step(0.25);
  check(equal(opt.get_initial_step_(x), std::vector<double>{0.25, 0.25}), "set_initial_step(double)");

  opt.set_initial_step(std::vector<double>{0.1, 0.2});
  check(equal(opt.get_initial_step_(x), std::vector<double>{0.1, 0.2}), "set_initial_step(vector)");

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
