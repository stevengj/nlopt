#include <iostream>
#include <vector>
#include <iomanip> // setprecision
#include <string> // to_string
#include <numeric> // inner_product
#include <utility> // move
#include <algorithm> // copy
#include <cassert> // assert
#include <cmath> // sin, cos

#include <nlopt.hpp>

typedef std::vector<double> Vector;

class QuadraticForm {
  /*
   * This class represents a quadratic form (½ xᵀAx + bᵀx) and provides two
   * methods to compute itself and its gradient.
   *
   * Input is always a (const double*) which has appropriate dimension.
   */
  private:
    Vector A; // symmetric matrix
    Vector b; // bias
    size_t dimension;

    Vector compute_A_times_x(const double* x) const
    {
      Vector A_times_x(dimension);
      for (size_t row = 0; row < dimension; ++row) {
        A_times_x[row] = dot(x, A, row * dimension);
      }
      return A_times_x;
    }

    double compute_quadratic_term(const double* x) const
    {
      auto A_times_x = compute_A_times_x(x);
      return 0.5 * dot(x, A_times_x);
    }

    double compute_linear_term(const double* x) const
    {
      return dot(x, b);
    }

    double dot(const double* x, const Vector& y, const size_t offset = 0) const
    {
      return std::inner_product(
          x,
          x + dimension,
          std::begin(y) + offset,
          double(0));
    }

  public:
    QuadraticForm() = delete;
    QuadraticForm(Vector A_, Vector b_) :
      A(std::move(A_)),
      b(std::move(b_))
    {
      if (A.size() != b.size() * b.size()) {
        throw std::runtime_error(
            "[QuadraticForm] matrix and bias dimension mismatch: " +
            std::to_string(A.size()) +
            " != " +
            std::to_string(b.size() * b.size()));
      }
      dimension = b.size();
    }

    double compute_form(const double* x) const
    {
      // ½ xᵀAx + bᵀx
      return compute_quadratic_term(x) + compute_linear_term(x);
    }

    void compute_gradient_in_place(const double* x, double* grad_array) const
    {
      // matrix is assumed symmetric, hence gradient is (Ax + b)
      Vector grad_vec = compute_A_times_x(x);

      auto b_iterator = std::begin(b);
      for (auto it = std::begin(grad_vec); it != std::end(grad_vec); ) {
        *it++ += *b_iterator++;
      }

      std::copy(std::begin(grad_vec), std::end(grad_vec), grad_array);
    }

}; // QuadraticForm


class LinearRegression {
  private:
    QuadraticForm quadratic_form;

  public:
    LinearRegression() = delete;
    LinearRegression(QuadraticForm quadratic_form_) :
      quadratic_form(std::move(quadratic_form_)) {}

    double operator()(unsigned n, const double* x, double* grad) const
    {
      const double result = quadratic_form.compute_form(x);
      if (!!grad) {
        quadratic_form.compute_gradient_in_place(x, grad);
      }

      return result;
    }
}; // LinearRegression


class SineRegression {
  private:
    QuadraticForm quadratic_form;

  public:
    SineRegression() = delete;
    SineRegression(QuadraticForm quadratic_form_) :
      quadratic_form(std::move(quadratic_form_)) {}

    double operator()(unsigned n, const double* x, double* grad) const
    {
      const double form_result = quadratic_form.compute_form(x);
      const double result = sin(form_result);

      if (!!grad) {
        const double grad_coefficient = cos(form_result);
        quadratic_form.compute_gradient_in_place(x, grad);

        for (size_t i = 0; i < n; ++i) {
          grad[i] *= grad_coefficient;
        }
      }

      return result;
    }
}; // SineRegression



int main()
{
  // quadratic form is positive-definite for tau ∈ [-1, 2]
  constexpr double tau = 1.5;
  constexpr size_t dimension = 3;

  Vector A = {
      2,   -1,   tau,
     -1,    2,    -1,
    tau,   -1,     2 };
  Vector b = { -3, 8, 1 };

  assert(b.size() == dimension);
  assert(A.size() == dimension * dimension);

  QuadraticForm    form(std::move(A), std::move(b));
  LinearRegression objective_1(form);
  SineRegression   objective_2(std::move(form));

  nlopt::opt optimizer("LD_MMA", dimension);
  optimizer.set_xtol_rel(1e-4);
  optimizer.set_maxeval(1000);


  // linear regression optimization
  optimizer.set_min_objective(std::move(objective_1));
  double minimum = 0.0;
  Vector x0 = {-1.0, 0.0, 0.1};

  try {
    optimizer.optimize(x0, minimum);
  } catch (const std::exception& e) {
    std::cerr << "NLopt failed: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "-- Linear regression --" << std::endl;
  std::cout << "found minimum at f("
    << x0[0] << ","
    << x0[1] << ","
    << x0[2] << ") = "
    << std::setprecision(10) << minimum << std::endl << std::endl;


  // sine regression optimization
  optimizer.set_min_objective(std::move(objective_2));
  minimum = 0.0;
  x0 = {-1.0, 0.0, 0.1};

  try {
    optimizer.optimize(x0, minimum);
  } catch (const std::exception& e) {
    std::cerr << "NLopt failed: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "-- Sine regression --" << std::endl;
  std::cout << "found minimum at f("
    << x0[0] << ","
    << x0[1] << ","
    << x0[2] << ") = "
    << std::setprecision(10) << minimum << std::endl;

  return 0;
}

