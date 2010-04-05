/* Copyright (c) 2007-2009 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

// C++ style wrapper around NLopt API
// nlopt.hpp is AUTOMATICALLY GENERATED from nlopt-in.hpp - edit the latter!

#include <nlopt.h>

#include <vector>
#include <stdexcept>
#include <new>

// convenience overloading for below (not in nlopt:: since has nlopt_ prefix)
inline nlopt_result nlopt_get_initial_step(const nlopt_opt opt, double *dx) {
      return nlopt_get_initial_step(opt, (const double *) NULL, dx);
}

namespace nlopt {

  //////////////////////////////////////////////////////////////////////
  // nlopt::* namespace versions of the C enumerated types
  //          AUTOMATICALLY GENERATED, DO NOT EDIT
  // GEN_ENUMS_HERE
  //////////////////////////////////////////////////////////////////////

  // virtual base class for objective function and constraints:
  class func {
  public:
    // should return function value, and set grad to gradient
    // (x and grad are length n)
    virtual double operator()(int n, const double *x, double *grad) = 0;

    // should return function value (x is length n)
    virtual double operator()(int n, const double *x) = 0;
  };

  // (Note: it is inefficient to use std::vector<double> for the arguments,
  //  since that would require a copy to be made of NLopt's double* data.)

  //////////////////////////////////////////////////////////////////////

  class opt {
  private:
    nlopt_opt o;
    
    void mythrow(nlopt_result ret) const {
      switch (ret) {
      case NLOPT_FAILURE: throw std::runtime_error("nlopt failure");
      case NLOPT_OUT_OF_MEMORY: throw std::bad_alloc();
      case NLOPT_INVALID_ARGS: throw std::invalid_argument("nlopt");
      case NLOPT_ROUNDOFF_LIMITED: throw std::runtime_error("nlopt roundoff");
      default: break;
      }
    }

    // nlopt_func wrapper around C++ "functional"
    static double myfunc(int n, const double *x, double *grad, void *f_) {
      func *f = reinterpret_cast<func*>(f_);
      return grad ? (*f)(n, x, grad) : (*f)(n, x);
    }

  public:
    // Constructors etc.
    opt() : o(NULL) {}
    opt(nlopt_algorithm a, int n) : o(nlopt_create(a, n)) {
      if (!o) throw std::bad_alloc();
    }
    opt(algorithm a, int n) : o(nlopt_create(nlopt_algorithm(a), n)) {
      if (!o) throw std::bad_alloc();
    }
    opt(const nlopt_opt o0) : o(nlopt_copy(o0)) {
      if (o0 && !o) throw std::bad_alloc();
    }
    ~opt() { nlopt_destroy(o); }
    opt(const opt& from) : o(nlopt_copy(from.o)) {
      if (from.o && !o) throw std::bad_alloc();
    }
    opt& operator=(opt const& f) {
      if (this == &f) return *this; // self-assignment
      nlopt_destroy(o);
      o = nlopt_copy(f.o);
      if (f.o && !o) throw std::bad_alloc();
      return *this;
    }

    // Do the optimization:
    result optimize(double *x, double &opt_f) {
      nlopt_result ret = nlopt_optimize(o, x, &opt_f);
      mythrow(ret);
      return result(ret);
    }
    result optimize(std::vector<double> &x, double &opt_f) {
      return optimize(x.empty() ? NULL : &x[0], opt_f);
    }

    // accessors:
    algorithm get_algorithm() const {
      if (!o) throw std::invalid_argument("uninitialized nlopt::opt");
      return algorithm(nlopt_get_algorithm(o));
    }
    int get_dimension() const {
      if (!o) throw std::invalid_argument("uninitialized nlopt::opt");
      return nlopt_get_dimension(o);
    }

    // Set the objective function
    void set_min_objective(nlopt_func f, void *f_data) {
      nlopt_result ret = nlopt_set_min_objective(o, f, f_data);
      mythrow(ret);
    }
    void set_min_objective(func *f) {
      set_min_objective(myfunc, f);
    }

    // Nonlinear constraints:

    void remove_inequality_constraints(void) {
      nlopt_result ret = nlopt_remove_inequality_constraints(o);
      mythrow(ret);
    }
    void add_inequality_constraint(nlopt_func f, void *f_data, double tol=0) {
      nlopt_result ret = nlopt_add_inequality_constraint(o, f, f_data, tol);
      mythrow(ret);
    }
    void add_inequality_constraint(func *f, double tol=0) {
      add_inequality_constraint(myfunc, f, tol);
    }

    void remove_equality_constraints(void) {
      nlopt_result ret = nlopt_remove_equality_constraints(o);
      mythrow(ret);
    }
    void add_equality_constraint(nlopt_func f, void *f_data, double tol=0) {
      nlopt_result ret = nlopt_add_equality_constraint(o, f, f_data, tol);
      mythrow(ret);
    }
    void add_equality_constraint(func *f, double tol=0) {
      add_equality_constraint(myfunc, f, tol);
    }

#define NLOPT_GETSET_VEC(name)						\
    void get_##name(double *v) const {					\
      nlopt_result ret = nlopt_get_##name(o, v);			\
      mythrow(ret);							\
    }									\
    void set_##name(const double *v) {					\
      nlopt_result ret = nlopt_set_##name(o, v);			\
      mythrow(ret);							\
    }									\
    void set_##name(double val) {					\
      nlopt_result ret = nlopt_set_##name##1(o, val);			\
      mythrow(ret);							\
    }									\
    void get_##name(std::vector<double> &v) const {			\
      if (o && unsigned(nlopt_get_dimension(o)) != v.size())		\
        throw std::invalid_argument("dimension mismatch");		\
      get_##name(v.empty() ? NULL : &v[0]);				\
    }									\
    std::vector<double> get_##name(void) const {			\
      if (!o) throw std::invalid_argument("uninitialized nlopt::opt");	\
      std::vector<double> v(unsigned(nlopt_get_dimension(o)));		\
      get_##name(v);							\
      return v;								\
    }			 						\
    void set_##name(const std::vector<double> &v) {			\
      if (o && unsigned(nlopt_get_dimension(o)) != v.size())		\
        throw std::invalid_argument("dimension mismatch");		\
      set_##name(v.empty() ? NULL : &v[0]);				\
    }

    NLOPT_GETSET_VEC(lower_bounds)
    NLOPT_GETSET_VEC(upper_bounds)

    // stopping criteria:

#define NLOPT_GETSET(T, name)						\
    T get_##name() const {						\
      if (!o) throw std::invalid_argument("uninitialized nlopt::opt");	\
      return nlopt_get_##name(o);					\
    }									\
    void set_##name(T name) {						\
      nlopt_result ret = nlopt_set_##name(o, name);			\
      mythrow(ret);							\
    }
    NLOPT_GETSET(double, stopval)
    NLOPT_GETSET(double, ftol_rel)
    NLOPT_GETSET(double, ftol_abs)
    NLOPT_GETSET(double, xtol_rel)
    NLOPT_GETSET_VEC(xtol_abs)
    NLOPT_GETSET(int, maxeval)
    NLOPT_GETSET(double, maxtime)

    // algorithm-specific parameters:

    void set_local_optimizer(const nlopt_opt lo) {
      nlopt_result ret = nlopt_set_local_optimizer(o, lo);
      mythrow(ret);
    }
    void set_local_optimizer(const opt &lo) {
      set_local_optimizer(lo.o);
    }

    NLOPT_GETSET(int, population)
    NLOPT_GETSET_VEC(initial_step)

    void set_default_initial_step(const double *x) {
      nlopt_result ret = nlopt_set_default_initial_step(o, x);
      mythrow(ret);
    }
    void set_default_initial_step(const std::vector<double> &x) {
      set_default_initial_step(x.empty() ? NULL : &x[0]);
    }
    void get_initial_step(const double *x, double *dx) const {
      nlopt_result ret = nlopt_get_initial_step(o, x, dx);
      mythrow(ret);
    }
    void get_initial_step(const std::vector<double> &x, std::vector<double> &dx) const {
      if (o && (unsigned(nlopt_get_dimension(o)) != x.size()
		|| unsigned(nlopt_get_dimension(o)) != dx.size()))
        throw std::invalid_argument("dimension mismatch");
      get_initial_step(x.empty() ? NULL : &x[0],
		       dx.empty() ? NULL : &dx[0]);
    }
    std::vector<double> get_initial_step(const std::vector<double> &x) const {
      if (!o) throw std::invalid_argument("uninitialized nlopt::opt");
      std::vector<double> v(unsigned(nlopt_get_dimension(o)));
      get_initial_step(x, v);
      return v;
    }
  };

  //////////////////////////////////////////////////////////////////////

} // namespace nlopt
