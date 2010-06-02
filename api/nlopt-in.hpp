/* Copyright (c) 2007-2010 Massachusetts Institute of Technology
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
#include <cstdlib>
#include <cstring>
#include <cmath>

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

  typedef nlopt_func func; // nlopt::func synoynm

  // alternative to nlopt_func that takes std::vector<double>
  // ... unfortunately requires a data copy
  typedef double (*vfunc)(const std::vector<double> &x,
			  std::vector<double> &grad, void *data);

  //////////////////////////////////////////////////////////////////////
  
  // NLopt-specific exceptions (corresponding to error codes):
  class roundoff_limited : public std::runtime_error {
  public:
    roundoff_limited() : std::runtime_error("nlopt roundoff-limited") {}
  };

  class forced_stop : public std::runtime_error {
  public:
    forced_stop() : std::runtime_error("nlopt forced stop") {}
  };

  //////////////////////////////////////////////////////////////////////

  class opt {
  private:
    nlopt_opt o;
    
    void mythrow(nlopt_result ret) const {
      switch (ret) {
      case NLOPT_FAILURE: throw std::runtime_error("nlopt failure");
      case NLOPT_OUT_OF_MEMORY: throw std::bad_alloc();
      case NLOPT_INVALID_ARGS: throw std::invalid_argument("nlopt invalid argument");
      case NLOPT_ROUNDOFF_LIMITED: throw roundoff_limited();
      case NLOPT_FORCED_STOP: throw forced_stop();
      default: break;
      }
    }

    typedef struct {
      opt *o;
      func f; void *f_data;
      vfunc vf;
    } myfunc_data;

    // nlopt_func wrapper that catches exceptions
    static double myfunc(unsigned n, const double *x, double *grad, void *d_) {
      myfunc_data *d = reinterpret_cast<myfunc_data*>(d_);
      try {
	return d->f(n, x, grad, d->f_data);
      }
      catch (...) {
	d->o->force_stop(); // stop gracefully, opt::optimize will re-throw
	return HUGE_VAL;
      }
    }

    std::vector<double> xtmp, gradtmp, gradtmp0; // scratch for myvfunc

    // nlopt_func wrapper, using std::vector<double>
    static double myvfunc(unsigned n, const double *x, double *grad, void *d_){
      myfunc_data *d = reinterpret_cast<myfunc_data*>(d_);
      try {
	std::vector<double> &xv = d->o->xtmp;
	if (n) std::memcpy(&xv[0], x, n * sizeof(double));
	double val=d->vf(xv, grad ? d->o->gradtmp : d->o->gradtmp0, d->f_data);
	if (grad && n) {
	  std::vector<double> &gradv = d->o->gradtmp;
	  std::memcpy(grad, &gradv[0], n * sizeof(double));
	}
	return val;
      }
      catch (...) {
	d->o->force_stop(); // stop gracefully, opt::optimize will re-throw
	return HUGE_VAL;
      }
    }

    void alloc_tmp() {
      if (xtmp.size() != nlopt_get_dimension(o)) {
	xtmp = std::vector<double>(nlopt_get_dimension(o));
	gradtmp = std::vector<double>(nlopt_get_dimension(o));
      }
    }

  public:
    // Constructors etc.
    opt() : o(NULL), xtmp(0), gradtmp(0), gradtmp0(0) {}
    ~opt() { nlopt_destroy(o); }
    opt(algorithm a, unsigned n) : 
      o(nlopt_create(nlopt_algorithm(a), n)), 
      xtmp(0), gradtmp(0), gradtmp0(0) {
      if (!o) throw std::bad_alloc();
      nlopt_set_free_f_data(o, 1);
    }
    opt(const opt& from) : o(nlopt_copy(from.o)) {
      if (from.o && !o) throw std::bad_alloc();
      mythrow(nlopt_dup_f_data(o, sizeof(myfunc_data)));
    }
    opt& operator=(opt const& f) {
      if (this == &f) return *this; // self-assignment
      nlopt_destroy(o);
      o = nlopt_copy(f.o);
      if (f.o && !o) throw std::bad_alloc();
      mythrow(nlopt_dup_f_data(o, sizeof(myfunc_data)));
      return *this;
    }

    // Do the optimization:
    result optimize(std::vector<double> &x, double &opt_f) {
      if (o && nlopt_get_dimension(o) != x.size())
        throw std::invalid_argument("dimension mismatch");
      nlopt_result ret = nlopt_optimize(o, x.empty() ? NULL : &x[0], &opt_f);
      mythrow(ret);
      return result(ret);
    }

    // accessors:
    algorithm get_algorithm() const {
      if (!o) throw std::runtime_error("uninitialized nlopt::opt");
      return algorithm(nlopt_get_algorithm(o));
    }
    const char *get_algorithm_name() const {
      if (!o) throw std::runtime_error("uninitialized nlopt::opt");
      return nlopt_algorithm_name(nlopt_get_algorithm(o));
    }
    unsigned get_dimension() const {
      if (!o) throw std::runtime_error("uninitialized nlopt::opt");
      return nlopt_get_dimension(o);
    }

    // Set the objective function
    void set_min_objective(func f, void *f_data) {
      myfunc_data *d = (myfunc_data *) std::malloc(sizeof(myfunc_data));
      if (!d) throw std::bad_alloc();
      d->o = this; d->f = f; d->f_data = f_data; d->vf = NULL;
      mythrow(nlopt_set_min_objective(o, myfunc, d)); // d freed via o
    }
    void set_min_objective(vfunc vf, void *f_data) {
      myfunc_data *d = (myfunc_data *) std::malloc(sizeof(myfunc_data));
      if (!d) throw std::bad_alloc();
      d->o = this; d->f = NULL; d->f_data = f_data; d->vf = vf;
      mythrow(nlopt_set_min_objective(o, myvfunc, d)); // d freed via o
      alloc_tmp();
    }
    void set_max_objective(func f, void *f_data) {
      myfunc_data *d = (myfunc_data *) std::malloc(sizeof(myfunc_data));
      if (!d) throw std::bad_alloc();
      d->o = this; d->f = f; d->f_data = f_data; d->vf = NULL;
      mythrow(nlopt_set_max_objective(o, myfunc, d)); // d freed via o
    }
    void set_max_objective(vfunc vf, void *f_data) {
      myfunc_data *d = (myfunc_data *) std::malloc(sizeof(myfunc_data));
      if (!d) throw std::bad_alloc();
      d->o = this; d->f = NULL; d->f_data = f_data; d->vf = vf;
      mythrow(nlopt_set_max_objective(o, myvfunc, d)); // d freed via o
      alloc_tmp();
    }

    // Nonlinear constraints:

    void remove_inequality_constraints(void) {
      nlopt_result ret = nlopt_remove_inequality_constraints(o);
      mythrow(ret);
    }
    void add_inequality_constraint(func f, void *f_data, double tol=0) {
      myfunc_data *d = (myfunc_data *) std::malloc(sizeof(myfunc_data));
      if (!d) throw std::bad_alloc();
      d->o = this; d->f = f; d->f_data = f_data; d->vf = NULL;
      mythrow(nlopt_add_inequality_constraint(o, myfunc, d, tol));
    }
    void add_inequality_constraint(vfunc vf, void *f_data, double tol=0) {
      myfunc_data *d = (myfunc_data *) std::malloc(sizeof(myfunc_data));
      if (!d) throw std::bad_alloc();
      d->o = this; d->f = NULL; d->f_data = f_data; d->vf = vf;
      mythrow(nlopt_add_inequality_constraint(o, myvfunc, d, tol));
      alloc_tmp();
    }

    void remove_equality_constraints(void) {
      nlopt_result ret = nlopt_remove_equality_constraints(o);
      mythrow(ret);
    }
    void add_equality_constraint(func f, void *f_data, double tol=0) {
      myfunc_data *d = (myfunc_data *) std::malloc(sizeof(myfunc_data));
      if (!d) throw std::bad_alloc();
      d->o = this; d->f = f; d->f_data = f_data; d->vf = NULL;
      mythrow(nlopt_add_equality_constraint(o, myfunc, d, tol));
    }
    void add_equality_constraint(vfunc vf, void *f_data, double tol=0) {
      myfunc_data *d = (myfunc_data *) std::malloc(sizeof(myfunc_data));
      if (!d) throw std::bad_alloc();
      d->o = this; d->f = NULL; d->f_data = f_data; d->vf = vf;
      mythrow(nlopt_add_equality_constraint(o, myvfunc, d, tol));
      alloc_tmp();
    }

#define NLOPT_GETSET_VEC(name)						\
    void set_##name(double val) {					\
      nlopt_result ret = nlopt_set_##name##1(o, val);			\
      mythrow(ret);							\
    }									\
    void get_##name(std::vector<double> &v) const {			\
      if (o && nlopt_get_dimension(o) != v.size())			\
        throw std::invalid_argument("dimension mismatch");		\
      nlopt_result ret = nlopt_get_##name(o, v.empty() ? NULL : &v[0]);	\
      mythrow(ret);							\
    }									\
    std::vector<double> get_##name(void) const {			\
      if (!o) throw std::runtime_error("uninitialized nlopt::opt");	\
      std::vector<double> v(nlopt_get_dimension(o));			\
      get_##name(v);							\
      return v;								\
    }			 						\
    void set_##name(const std::vector<double> &v) {			\
      if (o && nlopt_get_dimension(o) != v.size())			\
        throw std::invalid_argument("dimension mismatch");		\
      nlopt_result ret = nlopt_set_##name(o, v.empty() ? NULL : &v[0]);	\
      mythrow(ret);							\
    }

    NLOPT_GETSET_VEC(lower_bounds)
    NLOPT_GETSET_VEC(upper_bounds)

    // stopping criteria:

#define NLOPT_GETSET(T, name)						\
    T get_##name() const {						\
      if (!o) throw std::runtime_error("uninitialized nlopt::opt");	\
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

    NLOPT_GETSET(int, force_stop)
    void force_stop() { set_force_stop(1); }

    // algorithm-specific parameters:

    void set_local_optimizer(const opt &lo) {
      nlopt_result ret = nlopt_set_local_optimizer(o, lo.o);
      mythrow(ret);
    }

    NLOPT_GETSET(unsigned, population)
    NLOPT_GETSET_VEC(initial_step)

    void set_default_initial_step(const std::vector<double> &x) {
      nlopt_result ret 
	= nlopt_set_default_initial_step(o, x.empty() ? NULL : &x[0]);
      mythrow(ret);
    }
    void get_initial_step(const std::vector<double> &x, std::vector<double> &dx) const {
      if (o && (nlopt_get_dimension(o) != x.size()
		|| nlopt_get_dimension(o) != dx.size()))
        throw std::invalid_argument("dimension mismatch");
      nlopt_result ret = nlopt_get_initial_step(o, x.empty() ? NULL : &x[0],
						dx.empty() ? NULL : &dx[0]);
      mythrow(ret);
    }
    std::vector<double> get_initial_step(const std::vector<double> &x) const {
      if (!o) throw std::runtime_error("uninitialized nlopt::opt");
      std::vector<double> v(nlopt_get_dimension(o));
      get_initial_step(x, v);
      return v;
    }
  };

#undef NLOPT_GETSET
#undef NLOPT_GETSET_VEC

  //////////////////////////////////////////////////////////////////////

  inline void srand(unsigned long seed) { nlopt_srand(seed); }
  inline void srand_time(void) { nlopt_srand_time(); }
  inline void version(int &major, int &minor, int &bugfix) {
    nlopt_version(&major, &minor, &bugfix);
  }
  inline const char *algorithm_name(algorithm a) {
    return nlopt_algorithm_name(nlopt_algorithm(a));
  }

  //////////////////////////////////////////////////////////////////////

} // namespace nlopt
