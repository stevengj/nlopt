#include <octave/oct.h>
#include <octave/oct-map.h>
#include <octave/ov.h>
#include <math.h>
#include <stdio.h>

#include "nlopt.h"
#include "nlopt_minimize_usage.h"

static double struct_val_default(Octave_map &m, const std::string& k,
				 double dflt)
{
  if (m.contains(k)) {
    if (m.contents(k).length() == 1 && (m.contents(k))(0).is_real_scalar())
      return (m.contents(k))(0).double_value();
  }
  return dflt;
}

static Matrix struct_val_default(Octave_map &m, const std::string& k,
				 Matrix &dflt)
{
  if (m.contains(k)) {
    if ((m.contents(k)).length() == 1) {
      if ((m.contents(k))(0).is_real_scalar())
	return Matrix(1, dflt.length(), (m.contents(k))(0).double_value());
      else if ((m.contents(k))(0).is_real_matrix())
	return (m.contents(k))(0).matrix_value();
    }
  }
  return dflt;
}

typedef struct {
  octave_function *f;
  Cell f_data;
} user_function_data;

static double user_function(int n, const double *x,
			    double *gradient, /* NULL if not needed */
			    void *data_)
{
  user_function_data *data = (user_function_data *) data_;
  octave_value_list args(1 + data->f_data.length(), 0);
  Matrix xm(1,n);
  for (int i = 0; i < n; ++i)
    xm(i) = x[i];
  args(0) = xm;
  for (int i = 0; i < data->f_data.length(); ++i)
    args(1 + i) = data->f_data(i);
  octave_value_list res = data->f->do_multi_index_op(gradient ? 2 : 1, args); 
  if (res.length() < (gradient ? 2 : 1))
    gripe_user_supplied_eval("nlopt_minimize");
  else if (!res(0).is_real_scalar()
	   || (gradient && !res(1).is_real_matrix()
	       && !(n == 1 && res(1).is_real_scalar())))
    gripe_user_returned_invalid("nlopt_minimize");
  else {
    if (gradient) {
      if (n == 1 && res(1).is_real_scalar())
	gradient[0] = res(1).double_value();
      else {
	Matrix grad = res(1).matrix_value();
	for (int i = 0; i < n; ++i)
	  gradient[i] = grad(i);
      }
    }
    return res(0).double_value();
  }
  return 0;
}				 

#define CHECK(cond, msg) if (!(cond)) { fprintf(stderr, msg "\n\n"); print_usage("nlopt_minimize"); return retval; }

DEFUN_DLD(nlopt_minimize, args, nargout, NLOPT_MINIMIZE_USAGE)
{
  octave_value_list retval;
  double A;

  CHECK(args.length() == 7 && nargout <= 3, "wrong number of args");

  CHECK(args(0).is_real_scalar(), "n must be real scalar");
  nlopt_algorithm algorithm = nlopt_algorithm(args(0).int_value());

  user_function_data d;
  CHECK(args(1).is_function() || args(1).is_function_handle(), 
	"f must be function");
  d.f = args(1).function_value();
  CHECK(args(2).is_cell(), "f_data must be cell array");
  d.f_data = args(2).cell_value();

  CHECK(args(3).is_real_matrix() || args(3).is_real_scalar(),
	"lb must be real vector");
  Matrix lb = args(3).is_real_scalar() ?
    Matrix(1, 1, args(3).double_value()) : args(3).matrix_value();
  int n = lb.length();
  
  CHECK(args(4).is_real_matrix() || args(4).is_real_scalar(),
	"ub must be real vector");
  Matrix ub = args(4).is_real_scalar() ?
    Matrix(1, 1, args(4).double_value()) : args(4).matrix_value();
  CHECK(n == ub.length(), "lb and ub must have same length");

  CHECK(args(5).is_real_matrix() || args(5).is_real_scalar(),
	"x must be real vector");
  Matrix x = args(5).is_real_scalar() ?
    Matrix(1, 1, args(5).double_value()) : args(5).matrix_value();
  CHECK(n == x.length(), "x and lb/ub must have same length");

  CHECK(args(6).is_map(), "stop must be structure");
  Octave_map stop = args(6).map_value();
  double fmin_max = struct_val_default(stop, "fmin_max", -HUGE_VAL);
  double ftol_rel = struct_val_default(stop, "ftol_rel", 0);
  double ftol_abs = struct_val_default(stop, "ftol_abs", 0);
  double xtol_rel = struct_val_default(stop, "xtol_rel", 0);
  Matrix zeros(1, n, 0.0);
  Matrix xtol_abs = struct_val_default(stop, "xtol_abs", zeros);
  CHECK(n == xtol_abs.length(), "stop.xtol_abs must have same length as x");
  int maxeval = int(struct_val_default(stop, "maxeval", -1));
  double maxtime = struct_val_default(stop, "maxtime", -1);
  
  double fmin = HUGE_VAL;
  nlopt_result ret = nlopt_minimize(algorithm,
				    n,
				    user_function, &d,
				    lb.data(), ub.data(),
				    x.fortran_vec(), &fmin,
				    fmin_max, ftol_rel, ftol_abs,
				    xtol_rel, xtol_abs.data(),
				    maxeval, maxtime);
				    
  retval(0) = x;
  if (nargout > 1)
    retval(1) = fmin;
  if (nargout > 2)
    retval(2) = int(ret);

  return retval;
}
