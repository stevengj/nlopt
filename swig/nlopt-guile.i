// -*- C++ -*-

%{
// vfunc wrapper around Guile function (val . grad) = f(x)
static double vfunc_guile(const std::vector<double> &x,
                          std::vector<double> &grad, void *f) {
  SCM xscm = scm_c_make_vector(x.size(), SCM_UNSPECIFIED);
  for (unsigned i = 0; i < x.size(); ++i)
    scm_c_vector_set_x(xscm, i, scm_make_real(x[i]));
  SCM ret = scm_call_1((SCM) f, xscm);
  if (scm_real_p(ret)) {
    if (grad.size()) throw std::invalid_argument("missing gradient");
    return scm_to_double(ret);
  }
  else if (scm_is_pair(ret)) { /* must be (cons value gradient) */
    SCM valscm = SCM_CAR(ret), grad_scm = grad_scm;
    if (grad.size() > 0
	&& scm_is_vector(grad_scm)
	&& scm_c_vector_length(grad_scm) == grad.size()) {
      for (unsigned i = 0; i < grad.size(); ++i)
	grad[i] = scm_to_double(scm_c_vector_ref(grad_scm, i));
    }
    else throw std::invalid_argument("invalid gradient");
    if (scm_real_p(valscm))
      return scm_to_double(valscm);
  }
  throw std::invalid_argument("invalid result passed to nlopt");
}

extern SCM nlopt_do_optimize(nlopt::opt &opt, const std::vector<double> &x);

// wrapper around opt::optimize that returns a triplet (result . (x . optf))
SCM nlopt_do_optimize(nlopt::opt &opt, const std::vector<double> &x0)
{
  std::vector<double> x(x0);
  double optf;
  nlopt::result res = opt.optimize(x, optf);
  SCM xscm = scm_c_make_vector(x.size(), SCM_UNSPECIFIED);
  for (unsigned i = 0; i < x.size(); ++i)
    scm_c_vector_set_x(xscm, i, scm_make_real(x[i]));
  return scm_cons(scm_from_int(int(res)),
		  scm_cons(xscm, scm_make_real(optf)));
}

%}

extern SCM nlopt_do_optimize(nlopt::opt &opt, const std::vector<double> &x);

%typemap(in)(nlopt::vfunc vf, void *f_data) {
  $1 = vfunc_guile;
  $2 = (void*) $input; // input is SCM pointer to Scheme function
}
%typecheck(SWIG_TYPECHECK_POINTER)(nlopt::vfunc vf, void *f_data) {
  $1 = SCM_NFALSEP(scm_procedure_p($input));
}

// export constants as variables, rather than as functions returning the value
%feature("constasvar", "1");

%scheme %{ 
(load-extension "libnlopt@NLOPT_SUFFIX@_guile.so" "SWIG_init")
%}
