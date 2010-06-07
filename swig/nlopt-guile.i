// -*- C++ -*-

%{
// because our f_data pointer to the Scheme function is stored on the
// heap, rather than the stack, it may be missed by the Guile garbage
// collection and be accidentally freed.  Hence, use NLopts munge
// feature to prevent this, by incrementing Guile's reference count.
static void *free_guilefunc(void *p) { 
  scm_gc_unprotect_object((SCM) p); return p; }
static void *dup_guilefunc(void *p) { 
  scm_gc_protect_object((SCM) p); return p; }

// vfunc wrapper around Guile function (val . grad) = f(x)
static double func_guile(unsigned n, const double *x, double *grad, void *f) {
  SCM xscm = scm_c_make_vector(n, SCM_UNSPECIFIED);
  for (unsigned i = 0; i < n; ++i)
    scm_c_vector_set_x(xscm, i, scm_make_real(x[i]));
  SCM ret = scm_call_1((SCM) f, xscm);
  if (scm_real_p(ret)) {
    if (grad) throw std::invalid_argument("missing gradient");
    return scm_to_double(ret);
  }
  else if (scm_is_pair(ret)) { /* must be (cons value gradient) */
    SCM valscm = SCM_CAR(ret), grad_scm = grad_scm;
    if (grad
	&& scm_is_vector(grad_scm)
	&& scm_c_vector_length(grad_scm) == n) {
      for (unsigned i = 0; i < n; ++i)
	grad[i] = scm_to_double(scm_c_vector_ref(grad_scm, i));
    }
    else throw std::invalid_argument("invalid gradient");
    if (scm_real_p(valscm))
      return scm_to_double(valscm);
  }
  throw std::invalid_argument("invalid result passed to nlopt");
}
%}

%typemap(in)(nlopt::vfunc vf, void *f_data, nlopt_munge md, nlopt_munge mc) {
  $1 = vfunc_guile;
  $2 = dup_guilefunc((void*) $input); // input = SCM pointer to Scheme function
  $3 = free_guilefunc;
  $4 = dup_guilefunc;
}
%typecheck(SWIG_TYPECHECK_POINTER)(nlopt::vfunc vf, void *f_data, nlopt_munge md, nlopt_munge mc) {
  $1 = SCM_NFALSEP(scm_procedure_p($input));
}

// export constants as variables, rather than as functions returning the value
%feature("constasvar", "1");

%scheme %{ 
(load-extension "libnlopt@NLOPT_SUFFIX@_guile.so" "SWIG_init")
%}
