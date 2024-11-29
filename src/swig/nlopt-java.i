// -*- C++ -*-
// kate: hl c++

// use Java naming conventions
%rename("%(camelcase)s", match="enum") "";
%rename("%(camelcase)s", match="class") "";
%rename("%(lowercamelcase)s", %$isfunction) "";

// use proper Java enums
%include "enums.swg"
// use Java code for the constants in the enums instead of calling a C function
%javaconst(1);

// pointer-based API not supported, use version_{major,minor,bugfix} instead
%ignore version;
// pointer-based API not supported, use the other overload instead
%ignore optimize(std::vector<double> &, double &);
// unsupported function APIs, use the ones with nlopt_munge instead
%ignore set_min_objective(func, void *);
%ignore set_min_objective(vfunc, void *);
%ignore set_min_objective(functor_type);
%ignore set_max_objective(func, void *);
%ignore set_max_objective(vfunc, void *);
%ignore set_max_objective(functor_type);
%ignore add_inequality_constraint(func, void *);
%ignore add_inequality_constraint(func, void *, double);
%ignore add_inequality_constraint(vfunc, void *);
%ignore add_inequality_constraint(vfunc, void *, double);
%ignore add_inequality_mconstraint(mfunc, void *, const std::vector<double> &);
%ignore add_equality_constraint(func, void *);
%ignore add_equality_constraint(func, void *, double);
%ignore add_equality_constraint(vfunc, void *);
%ignore add_equality_constraint(vfunc, void *, double);
%ignore add_equality_mconstraint(mfunc, void *, const std::vector<double> &);

// Munge function types
%extend nlopt::opt {
  %proxycode {
    public static interface Func {
      public double apply(double[] x, double[] gradient);
    }

    public static interface MFunc {
      public double[] apply(double[] x, double[] gradient);
    }
  }
}

%{
struct jfunc {
  JNIEnv *jenv;
  jobject func;
  jmethodID method;
};

static void *free_jfunc(void *p) {
  ((jfunc *) p)->jenv->DeleteGlobalRef(((jfunc *) p)->func);
  delete (jfunc *) p;
  return (void *) 0;
}

static void *dup_jfunc(void *p) {
  jfunc *q = new jfunc;
  q->jenv = ((jfunc *) p)->jenv;
  q->func = q->jenv->NewGlobalRef(((jfunc *) p)->func);
  q->method = ((jfunc *) p)->method;
  return (void *) q;
}

static double func_java(unsigned n, const double *x, double *grad, void *f)
{
  JNIEnv *jenv = ((jfunc *) f)->jenv;
  jobject func = ((jfunc *) f)->func;
  jmethodID method = ((jfunc *) f)->method;

  jdoubleArray jx = jenv->NewDoubleArray(n);
  if (!jx || jenv->ExceptionCheck()) {
    throw nlopt::forced_stop();
  }
  jenv->SetDoubleArrayRegion(jx, 0, n, x);
  jdoubleArray jgrad = (jdoubleArray) 0;
  if (grad) {
    jgrad = jenv->NewDoubleArray(n);
    if (!jgrad || jenv->ExceptionCheck()) {
      jenv->DeleteLocalRef(jx);
      throw nlopt::forced_stop();
    }
    jenv->SetDoubleArrayRegion(jgrad, 0, n, grad);
  }

  jdouble res = jenv->CallDoubleMethod(func, method, jx, jgrad);
  jenv->DeleteLocalRef(jx);

  if (jenv->ExceptionCheck()) {
    if (jgrad) {
      jenv->DeleteLocalRef(jgrad);
    }
    throw nlopt::forced_stop();
  }

  if (grad) {
    jenv->GetDoubleArrayRegion(jgrad, 0, n, grad);
    jenv->DeleteLocalRef(jgrad);
  }

  return res;
}

static void mfunc_java(unsigned m, double *result,
			 unsigned n, const double *x, double *grad, void *f)
{
  JNIEnv *jenv = ((jfunc *) f)->jenv;
  jobject func = ((jfunc *) f)->func;
  jmethodID method = ((jfunc *) f)->method;

  jdoubleArray jx = jenv->NewDoubleArray(n);
  if (!jx || jenv->ExceptionCheck()) {
    throw nlopt::forced_stop();
  }
  jenv->SetDoubleArrayRegion(jx, 0, n, x);
  jdoubleArray jgrad = (jdoubleArray) 0;
  if (grad) {
    jgrad = jenv->NewDoubleArray(m * n);
    if (!jgrad || jenv->ExceptionCheck()) {
      jenv->DeleteLocalRef(jx);
      throw nlopt::forced_stop();
    }
    jenv->SetDoubleArrayRegion(jgrad, 0, m * n, grad);
  }

  jdoubleArray res = (jdoubleArray) jenv->CallObjectMethod(func, method, jx, jgrad);
  jenv->DeleteLocalRef(jx);

  if (!res || jenv->ExceptionCheck()) {
    if (jgrad) {
      jenv->DeleteLocalRef(jgrad);
    }
    if (res) {
      jenv->DeleteLocalRef(res);
    }
    throw nlopt::forced_stop();
  }

  jenv->GetDoubleArrayRegion(res, 0, m, result);
  jenv->DeleteLocalRef(res);

  if (grad) {
    jenv->GetDoubleArrayRegion(jgrad, 0, m * n, grad);
    jenv->DeleteLocalRef(jgrad);
  }
}
%}

%typemap(jni)(nlopt::func f, void *f_data, nlopt_munge md, nlopt_munge mc) "jobject"
%typemap(jtype)(nlopt::func f, void *f_data, nlopt_munge md, nlopt_munge mc) "java.lang.Object"
%typemap(jstype)(nlopt::func f, void *f_data, nlopt_munge md, nlopt_munge mc) "Func"
%typemap(in)(nlopt::func f, void *f_data, nlopt_munge md, nlopt_munge mc) {
  $1 = func_java;
  jfunc jf = {jenv, $input, jenv->GetMethodID(jenv->FindClass("nlopt/Opt$Func"), "apply", "([D[D)D")};
  $2 = dup_jfunc((void *) &jf);
  $3 = free_jfunc;
  $4 = dup_jfunc;
}
%typemap(javain)(nlopt::func f, void *f_data, nlopt_munge md, nlopt_munge mc) "$javainput"

%typemap(jni)(nlopt::mfunc mf, void *f_data, nlopt_munge md, nlopt_munge mc) "jobject"
%typemap(jtype)(nlopt::mfunc mf, void *f_data, nlopt_munge md, nlopt_munge mc) "java.lang.Object"
%typemap(jstype)(nlopt::mfunc mf, void *f_data, nlopt_munge md, nlopt_munge mc) "MFunc"
%typemap(in)(nlopt::mfunc mf, void *f_data, nlopt_munge md, nlopt_munge mc) {
  $1 = mfunc_java;
  jfunc jf = {jenv, $input, jenv->GetMethodID(jenv->FindClass("nlopt/Opt$MFunc"), "apply", "([D[D)[D")};
  $2 = dup_jfunc((void *) &jf);
  $3 = free_jfunc;
  $4 = dup_jfunc;
}
%typemap(javain)(nlopt::mfunc mf, void *f_data, nlopt_munge md, nlopt_munge mc) "$javainput"

// Make exception classes Java-compliant
%rename(ForcedStopException) nlopt::forced_stop;
%typemap(javabase) nlopt::forced_stop "java.lang.RuntimeException"
%typemap(javabody) nlopt::forced_stop ""
%typemap(javadestruct) nlopt::forced_stop ""
%typemap(javafinalize) nlopt::forced_stop ""
%ignore nlopt::forced_stop::forced_stop;
%extend nlopt::forced_stop {
  %proxycode {
    public ForcedStopException(String message) {
      super(message);
    }
  }
}
%rename(RoundoffLimitedException) nlopt::roundoff_limited;
%typemap(javabase) nlopt::roundoff_limited "java.lang.RuntimeException"
%typemap(javabody) nlopt::roundoff_limited ""
%typemap(javadestruct) nlopt::roundoff_limited ""
%typemap(javafinalize) nlopt::roundoff_limited ""
%ignore nlopt::roundoff_limited::roundoff_limited;
%extend nlopt::roundoff_limited {
  %proxycode {
    public RoundoffLimitedException(String message) {
      super(message);
    }
  }
}

// Map exceptions
%typemap(throws) std::bad_alloc %{
  SWIG_JavaThrowException(jenv, SWIG_JavaOutOfMemoryError, $1.what());
  return $null;
%}

%typemap(throws) nlopt::forced_stop %{
  if (!jenv->ExceptionCheck()) {
    jclass excep = jenv->FindClass("nlopt/ForcedStopException");
    if (excep)
      jenv->ThrowNew(excep, $1.what());
  }
  return $null;
%}

%typemap(throws) nlopt::roundoff_limited %{
  if (!jenv->ExceptionCheck()) {
    jclass excep = jenv->FindClass("nlopt/RoundoffLimitedException");
    if (excep)
      jenv->ThrowNew(excep, $1.what());
  }
  return $null;
%}

