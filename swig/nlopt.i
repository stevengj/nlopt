// -*- C++ -*-
%module nlopt
%{
#include "nlopt.hpp"
%}

%include "std_vector.i"
namespace std {
  %template(nlopt_doublevector) vector<double>;
};

%rename(nlopt_opt) nlopt::opt;
%rename(nlopt_roundoff_limited) nlopt::roundoff_limited;
%rename(nlopt_forced_stop) nlopt::forced_stop;
%rename(nlopt_srand) nlopt::srand;
%rename(nlopt_srand_time) nlopt::srand_time;
%rename(nlopt_version) nlopt::version;
%rename(nlopt_algorithm_name) nlopt::algorithm_name;
%include "nlopt-enum-renames.i"

%include "nlopt-exceptions.i"

#ifdef SWIGGUILE
%include "nlopt-guile.i"
#endif

%include "nlopt.hpp"
