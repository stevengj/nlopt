function [x, fval, exitflag, output] = nlopt_minimize(algorithm, f, f_data, lb, ub, x0, opt)
    %NLOPT_MINIMIZE  Minimizes a nonlinear multivariable function
    %
    % Minimizes a nonlinear multivariable function f(x, f_data{:}),
    % where x is a row vector, returning the optimal x found (xopt) along with
    % the minimum function value (fval = f(xopt)) and a return code (exitflag).
    % A variety of local and global optimization algorithms can be used,
    % as specified by the algorithm parameter described below.  lb and ub
    % are row vectors giving the upper and lower bounds on x, x0 is
    % a row vector giving the initial guess for x, and opt is a struct
    % containing termination conditions (see below).
    %
    % This function is a front-end for the external routine nlopt_minimize
    % in the free NLopt nonlinear-optimization library, which is a wrapper
    % around a number of free/open-source optimization subroutines.  More
    % details can be found on the NLopt web page (ab-initio.mit.edu/nlopt)
    % and also under 'man nlopt_minimize' on Unix.
    %
    % OBJECTIVE FUNCTION:
    %
    % f should be a handle (@) to a function of the form:
    %
    %    [val, gradient] = f(x, ...)
    %
    % where x is a row vector, val is the function value f(x), and gradient
    % is a row vector giving the gradient of the function with respect to x.
    % The gradient is only used for gradient-based optimization algorithms;
    % some of the algorithms (below) are derivative-free and only require
    % f to return val (its value).
    %
    % f can take additional arguments (...) which are passed via the argument
    % f_data: f_data is a cell array of the additional arguments to pass to f.
    % (Recall that cell arrays are specified by curly brackets { ... }.  For
    % example, pass f_data={} for functions that require no additional arguments.)
    %
    % BOUND CONSTRAINTS:
    %
    % Lower and/or upper bounds for the design variables x are specified
    % via lb and/or ub, respectively: these
    % are vectors (of the same length as x0, above) giving the bounds
    % in each component. An unbounded component may be specified by a
    % lower/upper bound of -inf/+inf, respectively.
    %
    %
    % ALGORITHMS:
    %
    % The algorithm should be one of the following constants (name and
    % interpretation are the same as for the C language interface).  Names
    % with G*_ are global optimization, and names with L*_ are local
    % optimization.  Names with *N_ are derivative-free, while names
    % with *D_ are gradient-based algorithms.  See nlopt_algorithms for the
    % full list.
    %
    % For more information on individual algorithms, see their individual
    % help pages.
    %
    % STOPPING CRITERIA:
    %
    % opt describes the termination criteria, and is a struct with a
    % number of optional fields:
    %
    %     opt.minf_max = stop when f < minf_max is found
    %     opt.ftol_rel = fractional tolerance on function value
    %     opt.ftol_abs = absolute tolerance on function value
    %     opt.xtol_rel = fractional tolerance on x
    %     opt.xtol_abs = row vector of absolute tolerances on x components
    %     opt.maxeval = maximum number of function evaluations
    %     opt.maxtime = maximum run time in seconds
    %     opt.verbose = > 0 indicates verbose output
    %
    % Minimization stops when any one of these conditions is met; any
    % condition that is omitted from opt will be ignored.  WARNING:
    % not all algorithms interpret the stopping criteria in exactly the
    % same way, and in any case ftol/xtol specify only a crude estimate
    % for the accuracy of the minimum function value/x.
    %
    % MORE INFORMATION:
    %
    % For more documentation, such as a detailed description of all the
    % algorithms, see the NLopt home page: http://ab-initio.mit.edu/nlopt
    %
    % See also: nlopt_optimize, nlopt_minimize_constrained, fminunc
    %

    if isfield(opt, 'minf_max')
        opt.stopval = opt.minf_max;
    end
    opt.algorithm = algorithm;
    opt.min_objective = @(x) f(x, f_data{:});
    opt.lower_bounds = lb;
    opt.upper_bounds = ub;
    [x, fval, exitflag, output] = nlopt_optimize(opt, x0);

end
