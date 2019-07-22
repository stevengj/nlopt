#define NLOPT_OPTIMIZE_USAGE \
"Usage: [xopt, fopt, retcode] = nlopt_optimize(opt, xinit)\n" \
"\n" \
"Optimizes (minimizes or maximizes) a nonlinear function under\n" \
"nonlinear constraints from the starting guess xinit, where the\n" \
"objective, constraints, stopping criteria, and other options are \n" \
"specified in the structure opt described below.  A variety of local\n" \
"and global optimization algorithms can be used, as specified by the \n" \
"opt.algorithm parameter described below.  Returns the optimum\n" \
"function value fopt, the location xopt of the optimum, and a\n" \
"return code retcode described below (> 0 on success).\n" \
"\n" \
"The dimension (n) of the problem, i.e. the number of design variables,\n" \
"is specified implicitly via the length of xinit.\n" \
"\n" \
"This function is a front-end for the external routine nlopt_optimize\n" \
"in the free NLopt nonlinear-optimization library, which is a wrapper\n" \
"around a number of free/open-source optimization subroutines.  More\n" \
"details can be found on the NLopt web page (ab-initio.mit.edu/nlopt)\n" \
"and also under 'man nlopt_minimize' on Unix.\n" \
"\n" \
"OBJECTIVE FUNCTION:\n" \
"\n" \
"The objective function f is specified via opt.min_objective or\n" \
"opt.max_objective for minimization or maximization, respectively.\n" \
"opt.min/max_objective should be a handle (@) to a function of the form:\n" \
"\n" \
"   [val, gradient] = f(x)\n" \
"\n" \
"where x is a row vector, val is the function value f(x), and gradient\n" \
"is a row vector giving the gradient of the function with respect to x.\n" \
"The gradient is only used for gradient-based optimization algorithms;\n" \
"some of the algorithms (below) are derivative-free and only require\n" \
"f to return val (its value).\n" \
"\n" \
"BOUND CONSTRAINTS:\n" \
"\n" \
"Lower and/or upper bounds for the design variables x are specified\n" \
"via opt.lower_bounds and/or opt.upper_bounds, respectively: these\n" \
"are vectors (of the same length as xinit, above) giving the bounds\n" \
"in each component. An unbounded component may be specified by a\n" \
"lower/upper bound of -inf/+inf, respectively.  If opt.lower_bounds\n" \
"and/or opt.upper_bounds are not specified, the default bounds are\n" \
"-inf/+inf (i.e. unbounded), respectively.\n" \
"\n" \
"NONLINEAR CONSTRAINTS:\n" \
"\n" \
"Several of the algorithms in NLopt (MMA, COBYLA, and ORIG_DIRECT) also\n" \
"support arbitrary nonlinear inequality constraints, and some also allow\n" \
"nonlinear equality constraints (ISRES and AUGLAG). For these \n" \
"algorithms, you can specify as many nonlinear constraints as you wish.\n" \
"(The default is no nonlinear constraints.)\n" \
"\n" \
"Inequality constraints of the form fc{i}(x) <= 0 are specified via opt.fc,\n" \
"which is a cell array of function handles (@) of the same form as\n" \
"the objective function above (i.e., returning the value and optionally\n" \
"the gradient of the constraint function fc{i}, where the gradient is\n" \
"only needed for gradient-based algorithms).\n" \
"\n" \
"Equality constraints of the form h{i}(x) = 0 are specified via opt.h,\n" \
"which is a cell array of function handles (@) of the same form as\n" \
"the objective function above (i.e., returning the value and optionally\n" \
"the gradient of the constraint function h{i}, where the gradient is\n" \
"only needed for gradient-based algorithms).\n" \
"\n" \
"For both inequality and equality constraints, you can supply a\n" \
"\"tolerance\" for each constraint: this tolerance is used for convergence\n" \
"tests only, and a point x is considered feasible for purposes of\n" \
"convergence if the constraint is violated by the given tolerance.\n" \
"The tolerances are specified via opt.fc_tol and opt.h_tol, respectively,\n" \
"which must be vectors of the same length as opt.fc and opt.h, so\n" \
"that opt.fc_tol(i) is the tolerance for opt.fc{i} and opt.h_tol(i)\n" \
"is the tolerance for opt.h{i}.  These tolerances default to zero; a\n" \
"small nonzero tolerance is recommended, however, especially for h_tol.\n" \
"\n" \
"ALGORITHMS\n" \
"\n" \
"The optimization algorithm must be specified via opt.algorithm.\n" \
"\n" \
"The algorithm should be one of the following constants (name and\n" \
"interpretation are the same as for the C language interface).  Names\n" \
"with _G*_ are global optimization, and names with _L*_ are local\n" \
"optimization.  Names with _*N_ are derivative-free, while names\n" \
"with _*D_ are gradient-based algorithms.  Algorithms:\n" \
"\n" \
"NLOPT_GD_MLSL_LDS, NLOPT_GD_MLSL, NLOPT_GD_STOGO, NLOPT_GD_STOGO_RAND, \n" \
"NLOPT_GN_CRS2_LM, NLOPT_GN_DIRECT_L, NLOPT_GN_DIRECT_L_NOSCAL, \n" \
"NLOPT_GN_DIRECT_L_RAND, NLOPT_GN_DIRECT_L_RAND_NOSCAL, NLOPT_GN_DIRECT, \n" \
"NLOPT_GN_DIRECT_NOSCAL, NLOPT_GN_ISRES, NLOPT_GN_MLSL_LDS, NLOPT_GN_MLSL, \n" \
"NLOPT_GN_ORIG_DIRECT_L, NLOPT_GN_ORIG_DIRECT, NLOPT_LD_AUGLAG_EQ, \n" \
"NLOPT_LD_AUGLAG, NLOPT_LD_LBFGS, NLOPT_LD_LBFGS_NOCEDAL, NLOPT_LD_MMA, \n" \
"NLOPT_LD_TNEWTON, NLOPT_LD_TNEWTON_PRECOND, \n" \
"NLOPT_LD_TNEWTON_PRECOND_RESTART, NLOPT_LD_TNEWTON_RESTART, \n" \
"NLOPT_LD_VAR1, NLOPT_LD_VAR2, NLOPT_LN_AUGLAG_EQ, NLOPT_LN_AUGLAG, \n" \
"NLOPT_LN_BOBYQA, NLOPT_LN_COBYLA, NLOPT_LN_NELDERMEAD, \n" \
"NLOPT_LN_NEWUOA_BOUND, NLOPT_LN_NEWUOA, NLOPT_LN_PRAXIS, NLOPT_LN_SBPLX\n" \
"\n" \
"For more information on individual algorithms, see their individual\n" \
"help pages (e.g. \"help NLOPT_LN_SBPLX\").\n" \
"\n" \
"STOPPING CRITERIA:\n" \
"\n" \
"Multiple stopping criteria can be specified by setting one or more of\n" \
"the following fields of opt.  The optimization halts whenever any\n" \
"one of the given criteria is satisfied.\n" \
"\n" \
"opt.stopval: Stop when an objective value of at least stopval is found.\n" \
"   That is, stop minimizing when a value <= stopval is found, or stop\n" \
"   maximizing when a value >= stopval is found.\n" \
"\n" \
"opt.ftol_rel: Relative tolerance on function value, to stop when\n" \
"   an optimization step (or an estimate of the optimum) changes\n" \
"   the function value by less than opt.ftol_rel multiplied by\n" \
"   the absolute value of the function.\n" \
"\n" \
"opt.ftol_abs: Absolute tolerance on function value, to stop when\n" \
"   an optimization step (or an estimate of the optimum) changes\n" \
"   the function value by less than opt.ftol_abs.\n" \
"\n" \
"opt.xtol_rel: Relative tolerance on function value, to stop when\n" \
"   an optimization step (or an estimate of the optimum) changes\n" \
"   every component of x by less than opt.xtol_rel multiplied by\n" \
"   the absolute value of that component of x.\n" \
"\n" \
"opt.xtol_abs: Absolute tolerance on function value, to stop when\n" \
"   an optimization step (or an estimate of the optimum) changes\n" \
"   every component of x by less than that component of opt.xtol_abs\n" \
"   -- should be a vector of same length as x.\n" \
"\n" \
"opt.maxeval: Maximum number of function evaluations.\n" \
"\n" \
"opt.maxtime: Maximum runtime (in seconds) for the optimization.\n" \
"\n" \
"RETURN CODE:\n" \
"\n" \
"The retcode result is positive upon successful completion, and\n" \
"negative for an error.  The specific values are:\n" \
"\n" \
"generic success code: +1\n" \
"     stopval reached: +2\n" \
"        ftol reached: +3\n" \
"        xtol reached: +4\n" \
"     maxeval reached: +5\n" \
"     maxtime reached: +6\n" \
"generic failure code: -1\n" \
"   invalid arguments: -2\n" \
"       out of memory: -3\n" \
"    roundoff-limited: -4\n" \
"\n" \
"LOCAL OPTIMIZER:\n" \
"\n" \
"Some of the algorithms, especially MLSL and AUGLAG, use a different\n" \
"optimization algorithm as a subroutine, typically for local optimization.\n" \
"By default, they use MMA or COBYLA for gradient-based or derivative-free\n" \
"searching, respectively.  However, you can change this by specifying\n" \
"opt.local_optimizer: this is a structure with the same types of fields as opt\n" \
"(stopping criteria, algorithm, etcetera).  The objective function\n" \
"and nonlinear constraint parameters of opt.local_optimizer are ignored.\n" \
"\n" \
"INITIAL STEP SIZE:\n" \
"\n" \
"For derivative-free local-optimization algorithms, the optimizer must\n" \
"somehow decide on some initial step size to perturb x by when it begins\n" \
"the optimization. This step size should be big enough that the value\n" \
"of the objective changes significantly, but not too big if you want to\n" \
"find the local optimum nearest to x. By default, NLopt chooses this\n" \
"initial step size heuristically from the bounds, tolerances, and other\n" \
"information, but this may not always be the best choice.\n" \
"\n" \
"You can modify the initial step by setting opt.initial_step, which\n" \
"is a vector of the same length as x containing the (nonzero) initial\n" \
"step size for each component of x.\n" \
"\n" \
"STOCHASTIC POPULATION:\n" \
"\n" \
"Several of the stochastic search algorithms (e.g., CRS, MLSL, and\n" \
"ISRES) start by generating some initial \"population\" of random points\n" \
"x. By default, this initial population size is chosen heuristically in\n" \
"some algorithm-specific way, but the initial population can by changed\n" \
"by setting opt.population to the desired initial population size.\n" \
"\n" \
"VERBOSE OUTPUT:\n" \
"\n" \
"If opt.verbose is set to a nonzero value, then nlopt_optimize\n" \
"will print out verbose output; for example, it will print the\n" \
"value of the objective function after each evaluation.\n" \
"\n" \
"MORE INFORMATION:\n" \
"\n" \
"For more documentation, such as a detailed description of all the\n" \
"algorithms, see the NLopt home page: http://ab-initio.mit.edu/nlopt\n" \

