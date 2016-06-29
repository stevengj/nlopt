function val = nlopt_algorithm(name)
    %NLOPT_ALGORITHM  Get algorithm constant by name
    %
    % Algorithms supported:
    %
    % * NLOPT_GN_DIRECT: DIRECT (global, no-derivative)
    % * NLOPT_GN_DIRECT_L: DIRECT-L (global, no-derivative)
    % * NLOPT_GN_DIRECT_L_RAND: Randomized DIRECT-L (global, no-derivative)
    % * NLOPT_GN_DIRECT_NOSCAL: Unscaled DIRECT (global, no-derivative)
    % * NLOPT_GN_DIRECT_L_NOSCAL: Unscaled DIRECT-L (global, no-derivative)
    % * NLOPT_GN_DIRECT_L_RAND_NOSCAL: Unscaled Randomized DIRECT-L (global, no-derivative)
    % * NLOPT_GN_ORIG_DIRECT: Original DIRECT version (global, no-derivative)
    % * NLOPT_GN_ORIG_DIRECT_L: Original DIRECT-L version (global, no-derivative)
    % * NLOPT_GD_STOGO: StoGO (global, derivative-based)
    % * NLOPT_GD_STOGO_RAND: StoGO with randomized search (global, derivative-based)
    % * NLOPT_LD_LBFGS_NOCEDAL: original L-BFGS code by Nocedal et al. (NOT COMPILED)
    % * NLOPT_LD_LBFGS: Limited-memory BFGS (L-BFGS) (local, derivative-based)
    % * NLOPT_LN_PRAXIS: Principal-axis, praxis (local, no-derivative)
    % * NLOPT_LD_VAR1: Limited-memory variable-metric, rank 1 (local, derivative-based)
    % * NLOPT_LD_VAR2: Limited-memory variable-metric, rank 2 (local, derivative-based)
    % * NLOPT_LD_TNEWTON: Truncated Newton (local, derivative-based)
    % * NLOPT_LD_TNEWTON_RESTART: Truncated Newton with restarting (local, derivative-based)
    % * NLOPT_LD_TNEWTON_PRECOND: Preconditioned truncated Newton (local, derivative-based)
    % * NLOPT_LD_TNEWTON_PRECOND_RESTART: Preconditioned truncated Newton with restarting (local, derivative-based)
    % * NLOPT_GN_CRS2_LM: Controlled random search (CRS2) with local mutation (global, no-derivative)
    % * NLOPT_GN_MLSL: Multi-level single-linkage (MLSL), random (global, no-derivative)
    % * NLOPT_GD_MLSL: Multi-level single-linkage (MLSL), random (global, derivative)
    % * NLOPT_GN_MLSL_LDS: Multi-level single-linkage (MLSL), quasi-random (global, no-derivative)
    % * NLOPT_GD_MLSL_LDS: Multi-level single-linkage (MLSL), quasi-random (global, derivative)
    % * NLOPT_LD_MMA: Method of Moving Asymptotes (MMA) (local, derivative)
    % * NLOPT_LN_COBYLA: COBYLA (Constrained Optimization BY Linear Approximations) (local, no-derivative)
    % * NLOPT_LN_NEWUOA: NEWUOA unconstrained optimization via quadratic models (local, no-derivative)
    % * NLOPT_LN_NEWUOA_BOUND: Bound-constrained optimization via NEWUOA-based quadratic models (local, no-derivative)
    % * NLOPT_LN_NELDERMEAD: Nelder-Mead simplex algorithm (local, no-derivative)
    % * NLOPT_LN_SBPLX: Sbplx variant of Nelder-Mead (re-implementation of Rowan's Subplex) (local, no-derivative)
    % * NLOPT_LN_AUGLAG: Augmented Lagrangian method (local, no-derivative)
    % * NLOPT_LD_AUGLAG: Augmented Lagrangian method (local, derivative)
    % * NLOPT_LN_AUGLAG_EQ: Augmented Lagrangian method for equality constraints (local, no-derivative)
    % * NLOPT_LD_AUGLAG_EQ: Augmented Lagrangian method for equality constraints (local, derivative)
    % * NLOPT_LN_BOBYQA: BOBYQA bound-constrained optimization via quadratic models (local, no-derivative)
    % * NLOPT_GN_ISRES: ISRES evolutionary constrained optimization (global, no-derivative)
    % * NLOPT_AUGLAG: Augmented Lagrangian method (needs sub-algorithm)
    % * NLOPT_AUGLAG_EQ: Augmented Lagrangian method for equality constraints (needs sub-algorithm)
    % * NLOPT_G_MLSL: Multi-level single-linkage (MLSL), random (global, needs sub-algorithm)
    % * NLOPT_G_MLSL_LDS: Multi-level single-linkage (MLSL), quasi-random (global, needs sub-algorithm)
    % * NLOPT_LD_SLSQP: Sequential Quadratic Programming (SQP) (local, derivative)
    % * NLOPT_LD_CCSAQ: CCSA (Conservative Convex Separable Approximations) with simple quadratic approximations (local, derivative)
    % * NLOPT_GN_ESCH: ESCH evolutionary strategy
    %
    % See also: nlopt_optimize
    %

    % same order as enum definition in header file
    algorithms = {
        'NLOPT_GN_DIRECT'
        'NLOPT_GN_DIRECT_L'
        'NLOPT_GN_DIRECT_L_RAND'
        'NLOPT_GN_DIRECT_NOSCAL'
        'NLOPT_GN_DIRECT_L_NOSCAL'
        'NLOPT_GN_DIRECT_L_RAND_NOSCAL'
        'NLOPT_GN_ORIG_DIRECT'
        'NLOPT_GN_ORIG_DIRECT_L'
        'NLOPT_GD_STOGO'
        'NLOPT_GD_STOGO_RAND'
        'NLOPT_LD_LBFGS_NOCEDAL'
        'NLOPT_LD_LBFGS'
        'NLOPT_LN_PRAXIS'
        'NLOPT_LD_VAR1'
        'NLOPT_LD_VAR2'
        'NLOPT_LD_TNEWTON'
        'NLOPT_LD_TNEWTON_RESTART'
        'NLOPT_LD_TNEWTON_PRECOND'
        'NLOPT_LD_TNEWTON_PRECOND_RESTART'
        'NLOPT_GN_CRS2_LM'
        'NLOPT_GN_MLSL'
        'NLOPT_GD_MLSL'
        'NLOPT_GN_MLSL_LDS'
        'NLOPT_GD_MLSL_LDS'
        'NLOPT_LD_MMA'
        'NLOPT_LN_COBYLA'
        'NLOPT_LN_NEWUOA'
        'NLOPT_LN_NEWUOA_BOUND'
        'NLOPT_LN_NELDERMEAD'
        'NLOPT_LN_SBPLX'
        'NLOPT_LN_AUGLAG'
        'NLOPT_LD_AUGLAG'
        'NLOPT_LN_AUGLAG_EQ'
        'NLOPT_LD_AUGLAG_EQ'
        'NLOPT_LN_BOBYQA'
        'NLOPT_GN_ISRES'
        'NLOPT_AUGLAG'
        'NLOPT_AUGLAG_EQ'
        'NLOPT_G_MLSL'
        'NLOPT_G_MLSL_LDS'
        'NLOPT_LD_SLSQP'
        'NLOPT_LD_CCSAQ'
        'NLOPT_GN_ESCH'
    };
    if nargin == 0
        val = algorithms;
    else
        % fuzzy matching
        if ~strncmpi(name, 'NLOPT_', 6)
            name = ['NLOPT_' name];
        end
        name = validatestring(name, algorithms);
        % corresponding index (0-based)
        val = find(strcmpi(name, algorithms), 1, 'first') - 1;
    end
end
