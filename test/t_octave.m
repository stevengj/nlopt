function t_octave()
    opt = struct();
    opt.algorithm = NLOPT_LD_MMA;
    %opt.algorithm = NLOPT_LN_COBYLA;
    opt.lower_bounds = [-inf, 0];
    opt.min_objective = @myfunc;
    opt.fc = { (@(x) myconstraint(x,2,0)), (@(x) myconstraint(x,-1,1)) };
    opt.fc_tol = [1e-8, 1e-8];
    opt.xtol_rel = 1e-4;
    display(opt)
    [xopt, fmin, retcode] = nlopt_optimize(opt, [1.234 5.678])
end

function [val, grad] = myfunc(x)
    val = sqrt(x(2));
    if (nargout > 1)
        grad = [0, 0.5 / val];
    end
end

function [val, grad] = myconstraint(x,a,b)
    val = (a*x(1) + b)^3 - x(2);
    if (nargout > 1)
        grad = [3*a*(a*x(1) + b)^2, -1];
    end
end
