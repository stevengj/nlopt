import nlopt
import math as m

numFuncEval = [0]
minf_ext = [1e300]

lb = [-10., -10.]
ub = [+10., +10.]


def myfunc(x, grad):
    for j in range(2):
        if m.isnan(x[j]) or x[j] < lb[j] or x[j] > ub[j]:
            return 1.e10
    numFuncEval[0] += 1

    x1 = x[0]
    x2 = x[1]

    # http://al-roomi.org/benchmarks/unconstrained/2-dimensions/65-powell-s-badly-scaled-function
    # Powell's Badly Scaled Function
    # Range of initial points: -10 < xj < 10 , j=1,2
    # Global minima: (x1,x2)=(1.098...e-5, 9.106...)
    # f(x1,x2)=0
    f1 = (10000. * x1 * x2 - 1.)**2
    f2 = (m.exp(-x1) + m.exp(-x2) - 1.0001)**2
    retval = f1 + f2

    if grad.size > 0:
        # raise ValueError('Cannot suppply gradient values')
        grad[0] = 2.0 * (10000. * x1 * x2 - 1.) * 10000. * x2 + 2.0 * (m.exp(-x1) + m.exp(-x2) - 1.0001) * -1.0
        grad[1] = 2.0 * (10000. * x1 * x2 - 1.) * 10000. * x1 + 2.0 * (m.exp(-x1) + m.exp(-x2) - 1.0001) * -1.0

    # print("myfunc: x:", x, ", val:", retval)

    if retval < minf_ext[0]:
        minf_ext[0] = retval

    return retval


for algo in range(nlopt.NUM_ALGORITHMS):
    if algo in [nlopt.LD_LBFGS, nlopt.LD_VAR1, nlopt.LD_VAR2]:
        continue
    opt = nlopt.opt(algo, 2)

    print('-'*40)
    print("Algo:", opt.get_algorithm_name(), algo)
    numFuncEval[0] = 0
    minf_ext = [1e300]

    opt.set_min_objective(myfunc)
    opt.set_lower_bounds(lb)
    opt.set_upper_bounds(ub)
    opt.set_maxeval(int(1e4))
    opt.set_xtol_rel(1e-4)
    local_opt = nlopt.opt(nlopt.LD_MMA, 2)
    opt.set_local_optimizer(local_opt)
    x0 = [0.0, 0.0]
    print("x0:", x0)
    try:
        x = opt.optimize(x0)
        minf = opt.last_optimum_value()
        print("optimum at ", x[0], x[1])
        print("minimum value = ", minf)
        print("result code = ", opt.last_optimize_result())
        print("num function evaluations:", numFuncEval[0])
        if minf_ext[0] < minf:
            raise ValueError(f"minimum value {minf} is not the true minimum {minf_ext[0]}")
    except nlopt.invalid_argument:
        # stogo/ags might not be enabled
        pass
