#!/usr/bin/env python

import nlopt
import numpy as np

print ('nlopt version='+nlopt.__version__)

def f(x, grad):
    F=x[0]
    L=x[1]
    E=x[2]
    I=x[3]
    D=F*L**3/(3.*E*I)
    return D

n = 4
opt = nlopt.opt(nlopt.LN_COBYLA, n)
opt.set_min_objective(f)
lb = np.array([40., 50., 30e3, 1.])
ub = np.array([60., 60., 40e3, 10.])
x = (lb+ub)/2.
opt.set_lower_bounds(lb)
opt.set_upper_bounds(ub)
opt.set_xtol_rel(1e-3)
opt.set_ftol_rel(1e-3)
xopt = opt.optimize(x)

opt_val = opt.last_optimum_value()
result = opt.last_optimize_result()
print ('opt_result='+str(result))
print ('optimizer='+str(xopt))
print ('opt_val='+str(opt_val))
