#!/usr/bin/env python

import nlopt
import math as m
import sys


def myfunc(x, grad):
    if grad.size > 0:
        grad[0] = 0.0
        grad[1] = 0.5 / m.sqrt(x[1])
    return m.sqrt(x[1])


def myconstraint(x, grad, a, b):
    if grad.size > 0:
        grad[0] = 3 * a * (a*x[0] + b)**2
        grad[1] = -1.0
    return (a*x[0] + b)**3 - x[1]


algo = nlopt.LD_MMA if len(sys.argv) < 2 else int(sys.argv[1])
opt = nlopt.opt(algo, 2)
print(f"algo: {opt.get_algorithm_name()}")
opt.set_lower_bounds([-float('inf'), 1e-6])
opt.set_min_objective(myfunc)
opt.add_inequality_constraint(lambda x, grad: myconstraint(x, grad, 2, 0), 1e-8)
opt.add_inequality_constraint(lambda x, grad: myconstraint(x, grad, -1, 1), 1e-8)
opt.set_xtol_rel(1e-4)
x0 = [1.234, 5.678]
x = opt.optimize(x0)
minf = opt.last_optimum_value()
print(f"optimum at {x}")
print(f"minimum value: {minf:.7g}")
print(f"result code: {opt.last_optimize_result()}")
print(f"nevals: {opt.get_numevals()}")
print(f"initial step: {opt.get_initial_step(x0)}")
assert m.fabs(minf - 0.544331) < 1e-3, "wrong optimum"
