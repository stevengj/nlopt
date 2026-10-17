#ifndef DIRECT_PYTHON_H
#define DIRECT_PYTHON_H

/* Thin layer translating the C API in direct.h to a PyObject* API,
   for use by scipy.optimize.direct.  Not part of the NLopt build. */

#include <Python.h>

#include "direct.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/* f(x_seq) or f(x_seq, args) is called for each evaluation, and
   callback(x_seq) after each iteration (unless callback is Py_None),
   where x_seq is a list of length dimension holding the point.  On
   return, x and x_seq hold the minimizer.  Returns NULL if a Python
   exception was raised, otherwise a new reference to Py_None. */
extern PyObject* direct_optimize_py(
    PyObject* f, double *x, PyObject *x_seq, PyObject* args,
    int dimension,
    const double *lower_bounds, const double *upper_bounds,
    double *minf,
    int max_feval, int max_iter,
    double magic_eps, double magic_eps_abs,
    double volume_reltol, double sigma_reltol,
    int *force_stop,
    double fglobal,
    double fglobal_reltol,
    FILE *logfile,
    direct_algorithm algorithm,
    direct_return_info *info,
    direct_return_code *ret_code,
    PyObject* callback);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* DIRECT_PYTHON_H */
