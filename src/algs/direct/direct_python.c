#include "direct_python.h"

typedef struct {
    PyObject *f, *x_seq, *args, *callback;
    int *force_stop;
    int error; /* set if a Python exception is pending */
} direct_py_data;

/* copy x[n] into the list x_seq; returns 0 if a Python exception was raised */
static int set_x_seq(PyObject *x_seq, int n, const double *x)
{
    int i;
    for (i = 0; i < n; ++i) {
        PyObject *xi = PyFloat_FromDouble(x[i]);
        if (!xi || PyList_SetItem(x_seq, i, xi) < 0)
            return 0;
    }
    return 1;
}

static double direct_py_objective(int n, const double *x, int *undefined_flag,
                                  void *data)
{
    direct_py_data *d = (direct_py_data *) data;
    PyObject *f_py;
    double f;

    (void) undefined_flag;
    // TODO: PyArray_SimpleNewFromData gives segmentation fault
    // and therefore using list to pass to the user function.
    // Once the above function works, replace with NumPy arrays.
    if (!set_x_seq(d->x_seq, n, x))
        goto fail;
    if (PyObject_IsTrue(d->args))
        f_py = PyObject_CallFunctionObjArgs(d->f, d->x_seq, d->args, NULL);
    else
        f_py = PyObject_CallFunctionObjArgs(d->f, d->x_seq, NULL);
    if (!f_py)
        goto fail;
    f = PyFloat_AsDouble(f_py);
    Py_DECREF(f_py);
    if (f == -1.0 && PyErr_Occurred())
        goto fail;
    return f;

 fail:
    /* stop DIRECT so that the exception can be propagated */
    d->error = 1;
    *d->force_stop = 1;
    return HUGE_VAL;
}

static int direct_py_callback(int n, const double *x, void *data)
{
    direct_py_data *d = (direct_py_data *) data;
    PyObject *ret;

    if (!set_x_seq(d->x_seq, n, x))
        goto fail;
    ret = PyObject_CallFunctionObjArgs(d->callback, d->x_seq, NULL);
    if (!ret)
        goto fail;
    Py_DECREF(ret);
    return 0;

 fail:
    d->error = 1;
    return 1;
}

PyObject* direct_optimize_py(
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
    PyObject* callback)
{
    int stop = 0;
    direct_py_data d;

    if (!force_stop) force_stop = &stop;
    d.f = f;
    d.x_seq = x_seq;
    d.args = args;
    d.callback = callback;
    d.force_stop = force_stop;
    d.error = 0;

    *ret_code = direct_optimize(direct_py_objective, &d, dimension,
                                lower_bounds, upper_bounds, x, minf,
                                max_feval, max_iter,
                                magic_eps, magic_eps_abs,
                                volume_reltol, sigma_reltol,
                                force_stop,
                                fglobal, fglobal_reltol,
                                logfile, algorithm, info,
                                callback == Py_None ? NULL : direct_py_callback);

    if (d.error || !set_x_seq(x_seq, dimension, x))
        return NULL;
    Py_RETURN_NONE;
}
