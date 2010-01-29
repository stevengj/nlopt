/* cobyla : contrained optimization by linear approximation */

/*
 * Copyright (c) 1992, Michael J. D. Powell (M.J.D.Powell@damtp.cam.ac.uk)
 * Copyright (c) 2004, Jean-Sebastien Roy (js@jeannot.org)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * This software is a C version of COBYLA2, a contrained optimization by linear
 * approximation package developed by Michael J. D. Powell in Fortran.
 * 
 * The original source code can be found at :
 * http://plato.la.asu.edu/topics/problems/nlores.html
 */

/* $Jeannot: cobyla.h,v 1.10 2004/04/18 09:51:37 js Exp $ */

#ifndef _COBYLA_
#define _COBYLA_

#include "nlopt.h"
#include "nlopt-util.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/*
 * Verbosity level
 */
typedef enum {
  COBYLA_MSG_NONE = 0, /* No messages */
  COBYLA_MSG_EXIT = 1, /* Exit reasons */
  COBYLA_MSG_ITER = 2, /* Rho and Sigma changes */
  COBYLA_MSG_INFO = 3, /* Informational messages */
} cobyla_message;

/*
 * A function as required by cobyla
 * state is a void pointer provided to the function at each call
 *
 * n     : the number of variables
 * m     : the number of constraints
 * x     : on input, then vector of variables (should not be modified)
 * f     : on output, the value of the function
 * con   : on output, the value of the constraints (vector of size m)
 * state : on input, the value of the state variable as provided to cobyla
 *
 * COBYLA will try to make all the values of the constraints positive.
 * So if you want to input a constraint j such as x[i] <= MAX, set:
 *   con[j] = MAX - x[i]
 * The function must returns 0 if no error occurs or 1 to immediately end the
 * minimization.
 *
 */
typedef int cobyla_function(int n, int m, double *x, double *f, double *con,
  void *state);

/*
 * cobyla : minimize a function subject to constraints
 *
 * n         : number of variables (>=0)
 * m         : number of constraints (>=0)
 * x         : on input, initial estimate ; on output, the solution
 * minf      : on output, minimum objective function
 * rhobeg    : a reasonable initial change to the variables
 * stop      : the NLopt stopping criteria
 * lb, ub    : lower and upper bounds on x
 * message   : see the cobyla_message enum
 * calcfc    : the function to minimize (see cobyla_function)
 * state     : used by function (see cobyla_function)
 *
 * The cobyla function returns the usual nlopt_result codes.
 *
 */
extern nlopt_result cobyla(int n, int m, double *x, double *minf, double rhobeg, nlopt_stopping *stop, const double *lb, const double *ub,
  int message, cobyla_function *calcfc, void *state);

/* NLopt-style interface function */
nlopt_result cobyla_minimize(int n, nlopt_func f, void *f_data,
                             int m, nlopt_func fc,
                             void *fc_data_, ptrdiff_t fc_datum_size,
                             const double *lb, const double *ub, /* bounds */
                             double *x, /* in: initial guess, out: minimizer */
                             double *minf,
                             nlopt_stopping *stop,
                             double rhobegin);

#ifdef __cplusplus
}
#endif

#endif /* _COBYLA_ */
