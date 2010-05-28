/* Copyright (c) 2007-2010 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

#ifndef CDIRECT_H
#define CDIRECT_H

#include "nlopt-util.h"
#include "nlopt.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

extern nlopt_result cdirect_unscaled(int n, nlopt_func f, void *f_data,
				     const double *lb, const double *ub,
				     double *x,
				     double *minf,
				     nlopt_stopping *stop,
				     double magic_eps, int which_alg);

extern nlopt_result cdirect(int n, nlopt_func f, void *f_data,
			    const double *lb, const double *ub,
			    double *x,
			    double *minf,
			    nlopt_stopping *stop,
			    double magic_eps, int which_alg);

extern nlopt_result cdirect_hybrid(int n, nlopt_func f, void *f_data,
			    const double *lb, const double *ub,
			    double *x,
			    double *minf,
			    nlopt_stopping *stop,
			    nlopt_algorithm local_alg,
			    int local_maxeval,
			    int randomized_div);

extern nlopt_result cdirect_hybrid_unscaled(int n, nlopt_func f, void *f_data,
			    const double *lb, const double *ub,
			    double *x,
			    double *minf,
			    nlopt_stopping *stop,
			    nlopt_algorithm local_alg,
			    int local_maxeval,
			    int randomized_div);

/* internal routines and data structures: */
extern int cdirect_hyperrect_compare(double *a, double *b);
typedef struct {
     nlopt_func f;
     void *f_data;
     double *x;
     const double *lb, *ub;
} cdirect_uf_data;
extern double cdirect_uf(unsigned n, const double *xu, double *grad, void *d_);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* DIRECT_H */
