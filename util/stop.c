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

#include <math.h>
#include "nlopt-util.h"

/* utility routines to implement the various stopping criteria */

static int relstop(double vold, double vnew, double reltol, double abstol)
{
     if (nlopt_isinf(vold)) return 0;
     return(fabs(vnew - vold) < abstol 
	    || fabs(vnew - vold) < reltol * (fabs(vnew) + fabs(vold)) * 0.5
	    || (reltol > 0 && vnew == vold)); /* catch vnew == vold == 0 */
}

int nlopt_stop_ftol(const nlopt_stopping *s, double f, double oldf)
{
     return (relstop(oldf, f, s->ftol_rel, s->ftol_abs));
}

int nlopt_stop_f(const nlopt_stopping *s, double f, double oldf)
{
     return (f <= s->minf_max || nlopt_stop_ftol(s, f, oldf));
}

int nlopt_stop_x(const nlopt_stopping *s, const double *x, const double *oldx)
{
     unsigned i;
     for (i = 0; i < s->n; ++i)
	  if (!relstop(oldx[i], x[i], s->xtol_rel, s->xtol_abs[i]))
	       return 0;
     return 1;
}

int nlopt_stop_dx(const nlopt_stopping *s, const double *x, const double *dx)
{
     unsigned i;
     for (i = 0; i < s->n; ++i)
	  if (!relstop(x[i] - dx[i], x[i], s->xtol_rel, s->xtol_abs[i]))
	       return 0;
     return 1;
}

static double sc(double x, double smin, double smax)
{
     return smin + x * (smax - smin);
}

/* some of the algorithms rescale x to a unit hypercube, so we need to
   scale back before we can compare to the tolerances */
int nlopt_stop_xs(const nlopt_stopping *s,
		  const double *xs, const double *oldxs,
		  const double *scale_min, const double *scale_max)
{
     unsigned i;
     for (i = 0; i < s->n; ++i)
	  if (relstop(sc(oldxs[i], scale_min[i], scale_max[i]), 
		      sc(xs[i], scale_min[i], scale_max[i]),
		      s->xtol_rel, s->xtol_abs[i]))
	       return 1;
     return 0;
}

int nlopt_stop_evals(const nlopt_stopping *s)
{
     return (s->maxeval > 0 && s->nevals >= s->maxeval);
}

int nlopt_stop_time(const nlopt_stopping *s)
{
     return (s->maxtime > 0 && nlopt_seconds() - s->start >= s->maxtime);
}

int nlopt_stop_evalstime(const nlopt_stopping *stop)
{
     return nlopt_stop_evals(stop) || nlopt_stop_time(stop);
}

int nlopt_stop_forced(const nlopt_stopping *stop)
{
     return stop->force_stop && *(stop->force_stop);
}

unsigned nlopt_count_constraints(unsigned p, const nlopt_constraint *c)
{
     unsigned i, count = 0;
     for (i = 0; i < p; ++i)
	  count += c[i].m;
     return count;
}

unsigned nlopt_max_constraint_dim(unsigned p, const nlopt_constraint *c)
{
     unsigned i, max_dim = 0;
     for (i = 0; i < p; ++i)
	  if (c[i].m > max_dim)
	       max_dim = c[i].m;
     return max_dim;
}

void nlopt_eval_constraint(double *result, double *grad,
			   const nlopt_constraint *c,
			   unsigned n, const double *x)
{
     if (c->f)
	  result[0] = c->f(n, x, grad, c->f_data);
     else
	  c->mf(c->m, result, n, x, grad, c->f_data);
}
