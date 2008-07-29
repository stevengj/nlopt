#include <math.h>
#include "nlopt-util.h"

/* utility routines to implement the various stopping criteria */

static int relstop(double old, double new, double reltol, double abstol)
{
     return(fabs(new - old) < abstol 
	    || fabs(new - old) < reltol * (fabs(new) + fabs(old)) * 0.5);
}

int nlopt_stop_ftol(const nlopt_stopping *s, const double f, double oldf)
{
     return (relstop(oldf, f, s->ftol_rel, s->ftol_abs));
}

int nlopt_stop_f(const nlopt_stopping *s, const double f, double oldf)
{
     return (f <= s->minf_max || nlopt_stop_ftol(s, f, oldf));
}

int nlopt_stop_x(const nlopt_stopping *s, const double *x, const double *oldx)
{
     int i;
     for (i = 0; i < s->n; ++i)
	  if (!relstop(oldx[i], x[i], s->xtol_rel, s->xtol_abs[i]))
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
     int i;
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
