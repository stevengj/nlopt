#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "nlopt-util.h"
#include "nlopt.h"
#include "cdirect.h"
#include "config.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/***************************************************************************/
/* basic data structure:
 *
 * a hyper-rectangle is stored as an array of length 2n+2, where [0]
 * is the value of the function at the center, [1] is the "size"
 * measure of the rectangle, [2..n+1] are the coordinates of the
 * center, and [n+2..2n+1] are the widths of the sides.
 *
 * a list of rectangles is just an array of N hyper-rectangles
 * stored as an N x (2n+1) in row-major order.  Generally,
 * we allocate more than we need, allocating Na hyper-rectangles.
 *
 * n > 0 always
 */

#define RECT_LEN(n) (2*(n)+2) /* number of double values in a hyperrect */

/* parameters of the search algorithm and various information that
   needs to be passed around */
typedef struct {
     double magic_eps; /* Jones' epsilon parameter (1e-4 is recommended) */
     int which_diam; /* which measure of hyper-rectangle diam to use:
			0 = Jones, 1 = Gablonsky */
     const double *lb, *ub;
     nlopt_stopping *stop; /* stopping criteria */
     int n;
     nlopt_func f; void *f_data;
     double *work; /* workspace, of length >= 2*n */
     int *iwork, iwork_len; /* workspace, of length iwork_len >= n */
     double fmin, *xmin; /* minimum so far */
} params;

/***************************************************************************/

/* evaluate the "diameter" (d) of a rectangle of widths w[n] */
static double rect_diameter(int n, const double *w, const params *p)
{
     int i;
     if (p->which_diam == 0) { /* Jones measure */
	  double sum = 0;
	  for (i = 0; i < n; ++i)
	       sum += w[i] * w[i];
	  return sqrt(sum) * 0.5; /* distance from center to a vertex */
     }
     else { /* Gablonsky measure */
	  double maxw = 0;
	  for (i = 0; i < n; ++i)
	       if (w[i] > maxw)
		    maxw = w[i];
	  return w[i] * 0.5; /* half-width of longest side */
     }
}

#define CUBE_TOL 5e-2 /* fractional tolerance to call something a "cube" */

/* return true if the elements of w[n] (> 0) are all equal to within a
   fractional tolerance of tol (i.e. they are the widths of a hypercube) */
static int all_equal(int n, const double *w, double tol)
{
     double wmin, wmax;
     int i;
     wmin = wmax = w[0];
     for (i = 1; i < n; ++i) {
	  if (w[i] < wmin) wmin = w[i];
	  if (w[i] > wmax) wmax = w[i];
     }
     return (wmax - wmin) < tol * wmax;
}

static double *alloc_rects(int n, int *Na, double *rects, int newN)
{
     if (newN <= *Na)
	  return rects;
     else {
	  (*Na) += newN;
	  return realloc(rects, sizeof(double) * RECT_LEN(n) * (*Na));
     }
}
#define ALLOC_RECTS(n, Nap, rects, newN) if (!(rects = alloc_rects(n, Nap, rects, newN))) return NLOPT_OUT_OF_MEMORY

static double *fv_qsort = 0;
static int sort_fv_compare(const void *a_, const void *b_)
{
     int a = *((const int *) a_), b = *((const int *) b_);
     double fa = MIN(fv_qsort[2*a], fv_qsort[2*a+1]);
     double fb = MIN(fv_qsort[2*b], fv_qsort[2*b+1]);
     if (fa < fb)
	  return -1;
     else if (fa > fb)
	  return +1;
     else
	  return 0;
}
static void sort_fv(int n, double *fv, int *isort)
{
     int i;
     for (i = 0; i < n; ++i) isort[i] = i;
     fv_qsort = fv; /* not re-entrant, sigh... */
     qsort(isort, (unsigned) n, sizeof(int), sort_fv_compare);
     fv_qsort = 0;
}

static double function_eval(const double *x, params *p) {
     double f = p->f(p->n, x, NULL, p->f_data);
     if (f < p->fmin) {
	  p->fmin = f;
	  memcpy(p->xmin, x, sizeof(double) * p->n);
     }
     p->stop->nevals++;
     return f;
}
#define FUNCTION_EVAL(fv,x,p) fv = function_eval(x, p); if (p->fmin < p->stop->fmin_max) return NLOPT_FMIN_MAX_REACHED; else if (nlopt_stop_evals((p)->stop)) return NLOPT_MAXEVAL_REACHED; else if (nlopt_stop_time((p)->stop)) return NLOPT_MAXTIME_REACHED

#define THIRD (0.3333333333333333333333)


/* divide rectangle idiv in the list rects */
static nlopt_result divide_rect(int *N, int *Na, double **rects, int idiv,
				params *p)
{
     int i;
     const const int n = p->n;
     const int L = RECT_LEN(n);
     double *r = *rects;
     double *c = r + L*idiv + 2; /* center of rect to divide */
     double *w = c + n; /* widths of rect to divide */

     if (all_equal(n, w, CUBE_TOL)) { /* divide all dimensions */
	  double *fv = p->work;
	  int *isort = p->iwork;
	  for (i = 0; i < n; ++i) {
	       double csave = c[i];
	       c[i] = csave - w[i] * THIRD;
	       FUNCTION_EVAL(fv[2*i], c, p);
	       c[i] = csave + w[i] * THIRD;
	       FUNCTION_EVAL(fv[2*i+1], c, p);
	       c[i] = csave;
	  }
	  sort_fv(n, fv, isort);
	  ALLOC_RECTS(n, Na, r, (*N)+2*n); 
	  *rects = r; c = r + L*idiv + 2; w = c + n;
	  for (i = 0; i < n; ++i) {
	       int k;
	       w[isort[i]] *= THIRD;
	       r[L*idiv + 1] = rect_diameter(n, w, p);
	       for (k = 0; k <= 1; ++k) {
		    memcpy(r + L*(*N) + 1, c-1, sizeof(double) * (2*n+1));
		    r[L*(*N) + 2 + isort[i]] += w[isort[i]] * (2*k-1);
		    r[L*(*N)] = fv[2*isort[i]+k];
		    ++(*N);
	       }
	  }
     }
     else { /* divide longest side by 3 and split off 2 new rectangles */
	  int imax = 0;
	  double wmax = w[0];
	  for (i = 1; i < n; ++i)
	       if (w[i] > wmax)
		    wmax = w[imax = i];
	  ALLOC_RECTS(n, Na, r, (*N)+2);
	  *rects = r; c = r + L*idiv + 2; w = c + n;
	  w[imax] *= THIRD;
	  r[L*idiv + 1] = rect_diameter(n, w, p);
	  for (i = 0; i <= 1; ++i) {
	       memcpy(r + L*(*N) + 1, c-1, sizeof(double) * (2*n+1));
	       r[L*(*N) + 2 + imax] += w[imax] * (2*i-1); /* move center */
	       ++(*N);
	       FUNCTION_EVAL(r[L*((*N)-1)], r + L*((*N)-1) + 2, p);
	  }
     }
     return NLOPT_SUCCESS;
}

/***************************************************************************/
/* O(N log N) convex hull algorithm, used later to find the potentially
   optimal points */

/* sort ihull by xy in lexographic order by x,y */
static int s_qsort = 1; static double *xy_qsort = 0;
static int sort_xy_compare(const void *a_, const void *b_)
{
     int a = *((const int *) a_), b = *((const int *) b_);
     double xa = xy_qsort[a*s_qsort+1], xb = xy_qsort[b*s_qsort+1];
     if (xa < xb) return -1;
     else if (xb < xa) return +1;
     else {
	  double ya = xy_qsort[a*s_qsort], yb = xy_qsort[b*s_qsort];
	  if (ya < yb) return -1;
	  else if (ya > yb) return +1;
	  else return 0;
     }
}
static void sort_xy(int N, double *xy, int s, int *isort)
{
     int i;

     for (i = 0; i < N; ++i) isort[i] = i;
     s_qsort = s; xy_qsort = xy;
     qsort(isort, (unsigned) N, sizeof(int), sort_xy_compare);
     xy_qsort = 0;
}

/* Find the lower convex hull of a set of points (xy[s*i+1], xy[s*i]), where
   0 <= i < N and s >= 2.

   The return value is the number of points in the hull, with indices
   stored in ihull.  ihull and is should point to arrays of length >= N.

   Note that we don't allow redundant points along the same line in the
   hull, similar to Gablonsky's version of DIRECT and differing from
   Jones'. */
static int convex_hull(int N, double *xy, int s, int *ihull, int *is)
{
     int minmin; /* first index (smallest y) with min x */
     int minmax; /* last index (largest y) with min x */
     int maxmin; /* first index (smallest y) with max x */
     int maxmax; /* last index (largest y) with max x */
     int i, nhull = 0;
     double minslope;
     double xmin, xmax;

     /* Monotone chain algorithm [Andrew, 1979]. */

     sort_xy(N, xy, s, is);

     xmin = xy[s*is[minmin=0]+1]; xmax = xy[s*is[maxmax=N-1]+1];

     if (xmin == xmax) { /* degenerate case */
	  ihull[nhull++] = is[minmin];
	  return nhull;
     }

     for (minmax = minmin; minmax+1 < N && xy[s*is[minmax+1]+1]==xmin; 
	  ++minmax);
     for (maxmin = maxmax; maxmin-1>=0 && xy[s*is[maxmin-1]+1]==xmax;
	  --maxmin);

     minslope = (xy[s*is[maxmin]] - xy[s*is[minmin]]) / (xmax - xmin);

     ihull[nhull++] = is[minmin];
     for (i = minmax + 1; i < maxmin; ++i) {
	  int k = is[i];
	  if (xy[s*k] > xy[s*is[minmin]] + (xy[s*k+1] - xmin) * minslope)
	       continue;
	  /* remove points until we are making a "left turn" to k */
	  while (nhull > 1) {
	       int t1 = ihull[nhull - 1], t2 = ihull[nhull - 2];
	       /* cross product (t1-t2) x (k-t2) > 0 for a left turn: */
	       if ((xy[s*t1+1]-xy[s*t2+1]) * (xy[s*k]-xy[s*t2])
		   - (xy[s*t1]-xy[s*t2]) * (xy[s*k+1]-xy[s*t2+1]) > 0)
		    break;
	       --nhull;
	  }
	  ihull[nhull++] = k;
     }
     ihull[nhull++] = is[maxmin];
     return nhull;
}

/***************************************************************************/

static int small(double *w, params *p)
{
     int i;
     for (i = 0; i < p->n; ++i)
	  if (w[i] > p->stop->xtol_abs[i] &&
	      w[i] > (p->ub[i] - p->lb[i]) * p->stop->xtol_rel)
	       return 0;
     return 1;
}

static nlopt_result divide_good_rects(int *N, int *Na, double **rects, 
				      params *p)
{
     const int n = p->n;
     const int L = RECT_LEN(n);
     int *ihull, nhull, i, xtol_reached = 1, divided_some = 0;
     double *r = *rects;
     double magic_eps = p->magic_eps;

     if (p->iwork_len < n + 2*(*N)) {
	  p->iwork_len = p->iwork_len + n + 2*(*N);
	  p->iwork = (int *) realloc(p->iwork, sizeof(int) * p->iwork_len);
	  if (!p->iwork)
	       return NLOPT_OUT_OF_MEMORY;
     }
     ihull = p->iwork;
     nhull = convex_hull(*N, r, L, ihull, ihull + *N);
 divisions:
     for (i = 0; i < nhull; ++i) {
	  double K1 = -HUGE_VAL, K2 = -HUGE_VAL, K;
	  if (i > 0)
	       K1 = (r[L*ihull[i]] - r[L*ihull[i-1]]) /
		    (r[L*ihull[i]+1] - r[L*ihull[i-1]+1]);
	  if (i < nhull-1)
	       K1 = (r[L*ihull[i]] - r[L*ihull[i+1]]) /
		    (r[L*ihull[i]+1] - r[L*ihull[i+1]+1]);
	  K = MAX(K1, K2);
	  if (r[L*ihull[i]] - K * r[L*ihull[i]+1]
	      <= p->fmin - magic_eps * fabs(p->fmin)) {
	       /* "potentially optimal" rectangle, so subdivide */
	       divided_some = 1;
	       nlopt_result ret;
	       ret = divide_rect(N, Na, rects, ihull[i], p);
	       r = *rects; /* may have grown */
	       if (ret != NLOPT_SUCCESS) return ret;
	       xtol_reached = xtol_reached && small(r + L*ihull[i] + 2+n, p);
	  }
     }
     if (!divided_some) {
	  if (magic_eps != 0) {
	       magic_eps = 0;
	       goto divisions; /* try again */
	  }
	  else { /* WTF? divide largest rectangle */
	       double wmax = r[1];
	       int imax = 0;
	       for (i = 1; i < *N; ++i)
		    if (r[L*i+1] > wmax)
			 wmax = r[L*(imax=i)+1];
	       return divide_rect(N, Na, rects, imax, p);
	  }
     }
     return xtol_reached ? NLOPT_XTOL_REACHED : NLOPT_SUCCESS;
}

/***************************************************************************/

nlopt_result cdirect_unscaled(int n, nlopt_func f, void *f_data,
			      const double *lb, const double *ub,
			      double *x,
			      double *fmin,
			      nlopt_stopping *stop,
			      double magic_eps, int which_alg)
{
     params p;
     double *rects;
     int Na = 100, N = 1, i, x_center = 1;
     nlopt_result ret = NLOPT_OUT_OF_MEMORY;

     p.magic_eps = magic_eps;
     p.which_diam = which_alg & 1;
     p.lb = lb; p.ub = ub;
     p.stop = stop;
     p.n = n;
     p.f = f;
     p.f_data = f_data;
     p.xmin = x;
     p.fmin = f(n, x, NULL, f_data); stop->nevals++;
     p.work = 0;
     p.iwork = 0;
     rects = 0;
     p.work = (double *) malloc(sizeof(double) * 2*n);
     if (!p.work) goto done;
     rects = (double *) malloc(sizeof(double) * Na * RECT_LEN(n));
     if (!rects) goto done;
     p.iwork = (int *) malloc(sizeof(int) * (p.iwork_len = 2*Na + n));
     if (!p.iwork) goto done;

     for (i = 0; i < n; ++i) {
	  rects[2+i] = 0.5 * (lb[i] + ub[i]);
	  x_center = x_center
	       && (fabs(rects[2+i]-x[i]) < 1e-13*(1+fabs(x[i])));
	  rects[2+n+i] = ub[i] - lb[i];
     }
     rects[1] = rect_diameter(n, rects+2+n, &p);
     if (x_center)
	  rects[0] = p.fmin; /* avoid computing f(center) twice */
     else
	  rects[0] = function_eval(rects+2, &p);

     ret = divide_rect(&N, &Na, &rects, 0, &p);
     if (ret != NLOPT_SUCCESS) goto done;

     while (1) {
	  double fmin0 = p.fmin;
	  ret = divide_good_rects(&N, &Na, &rects, &p);
	  if (ret != NLOPT_SUCCESS) goto done;
	  if (nlopt_stop_f(p.stop, p.fmin, fmin0)) {
	       ret = NLOPT_FTOL_REACHED;
	       goto done;
	  }
     }

 done:
     free(p.iwork);
     free(rects);
     free(p.work);
	      
     *fmin = p.fmin;
     return ret;
}

/* in the conventional DIRECT-type algorithm, we first rescale our
   coordinates to a unit hypercube ... we do this simply by
   wrapping cdirect() around cdirect_unscaled(). */

typedef struct {
     nlopt_func f;
     void *f_data;
     double *x;
     const double *lb, *ub;
} uf_data;
static double uf(int n, const double *xu, double *grad, void *d_)
{
     uf_data *d = (uf_data *) d_;
     double f;
     int i;
     for (i = 0; i < n; ++i)
	  d->x[i] = d->lb[i] + xu[i] * (d->ub[i] - d->lb[i]);
     f = d->f(n, d->x, grad, d->f_data);
     if (grad)
	  for (i = 0; i < n; ++i)
	       grad[i] *= d->ub[i] - d->lb[i];
     return f;
}

nlopt_result cdirect(int n, nlopt_func f, void *f_data,
                     const double *lb, const double *ub,
                     double *x,
                     double *fmin,
                     nlopt_stopping *stop,
                     double magic_eps, int which_alg)
{
     uf_data d;
     nlopt_result ret;
     const double *xtol_abs_save;
     int i;

     d.f = f; d.f_data = f_data; d.lb = lb; d.ub = ub;
     d.x = (double *) malloc(sizeof(double) * n*4);
     if (!d.x) return NLOPT_OUT_OF_MEMORY;
     
     for (i = 0; i < n; ++i) {
	  x[i] = (x[i] - lb[i]) / (ub[i] - lb[i]);
	  d.x[n+i] = 0;
	  d.x[2*n+i] = 1;
	  d.x[3*n+i] = stop->xtol_abs[i] / (ub[i] - lb[i]);
     }
     xtol_abs_save = stop->xtol_abs;
     stop->xtol_abs = d.x + 3*n;
     ret = cdirect_unscaled(n, uf, &d, d.x+n, d.x+2*n, x, fmin, stop,
			    magic_eps, which_alg);
     stop->xtol_abs = xtol_abs_save;
     for (i = 0; i < n; ++i)
	  x[i] = lb[i]+ x[i] * (ub[i] - lb[i]);
     free(d.x);
     return ret;
}
