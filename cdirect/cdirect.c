#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "nlopt-util.h"
#include "nlopt.h"
#include "cdirect.h"
#include "redblack.h"
#include "config.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/***************************************************************************/
/* basic data structure:
 *
 * a hyper-rectangle is stored as an array of length L = 2n+2, where [0]
 * is the value (f) of the function at the center, [1] is the "size"
 * measure (d) of the rectangle, [2..n+1] are the coordinates of the
 * center (c), and [n+2..2n+1] are the widths of the sides (w).
 *
 * a list of rectangles is just an array of N hyper-rectangles
 * stored as an N x L in row-major order.  Generally,
 * we allocate more than we need, allocating Na hyper-rectangles.
 *
 * n > 0 always
 */

#define RECT_LEN(n) (2*(n)+2) /* number of double values in a hyperrect */

/* parameters of the search algorithm and various information that
   needs to be passed around */
typedef struct {
     int n; /* dimension */
     int L; /* RECT_LEN(n) */
     double *rects; /* the hyper-rectangles */
     double magic_eps; /* Jones' epsilon parameter (1e-4 is recommended) */
     int which_diam; /* which measure of hyper-rectangle diam to use:
			0 = Jones, 1 = Gablonsky */
     int which_div; /* which way to divide rects:
		       0: Gablonsky (cubes divide all, rects longest)
		       1: orig. Jones (divide all longest sides)
		       2: Jones Encyc. Opt.: pick random longest side */
     const double *lb, *ub;
     nlopt_stopping *stop; /* stopping criteria */
     nlopt_func f; void *f_data;
     double *work; /* workspace, of length >= 2*n */
     int *iwork, iwork_len; /* workspace, of length iwork_len >= n */
     double fmin, *xmin; /* minimum so far */
     
     /* red-black tree of hyperrect indices, sorted by (d,f) in
	lexographical order */
     rb_tree rtree;
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
	  return maxw * 0.5; /* half-width of longest side */
     }
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

#define EQUAL_SIDE_TOL 5e-2 /* tolerance to equate side sizes */

/* divide rectangle idiv in the list p->rects */
static nlopt_result divide_rect(int *N, int *Na, int idiv, params *p)
{
     int i;
     const const int n = p->n;
     const int L = p->L;
     double *r = p->rects;
     double *c = r + L*idiv + 2; /* center of rect to divide */
     double *w = c + n; /* widths of rect to divide */
     double wmax = w[0];
     int imax = 0, nlongest = 0;
     rb_node *node;

     for (i = 1; i < n; ++i)
	  if (w[i] > wmax)
	       wmax = w[imax = i];
     for (i = 0; i < n; ++i)
	  if (wmax - w[i] <= wmax * EQUAL_SIDE_TOL)
	       ++nlongest;
     if (p->which_div == 1 || (p->which_div == 0 && nlongest == n)) {
	  /* trisect all longest sides, in increasing order of the average
	     function value along that direction */
	  double *fv = p->work;
	  int *isort = p->iwork;
	  for (i = 0; i < n; ++i) {
	       if (wmax - w[i] <= wmax * EQUAL_SIDE_TOL) {
		    double csave = c[i];
		    c[i] = csave - w[i] * THIRD;
		    FUNCTION_EVAL(fv[2*i], c, p);
		    c[i] = csave + w[i] * THIRD;
		    FUNCTION_EVAL(fv[2*i+1], c, p);
		    c[i] = csave;
	       }
	       else {
		    fv[2*i] = fv[2*i+1] = HUGE_VAL;
	       }
	  }
	  sort_fv(n, fv, isort);
	  ALLOC_RECTS(n, Na, r, (*N)+2*nlongest); 
	  p->rects = r; c = r + L*idiv + 2; w = c + n;
	  for (i = 0; i < nlongest; ++i) {
	       int k;
	       if (!(node = rb_tree_find(&p->rtree, idiv)))
		    return NLOPT_FAILURE;
	       w[isort[i]] *= THIRD;
	       r[L*idiv + 1] = rect_diameter(n, w, p);
	       rb_tree_resort(&p->rtree, node);
	       for (k = 0; k <= 1; ++k) {
		    memcpy(r + L*(*N) + 1, c-1, sizeof(double) * (2*n+1));
		    r[L*(*N) + 2 + isort[i]] += w[isort[i]] * (2*k-1);
		    r[L*(*N)] = fv[2*isort[i]+k];
		    if (!rb_tree_insert(&p->rtree, *N))
			 return NLOPT_FAILURE;
		    ++(*N);
	       }
	  }
     }
     else {
	  int k;
	  if (nlongest > 1 && p->which_div == 2) { 
               /* randomly choose longest side */
	       i = nlopt_iurand(nlongest);
	       for (k = 0; k < n; ++k)
		    if (wmax - w[k] <= wmax * EQUAL_SIDE_TOL) {
			 if (!i) { i = k; break; }
			 --i;
		    }
	  }
	  else
	       i = imax; /* trisect longest side */
	  ALLOC_RECTS(n, Na, r, (*N)+2);
          p->rects = r; c = r + L*idiv + 2; w = c + n;
	  if (!(node = rb_tree_find(&p->rtree, idiv)))
	       return NLOPT_FAILURE;
	  w[i] *= THIRD;
	  r[L*idiv + 1] = rect_diameter(n, w, p);
	  rb_tree_resort(&p->rtree, node);
	  for (k = 0; k <= 1; ++k) {
	       memcpy(r + L*(*N) + 1, c-1, sizeof(double) * (2*n+1));
	       r[L*(*N) + 2 + i] += w[i] * (2*k-1);
	       FUNCTION_EVAL(r[L*(*N)], r + L*(*N) + 2, p);
	       if (!rb_tree_insert(&p->rtree, *N))
		    return NLOPT_FAILURE;
	       ++(*N);
	  }
     }
     return NLOPT_SUCCESS;
}

/***************************************************************************/
/* O(N log N) convex hull algorithm, used later to find the potentially
   optimal points */

/* Find the lower convex hull of a set of points (xy[s*i+1], xy[s*i]), where
   0 <= i < N and s >= 2.

   The return value is the number of points in the hull, with indices
   stored in ihull.  ihull should point to arrays of length >= N.
   rb_tree should be a red-black tree of indices (keys == i) sorted
   in lexographic order by (xy[s*i+1], xy[s*i]).
*/
static int convex_hull(int N, double *xy, int s, int *ihull, rb_tree *t)
{
     int nhull;
     double minslope;
     double xmin, xmax, yminmin, ymaxmin;
     rb_node *n, *nmax;

     if (N == 0) return 0;
     
     /* Monotone chain algorithm [Andrew, 1979]. */

     n = rb_tree_min(t);
     nmax = rb_tree_max(t);

     xmin = xy[s*(n->k)+1];
     yminmin = xy[s*(n->k)];
     xmax = xy[s*(nmax->k)+1];

     ihull[nhull = 1] = n->k;
     if (xmin == xmax) return nhull;

     /* set nmax = min mode with x == xmax */
     while (xy[s*(nmax->k)+1] == xmax)
	  nmax = rb_tree_pred(t, nmax); /* non-NULL since xmin != xmax */
     nmax = rb_tree_succ(t, nmax);

     ymaxmin = xy[s*(nmax->k)];
     minslope = (ymaxmin - yminmin) / (xmax - xmin);

     /* set n = first node with x != xmin */
     while (xy[s*(n->k)+1] == xmin)
	  n = rb_tree_succ(t, n); /* non-NULL since xmin != xmax */

     for (; n != nmax; n = rb_tree_succ(t, n)) { 
	  int k = n->k;
	  if (xy[s*k] > yminmin + (xy[s*k+1] - xmin) * minslope)
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
     ihull[nhull++] = nmax->k;
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

static nlopt_result divide_good_rects(int *N, int *Na, params *p)
{
     const int n = p->n;
     const int L = p->L;
     int *ihull, nhull, i, xtol_reached = 1, divided_some = 0;
     double *r = p->rects;
     double magic_eps = p->magic_eps;

     if (p->iwork_len < *N) {
	  p->iwork_len = p->iwork_len + *N;
	  p->iwork = (int *) realloc(p->iwork, sizeof(int) * p->iwork_len);
	  if (!p->iwork)
	       return NLOPT_OUT_OF_MEMORY;
     }
     ihull = p->iwork;
     nhull = convex_hull(*N, r, L, ihull, &p->rtree);
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
	       ret = divide_rect(N, Na, ihull[i], p);
	       r = p->rects; /* may have grown */
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
	       return divide_rect(N, Na, imax, p);
	  }
     }
     return xtol_reached ? NLOPT_XTOL_REACHED : NLOPT_SUCCESS;
}

/***************************************************************************/

/* lexographic sort order (d,f) of hyper-rects, for red-black tree */
static int hyperrect_compare(int i, int j, void *p_)
{
     params *p = (params *) p_;
     int L = p->L;
     double *r = p->rects;
     double di = r[i*L+1], dj = r[j*L+1], fi, fj;
     if (di < dj) return -1;
     if (dj < di) return +1;
     fi = r[i*L]; fj = r[j*L];
     if (fi < fj) return -1;
     if (fj < fi) return +1;
     return 0;
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
     int Na = 100, N = 1, i, x_center = 1;
     nlopt_result ret = NLOPT_OUT_OF_MEMORY;

     p.magic_eps = magic_eps;
     p.which_diam = which_alg % 2;
     p.which_div = 0;
     p.lb = lb; p.ub = ub;
     p.stop = stop;
     p.n = n;
     p.L = RECT_LEN(n);
     p.f = f;
     p.f_data = f_data;
     p.xmin = x;
     p.fmin = f(n, x, NULL, f_data); stop->nevals++;
     p.work = 0;
     p.iwork = 0;
     p.rects = 0;

     if (!rb_tree_init(&p.rtree, hyperrect_compare, &p)) goto done;
     p.work = (double *) malloc(sizeof(double) * 2*n);
     if (!p.work) goto done;
     p.rects = (double *) malloc(sizeof(double) * Na * RECT_LEN(n));
     if (!p.rects) goto done;
     p.iwork = (int *) malloc(sizeof(int) * (p.iwork_len = Na + n));
     if (!p.iwork) goto done;

     for (i = 0; i < n; ++i) {
	  p.rects[2+i] = 0.5 * (lb[i] + ub[i]);
	  x_center = x_center
	       && (fabs(p.rects[2+i]-x[i]) < 1e-13*(1+fabs(x[i])));
	  p.rects[2+n+i] = ub[i] - lb[i];
     }
     p.rects[1] = rect_diameter(n, p.rects+2+n, &p);
     if (x_center)
	  p.rects[0] = p.fmin; /* avoid computing f(center) twice */
     else
	  p.rects[0] = function_eval(p.rects+2, &p);
     if (!rb_tree_insert(&p.rtree, 0)) {
	  ret = NLOPT_FAILURE;
	  goto done;
     }

     ret = divide_rect(&N, &Na, 0, &p);
     if (ret != NLOPT_SUCCESS) goto done;

     while (1) {
	  double fmin0 = p.fmin;
	  ret = divide_good_rects(&N, &Na, &p);
	  if (ret != NLOPT_SUCCESS) goto done;
	  if (nlopt_stop_f(p.stop, p.fmin, fmin0)) {
	       ret = NLOPT_FTOL_REACHED;
	       goto done;
	  }
     }

 done:
     rb_tree_destroy(&p.rtree);
     free(p.iwork);
     free(p.rects);
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
