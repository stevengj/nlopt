/* SLSQP: Sequentional Least Squares Programming (aka sequential quadratic programming SQP)
   method for nonlinearly constrained nonlinear optimization, by Dieter Kraft (1991).
   Fortran released under a free (BSD) license by ACM to the SciPy project and used there.
   C translation via f2c + hand-cleanup and incorporation into NLopt by S. G. Johnson (2009).
   The f2c translation of the algorithm (main iteration, LSQ, LSEI, LSI, LDP, NNLS and the
   LDL' update) was later replaced by SciPy's C rewrite of the same routines
   (Ilhan Polat, scipy/scipy#22524, BSD license), which delegates the dense linear
   algebra to BLAS/LAPACK (or to the built-in versions in util/nlopt-blas-builtin.h
   when NLopt is configured without one).  The NLopt-specific changes of the original
   translation (bound clamping, non-finite line-search steps, no iteration limit,
   gradient at the first line-search trial point) were carried over. */

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "slsqp.h"

/*      ALGORITHM 733, COLLECTED ALGORITHMS FROM ACM. */
/*      TRANSACTIONS ON MATHEMATICAL SOFTWARE, */
/*      VOL. 20, NO. 3, SEPTEMBER, 1994, PP. 262-281. */
/*      http://doi.acm.org/10.1145/192115.192124 */


/*      http://permalink.gmane.org/gmane.comp.python.scientific.devel/6725 */
/*      ------ */
/*      From: Deborah Cotton <cotton@hq.acm.org> */
/*      Date: Fri, 14 Sep 2007 12:35:55 -0500 */
/*      Subject: RE: Algorithm License requested */
/*      To: Alan Isaac */

/*      Prof. Issac, */

/*      In that case, then because the author consents to [the ACM] releasing */
/*      the code currently archived at http://www.netlib.org/toms/733 under the */
/*      BSD license, the ACM hereby releases this code under the BSD license. */

/*      Regards, */

/*      Deborah Cotton, Copyright & Permissions */
/*      ACM Publications */
/*      2 Penn Plaza, Suite 701** */
/*      New York, NY 10121-0701 */
/*      permissions@acm.org */
/*      212.869.7440 ext. 652 */
/*      Fax. 212.869.0481 */
/*      ------ */

#include "nlopt-blas.h"

/********************************* NNLS ***********************************/

/* Lawson & Hanson nonnegative least squares, min |Ax - b| s.t. x >= 0.
   a is (m x n), b is (m), x, w, work are (n), zz is (m), indices is (n). */
static void
nnls(const blas_int m, const blas_int n, double* restrict a, double* restrict b,
     double* restrict x, double* restrict w, double* restrict zz, double* restrict work,
     blas_int* restrict indices, const int64_t maxiter, double* rnorm, int64_t* info)
{
    blas_int i = 0, ii = 0, ip = 0, indz = 0, iteration = 0, iz = 0, izmax = 0;
    blas_int j = 0, jj = 0, k = 0, one = 1, tmpint = 0;
    double tau = 0.0, unorm = 0.0, ztest, alpha, cc, ss, wmax, T, done = 1.0, dzero = 0.0;
    double pivot = 1.0, tauvtb = 0.0, tmp = 0.0, spacing = 0.0;
    *info = 1;
    if (m <= 0 || n <= 0)
    {
        *info = 2;
        return;
    }

    // Initialize the indices and the solution vector x.
    for (i = 0; i < n; i++) { indices[i] = i; }
    for (i = 0; i < n; i++) { x[i] = 0.0; }

    // Outer loop
    while (indz < (m < n ? m : n))
    {
        // Compute the dual vector components in set Z, in NumPy notation;
        // w[indices[indz:]] = A[indz:m, indices[indz:]].T @ b[indz:m]
        // as one GEMV over all columns (rather than one DOT per column of Z
        // as in scipy).  The columns in set P are zero in these rows; their dual
        // components stay zero.
        tmpint = m - indz;
        NLOPT_BLAS(dgemv)("T", &tmpint, &n, &done, &a[indz], &m, &b[indz], &one, &dzero, w, &one);
        for (i = 0; i < indz; i++) { w[indices[i]] = 0.0; }

        // Find the next linearly independent column that corresponds to the
        // largest entry in the dual vector w.
        // ====================================================================
        while (1)
        {
            // Finding the largest w[j] and its index
            // izmax, wmax = argmax(w[indices[indz:]])
            wmax = 0.0;
            for (k = indz; k < n; k++)
            {
                j = indices[k];
                if (w[j] > wmax) { wmax = w[j]; izmax = k; }
            }
            // If wmax <= 0.0, terminate since this is a KKT certificate.
            if (wmax <= 0.0) { goto END; }
            iz = izmax;
            j = indices[iz];

            // The sign of wmax is OK for j to be moved to set p. Column j and b
            // are not modified until j is accepted.

            // Compute the norm of a[0:indz, j] to check for linear dependence.
            unorm = (indz > 0 ? NLOPT_BLAS(dnrm2)(&indz, &a[j*m], &one) : 0.0);
            // unorm is nonnegative
            spacing = (unorm > 0.0 ? nextafter(unorm, 2*unorm) - unorm : 0.0);

            // Test for independence by checking the new diagonal element for
            // zero. Its magnitude is the norm of a[indz:m, j].
            tmpint = m - indz;
            if (NLOPT_BLAS(dnrm2)(&tmpint, &a[indz + j*m], &one) > 100.0*spacing)
            {
                // Column j is sufficiently independent. Form the transformation
                // on a copy in zz[indz:m] and solve for ztest which is the new
                // prospective value for x[j].
                for (i = indz; i < m; i++) { zz[i] = a[i + j*m]; }
                NLOPT_BLAS(dlarfgp)(&tmpint, &zz[indz], &zz[indz + 1], &one, &tau);
                pivot = zz[indz];
                zz[indz] = 1.0;
                // (H*b)[indz] = b[indz] - tau*v'b, with v = zz[indz:m].
                tmpint = m - indz - 1;
                tauvtb = tau*(b[indz] + NLOPT_BLAS(ddot)(&tmpint, &zz[indz + 1], &one, &b[indz + 1], &one));
                // See if ztest is positive. This is from the original F77 code and
                // guards against round-off making the new x[j] nonpositive.
                ztest = (b[indz] - tauvtb) / pivot;
                if (ztest > 0.0) { break; }
            }
            // Reject j as a candidate to be moved from set z to set p.
            // Set w(j)=0., and loop back to test dual coeffs again.
            w[j] = 0.0;
        }
        // ====================================================================

        // the index j=indices[iz]  has been selected to be moved from set z to
        // set p. Store the transformation in col j, apply it to b, update
        // indices, apply householder transformations to cols in new set z,
        // zero subdiagonal elements in col j,  set w(j)=0.
        for (i = indz; i < m; i++) { a[i + j*m] = zz[i]; }
        b[indz] -= tauvtb;
        tmpint = m - indz - 1;
        tmp = -tauvtb;
        NLOPT_BLAS(daxpy)(&tmpint, &tmp, &zz[indz + 1], &one, &b[indz + 1], &one);
        indices[iz] = indices[indz];
        indices[indz] = j;
        indz++;
        // Apply the householder transformation to the remaining columns.
        // The columns in set P other than j are zero in rows indz-1:m, which
        // the transformation leaves unchanged, so it is applied to all columns
        // left and right of column j as two blocks rather than to each column
        // of Z separately (as in scipy).
        if (indz < n)
        {
            blas_int nleft = j, nright = n - j - 1;
            tmpint = m - indz + 1;
            if (nleft > 0)
            {
                NLOPT_BLAS(dlarf)("L", &tmpint, &nleft, &a[indz - 1 + j*m], &one, &tau, &a[indz - 1], &m, work);
            }
            if (nright > 0)
            {
                NLOPT_BLAS(dlarf)("L", &tmpint, &nright, &a[indz - 1 + j*m], &one, &tau, &a[indz - 1 + (j + 1)*m], &m, work);
            }
        }
        // Restore the pivot element into a, zero the subdiagonal elements in col j
        a[indz - 1 + j*m] = pivot;
        if (indz < m) { for (i = indz; i < m; i++) { a[j*m + i] = 0.0; } }
        // Zero the dual coefficient for the column.
        w[j] = 0.0;

        // Solve the permuted triangular system, store in zz.
        for (i = 0; i < indz; i++) { zz[i] = b[i]; }
        for (k = 0; k < indz; k++)
        {
            // ip traverses the indices of P set in reverse
            ip = indz - 1 - k;
            if (k != 0)
            {
                for (i = 0; i <= ip; i++)
                {
                    zz[i] = zz[i] - a[i + jj*m] * zz[ip + 1];
                }
            }
            jj = indices[ip];
            zz[ip] = zz[ip] / a[ip + jj*m];
        }

        // ****** Inner loop ******
        while (1)
        {
            iteration++;
            if (iteration >= maxiter) { *info = 3; goto END; }

            // See if all new constrained coefficients are feasible,
            // if not compute alpha
            alpha = 2.0;
            for (ip = 0; ip < indz; ip++)
            {
                k = indices[ip];
                if (zz[ip] <= 0.0)
                {
                    T = -x[k] / (zz[ip] - x[k]);
                    if (alpha > T)
                    {
                        alpha = T;
                        jj = ip;
                    }
                }
            }
            // If all new constrained coefficients are feasible, alpha is still
            // 2.0. If so exit from secondary loop to main loop.
            if (alpha == 2.0) { break; }  // Get back to outer loop

            // Otherwise interpolate between old x and zz.
            for (ip = 0; ip < indz; ip++)
            {
                k = indices[ip];
                x[k] = x[k] + alpha*(zz[ip] - x[k]);
            }

            // Modify a, b, and the indices to move coefficient i from set p
            // to set z. While loop simulates a goto in the original F77 code.
            i = indices[jj];
            while (1)
            {
                x[i] = 0.0;

                if (jj != indz-1)
                {
                    jj++;
                    for (j = jj; j < indz; j++)
                    {
                        ii = indices[j];
                        indices[j-1] = ii;
                        NLOPT_BLAS(dlartgp)(&a[j-1 + ii*m], &a[j + ii*m], &cc, &ss, &a[j-1 + ii*m]);
                        a[j + ii*m] = 0.0;
                        // Apply the Givens rotation to all columns except ii.
                        // Because the columns are not ordered we do it manually.
                        for (k = 0; k < n; k++)
                        {
                            if (k != ii)
                            {
                                tmp = a[j-1 + k*m];
                                a[j-1 + k*m] =  cc*tmp + ss*a[j + k*m];
                                a[j   + k*m] = -ss*tmp + cc*a[j + k*m];
                            }
                        }
                        tmp = b[j-1];
                        b[j-1] =  cc*tmp + ss*b[j];
                        b[j]   = -ss*tmp + cc*b[j];
                    }
                }
                indz--;
                indices[indz] = i;

                // See if remaining coefficients in set P are feasible
                // since determination of alpha guarantees it. If still
                // there are infeasible ones, they are due to numerical
                // noise. Any that are nonpositive will be set to zero
                // and moved from set p to set z. If set p is now empty there is
                // nothing left to check (scipy loops on with indz = -1 here).
                int nobreak = (indz == 0);
                for (jj = 0; jj < indz; jj++)
                {
                    i = indices[jj];
                    if (x[i] <= 0.0) { break; }
                    if (jj == indz - 1) { nobreak = 1; }
                }
                // If for loop completes without break, then leave the while loop
                if (nobreak) { break; }
            }

            for (i = 0; i < m; i++) { zz[i] = b[i]; }
            for (k = 0; k < indz; k++)
            {
                ip = indz - 1 - k;
                if (k != 0)
                {
                    for (i = 0; i <= ip; i++)
                    {
                        zz[i] = zz[i] - a[i + jj*m] * zz[ip + 1];

                    }
                }
                jj = indices[ip];
                zz[ip] = zz[ip] / a[ip + jj*m];
            }
            // ****** end of inner loop ******
        }

        // Back in the outer loop
        for (k = 0; k < indz; k++)
        {
            i = indices[k];
            x[i] = zz[k];
        }
        // ****** end of outer loop ******
    }
END:
    // Compute the residual vector and its norm.
    if (indz < m)
    {
        tmpint = m - indz;
        *rnorm = NLOPT_BLAS(dnrm2)(&tmpint, &b[indz], &one);
    } else {
        for (i = 0; i < n; i++) { w[i] = 0.0; }
        *rnorm = 0.0;
    }
    return;
}

/********************************* SLSQP **********************************/

/* State of the reverse-communication SLSQP iteration, kept across calls. */
typedef struct {
    double acc, alpha, f0, gs, h1, h2, h3, h4, t, t0, tol;
    int64_t inconsistent, reset, iter, itermax, line, m, meq, mode, n;
    /* whether the current QP subproblem was inconsistent (augmented) */
    int badlin;
    /* SGJ 2010: evaluate the gradient together with the function at the
       first trial point of each line search (mode = -2) */
    int eager_grad;
} slsqp_state;

static void ldp(blas_int m, blas_int n, double* restrict g, double* restrict h, double* restrict x, double* restrict buffer, blas_int* indices, double* xnorm, int64_t* mode);
static void lsi(blas_int ma, blas_int mg, blas_int n, double* restrict a, double* restrict b, double* restrict g, double* restrict h, double* restrict x, double* restrict buffer, blas_int* jw, double* xnorm, int64_t* mode);
static void lsei(blas_int ma, blas_int me, blas_int mg, blas_int n, double* restrict a, double* restrict b, double* restrict e, double* restrict f, double* restrict g, double* restrict h, double* restrict x, double* restrict buffer, blas_int* jw, double* xnorm, int64_t* mode);
static void lsq(blas_int m, blas_int meq, blas_int n, blas_int augment, double aug_weight, const double* restrict Lf, const double* restrict gradx, const double* restrict C, const double* restrict d, double* restrict xl, double* restrict xu, double* restrict x, double* restrict y, double* buffer, blas_int* jw, int64_t* mode);
static void ldl_update(blas_int n, double* restrict a, double* restrict z, double sigma, double* restrict w);

/*
 * The main SLSQP function, called repeatedly in reverse-communication style.
 * Naming convention of the arguments:
 *
 *  - funx: The function value at the current point.
 *  - gradx: The gradient of the function at the current point. (n)
 *  - C: The equality and inequality constraint normals. (m x n, column major,
 *       leading dimension max(m, 1))
 *  - d: The  equality and inequality constraints, (m)
 *  - xl: The lower bounds on x, NaN if unbounded (n)
 *  - xu: The upper bounds on x, NaN if unbounded (n)
 *  - sol: The solution vector, (n)
 *  - mult: The Lagrange multipliers, (m + 2*n + 2)
 *  - buffer: A buffer to hold various intermediate arrays, see slsqp_buffer_size.
 *  - indices: An array to hold the indices of the active constraints. (m + 2*n + 2)
 *
 *  Constraints follow Kraft's convention: d[i] == 0 for i < meq, d[i] >= 0 otherwise.
 *
 *  On return S->mode is
 *   1: evaluate f and the constraints at sol, then call again;
 *  -2: (only with S->eager_grad) same as 1, but the gradients may be evaluated too;
 *  -1: evaluate the gradients at sol (f and the constraints are already known
 *      there), then call again;
 *   0: converged (never happens with S->acc == 0);
 *  2-9: failure modes, see nlopt_slsqp.
 *
 *  The first n*(n+1)/2 entries of buffer hold the LDL' factor of the BFGS
 *  matrix and the following n entries hold x0, the start point of the
 *  current line search.
 *
 *  If applicable, the following are the problem matrix naming convention:
 *
 *  - A: The coefficient matrix of cost function |Ax - b|
 *  - b: The RHS of cost function |Ax - b|
 *  - E: The (E)quality constraint matrix of Ex = f
 *  - f: The equality constraint RHS of Ex = f
 *  - G: The inequality constraint matrix of Gx >= h
 *  - h: The inequality constraint RHS of Gx >= h
 *
 */
static void
slsqp_body(
    slsqp_state* S, double funx, const double* restrict gradx,
    const double* restrict C, const double* restrict d, double* restrict sol,
    double* restrict mult, const double* restrict xl, const double* restrict xu,
    double* buffer, blas_int* indices)
{

    blas_int one = 1, lda = (S->m > 0 ? S->m : 1);
    blas_int j;
    double done = 1.0, dmone = -1.0, alfmin = 0.1;
    blas_int n = S->n;
    blas_int m = S->m;
    blas_int n1 = n + 1;
    blas_int n2 = n1*n/2;

    // Chop the buffer for various array pointers.
    double* restrict bfgs       = &buffer[0];
    double* restrict x0         = &buffer[n2];
    double* restrict mu         = &buffer[n2 + n];
    double* restrict s          = &buffer[n2 + n + m];
    double* restrict u          = &buffer[n2 + n + m + n1];
    double* restrict v          = &buffer[n2 + n + m + n1 + n1];
    double* restrict lsq_buffer = &buffer[n2 + n + m + n1 + n1 + n1];

    // Fortran code uses reverse communication for the iterations hence it
    // needs to jump back to where it left off. Thus the goto statements are
    // kept as is. Fortunately, they do not overlap too much and have a relatively
    // clean separation.
    if (S->mode == -1) { goto MODEM1; }
    if ((S->mode == 1) || (S->mode == -2)) { goto MODE1; }

    // MODE0:
    // We always use inexact line search, since exact search is broken in the
    // original Fortran code.
    S->acc = fabs(S->acc);
    S->tol = 10*S->acc;
    S->iter = 0;
    S->reset = 0;
    S->badlin = 0;
    for (int64_t i = 0; i < n; i++) { s[i] = 0.0; }
    for (int64_t i = 0; i < m; i++) { mu[i] = 0.0; }

RESET_BFGS:
    // Reset the BFGS matrix stored in packed format
    S->reset++;
    if (S->reset > 5) { goto LABEL255;}
    for (int64_t i = 0; i < n2; i++) { bfgs[i] = 0.0; }
    j = 0;
    for (int64_t i = 0; i < n; i++)
    {
        bfgs[j] = 1.0;
        j += n - i;
    }
    // 120

ITER_START:
    // Main iteration: Search direction, steplength, LDL'-update
    // 130
    S->mode = 9;
    // SGJ 2010: ignore if itermax <= 0
    if ((S->itermax > 0) && (S->iter >= S->itermax)) { return; }
    S->iter++;

    // Search direction as solution of the QP-problem
    for (int64_t i = 0; i < n; i++)
    {
        u[i] = -sol[i] + xl[i] ;
        v[i] = -sol[i] + xu[i] ;
    }

    S->h4 = 1.0;
    // augment and aug_weight are not used and hence 0.
    lsq(m, S->meq, n, 0, 0, bfgs, gradx, C, d, u, v, s, mult, lsq_buffer, indices, &S->mode);

    // Augmented problem for inconsistent linearization

    // If it turns out that the original SQP problem is inconsistent,
    // disallow termination with convergence on this iteration,
    // even if the augmented problem was solved.
    S->badlin = 0;

    // If equality constraints are not full rank and all are equality constrained
    // then the problem is inconsistent.
    if ((S->mode == 6) && (n == S->meq)) { S->mode = 4;}

    // If inconsistency detected, we augment the problem and try again.
    // Fortran code augments the problem matrices by embedding them in larger
    // buffers then calls lsq. However, these matrices are then copied into
    // another buffer inside lsq hence we can let lsq insert into the second
    // buffer without modifying the original matrices. The only change lsq needs
    // is the weightvalue of the augmented variable which starts at 100 and
    // being multiplied by 10 on each iteration. Hence we only pass that value
    // with "aug_weight".
    if (S->mode == 4)
    {
        S->badlin = 1;
        // Reset the RHS of the constraints to zero of the augmented system.
        for (int64_t i = 0; i < n; i++) { s[i] = 0.0; }
        S->h3 = 0.0;
        double rho = 100.0;
        S->inconsistent = 0;
        while (1)
        {
            lsq(m, S->meq, n, 1, rho, bfgs, gradx, C, d, u, v, s, mult, lsq_buffer, indices, &S->mode);
            S->h4 = 1.0 - s[n];
            if (S->mode == 4)
            {
                rho *= 10.0;
                S->inconsistent++;
                if (S->inconsistent > 5) { return; }
                continue;
            } else if (S->mode != 1) {
                return;
            }
            break;
        }
    } else if (S->mode != 1) {
        return;
    }

    // Update multipliers for L1-test
    for (int64_t i = 0; i < n; i++) { v[i] = gradx[i]; }
    NLOPT_BLAS(dgemv)("T", &m, &n, &dmone, C, &lda, mult, &one, &done, v, &one);

    S->f0 = funx;
    for (int64_t i = 0; i < n; i++) { x0[i] = sol[i]; }
    S->gs = NLOPT_BLAS(ddot)(&n, gradx, &one, s, &one);
    S->h1 = fabs(S->gs);
    S->h2 = 0.0;
    for (int64_t j = 0; j < m; j++)
    {
        if (j < S->meq)
        {
            S->h3 = d[j];
        } else {
            S->h3 = 0.0;
        }
        S->h2 = S->h2 + fmax(-d[j], S->h3);
        S->h3 = fabs(mult[j]);
        mu[j] = fmax(S->h3, (mu[j] + S->h3)/2.0);
        S->h1 = S->h1 + S->h3*fabs(d[j]);
    }

    // Check convergence
    S->mode = 0;
    if ((S->h1 < S->acc) && (S->h2 < S->acc) && (!S->badlin) && (funx == funx)) { return; }
    S->h1 = 0.0;
    for (int64_t j = 0; j < m; j++)
    {
        if (j < S->meq)
        {
            S->h3 = d[j];
        } else {
            S->h3 = 0.0;
        }
        S->h1 += mu[j]*fmax(-d[j], S->h3);
    }
    // 180
    S->t0 = funx + S->h1;
    S->h3 = S->gs - S->h1*S->h4;
    S->mode = 8;
    if (S->h3 >= 0.0) { goto RESET_BFGS; }

    // Line search with an L1 test function
    S->line = 0;
    S->alpha = 1.0;

    // Inexact line search
LINE_SEARCH:

    S->line++;
    S->h3 = (S->alpha) * (S->h3);
    NLOPT_BLAS(dscal)(&n, &S->alpha, s, &one);
    for (int64_t i = 0; i < n; i++) { sol[i] = x0[i]; }
    NLOPT_BLAS(daxpy)(&n, &done, s, &one, sol, &one);

    /* SGJ 2010: optimizing for the common case where the inexact line
       search succeeds in one step, use special mode = -2 here to
       eliminate a a subsequent unnecessary mode = -1 call, at the
       expense of extra gradient evaluations when more than one inexact
       line-search step is required */
    S->mode = ((S->eager_grad) && (S->line == 1)) ? -2 : 1;
    return;

MODE1:

    S->t = funx;
    for (int64_t j = 0; j < m; j++)
    {
        if (j < S->meq)
        {
            S->h1 = d[j];
        } else {
            S->h1 = 0.0;
        }
        S->t = S->t + mu[j]*fmax(-d[j], S->h1);
    }
    S->h1 = S->t - S->t0;

    if (nlopt_isfinite(S->h1))
    {
        if (!((S->h1 <= (S->h3 / 10.0)) || (S->line > 10)))
        {
            S->alpha = fmax(S->h3/(2.0*(S->h3 - S->h1)), alfmin);
            goto LINE_SEARCH;
        }
    } else {
        // SGJ: back off from a non-finite trial point
        S->alpha = fmax(S->alpha*0.5, alfmin);
        goto LINE_SEARCH;
    }

    // Check convergence
    S->h3 = 0.0;
    for (int64_t j = 0; j < m; j++)
    {
        if (j < S->meq)
        {
            S->h1 = d[j];
        } else {
            S->h1 = 0.0;
        }
        S->h3 = S->h3 + fmax(-d[j], S->h1);
    }
    if (
        ((fabs(funx - S->f0) < S->acc) || (NLOPT_BLAS(dnrm2)(&n, s, &one) < S->acc)) &&
        (S->h3 < S->acc) &&
        (!S->badlin) &&
        (funx == funx) // To filter for finite entries
    )
    {
        S->mode = 0;
        return;
    } else {
        S->mode = -1;
    }
    return;

LABEL255:
    // Check relaxed convergence in case of positive directional derivative
    S->h3 = 0.0;
    for (int64_t j = 0; j < m; j++)
    {
        if (j < S->meq)
        {
            S->h1 = d[j];
        } else {
            S->h1 = 0.0;
        }
        S->h3 = S->h3 + fmax(-d[j], S->h1);
    }
    if (((fabs(funx - S->f0) < S->tol) || (NLOPT_BLAS(dnrm2)(&n, s, &one) < S->tol)) &&
        (S->h3 < S->tol) &&
        (!S->badlin) &&
        (funx == funx)
    )
    {
        S->mode = 0;
    } else {
        S->mode = 8;
    }
    return;

MODEM1:

    // Call Jacobian at current x

    // Update Cholesky factors of Hessian matrix modified by BFGS formula
    // u[i] = gradx[i] - C.T @ mult - v[i]

    for (int64_t i = 0; i < n; i++) { u[i] = gradx[i]; }
    NLOPT_BLAS(dgemv)("T", &m, &n, &dmone, C, &lda, mult, &one, &done, u, &one);
    for (int64_t i = 0; i < n; i++)
    {
        u[i] = u[i] - v[i];
    }

    // L'*S
    for (int64_t i = 0; i < n; i++) { v[i] = s[i]; }
    NLOPT_BLAS(dtpmv)("L", "T", "U", &n, bfgs, v, &one);

    // D*L'*S
    j = 0;
    for (int64_t i = 0; i < n; i++) {
        v[i] = bfgs[j]*v[i];
        j += n - i;
    }

    // L*D*L'*S
    NLOPT_BLAS(dtpmv)("L", "N", "U", &n, bfgs, v, &one);

    S->h1 = NLOPT_BLAS(ddot)(&n, s, &one, u, &one);
    S->h2 = NLOPT_BLAS(ddot)(&n, s, &one, v, &one);
    S->h3 = 0.2*(S->h2);
    if (S->h1 < S->h3)
    {
        S->h4 = (S->h2 - S->h3) / (S->h2 - S->h1);
        S->h1 = S->h3;
        double tmp_dbl = 1.0 - S->h4;
        NLOPT_BLAS(dscal)(&n, &S->h4, u, &one);
        NLOPT_BLAS(daxpy)(&n, &tmp_dbl, v, &one, u, &one);
    }

    // Test for singular update, and reset hessian if so
    if ((S->h1 == 0.0) || (S->h2 == 0.0)) { goto RESET_BFGS; }

    ldl_update(n, bfgs, u, 1.0 / S->h1, v);
    ldl_update(n, bfgs, v, -1.0 / S->h2, u);

    // End of main iteration
    goto ITER_START;
}


/*
 *          min     |A*x - b|
 *        E*x = f
 *        G*x >= h
 *      xl <= x <= xu
 *
 * Problem data is kept in Lf, gradx, C, d, xl, xu arrays in a rather tedious
 * format. C(m, n) is the constraint normals, d(n) is the constraint bounds.
 * xl(n) and xu(n) are the lower and upper bounds on x.
 *
 * Lf is the LDL' factor of the BFGS matrix also holding the diagonal entries.
 *
 * NaN entries in xl, xu, signify unconstrained variables and hence not included.
 *
 * The C matrix, for a problem with all x bounds are given and finite,
 * broken into E and G as follows:
 *
 *                      ┌────┐    ┌────┐    ┌┐
 *                  meq │    │    │ E  │  = ││ f
 *                      │   ─┼────┼>   │    ││
 *                      ┼────┼    └────┘    └┘
 *                      │    │    ┌────┐    ┌┐
 *                      │    │    │    │    ││
 *      mineq = m - meq │   ─┼────┼>   │    ││
 *                      │    │    │    │    ││
 *                      │    │    │    │    ││
 *                      └────┘    │    │ >= ││
 *                        C       ┼────┼    ┼┼
 *                              n │  I │    ││  xl
 *                                ┼────┼    ┼┼
 *                              n │ -I │    ││ -xu
 *                                └────┘    └┘
 *                                   G       h
 *
 * A and b are stored in Lf[] in LAPACK packed format where Lf holds a unit, lower
 * triangular matrix with diagonal entries are overwritten by the entries of d[]
 * and vector and gradx[].
 *
 *  Lf[] = [d[0], s[1], s[2], . , d[1], s[n + 2], d[2], ...]
 *
 *  interpreted as:
 *
 *         [d[ 0 ],                          ]
 *         [s[ 1 ], d[ 1 ], .  ,             ]
 *  Lf[] = [s[ 2 ], s[n+2], .  ,             ]
 *         [ .    ,   .   , .  , d[n-1]      ]
 *         [s[ n ], s[2*n], .  ,   .   , d[n]]
 *
 * Then, the following relations recover the A and b
 *
 *          A = sqrt(d[]) * Lf[]^T
 *          b = - inv( Lf[] * sqrt(d[]) ) * gradx[]
 *
 * The solution is returned in x() and the Lagrange multipliers are returned in y().
 *
 * For solving the problem in case of a detection of inconsistent linearization,
 * see D. Kraft, "A software package for Sequential Quadratic Programming"
 * Section 2.2.3
 *
 * In the original code, the augmented system is detected by mismatch of certain
 * integers which is making things quite unreadable. Here we explicitly pass a
 * flag.
 *
 * Inconsistent linearization augments all arrays to accommodate for the dummy
 * variable. The function is still called with the original sizes but the flag
 * allows for enlarging the problem and hence the supplied buffer should accommodate
 * for this extra space.
 *
 * The required buffer size is given by:
 * (2*(m - meq)*(n + 1)+2)*(n - meq +1) + 2*2*(m-meq)*(n + 1) + 2*(m-meq)*(n + 1)
 *  + 2*meq + ld + (ld + 2*(m-meq)*(n + 1))*(n - meq)
 *
 */
static void
lsq(
    blas_int m, blas_int meq, blas_int n, blas_int augment, double aug_weight, const double* restrict Lf,
    const double* restrict gradx, const double* restrict C, const double* restrict d,
    double* restrict xl, double* restrict xu, double* restrict x,
    double* restrict y, double* buffer, blas_int* jw, int64_t* mode)
{
    blas_int one = 1, orign = n;
    blas_int mineq = m - meq;
    double xnorm = 0.0;
    blas_int cursor = 0;
    blas_int ld = n;
    blas_int n_wG_rows = 0;

    if (augment) {
        ld = n + 1;
        x[n]     = 1.0;
        xl[n]    = 0.0;
        xu[n]    = 1.0;
    }

    // Recover A and b from Lf and gradx
    for (int64_t i = 0; i < (ld+2)*ld; i++) { buffer[i] = 0.0; }
    double* restrict wA = buffer;
    double* restrict wb = &buffer[ld*(ld+1)];

    // Depending on augmented, wA is either the full array or the top-left block.

    for (int64_t j = 0; j < n; j++)
    {
        double diag = sqrt(Lf[cursor++]);      // Extract the diagonal value from Lf.
        wA[j + j * ld] = diag;                 // Place the sqrt diagonal.
        for (int64_t i = j + 1; i < n; i++)
        {
            wA[j + i * ld] = Lf[cursor++] * diag;
        }
    }

    // Compute b = - 1/sqrt(d[]) * inv(Lf[]) * gradx[]. Lf is already in packed format.
    for (int64_t i = 0; i < n; i++) { wb[i] = gradx[i]; }
    NLOPT_BLAS(dtpsv)("L", "N", "U", &n, Lf, wb, &one);
    cursor = 0;
    for (int64_t i = 0; i < n; i++)
    {
        wb[i] /= -sqrt(Lf[cursor]);
        cursor += n - i;
    }

    // If augmented, fill in the extra entry in the bottom right corner.
    if (augment) { wA[ld*ld - 1] = aug_weight; }

    // If augmented, also increase the number of variables by 1.
    if (augment) { n++; }

    // Get the equality constraints if given.
    double* restrict wE = &buffer[n*(n+1) + n];
    double* restrict wf = &buffer[n*(n+1) + n + n*meq];
    if (meq > 0)
    {
        for (int64_t j = 0; j < n-1; j++)
        {
            for (int64_t i = 0; i < meq; i++)
            {
                wE[i + j*meq] = C[i + j*m];
            }
        }
        if (augment)
        {
            // n is incremented hence all Ceq is now in wE. Add the extra column.
            for (int64_t i = 0; i < meq; i++) { wE[i + (n-1)*meq] = -d[i]; }

        } else  {
            // If not augmented then handle j = n - 1 that is skipped.
            for (int64_t i = 0; i < meq; i++) { wE[i + (n-1)*meq] = C[i + (n-1)*m]; }

        }
        for (int64_t i = 0; i < meq; i++) { wf[i] = -d[i]; }
    }

    // Get the inequality constraints if given. First zero out wG and wh.
    double* restrict wG = &buffer[n*(n+1) + n + n*meq + meq];
    double* restrict wh = &buffer[n*(n+1) + n + n*meq + meq + (mineq + 2*n)*ld];
    // Zero out wG and wh
    for (int64_t i = 0; i < (mineq + 2*n)*(ld + 1); i++) { wG[i] = 0.0; }

    // Convert the bounds on x to +I and -I blocks in G.
    // Augment h by xl and -xu.
    // Unbounded constraints are signified by NaN values and they do not appear
    // in G and h. Hence there is a nancount tab to keep track of them.

    // We first populate "wh" to get the number of unbounded constraints. That will
    // define the unskipped row number of wG. This is different than the original
    // Fortran code where the max allocated row number and the actual row number
    // of wG has been kept separate and it causes to be sent to every nested
    // function call. Instead we form wG and wh once with fixed size.

    blas_int nancount = 0;
    blas_int nrow = mineq;
    if (m > meq)
    {
        for (int64_t i = 0; i < mineq; i++) { wh[i] = -d[meq + i]; }
    }
    for (int64_t i = 0; i < n; i++)
    {
        if (isnan(xl[i]))
        {
            nancount++;
        } else {
            wh[nrow++] = xl[i];
        }
    }
    for (int64_t i = 0; i < n; i++)
    {
        if (isnan(xu[i]))
        {
            nancount++;
        } else {
            wh[nrow++] = -xu[i];
        }
    }

    n_wG_rows = mineq + 2*n - nancount;

    // Now that we know the actual row number of wG, we can finally populate
    // the top part with C.
    if (m > meq)
    {
        for (int64_t j = 0; j < orign; j++)
        {
            for (int64_t i = 0; i < mineq; i++)
            {
                wG[i + j*n_wG_rows] = C[meq + i + j*m];
            }
        }
    }

    // If augmented add the extra column.
    if (augment)
    {
        for (int64_t i = 0; i < mineq; i++)
        {
            wG[i + orign*n_wG_rows] = fmax(-d[meq + i], 0.0);
        }
    }

    // Reset counter
    nrow = mineq;
    for (int64_t i = 0; i < n; i++)
    {
        if (!isnan(xl[i]))
        {
            wG[nrow + i*n_wG_rows] = 1.0;
            nrow++;
        }
    }
    for (int64_t i = 0; i < n; i++)
    {
        if (!isnan(xu[i]))
        {
            wG[nrow + i*n_wG_rows] = -1.0;
            nrow++;
        }
    }

    // Assign the remaining part of the buffer to the LSEI problem.
    double* restrict lsei_scratch = &wh[mineq + 2*n];

    lsei(ld, meq, n_wG_rows, n, wA, wb, wE, wf, wG, wh, x, lsei_scratch, jw, &xnorm, mode);

    if (*mode == 1)
    {
        // Restore the Lagrange multipliers, first equality, then inequality.
        for (int64_t i = 0; i < meq; i++) { y[i] = lsei_scratch[i+n_wG_rows]; }
        for (int64_t i = 0; i < mineq; i++) { y[meq + i] = lsei_scratch[i]; }

        // Set the user-defined bounds on x to NaN
        for (int64_t i = 0; i < 2*n; i++) { y[m + i] = NAN; }
    }

    // Clamp the solution, if given, to the finite bound interval
    for (int64_t i = 0; i < n; i++)
    {
        if ((!isnan(xl[i])) && (x[i] < xl[i])) { x[i] = xl[i]; }
        else if ((!isnan(xu[i])) && (x[i] > xu[i])) { x[i] = xu[i]; }
    }

    return;
}


/*
 * Least squares problem min |A*x - b| for an A that is upper triangular (zero
 * below the diagonal), as the LSQ matrix is when there are no equality
 * constraints.  DGELSY with all columns fixed leaves such an A unchanged in its
 * QR factorization (every reflector is the identity), but still spends O(m*n^2)
 * on it, and LAPACK's blocked QR does not skip the zeros.  This does the
 * remaining steps of DGELSY directly, with the same rank estimation and the same
 * triangular solve, so the result is the same.  Returns 0 without modifying b
 * if DGELSY would rescale A or b, or if A is rank deficient (for the
 * minimum-norm solution); the caller then calls DGELSY.
 *
 * A is (m x n) with leading dimension m, b is (m), work is at least 2*n.  Only
 * the upper triangle of A is read (the rest is checked in debug builds).
 */
static int
lsq_triangular(blas_int m, blas_int n, const double* restrict a, double* restrict b,
               double rcond, blas_int* rank, double* restrict work)
{
    // DGELSY's scaling thresholds: DLAMCH('S')/DLAMCH('P') and its inverse.
    const double smlnum = 2.2250738585072014e-308/2.220446049250313e-16;
    const double bignum = 1.0/smlnum;
    blas_int one = 1, imax = 1, imin = 2;
    double anrm = 0.0, bnrm = 0.0, smax, smin, done = 1.0;
    double* restrict xmin = work;
    double* restrict xmax = &work[n];

    if ((m < n) || (n == 0)) { return 0; }
#ifndef NDEBUG
    for (blas_int j = 0; j < n; j++)
    {
        for (blas_int i = j + 1; i < m; i++) { assert(a[i + j*m] == 0.0); }
    }
#endif
    for (blas_int j = 0; j < n; j++)
    {
        for (blas_int i = 0; i <= j; i++)
        {
            double v = fabs(a[i + j*m]);
            if (!(v <= bignum)) { return 0; }
            if (v > anrm) { anrm = v; }
        }
    }
    if (anrm < smlnum) { return 0; }
    for (blas_int i = 0; i < m; i++)
    {
        double v = fabs(b[i]);
        if (!(v <= bignum)) { return 0; }
        if (v > bnrm) { bnrm = v; }
    }
    if ((bnrm > 0.0) && (bnrm < smlnum)) { return 0; }

    // Incremental condition estimation, as in DGELSY.
    smax = fabs(a[0]);
    smin = smax;
    if (smax == 0.0) { return 0; }
    xmin[0] = 1.0;
    xmax[0] = 1.0;
    *rank = 1;
    while (*rank < n)
    {
        blas_int i = *rank;
        double sminpr, smaxpr, s1, c1, s2, c2;
        NLOPT_BLAS(dlaic1)(&imin, &i, xmin, &smin, &a[i*m], &a[i + i*m], &sminpr, &s1, &c1);
        NLOPT_BLAS(dlaic1)(&imax, &i, xmax, &smax, &a[i*m], &a[i + i*m], &smaxpr, &s2, &c2);
        if (!(smaxpr*rcond <= sminpr)) { return 0; }
        for (blas_int k = 0; k < i; k++)
        {
            xmin[k] *= s1;
            xmax[k] *= s2;
        }
        xmin[i] = c1;
        xmax[i] = c2;
        smin = sminpr;
        smax = smaxpr;
        (*rank)++;
    }

    // Q^T*b = b, solve R*x = b.
    NLOPT_BLAS(dtrsm)("L", "U", "N", "N", &n, &one, &done, a, &m, b, &m);
    return 1;
}


/*
 * Solve equality and inequality constrained least squares problem (LSEI)
 *      min |A*x - b|, subject to E*x = f, G*x >= h.
 *
 *  ma, me, mg : number of rows in A, E, G
 *  n          : number of columns in A, x
 *  a          : matrix A (ma x n)
 *  b          : vector b (ma)
 *  e          : matrix E (me x n)
 *  f          : vector f (me)
 *  g          : matrix G (mg x n)
 *  h          : vector h (mg)
 *  x          : solution vector x (n)
 *  buffer     : work buffer (mg + 2)*(n - me +1) + 3*mg + 2*me + ma + (ma + mg)*(n - me)
 *  jw         : integer work array
 *  xnorm      : norm of the solution
 *  mode       : return code
 *
 *  The buffer pointers that will be used:
 *  buffer[0]                : Lagrange multipliers (mg + me)
 *  buffer[mg + me]          : wb, Modified b vector (ma)
 *  buffer[mg + me + ma]     : tau, Pivots for the RQ decomposition of E (me)
 *  buffer[mg + 2*me + ma]   : Scratch space
 *
 */
static void
lsei(blas_int ma, blas_int me, blas_int mg, blas_int n,
     double* restrict a, double* restrict b, double* restrict e,
     double* restrict f, double* restrict g, double* restrict h,
     double* restrict x, double* restrict buffer, blas_int* jw,
     double* xnorm, int64_t* mode)
{
    blas_int one = 1, nvars = 0, info = 0, lde = 0, ldg = 0;
    double done = 1.0, dmone = -1.0, dzero = 0.0, t= 0.0;
    const double epsmach = 2.220446049250313e-16;

    for (int64_t i = 0; i < n; i++) { x[i] = 0.0; }
    // Return if the problem is over-constrained.
    if (me > n) { *mode = 2; return; }

    //    [E]         [E2 |  R]                                [x ]
    //    [A] @ Q.T = [A2 | A1]  ,and, x is partitioned as x = [--]
    //    [G]         [G2 | G1]                                [xe]

    // me = 0 skips the equality constraint related computations even though it
    // causes aliasing below. The aliased arrays are not referenced in that case.
    // Use at least 1 for the leading dimension of E even when me = 0 for LAPACK
    // calls.
    nvars = (n - me);
    double* restrict gmults      = &buffer[0];
    double* restrict emults      = &buffer[mg];
    double* restrict wb          = &buffer[me + mg];
    double* restrict tau         = &buffer[me + mg + ma];
    double* restrict a2          = &buffer[mg + 2*me + ma];
    double* restrict g2          = &buffer[mg + 2*me + ma + ma*nvars];
    double* restrict lsi_scratch = &buffer[mg + 2*me + ma + (ma + mg)*nvars];

    // RQ decomposition of equality constraint data E and application to A, G.
    // LAPACK RQ routine dgerq2 forms R on the right.
    // dgeqr2 is the unblocked versions of dgeqrf without the memory allocation.
    // Use top of the yet unutilized scratch space for throw-away work.
    lde = (me > 0 ? me : 1);
    ldg = (mg > 0 ? mg : 1);
    NLOPT_BLAS(dgerq2)(&me, &n, e, &lde, tau, lsi_scratch, &info);

    // Right triangularize E and apply Q.T to A and G from the right.
    NLOPT_BLAS(dormr2)("R", "T", &ma, &n, &me, e, &lde, tau, a, &ma, lsi_scratch, &info);
    NLOPT_BLAS(dormr2)("R", "T", &mg, &n, &me, e, &lde, tau, g, &ldg, lsi_scratch, &info);

    // Check the diagonal elements of E for rank deficiency.
    for (int64_t i = 0; i < me; i++)
    {
        if (!(fabs(e[i + (nvars + i)*me]) >= epsmach)) { *mode = 6;return; }
    }
    // Solve E*x = f and modify b.
    // Note: RQ forms R at the right of E instead of [0, 0] position.
    for (int64_t i = 0; i < me; i++) { x[nvars + i] = f[i]; }
    NLOPT_BLAS(dtrsv)("U", "N", "N", &me, &e[(nvars)*me], &lde, &x[nvars], &one);

    *mode = 1;
    // Zero out the inequality multiplier.
    for (int64_t i = 0; i < mg; i++) { gmults[i] = 0.0; }

    // If the problem is fully equality-constrained, revert the basis and return.
    if (me == n) { goto ORIGINAL_BASIS; }

    // Compute the modified RHS wb = b - A1*x
    // Copy b into wb
    for (int64_t i = 0; i < ma; i++) { wb[i] = b[i]; }
    // Compute wb -= A1*xe
    NLOPT_BLAS(dgemv)("N", &ma, &me, &dmone, &a[ma*nvars], &ma, &x[nvars], &one, &done, wb, &one);

    // Store the transformed A2 and G2 in the buffer
    for (int64_t j = 0; j < nvars; j++)
    {
        for (int64_t i = 0; i < ma; i++)
        {
            a2[i + j*ma] = a[i + j*ma];
        }
        for (int64_t i = 0; i < mg; i++)
        {
            g2[i + j*mg] = g[i + j*mg];
        }
    }

    if (mg == 0)
    {
        // No inequality constraints, solve the least squares problem directly.
        // We deliberately use the unblocked algorithm to avoid allocation.
        blas_int lwork = ma*nvars + 3*nvars + 1;
        // Save the RHS for residual computation
        double* restrict wb_orig = &lsi_scratch[lwork];
        for (int64_t i = 0; i < ma; i++) { wb_orig[i] = wb[i]; }

        blas_int krank = 0;
        t = sqrt(epsmach);
        // Without equality constraints A is still the upper triangular LSQ
        // matrix (built by lsq; the RQ transformation above is empty).
        if ((me == 0) && lsq_triangular(ma, nvars, a2, wb, t, &krank, lsi_scratch))
        {
            // Same result as dgelsy below, including the permutation.
            for (int64_t i = 0; i < nvars; i++) { jw[i] = i + 1; }
        } else {
            // Mark all columns as fixed so that dgelsy uses the (blocked) QR
            // without column pivoting, which is about twice as fast here. scipy
            // gets the same behavior implicitly, except on the first call, since
            // jw still holds the permutation from the previous call.
            for (int64_t i = 0; i < nvars; i++) { jw[i] = 1; }
            NLOPT_BLAS(dgelsy)(&ma, &nvars, &one, a2, &ma, wb, &ma, jw, &t, &krank, lsi_scratch, &lwork, &info);
        }

        // Copy the solution to x
        for (int64_t i = 0; i < nvars; i++) { x[i] = wb[i]; }

        // Compute the residual and its norm, use a since a2 is overwritten.
        NLOPT_BLAS(dgemv)("N", &ma, &nvars, &done, a, &ma, x, &one, &dmone, wb_orig, &one);
        *xnorm = NLOPT_BLAS(dnrm2)(&ma, wb_orig, &one);

        *mode = 7;
        if (krank < nvars) { return; }
        *mode = 1;
        goto ORIGINAL_BASIS;
    }

    // Modify h, and solve the inequality constrained least squares problem.
    // h -= G1*xe
    NLOPT_BLAS(dgemv)("N", &mg, &me, &dmone, &g[mg*nvars], &ldg, &x[nvars], &one, &done, h, &one);

    lsi(ma, mg, nvars, a2, wb, g2, h, x, lsi_scratch, jw, xnorm, mode);

    // Copy multipliers from scratch to gmults
    for (int64_t i = 0; i < mg; i++) { gmults[i] = lsi_scratch[i]; }

    // If no equality constraints this was an LSI problem all along.
    if (me == 0) { return; }

    t = NLOPT_BLAS(dnrm2)(&me, &x[nvars], &one);
    // Modify the norm by adding the equality solution.
    *xnorm = hypot(*xnorm, t);
    if (*mode != 1) { return; }

ORIGINAL_BASIS:
    // Convert the solution and multipliers to the original basis.
    // b = A*x - b (residuals)
    NLOPT_BLAS(dgemv)("N", &ma, &n, &done, a, &ma, x, &one, &dmone, b, &one);
    // f = A1^T*b - G1^T*w
    NLOPT_BLAS(dgemv)("T", &ma, &me, &done, &a[nvars*ma], &ma, b, &one, &dzero, f, &one);
    NLOPT_BLAS(dgemv)("T", &mg, &me, &dmone, &g[nvars*mg], &ldg, gmults, &one, &done, f, &one);

    // x = Q.T*x
    NLOPT_BLAS(dormr2)("L", "T", &n, &one, &me, e, &lde, tau, x, &n, lsi_scratch, &info);

    // Solve the triangular system for the equality multipliers, emults.
    for (int64_t i = 0; i < me; i++) { emults[i] = f[i]; }
    NLOPT_BLAS(dtrsv)("U", "T", "N", &me, &e[(n - me)*me], &lde, emults, &one);

    return;
}


/*
 * Solve inequality constrained least squares problem
 *      min |Ax - b|  subject to Gx >= h
 *
 * A is (ma x n), b is (ma), G is (mg x n), h is (mg), x is (n)
 * buffer is at least (mg+2)*(n+1) + 3*mg
 * jw is at least (mg)
 * xnorm is the 2-norm of the residual vector
 * mode is the integer return code
 *
 * Return codes for mode
 *  1: successful computation
 *  2: error return because of wrong dimensions
 *  3: iteration count exceeded by nnls
 *  4: inequality constraints incompatible
 *  5: matrix A is not rank n
 *
*/
static void
lsi(blas_int ma, blas_int mg, blas_int n, double* restrict a, double* restrict b, double* restrict g,
    double* restrict h, double* restrict x, double* restrict buffer, blas_int* jw,
    double* xnorm, int64_t* mode)
{
    blas_int one = 1, tmp_int = 0, info = 0;
    double done = 1.0, dmone = -1.0, tmp_dbl = 0.0;
    const double epsmach = 2.220446049250313e-16;

    // QR decomposition of A and application to b.
    // We use the unblocked versions of the LAPACK routines to avoid
    // allocating extra "work" memory for the blocked versions.
    tmp_int = (ma < n ? ma : n);
    NLOPT_BLAS(dgeqr2)(&ma, &n, a, &ma, buffer, &buffer[tmp_int], &info);

    // Compute Q^T b
    NLOPT_BLAS(dorm2r)("L", "T", &ma, &one, &tmp_int, a, &ma, buffer, b, &ma, &buffer[tmp_int], &info);

    // Check the diagonal elements of R for rank deficiency.
    *mode = 5;
    *xnorm = 0.0;
    for (int64_t i = 0; i < tmp_int; i++) {
        if (!(fabs(a[i + i*ma]) >= epsmach)) { return; }
    }
    // Transform G and h to form the LDP problem.
    // Solve XR = G where R is the upper triangular matrix from the QR.
    // The result is stored in G.
    // Note: There is an inherent assumption that ma >= n. This is a bug carried
    // over here from the original slsqp implementation.
    NLOPT_BLAS(dtrsm)("R", "U", "N", "N", &mg, &n, &done, a, &ma, g, &mg);
    // h = h - Xf
    NLOPT_BLAS(dgemv)("N", &mg, &n, &dmone, g, &mg, b, &one, &done, h, &one);

    // Solve the LDP problem.
    ldp(mg, n, g, h, x, buffer, jw, xnorm, mode);
    if (*mode != 1) { return; }

    // Convert to the solution of the original problem.
    NLOPT_BLAS(daxpy)(&n, &done, b, &one, x, &one);
    NLOPT_BLAS(dtrsv)("U", "N", "N", &n, a, &ma, x, &one);

    // If any, compute the norm of the tail of b and add to xnorm
    tmp_int = ma - n;
    tmp_dbl = NLOPT_BLAS(dnrm2)(&tmp_int, &b[(n + 1 > ma ? ma : n + 1) - 1], &one);
    *xnorm = hypot(*xnorm, tmp_dbl);

    return;
}

/*
 * Solve least distance problem
 *  min (1/2)|x|^2  subject to  Gx >= h
 *
 * G is (m x n), h is (m)
 * buffer is at least (m+2)*(n+1) + 3*m
 * indices is int(n)
 * x is (n)
 * xnorm is the norm of the solution if succeeded
 * mode is the return code integer
 *
 * Mode return values
 *  1  : solution found
 *  2  : bad input dimensions
 *  3  : iteration count exceeded by nnls
 *  4  : inequality constraints incompatible
 *
*/
static void
ldp(blas_int m, blas_int n, double* restrict g, double* restrict h, double* restrict x,
    double* restrict buffer, blas_int* indices, double* xnorm, int64_t* mode)
{
    blas_int one = 1;
    double dzero = 0.0, rnorm = 0.0;
    // Check for inputs and initialize x
    if (n <= 0) { *mode = 2; return; }
    for (int64_t i = 0; i < n; i++) { x[i] = 0.0; }
    if (m == 0) { *mode = 1; return; }

    // Define pointers for the variables on buffer
    double* restrict a    = &buffer[0];
    double* restrict b    = &buffer[m*(n+1)];
    double* restrict zz   = &buffer[(m+1)*(n+1)];
    double* restrict y    = &buffer[(m+2)*(n+1)];
    double* restrict w    = &buffer[(m+2)*(n+1) + m];
    double* restrict work = &buffer[(m+2)*(n+1) + 2*m];

    // Save the dual problem data into buffer
    //       dual problem [G^T] [x] = [0]
    //                    [h^T]       [1]

    // LHS, G is (m x n), h is (m). Both transposed and stacked into (n+1) x m.
    for (int64_t j = 0; j < m; j++)
    {
        for (int64_t i = 0; i < n; i++)
        {
            a[i + j*(n+1)] = g[j + i*m];
        }
        // Place h in the last row.
        a[n + j*(n+1)] = h[j];
    }
    // RHS is (n+1)
    for (int64_t i = 0; i < n; i++) { b[i] = 0.0; }
    b[n] = 1.0;

    // Solve the dual problem
    nnls(n+1, m, a, b, y, w, zz, work, indices, 3*m, &rnorm, mode);
    if (*mode != 1) { return; }
    *mode = 4;
    if (rnorm <= 0.0) { return; }

    // Solve the primal problem
    double fac = 1.0 - NLOPT_BLAS(ddot)(&m, h, &one, y, &one);
    if (!((1.0 + fac) - 1.0 > 0.0)) { return; }
    *mode = 1;
    fac = 1.0 / fac;
    NLOPT_BLAS(dgemv)("T", &m, &n, &fac, g, &m, y, &one, &dzero, x, &one);
    *xnorm = NLOPT_BLAS(dnrm2)(&n, x, &one);

    // Compute the lagrange multipliers for the primal problem
    for (int64_t i = 0; i < m; i++) { buffer[i] = fac*y[i]; }
    return;
}


/*
 *
 * Updates the LDL' factors of matrix a by rank-one matrix sigma*z*z'
 * n     : order of the coefficient matrix a
 * a     : positive definite matrix of dimension n; only the lower triangle is
 *         used and is stored column by column as one dimensional array of
 *         dimension n*(n+1)/2.
 * z     : vector of dimension n of updating elements
 * sigma : scalar factor by which the modifying dyade z*z' is multiplied
 * w     : working array of dimension n
 *
 * Uses the composite-t method of fletcher and powell as described in "On the
 * modification of LDL' factorizations", DOI:10.1090/S0025-5718-1974-0359297-1
 *
 * Implemented by: Dieter Kraft, dfvlr - Institut für Dynamik der Flugsysteme
 *                 D-8031  Oberpfaffenhofen
 *
 */
static void
ldl_update(blas_int n, double* restrict a, double* restrict z, double sigma, double* restrict w)
{
    blas_int j, ij = 0;
    const double epsmach = 2.220446049250313e-16;
    if (sigma == 0.0) { return; }
    double alpha, beta, delta, gamma, u, v, tp, t = 1.0 / sigma;

    if (sigma <= 0.0)
    {
        // Negative update
        for (int64_t i = 0; i < n; i++) { w[i] = z[i]; }
        for (int64_t i = 0; i < n; i++)
        {
            v = w[i];
            t = t + v*v/a[ij];
            for (int64_t j = i + 1; j < n; j++)
            {
                ij++;
                w[j] = w[j] - v*a[ij];
            }
            ij++;
        }
        if (t >= 0.0) { t = epsmach / sigma; }

        for (int64_t i = 0; i < n; i++)
        {
            j = n - i - 1;
            ij -= i + 1;
            u = w[j];
            w[j] = t;
            t = t - u*u / a[ij];
        }
    }

    // Positive update
    for (int64_t i = 0; i < n; i++)
    {
        v = z[i];
        delta = v / a[ij];
        // sigma == 0.0 is handled at the beginning.
        tp = (sigma < 0.0 ? w[i] : t + delta*v);
        alpha = tp / t;
        a[ij] = alpha*a[ij];
        if (i == n - 1) { return; }
        beta = delta / tp;
        if (alpha <= 4.0)
        {
            for (int64_t j = i + 1; j < n; j++)
            {
                ij++;
                z[j] = z[j] - v * a[ij];
                a[ij] = a[ij] + beta * z[j];
            }
        } else {
            gamma = t / tp;
            for (int64_t j = i + 1; j < n; j++)
            {
                ij++;
                u = a[ij];
                a[ij] = gamma * u + beta * z[j];
                z[j] = z[j] - v * u;
            }
        }
        ij++;
        t = tp;
    }

    return;
}

#define MAX2(a,b) ((a) >= (b) ? (a) : (b))

/* Worst-case size of the slsqp_body buffer, summed over
     n*(n+1)//2 + m + 4*n + 3                                           # SLSQP
     (n+1)*(n+2) + (n+1)*meq + m + (mineq + 2*n + 2)*(n+1) +  3*n + 3   # LSQ
     mineq + 2n + 2 + 2*meq + (n+1) + (mineq + 3n + 3)*(n + 1 - meq)    # LSEI
     (mineq + 2n + 2 + 2)*(n + 2) + mineq + 2n + 2                      # LDP
     2*(mineq + 2n + 2)                                                 # NNLS
   (from scipy/optimize/_slsqp_py.py, plus the NNLS work array for the
   blocked Householder application) */
static size_t slsqp_buffer_size(int64_t m, int64_t meq, int64_t n)
{
     int64_t size = n*(n+1)/2 + 3*m*n - (m + 5*n + 7)*meq + 9*m + 8*n*n + 35*n + meq*meq + 28
	  + (m - meq) + 2*n + 2;
     /* if no inequality constraints are given, top up workspace for
        the dgelsy path in lsei */
     if (m == meq)
	  size += 2*n*(n + 1);
     return (size_t) (size > 1 ? size : 1);
}

nlopt_result nlopt_slsqp(unsigned n, nlopt_func f, void *f_data,
			 unsigned m, nlopt_constraint *fc,
			 unsigned p, nlopt_constraint *h,
			 const double *lb, const double *ub,
			 double *x, double *minf,
			 nlopt_stopping *stop, int eager_grad, double acc)
{
     slsqp_state state;
     unsigned mtot = nlopt_count_constraints(m, fc);
     unsigned ptot = nlopt_count_constraints(p, h);
     double *work, *cgrad, *c, *grad, *xl, *xu, *mult, *buffer,
	  fcur, *xcur, fprev, *xprev, *cgradtmp;
     blas_int *indices;
     int mpi = (int) (mtot + ptot), pi = (int) ptot, mpi1 = mpi > 0 ? mpi : 1;
     size_t len_buffer, len_mult;
     int mode = 0, prev_mode = 0;
     unsigned i, ii;
     nlopt_result ret = NLOPT_SUCCESS;
     int feasible, feasible_cur;
     double infeasibility = HUGE_VAL, infeasibility_cur = HUGE_VAL;
     unsigned max_cdim;
     int want_grad = 1;

     if (ptot > n) {
       nlopt_stop_msg(stop, "slsqp: more equality constraints than variables");
       return NLOPT_INVALID_ARGS;
     }

     max_cdim = MAX2(nlopt_max_constraint_dim(m, fc),
		    nlopt_max_constraint_dim(p, h));
     len_buffer = slsqp_buffer_size(mpi, pi, n);
#define U(n) ((unsigned) (n))
     len_mult = U(mpi) + 2*n + 2;

     work = (double *) malloc(sizeof(double) * (U(mpi1) * n
						+ U(mpi1)
						+ n + n + n + max_cdim*n
						+ n + n + len_mult
						+ len_buffer));
     if (!work) return NLOPT_OUT_OF_MEMORY;
     /* indices doubles as the dgelsy pivot array, must start zeroed */
     indices = (blas_int *) calloc(len_mult, sizeof(blas_int));
     if (!indices) { free(work); return NLOPT_OUT_OF_MEMORY; }
     cgrad = work;
     c = cgrad + U(mpi1) * n;
     grad = c + mpi1;
     xcur = grad + n;
     xprev = xcur + n;
     cgradtmp = xprev + n;
     xl = cgradtmp + max_cdim*n;
     xu = xl + n;
     mult = xu + n;
     buffer = mult + len_mult;
     memset(mult, 0, sizeof(double) * len_mult);
     memset(buffer, 0, sizeof(double) * len_buffer);

     /* the subproblem solver marks unbounded variables with NaN bounds */
     for (i = 0; i < n; ++i) {
	  xl[i] = nlopt_isinf(lb[i]) ? NAN : lb[i];
	  xu[i] = nlopt_isinf(ub[i]) ? NAN : ub[i];
     }

     memset(&state, 0, sizeof(state));
     /* acc = 0 disables the convergence tests of the original code, so that
        only NLopt's own tests below apply */
     state.acc = fabs(acc);
     state.itermax = 0; /* tell slsqp to ignore this check, since we check evaluation counts ourselves */
     state.m = mpi;
     state.meq = pi;
     state.n = n;
     state.mode = 0;
     state.eager_grad = eager_grad;

     memcpy(xcur, x, sizeof(double) * n);
     memcpy(xprev, x, sizeof(double) * n);
     fprev = fcur = *minf = HUGE_VAL;
     feasible = feasible_cur = 0;

     goto eval_f_and_grad; /* eval before calling slsqp the first time */

     do {
	  slsqp_body(&state, fcur, grad, cgrad, c, xcur, mult, xl, xu,
		     buffer, indices);
	  mode = (int) state.mode;

	  if (mode == 1 || mode == -1 || mode == -2) {
	       /* SGJ 2010: ensure roundoff doesn't push us past bound constraints */
	       for (i = 0; i < n; ++i) {
		    if (xcur[i] < lb[i]) xcur[i] = lb[i];
		    else if (xcur[i] > ub[i]) xcur[i] = ub[i];
	       }
	  }

	  switch (mode) {
	  case -1:  /* objective & gradient evaluation */
	      if (prev_mode == -2 && !want_grad) break; /* just evaluated this point */
	      /* fall through */
	  case -2:
	      eval_f_and_grad:
	      want_grad = 1;
	      /* fall through */
	  case 1:{ /* don't need grad unless we don't have it yet */
	      double *newgrad = 0;
	      double *newcgrad = 0;
	      if (want_grad) {
		  newgrad = grad;
		  newcgrad = cgradtmp;
	      }
	      feasible_cur = 1; infeasibility_cur = 0;
	      fcur = f(n, xcur, newgrad, f_data);
	      ++ *(stop->nevals_p);
	      if (nlopt_stop_forced(stop)) {
		  fcur = HUGE_VAL; ret = NLOPT_FORCED_STOP; goto done; }
	      if (nlopt_isfinite(fcur)) {
		  want_grad = 0;
		  ii = 0;
		  for (i = 0; i < p; ++i) {
		      unsigned j, k;
		      nlopt_eval_constraint(c+ii, newcgrad, h+i, n, xcur);
		      if (nlopt_stop_forced(stop)) {
			  ret = NLOPT_FORCED_STOP; goto done; }
		      for (k = 0; k < h[i].m; ++k, ++ii) {
			  infeasibility_cur =
			      MAX2(infeasibility_cur, fabs(c[ii]));
			  feasible_cur =
			      feasible_cur && fabs(c[ii]) <= h[i].tol[k];
			  if (newcgrad) {
			      for (j = 0; j < n; ++ j)
				  cgrad[j*U(mpi1) + ii] = cgradtmp[k*n + j];
			  }
		      }
		  }
		  for (i = 0; i < m; ++i) {
		      unsigned j, k;
		      nlopt_eval_constraint(c+ii, newcgrad, fc+i, n, xcur);
		      if (nlopt_stop_forced(stop)) {
			  ret = NLOPT_FORCED_STOP; goto done; }
		      for (k = 0; k < fc[i].m; ++k, ++ii) {
			  infeasibility_cur =
			      MAX2(infeasibility_cur, c[ii]);
			  feasible_cur =
			      feasible_cur && c[ii] <= fc[i].tol[k];
			  if (newcgrad) {
			      for (j = 0; j < n; ++ j)
				  cgrad[j*U(mpi1) + ii] = -cgradtmp[k*n + j];
			  }
			  c[ii] = -c[ii]; /* slsqp sign convention */
		      }
		  }
	      }
	      break;}
	      case 0: /* required accuracy for solution obtained */
		  goto done;
	      case 8: /* positive directional derivative for linesearch */
		  /* relaxed convergence check for a feasible_cur point,
		     as in the SLSQP code (except xtol as well as ftol) */
		  ret = NLOPT_ROUNDOFF_LIMITED; /* usually why deriv>0 */
		  if (feasible_cur) {
		      double save_ftol_rel = stop->ftol_rel;
		      double save_xtol_rel = stop->xtol_rel;
		      double save_ftol_abs = stop->ftol_abs;
		      stop->ftol_rel *= 10;
		      stop->ftol_abs *= 10;
		      stop->xtol_rel *= 10;
		      if (nlopt_stop_ftol(stop, fcur, state.f0))
			  ret = NLOPT_FTOL_REACHED;
		      else if (nlopt_stop_x(stop, xcur, buffer + n*(n+1)/2 /* x0 */))
			  ret = NLOPT_XTOL_REACHED;
		      stop->ftol_rel = save_ftol_rel;
		      stop->ftol_abs = save_ftol_abs;
		      stop->xtol_rel = save_xtol_rel;
		  }
		  goto done;
	      case 5: /* singular matrix E in LSQ subproblem */
	      case 6: /* singular matrix C in LSQ subproblem */
	      case 7: /* rank-deficient equality constraint subproblem */
		  ret = NLOPT_ROUNDOFF_LIMITED;
		  goto done;
	      case 4: /* inequality constraints incompatible */
	      case 3: /* more than 3*n iterations in LSQ subproblem */
	      case 9: /* more than iter iterations in SQP */
                  nlopt_stop_msg(stop, "bug: more than iter SQP iterations");
		  ret = NLOPT_FAILURE;
		  goto done;
	      case 2: /* number of equality constraints > n */
	      default:
                  nlopt_stop_msg(stop, "bug: invalid slsqp mode %d", mode);
		  ret = NLOPT_INVALID_ARGS;
		  goto done;
	  }
	  prev_mode = mode;

	  /* update best point so far */
	  if (nlopt_isfinite(fcur) && ((fcur < *minf && (feasible_cur || !feasible))
				       || (!feasible && infeasibility_cur < infeasibility))) {
	       *minf = fcur;
	       feasible = feasible_cur;
	       infeasibility = infeasibility_cur;
	       memcpy(x, xcur, sizeof(double)*n);
	  }

	  /* note: mode == -1 corresponds to the completion of a line search,
	     and is the only time we should check convergence (as in original slsqp code) */
	  if (mode == -1) {
	       if (!nlopt_isinf(fprev) && feasible_cur) {
		    if (nlopt_stop_ftol(stop, fcur, fprev))
			 ret = NLOPT_FTOL_REACHED;
		    else if (nlopt_stop_x(stop, xcur, xprev))
			 ret = NLOPT_XTOL_REACHED;
	       }
	       fprev = fcur;
	       memcpy(xprev, xcur, sizeof(double)*n);
	  }

	  /* do some additional termination tests */
	  if (nlopt_stop_evals(stop)) ret = NLOPT_MAXEVAL_REACHED;
	  else if (nlopt_stop_time(stop)) ret = NLOPT_MAXTIME_REACHED;
	  else if (feasible && *minf < stop->minf_max) ret = NLOPT_MINF_MAX_REACHED;
     } while (ret == NLOPT_SUCCESS);

done:
     if (nlopt_isinf(*minf)) { /* didn't find any feasible points, just return last point evaluated */
	  if (nlopt_isinf(fcur)) { /* invalid cur. point, use previous pt. */
	       *minf = fprev;
	       memcpy(x, xprev, sizeof(double)*n);
	  }
	  else {
	       *minf = fcur;
	       memcpy(x, xcur, sizeof(double)*n);
	  }
     }

     free(indices);
     free(work);
     return ret;
}
