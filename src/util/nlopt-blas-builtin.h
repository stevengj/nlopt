/* Built-in implementations of the BLAS/LAPACK routines in nlopt-blas.h, for
   builds without an external BLAS/LAPACK (CMake option NLOPT_BLAS).  They
   are static inline functions, included by nlopt-blas.h, so that each
   algorithm can inline them.

   These are plain C ports of the reference BLAS and LAPACK 3.12 routines,
   restricted to the arguments that NLopt's algorithms pass (noted per
   routine).  Vector increments must be positive.

   Reference BLAS/LAPACK, distributed under the following license:

   Copyright (c) 1992-2023 The University of Tennessee and The University
                           of Tennessee Research Foundation.  All rights
                           reserved.
   Copyright (c) 2000-2023 The University of California Berkeley. All
                           rights reserved.
   Copyright (c) 2006-2023 The University of Colorado Denver.  All rights
                           reserved.

   $COPYRIGHT$

   Additional copyrights may follow

   $HEADER$

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

   - Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer listed
     in this license in the documentation and/or other materials
     provided with the distribution.

   - Neither the name of the copyright holders nor the names of its
     contributors may be used to endorse or promote products derived from
     this software without specific prior written permission.

   The copyright holders provide no reassurances that the source code
   provided does not infringe any patent, copyright, or any other
   intellectual property rights of third parties.  The copyright holders
   disclaim any liability to any recipient for claims brought against
   recipient by any third party for infringement of that parties
   intellectual property rights.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef NLOPT_BLAS_BUILTIN_H
#define NLOPT_BLAS_BUILTIN_H

#include <math.h>

/* Machine constants as returned by DLAMCH. */
#define NLOPT_DLAMCH_E 1.1102230246251565e-16   /* eps (rounding) = 2^-53 */
#define NLOPT_DLAMCH_P 2.220446049250313e-16    /* eps*base = 2^-52 */
#define NLOPT_DLAMCH_S 2.2250738585072014e-308  /* safe minimum = 2^-1022 */

static inline int nlopt_lsame(const char* ca, char cb) { return (ca[0] == cb) || (ca[0] == cb + ('a' - 'A')); }

/******************************** BLAS 1 *********************************/

static inline void
nlopt_daxpy(const blas_int* n, const double* alpha, const double* x, const blas_int* incx, double* y, const blas_int* incy)
{
    if ((*n <= 0) || (*alpha == 0.0)) { return; }
    for (blas_int i = 0; i < *n; i++) { y[i*(*incy)] += (*alpha)*x[i*(*incx)]; }
}

static inline void
nlopt_dscal(const blas_int* n, const double* alpha, double* x, const blas_int* incx)
{
    for (blas_int i = 0; i < *n; i++) { x[i*(*incx)] *= *alpha; }
}

static inline double
nlopt_ddot(const blas_int* n, const double* x, const blas_int* incx, const double* y, const blas_int* incy)
{
    double dtemp = 0.0;
    for (blas_int i = 0; i < *n; i++) { dtemp += x[i*(*incx)]*y[i*(*incy)]; }
    return dtemp;
}

/* Blue's algorithm, as in the reference DNRM2 (BLAS 3.10). */
static inline double
nlopt_dnrm2(const blas_int* n, const double* x, const blas_int* incx)
{
    const double maxn = 1.7976931348623157e308;
    const double tsml = ldexp(1.0, -511), tbig = ldexp(1.0, 486);
    const double ssml = ldexp(1.0, 537), sbig = ldexp(1.0, -538);
    double abig = 0.0, amed = 0.0, asml = 0.0, scl = 1.0, sumsq = 0.0;
    int notbig = 1;

    if (*n <= 0) { return 0.0; }
    for (blas_int i = 0; i < *n; i++)
    {
        double ax = fabs(x[i*(*incx)]);
        if (ax > tbig) {
            abig += (ax*sbig)*(ax*sbig);
            notbig = 0;
        } else if (ax < tsml) {
            if (notbig) { asml += (ax*ssml)*(ax*ssml); }
        } else {
            amed += ax*ax;
        }
    }
    if (abig > 0.0) {
        if ((amed > 0.0) || (amed > maxn) || (amed != amed)) { abig += (amed*sbig)*sbig; }
        scl = 1.0/sbig;
        sumsq = abig;
    } else if (asml > 0.0) {
        if ((amed > 0.0) || (amed > maxn) || (amed != amed)) {
            double ymin, ymax;
            amed = sqrt(amed);
            asml = sqrt(asml)/ssml;
            if (asml > amed) { ymin = amed; ymax = asml; } else { ymin = asml; ymax = amed; }
            scl = 1.0;
            sumsq = ymax*ymax*(1.0 + (ymin/ymax)*(ymin/ymax));
        } else {
            scl = 1.0/ssml;
            sumsq = asml;
        }
    } else {
        scl = 1.0;
        sumsq = amed;
    }
    return scl*sqrt(sumsq);
}

/******************************** BLAS 2/3 *******************************/

/* trans = "N" or "T" */
static inline void
nlopt_dgemv(const char* trans, const blas_int* m, const blas_int* n, const double* alpha,
            const double* a, const blas_int* lda, const double* x, const blas_int* incx,
            const double* beta, double* y, const blas_int* incy)
{
    blas_int leny = nlopt_lsame(trans, 'N') ? *m : *n;

    if ((*m == 0) || (*n == 0) || ((*alpha == 0.0) && (*beta == 1.0))) { return; }
    if (*beta != 1.0)
    {
        if (*beta == 0.0) {
            for (blas_int i = 0; i < leny; i++) { y[i*(*incy)] = 0.0; }
        } else {
            for (blas_int i = 0; i < leny; i++) { y[i*(*incy)] *= *beta; }
        }
    }
    if (*alpha == 0.0) { return; }
    if (nlopt_lsame(trans, 'N')) {
        // y += alpha*A*x
        for (blas_int j = 0; j < *n; j++)
        {
            double temp = (*alpha)*x[j*(*incx)];
            for (blas_int i = 0; i < *m; i++) { y[i*(*incy)] += temp*a[i + j*(*lda)]; }
        }
    } else {
        // y += alpha*A^T*x
        for (blas_int j = 0; j < *n; j++)
        {
            double temp = 0.0;
            for (blas_int i = 0; i < *m; i++) { temp += a[i + j*(*lda)]*x[i*(*incx)]; }
            y[j*(*incy)] += (*alpha)*temp;
        }
    }
}

/* Packed lower triangular matrix, column-major: A(i,j), i >= j. */
#define NLOPT_BLAS_TP_L(ap, n, i, j) ((ap)[(j)*(n) - ((j)*((j) - 1))/2 + (i) - (j)])

/* uplo = "L", trans = "N" or "T" */
static inline void
nlopt_dtpmv(const char* uplo, const char* trans, const char* diag, const blas_int* n,
            const double* ap, double* x, const blas_int* incx)
{
    blas_int nn = *n, inc = *incx;
    int nounit = nlopt_lsame(diag, 'N');
    (void)uplo;

    if (nlopt_lsame(trans, 'N')) {
        // x := L*x
        for (blas_int j = nn - 1; j >= 0; j--)
        {
            double temp = x[j*inc];
            if (temp != 0.0)
            {
                for (blas_int i = nn - 1; i > j; i--) { x[i*inc] += temp*NLOPT_BLAS_TP_L(ap, nn, i, j); }
                if (nounit) { x[j*inc] *= NLOPT_BLAS_TP_L(ap, nn, j, j); }
            }
        }
    } else {
        // x := L^T*x
        for (blas_int j = 0; j < nn; j++)
        {
            double temp = x[j*inc];
            if (nounit) { temp *= NLOPT_BLAS_TP_L(ap, nn, j, j); }
            for (blas_int i = j + 1; i < nn; i++) { temp += NLOPT_BLAS_TP_L(ap, nn, i, j)*x[i*inc]; }
            x[j*inc] = temp;
        }
    }
}

/* uplo = "L", trans = "N" */
static inline void
nlopt_dtpsv(const char* uplo, const char* trans, const char* diag, const blas_int* n,
            const double* ap, double* x, const blas_int* incx)
{
    blas_int nn = *n, inc = *incx;
    int nounit = nlopt_lsame(diag, 'N');
    (void)uplo; (void)trans;

    // Solve L*x = b
    for (blas_int j = 0; j < nn; j++)
    {
        if (x[j*inc] != 0.0)
        {
            if (nounit) { x[j*inc] /= NLOPT_BLAS_TP_L(ap, nn, j, j); }
            double temp = x[j*inc];
            for (blas_int i = j + 1; i < nn; i++) { x[i*inc] -= temp*NLOPT_BLAS_TP_L(ap, nn, i, j); }
        }
    }
}

/* uplo = "U", trans = "N" or "T" */
static inline void
nlopt_dtrsv(const char* uplo, const char* trans, const char* diag, const blas_int* n,
            const double* a, const blas_int* lda, double* x, const blas_int* incx)
{
    blas_int nn = *n, ld = *lda, inc = *incx;
    int nounit = nlopt_lsame(diag, 'N');
    (void)uplo;

    if (nlopt_lsame(trans, 'N')) {
        // Solve U*x = b
        for (blas_int j = nn - 1; j >= 0; j--)
        {
            if (x[j*inc] != 0.0)
            {
                if (nounit) { x[j*inc] /= a[j + j*ld]; }
                double temp = x[j*inc];
                for (blas_int i = j - 1; i >= 0; i--) { x[i*inc] -= temp*a[i + j*ld]; }
            }
        }
    } else {
        // Solve U^T*x = b
        for (blas_int j = 0; j < nn; j++)
        {
            double temp = x[j*inc];
            for (blas_int i = 0; i < j; i++) { temp -= a[i + j*ld]*x[i*inc]; }
            if (nounit) { temp /= a[j + j*ld]; }
            x[j*inc] = temp;
        }
    }
}

/* uplo = "U", transa = "N", side = "L" or "R": solve A*X = alpha*B or X*A = alpha*B */
static inline void
nlopt_dtrsm(const char* side, const char* uplo, const char* transa, const char* diag,
            const blas_int* m, const blas_int* n, const double* alpha, const double* a,
            const blas_int* lda, double* b, const blas_int* ldb)
{
    blas_int mm = *m, nn = *n, la = *lda, lb = *ldb;
    int nounit = nlopt_lsame(diag, 'N');
    (void)uplo; (void)transa;

    if ((mm == 0) || (nn == 0)) { return; }
    if (nlopt_lsame(side, 'L'))
    {
        for (blas_int j = 0; j < nn; j++)
        {
            if (*alpha != 1.0)
            {
                for (blas_int i = 0; i < mm; i++) { b[i + j*lb] *= *alpha; }
            }
            for (blas_int k = mm - 1; k >= 0; k--)
            {
                if (b[k + j*lb] != 0.0)
                {
                    if (nounit) { b[k + j*lb] /= a[k + k*la]; }
                    double temp = b[k + j*lb];
                    for (blas_int i = 0; i < k; i++) { b[i + j*lb] -= temp*a[i + k*la]; }
                }
            }
        }
        return;
    }
    for (blas_int j = 0; j < nn; j++)
    {
        if (*alpha != 1.0)
        {
            for (blas_int i = 0; i < mm; i++) { b[i + j*lb] *= *alpha; }
        }
        for (blas_int k = 0; k < j; k++)
        {
            double akj = a[k + j*la];
            if (akj != 0.0)
            {
                for (blas_int i = 0; i < mm; i++) { b[i + j*lb] -= akj*b[i + k*lb]; }
            }
        }
        if (nounit)
        {
            double temp = 1.0/a[j + j*la];
            for (blas_int i = 0; i < mm; i++) { b[i + j*lb] *= temp; }
        }
    }
}

/******************************** LAPACK *********************************/

static inline double
nlopt_dlapy2(double x, double y)
{
    int x_is_nan = (x != x), y_is_nan = (y != y);
    if (x_is_nan) { return x; }
    if (y_is_nan) { return y; }
    double xabs = fabs(x), yabs = fabs(y);
    double w = (xabs > yabs ? xabs : yabs), z = (xabs < yabs ? xabs : yabs);
    if ((z == 0.0) || (w > 1.7976931348623157e308)) { return w; }
    return w*sqrt(1.0 + (z/w)*(z/w));
}

/* Apply H = I - tau*v*v^T to the (m x n) matrix C from the left (side "L") or
   the right (side "R"). v[one_at*incv], if one_at >= 0, is taken to be 1. */
static inline void
nlopt_larf(int left, blas_int m, blas_int n, const double* v, blas_int incv, blas_int one_at,
           double tau, double* c, blas_int ldc, double* work)
{
    blas_int lastv;

    if (tau == 0.0) { return; }
    // Look for the last non-zero entry of v.
    lastv = (left ? m : n);
    while ((lastv > 0) && (lastv - 1 != one_at) && (v[(lastv - 1)*incv] == 0.0)) { lastv--; }
#define NLOPT_BLAS_V(i) ((i) == one_at ? 1.0 : v[(i)*incv])
    if (left) {
        // work = C(0:lastv, :)^T*v; C -= tau*v*work^T
        for (blas_int j = 0; j < n; j++)
        {
            double temp = 0.0;
            for (blas_int i = 0; i < lastv; i++) { temp += c[i + j*ldc]*NLOPT_BLAS_V(i); }
            work[j] = temp;
        }
        for (blas_int j = 0; j < n; j++)
        {
            if (work[j] != 0.0)
            {
                double temp = -tau*work[j];
                for (blas_int i = 0; i < lastv; i++) { c[i + j*ldc] += NLOPT_BLAS_V(i)*temp; }
            }
        }
    } else {
        // work = C(:, 0:lastv)*v; C -= tau*work*v^T
        for (blas_int i = 0; i < m; i++) { work[i] = 0.0; }
        for (blas_int j = 0; j < lastv; j++)
        {
            double temp = NLOPT_BLAS_V(j);
            for (blas_int i = 0; i < m; i++) { work[i] += temp*c[i + j*ldc]; }
        }
        for (blas_int j = 0; j < lastv; j++)
        {
            double vj = NLOPT_BLAS_V(j);
            if (vj != 0.0)
            {
                double temp = -tau*vj;
                for (blas_int i = 0; i < m; i++) { c[i + j*ldc] += work[i]*temp; }
            }
        }
    }
#undef NLOPT_BLAS_V
}

static inline void
nlopt_dlarf(const char* side, const blas_int* m, const blas_int* n, const double* v,
            const blas_int* incv, const double* tau, double* c, const blas_int* ldc, double* work)
{
    nlopt_larf(nlopt_lsame(side, 'L'), *m, *n, v, *incv, -1, *tau, c, *ldc, work);
}

/* Generate an elementary reflector H with H*(alpha; x) = (beta; 0). */
static inline void
nlopt_dlarfg(blas_int n, double* alpha, double* x, blas_int incx, double* tau)
{
    blas_int nm1 = n - 1, knt = 0;
    double xnorm, beta, safmin, scal;

    if (n <= 1) { *tau = 0.0; return; }
    xnorm = nlopt_dnrm2(&nm1, x, &incx);
    if (xnorm == 0.0) { *tau = 0.0; return; }
    beta = -copysign(nlopt_dlapy2(*alpha, xnorm), *alpha);
    safmin = NLOPT_DLAMCH_S/NLOPT_DLAMCH_E;
    if (fabs(beta) < safmin)
    {
        double rsafmn = 1.0/safmin;
        do {
            knt++;
            nlopt_dscal(&nm1, &rsafmn, x, &incx);
            beta *= rsafmn;
            *alpha *= rsafmn;
        } while ((fabs(beta) < safmin) && (knt < 20));
        xnorm = nlopt_dnrm2(&nm1, x, &incx);
        beta = -copysign(nlopt_dlapy2(*alpha, xnorm), *alpha);
    }
    *tau = (beta - *alpha)/beta;
    scal = 1.0/(*alpha - beta);
    nlopt_dscal(&nm1, &scal, x, &incx);
    for (blas_int j = 0; j < knt; j++) { beta *= safmin; }
    *alpha = beta;
}

/* As dlarfg, but with beta >= 0. */
static inline void
nlopt_dlarfgp(const blas_int* n, double* alpha, double* x, const blas_int* incx, double* tau)
{
    blas_int nm1 = *n - 1, knt = 0;
    double xnorm, beta, smlnum, savealpha, scal;

    if (*n <= 0) { *tau = 0.0; return; }
    xnorm = nlopt_dnrm2(&nm1, x, incx);
    if (xnorm <= NLOPT_DLAMCH_P*fabs(*alpha)) {
        if (*alpha >= 0.0) {
            *tau = 0.0;
        } else {
            *tau = 2.0;
            for (blas_int j = 0; j < nm1; j++) { x[j*(*incx)] = 0.0; }
            *alpha = -*alpha;
        }
        return;
    }
    beta = copysign(nlopt_dlapy2(*alpha, xnorm), *alpha);
    smlnum = NLOPT_DLAMCH_S/NLOPT_DLAMCH_E;
    if (fabs(beta) < smlnum)
    {
        double bignum = 1.0/smlnum;
        do {
            knt++;
            nlopt_dscal(&nm1, &bignum, x, incx);
            beta *= bignum;
            *alpha *= bignum;
        } while ((fabs(beta) < smlnum) && (knt < 20));
        xnorm = nlopt_dnrm2(&nm1, x, incx);
        beta = copysign(nlopt_dlapy2(*alpha, xnorm), *alpha);
    }
    savealpha = *alpha;
    *alpha += beta;
    if (beta < 0.0) {
        beta = -beta;
        *tau = -*alpha/beta;
    } else {
        *alpha = xnorm*(xnorm/(*alpha));
        *tau = *alpha/beta;
        *alpha = -*alpha;
    }
    if (fabs(*tau) <= smlnum) {
        if (savealpha >= 0.0) {
            *tau = 0.0;
        } else {
            *tau = 2.0;
            for (blas_int j = 0; j < nm1; j++) { x[j*(*incx)] = 0.0; }
            beta = -savealpha;
        }
    } else {
        scal = 1.0/(*alpha);
        nlopt_dscal(&nm1, &scal, x, incx);
    }
    for (blas_int j = 0; j < knt; j++) { beta *= smlnum; }
    *alpha = beta;
}

/* Plane rotation with r >= 0. */
static inline void
nlopt_dlartgp(const double* f, const double* g, double* cs, double* sn, double* r)
{
    const double safmn2 = ldexp(1.0, -484), safmx2 = ldexp(1.0, 484);
    double f1 = *f, g1 = *g, scale;
    int count = 0;

    if (*g == 0.0) {
        *cs = copysign(1.0, *f);
        *sn = 0.0;
        *r = fabs(*f);
        return;
    }
    if (*f == 0.0) {
        *cs = 0.0;
        *sn = copysign(1.0, *g);
        *r = fabs(*g);
        return;
    }
    scale = fmax(fabs(f1), fabs(g1));
    if (scale >= safmx2) {
        do {
            count++;
            f1 *= safmn2;
            g1 *= safmn2;
            scale = fmax(fabs(f1), fabs(g1));
        } while ((scale >= safmx2) && (count < 20));
        *r = sqrt(f1*f1 + g1*g1);
        *cs = f1/(*r);
        *sn = g1/(*r);
        for (int i = 0; i < count; i++) { *r *= safmx2; }
    } else if (scale <= safmn2) {
        do {
            count++;
            f1 *= safmx2;
            g1 *= safmx2;
            scale = fmax(fabs(f1), fabs(g1));
        } while (scale <= safmn2);
        *r = sqrt(f1*f1 + g1*g1);
        *cs = f1/(*r);
        *sn = g1/(*r);
        for (int i = 0; i < count; i++) { *r *= safmn2; }
    } else {
        *r = sqrt(f1*f1 + g1*g1);
        *cs = f1/(*r);
        *sn = g1/(*r);
    }
    if (*r < 0.0)
    {
        *cs = -*cs;
        *sn = -*sn;
        *r = -*r;
    }
}

/* QR factorization, A = Q*R. */
static inline void
nlopt_dgeqr2(const blas_int* m, const blas_int* n, double* a, const blas_int* lda,
             double* tau, double* work, blas_int* info)
{
    blas_int mm = *m, nn = *n, ld = *lda, k = (mm < nn ? mm : nn);

    *info = 0;
    for (blas_int i = 0; i < k; i++)
    {
        nlopt_dlarfg(mm - i, &a[i + i*ld], &a[(i + 1 < mm ? i + 1 : mm - 1) + i*ld], 1, &tau[i]);
        if (i < nn - 1)
        {
            nlopt_larf(1, mm - i, nn - i - 1, &a[i + i*ld], 1, 0, tau[i], &a[i + (i + 1)*ld], ld, work);
        }
    }
}

/* RQ factorization, A = R*Q. */
static inline void
nlopt_dgerq2(const blas_int* m, const blas_int* n, double* a, const blas_int* lda,
             double* tau, double* work, blas_int* info)
{
    blas_int mm = *m, nn = *n, ld = *lda, k = (mm < nn ? mm : nn);

    *info = 0;
    for (blas_int i = k - 1; i >= 0; i--)
    {
        // Row mm-k+i, reflector over columns 0 .. nn-k+i
        blas_int r = mm - k + i, c = nn - k + i;
        nlopt_dlarfg(c + 1, &a[r + c*ld], &a[r], ld, &tau[i]);
        nlopt_larf(0, r, c + 1, &a[r], ld, c, tau[i], a, ld, work);
    }
}

/* side = "L", trans = "T": C := Q^T*C with Q from dgeqr2. */
static inline void
nlopt_dorm2r(const char* side, const char* trans, const blas_int* m, const blas_int* n,
             const blas_int* k, const double* a, const blas_int* lda, const double* tau,
             double* c, const blas_int* ldc, double* work, blas_int* info)
{
    blas_int ld = *lda;
    (void)side; (void)trans;

    *info = 0;
    if ((*m == 0) || (*n == 0) || (*k == 0)) { return; }
    for (blas_int i = 0; i < *k; i++)
    {
        nlopt_larf(1, *m - i, *n, &a[i + i*ld], 1, 0, tau[i], &c[i], *ldc, work);
    }
}

/* trans = "T", side = "L" or "R": C := Q^T*C or C*Q^T with Q from dgerq2. */
static inline void
nlopt_dormr2(const char* side, const char* trans, const blas_int* m, const blas_int* n,
             const blas_int* k, const double* a, const blas_int* lda, const double* tau,
             double* c, const blas_int* ldc, double* work, blas_int* info)
{
    blas_int kk = *k, ld = *lda;
    (void)trans;

    *info = 0;
    if ((*m == 0) || (*n == 0) || (kk == 0)) { return; }
    if (nlopt_lsame(side, 'L')) {
        for (blas_int i = 0; i < kk; i++)
        {
            blas_int mi = *m - kk + i + 1;
            nlopt_larf(1, mi, *n, &a[i], ld, mi - 1, tau[i], c, *ldc, work);
        }
    } else {
        for (blas_int i = kk - 1; i >= 0; i--)
        {
            blas_int ni = *n - kk + i + 1;
            nlopt_larf(0, *m, ni, &a[i], ld, ni - 1, tau[i], c, *ldc, work);
        }
    }
}

/* Incremental condition estimation (DLAIC1); job 1 = largest, 2 = smallest
   singular value. */
static inline void
nlopt_dlaic1(const blas_int* job, const blas_int* j, const double* x, const double* sestp,
             const double* w, const double* gammap, double* sestpr, double* s, double* c)
{
    const double eps = NLOPT_DLAMCH_E;
    const double sest = *sestp, gamma = *gammap;
    blas_int one = 1;
    double alpha = nlopt_ddot(j, x, &one, w, &one);
    double absalp = fabs(alpha), absgam = fabs(gamma), absest = fabs(sest);
    double s1, s2, tmp, b, t, zeta1, zeta2, sine, cosine, norma, test;

    if (*job == 1) {
        if (sest == 0.0) {
            s1 = fmax(absgam, absalp);
            if (s1 == 0.0) {
                *s = 0.0; *c = 1.0; *sestpr = 0.0;
            } else {
                *s = alpha/s1;
                *c = gamma/s1;
                tmp = sqrt((*s)*(*s) + (*c)*(*c));
                *s /= tmp;
                *c /= tmp;
                *sestpr = s1*tmp;
            }
        } else if (absgam <= eps*absest) {
            *s = 1.0; *c = 0.0;
            tmp = fmax(absest, absalp);
            s1 = absest/tmp;
            s2 = absalp/tmp;
            *sestpr = tmp*sqrt(s1*s1 + s2*s2);
        } else if (absalp <= eps*absest) {
            s1 = absgam; s2 = absest;
            if (s1 <= s2) { *s = 1.0; *c = 0.0; *sestpr = s2; }
            else { *s = 0.0; *c = 1.0; *sestpr = s1; }
        } else if ((absest <= eps*absalp) || (absest <= eps*absgam)) {
            s1 = absgam; s2 = absalp;
            if (s1 <= s2) {
                tmp = s1/s2;
                *s = sqrt(1.0 + tmp*tmp);
                *sestpr = s2*(*s);
                *c = (gamma/s2)/(*s);
                *s = copysign(1.0, alpha)/(*s);
            } else {
                tmp = s2/s1;
                *c = sqrt(1.0 + tmp*tmp);
                *sestpr = s1*(*c);
                *s = (alpha/s1)/(*c);
                *c = copysign(1.0, gamma)/(*c);
            }
        } else {
            zeta1 = alpha/absest;
            zeta2 = gamma/absest;
            b = (1.0 - zeta1*zeta1 - zeta2*zeta2)*0.5;
            *c = zeta1*zeta1;
            if (b > 0.0) { t = (*c)/(b + sqrt(b*b + (*c))); }
            else { t = sqrt(b*b + (*c)) - b; }
            sine = -zeta1/t;
            cosine = -zeta2/(1.0 + t);
            tmp = sqrt(sine*sine + cosine*cosine);
            *s = sine/tmp;
            *c = cosine/tmp;
            *sestpr = sqrt(t + 1.0)*absest;
        }
    } else {
        if (sest == 0.0) {
            *sestpr = 0.0;
            if (fmax(absgam, absalp) == 0.0) { sine = 1.0; cosine = 0.0; }
            else { sine = -gamma; cosine = alpha; }
            s1 = fmax(fabs(sine), fabs(cosine));
            *s = sine/s1;
            *c = cosine/s1;
            tmp = sqrt((*s)*(*s) + (*c)*(*c));
            *s /= tmp;
            *c /= tmp;
        } else if (absgam <= eps*absest) {
            *s = 0.0; *c = 1.0; *sestpr = absgam;
        } else if (absalp <= eps*absest) {
            s1 = absgam; s2 = absest;
            if (s1 <= s2) { *s = 0.0; *c = 1.0; *sestpr = s1; }
            else { *s = 1.0; *c = 0.0; *sestpr = s2; }
        } else if ((absest <= eps*absalp) || (absest <= eps*absgam)) {
            s1 = absgam; s2 = absalp;
            if (s1 <= s2) {
                tmp = s1/s2;
                *c = sqrt(1.0 + tmp*tmp);
                *sestpr = absest*(tmp/(*c));
                *s = -(gamma/s2)/(*c);
                *c = copysign(1.0, alpha)/(*c);
            } else {
                tmp = s2/s1;
                *s = sqrt(1.0 + tmp*tmp);
                *sestpr = absest/(*s);
                *c = (alpha/s1)/(*s);
                *s = -copysign(1.0, gamma)/(*s);
            }
        } else {
            zeta1 = alpha/absest;
            zeta2 = gamma/absest;
            norma = fmax(1.0 + zeta1*zeta1 + fabs(zeta1*zeta2), fabs(zeta1*zeta2) + zeta2*zeta2);
            test = 1.0 + 2.0*(zeta1 - zeta2)*(zeta1 + zeta2);
            if (test >= 0.0) {
                b = (zeta1*zeta1 + zeta2*zeta2 + 1.0)*0.5;
                *c = zeta2*zeta2;
                t = (*c)/(b + sqrt(fabs(b*b - (*c))));
                sine = zeta1/(1.0 - t);
                cosine = -zeta2/t;
                *sestpr = sqrt(t + 4.0*eps*eps*norma)*absest;
            } else {
                b = (zeta2*zeta2 + zeta1*zeta1 - 1.0)*0.5;
                *c = zeta1*zeta1;
                if (b >= 0.0) { t = -(*c)/(b + sqrt(b*b + (*c))); }
                else { t = b - sqrt(b*b + (*c)); }
                sine = -zeta1/t;
                cosine = -zeta2/(1.0 + t);
                *sestpr = sqrt(1.0 + t + 4.0*eps*eps*norma)*absest;
            }
            tmp = sqrt(sine*sine + cosine*cosine);
            *s = sine/tmp;
            *c = cosine/tmp;
        }
    }
}

/* Least squares solution of A*X = B via QR with the rank determined as in
   DGELSY.  slsqp.c marks all columns as fixed (jpvt != 0), so there is no
   column pivoting.  Unlike DGELSY, A and B are not rescaled (Householder QR
   is scale invariant short of overflow), and for rank < n the basic solution
   (zeros in the trailing components) is returned instead of the minimum norm
   one; slsqp.c only uses the rank in that case (mode 7). */
static inline void
nlopt_dgelsy(const blas_int* m, const blas_int* n, const blas_int* nrhs, double* a,
             const blas_int* lda, double* b, const blas_int* ldb, blas_int* jpvt,
             const double* rcond, blas_int* rank, double* work, const blas_int* lwork,
             blas_int* info)
{
    blas_int mm = *m, nn = *n, ld = *lda, mn = (mm < nn ? mm : nn);
    double* tau = work;
    double* xmin = &work[mn];
    double* xmax = &work[2*mn];
    double smax, smin;
    (void)lwork;

    *info = 0;
    *rank = 0;
    if ((mn == 0) || (*nrhs == 0)) { return; }
    for (blas_int i = 0; i < nn; i++) { jpvt[i] = i + 1; }

    nlopt_dgeqr2(m, n, a, lda, tau, &work[mn], info);

    // Determine the rank using incremental condition estimation.
    xmin[0] = 1.0;
    xmax[0] = 1.0;
    smax = fabs(a[0]);
    smin = smax;
    if (smax == 0.0)
    {
        for (blas_int j = 0; j < *nrhs; j++)
        {
            for (blas_int i = 0; i < (mm > nn ? mm : nn); i++) { b[i + j*(*ldb)] = 0.0; }
        }
        return;
    }
    *rank = 1;
    while (*rank < mn)
    {
        blas_int i = *rank, imax = 1, imin = 2;
        double sminpr, smaxpr, s1, c1, s2, c2;
        nlopt_dlaic1(&imin, &i, xmin, &smin, &a[i*ld], &a[i + i*ld], &sminpr, &s1, &c1);
        nlopt_dlaic1(&imax, &i, xmax, &smax, &a[i*ld], &a[i + i*ld], &smaxpr, &s2, &c2);
        if (!(smaxpr*(*rcond) <= sminpr)) { break; }
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

    // B := Q^T*B, solve R(0:rank, 0:rank)*X = B(0:rank) and zero the rest.
    nlopt_dorm2r("L", "T", m, nrhs, &mn, a, lda, tau, b, ldb, &work[mn], info);
    for (blas_int j = 0; j < *nrhs; j++)
    {
        blas_int one = 1;
        nlopt_dtrsv("U", "N", "N", rank, a, lda, &b[j*(*ldb)], &one);
        for (blas_int i = *rank; i < nn; i++) { b[i + j*(*ldb)] = 0.0; }
    }
}

#undef NLOPT_BLAS_TP_L

#endif /* NLOPT_BLAS_BUILTIN_H */
