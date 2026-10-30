/* BLAS/LAPACK routines for NLopt's algorithms.

   They come either from an external BLAS/LAPACK found at configure time
   (CMake option NLOPT_BLAS), as the Fortran symbols name##NLOPT_BLAS_SUFFIX
   with NLOPT_BLAS_INT integers, or, without one, from the built-in versions
   in nlopt-blas-builtin.h (plain C ports of the reference routines, named
   nlopt_##name, with pointer-sized integers, as static inline functions so
   that they can be inlined).  Either way they are called as
   NLOPT_BLAS(name)(...), with the Fortran calling convention (every argument
   by pointer). */

#ifndef NLOPT_BLAS_H
#define NLOPT_BLAS_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef NLOPT_HAVE_BLAS
typedef NLOPT_BLAS_INT blas_int;
#define NLOPT_BLAS_PASTE2(a, b) a##b
#define NLOPT_BLAS_PASTE(a, b) NLOPT_BLAS_PASTE2(a, b)
#define NLOPT_BLAS(name) NLOPT_BLAS_PASTE(name, NLOPT_BLAS_SUFFIX)

void   NLOPT_BLAS(daxpy)(const blas_int* n, const double* alpha, const double* x, const blas_int* incx, double* y, const blas_int* incy);
void   NLOPT_BLAS(dscal)(const blas_int* n, const double* alpha, double* x, const blas_int* incx);
double NLOPT_BLAS(dnrm2)(const blas_int* n, const double* x, const blas_int* incx);
double NLOPT_BLAS(ddot)(const blas_int* n, const double* x, const blas_int* incx, const double* y, const blas_int* incy);
void   NLOPT_BLAS(dgemv)(const char* trans, const blas_int* m, const blas_int* n, const double* alpha, const double* a, const blas_int* lda, const double* x, const blas_int* incx, const double* beta, double* y, const blas_int* incy);
void   NLOPT_BLAS(dtpmv)(const char* uplo, const char* trans, const char* diag, const blas_int* n, const double* ap, double* x, const blas_int* incx);
void   NLOPT_BLAS(dtpsv)(const char* uplo, const char* trans, const char* diag, const blas_int* n, const double* ap, double* x, const blas_int* incx);
void   NLOPT_BLAS(dtrsm)(const char* side, const char* uplo, const char* transa, const char* diag, const blas_int* m, const blas_int* n, const double* alpha, const double* a, const blas_int* lda, double* b, const blas_int* ldb);
void   NLOPT_BLAS(dtrsv)(const char* uplo, const char* trans, const char* diag, const blas_int* n, const double* a, const blas_int* lda, double* x, const blas_int* incx);
void   NLOPT_BLAS(dgelsy)(const blas_int* m, const blas_int* n, const blas_int* nrhs, double* a, const blas_int* lda, double* b, const blas_int* ldb, blas_int* jpvt, const double* rcond, blas_int* rank, double* work, const blas_int* lwork, blas_int* info);
void   NLOPT_BLAS(dgeqr2)(const blas_int* m, const blas_int* n, double* a, const blas_int* lda, double* tau, double* work, blas_int* info);
void   NLOPT_BLAS(dgerq2)(const blas_int* m, const blas_int* n, double* a, const blas_int* lda, double* tau, double* work, blas_int* info);
void   NLOPT_BLAS(dlaic1)(const blas_int* job, const blas_int* j, const double* x, const double* sest, const double* w, const double* gamma, double* sestpr, double* s, double* c);
void   NLOPT_BLAS(dlarf)(const char* side, const blas_int* m, const blas_int* n, const double* v, const blas_int* incv, const double* tau, double* c, const blas_int* ldc, double* work);
void   NLOPT_BLAS(dlarfgp)(const blas_int* n, double* alpha, double* x, const blas_int* incx, double* tau);
void   NLOPT_BLAS(dlartgp)(const double* f, const double* g, double* cs, double* sn, double* r);
void   NLOPT_BLAS(dorm2r)(const char* side, const char* trans, const blas_int* m, const blas_int* n, const blas_int* k, const double* a, const blas_int* lda, const double* tau, double* c, const blas_int* ldc, double* work, blas_int* info);
void   NLOPT_BLAS(dormr2)(const char* side, const char* trans, const blas_int* m, const blas_int* n, const blas_int* k, const double* a, const blas_int* lda, const double* tau, double* c, const blas_int* ldc, double* work, blas_int* info);
#else
typedef ptrdiff_t blas_int;
#define NLOPT_BLAS(name) nlopt_##name
#include "nlopt-blas-builtin.h"
#endif

#ifdef __cplusplus
}
#endif

#endif /* NLOPT_BLAS_H */
