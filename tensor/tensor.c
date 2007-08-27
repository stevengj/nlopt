/* tensor.f -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__5 = 5;
static integer c__1 = 1;
static integer c__9 = 9;
static integer c__3 = 3;
static doublereal c_b111 = 2.;
static doublereal c_b134 = 10.;
static doublereal c_b250 = .33333333333333331;
static doublereal c_b324 = 1.;
static doublereal c_b384 = 3.;

/*      ALGORITHM 739, COLLECTED ALGORITHMS FROM ACM. */
/*      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE, */
/*      VOL. 20, NO. 4, DECEMBER, 1994, P.518-530. */
/* *** driver.f */
/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_211[] = "(\002    G=\002,999e20.13)";

    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;
    alist al__1, al__2;

    /* Builtin functions */
    integer f_open(olist *), s_rsle(cilist *), do_lio(integer *, integer *, 
	    char *, ftnlen), e_rsle(void), s_wsle(cilist *), e_wsle(void), 
	    s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     f_rew(alist *), f_end(alist *), f_clos(cllist *);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    doublereal h__[200]	/* was [50][4] */;
    integer i__, n;
    doublereal x[50];
    integer nr;
    extern /* Subroutine */ int fcn_(), grd_();
    integer msg;
    extern /* Subroutine */ int hsn_();
    integer ipr;
    doublereal wrk[400]	/* was [50][8] */, fpls, gpls[50];
    integer iwrk[4];
    doublereal xpls[50], typx[50];
    integer itnno;
    doublereal fscale;
    integer grdflg, hesflg;
    doublereal gradtl;
    integer ndigit;
    extern /* Subroutine */ int dfault_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    integer method;
    extern /* Subroutine */ int prtfcn_(integer *);
    integer itnlim;
    extern /* Subroutine */ int tensrd_(integer *, integer *, doublereal *, 
	    U_fp, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), tensor_(
	    integer *, integer *, doublereal *, U_fp, U_fp, U_fp, doublereal *
	    , doublereal *, doublereal *, doublereal *, integer *, doublereal 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    doublereal steptl, stepmx;

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 5, 0, 0, 0 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };
    static cilist io___7 = { 0, 6, 0, 0, 0 };
    static cilist io___8 = { 0, 6, 0, 0, 0 };
    static cilist io___9 = { 0, 6, 0, 0, 0 };
    static cilist io___10 = { 0, 6, 0, 0, 0 };
    static cilist io___11 = { 0, 6, 0, 0, 0 };
    static cilist io___12 = { 0, 6, 0, 0, 0 };
    static cilist io___21 = { 0, 6, 0, 0, 0 };
    static cilist io___22 = { 0, 6, 0, fmt_211, 0 };
    static cilist io___23 = { 0, 6, 0, 0, 0 };
    static cilist io___24 = { 0, 6, 0, 0, 0 };
    static cilist io___25 = { 0, 5, 0, 0, 0 };
    static cilist io___26 = { 0, 6, 0, 0, 0 };
    static cilist io___27 = { 0, 6, 0, 0, 0 };
    static cilist io___28 = { 0, 6, 0, 0, 0 };
    static cilist io___29 = { 0, 6, 0, 0, 0 };
    static cilist io___30 = { 0, 6, 0, 0, 0 };
    static cilist io___31 = { 0, 6, 0, 0, 0 };
    static cilist io___32 = { 0, 6, 0, 0, 0 };
    static cilist io___33 = { 0, 6, 0, 0, 0 };
    static cilist io___45 = { 0, 6, 0, 0, 0 };
    static cilist io___46 = { 0, 6, 0, fmt_211, 0 };
    static cilist io___47 = { 0, 6, 0, 0, 0 };
    static cilist io___48 = { 0, 6, 0, 0, 0 };


    nr = 50;
    o__1.oerr = 0;
    o__1.ounit = 5;
    o__1.ofnmlen = 8;
    o__1.ofnm = "wood.inp";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    prtfcn_(&n);
    s_rsle(&io___3);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_lio(&c__5, &c__1, (char *)&x[i__ - 1], (ftnlen)sizeof(doublereal));
    }
    e_rsle();
    s_wsle(&io___6);
    do_lio(&c__9, &c__1, "INITIAL APPROXIMATION TO THE SOLUTION X*:", (ftnlen)
	    41);
    e_wsle();
    s_wsle(&io___7);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_lio(&c__5, &c__1, (char *)&x[i__ - 1], (ftnlen)sizeof(doublereal));
    }
    e_wsle();
/*  SHORT FORM */
    s_wsle(&io___8);
    do_lio(&c__9, &c__1, " ", (ftnlen)1);
    e_wsle();
    s_wsle(&io___9);
    do_lio(&c__9, &c__1, "CALLING TENSRD - ALL INPUT PARAMETERS  ARE SET", (
	    ftnlen)46);
    e_wsle();
    s_wsle(&io___10);
    do_lio(&c__9, &c__1, "                 TO THEIR DEFAULT VALUES.", (ftnlen)
	    41);
    e_wsle();
    s_wsle(&io___11);
    do_lio(&c__9, &c__1, " OUTPUT WILL BE WRITTEN TO THE STANDARD OUTPUT.", (
	    ftnlen)47);
    e_wsle();
    s_wsle(&io___12);
    do_lio(&c__9, &c__1, " ", (ftnlen)1);
    e_wsle();
    tensrd_(&nr, &n, x, (U_fp)fcn_, &msg, xpls, &fpls, gpls, h__, &itnno, wrk,
	     iwrk);
    s_wsle(&io___21);
    do_lio(&c__9, &c__1, "RESULTS OF TENSRD:", (ftnlen)18);
    e_wsle();
    s_wsfe(&io___22);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&gpls[i__ - 1], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    s_wsle(&io___23);
    do_lio(&c__9, &c__1, "XPLS=", (ftnlen)5);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_lio(&c__5, &c__1, (char *)&xpls[i__ - 1], (ftnlen)sizeof(
		doublereal));
    }
    e_wsle();
    s_wsle(&io___24);
    do_lio(&c__9, &c__1, "FPLS=", (ftnlen)5);
    do_lio(&c__5, &c__1, (char *)&fpls, (ftnlen)sizeof(doublereal));
    do_lio(&c__9, &c__1, "  ITNNO=", (ftnlen)8);
    do_lio(&c__3, &c__1, (char *)&itnno, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, "  MSG=", (ftnlen)6);
    do_lio(&c__3, &c__1, (char *)&msg, (ftnlen)sizeof(integer));
    e_wsle();
/*  LONG FORM */
    prtfcn_(&n);
    al__1.aerr = 0;
    al__1.aunit = 5;
    f_rew(&al__1);
    s_rsle(&io___25);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_lio(&c__5, &c__1, (char *)&x[i__ - 1], (ftnlen)sizeof(doublereal));
    }
    e_rsle();
    s_wsle(&io___26);
    do_lio(&c__9, &c__1, "INITIAL APPROXIMATION TO THE SOLUTION X*:", (ftnlen)
	    41);
    e_wsle();
    s_wsle(&io___27);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_lio(&c__5, &c__1, (char *)&x[i__ - 1], (ftnlen)sizeof(doublereal));
    }
    e_wsle();
    s_wsle(&io___28);
    do_lio(&c__9, &c__1, " ", (ftnlen)1);
    e_wsle();
    s_wsle(&io___29);
    do_lio(&c__9, &c__1, "CALLING TENSOR AFTER RESETTING INPUT PARAMETERS", (
	    ftnlen)47);
    e_wsle();
    s_wsle(&io___30);
    do_lio(&c__9, &c__1, "IPR, MSG AND METHOD.", (ftnlen)20);
    e_wsle();
    s_wsle(&io___31);
    do_lio(&c__9, &c__1, "THE STANDARD METHOD IS USED IN THIS EXAMPLE.", (
	    ftnlen)44);
    e_wsle();
    s_wsle(&io___32);
    do_lio(&c__9, &c__1, "ADDITIONAL OUTPUT WILL BE WRITTEN TO FILE 'FORT8'.",
	     (ftnlen)50);
    e_wsle();
    s_wsle(&io___33);
    do_lio(&c__9, &c__1, " ", (ftnlen)1);
    e_wsle();
/* L997: */
    dfault_(&n, typx, &fscale, &gradtl, &steptl, &itnlim, &stepmx, &ipr, &
	    method, &grdflg, &hesflg, &ndigit, &msg);
    o__1.oerr = 0;
    o__1.ounit = 8;
    o__1.ofnmlen = 5;
    o__1.ofnm = "FORT8";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    ipr = 8;
    msg = 2;
    method = 0;
    tensor_(&nr, &n, x, (U_fp)fcn_, (U_fp)grd_, (U_fp)hsn_, typx, &fscale, &
	    gradtl, &steptl, &itnlim, &stepmx, &ipr, &method, &grdflg, &
	    hesflg, &ndigit, &msg, xpls, &fpls, gpls, h__, &itnno, wrk, iwrk);
    s_wsle(&io___45);
    do_lio(&c__9, &c__1, "RESULTS OF TENSOR:", (ftnlen)18);
    e_wsle();
    s_wsfe(&io___46);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&gpls[i__ - 1], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    s_wsle(&io___47);
    do_lio(&c__9, &c__1, "XPLS=", (ftnlen)5);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_lio(&c__5, &c__1, (char *)&xpls[i__ - 1], (ftnlen)sizeof(
		doublereal));
    }
    e_wsle();
    s_wsle(&io___48);
    do_lio(&c__9, &c__1, "FPLS=", (ftnlen)5);
    do_lio(&c__5, &c__1, (char *)&fpls, (ftnlen)sizeof(doublereal));
    do_lio(&c__9, &c__1, "  ITNNO=", (ftnlen)8);
    do_lio(&c__3, &c__1, (char *)&itnno, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, "  MSG=", (ftnlen)6);
    do_lio(&c__3, &c__1, (char *)&msg, (ftnlen)sizeof(integer));
    e_wsle();
    al__2.aerr = 0;
    al__2.aunit = 8;
    f_end(&al__2);
    cl__1.cerr = 0;
    cl__1.cunit = 8;
    cl__1.csta = 0;
    f_clos(&cl__1);
    cl__1.cerr = 0;
    cl__1.cunit = 5;
    cl__1.csta = 0;
    f_clos(&cl__1);
    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

/* Subroutine */ int fcn_(integer *n, doublereal *x, doublereal *f)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*n <= 0) {
	*n = 4;
    } else {
	d__1 = x[1] * x[1] - x[2];
	d__2 = 1. - x[1];
	d__3 = x[3] * x[3] - x[4];
	d__4 = 1. - x[3];
	d__5 = 1. - x[2];
	d__6 = 1. - x[4];
	*f = pow_dd(&d__1, &c_b111) * 100. + pow_dd(&d__2, &c_b111) + pow_dd(&
		d__3, &c_b111) * 90. + pow_dd(&d__4, &c_b111) + (pow_dd(&d__5,
		 &c_b111) + pow_dd(&d__6, &c_b111)) * 10.1 + (1. - x[2]) * 
		19.8 * (1. - x[4]);
    }
    return 0;
} /* fcn_ */

/* Subroutine */ int prtfcn_(integer *n)
{
    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Fortran I/O blocks */
    static cilist io___49 = { 0, 6, 0, 0, 0 };
    static cilist io___50 = { 0, 6, 0, 0, 0 };
    static cilist io___51 = { 0, 6, 0, 0, 0 };


    *n = 4;
    s_wsle(&io___49);
    do_lio(&c__9, &c__1, "__________________________________________", (
	    ftnlen)42);
    e_wsle();
    s_wsle(&io___50);
    do_lio(&c__9, &c__1, "f(x)= Wood function", (ftnlen)19);
    e_wsle();
    s_wsle(&io___51);
    do_lio(&c__9, &c__1, "__________________________________________", (
	    ftnlen)42);
    e_wsle();
    return 0;
} /* prtfcn_ */

/*  ---------------------- */
/*  |  G R D              | */
/*  ---------------------- */
/* Subroutine */ int grd_(integer *n, real *x, real *g)
{
/*     DUMMY ROUTINE */
    return 0;
} /* grd_ */

/*  ---------------------- */
/*  |  H S N              | */
/*  ---------------------- */
/* Subroutine */ int hsn_(integer *nr, integer *n, real *x, real *h__)
{
/*     DUMMY ROUTINE */
    return 0;
} /* hsn_ */

/* *** tensrd.f */
/*  ---------------------- */
/*  |  T E N S O R       | */
/*  ---------------------- */
/* Subroutine */ int tensor_(integer *nr, integer *n, doublereal *x, U_fp fcn,
	 U_fp grd, U_fp hsn, doublereal *typx, doublereal *fscale, doublereal 
	*gradtl, doublereal *steptl, integer *itnlim, doublereal *stepmx, 
	integer *ipr, integer *method, integer *grdflg, integer *hesflg, 
	integer *ndigit, integer *msg, doublereal *xpls, doublereal *fpls, 
	doublereal *gpls, doublereal *h__, integer *itnno, doublereal *wrk, 
	integer *iwrk)
{
    /* System generated locals */
    integer h_dim1, h_offset, wrk_dim1, wrk_offset;

    /* Local variables */
    extern /* Subroutine */ int opt_(integer *, integer *, doublereal *, U_fp,
	     U_fp, U_fp, doublereal *, doublereal *, doublereal *, doublereal 
	    *, integer *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *);


/* AUTHORS:   TA-TUNG CHOW, ELIZABETH ESKOW AND ROBERT B. SCHNABEL */

/* DATE:      MARCH 29, 1991 */

/* PURPOSE: */
/*   DERIVATIVE TENSOR METHOD FOR UNCONSTRAINED OPTIMIZATION */
/*     THE METHOD BASES EACH ITERATION ON A SPECIALLY CONSTRUCTED */
/*     FOURTH ORDER MODEL OF THE OBJECTIVE FUNCTION.  THE MODEL */
/*     INTERPOLATES THE FUNCTION VALUE AND GRADIENT FROM THE PREVIOUS */
/*     ITERATE AND THE CURRENT FUNCTION VALUE, GRADIENT AND HESSIAN */
/*     MATRIX. */

/* BLAS SUBROUTINES: DCOPY,DDOT,DSCAL */
/* UNCMIN SUBROUTINES: DFAULT, OPTCHK, GRDCHK, HESCHK, LNSRCH, FSTOFD, */
/*                     SNDOFD, BAKSLV, FORSLV, OPTSTP */
/* MODCHL SUBROUTINES: MODCHL, INIT, GERCH, FIN2X2 */

/* ----------------------------------------------------------------------- */

/* PARAMETERS: */

/*   NR      --> ROW DIMENSION OF MATRIX */
/*   N       --> DIMENSION OF PROBLEM */
/*   X(N)    --> INITIAL GUESS (INPUT) AND CURRENT POINT */
/*   FCN     --> NAME OF SUBROUTINE TO EVALUATE FUNCTION VALUE */
/*   GRD     --> NAME OF SUBROUTINE TO EVALUATE ANALYTICAL GRADIENT */
/*   HSN     --> NAME OF SUBROUTINE TO EVALUATE ANALYTICAL HESSIAN */
/*               HSN SHOULD FILL ONLY THE LOWER TRIANGULAR PART */
/*               AND DIAGONAL OF H */
/*   TYPX(N) --> TYPICAL SIZE OF EACH COMPONENT OF X */
/*   FSCALE  --> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION FCN */
/*   GRADTL  --> GRADIENT TOLERANCE */
/*   STEPTL  --> STEP TOLERANCE */
/*   ITNLIM  --> ITERATION LIMIT */
/*   STEPMX  --> MAXIMUM STEP LENGTH ALLOWED */
/*   IPR     --> OUTPUT UNIT NUMBER */
/*   METHOD  --> IF VALUE IS 0 THEN USE ONLY NEWTON STEP AT */
/*               EACH ITERATION */
/*               IF VALUE IS 1 THEN TRY BOTH TENSOR AND NEWTON */
/*               STEPS AT EACH ITERATION */
/*   GRDFLG  --> = 1 OR 2 IF ANALYTICAL GRADIENT IS SUPPLIED */
/*   HESFLG  --> = 1 OR 2 IF ANALYTICAL HESSIAN IS SUPPLIED */
/*   NDIGIT  --> NUMBER OF GOOD DIGITS IN OPTIMIZATION FUNCTION FCN */
/*   MSG     --> OUTPUT MESSAGE CONTROL */
/*   XPLS(N) <--  NEW POINT AND FINAL POINT (OUTPUT) */
/*   FPLS    <--  NEW FUNCTION VALUE AND FINAL FUNCTION VALUE (OUTPUT) */
/*   GPLS(N) <--  CURRENT GRADIENT AND GRADIENT AT FINAL POINT (OUTPUT) */
/*   H(N,N)  --> HESSIAN */
/*   ITNNO   <--  NUMBER OF ITERATIONS */
/*   WRK(N,8)--> WORK SPACE */
/*   IWRK(N) --> WORK SPACE */


/*   EQUIVALENCE WRK(N,1) = G(N) */
/*               WRK(N,2) = S(N) */
/*               WRK(N,3) = D(N) */
/*               WRK(N,4) = DN(N) */
/*               WRK(N,6) = E(N) */
/*               WRK(N,7) = WK1(N) */
/*               WRK(N,8) = WK2(N) */

    /* Parameter adjustments */
    wrk_dim1 = *nr;
    wrk_offset = 1 + wrk_dim1;
    wrk -= wrk_offset;
    --iwrk;
    h_dim1 = *nr;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --gpls;
    --xpls;
    --typx;
    --x;

    /* Function Body */
    opt_(nr, n, &x[1], (U_fp)fcn, (U_fp)grd, (U_fp)hsn, &typx[1], fscale, 
	    gradtl, steptl, itnlim, stepmx, ipr, method, grdflg, hesflg, 
	    ndigit, msg, &xpls[1], fpls, &gpls[1], &h__[h_offset], itnno, &
	    wrk[wrk_dim1 + 1], &wrk[(wrk_dim1 << 1) + 1], &wrk[wrk_dim1 * 3 + 
	    1], &wrk[(wrk_dim1 << 2) + 1], &wrk[wrk_dim1 * 6 + 1], &wrk[
	    wrk_dim1 * 7 + 1], &wrk[(wrk_dim1 << 3) + 1], &iwrk[1]);
    return 0;
} /* tensor_ */

/*  ---------------------- */
/*  |  T E N S R D       | */
/*  ---------------------- */
/* Subroutine */ int tensrd_(integer *nr, integer *n, doublereal *x, U_fp fcn,
	 integer *msg, doublereal *xpls, doublereal *fpls, doublereal *gpls, 
	doublereal *h__, integer *itnno, doublereal *wrk, integer *iwrk)
{
    /* System generated locals */
    integer h_dim1, h_offset, wrk_dim1, wrk_offset;

    /* Local variables */
    extern /* Subroutine */ int grd_(integer *, real *, real *), hsn_(integer 
	    *, integer *, real *, real *);
    integer ipr;
    doublereal fscale;
    integer grdflg, hesflg;
    doublereal gradtl;
    integer ndigit;
    extern /* Subroutine */ int dfault_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    integer method, itnlim;
    extern /* Subroutine */ int tensor_(integer *, integer *, doublereal *, 
	    U_fp, S_fp, S_fp, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *);
    doublereal steptl, stepmx;


/* PURPOSE: */
/*   SHORT CALLING SEQUENCE SUBROUTINE */

/* ----------------------------------------------------------------------- */

/* PARAMETERS: */

/*   NR         --> ROW DIMENSION OF MATRIX */
/*   N          --> DIMENSION OF PROBLEM */
/*   X(N)       --> INITIAL GUESS (INPUT) AND CURRENT POINT */
/*   FCN        --> NAME OF SUBROUTINE TO EVALUATE FUNCTION VALUE */
/*   MSG        --> OUTPUT MESSAGE CONTROL */
/*   XPLS(N)    <--  NEW POINT AND FINAL POINT (OUTPUT) */
/*   FPLS       <--  NEW FUNCTION VALUE AND FINAL FUNCTION VALUE (OUTPUT) */
/*   GPLS(N)    <--  GRADIENT AT FINAL POINT (OUTPUT) */
/*   H(N,N)     --> HESSIAN */
/*   ITNNO      <--  NUMBER OF ITERATIONS */
/*   WRK(N,8)   --> WORK SPACE */
/*   IWRK(N)    --> WORK SPACE */


/*   EQUIVALENCE WRK(N,1) = G(N) */
/*               WRK(N,2) = S(N) */
/*               WRK(N,3) = D(N) */
/*               WRK(N,4) = DN(N) */
/*               WRK(N,5) = TYPX(N) */
/*               WRK(N,6) = E(N) */
/*               WRK(N,7) = WK1(N) */
/*               WRK(N,8) = WK2(N) */

    /* Parameter adjustments */
    wrk_dim1 = *nr;
    wrk_offset = 1 + wrk_dim1;
    wrk -= wrk_offset;
    --iwrk;
    h_dim1 = *nr;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --gpls;
    --xpls;
    --x;

    /* Function Body */
    dfault_(n, &wrk[wrk_dim1 * 5 + 1], &fscale, &gradtl, &steptl, &itnlim, &
	    stepmx, &ipr, &method, &grdflg, &hesflg, &ndigit, msg);
    tensor_(nr, n, &x[1], (U_fp)fcn, (S_fp)grd_, (S_fp)hsn_, &wrk[wrk_dim1 * 
	    5 + 1], &fscale, &gradtl, &steptl, &itnlim, &stepmx, &ipr, &
	    method, &grdflg, &hesflg, &ndigit, msg, &xpls[1], fpls, &gpls[1], 
	    &h__[h_offset], itnno, &wrk[wrk_offset], &iwrk[1]);
    return 0;
} /* tensrd_ */

/*  ---------------------- */
/*  |       O P T        | */
/*  ---------------------- */
/* Subroutine */ int opt_(integer *nr, integer *n, doublereal *x, S_fp fcn, 
	S_fp grd, S_fp hsn, doublereal *typx, doublereal *fscale, doublereal *
	gradtl, doublereal *steptl, integer *itnlim, doublereal *stepmx, 
	integer *ipr, integer *method, integer *grdflg, integer *hesflg, 
	integer *ndigit, integer *msg, doublereal *xpls, doublereal *fpls, 
	doublereal *gpls, doublereal *h__, integer *itnno, doublereal *g, 
	doublereal *s, doublereal *d__, doublereal *dn, doublereal *e, 
	doublereal *wk1, doublereal *wk2, integer *pivot)
{
    /* Format strings */
    static char fmt_25[] = "(\002 INITIAL FUNCTION VALUE F=\002,e20.13)";
    static char fmt_30[] = "(\002 INITIAL GRADIENT G=\002,999e20.13)";
    static char fmt_103[] = "(\002 -----------    ITERATION \002,i4,\002 ---"
	    "-------------\002)";
    static char fmt_104[] = "(\002 X=\002,999e20.13)";
    static char fmt_105[] = "(\002 F(X)=\002,e20.13)";
    static char fmt_106[] = "(\002 G(X)=\002,999e20.13)";
    static char fmt_108[] = "(\002FUNCTION EVAL COUNT:\002,i6,\002 REL. GRAD"
	    ". MAX:\002,e20.13)";
    static char fmt_110[] = "(\002REL. STEP MAX :\002,e20.13)";
    static char fmt_370[] = "(\002 FINAL X=\002,999e20.13)";
    static char fmt_380[] = "(\002 GRADIENT G=\002,999e20.13)";
    static char fmt_390[] = "(\002 FUNCTION VALUE F(X)=\002,e20.13,\002 AT I"
	    "TERATION \002,i4)";
    static char fmt_400[] = "(\002FUNCTION EVAL COUNT:\002,i6,\002 REL. GRAD"
	    ". MAX:\002,e20.13)";
    static char fmt_410[] = "(\002REL. STEP MAX :\002,e20.13)";
    static char fmt_560[] = "(\002 -----------    ITERATION \002,i4,\002 ---"
	    "-------------\002)";
    static char fmt_570[] = "(\002 X=\002,999e20.13)";
    static char fmt_580[] = "(\002 F(X)=\002,e20.13)";
    static char fmt_590[] = "(\002 G(X)=\002,999e20.13)";
    static char fmt_600[] = "(\002FUNCTION EVAL COUNT:\002,i6,\002 REL. GRAD"
	    ". MAX:\002,e20.13)";
    static char fmt_610[] = "(\002REL. STEP MAX :\002,e20.13)";

    /* System generated locals */
    integer h_dim1, h_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_di(doublereal *, integer *), sqrt(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    doublereal f;
    integer i__;
    doublereal gd, fn, fp, rnf, eps, rgx, rsx, beta;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    integer imsg;
    doublereal temp, alpha;
    extern /* Subroutine */ int mkmdl_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    integer nfcnt;
    doublereal finit;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    logical nomin;
    doublereal gnorm, addmax;
    extern /* Subroutine */ int grdchk_(integer *, doublereal *, S_fp, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *), heschk_(integer *, integer *, doublereal *
	    , S_fp, S_fp, S_fp, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *);
    integer iretcd;
    doublereal analtl;
    extern /* Subroutine */ int choldr_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *), mcheps_(doublereal *), sndofd_(integer *, integer 
	    *, doublereal *, S_fp, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    integer itrmcd;
    extern /* Subroutine */ int bakslv_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), fstofd_(integer *, integer *, 
	    integer *, doublereal *, S_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *);
    integer icscmx;
    extern /* Subroutine */ int optchk_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *);
    logical mxtake;
    extern /* Subroutine */ int lnsrch_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, logical *, integer *, doublereal *, doublereal *, 
	    doublereal *, S_fp, doublereal *, integer *), slvmdl_(integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, logical 
	    *, doublereal *), forslv_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    extern doublereal twonrm_(integer *, doublereal *);
    extern /* Subroutine */ int optstp_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    logical *, integer *, integer *, doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___70 = { 0, 0, 0, fmt_25, 0 };
    static cilist io___71 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___81 = { 0, 0, 0, fmt_103, 0 };
    static cilist io___82 = { 0, 0, 0, fmt_104, 0 };
    static cilist io___83 = { 0, 0, 0, fmt_105, 0 };
    static cilist io___84 = { 0, 0, 0, fmt_106, 0 };
    static cilist io___85 = { 0, 0, 0, fmt_108, 0 };
    static cilist io___86 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___93 = { 0, 0, 0, fmt_370, 0 };
    static cilist io___94 = { 0, 0, 0, fmt_380, 0 };
    static cilist io___95 = { 0, 0, 0, fmt_390, 0 };
    static cilist io___96 = { 0, 0, 0, fmt_400, 0 };
    static cilist io___97 = { 0, 0, 0, fmt_410, 0 };
    static cilist io___98 = { 0, 0, 0, fmt_560, 0 };
    static cilist io___99 = { 0, 0, 0, fmt_570, 0 };
    static cilist io___100 = { 0, 0, 0, fmt_580, 0 };
    static cilist io___101 = { 0, 0, 0, fmt_590, 0 };
    static cilist io___102 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___103 = { 0, 0, 0, fmt_610, 0 };



/* PURPOSE: */
/*   DERIVATIVE TENSOR METHODS FOR UNCONSTRAINED OPTIMIZATION */

/* ----------------------------------------------------------------------- */

/* PARAMETERS: */

/*   NR      --> ROW DIMENSION OF MATRIX */
/*   N       --> DIMENSION OF PROBLEM */
/*   X(N)    --> INITIAL GUESS (INPUT) AND CURRENT POINT */
/*   FCN     --> NAME OF SUBROUTINE TO EVALUATE FUNCTION VALUE */
/*               FCN: R(N) --> R(1) */
/*   GRD     --> NAME OF SUBROUTINE TO EVALUATE ANALYTICAL GRADIENT */
/*               FCN: R(N) --> R(N) */
/*   HSN     --> NAME OF SUBROUTINE TO EVALUATE ANALYTICAL HESSIAN */
/*               FCN: R(N) --> R(N X N) */
/*               HSN SHOULD FILL ONLY THE LOWER TRIANGULAR PART */
/*               AND DIAGONAL OF H */
/*   TYPX(N) --> TYPICAL SIZE OF EACH COMPONENT OF X */
/*   FSCALE  --> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION FCN */
/*   GRADTL  --> GRADIENT TOLERANCE */
/*   STEPTL  --> STEP TOLERANCE */
/*   ITNLIM  --> ITERATION LIMIT */
/*   STEPMX  --> MAXIMUM STEP LENGTH ALLOWED */
/*   IPR     --> OUTPUT UNIT NUMBER */
/*   METHOD  --> IF VALUE IS 0 THEN USE ONLY NEWTON STEP AT */
/*               EACH ITERATION */
/*   GRDFLG  --> = 1 OR 2 IF ANALYTICAL GRADIENT IS SUPPLIED */
/*   HESFLG  --> = 1 OR 2 IF ANALYTICAL HESSIAN IS SUPPLIED */
/*   NDIGIT  --> NUMBER OF GOOD DIGITS IN OPTIMIZATION FUNCTION FCN */
/*   MSG     --> OUTPUT MESSAGE CONTRO */
/*   XPLS(N) <--  NEW POINT AND FINAL POINT (OUTPUT) */
/*   FPLS    <--  NEW FUNCTION VALUE AND FINAL FUNCTION VALUE (OUTPUT) */
/*   GPLS(N) <--  CURRENT GRADIENT AND GRADIENT AT FINAL POINT (OUTPUT) */
/*   H(N,N)  --> HESSIAN */
/*   ITNNO   <--  NUMBER OF ITERATIONS */
/*   G(N)    --> PREVIOUS GRADIENT */
/*   S       --> STEP TO PREVIOUS ITERATE (FOR TENSOR MODEL) */
/*   D       --> CURRENT STEP (TENSOR OR NEWTON) */
/*   DN      --> NEWTON STEP */
/*   E(N)    --> DIAGONAL ADDED TO HESSIAN IN CHOLESKY DECOMPOSITION */
/*   WK1(N)  --> WORKSPACE */
/*   WK2(N)  --> WORKSPACE */
/*   PIVOT(N)--> PIVOT VECTOR FOR CHOLESKY DECOMPOSITION */

/* -------------------------------- */
/*     INITIALIZATION           | */
/* -------------------------------- */

    /* Parameter adjustments */
    --pivot;
    --wk2;
    --wk1;
    --e;
    --dn;
    --d__;
    --s;
    --g;
    h_dim1 = *nr;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --gpls;
    --xpls;
    --typx;
    --x;

    /* Function Body */
    *itnno = 0;
    icscmx = 0;
    nfcnt = 0;
    mcheps_(&eps);
    optchk_(nr, n, msg, &x[1], &typx[1], fscale, gradtl, steptl, itnlim, 
	    ndigit, &eps, method, grdflg, hesflg, stepmx, ipr);
    if (*msg < 0) {
	return 0;
    }

/*     SCALE X */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] /= typx[i__];
/* L10: */
    }

/* -------------------------------- */
/*     INITIAL ITERATION | */
/* -------------------------------- */

/*     COMPUTE TYPICAL SIZE OF X */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dn[i__] = 1. / typx[i__];
/* L15: */
    }
/* Computing MAX */
    i__1 = -(*ndigit);
    d__1 = pow_di(&c_b134, &i__1);
    rnf = max(d__1,eps);
/* Computing MAX */
    d__1 = .01, d__2 = sqrt(rnf);
    analtl = max(d__1,d__2);
/*     UNSCALE X AND COMPUTE F AND G */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wk1[i__] = x[i__] * typx[i__];
/* L20: */
    }
    (*fcn)(n, &wk1[1], &f);
    ++nfcnt;
    if (*grdflg >= 1) {
	(*grd)(n, &wk1[1], &g[1]);
	if (*grdflg == 1) {
	    *fscale = 1.;
	    grdchk_(n, &wk1[1], (S_fp)fcn, &f, &g[1], &dn[1], &typx[1], 
		    fscale, &rnf, &analtl, &wk2[1], msg, ipr, &nfcnt);
	}
    } else {
	fstofd_(&c__1, &c__1, n, &wk1[1], (S_fp)fcn, &f, &g[1], &typx[1], &
		rnf, &wk2[1], &c__1, &nfcnt);
    }
    if (*msg < -20) {
	return 0;
    }

    gnorm = twonrm_(n, &g[1]);

/*     PRINT OUT 1ST ITERATION? */
    if (*msg >= 1) {
	io___70.ciunit = *ipr;
	s_wsfe(&io___70);
	do_fio(&c__1, (char *)&f, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___71.ciunit = *ipr;
	s_wsfe(&io___71);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&g[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }

/*     TEST WHETHER THE INITIAL GUESS SATISFIES THE STOPPING CRITERIA */
    if (gnorm <= *gradtl) {
	*fpls = f;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xpls[i__] = x[i__];
	    gpls[i__] = g[i__];
/* L40: */
	}
	optstp_(n, &xpls[1], fpls, &gpls[1], &x[1], itnno, &icscmx, &itrmcd, 
		gradtl, steptl, fscale, itnlim, &iretcd, &mxtake, ipr, msg, &
		rgx, &rsx);
	goto L350;
    }
    finit = f;

/* ------------------------ */
/*     ITERATION 1 | */
/* ------------------------ */

    ++(*itnno);

/*     COMPUTE HESSIAN */
    if (*hesflg == 0) {
	if (*grdflg == 1) {
	    fstofd_(nr, n, n, &wk1[1], (S_fp)grd, &g[1], &h__[h_offset], &
		    typx[1], &rnf, &wk2[1], &c__3, &nfcnt);
	} else {
	    sndofd_(nr, n, &wk1[1], (S_fp)fcn, &f, &h__[h_offset], &typx[1], &
		    rnf, &wk2[1], &d__[1], &nfcnt);
	}
    } else {
	if (*hesflg == 2) {
	    (*hsn)(nr, n, &wk1[1], &h__[h_offset]);
	} else {
	    if (*hesflg == 1) {
/*           IN HESCHK GPLS,XPLS AND E ARE USED AS WORK SPACE */
		heschk_(nr, n, &wk1[1], (S_fp)fcn, (S_fp)grd, (S_fp)hsn, &f, &
			g[1], &h__[h_offset], &dn[1], &typx[1], &rnf, &analtl,
			 grdflg, &gpls[1], &xpls[1], &e[1], msg, ipr, &nfcnt);
	    }
	}
    }
    if (*msg < -20) {
	return 0;
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = i__ - 1;
	dcopy_(&i__2, &h__[i__ + h_dim1], nr, &h__[i__ * h_dim1 + 1], &c__1);
/* L50: */
    }

/*     CHOLESKY DECOMPOSITION FOR H (H=LLT) */
    choldr_(nr, n, &h__[h_offset], &d__[1], &eps, &pivot[1], &e[1], &wk1[1], &
	    addmax);

/*     SOLVE FOR NEWTON STEP D */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L60: */
	wk2[i__] = -g[i__];
    }
    forslv_(nr, n, &h__[h_offset], &wk1[1], &wk2[1]);
    bakslv_(nr, n, &h__[h_offset], &d__[1], &wk1[1]);

/*     APPLY LINESEARCH TO THE NEWTON STEP */
    lnsrch_(nr, n, &x[1], &f, &g[1], &d__[1], &xpls[1], fpls, &mxtake, &
	    iretcd, stepmx, steptl, &typx[1], (S_fp)fcn, &wk1[1], &nfcnt);

/*     UPDATE G */
/*      CALL DCOPY(N,GPLS(1),1,GP(1),1) */

/*     UNSCALE XPLS AND COMPUTE GPLS */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wk1[i__] = xpls[i__] * typx[i__];
/* L80: */
    }
    if (*grdflg >= 1) {
	(*grd)(n, &wk1[1], &gpls[1]);
    } else {
	fstofd_(&c__1, &c__1, n, &wk1[1], (S_fp)fcn, fpls, &gpls[1], &typx[1],
		 &rnf, &wk2[1], &c__1, &nfcnt);
    }

/*     CHECK STOPPING CONDITIONS */
    optstp_(n, &xpls[1], fpls, &gpls[1], &x[1], itnno, &icscmx, &itrmcd, 
	    gradtl, steptl, fscale, itnlim, &iretcd, &mxtake, ipr, msg, &rgx, 
	    &rsx);

/*     IF ITRMCD > 0 THEN STOPPING CONDITIONS SATISFIED */
    if (itrmcd > 0) {
	goto L350;
    }

/*     UPDATE X,F AND S FOR TENSOR MODEL */
    fp = f;
    f = *fpls;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp = xpls[i__];
	s[i__] = x[i__] - temp;
	x[i__] = temp;
/* L90: */
    }

/*     IF MSG >= 2 THEN PRINT OUT EACH ITERATION */
    if (*msg >= 2) {
	io___81.ciunit = *ipr;
	s_wsfe(&io___81);
	do_fio(&c__1, (char *)&(*itnno), (ftnlen)sizeof(integer));
	e_wsfe();
	io___82.ciunit = *ipr;
	s_wsfe(&io___82);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
	io___83.ciunit = *ipr;
	s_wsfe(&io___83);
	do_fio(&c__1, (char *)&(*fpls), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___84.ciunit = *ipr;
	s_wsfe(&io___84);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&gpls[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    if (*msg >= 3) {
	io___85.ciunit = *ipr;
	s_wsfe(&io___85);
	do_fio(&c__1, (char *)&nfcnt, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&rgx, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___86.ciunit = *ipr;
	s_wsfe(&io___86);
	do_fio(&c__1, (char *)&rsx, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/* ------------------------ */
/*     ITERATION > 1     | */
/* ------------------------ */

/*     UNSCALE X AND COMPUTE H */
L200:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wk1[i__] = x[i__] * typx[i__];
/* L210: */
    }
    if (*hesflg == 0) {
	if (*grdflg == 1) {
	    fstofd_(nr, n, n, &wk1[1], (S_fp)grd, &g[1], &h__[h_offset], &
		    typx[1], &rnf, &wk2[1], &c__3, &nfcnt);
	} else {
	    sndofd_(nr, n, &wk1[1], (S_fp)fcn, &f, &h__[h_offset], &typx[1], &
		    rnf, &wk2[1], &d__[1], &nfcnt);
	}
    } else {
	(*hsn)(nr, n, &wk1[1], &h__[h_offset]);
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = i__ - 1;
	dcopy_(&i__2, &h__[i__ + h_dim1], nr, &h__[i__ * h_dim1 + 1], &c__1);
/* L230: */
    }

/*     IF METHOD = 0 THEN USE NEWTON STEP ONLY */
    if (*method == 0) {

/*       CHOLESKY DECOMPOSITION FOR H */
	choldr_(nr, n, &h__[h_offset], &wk2[1], &eps, &pivot[1], &e[1], &wk1[
		1], &addmax);

/*       COMPUTE NEWTON STEP */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    wk1[i__] = -gpls[i__];
/* L240: */
	}
	forslv_(nr, n, &h__[h_offset], &wk2[1], &wk1[1]);
	bakslv_(nr, n, &h__[h_offset], &d__[1], &wk2[1]);

/*       NO TENSOR STEP */
	nomin = TRUE_;
	goto L300;

    }

/*     FORM TENSOR MODEL */
    mkmdl_(nr, n, &f, &fp, &gpls[1], &g[1], &s[1], &h__[h_offset], &alpha, &
	    beta, &wk1[1], &d__[1]);

/*   SOLVE TENSOR MODEL AND COMPUTE NEWTON STEP */
/*   ON INPUT : SH IS STORED IN WK1 */
/*            A=(G-GPLS-SH-S*BETA/(6*STS)) IS STORED IN D */
/*   ON OUTPUT: NEWTON STEP IS STORED IN DN */
/*            TENSOR STEP IS STORED IN D */
    slvmdl_(nr, n, &h__[h_offset], &xpls[1], &wk2[1], &e[1], &g[1], &s[1], &
	    gpls[1], &pivot[1], &d__[1], &wk1[1], &dn[1], &alpha, &beta, &
	    nomin, &eps);

/*     IF TENSOR MODEL HAS NO MINIMIZER THEN USE NEWTON STEP */
    if (nomin) {
	dcopy_(n, &dn[1], &c__1, &d__[1], &c__1);
	goto L300;
    }

/*     IF TENSOR STEP IS NOT IN DESCENT DIRECTION THEN USE NEWTON STEP */
    gd = ddot_(n, &gpls[1], &c__1, &d__[1], &c__1);
    if (gd > 0.) {
	dcopy_(n, &dn[1], &c__1, &d__[1], &c__1);
	nomin = TRUE_;
    }

L300:
    ++(*itnno);
    dcopy_(n, &gpls[1], &c__1, &g[1], &c__1);

/*     APPLY LINESEARCH TO TENSOR (OR NEWTON) STEP */
    lnsrch_(nr, n, &x[1], &f, &g[1], &d__[1], &xpls[1], fpls, &mxtake, &
	    iretcd, stepmx, steptl, &typx[1], (S_fp)fcn, &wk1[1], &nfcnt);

    if (! nomin) {
/*       TENSOR STEP IS FOUND AND IN DESCENT DIRECTION, */
/*       APPLY LINESEARCH TO NEWTON STEP */
/*       NEW NEWTON POINT IN WK2 */
	lnsrch_(nr, n, &x[1], &f, &g[1], &dn[1], &wk2[1], &fn, &mxtake, &
		iretcd, stepmx, steptl, &typx[1], (S_fp)fcn, &wk1[1], &nfcnt);

/*       COMPARE TENSOR STEP TO NEWTON STEP */
/*       IF NEWTON STEP IS BETTER, SET NEXT ITERATE TO NEW NEWTON POINT */
	if (fn < *fpls) {
	    *fpls = fn;
	    dcopy_(n, &dn[1], &c__1, &d__[1], &c__1);
	    dcopy_(n, &wk2[1], &c__1, &xpls[1], &c__1);
	}
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__[i__] = xpls[i__] - x[i__];
/* L320: */
    }

/*     UNSCALE XPLS, AND COMPUTE FPLS AND GPLS */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wk1[i__] = xpls[i__] * typx[i__];
/* L330: */
    }
    (*fcn)(n, &wk1[1], fpls);
    ++nfcnt;
    if (*grdflg >= 1) {
	(*grd)(n, &wk1[1], &gpls[1]);
    } else {
	fstofd_(&c__1, &c__1, n, &wk1[1], (S_fp)fcn, fpls, &gpls[1], &typx[1],
		 &rnf, &wk2[1], &c__1, &nfcnt);
    }

/*     CHECK STOPPING CONDITIONS */
    imsg = *msg;
    optstp_(n, &xpls[1], fpls, &gpls[1], &x[1], itnno, &icscmx, &itrmcd, 
	    gradtl, steptl, fscale, itnlim, &iretcd, &mxtake, ipr, msg, &rgx, 
	    &rsx);

/*     IF ITRMCD = 0 THEN NOT OVER YET */
    if (itrmcd == 0) {
	goto L500;
    }

/*     IF MSG >= 1 THEN PRINT OUT FINAL ITERATION */
L350:
    if (imsg >= 1) {

/*       TRANSFORM X BACK TO ORIGINAL SPACE */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xpls[i__] *= typx[i__];
/* L360: */
	}
	io___93.ciunit = *ipr;
	s_wsfe(&io___93);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&xpls[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
	io___94.ciunit = *ipr;
	s_wsfe(&io___94);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&gpls[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
	io___95.ciunit = *ipr;
	s_wsfe(&io___95);
	do_fio(&c__1, (char *)&(*fpls), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*itnno), (ftnlen)sizeof(integer));
	e_wsfe();
	if (imsg >= 3) {
	    io___96.ciunit = *ipr;
	    s_wsfe(&io___96);
	    do_fio(&c__1, (char *)&nfcnt, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&rgx, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    io___97.ciunit = *ipr;
	    s_wsfe(&io___97);
	    do_fio(&c__1, (char *)&rsx, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
/*       UPDATE THE HESSIAN */
	if (*hesflg == 0) {
	    if (*grdflg == 1) {
		fstofd_(nr, n, n, &xpls[1], (S_fp)grd, &gpls[1], &h__[
			h_offset], &typx[1], &rnf, &wk2[1], &c__3, &nfcnt);
	    } else {
		sndofd_(nr, n, &xpls[1], (S_fp)fcn, fpls, &h__[h_offset], &
			typx[1], &rnf, &wk2[1], &d__[1], &nfcnt);
	    }
	} else {
	    (*hsn)(nr, n, &xpls[1], &h__[h_offset]);
	}
    }
    return 0;

/*     UPDATE INFORMATION AT THE CURRENT POINT */
L500:
    dcopy_(n, &xpls[1], &c__1, &x[1], &c__1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s[i__] = -d__[i__];
/* L550: */
    }

/*     IF TOO MANY ITERATIONS THEN RETURN */
    if (*itnno > *itnlim) {
	goto L350;
    }

/*     IF MSG >= 2 THEN PRINT OUT EACH ITERATION */
    if (*msg >= 2) {
	io___98.ciunit = *ipr;
	s_wsfe(&io___98);
	do_fio(&c__1, (char *)&(*itnno), (ftnlen)sizeof(integer));
	e_wsfe();
	io___99.ciunit = *ipr;
	s_wsfe(&io___99);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
	io___100.ciunit = *ipr;
	s_wsfe(&io___100);
	do_fio(&c__1, (char *)&(*fpls), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___101.ciunit = *ipr;
	s_wsfe(&io___101);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&gpls[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    if (*msg >= 3) {
	io___102.ciunit = *ipr;
	s_wsfe(&io___102);
	do_fio(&c__1, (char *)&nfcnt, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&rgx, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___103.ciunit = *ipr;
	s_wsfe(&io___103);
	do_fio(&c__1, (char *)&rsx, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
/*     UPDATE F */
    fp = f;
    f = *fpls;

/*     PERFORM NEXT ITERATION */
    goto L200;

/*     END OF ITERATION > 1 */

} /* opt_ */

/*  ---------------------- */
/*  |  D F A U L T      | */
/*  ---------------------- */
/* Subroutine */ int dfault_(integer *n, doublereal *typx, doublereal *fscale,
	 doublereal *gradtl, doublereal *steptl, integer *itnlim, doublereal *
	stepmx, integer *ipr, integer *method, integer *grdflg, integer *
	hesflg, integer *ndigit, integer *msg)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), d_lg10(doublereal *);

    /* Local variables */
    integer i__;
    doublereal eps, temp;
    extern /* Subroutine */ int mcheps_(doublereal *);

    /* Parameter adjustments */
    --typx;

    /* Function Body */
    mcheps_(&eps);
    *method = 1;
    *fscale = 1.;
    *grdflg = 0;
    *hesflg = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	typx[i__] = 1.;
/* L1: */
    }
    temp = pow_dd(&eps, &c_b250);
    *gradtl = temp;
    *steptl = temp * temp;
    *ndigit = (integer) (-d_lg10(&eps));
/*     SET ACTUAL DFAULT VALUE OF STEPMX IN OPTCHK */
    *stepmx = 0.;
    *itnlim = 100;
    *ipr = 6;
    *msg = 1;
    return 0;
} /* dfault_ */

/*  ---------------------- */
/*  |  O P T C H K       | */
/*  ---------------------- */
/* Subroutine */ int optchk_(integer *nr, integer *n, integer *msg, 
	doublereal *x, doublereal *typx, doublereal *fscale, doublereal *
	gradtl, doublereal *steptl, integer *itnlim, integer *ndigit, 
	doublereal *eps, integer *method, integer *grdflg, integer *hesflg, 
	doublereal *stepmx, integer *ipr)
{
    /* Format strings */
    static char fmt_901[] = "(\0020OPTCHK    ILLEGAL DIMENSION, N=\002,i5)";
    static char fmt_902[] = "(\002OPTCHK    ILLEGAL INPUT VALUES OF: NR=\002"
	    ",i5,\002, N=\002,i5,\002, NR MUST BE <=  N.\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *), d_lg10(
	    doublereal *), pow_di(doublereal *, integer *);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    integer i__;
    doublereal temp, stpsiz;

    /* Fortran I/O blocks */
    static cilist io___110 = { 0, 0, 0, fmt_901, 0 };
    static cilist io___111 = { 0, 0, 0, fmt_902, 0 };



/* PURPOSE */
/* ------- */
/* CHECK INPUT FOR REASONABLENESS */

/* PARAMETERS */
/* ---------- */
/* NR          <--> ROW DIMENSION OF H AND WRK */
/* N            --> DIMENSION OF PROBLEM */
/* X(N)         --> ON ENTRY, ESTIMATE TO ROOT OF FCN */
/* TYPX(N)     <--> TYPICAL SIZE OF EACH COMPONENT OF X */
/* FSCALE      <--> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION FCN */
/* GRADTL      <--> TOLERANCE AT WHICH GRADIENT CONSIDERED CLOSE */
/*                  ENOUGH TO ZERO TO TERMINATE ALGORITHM */
/* STEPTL      <--> TOLERANCE AT WHICH STEP LENGTH CONSIDERED CLOSE */
/*                  ENOUGH TO ZERO TO TERMINATE ALGORITHM */
/* ITNLIM      <--> MAXIMUM NUMBER OF ALLOWABLE ITERATIONS */
/* NDIGIT      <--> NUMBER OF GOOD DIGITS IN OPTIMIZATION FUNCTION FCN */
/* EPS          --> MACHINE PRECISION */
/* METHOD      <--> ALGORITHM INDICATOR */
/* GRDFLG      <--> =1 IF ANALYTIC GRADIENT SUPPLIED */
/* HESFLG      <--> =1 IF ANALYTIC HESSIAN SUPPLIED */
/* STEPMX      <--> MAXIMUM STEP SIZE */
/* MSG         <--> MESSAGE AND ERROR CODE */
/* IPR          --> DEVICE TO WHICH TO SEND OUTPUT */


/* CHECK DIMENSION OF PROBLEM */

    /* Parameter adjustments */
    --typx;
    --x;

    /* Function Body */
    if (*n <= 0) {
	goto L805;
    }
    if (*nr < *n) {
	goto L806;
    }

/* CHECK THAT PARAMETERS ONLY TAKE ON ACCEPTABLE VALUES. */
/* IF NOT, SET THEM TO DEFAULT VALUES. */
    if (*method != 0) {
	*method = 1;
    }
    if (*grdflg != 2 && *grdflg != 1) {
	*grdflg = 0;
    }
    if (*hesflg != 2 && *hesflg != 1) {
	*hesflg = 0;
    }
    if (*msg > 3 || (real) (*msg) < 0.f) {
	*msg = 1;
    }

/* COMPUTE SCALE MATRIX */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (typx[i__] == 0.f) {
	    typx[i__] = 1.;
	}
	if (typx[i__] < 0.f) {
	    typx[i__] = -typx[i__];
	}
/* L10: */
    }

/* CHECK MAXIMUM STEP SIZE */

    if (*stepmx > 0.) {
	goto L20;
    }
    stpsiz = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stpsiz += x[i__] * x[i__] / typx[i__] * typx[i__];
/* L15: */
    }
    stpsiz = sqrt(stpsiz);
/* Computing MAX */
    d__1 = stpsiz * 1e3;
    *stepmx = max(d__1,1e3);
L20:
/* CHECK FUNCTION SCALE */
    if (*fscale == 0.f) {
	*fscale = 1.;
    }
    if (*fscale < 0.f) {
	*fscale = -(*fscale);
    }
/* CHECK GRADIENT TOLERANCE */
    if (*gradtl < 0.f) {
	*gradtl = pow_dd(eps, &c_b250);
    }
/* CHECK STEP TOLERANCE */
    if (*steptl < 0.f) {
	temp = pow_dd(eps, &c_b250);
	*steptl = temp * temp;
    }

/* CHECK ITERATION LIMIT */
    if (*itnlim <= 0) {
	*itnlim = 100;
    }

/* CHECK NUMBER OF DIGITS OF ACCURACY IN FUNCTION FCN */
    if (*ndigit <= 0) {
	*ndigit = (integer) (-d_lg10(eps));
    }
    i__1 = -(*ndigit);
    if (pow_di(&c_b134, &i__1) <= *eps) {
	*ndigit = (integer) (-d_lg10(eps));
    }
    return 0;

/* ERROR EXITS */

L805:
    io___110.ciunit = *ipr;
    s_wsfe(&io___110);
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    e_wsfe();
    *msg = -20;
    return 0;
L806:
    io___111.ciunit = *ipr;
    s_wsfe(&io___111);
    do_fio(&c__1, (char *)&(*nr), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    e_wsfe();
    *msg = -20;
    return 0;
} /* optchk_ */

/*  ---------------------- */
/*  |  G R D C H K       | */
/*  ---------------------- */
/* Subroutine */ int grdchk_(integer *n, doublereal *x, S_fp fcn, doublereal *
	f, doublereal *g, doublereal *typsiz, doublereal *typx, doublereal *
	fscale, doublereal *rnf, doublereal *analtl, doublereal *wrk1, 
	integer *msg, integer *ipr, integer *ifn)
{
    /* Format strings */
    static char fmt_901[] = "(\0020GRDCHK    PROBABLE ERROR IN CODING OF ANA"
	    "LYTIC\002,\002 GRADIENT FUNCTION.\002/\002 GRDCHK     COMP\002,1"
	    "2x,\002ANALYTIC\002,12x,\002ESTIMATE\002)";
    static char fmt_902[] = "(\002 GRDCHK    \002,i5,3x,e20.13,3x,e20.13)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    integer i__;
    doublereal gs;
    integer ker;
    doublereal wrk;
    extern /* Subroutine */ int fstofd_(integer *, integer *, integer *, 
	    doublereal *, S_fp, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___116 = { 0, 0, 0, fmt_901, 0 };
    static cilist io___117 = { 0, 0, 0, fmt_902, 0 };



/* PURPOSE */
/* ------- */
/* CHECK ANALYTIC GRADIENT AGAINST ESTIMATED GRADIENT */

/* PARAMETERS */
/* ---------- */
/* N            --> DIMENSION OF PROBLEM */
/* X(N)         --> ESTIMATE TO A ROOT OF FCN */
/* FCN          --> NAME OF SUBROUTINE TO EVALUATE OPTIMIZATION FUNCTION */
/*                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE */
/*                       FCN:  R(N) --> R(1) */
/* F            --> FUNCTION VALUE:  FCN(X) */
/* G(N)         --> GRADIENT:  G(X) */
/* TYPSIZ(N)    --> TYPICAL SIZE FOR EACH COMPONENT OF X */
/* FSCALE       --> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION FCN */
/* RNF          --> RELATIVE NOISE IN OPTIMIZATION FUNCTION FCN */
/* ANALTL       --> TOLERANCE FOR COMPARISON OF ESTIMATED AND */
/*                  ANALYTICAL GRADIENTS */
/* WRK1(N)      --> WORKSPACE */
/* MSG         <--  MESSAGE OR ERROR CODE */
/*                    ON OUTPUT: =-21, PROBABLE CODING ERROR OF GRADIENT */
/* IPR          --> DEVICE TO WHICH TO SEND OUTPUT */
/* IFN         <--> NUMBER OF FUNCTION EVALUATIONS */

/* COMPUTE FIRST ORDER FINITE DIFFERENCE GRADIENT AND COMPARE TO */
/* ANALYTIC GRADIENT. */

    /* Parameter adjustments */
    --wrk1;
    --typx;
    --typsiz;
    --g;
    --x;

    /* Function Body */
    fstofd_(&c__1, &c__1, n, &x[1], (S_fp)fcn, f, &wrk1[1], &typx[1], rnf, &
	    wrk, &c__1, ifn);
    ker = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__2 = abs(*f);
/* Computing MAX */
	d__3 = (d__1 = x[i__], abs(d__1)), d__4 = typsiz[i__];
	gs = max(d__2,*fscale) / max(d__3,d__4);
/* Computing MAX */
	d__3 = (d__1 = g[i__], abs(d__1));
	if ((d__2 = g[i__] - wrk1[i__], abs(d__2)) > max(d__3,gs) * *analtl) {
	    ker = 1;
	}
/* L5: */
    }
    if (ker == 0) {
	goto L20;
    }
    io___116.ciunit = *ipr;
    s_wsfe(&io___116);
    e_wsfe();
    io___117.ciunit = *ipr;
    s_wsfe(&io___117);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&g[i__], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&wrk1[i__], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    *msg = -21;
L20:
    return 0;
} /* grdchk_ */

/*  ---------------------- */
/*  |  H E S C H K       | */
/*  ---------------------- */
/* Subroutine */ int heschk_(integer *nr, integer *n, doublereal *x, S_fp fcn,
	 S_fp grd, S_fp hsn, doublereal *f, doublereal *g, doublereal *a, 
	doublereal *typsiz, doublereal *typx, doublereal *rnf, doublereal *
	analtl, integer *iagflg, doublereal *udiag, doublereal *wrk1, 
	doublereal *wrk2, integer *msg, integer *ipr, integer *ifn)
{
    /* Format strings */
    static char fmt_901[] = "(\002 HESCHK    PROBABLE ERROR IN CODING OF ANA"
	    "LYTIC\002,\002 HESSIAN FUNCTION.\002/\002 HESCHK      ROW  CO"
	    "L\002,14x,\002ANALYTIC\002,14x,\002(ESTIMATE)\002)";
    static char fmt_902[] = "(\002 HESCHK    \002,2i5,2x,e20.13,2x,\002(\002"
	    ",e20.13,\002)\002)";

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    integer i__, j;
    doublereal hs;
    integer im1, jp1, ker;
    extern /* Subroutine */ int sndofd_(integer *, integer *, doublereal *, 
	    S_fp, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *), fstofd_(integer *, 
	    integer *, integer *, doublereal *, S_fp, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *);

    /* Fortran I/O blocks */
    static cilist io___123 = { 0, 0, 0, fmt_901, 0 };
    static cilist io___125 = { 0, 0, 0, fmt_902, 0 };
    static cilist io___126 = { 0, 0, 0, fmt_902, 0 };



/* PURPOSE */
/* ------- */
/* CHECK ANALYTIC HESSIAN AGAINST ESTIMATED HESSIAN */
/*  (THIS MAY BE DONE ONLY IF THE USER SUPPLIED ANALYTIC HESSIAN */
/*   HSN FILLS ONLY THE LOWER TRIANGULAR PART AND DIAGONAL OF A) */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* X(N)         --> ESTIMATE TO A ROOT OF FCN */
/* FCN          --> NAME OF SUBROUTINE TO EVALUATE OPTIMIZATION FUNCTION */
/*                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE */
/*                       FCN:  R(N) --> R(1) */
/* GRD          --> NAME OF SUBROUTINE TO EVALUATE GRADIENT OF FCN. */
/*                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE */
/* HSN          --> NAME OF SUBROUTINE TO EVALUATE HESSIAN OF FCN. */
/*                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE */
/*                  HSN SHOULD FILL ONLY THE LOWER TRIANGULAR PART */
/*                  AND DIAGONAL OF A */
/* F            --> FUNCTION VALUE:  FCN(X) */
/* G(N)        <--  GRADIENT:  G(X) */
/* A(N,N)      <--  ON EXIT:  HESSIAN IN LOWER TRIANGULAR PART AND DIAG */
/* TYPSIZ(N)    --> TYPICAL SIZE FOR EACH COMPONENT OF X */
/* RNF          --> RELATIVE NOISE IN OPTIMIZATION FUNCTION FCN */
/* ANALTL       --> TOLERANCE FOR COMPARISON OF ESTIMATED AND */
/*                  ANALYTICAL GRADIENTS */
/* IAGFLG       --> =1 IF ANALYTIC GRADIENT SUPPLIED */
/* UDIAG(N)     --> WORKSPACE */
/* WRK1(N)      --> WORKSPACE */
/* WRK2(N)      --> WORKSPACE */
/* MSG         <--> MESSAGE OR ERROR CODE */
/*                    ON INPUT : IF =1XX DO NOT COMPARE ANAL + EST HESS */
/*                    ON OUTPUT: =-22, PROBABLE CODING ERROR OF HESSIAN */
/* IPR          --> DEVICE TO WHICH TO SEND OUTPUT */
/* IFN         <--> NUMBER OF FUNCTION EVALUTATIONS */

/* COMPUTE FINITE DIFFERENCE APPROXIMATION A TO THE HESSIAN. */

    /* Parameter adjustments */
    a_dim1 = *nr;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --wrk2;
    --wrk1;
    --udiag;
    --typx;
    --typsiz;
    --g;
    --x;

    /* Function Body */
    if (*iagflg == 1) {
	fstofd_(nr, n, n, &x[1], (S_fp)grd, &g[1], &a[a_offset], &typx[1], 
		rnf, &wrk1[1], &c__3, ifn);
    }
    if (*iagflg != 1) {
	sndofd_(nr, n, &x[1], (S_fp)fcn, f, &a[a_offset], &typx[1], rnf, &
		wrk1[1], &wrk2[1], ifn);
    }

    ker = 0;

/* COPY LOWER TRIANGULAR PART OF "A" TO UPPER TRIANGULAR PART */
/* AND DIAGONAL OF "A" TO UDIAG */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	udiag[j] = a[j + j * a_dim1];
	if (j == *n) {
	    goto L30;
	}
	jp1 = j + 1;
	i__2 = *n;
	for (i__ = jp1; i__ <= i__2; ++i__) {
	    a[j + i__ * a_dim1] = a[i__ + j * a_dim1];
/* L25: */
	}
L30:
	;
    }

/* COMPUTE ANALYTIC HESSIAN AND COMPARE TO FINITE DIFFERENCE */
/* APPROXIMATION. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
	    a[j + i__ * a_dim1] = 0.;
/* L32: */
	}
    }

    (*hsn)(nr, n, &x[1], &a[a_offset]);
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
/* Computing MAX */
	d__3 = (d__1 = g[j], abs(d__1));
/* Computing MAX */
	d__4 = (d__2 = x[j], abs(d__2)), d__5 = typsiz[j];
	hs = max(d__3,1.) / max(d__4,d__5);
/* Computing MAX */
	d__3 = (d__1 = udiag[j], abs(d__1));
	if ((d__2 = a[j + j * a_dim1] - udiag[j], abs(d__2)) > max(d__3,hs) * 
		*analtl) {
	    ker = 1;
	}
	if (j == *n) {
	    goto L40;
	}
	jp1 = j + 1;
	i__1 = *n;
	for (i__ = jp1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
	    if ((d__2 = a[i__ + j * a_dim1] - a[j + i__ * a_dim1], abs(d__2)) 
		    > max(d__3,hs) * *analtl) {
		ker = 1;
	    }
/* L35: */
	}
L40:
	;
    }

    if (ker == 0) {
	goto L90;
    }
    io___123.ciunit = *ipr;
    s_wsfe(&io___123);
    e_wsfe();
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (i__ == 1) {
	    goto L45;
	}
	im1 = i__ - 1;
	i__1 = im1;
	for (j = 1; j <= i__1; ++j) {
	    io___125.ciunit = *ipr;
	    s_wsfe(&io___125);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&a[i__ + j * a_dim1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&a[j + i__ * a_dim1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
/* L43: */
	}
L45:
	io___126.ciunit = *ipr;
	s_wsfe(&io___126);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&a[i__ + i__ * a_dim1], (ftnlen)sizeof(
		doublereal));
	do_fio(&c__1, (char *)&udiag[i__], (ftnlen)sizeof(doublereal));
	e_wsfe();
/* L50: */
    }
    *msg = -22;
/*     ENDIF */
L90:
    return 0;
} /* heschk_ */

/*  ----------------------- */
/*  |    M C H E P S    | */
/*  ----------------------- */
/* Subroutine */ int mcheps_(doublereal *eps)
{
    doublereal temp, temp1;


/* PURPOSE: */
/*   COMPUTE MACHINE PRECISION */

/* ------------------------------------------------------------------------- */

/* PARAMETERS: */

/*   EPS <-- MACHINE PRECISION */

    temp = 1.;
L20:
    temp /= 2.;
    temp1 = temp + 1.;
    if (temp1 != 1.) {
	goto L20;
    }
    *eps = temp * 2.;
    return 0;
} /* mcheps_ */


doublereal twonrm_(integer *n, doublereal *v)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    doublereal temp;


/* PURPOSE: */
/*   COMPUTER L-2 NORM */

/* -------------------------------------------------------------------------- */

/* PARAMETER: */

/*   N       --> DIMENSION OF PROBLEM */
/*   V(N)    --> VECTOR WHICH L-2 NORM IS EVALUATED */

    /* Parameter adjustments */
    --v;

    /* Function Body */
    temp = ddot_(n, &v[1], &c__1, &v[1], &c__1);
    ret_val = sqrt(temp);
    return ret_val;
} /* twonrm_ */

/*  ---------------------- */
/*  |  L N S R C H       | */
/*  ---------------------- */
/* Subroutine */ int lnsrch_(integer *nr, integer *n, doublereal *x, 
	doublereal *f, doublereal *g, doublereal *p, doublereal *xpls, 
	doublereal *fpls, logical *mxtake, integer *iretcd, doublereal *
	stepmx, doublereal *steptl, doublereal *typx, S_fp fcn, doublereal *
	w2, integer *nfcnt)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    doublereal a, b;
    integer i__, k;
    doublereal t1, t2, t3, scl, rln, sln, slp, tmp, disc;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    doublereal temp, temp1, temp2, alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    doublereal pfpls, almbda, plmbda, tlmbda, almbmn;


/*     THE ALPHA CONDITION ONLY LINE SEARCH */

/* PURPOSE */
/* ------- */
/* FIND A NEXT NEWTON ITERATE BY LINE SEARCH. */

/* PARAMETERS */
/* ---------- */
/* N            --> DIMENSION OF PROBLEM */
/* X(N)         --> OLD ITERATE:   X[K-1] */
/* F            --> FUNCTION VALUE AT OLD ITERATE, F(X) */
/* G(N)         --> GRADIENT AT OLD ITERATE, G(X), OR APPROXIMATE */
/* P(N)         --> NON-ZERO NEWTON STEP */
/* XPLS(N)     <--  NEW ITERATE X[K] */
/* FPLS        <--  FUNCTION VALUE AT NEW ITERATE, F(XPLS) */
/* FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION */
/* IRETCD      <--  RETURN CODE */
/* MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED */
/* STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE */
/* STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES */
/*                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM */
/* TYPX(N)            --> DIAGONAL SCALING MATRIX FOR X (NOT IN UNCMIN) */
/* IPR          --> DEVICE TO WHICH TO SEND OUTPUT */
/* W2         --> WORKING SPACE */
/* NFCNT      <--> NUMBER OF FUNCTION EVALUTIONS */

/* INTERNAL VARIABLES */
/* ------------------ */
/* SLN              NEWTON LENGTH */
/* RLN              RELATIVE LENGTH OF NEWTON STEP */


    /* Parameter adjustments */
    --w2;
    --typx;
    --xpls;
    --p;
    --g;
    --x;

    /* Function Body */
    *mxtake = FALSE_;
    *iretcd = 2;
    alpha = 1e-4;
/* $    WRITE(IPR,954) */
/* $    WRITE(IPR,955) (P(I),I=1,N) */
    tmp = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tmp += p[i__] * p[i__];
/* L5: */
    }
    sln = sqrt(tmp);
    if (sln <= *stepmx) {
	goto L10;
    }

/* NEWTON STEP LONGER THAN MAXIMUM ALLOWED */
    scl = *stepmx / sln;
    dscal_(n, &scl, &p[1], &c__1);
    sln = *stepmx;
/* $     WRITE(IPR,954) */
/* $     WRITE(IPR,955) (P(I),I=1,N) */
L10:
    slp = ddot_(n, &g[1], &c__1, &p[1], &c__1);
    rln = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp = 1.;
	temp1 = (d__1 = x[i__], abs(d__1));
	temp2 = max(temp1,temp);
	temp1 = (d__1 = p[i__], abs(d__1));
/* Computing MAX */
	d__1 = rln, d__2 = temp1 / temp2;
	rln = max(d__1,d__2);
/* L15: */
    }
    almbmn = *steptl / rln;
    almbda = 1.;
/* $    WRITE(IPR,952) SLN,SLP,RMNLMB,STEPMX,STEPTL */

/* LOOP */
/* CHECK IF NEW ITERATE SATISFACTORY.  GENERATE NEW LAMBDA IF NECESSARY. */

L100:
    if (*iretcd < 2) {
	return 0;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xpls[i__] = x[i__] + almbda * p[i__];
/* L105: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	w2[k] = xpls[k] * typx[k];
/* L101: */
    }
    (*fcn)(n, &w2[1], fpls);
    ++(*nfcnt);
/* $    WRITE(IPR,956) ALMBDA */
/* $    WRITE(IPR,951) */
/* $    WRITE(IPR,955) (XPLS(I),I=1,N) */
/* $    WRITE(IPR,953) FPLS */
    if (*fpls > *f + slp * alpha * almbda) {
	goto L130;
    }
/*     IF(FPLS.LE. F+SLP*1.E-4*ALMBDA) */
/*     THEN */

/* SOLUTION FOUND */

    *iretcd = 0;
    if (almbda == 1. && sln > *stepmx * .99) {
	*mxtake = TRUE_;
    }
    goto L100;

/* SOLUTION NOT (YET) FOUND */

/*     ELSE */
L130:
    if (almbda >= almbmn) {
	goto L140;
    }
/*       IF(ALMBDA .LT. ALMBMN) */
/*       THEN */

/* NO SATISFACTORY XPLS FOUND SUFFICIENTLY DISTINCT FROM X */

    *iretcd = 1;
    goto L100;
/*       ELSE */

/* CALCULATE NEW LAMBDA */

L140:
    if (almbda != 1.) {
	goto L150;
    }
/*         IF(ALMBDA.EQ.1.0) */
/*         THEN */

/* FIRST BACKTRACK: QUADRATIC FIT */

    tlmbda = -slp / ((*fpls - *f - slp) * 2.);
    goto L170;
/*         ELSE */

/* ALL SUBSEQUENT BACKTRACKS: CUBIC FIT */

L150:
    t1 = *fpls - *f - almbda * slp;
    t2 = pfpls - *f - plmbda * slp;
    t3 = 1. / (almbda - plmbda);
    a = t3 * (t1 / (almbda * almbda) - t2 / (plmbda * plmbda));
    b = t3 * (t2 * almbda / (plmbda * plmbda) - t1 * plmbda / (almbda * 
	    almbda));
    disc = b * b - a * 3.f * slp;
    if (disc <= b * b) {
	goto L160;
    }
/*           IF(DISC.GT. B*B) */
/*           THEN */

/* ONLY ONE POSITIVE CRITICAL POINT, MUST BE MINIMUM */

    tlmbda = (-b + d_sign(&c_b324, &a) * sqrt(disc)) / (a * 3.f);
    goto L165;
/*           ELSE */

/* BOTH CRITICAL POINTS POSITIVE, FIRST IS MINIMUM */

L160:
    tlmbda = (-b - d_sign(&c_b324, &a) * sqrt(disc)) / (a * 3.f);
/*           ENDIF */
L165:
    if (tlmbda > almbda * .5) {
	tlmbda = almbda * .5;
    }
/*         ENDIF */
L170:
    plmbda = almbda;
    pfpls = *fpls;
    if (tlmbda >= almbda * .1) {
	goto L180;
    }
/*         IF(TLMBDA.LT.ALMBDA/10.) */
/*         THEN */
    almbda *= .1f;
    goto L190;
/*         ELSE */
L180:
    almbda = tlmbda;
/*         ENDIF */
/*       ENDIF */
/*     ENDIF */
L190:
    goto L100;
/* L956: */
/* L951: */
/* L952: */
/* L953: */
/* L954: */
/* L955: */
} /* lnsrch_ */

/*  ---------------------- */
/*  |       Z H Z        | */
/*  ---------------------- */
/* Subroutine */ int zhz_(integer *nr, integer *n, doublereal *y, doublereal *
	h__, doublereal *u, doublereal *t)
{
    /* System generated locals */
    integer h_dim1, h_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal d__;
    integer i__, j;
    doublereal s, sgn;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    doublereal temp1, temp2;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    doublereal ynorm;


/* PURPOSE: */
/*   COMPUTE QTHQ(N,N) AND ZTHZ(N-1,N-1) = FIRST N-1 ROWS AND */
/*   FIRST N-1 COLUMNS OF QTHQ */

/* --------------------------------------------------------------------------- */

/* PARAMETERS: */

/*   NR            --> ROW DIMENSION OF MATRIX */
/*   N       --> DIMENSION OF PROBLEM */
/*   Y(N)    --> FIRST BASIS IN Q */
/*   H(N,N)     <--> ON INPUT : HESSIAN */
/*               ON OUTPUT: QTHQ (ZTHZ) */
/*   U(N)       <--  VECTOR TO FORM Q AND Z */
/*   T(N)        --> WORKSPACE */


/*     U=Y+SGN(Y(N))||Y||E(N) */
    /* Parameter adjustments */
    --t;
    --u;
    h_dim1 = *nr;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --y;

    /* Function Body */
    if (y[*n] != 0.) {
	sgn = y[*n] / (d__1 = y[*n], abs(d__1));
    } else {
	sgn = 1.;
    }
    ynorm = ddot_(n, &y[1], &c__1, &y[1], &c__1);
    ynorm = sqrt(ynorm);
    u[*n] = y[*n] + sgn * ynorm;
    i__1 = *n - 1;
    dcopy_(&i__1, &y[1], &c__1, &u[1], &c__1);

/*     D=UTU/2 */
    d__ = ddot_(n, &u[1], &c__1, &u[1], &c__1);
    d__ /= 2.;

/*     T=2HU/UTU */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t[i__] = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    t[i__] = h__[i__ + j * h_dim1] * u[j] + t[i__];
/* L30: */
	}
	t[i__] /= d__;
/* L40: */
    }

/*     S=4UHU/(UTU)**2 */
    s = ddot_(n, &u[1], &c__1, &t[1], &c__1);
    s /= d__;

/*     COMPUTE QTHQ (ZTHZ) */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp1 = u[i__];
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    temp2 = u[j];
	    h__[i__ + j * h_dim1] = h__[i__ + j * h_dim1] - u[i__] * t[j] - t[
		    i__] * u[j] + u[i__] * u[j] * s;
/* L60: */
	}
/* L70: */
    }
    return 0;
} /* zhz_ */

/*  ---------------------- */
/*  |    S O L V E W     | */
/*  ---------------------- */
/* Subroutine */ int solvew_(integer *nr, integer *n, doublereal *al, 
	doublereal *u, doublereal *w, doublereal *b)
{
    /* System generated locals */
    integer al_dim1, al_offset, i__1, i__2;

    /* Local variables */
    doublereal d__;
    integer i__, j;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int forslv_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);


/* PURPOSE: */
/*   SOLVE L*W=ZT*V */

/* ---------------------------------------------------------------------- */

/* PARAMETERS: */
/*   NR            --> ROW DIMENSION OF MATRIX */
/*   N       --> DIMENSION OF PROBLEM */
/*   AL(N-1,N-1)   --> LOWER TRIAGULAR MATRIX */
/*   U(N)    --> VECTOR TO FORM Z */
/*   W(N)    --> ON INPUT : VECTOR V IN SYSTEM OF LINEAR EQUATIONS */
/*               ON OUTPUT: SOLUTION OF SYSTEM OF LINEAR EQUATIONS */
/*   B(N)    --> WORKSPACE TO STORE ZT*V */

/*     FORM ZT*V (STORED IN B) */
    /* Parameter adjustments */
    al_dim1 = *nr;
    al_offset = 1 + al_dim1;
    al -= al_offset;
    --b;
    --w;
    --u;

    /* Function Body */
    d__ = ddot_(n, &u[1], &c__1, &u[1], &c__1);
    d__ /= 2.;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b[i__] = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    b[i__] += u[j] * u[i__] * w[j] / d__;
/* L15: */
	}
	b[i__] = w[i__] - b[i__];
/* L20: */
    }

/*     SOLVE LW=ZT*V */
    i__1 = *n - 1;
    forslv_(nr, &i__1, &al[al_offset], &w[1], &b[1]);
    return 0;
} /* solvew_ */

/*  ---------------------- */
/*  |     D S T A R      | */
/*  ---------------------- */
/* Subroutine */ int dstar_(integer *nr, integer *n, doublereal *u, 
	doublereal *s, doublereal *w1, doublereal *w2, doublereal *w3, 
	doublereal *sigma, doublereal *al, doublereal *d__)
{
    /* System generated locals */
    integer al_dim1, al_offset, i__1;

    /* Local variables */
    integer i__;
    doublereal utt, utu;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    doublereal temp;
    extern /* Subroutine */ int bakslv_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);


/* PURPOSE: */
/*   COMPUTE TENSOR STEP D=SIGMA*S+ZT*T(SIGMA) */

/* ------------------------------------------------------------------------ */

/* PARAMETERS: */
/*   NR            --> ROW DIMENSION OF MATRIX */
/*   N       --> DIMENSION OF PROBLEM */
/*   U(N)    --> VECTOR TO FORM Z */
/*   S(N)    --> PREVIOUS STEP */
/*   W1(N)   --> L**-1*ZT*A, WHERE A IS DESCRIBED IN SUBROUTINE SLVMDL */
/*   W2(N)   --> L**-1*ZT*SH, WHERE H IS CURRENT HESSIAN */
/*   W3(N)   --> L**-1*ZT*G, WHERE G IS CURRENT GRADIENT */
/*   SIGMA   --> SOLUTION FOR REDUCED ONE VARIABLE MODEL */
/*   AL(N-1,N-1) --> LOWER TRIANGULAR MATRIX L */
/*   D(N)    --> TENSOR STEP */

    /* Parameter adjustments */
    al_dim1 = *nr;
    al_offset = 1 + al_dim1;
    al -= al_offset;
    --d__;
    --w3;
    --w2;
    --w1;
    --s;
    --u;

    /* Function Body */
    if (*n == 1) {
	d__[1] = *sigma * s[1];
    } else {

/*     COMPUTE T(SIGMA)=-(ZTHZ)*ZT*(G+SIGMA*SH+SIGMA**2*A/2) (STORED IN D) */
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    w2[i__] = w3[i__] + *sigma * w2[i__] + w1[i__] * .5 * *sigma * *
		    sigma;
/* L10: */
	}
	i__1 = *n - 1;
	bakslv_(nr, &i__1, &al[al_offset], &d__[1], &w2[1]);
	d__[*n] = 0.;

/*     COMPUTE TENSOR STEP D=SIGMA*S+ZT*T(SIGMA) */
	utu = ddot_(n, &u[1], &c__1, &u[1], &c__1);
	utt = ddot_(n, &u[1], &c__1, &d__[1], &c__1);
	temp = utt / utu;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d__[i__] = *sigma * s[i__] - (d__[i__] - u[i__] * 2. * temp);
/* L50: */
	}
    }
    return 0;
} /* dstar_ */

/*  ---------------------- */
/*  |     M K M D L      | */
/*  ---------------------- */
/* Subroutine */ int mkmdl_(integer *nr, integer *n, doublereal *f, 
	doublereal *fp, doublereal *g, doublereal *gp, doublereal *s, 
	doublereal *h__, doublereal *alpha, doublereal *beta, doublereal *sh, 
	doublereal *a)
{
    /* System generated locals */
    integer h_dim1, h_offset, i__1, i__2;

    /* Local variables */
    integer i__, j;
    doublereal b1, b2, gs, gps, shs, sts;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);


/* PURPOSE: */
/*   FORM TENSOR MODEL */

/* ----------------------------------------------------------------------- */

/* PARAMETERS: */
/*   NR            --> ROW DIMENSION OF MATRIX */
/*   N       --> DIMENSION OF PROBLEM */
/*   F       --> CURRENT FUNCTION VALUE */
/*   FP            --> PREVIOUS FUNCTION VALUE */
/*   G(N)    --> CURRENT GRADIENT */
/*   GP(N)   --> PREVIOUS GRADIENT */
/*   S(N)    --> STEP TO PREVIOUS POINT */
/*   H(N,N)  --> HESSIAN */
/*   ALPHA      <--  SCALAR TO FORM 3RD ORDER TERM OF TENSOR MODEL */
/*   BETA       <--  SCALAR TO FORM 4TH ORDER TERM OF TENSOR MODEL */
/*   SH(N)      <--  SH */
/*   A(N)       <--  A=2*(GP-G-SH-S*BETA/(6*STS)) */


/*     COMPUTE SH */
    /* Parameter adjustments */
    --a;
    --sh;
    h_dim1 = *nr;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --s;
    --gp;
    --g;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sh[i__] = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    sh[i__] += s[j] * h__[j + i__ * h_dim1];
/* L10: */
	}
/* L20: */
    }
    gs = ddot_(n, &g[1], &c__1, &s[1], &c__1);
    gps = ddot_(n, &gp[1], &c__1, &s[1], &c__1);
    shs = ddot_(n, &sh[1], &c__1, &s[1], &c__1);
    b1 = gps - gs - shs;
    b2 = *fp - *f - gs - shs * .5;
    *alpha = b2 * 24. - b1 * 6.;
    *beta = b1 * 24. - b2 * 72.;

/*     COMPUTE A */
    sts = ddot_(n, &s[1], &c__1, &s[1], &c__1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__] = (gp[i__] - g[i__] - sh[i__] - s[i__] * *beta / (sts * 6.)) * 
		2.;
/* L50: */
    }
    return 0;
} /* mkmdl_ */

/*  ---------------------- */
/*  |     S I G M A      | */
/*  ---------------------- */
/* Subroutine */ int sigma_(doublereal *sgstar, doublereal *a, doublereal *b, 
	doublereal *c__, doublereal *d__)
{
    doublereal s1, s2, s3;
    extern /* Subroutine */ int roots_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *), 
	    sortrt_(doublereal *, doublereal *, doublereal *);


/* PURPOSE: */
/*   COMPUTE DESIRABLE ROOT OF REDUCED ONE VARIABLE EQUATION */

/* ------------------------------------------------------------------------- */

/* PARAMETERS: */
/*   SGSTAR     --> DESIRABLE ROOT */
/*   A       --> COEFFICIENT OF 3RD ORDER TERM */
/*   B       --> COEFFICIENT OF 2ND ORDER TERM */
/*   C       --> COEFFICIENT OF 1ST ORDER TERM */
/*   D       --> COEFFICIENT OF CONSTANT TERM */


/*     COMPUTE ALL THREE ROOTS */
    roots_(&s1, &s2, &s3, a, b, c__, d__);

/*     SORT ROOTS */
    sortrt_(&s1, &s2, &s3);

/*     CHOOSE DESIRABLE ROOT */
    if (*a > 0.) {
	*sgstar = s3;
	if (s2 >= 0.) {
	    *sgstar = s1;
	}
    } else {
/*       IF  A  <  0  THEN */
	*sgstar = s2;
	if (s1 > 0. || s3 < 0.) {
	    if (s1 > 0.) {
		*sgstar = s1;
	    } else {
		*sgstar = s3;
	    }
	    *a = 0.;
	}
    }
    return 0;
} /* sigma_ */

/*  ---------------------- */
/*  |     R O O T S      | */
/*  ---------------------- */
/* Subroutine */ int roots_(doublereal *s1, doublereal *s2, doublereal *s3, 
	doublereal *a, doublereal *b, doublereal *c__, doublereal *d__)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *), acos(
	    doublereal), cos(doublereal);

    /* Local variables */
    doublereal q, r__, s, t, v, a1, a2, a3, pi, temp, theta;


/* PURPOSE: */
/*   COMPUTE ROOT(S) OF 3RD ORDER EQUATION */

/* --------------------------------------------------------------------------- */

/* PARAMETERS: */
/*   S1             <--  ROOT   (IF THREE ROOTS ARE */
/*   S2             <--  ROOT    EQUAL, THEN S1=S2=S3) */
/*   S3             <--  ROOT */
/*   A       --> COEFFICIENT OF 3RD ORDER TERM */
/*   B       --> COEFFICIENT OF 2ND ORDER TERM */
/*   C       --> COEFFICIENT OF 1ST ORDER TERM */
/*   D       --> COEFFICIENT OF CONSTANT TERM */

/*     SET VALUE OF PI */
    pi = 3.141592653589793;
    a1 = *b / *a;
    a2 = *c__ / *a;
    a3 = *d__ / *a;
    q = (a2 * 3. - a1 * a1) / 9.;
    r__ = (a1 * 9. * a2 - a3 * 27. - a1 * 2. * a1 * a1) / 54.;
    v = q * q * q + r__ * r__;
    if (v > 0.) {
	s = r__ + sqrt(v);
	t = r__ - sqrt(v);
	if (t < 0.) {
	    d__1 = -t;
	    t = -pow_dd(&d__1, &c_b250);
	} else {
	    t = pow_dd(&t, &c_b250);
	}
	if (s < 0.) {
	    d__1 = -s;
	    s = -pow_dd(&d__1, &c_b250);
	} else {
	    s = pow_dd(&s, &c_b250);
	}
	*s1 = s + t - a1 / 3.;
	*s3 = *s1;
	*s2 = *s1;
    } else {
	temp = r__ / sqrt(-pow_dd(&q, &c_b384));
	theta = acos(temp);
	theta /= 3.;
	temp = sqrt(-q) * 2.;
	*s1 = temp * cos(theta) - a1 / 3.;
	*s2 = temp * cos(theta + pi * 2. / 3.) - a1 / 3.;
	*s3 = temp * cos(theta + pi * 4. / 3.) - a1 / 3.;
    }
    return 0;
} /* roots_ */

/*  ---------------------- */
/*  | S O R T R T         | */
/*  ---------------------- */
/* Subroutine */ int sortrt_(doublereal *s1, doublereal *s2, doublereal *s3)
{
    doublereal t;


/* PURPOSE: */
/*   SORT ROOTS INTO ASCENDING ORDER */

/* ----------------------------------------------------------------------------- */

/* PARAMETERS: */
/*   S1             <--> ROOT */
/*   S2             <--> ROOT */
/*   S3             <--> ROOT */

    if (*s1 > *s2) {
	t = *s1;
	*s1 = *s2;
	*s2 = t;
    }
    if (*s2 > *s3) {
	t = *s2;
	*s2 = *s3;
	*s3 = t;
    }
    if (*s1 > *s2) {
	t = *s1;
	*s1 = *s2;
	*s2 = t;
    }
    return 0;
} /* sortrt_ */

/*  ---------------------- */
/*  |  F S T O F D       | */
/*  ---------------------- */
/* Subroutine */ int fstofd_(integer *nr, integer *m, integer *n, doublereal *
	xpls, S_fp fcn, doublereal *fpls, doublereal *a, doublereal *typx, 
	doublereal *rnoise, doublereal *fhat, integer *icase, integer *nfcnt)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__, j, jp1, nm1;
    doublereal xtmpj, stepsz;

/* PURPOSE */
/* ------- */
/* FIND FIRST ORDER FORWARD FINITE DIFFERENCE APPROXIMATION "A" TO THE */
/* FIRST DERIVATIVE OF THE FUNCTION DEFINED BY THE SUBPROGRAM "FNAME" */
/* EVALUATED AT THE NEW ITERATE "XPLS". */


/* FOR OPTIMIZATION USE THIS ROUTINE TO ESTIMATE: */
/* 1) THE FIRST DERIVATIVE (GRADIENT) OF THE OPTIMIZATION FUNCTION "FCN */
/*    ANALYTIC USER ROUTINE HAS BEEN SUPPLIED; */
/* 2) THE SECOND DERIVATIVE (HESSIAN) OF THE OPTIMIZATION FUNCTION */
/*    IF NO ANALYTIC USER ROUTINE HAS BEEN SUPPLIED FOR THE HESSIAN BUT */
/*    ONE HAS BEEN SUPPLIED FOR THE GRADIENT ("FCN") AND IF THE */
/*    OPTIMIZATION FUNCTION IS INEXPENSIVE TO EVALUATE */

/* NOTE */
/* ---- */
/* _M=1 (OPTIMIZATION) ALGORITHM ESTIMATES THE GRADIENT OF THE FUNCTION */
/*      (FCN).   FCN(X) # F: R(N)-->R(1) */
/* _M=N (SYSTEMS) ALGORITHM ESTIMATES THE JACOBIAN OF THE FUNCTION */
/*      FCN(X) # F: R(N)-->R(N). */
/* _M=N (OPTIMIZATION) ALGORITHM ESTIMATES THE HESSIAN OF THE OPTIMIZATIO */
/*      FUNCTION, WHERE THE HESSIAN IS THE FIRST DERIVATIVE OF "FCN" */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* M            --> NUMBER OF ROWS IN A */
/* N            --> NUMBER OF COLUMNS IN A; DIMENSION OF PROBLEM */
/* XPLS(N)      --> NEW ITERATE:  X[K] */
/* FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION */
/* FPLS(M)      --> _M=1 (OPTIMIZATION) FUNCTION VALUE AT NEW ITERATE: */
/*                       FCN(XPLS) */
/*                  _M=N (OPTIMIZATION) VALUE OF FIRST DERIVATIVE */
/*                       (GRADIENT) GIVEN BY USER FUNCTION FCN */
/*                  _M=N (SYSTEMS)  FUNCTION VALUE OF ASSOCIATED */
/*                       MINIMIZATION FUNCTION */
/* A(NR,N)     <--  FINITE DIFFERENCE APPROXIMATION (SEE NOTE).  ONLY */
/*                  LOWER TRIANGULAR MATRIX AND DIAGONAL ARE RETURNED */
/* RNOISE       --> RELATIVE NOISE IN FCN [F(X)] */
/* FHAT(M)      --> WORKSPACE */
/* ICASE        --> =1 OPTIMIZATION (GRADIENT) */
/*                  =2 SYSTEMS */
/*                  =3 OPTIMIZATION (HESSIAN) */
/* NFCNT       <--> NUMBER OF FUNCTION EVALUTIONS */

/* INTERNAL VARIABLES */
/* ------------------ */
/* STEPSZ - STEPSIZE IN THE J-TH VARIABLE DIRECTION */


/* FIND J-TH COLUMN OF A */
/* EACH COLUMN IS DERIVATIVE OF F(FCN) WITH RESPECT TO XPLS(J) */

    /* Parameter adjustments */
    --fhat;
    --fpls;
    --typx;
    a_dim1 = *nr;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --xpls;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	xtmpj = xpls[j];
/* Computing MAX */
	d__2 = (d__1 = xpls[j], abs(d__1));
	stepsz = sqrt(*rnoise) * max(d__2,1.);
	xpls[j] = xtmpj + stepsz;
	(*fcn)(n, &xpls[1], &fhat[1]);
	++(*nfcnt);
	xpls[j] = xtmpj;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    a[i__ + j * a_dim1] = (fhat[i__] - fpls[i__]) / stepsz;
	    a[i__ + j * a_dim1] *= typx[j];
/* L20: */
	}
/* L30: */
    }
    if (*icase != 3) {
	return 0;
    }

/* IF COMPUTING HESSIAN, A MUST BE SYMMETRIC */

    if (*n == 1) {
	return 0;
    }
    nm1 = *n - 1;
    i__1 = nm1;
    for (j = 1; j <= i__1; ++j) {
	jp1 = j + 1;
	i__2 = *m;
	for (i__ = jp1; i__ <= i__2; ++i__) {
	    a[i__ + j * a_dim1] = (a[i__ + j * a_dim1] + a[j + i__ * a_dim1]) 
		    / 2.f;
/* L40: */
	}
/* L50: */
    }
    return 0;
} /* fstofd_ */

/*  ---------------------- */
/*  |  S N D O F D       | */
/*  ---------------------- */
/* Subroutine */ int sndofd_(integer *nr, integer *n, doublereal *xpls, S_fp 
	fcn, doublereal *fpls, doublereal *a, doublereal *typx, doublereal *
	rnoise, doublereal *stepsz, doublereal *anbr, integer *nfcnt)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    integer i__, j, ip1;
    doublereal ov3, fhat, xtmpi, xtmpj;

/* PURPOSE */
/* ------- */
/* FIND SECOND ORDER FORWARD FINITE DIFFERENCE APPROXIMATION "A" */
/* TO THE SECOND DERIVATIVE (HESSIAN) OF THE FUNCTION DEFINED BY THE SUBP */
/* "FCN" EVALUATED AT THE NEW ITERATE "XPLS" */

/* FOR OPTIMIZATION USE THIS ROUTINE TO ESTIMATE */
/* 1) THE SECOND DERIVATIVE (HESSIAN) OF THE OPTIMIZATION FUNCTION */
/*    IF NO ANALYTICAL USER FUNCTION HAS BEEN SUPPLIED FOR EITHER */
/*    THE GRADIENT OR THE HESSIAN AND IF THE OPTIMIZATION FUNCTION */
/*    "FCN" IS INEXPENSIVE TO EVALUATE. */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* XPLS(N)      --> NEW ITERATE:   X[K] */
/* FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION */
/* FPLS         --> FUNCTION VALUE AT NEW ITERATE, F(XPLS) */
/* A(N,N)      <--  FINITE DIFFERENCE APPROXIMATION TO HESSIAN */
/*                  ONLY LOWER TRIANGULAR MATRIX AND DIAGONAL */
/*                  ARE RETURNED */
/* RNOISE       --> RELATIVE NOISE IN FNAME [F(X)] */
/* STEPSZ(N)    --> WORKSPACE (STEPSIZE IN I-TH COMPONENT DIRECTION) */
/* ANBR(N)      --> WORKSPACE (NEIGHBOR IN I-TH DIRECTION) */
/* NFCNT       <--> NUMBER OF FUNCTION EVALUATIONS */



/* FIND I-TH STEPSIZE AND EVALUATE NEIGHBOR IN DIRECTION */
/* OF I-TH UNIT VECTOR. */

    /* Parameter adjustments */
    --anbr;
    --stepsz;
    --typx;
    a_dim1 = *nr;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --xpls;

    /* Function Body */
    ov3 = .33333333333333331;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xtmpi = xpls[i__];
/* Computing MAX */
	d__2 = (d__1 = xpls[i__], abs(d__1));
	stepsz[i__] = pow_dd(rnoise, &ov3) * max(d__2,1.);
	xpls[i__] = xtmpi + stepsz[i__];
	(*fcn)(n, &xpls[1], &anbr[i__]);
	++(*nfcnt);
	xpls[i__] = xtmpi;
/* L10: */
    }

/* CALCULATE COLUMN I OF A */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xtmpi = xpls[i__];
	xpls[i__] = xtmpi + stepsz[i__] * 2.f;
	(*fcn)(n, &xpls[1], &fhat);
	++(*nfcnt);
	a[i__ + i__ * a_dim1] = (*fpls - anbr[i__] + (fhat - anbr[i__])) / (
		stepsz[i__] * stepsz[i__]);
	a[i__ + i__ * a_dim1] *= typx[i__] * typx[i__];

/* CALCULATE SUB-DIAGONAL ELEMENTS OF COLUMN */
	if (i__ == *n) {
	    goto L25;
	}
	xpls[i__] = xtmpi + stepsz[i__];
	ip1 = i__ + 1;
	i__2 = *n;
	for (j = ip1; j <= i__2; ++j) {
	    xtmpj = xpls[j];
	    xpls[j] = xtmpj + stepsz[j];
	    (*fcn)(n, &xpls[1], &fhat);
	    ++(*nfcnt);
	    a[j + i__ * a_dim1] = (*fpls - anbr[i__] + (fhat - anbr[j])) / (
		    stepsz[i__] * stepsz[j]);
	    a[j + i__ * a_dim1] *= typx[i__] * typx[j];
	    xpls[j] = xtmpj;
/* L20: */
	}
L25:
	xpls[i__] = xtmpi;
/* L30: */
    }
    return 0;
} /* sndofd_ */

/*  ---------------------- */
/*  |  B A K S L V       | */
/*  ---------------------- */
/* Subroutine */ int bakslv_(integer *nr, integer *n, doublereal *a, 
	doublereal *x, doublereal *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    integer i__, j, ip1;
    doublereal sum;


/* PURPOSE */
/* ------- */
/* SOLVE  AX=B  WHERE A IS UPPER TRIANGULAR MATRIX. */
/* NOTE THAT A IS INPUT AS A LOWER TRIANGULAR MATRIX AND */
/* THAT THIS ROUTINE TAKES ITS TRANSPOSE IMPLICITLY. */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* A(N,N)       --> LOWER TRIANGULAR MATRIX (PRESERVED) */
/* X(N)        <--  SOLUTION VECTOR */
/* B(N)         --> RIGHT-HAND SIDE VECTOR */

/* NOTE */
/* ---- */
/* IF B IS NO LONGER REQUIRED BY CALLING ROUTINE, */
/* THEN VECTORS B AND X MAY SHARE THE SAME STORAGE. */


/* SOLVE (L-TRANSPOSE)X=B. (BACK SOLVE) */

    /* Parameter adjustments */
    --b;
    --x;
    a_dim1 = *nr;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__ = *n;
    x[i__] = b[i__] / a[i__ + i__ * a_dim1];
    if (*n == 1) {
	return 0;
    }
L30:
    ip1 = i__;
    --i__;
    sum = 0.f;
    i__1 = *n;
    for (j = ip1; j <= i__1; ++j) {
	sum += a[j + i__ * a_dim1] * x[j];
/* L40: */
    }
    x[i__] = (b[i__] - sum) / a[i__ + i__ * a_dim1];
    if (i__ > 1) {
	goto L30;
    }
    return 0;
} /* bakslv_ */

/*  ---------------------- */
/*  |  F O R S L V       | */
/*  ---------------------- */
/* Subroutine */ int forslv_(integer *nr, integer *n, doublereal *a, 
	doublereal *x, doublereal *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, im1;
    doublereal sum;


/* PURPOSE */
/* -------- */
/* SOLVE  AX=B  WHERE A  IS LOWER TRIANGULAR  MATRIX */

/* PARAMETERS */
/* --------- */

/* NR            -----> ROW DIMENSION OF MATRIX */
/* N             -----> DIMENSION OF PROBLEM */
/* A(N,N)        -----> LOWER TRIANGULAR MATRIX (PRESERVED) */
/* X(N)          <----  SOLUTION VECTOR */
/* B(N)           ----> RIGHT-HAND SIDE VECTOR */

/* NOTE */
/* ----- */
/* THEN VECTORS B AND X MAY SHARE THE SAME STORAGE */


/* SOLVE LX=B.  (FOREWARD  SOLVE) */

    /* Parameter adjustments */
    --b;
    --x;
    a_dim1 = *nr;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    x[1] = b[1] / a[a_dim1 + 1];
    if (*n == 1) {
	return 0;
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	sum = 0.f;
	im1 = i__ - 1;
	i__2 = im1;
	for (j = 1; j <= i__2; ++j) {
	    sum += a[i__ + j * a_dim1] * x[j];
/* L10: */
	}
	x[i__] = (b[i__] - sum) / a[i__ + i__ * a_dim1];
/* L20: */
    }
    return 0;
} /* forslv_ */

/*  ---------------------- */
/*  |  C H O L D R       | */
/*  ---------------------- */
/* Subroutine */ int choldr_(integer *nr, integer *n, doublereal *h__, 
	doublereal *g, doublereal *eps, integer *pivot, doublereal *e, 
	doublereal *diag, doublereal *addmax)
{
    /* System generated locals */
    integer h_dim1, h_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    integer i__, j, k;
    doublereal tau1, tau2;
    logical redo;
    doublereal temp;
    extern /* Subroutine */ int modchl_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *);


/* PURPOSE: */
/*   DRIVER FOR CHOLESKY DECOMPOSITION */

/* ---------------------------------------------------------------------- */

/* PARAMETERS: */

/*   NR            --> ROW DIMENSION */
/*   N       --> DIMENSION OF PROBLEM */
/*   H(N,N)  --> MATRIX */
/*   G(N)    --> WORK SPACE */
/*   EPS           --> MACHINE EPSILON */
/*   PIVOT(N)      --> PIVOTING VECTOR */
/*   E(N)    --> DIAGONAL MATRIX ADDED TO H FOR MAKING H P.D. */
/*   DIAG(N) --> DIAGONAL OF H */
/*   ADDMAX  --> ADDMAX * I  IS ADDED TO H */
    /* Parameter adjustments */
    --diag;
    --e;
    --pivot;
    --g;
    h_dim1 = *nr;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;

    /* Function Body */
    redo = FALSE_;

/*     SAVE DIAGONAL OF H */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	diag[i__] = h__[i__ + i__ * h_dim1];
/* L10: */
    }
    tau1 = pow_dd(eps, &c_b250);
    tau2 = tau1;
    modchl_(nr, n, &h__[h_offset], &g[1], eps, &tau1, &tau2, &pivot[1], &e[1])
	    ;
    *addmax = e[*n];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (pivot[i__] != i__) {
	    redo = TRUE_;
	}
/* L22: */
    }
    if (*addmax > 0. || redo) {
/* ******************************** */
/*                       * */
/*       H IS NOT P.D.         * */
/*                       * */
/* ******************************** */

/*     H=H+UI */
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		h__[i__ + j * h_dim1] = h__[j + i__ * h_dim1];
/* L32: */
	    }
/* L30: */
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    pivot[i__] = i__;
	    h__[i__ + i__ * h_dim1] = diag[i__] + *addmax;
/* L34: */
	}
/* ******************************** */
/*                       * */
/*        COMPUTE L            * */
/*                       * */
/* ******************************** */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {

/*     COMPUTE L(J,J) */
	    temp = 0.;
	    if (j > 1) {
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp += h__[j + i__ * h_dim1] * h__[j + i__ * h_dim1];
/* L41: */
		}
	    }
	    h__[j + j * h_dim1] -= temp;
	    h__[j + j * h_dim1] = sqrt(h__[j + j * h_dim1]);

/*     COMPUTE L(I,J) */
	    i__2 = *n;
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
		temp = 0.;
		if (j > 1) {
		    i__3 = j - 1;
		    for (k = 1; k <= i__3; ++k) {
			temp += h__[i__ + k * h_dim1] * h__[j + k * h_dim1];
/* L45: */
		    }
		}
		h__[i__ + j * h_dim1] = h__[j + i__ * h_dim1] - temp;
		h__[i__ + j * h_dim1] /= h__[j + j * h_dim1];
/* L43: */
	    }
/* L40: */
	}
    }
    return 0;
} /* choldr_ */

/*  ---------------------- */
/*  |  M O D C H L       | */
/*  ---------------------- */
/* ********************************************************************* */

/*       SUBROUTINE NAME: MODCHL */

/*       AUTHORS :  ELIZABETH ESKOW AND ROBERT B. SCHNABEL */

/*       DATE    : DECEMBER, 1988 */

/*       PURPOSE : PERFORM A MODIFIED CHOLESKY FACTORIZATION */
/*                 OF THE FORM (PTRANSPOSE)AP  + E = L(LTRANSPOSE), */
/*       WHERE L IS STORED IN THE LOWER TRIANGLE OF THE */
/*       ORIGINAL MATRIX A. */
/*       THE FACTORIZATION HAS 2 PHASES: */
/*        PHASE 1: PIVOT ON THE MAXIMUM DIAGONAL ELEMENT. */
/*            CHECK THAT THE NORMAL CHOLESKY UPDATE */
/*            WOULD RESULT IN A POSITIVE DIAGONAL */
/*            AT THE CURRENT ITERATION, AND */
/*            IF SO, DO THE NORMAL CHOLESKY UPDATE, */
/*            OTHERWISE SWITCH TO PHASE 2. */
/*        PHASE 2: PIVOT ON THE MINIMUM OF THE NEGATIVES */
/*            OF THE LOWER GERSCHGORIN BOUND */
/*            ESTIMATES. */
/*            COMPUTE THE AMOUNT TO ADD TO THE */
/*            PIVOT ELEMENT AND ADD THIS */
/*            TO THE PIVOT ELEMENT. */
/*            DO THE CHOLESKY UPDATE. */
/*            UPDATE THE ESTIMATES OF THE */
/*            GERSCHGORIN BOUNDS. */

/*       INPUT   : NDIM    - LARGEST DIMENSION OF MATRIX THAT */
/*                           WILL BE USED */

/*                 N       - DIMENSION OF MATRIX A */

/*                 A       - N*N SYMMETRIC MATRIX (ONLY LOWER TRIANGULAR */
/*            PORTION OF A, INCLUDING THE MAIN DIAGONAL, IS USED) */

/*                 G       - N*1 WORK ARRAY */

/*                 MCHEPS - MACHINE PRECISION */

/*                TAU1    - TOLERANCE USED FOR DETERMINING WHEN TO SWITCH TO */
/*                          PHASE 2 */

/*                TAU2    - TOLERANCE USED FOR DETERMINING THE MAXIMUM */
/*                          CONDITION NUMBER OF THE FINAL 2X2 SUBMATRIX. */


/*       OUTPUT  : L     - STORED IN THE MATRIX A (IN LOWER TRIANGULAR */
/*                           PORTION OF A, INCLUDING THE MAIN DIAGONAL) */

/*                 P     - A RECORD OF HOW THE ROWS AND COLUMNS */
/*                         OF THE MATRIX WERE PERMUTED WHILE */
/*                         PERFORMING THE DECOMPOSITION */

/*                 E     - N*1 ARRAY, THE ITH ELEMENT IS THE */
/*                         AMOUNT ADDED TO THE DIAGONAL OF A */
/*                         AT THE ITH ITERATION */


/* ************************************************************************ */
/* Subroutine */ int modchl_(integer *ndim, integer *n, doublereal *a, 
	doublereal *g, doublereal *mcheps, doublereal *tau1, doublereal *tau2,
	 integer *p, doublereal *e)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__, j, k, jp1;
    doublereal maxd, ming;
    extern /* Subroutine */ int init_(integer *, integer *, doublereal *, 
	    logical *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    doublereal temp, gamma, delta, jdmin;
    integer imaxd, iming;
    extern /* Subroutine */ int fin2x2_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *)
	    ;
    doublereal tdmin;
    integer itemp;
    doublereal normj, delta1;
    logical phase1;
    extern /* Subroutine */ int gersch_(integer *, integer *, doublereal *, 
	    integer *, doublereal *);
    doublereal taugam, tempjj;


/*  J        - CURRENT ITERATION NUMBER */
/*  IMING    - INDEX OF THE ROW WITH THE MIN. OF THE */
/*           NEG. LOWER GERSCH. BOUNDS */
/*  IMAXD    - INDEX OF THE ROW WITH THE MAXIMUM DIAG. */
/*           ELEMENT */
/*  I,ITEMP,JPL,K  - TEMPORARY INTEGER VARIABLES */
/*  DELTA    - AMOUNT TO ADD TO AJJ AT THE JTH ITERATION */
/*  GAMMA    - THE MAXIMUM DIAGONAL ELEMENT OF THE ORIGINAL */
/*           MATRIX A. */
/*  NORMJ    - THE 1 NORM OF A(COLJ), ROWS J+1 --> N. */
/*  MING     - THE MINIMUM OF THE NEG. LOWER GERSCH. BOUNDS */
/*  MAXD     - THE MAXIMUM DIAGONAL ELEMENT */
/*  TAUGAM - TAU1 * GAMMA */
/*  PHASE1      - LOGICAL, TRUE IF IN PHASE1, OTHERWISE FALSE */
/*  DELTA1,TEMP,JDMIN,TDMIN,TEMPJJ - TEMPORARY DOUBLE PRECISION VARS. */

    /* Parameter adjustments */
    --e;
    --p;
    --g;
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    init_(n, ndim, &a[a_offset], &phase1, &delta, &p[1], &g[1], &e[1], &ming, 
	    tau1, &gamma, &taugam);
/*     CHECK FOR N=1 */
    if (*n == 1) {
	delta = *tau2 * (d__1 = a[a_dim1 + 1], abs(d__1)) - a[a_dim1 + 1];
	if (delta > 0.) {
	    e[1] = delta;
	}
	if (a[a_dim1 + 1] == 0.) {
	    e[1] = *tau2;
	}
	a[a_dim1 + 1] = sqrt(a[a_dim1 + 1] + e[1]);
    }


    i__1 = *n - 1;
    for (j = 1; j <= i__1; ++j) {

/*        PHASE 1 */

	if (phase1) {

/*           FIND INDEX OF MAXIMUM DIAGONAL ELEMENT A(I,I) WHERE I>=J */

	    maxd = a[j + j * a_dim1];
	    imaxd = j;
	    i__2 = *n;
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
		if (maxd < a[i__ + i__ * a_dim1]) {
		    maxd = a[i__ + i__ * a_dim1];
		    imaxd = i__;
		}
/* L20: */
	    }

/*           PIVOT TO THE TOP THE ROW AND COLUMN WITH THE MAX DIAG */

	    if (imaxd != j) {

/*              SWAP ROW J WITH ROW OF MAX DIAG */

		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = a[j + i__ * a_dim1];
		    a[j + i__ * a_dim1] = a[imaxd + i__ * a_dim1];
		    a[imaxd + i__ * a_dim1] = temp;
/* L30: */
		}

/*              SWAP COLJ AND ROW MAXDIAG BETWEEN J AND MAXDIAG */

		i__2 = imaxd - 1;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    temp = a[i__ + j * a_dim1];
		    a[i__ + j * a_dim1] = a[imaxd + i__ * a_dim1];
		    a[imaxd + i__ * a_dim1] = temp;
/* L35: */
		}

/*              SWAP COLUMN J WITH COLUMN OF MAX DIAG */

		i__2 = *n;
		for (i__ = imaxd + 1; i__ <= i__2; ++i__) {
		    temp = a[i__ + j * a_dim1];
		    a[i__ + j * a_dim1] = a[i__ + imaxd * a_dim1];
		    a[i__ + imaxd * a_dim1] = temp;
/* L40: */
		}

/*              SWAP DIAG ELEMENTS */

		temp = a[j + j * a_dim1];
		a[j + j * a_dim1] = a[imaxd + imaxd * a_dim1];
		a[imaxd + imaxd * a_dim1] = temp;

/*              SWAP ELEMENTS OF THE PERMUTATION VECTOR */

		itemp = p[j];
		p[j] = p[imaxd];
		p[imaxd] = itemp;
	    }
/*           CHECK TO SEE WHETHER THE NORMAL CHOLESKY UPDATE FOR THIS */
/*           ITERATION WOULD RESULT IN A POSITIVE DIAGONAL, */
/*           AND IF NOT THEN SWITCH TO PHASE 2. */
	    jp1 = j + 1;
	    tempjj = a[j + j * a_dim1];
	    if (tempjj > 0.) {
		jdmin = a[jp1 + jp1 * a_dim1];
		i__2 = *n;
		for (i__ = jp1; i__ <= i__2; ++i__) {
		    temp = a[i__ + j * a_dim1] * a[i__ + j * a_dim1] / tempjj;
		    tdmin = a[i__ + i__ * a_dim1] - temp;
		    jdmin = min(jdmin,tdmin);
/* L60: */
		}
		if (jdmin < taugam) {
		    phase1 = FALSE_;
		}
	    } else {
		phase1 = FALSE_;
	    }
	    if (phase1) {

/*              DO THE NORMAL CHOLESKY UPDATE IF STILL IN PHASE 1 */

		a[j + j * a_dim1] = sqrt(a[j + j * a_dim1]);
		tempjj = a[j + j * a_dim1];
		i__2 = *n;
		for (i__ = jp1; i__ <= i__2; ++i__) {
		    a[i__ + j * a_dim1] /= tempjj;
/* L70: */
		}
		i__2 = *n;
		for (i__ = jp1; i__ <= i__2; ++i__) {
		    temp = a[i__ + j * a_dim1];
		    i__3 = i__;
		    for (k = jp1; k <= i__3; ++k) {
			a[i__ + k * a_dim1] -= temp * a[k + j * a_dim1];
/* L75: */
		    }
/* L80: */
		}
		if (j == *n - 1) {
		    a[*n + *n * a_dim1] = sqrt(a[*n + *n * a_dim1]);
		}
	    } else {

/*              CALCULATE THE NEGATIVES OF THE LOWER GERSCHGORIN BOUNDS */

		gersch_(ndim, n, &a[a_offset], &j, &g[1]);
	    }
	}

/*        PHASE 2 */

	if (! phase1) {
	    if (j != *n - 1) {

/*              FIND THE MINIMUM NEGATIVE GERSHGORIN BOUND */

		iming = j;
		ming = g[j];
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    if (ming > g[i__]) {
			ming = g[i__];
			iming = i__;
		    }
/* L90: */
		}

/*               PIVOT TO THE TOP THE ROW AND COLUMN WITH THE */
/*               MINIMUM NEGATIVE GERSCHGORIN BOUND */

		if (iming != j) {

/*                  SWAP ROW J WITH ROW OF MIN GERSCH BOUND */

		    i__2 = j - 1;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp = a[j + i__ * a_dim1];
			a[j + i__ * a_dim1] = a[iming + i__ * a_dim1];
			a[iming + i__ * a_dim1] = temp;
/* L100: */
		    }

/*                  SWAP COLJ WITH ROW IMING FROM J TO IMING */

		    i__2 = iming - 1;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			temp = a[i__ + j * a_dim1];
			a[i__ + j * a_dim1] = a[iming + i__ * a_dim1];
			a[iming + i__ * a_dim1] = temp;
/* L105: */
		    }

/*                 SWAP COLUMN J WITH COLUMN OF MIN GERSCH BOUND */

		    i__2 = *n;
		    for (i__ = iming + 1; i__ <= i__2; ++i__) {
			temp = a[i__ + j * a_dim1];
			a[i__ + j * a_dim1] = a[i__ + iming * a_dim1];
			a[i__ + iming * a_dim1] = temp;
/* L110: */
		    }

/*                 SWAP DIAGONAL ELEMENTS */

		    temp = a[j + j * a_dim1];
		    a[j + j * a_dim1] = a[iming + iming * a_dim1];
		    a[iming + iming * a_dim1] = temp;

/*                 SWAP ELEMENTS OF THE PERMUTATION VECTOR */

		    itemp = p[j];
		    p[j] = p[iming];
		    p[iming] = itemp;

/*                 SWAP ELEMENTS OF THE NEGATIVE GERSCHGORIN BOUNDS VECTOR */

		    temp = g[j];
		    g[j] = g[iming];
		    g[iming] = temp;
		}

/*              CALCULATE DELTA AND ADD TO THE DIAGONAL. */
/*              DELTA=MAX{0,-A(J,J) + MAX{NORMJ,TAUGAM},DELTA_PREVIOUS} */
/*              WHERE NORMJ=SUM OF |A(I,J)|,FOR I=1,N, */
/*              DELTA_PREVIOUS IS THE DELTA COMPUTED AT THE PREVIOUS ITERATION, */
/*              AND TAUGAM IS TAU1*GAMMA. */

		normj = 0.f;
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    normj += (d__1 = a[i__ + j * a_dim1], abs(d__1));
/* L140: */
		}
		temp = max(normj,taugam);
		delta1 = temp - a[j + j * a_dim1];
		temp = 0.f;
		delta1 = max(temp,delta1);
		delta = max(delta1,delta);
		e[j] = delta;
		a[j + j * a_dim1] += e[j];

/*              UPDATE THE GERSCHGORIN BOUND ESTIMATES */
/*              (NOTE: G(I) IS THE NEGATIVE OF THE */
/*               GERSCHGORIN LOWER BOUND.) */

		if (a[j + j * a_dim1] != normj) {
		    temp = normj / a[j + j * a_dim1] - 1.f;
		    i__2 = *n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			g[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1)) * 
				temp;
/* L150: */
		    }
		}

/*              DO THE CHOLESKY UPDATE */

		a[j + j * a_dim1] = sqrt(a[j + j * a_dim1]);
		tempjj = a[j + j * a_dim1];
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    a[i__ + j * a_dim1] /= tempjj;
/* L160: */
		}
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    temp = a[i__ + j * a_dim1];
		    i__3 = i__;
		    for (k = j + 1; k <= i__3; ++k) {
			a[i__ + k * a_dim1] -= temp * a[k + j * a_dim1];
/* L170: */
		    }
/* L180: */
		}
	    } else {
		fin2x2_(ndim, n, &a[a_offset], &e[1], &j, tau2, &delta, &
			gamma);
	    }
	}
/* L200: */
    }
    return 0;
} /* modchl_ */

/* ************************************************************************ */
/*       SUBROUTINE NAME : INIT */

/*       PURPOSE : SET UP FOR START OF CHOLESKY FACTORIZATION */

/*       INPUT : N, NDIM, A, TAU1 */

/*       OUTPUT : PHASE1    - BOOLEAN VALUE SET TO TRUE IF IN PHASE ONE, */
/*             OTHERWISE FALSE. */
/*      DELTA     - AMOUNT TO ADD TO AJJ AT ITERATION J */
/*      P,G,E - DESCRIBED ABOVE IN MODCHL */
/*      MING      - THE MINIMUM NEGATIVE GERSCHGORIN BOUND */
/*      GAMMA     - THE MAXIMUM DIAGONAL ELEMENT OF A */
/*      TAUGAM  - TAU1 * GAMMA */

/* ************************************************************************ */
/* Subroutine */ int init_(integer *n, integer *ndim, doublereal *a, logical *
	phase1, doublereal *delta, integer *p, doublereal *g, doublereal *e, 
	doublereal *ming, doublereal *tau1, doublereal *gamma, doublereal *
	taugam)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    integer i__;
    extern /* Subroutine */ int gersch_(integer *, integer *, doublereal *, 
	    integer *, doublereal *);

    /* Parameter adjustments */
    --e;
    --g;
    --p;
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *phase1 = TRUE_;
    *delta = 0.f;
    *ming = 0.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__] = i__;
	g[i__] = 0.f;
	e[i__] = 0.f;
/* L10: */
    }

/*     FIND THE MAXIMUM MAGNITUDE OF THE DIAGONAL ELEMENTS. */
/*     IF ANY DIAGONAL ELEMENT IS NEGATIVE, THEN PHASE1 IS FALSE. */

    *gamma = 0.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__2 = *gamma, d__3 = (d__1 = a[i__ + i__ * a_dim1], abs(d__1));
	*gamma = max(d__2,d__3);
	if (a[i__ + i__ * a_dim1] < 0.f) {
	    *phase1 = FALSE_;
	}
/* L20: */
    }
    *taugam = *tau1 * *gamma;

/*     IF NOT IN PHASE1, THEN CALCULATE THE INITIAL GERSCHGORIN BOUNDS */
/*     NEEDED FOR THE START OF PHASE2. */

    if (! (*phase1)) {
	gersch_(ndim, n, &a[a_offset], &c__1, &g[1]);
    }
    return 0;
} /* init_ */

/* ************************************************************************ */

/*       SUBROUTINE NAME : GERSCH */

/*       PURPOSE : CALCULATE THE NEGATIVE OF THE GERSCHGORIN BOUNDS */
/*                 CALLED ONCE AT THE START OF PHASE II. */

/*       INPUT   : NDIM, N, A, J */

/*       OUTPUT  : G - AN N VECTOR CONTAINING THE NEGATIVES OF THE */
/*           GERSCHGORIN BOUNDS. */

/* ************************************************************************ */
/* Subroutine */ int gersch_(integer *ndim, integer *n, doublereal *a, 
	integer *j, doublereal *g)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    integer i__, k;
    doublereal offrow;

    /* Parameter adjustments */
    --g;
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = *j; i__ <= i__1; ++i__) {
	offrow = 0.f;
	i__2 = i__ - 1;
	for (k = *j; k <= i__2; ++k) {
	    offrow += (d__1 = a[i__ + k * a_dim1], abs(d__1));
/* L10: */
	}
	i__2 = *n;
	for (k = i__ + 1; k <= i__2; ++k) {
	    offrow += (d__1 = a[k + i__ * a_dim1], abs(d__1));
/* L20: */
	}
	g[i__] = offrow - a[i__ + i__ * a_dim1];
/* L30: */
    }
    return 0;
} /* gersch_ */

/* ************************************************************************ */

/*  SUBROUTINE NAME : FIN2X2 */

/*  PURPOSE : HANDLES FINAL 2X2 SUBMATRIX IN PHASE II. */
/*            FINDS EIGENVALUES OF FINAL 2 BY 2 SUBMATRIX, */
/*            CALCULATES THE AMOUNT TO ADD TO THE DIAGONAL, */
/*            ADDS TO THE FINAL 2 DIAGONAL ELEMENTS, */
/*            AND DOES THE FINAL UPDATE. */

/*  INPUT : NDIM, N, A, E, J, TAU2, */
/*          DELTA - AMOUNT ADDED TO THE DIAGONAL IN THE */
/*                  PREVIOUS ITERATION */

/*  OUTPUT : A - MATRIX WITH COMPLETE L FACTOR IN THE LOWER TRIANGLE, */
/*           E - N*1 VECTOR CONTAINING THE AMOUNT ADDED TO THE DIAGONAL */
/*               AT EACH ITERATION, */
/*           DELTA - AMOUNT ADDED TO DIAGONAL ELEMENTS N-1 AND N. */

/* ************************************************************************ */
/* Subroutine */ int fin2x2_(integer *ndim, integer *n, doublereal *a, 
	doublereal *e, integer *j, doublereal *tau2, doublereal *delta, 
	doublereal *gamma)
{
    /* System generated locals */
    integer a_dim1, a_offset;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal t1, t2, t3, t1a, t2a, temp, lmbd1, lmbd2, delta1, lmbdhi, 
	    lmbdlo;


/*     FIND EIGENVALUES OF FINAL 2 BY 2 SUBMATRIX */

    /* Parameter adjustments */
    --e;
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    t1 = a[*n - 1 + (*n - 1) * a_dim1] + a[*n + *n * a_dim1];
    t2 = a[*n - 1 + (*n - 1) * a_dim1] - a[*n + *n * a_dim1];
    t1a = abs(t2);
    t2a = (d__1 = a[*n + (*n - 1) * a_dim1], abs(d__1)) * 2.;
    if (t1a >= t2a) {
	if (t1a > 0.) {
	    t2a /= t1a;
	}
/* Computing 2nd power */
	d__1 = t2a;
	t3 = t1a * sqrt(d__1 * d__1 + 1.);
    } else {
	t1a /= t2a;
/* Computing 2nd power */
	d__1 = t1a;
	t3 = t2a * sqrt(d__1 * d__1 + 1.);
    }
    lmbd1 = (t1 - t3) / 2.f;
    lmbd2 = (t1 + t3) / 2.f;
    lmbdhi = max(lmbd1,lmbd2);
    lmbdlo = min(lmbd1,lmbd2);

/*     FIND DELTA SUCH THAT: */
/*     1.  THE L2 CONDITION NUMBER OF THE FINAL */
/*     2X2 SUBMATRIX + DELTA*I <= TAU2 */
/*     2. DELTA >= PREVIOUS DELTA, */
/*     3. LMBDLO + DELTA >= TAU2 * GAMMA, */
/*     WHERE LMBDLO IS THE SMALLEST EIGENVALUE OF THE FINAL */
/*     2X2 SUBMATRIX */

    delta1 = (lmbdhi - lmbdlo) / (1.f - *tau2);
    delta1 = max(delta1,*gamma);
    delta1 = *tau2 * delta1 - lmbdlo;
    temp = 0.f;
    *delta = max(*delta,temp);
    *delta = max(*delta,delta1);
    if (*delta > 0.f) {
	a[*n - 1 + (*n - 1) * a_dim1] += *delta;
	a[*n + *n * a_dim1] += *delta;
	e[*n - 1] = *delta;
	e[*n] = *delta;
    }

/*     FINAL UPDATE */

    a[*n - 1 + (*n - 1) * a_dim1] = sqrt(a[*n - 1 + (*n - 1) * a_dim1]);
    a[*n + (*n - 1) * a_dim1] /= a[*n - 1 + (*n - 1) * a_dim1];
/* Computing 2nd power */
    d__1 = a[*n + (*n - 1) * a_dim1];
    a[*n + *n * a_dim1] -= d__1 * d__1;
    a[*n + *n * a_dim1] = sqrt(a[*n + *n * a_dim1]);
    return 0;
} /* fin2x2_ */

/*  ---------------------- */
/*  |    S L V M D L     | */
/*  ---------------------- */
/* Subroutine */ int slvmdl_(integer *nr, integer *n, doublereal *h__, 
	doublereal *u, doublereal *t, doublereal *e, doublereal *diag, 
	doublereal *s, doublereal *g, integer *pivot, doublereal *w1, 
	doublereal *w2, doublereal *w3, doublereal *alpha, doublereal *beta, 
	logical *nomin, doublereal *eps)
{
    /* System generated locals */
    integer h_dim1, h_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__, j;
    doublereal r__, r1, r2, ca, cb, cc, cd, w11, sg, w22, w12, w33, w13, w23, 
	    ss, uu, shs;
    extern /* Subroutine */ int zhz_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    doublereal temp;
    extern /* Subroutine */ int sigma_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *), dstar_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    doublereal addmax;
    extern /* Subroutine */ int choldr_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *), bakslv_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    doublereal sgstar;
    extern /* Subroutine */ int forslv_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), solvew_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);


/* PURPOSE: */
/*   COMPUTE TENSOR AND NEWTON STEPS */

/* ---------------------------------------------------------------------------- */

/* PARAMETERS: */

/*   NR            --> ROW DIMENSION OF MATRIX */
/*   N       --> DIMENSION OF PROBLEM */
/*   H(N,N)  --> HESSIAN */
/*   U(N)    --> VECTOR TO FORM Q IN QR */
/*   T(N)    --> WORKSPACE */
/*   E(N)    --> DIAGONAL ADDED TO HESSIAN IN CHOLESKY DECOMPOSITION */
/*   DIAG(N) --> DIAGONAL OF HESSIAN */
/*   S(N)    --> STEP TO PREVIOUS POINT (FOR TENSOR MODEL) */
/*   G(N)    --> CURRENT GRADIENT */
/*   PIVOT(N)      --> PIVOT VECTOR FOR CHOLESKY DECOMPOSITION */
/*   W1(N)      <--> ON INPUT: A=2*(GP-G-HS-S*BETA/(6*STS)) */
/*               ON OUTPUT: TENSOR STEP */
/*   W2(N)   --> SH */
/*   W3(N)      <--  NEWTON STEP */
/*   ALPHA   --> SCALAR FOR 3RD ORDER TERM OF TENSOR MODEL */
/*   BETA    --> SCALAR FOR 4TH ORDER TERM OF TENSOR MODEL */
/*   NOMIN      <--  =.TRUE. IF TENSOR MODEL HAS NO MINIMIZER */
/*   EPS           --> MACHINE EPSILON */


/*     S O L V E    M O D E L */

    /* Parameter adjustments */
    --w3;
    --w2;
    --w1;
    --pivot;
    --g;
    --s;
    --diag;
    --e;
    --t;
    --u;
    h_dim1 = *nr;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;

    /* Function Body */
    *nomin = FALSE_;

/*     COMPUTE QTHQ(N,N), ZTHZ(N-1,N-1) = FIRST N-1 ROWS AND N-1 */
/*     COLUMNS OF QTHQ */
    if (*n > 1) {
	zhz_(nr, n, &s[1], &h__[h_offset], &u[1], &t[1]);

/*     IN CHOLESKY DECOMPOSITION WILL STORE H(1,1) ... H(N-1,N-1) */
/*     IN DIAG(1) ... DIAG(N-1), STORE H(N,N) IN DIAG(N) FIRST */
	diag[*n] = h__[*n + *n * h_dim1];

/*     COLESKY DECOMPOSITION FOR FIRST N-1 ROWS AND N-1 COLUMNS OF ZTHZ */
/*     ZTHZ(N-1,N-1)=LLT */
	i__1 = *n - 1;
	choldr_(nr, &i__1, &h__[h_offset], &t[1], eps, &pivot[1], &e[1], &
		diag[1], &addmax);
    }

/*     ON INPUT: SH IS STORED IN W2 */
    shs = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	shs += w2[i__] * s[i__];
	w3[i__] = g[i__];
/* L100: */
    }

/*   COMPUTE W1,W2,W3 */
/*   W1=L**-1*ZT*A */
/*   W2=L**-1*ZT*SH */
/*   W3=L**-1*ZT*G */

    if (*n > 1) {
	solvew_(nr, n, &h__[h_offset], &u[1], &w1[1], &t[1]);
	solvew_(nr, n, &h__[h_offset], &u[1], &w2[1], &t[1]);
	solvew_(nr, n, &h__[h_offset], &u[1], &w3[1], &t[1]);
    }

/*     COMPUTE COEFFICIENTS CA,CB,CC AND CD FOR REDUCED ONE VARIABLE */
/*     3RD ORDER EQUATION */
    w11 = 0.;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w11 += w1[i__] * w1[i__];
/* L110: */
    }
    ca = *beta / 6. - w11 / 2.;
    w12 = 0.;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w12 += w1[i__] * w2[i__];
/* L120: */
    }
    cb = *alpha / 2. - w12 * 3. / 2.;
    w13 = 0.;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w13 += w1[i__] * w3[i__];
/* L130: */
    }
    w22 = 0.;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w22 += w2[i__] * w2[i__];
/* L133: */
    }
    cc = shs - w22 - w13;
    sg = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sg += s[i__] * g[i__];
/* L140: */
    }
    w23 = 0.;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w23 += w2[i__] * w3[i__];
/* L145: */
    }
    cd = sg - w23;
    w33 = 0.;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w33 += w3[i__] * w3[i__];
/* L147: */
    }

/*     COMPUTE DESIRABLE ROOT, SGSTAR, OF 3RD ORDER EQUATION */
    if (ca != 0.) {
	sigma_(&sgstar, &ca, &cb, &cc, &cd);
	if (ca == 0.) {
	    *nomin = TRUE_;
	    goto L200;
	}
    } else {
/*       2ND ORDER ( CA=0 ) */
	if (cb != 0.) {
	    r__ = cc * cc - cb * 4. * cd;
	    if (r__ < 0.) {
		*nomin = TRUE_;
		goto L200;
	    } else {
		r1 = (-cc + sqrt(r__)) / (cb * 2.);
		r2 = (-cc - sqrt(r__)) / (cb * 2.);
		if (r2 < r1) {
		    temp = r1;
		    r1 = r2;
		    r2 = temp;
		}
		if (cb > 0.) {
		    sgstar = r2;
		} else {
		    sgstar = r1;
		}
		if (r1 > 0. && sgstar == r2 || r2 < 0. && sgstar == r1) {
		    *nomin = TRUE_;
		    goto L200;
		}
	    }
	} else {
/*         1ST ORDER (CA=0,CB=0) */
	    if (cc > 0.) {
		sgstar = -cd / cc;
	    } else {
		*nomin = TRUE_;
		goto L200;
	    }
	}
    }

/*     FIND TENSOR STEP, W1 (FUNCTION OF SGSTAR) */
    dstar_(nr, n, &u[1], &s[1], &w1[1], &w2[1], &w3[1], &sgstar, &h__[
	    h_offset], &w1[1]);

/*     COMPUTE DN */
L200:
    uu = 0.;
    ss = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	uu += u[i__] * u[i__];
	ss += s[i__] * s[i__];
/* L202: */
    }
    uu /= 2.;
    ss = sqrt(ss);
    if (*n == 1) {
	choldr_(nr, n, &h__[h_offset], &t[1], eps, &pivot[1], &e[1], &diag[1],
		 &addmax);
    } else {

/* COMPUTE LAST ROW OF L(N,N) */
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    temp = 0.;
	    if (i__ > 1) {
		i__2 = i__ - 1;
		for (j = 1; j <= i__2; ++j) {
		    temp += h__[*n + j * h_dim1] * h__[i__ + j * h_dim1];
/* L210: */
		}
	    }
	    h__[*n + i__ * h_dim1] = (h__[i__ + *n * h_dim1] - temp) / h__[
		    i__ + i__ * h_dim1];
/* L220: */
	}
	temp = 0.;
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    temp += h__[*n + i__ * h_dim1] * h__[*n + i__ * h_dim1];
/* L224: */
	}
	h__[*n + *n * h_dim1] = diag[*n] - temp + addmax;
	if (h__[*n + *n * h_dim1] > 0.) {
	    h__[*n + *n * h_dim1] = sqrt(h__[*n + *n * h_dim1]);
	} else {
/*     AFTER ADDING THE LAST COLUMN AND ROW */
/*     QTHQ IS NOT P.D., NEED TO REDO CHOLESKY DECOMPOSITION */
	    i__1 = *n;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		i__2 = i__ - 1;
		for (j = 1; j <= i__2; ++j) {
		    h__[i__ + j * h_dim1] = h__[j + i__ * h_dim1];
/* L230: */
		}
		h__[i__ + i__ * h_dim1] = diag[i__];
/* L232: */
	    }
	    h__[h_dim1 + 1] = diag[1];
	    choldr_(nr, n, &h__[h_offset], &t[1], eps, &pivot[1], &e[1], &
		    diag[1], &addmax);
	}
    }

/*   SOLVE QTHQ*QT*W3=-QT*G, WHERE W3 IS NEWTON STEP */
/*   W2=-QT*G */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w2[i__] = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    w2[i__] += u[j] * u[i__] * g[j] / uu;
/* L300: */
	}
	w2[i__] -= g[i__];
/* L302: */
    }
    forslv_(nr, n, &h__[h_offset], &w3[1], &w2[1]);
    bakslv_(nr, n, &h__[h_offset], &w2[1], &w3[1]);
/*     W2=QT*W3 => W3=Q*W2 --- NEWTON STEP */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w3[i__] = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    w3[i__] += u[i__] * u[j] * w2[j] / uu;
/* L310: */
	}
	w3[i__] = w2[i__] - w3[i__];
/* L312: */
    }
    return 0;
} /* slvmdl_ */

/*  ---------------------- */
/*  |  O P T S T P       | */
/*  ---------------------- */
/* Subroutine */ int optstp_(integer *n, doublereal *xpls, doublereal *fpls, 
	doublereal *gpls, doublereal *x, integer *itncnt, integer *icscmx, 
	integer *itrmcd, doublereal *gradtl, doublereal *steptl, doublereal *
	fscale, integer *itnlim, integer *iretcd, logical *mxtake, integer *
	ipr, integer *msg, doublereal *rgx, doublereal *rsx)
{
    /* Format strings */
    static char fmt_900[] = "(\002 OPTSTP    STEP OF MAXIMUM LENGTH (STEPMX)"
	    " TAKEN\002)";
    static char fmt_901[] = "(\002 OPTSTP    RELATIVE GRADIENT CLOSE TO ZE"
	    "RO.\002/\002 OPTSTP    CURRENT ITERATE IS PROBABLY SOLUTION.\002)"
	    ;
    static char fmt_902[] = "(\002 OPTSTP    SUCCESSIVE ITERATES WITHIN TOLE"
	    "RANCE.\002/\002 OPTSTP    CURRENT ITERATE IS PROBABLY SOLUTION"
	    ".\002)";
    static char fmt_903[] = "(\002 OPTSTP    LAST GLOBAL STEP FAILED TO LOCA"
	    "TE A POINT\002,\002 LOWER THAN X.\002/\002 OPTSTP    EITHER X IS"
	    " AN APPROXIMATE LOCAL MINIMUM\002,\002 OF THE FUNCTION,\002/\002"
	    " OPTSTP    THE FUNCTION IS TOO NON-LINEAR FOR THIS\002,\002 ALGO"
	    "RITHM,\002/\002 OPTSTP    OR STEPTL IS TOO LARGE.\002)";
    static char fmt_904[] = "(\002 OPTSTP    ITERATION LIMIT EXCEEDED.\002"
	    "/\002 OPTSTP    ALGORITHM FAILED.\002)";
    static char fmt_905[] = "(\002 OPTSTP    MAXIMUM STEP SIZE EXCEEDED 5"
	    "\002,\002 CONSECUTIVE TIMES.\002/\002 OPTSTP    EITHER THE FUNCT"
	    "ION IS UNBOUNDED BELOW,\002/\002 OPTSTP    BECOMES ASYMPTOTIC TO"
	    " A FINITE VALUE\002,\002 FROM ABOVE IN SOME DIRECTION,\002/\002 "
	    "OPTSTP    OR STEPMX IS TOO SMALL\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    doublereal d__;
    integer i__, imsg;
    doublereal relgrd;
    integer jtrmcd;
    doublereal relstp;

    /* Fortran I/O blocks */
    static cilist io___280 = { 0, 0, 0, fmt_900, 0 };
    static cilist io___281 = { 0, 0, 0, fmt_901, 0 };
    static cilist io___282 = { 0, 0, 0, fmt_902, 0 };
    static cilist io___283 = { 0, 0, 0, fmt_903, 0 };
    static cilist io___284 = { 0, 0, 0, fmt_904, 0 };
    static cilist io___285 = { 0, 0, 0, fmt_905, 0 };



/* UNCONSTRAINED MINIMIZATION STOPPING CRITERIA */
/* -------------------------------------------- */
/* FIND WHETHER THE ALGORITHM SHOULD TERMINATE, DUE TO ANY */
/* OF THE FOLLOWING: */
/* 1) PROBLEM SOLVED WITHIN USER TOLERANCE */
/* 2) CONVERGENCE WITHIN USER TOLERANCE */
/* 3) ITERATION LIMIT REACHED */
/* 4) DIVERGENCE OR TOO RESTRICTIVE MAXIMUM STEP (STEPMX) SUSPECTED */


/* PARAMETERS */
/* ---------- */
/* N            --> DIMENSION OF PROBLEM */
/* XPLS(N)      --> NEW ITERATE X[K] */
/* FPLS         --> FUNCTION VALUE AT NEW ITERATE F(XPLS) */
/* GPLS(N)      --> GRADIENT AT NEW ITERATE, G(XPLS), OR APPROXIMATE */
/* X(N)         --> OLD ITERATE X[K-1] */
/* ITNCNT       --> CURRENT ITERATION K */
/* ICSCMX      <--> NUMBER CONSECUTIVE STEPS .GE. STEPMX */
/*                  [RETAIN VALUE BETWEEN SUCCESSIVE CALLS] */
/* ITRMCD      <--  TERMINATION CODE */
/* GRADTL       --> TOLERANCE AT WHICH RELATIVE GRADIENT CONSIDERED CLOSE */
/*                  ENOUGH TO ZERO TO TERMINATE ALGORITHM */
/* STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES */
/*                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM */
/* FSCALE       --> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION */
/* ITNLIM       --> MAXIMUM NUMBER OF ALLOWABLE ITERATIONS */
/* IRETCD       --> RETURN CODE */
/* MXTAKE       --> BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED */
/* IPR          --> DEVICE TO WHICH TO SEND OUTPUT */
/* MSG         <--> CONTROL OUTPUT ON INPUT AND CONTAIN STOPPING */
/*                  CONDITION ON OUTPUT */



    /* Parameter adjustments */
    --x;
    --gpls;
    --xpls;

    /* Function Body */
    *itrmcd = 0;
    imsg = *msg;
    *rgx = 0.;
    *rsx = 0.;

/* LAST GLOBAL STEP FAILED TO LOCATE A POINT LOWER THAN X */
    if (*iretcd != 1) {
	goto L50;
    }
/*     IF(IRETCD.EQ.1) */
/*     THEN */
    jtrmcd = 3;
    goto L600;
/*     ENDIF */
L50:

/* FIND DIRECTION IN WHICH RELATIVE GRADIENT MAXIMUM. */
/* CHECK WHETHER WITHIN TOLERANCE */

/* Computing MAX */
    d__1 = abs(*fpls);
    d__ = max(d__1,*fscale);
/*     D=1D0 */
    *rgx = 0.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__3 = (d__2 = xpls[i__], abs(d__2));
	relgrd = (d__1 = gpls[i__], abs(d__1)) * max(d__3,1.) / d__;
	*rgx = max(*rgx,relgrd);
/* L100: */
    }
    jtrmcd = 1;
    if (*rgx <= *gradtl) {
	goto L600;
    }

    if (*itncnt == 0) {
	return 0;
    }

/* FIND DIRECTION IN WHICH RELATIVE STEPSIZE MAXIMUM */
/* CHECK WHETHER WITHIN TOLERANCE. */

    *rsx = 0.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__3 = (d__2 = xpls[i__], abs(d__2));
	relstp = (d__1 = xpls[i__] - x[i__], abs(d__1)) / max(d__3,1.);
	*rsx = max(*rsx,relstp);
/* L120: */
    }
    jtrmcd = 2;
    if (*rsx <= *steptl) {
	goto L600;
    }

/* CHECK ITERATION LIMIT */

    jtrmcd = 4;
    if (*itncnt >= *itnlim) {
	goto L600;
    }

/* CHECK NUMBER OF CONSECUTIVE STEPS \ STEPMX */

    if (*mxtake) {
	goto L140;
    }
/*     IF(.NOT.MXTAKE) */
/*     THEN */
    *icscmx = 0;
    return 0;
/*     ELSE */
L140:
/*       IF (MOD(MSG/8,2) .EQ. 0) WRITE(IPR,900) */
    if (imsg >= 1) {
	io___280.ciunit = *ipr;
	s_wsfe(&io___280);
	e_wsfe();
    }
    ++(*icscmx);
    if (*icscmx < 5) {
	return 0;
    }
    jtrmcd = 5;
/*     ENDIF */


/* PRINT TERMINATION CODE */

L600:
    *itrmcd = jtrmcd;
/*     IF (MOD(MSG/8,2) .EQ. 0) GO TO(601,602,603,604,605), ITRMCD */
    if (imsg >= 1) {
	switch (*itrmcd) {
	    case 1:  goto L601;
	    case 2:  goto L602;
	    case 3:  goto L603;
	    case 4:  goto L604;
	    case 5:  goto L605;
	}
    }
    goto L700;
L601:
    io___281.ciunit = *ipr;
    s_wsfe(&io___281);
    e_wsfe();
    goto L700;
L602:
    io___282.ciunit = *ipr;
    s_wsfe(&io___282);
    e_wsfe();
    goto L700;
L603:
    io___283.ciunit = *ipr;
    s_wsfe(&io___283);
    e_wsfe();
    goto L700;
L604:
    io___284.ciunit = *ipr;
    s_wsfe(&io___284);
    e_wsfe();
    goto L700;
L605:
    io___285.ciunit = *ipr;
    s_wsfe(&io___285);
    e_wsfe();

L700:
    *msg = -(*itrmcd);
    return 0;

} /* optstp_ */


/* Subroutine */ int dcopy_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, m, ix, iy, mp1;


/*     COPIES A VECTOR, X, TO A VECTOR, Y. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[iy] = dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP */

L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[i__] = dx[i__];
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 7) {
	dy[i__] = dx[i__];
	dy[i__ + 1] = dx[i__ + 1];
	dy[i__ + 2] = dx[i__ + 2];
	dy[i__ + 3] = dx[i__ + 3];
	dy[i__ + 4] = dx[i__ + 4];
	dy[i__ + 5] = dx[i__ + 5];
	dy[i__ + 6] = dx[i__ + 6];
/* L50: */
    }
    return 0;
} /* dcopy_ */


doublereal ddot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, 
	integer *incy)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    integer i__, m, ix, iy, mp1;
    doublereal dtemp;


/*     FORMS THE DOT PRODUCT OF TWO VECTORS. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += dx[ix] * dy[iy];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += dx[i__] * dy[i__];
/* L30: */
    }
    if (*n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
	dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + dx[
		i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + dx[i__ + 
		4] * dy[i__ + 4];
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* ddot_ */


/* Subroutine */ int dscal_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, m, mp1, nincx;


/*     SCALES A VECTOR BY A CONSTANT. */
/*     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENT NOT EQUAL TO 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	dx[i__] = *da * dx[i__];
/* L10: */
    }
    return 0;

/*        CODE FOR INCREMENT EQUAL TO 1 */


/*        CLEAN-UP LOOP */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	dx[i__] = *da * dx[i__];
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i__ = mp1; i__ <= i__2; i__ += 5) {
	dx[i__] = *da * dx[i__];
	dx[i__ + 1] = *da * dx[i__ + 1];
	dx[i__ + 2] = *da * dx[i__ + 2];
	dx[i__ + 3] = *da * dx[i__ + 3];
	dx[i__ + 4] = *da * dx[i__ + 4];
/* L50: */
    }
    return 0;
} /* dscal_ */

/* Main program alias */ int driver_ () { MAIN__ (); return 0; }
