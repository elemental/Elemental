/* dlarnv.f -- translated by f2c (version 20061008) */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#define imin(a,b) ( (a) < (b) ? (a) : (b) )

/* Subroutine */ 
int odrnv(int *idist, int *iseed, int *n, double *x)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Builtin functions */
    // double log(double), sqrt(double), cos(double);

    /* Local variables */
    int i__;
    double u[128];
    int il, iv, il2;
    extern int odruv(int *, int *, double *);


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  ODRNV returns a vector of n random real numbers from a uniform or */
/*  normal distribution. */

/*  Arguments */
/*  ========= */

/*  IDIST   (input) INT */
/*          Specifies the distribution of the random numbers: */
/*          = 1:  uniform (0,1) */
/*          = 2:  uniform (-1,1) */
/*          = 3:  normal (0,1) */

/*  ISEED   (input/output) INT array, dimension (4) */
/*          On entry, the seed of the random number generator; the array */
/*          elements must be between 0 and 4095, and ISEED(4) must be */
/*          odd. */
/*          On exit, the seed is updated. */

/*  N       (input) INT */
/*          The number of random numbers to be generated. */

/*  X       (output) DOUBLE PRECISION array, dimension (N) */
/*          The generated random numbers. */

/*  Further Details */
/*  =============== */

/*  This routine calls the auxiliary routine ODRUV to generate random */
/*  real numbers from a uniform (0,1) distribution, in batches of up to */
/*  128 using vectorisable code. The Box-Muller method is used to */
/*  transform numbers from a uniform to a normal distribution. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --x;
    --iseed;

    /* Function Body */
    i__1 = *n;
    for (iv = 1; iv <= i__1; iv += 64) {
/* Computing MIN */
	i__2 = 64, i__3 = *n - iv + 1;
	il = imin(i__2,i__3);
	if (*idist == 3) {
	    il2 = il << 1;
	} else {
	    il2 = il;
	}

/*        Call ODRUV to generate IL2 numbers from a uniform (0,1) */
/*        distribution (IL2 <= LV) */

	odruv(&iseed[1], &il2, u);

	if (*idist == 1) {

/*           Copy generated numbers */

	    i__2 = il;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		x[iv + i__ - 1] = u[i__ - 1];
/* L10: */
	    }
	} else if (*idist == 2) {

/*           Convert generated numbers to uniform (-1,1) distribution */

	    i__2 = il;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		x[iv + i__ - 1] = u[i__ - 1] * 2. - 1.;
/* L20: */
	    }
	} else if (*idist == 3) {

/*           Convert generated numbers to normal (0,1) distribution */

	    i__2 = il;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		x[iv + i__ - 1] = sqrt(log(u[(i__ << 1) - 2]) * -2.) * cos(u[(
			i__ << 1) - 1] * 6.2831853071795864769252867663);
/* L30: */
	    }
	}
/* L40: */
    }
    return 0;

/*     End of ODRNV */

} /* odrnv_ */
