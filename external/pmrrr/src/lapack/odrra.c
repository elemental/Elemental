/* dlarra.f -- translated by f2c (version 20061008) */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

/* Subroutine */ 
int odrra(int *n, double *d__, double *e, 
	double *e2, double *spltol, double *tnrm, int *nsplit, 
	 int *isplit, int *info)
{
    /* System generated locals */
    int i__1;
    double d__1, d__2;

    /* Builtin functions */
    // double sqrt(double);

    /* Local variables */
    int i__;
    double tmp1, eabs;


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  Compute the splitting points with threshold SPLTOL. */
/*  ODRRA sets any "small" off-diagonal elements to zero. */

/*  Arguments */
/*  ========= */

/*  N       (input) INT */
/*          The order of the matrix. N > 0. */

/*  D       (input) DOUBLE PRECISION array, dimension (N) */
/*          On entry, the N diagonal elements of the tridiagonal */
/*          matrix T. */

/*  E       (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On entry, the first (N-1) entries contain the subdiagonal */
/*          elements of the tridiagonal matrix T; E(N) need not be set. */
/*          On exit, the entries E( ISPLIT( I ) ), 1 <= I <= NSPLIT, */
/*          are set to zero, the other entries of E are untouched. */

/*  E2      (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On entry, the first (N-1) entries contain the SQUARES of the */
/*          subdiagonal elements of the tridiagonal matrix T; */
/*          E2(N) need not be set. */
/*          On exit, the entries E2( ISPLIT( I ) ), */
/*          1 <= I <= NSPLIT, have been set to zero */

/*  SPLTOL (input) DOUBLE PRECISION */
/*          The threshold for splitting. Two criteria can be used: */
/*          SPLTOL<0 : criterion based on absolute off-diagonal value */
/*          SPLTOL>0 : criterion that preserves relative accuracy */

/*  TNRM (input) DOUBLE PRECISION */
/*          The norm of the matrix. */

/*  NSPLIT  (output) INT */
/*          The number of blocks T splits into. 1 <= NSPLIT <= N. */

/*  ISPLIT  (output) INT array, dimension (N) */
/*          The splitting points, at which T breaks up into blocks. */
/*          The first block consists of rows/columns 1 to ISPLIT(1), */
/*          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2), */
/*          etc., and the NSPLIT-th consists of rows/columns */
/*          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N. */


/*  INFO    (output) INT */
/*          = 0:  successful exit */

/*  Further Details */
/*  =============== */

/*  Based on contributions by */
/*     Beresford Parlett, University of California, Berkeley, USA */
/*     Jim Demmel, University of California, Berkeley, USA */
/*     Inderjit Dhillon, University of Texas, Austin, USA */
/*     Osni Marques, LBNL/NERSC, USA */
/*     Christof Voemel, University of California, Berkeley, USA */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --isplit;
    --e2;
    --e;
    --d__;

    /* Function Body */
    *info = 0;
/*     Compute splitting points */
    *nsplit = 1;
    if (*spltol < 0.) {
/*        Criterion based on absolute off-diagonal value */
	tmp1 = fabs(*spltol) * *tnrm;
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    eabs = (d__1 = e[i__], fabs(d__1));
	    if (eabs <= tmp1) {
		e[i__] = 0.;
		e2[i__] = 0.;
		isplit[*nsplit] = i__;
		++(*nsplit);
	    }
/* L9: */
	}
    } else {
/*        Criterion that guarantees relative accuracy */
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    eabs = (d__1 = e[i__], fabs(d__1));
	    if (eabs <= *spltol * sqrt((d__1 = d__[i__], fabs(d__1))) * sqrt((
		    d__2 = d__[i__ + 1], fabs(d__2)))) {
		e[i__] = 0.;
		e2[i__] = 0.;
		isplit[*nsplit] = i__;
		++(*nsplit);
	    }
/* L10: */
	}
    }
    isplit[*nsplit] = *n;
    return 0;

/*     End of ODRRA */

} /* odrra_ */
