/* dlarrc.f -- translated by f2c (version 20061008) */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>


/* Subroutine */ 
int odrrc(char *jobt, int *n, double *vl, 
	double *vu, double *d__, double *e, double *pivmin, 
	int *eigcnt, int *lcnt, int *rcnt, int *info)
{
    /* System generated locals */
    int i__1;
    double d__1;

    /* Local variables */
    int i__;
    double sl, su, tmp, tmp2;
    int matt;
    extern int olsame(char *, char *);
    double lpivot, rpivot;


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  Find the number of eigenvalues of the symmetric tridiagonal matrix T */
/*  that are in the interval (VL,VU] if JOBT = 'T', and of L D L^T */
/*  if JOBT = 'L'. */

/*  Arguments */
/*  ========= */

/*  JOBT    (input) CHARACTER*1 */
/*          = 'T':  Compute Sturm count for matrix T. */
/*          = 'L':  Compute Sturm count for matrix L D L^T. */

/*  N       (input) INT */
/*          The order of the matrix. N > 0. */

/*  VL      (input) DOUBLE PRECISION */
/*  VU      (input) DOUBLE PRECISION */
/*          The lower and upper bounds for the eigenvalues. */

/*  D       (input) DOUBLE PRECISION array, dimension (N) */
/*          JOBT = 'T': The N diagonal elements of the tridiagonal matrix T. */
/*          JOBT = 'L': The N diagonal elements of the diagonal matrix D. */

/*  E       (input) DOUBLE PRECISION array, dimension (N) */
/*          JOBT = 'T': The N-1 offdiagonal elements of the matrix T. */
/*          JOBT = 'L': The N-1 offdiagonal elements of the matrix L. */

/*  PIVMIN  (input) DOUBLE PRECISION */
/*          The minimum pivot in the Sturm sequence for T. */

/*  EIGCNT  (output) INT */
/*          The number of eigenvalues of the symmetric tridiagonal matrix T */
/*          that are in the interval (VL,VU] */

/*  LCNT    (output) INT */
/*  RCNT    (output) INT */
/*          The left and right negcounts of the interval. */

/*  INFO    (output) INT */

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
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --e;
    --d__;

    /* Function Body */
    *info = 0;
    *lcnt = 0;
    *rcnt = 0;
    *eigcnt = 0;
    matt = olsame(jobt, "T");
    if (matt) {
/*        Sturm sequence count on T */
	lpivot = d__[1] - *vl;
	rpivot = d__[1] - *vu;
	if (lpivot <= 0.) {
	    ++(*lcnt);
	}
	if (rpivot <= 0.) {
	    ++(*rcnt);
	}
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = e[i__];
	    tmp = d__1 * d__1;
	    lpivot = d__[i__ + 1] - *vl - tmp / lpivot;
	    rpivot = d__[i__ + 1] - *vu - tmp / rpivot;
	    if (lpivot <= 0.) {
		++(*lcnt);
	    }
	    if (rpivot <= 0.) {
		++(*rcnt);
	    }
/* L10: */
	}
    } else {
/*        Sturm sequence count on L D L^T */
	sl = -(*vl);
	su = -(*vu);
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    lpivot = d__[i__] + sl;
	    rpivot = d__[i__] + su;
	    if (lpivot <= 0.) {
		++(*lcnt);
	    }
	    if (rpivot <= 0.) {
		++(*rcnt);
	    }
	    tmp = e[i__] * d__[i__] * e[i__];

	    tmp2 = tmp / lpivot;
	    if (tmp2 == 0.) {
		sl = tmp - *vl;
	    } else {
		sl = sl * tmp2 - *vl;
	    }

	    tmp2 = tmp / rpivot;
	    if (tmp2 == 0.) {
		su = tmp - *vu;
	    } else {
		su = su * tmp2 - *vu;
	    }
/* L20: */
	}
	lpivot = d__[*n] + sl;
	rpivot = d__[*n] + su;
	if (lpivot <= 0.) {
	    ++(*lcnt);
	}
	if (rpivot <= 0.) {
	    ++(*rcnt);
	}
    }
    *eigcnt = *rcnt - *lcnt;
    return 0;

/*     end of ODRRC */

} /* odrrc_ */
