/* dlaisnan.f -- translated by f2c (version 20061008) */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

int odsnan(double *din1, double *din2)
{
    /* System generated locals */
    int ret_val;

/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  This routine is not for general use.  It exists solely to avoid */
/*  over-optimization in DISNAN. */

/*  ODSNAN checks for NaNs by comparing its two arguments for */
/*  inequality.  NaN is the only floating-point value where NaN != NaN */
/*  returns .TRUE.  To check for NaNs, pass the same variable as both */
/*  arguments. */

/*  A compiler must assume that the two arguments are */
/*  not the same variable, and the test will not be optimized away. */
/*  Interprocedural or whole-program optimization may delete this */
/*  test.  The ISNAN functions will be replaced by the correct */
/*  Fortran 03 intrinsic once the intrinsic is widely available. */

/*  Arguments */
/*  ========= */

/*  DIN1     (input) DOUBLE PRECISION */
/*  DIN2     (input) DOUBLE PRECISION */
/*          Two numbers to compare for inequality. */

/*  ===================================================================== */

/*  .. Executable Statements .. */
    ret_val = *din1 != *din2;
    return ret_val;
} /* odsnan_ */

