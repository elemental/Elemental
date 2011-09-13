      DOUBLE PRECISION FUNCTION XDLAMCH( CMACH )
*
*  -- LAPACK auxiliary routine (version 3.2.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     Based on LAPACK XDLAMCH but with Fortran 95 query functions
*     See: http://www.cs.utk.edu/~luszczek/lapack/lamch.html
*     and  http://www.netlib.org/lapack-dev/lapack-coding/program-style.html#id2537289
*     July 2010
*
*  -- Modified by Jack Poulson in December 2010 to avoid symbol conflicts with
*     a full LAPACK library.
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
*  Purpose
*  =======
*
*  XDLAMCH determines double precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by XDLAMCH:
*          = 'E' or 'e',   XDLAMCH := eps
*          = 'S' or 's ,   XDLAMCH := sfmin
*          = 'B' or 'b',   XDLAMCH := base
*          = 'P' or 'p',   XDLAMCH := eps*base
*          = 'N' or 'n',   XDLAMCH := t
*          = 'R' or 'r',   XDLAMCH := rnd
*          = 'M' or 'm',   XDLAMCH := emin
*          = 'U' or 'u',   XDLAMCH := rmin
*          = 'L' or 'l',   XDLAMCH := emax
*          = 'O' or 'o',   XDLAMCH := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   RND, EPS, SFMIN, SMALL, RMACH
*     ..
*     .. External Functions ..
      LOGICAL            XLSAME
      EXTERNAL           XLSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DIGITS, EPSILON, HUGE, MAXEXPONENT,
     $                   MINEXPONENT, RADIX, TINY
*     ..
*     .. Executable Statements ..
*
*
*     Assume rounding, not chopping. Always.
*
      RND = ONE
*
      IF( ONE.EQ.RND ) THEN
         EPS = EPSILON(ZERO) * 0.5
      ELSE
         EPS = EPSILON(ZERO)
      END IF
*
      IF( XLSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( XLSAME( CMACH, 'S' ) ) THEN
         SFMIN = TINY(ZERO)
         SMALL = ONE / HUGE(ZERO)
         IF( SMALL.GE.SFMIN ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            SFMIN = SMALL*( ONE+EPS )
         END IF
         RMACH = SFMIN
      ELSE IF( XLSAME( CMACH, 'B' ) ) THEN
         RMACH = RADIX(ZERO)
      ELSE IF( XLSAME( CMACH, 'P' ) ) THEN
         RMACH = EPS * RADIX(ZERO)
      ELSE IF( XLSAME( CMACH, 'N' ) ) THEN
         RMACH = DIGITS(ZERO)
      ELSE IF( XLSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( XLSAME( CMACH, 'M' ) ) THEN
         RMACH = MINEXPONENT(ZERO)
      ELSE IF( XLSAME( CMACH, 'U' ) ) THEN
         RMACH = tiny(zero)
      ELSE IF( XLSAME( CMACH, 'L' ) ) THEN
         RMACH = MAXEXPONENT(ZERO)
      ELSE IF( XLSAME( CMACH, 'O' ) ) THEN
         RMACH = HUGE(ZERO)
      ELSE
         RMACH = ZERO
      END IF
*
      XDLAMCH = RMACH
      RETURN
*
*     End of XDLAMCH
*
      END
************************************************************************
*
      DOUBLE PRECISION FUNCTION XDLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
*
*  Purpose
*  =======
*
*  XDLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*  B       (input) DOUBLE PRECISION
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      XDLAMC3 = A + B
*
      RETURN
*
*     End of XDLAMC3
*
      END
*
************************************************************************
