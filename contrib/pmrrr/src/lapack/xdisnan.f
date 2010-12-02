      LOGICAL FUNCTION XDISNAN( DIN )
*
*  -- LAPACK auxiliary routine (version 3.2.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2010
*
*  -- Modified by Jack Poulson in December 2010 to avoid symbol conflicts with
*     a full LAPACK library.
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   DIN
*     ..
*
*  Purpose
*  =======
*
*  XDISNAN returns .TRUE. if its argument is NaN, and .FALSE.
*  otherwise.  To be replaced by the Fortran 2003 intrinsic in the
*  future.
*
*  Arguments
*  =========
*
*  DIN     (input) DOUBLE PRECISION
*          Input to test for NaN.
*
*  =====================================================================
*
*  .. External Functions ..
      LOGICAL XDLAISNAN
      EXTERNAL XDLAISNAN
*  ..
*  .. Executable Statements ..
      XDISNAN = XDLAISNAN(DIN,DIN)
      RETURN
      END
