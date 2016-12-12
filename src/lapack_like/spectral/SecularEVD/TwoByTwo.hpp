/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SECULAR_EVD_TWOBYTWO_HPP
#define EL_SECULAR_EVD_TWOBYTWO_HPP

namespace El {
namespace secular_evd {

// Compute an eigenvalue of the diagonal plus rank-one matrix
//
//     D + rho z z^T,
//
// where || z ||_2 = 1, with D = diag( delta_0, delta_1 ) such that 
//
//     delta_0 < delta_1
//
// and rho > 0.
//
// This routine loosely corresponds to LAPACK's {s,d}laed5 [CITATION].
//

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real TwoByTwo
( Int whichValue,
  const Real& delta0,
  const Real& delta1,
  const Real& rho,
  const Real& zeta0,
  const Real& zeta1,
        Real& delta0MinusShift,
        Real& delta1MinusShift,
  FlipOrClip negativeFix=CLIP_NEGATIVES )
{
    EL_DEBUG_CSE
    const Real zero(0), one(1), two(2);
    EL_DEBUG_ONLY(
      if( whichValue < 0 || whichValue > 1 )
          LogicError("Invalid singular value request");
      if( delta1 <= delta0 )
          LogicError("Assumption that delta0 < delta1 was broken");
      if( rho <= zero )
          LogicError("Assumption that rho > 0 was broken");
      // TODO(poulson): Check the assumption that || z ||_2 = 1
    )

    const Real diagDiff = delta1 - delta0;    
    if( whichValue == 0 )
    {
        // Find the eigenvalue in (delta_0,delta_1)

        // Determine whether to shift the origin to delta0 or delta1 by 
        // testing whether the eigenvalue occurs left or right of 
        // (delta0+delta1)/2 via testing the sign of the secular equation
        // at ((delta0+delta1)/2).
        //
        // Since  
        //
        //   omega(x) = 1 + rho (zeta_0^2/(delta_0-x) +
        //                       zeta_1^2/(delta_1-x)),
        //
        // substituting x = ((delta_0+delta_1)/2) yields
        //
        //    1 + 2 rho (zeta_1^2 - zeta_0^2) / (delta_1 - delta_0).
        //
        const Real omega = one + two*rho*(zeta1*zeta1-zeta0*zeta0)/diagDiff;

        if( omega > zero )
        {
            // The eigenvalue is closer to delta_0 than delta_1, so shifting
            // the origin to delta_0 and setting
            //
            //   eta = lambda - delta_0,
            //
            // the roots of the secular equation occur where
            //
            //   1 + rho (-zeta_0^2/eta +
            //                zeta_1^2/((delta_1 - delta_0) - eta)) = 0.
            //
            // Multiplying by eta*((delta_1 - delta_0) - eta) yield a
            // quadratic equation x^2 + b x + c = 0, with
            //
            //   b = -(delta_1 - delta_0) - rho (zeta_0^2 + zeta_1^2)
            //     = -(delta_1 - delta_0) - rho, and
            //
            //   c = rho zeta_0^2 (delta_1 - delta_0).
            //

            // Note that LAPACK's {s,d}laed5 [CITATION] claims to assume
            // that || z ||_2 = 1, but it does not actively exploit this fact
            // when computing the analogue of our 'bNeg'.
            const Real bNeg = diagDiff + rho;
            const Real c = rho*zeta0*zeta0*diagDiff;

            // We inline SolveQuadraticMinus to avoid a branch;
            // we do so to respect LAPACK's strategy, but the gain for the
            // complexity is questionable.
            Real eta;
            {
                Real discrim = bNeg*bNeg - 4*c;
                if( negativeFix == CLIP_NEGATIVES )
                    discrim = Max( discrim, zero );
                else
                    discrim = Abs( discrim );

                // Clearly b is always negative, and so we avoid cancellation in
                // the formula
                //
                //   eta = (-b - sqrt(b^2 - 4c)) / 2,
                //
                // by using the "inverted" quadratic formula,
                // 
                //   eta = 2c / (-b + sqrt(b^2 - 4c)).
                //
                eta = (2*c) / (bNeg + Sqrt(discrim));
            }

            delta0MinusShift = -eta;
            delta1MinusShift = diagDiff - eta;
            return eta + delta0;
        }
        else
        {
            // The eigenvalue is at least as close to delta_1 as to delta_0,
            // so shifting the origin to delta_1 and setting
            //
            //   eta = lambda - delta_1,
            //
            // the secular equation becomes a quadratic x^2 + b x + c = 0, with
            //
            //   b = (delta_1 - delta_0) - rho, and
            //   c = -rho zeta_1^2 (delta_1 - delta_0).
            //
            const Real bNeg = -diagDiff + rho;
            const Real c = -rho*zeta1*zeta1*diagDiff;

            const Real eta = SolveQuadraticMinus( bNeg, c, negativeFix );

            delta0MinusShift = -(diagDiff + eta);
            delta1MinusShift = -eta;
            return eta + delta1;
        }
    }
    else
    {
        // Find the eigenvalue above delta_1 by shifting the origin to 
        // delta_1 (similar to above, but with the '+' branch).
        const Real bNeg = -diagDiff + rho;
        const Real c = -rho*zeta1*zeta1*diagDiff;

        const Real eta = SolveQuadraticPlus( bNeg, c, negativeFix );

        delta0MinusShift = -(diagDiff + eta);
        delta1MinusShift = -eta;
        return eta + delta1;
    }
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real TwoByTwo
( Int whichValue,
  const Real& delta0,
  const Real& delta1,
  const Real& rho,
  const Real& zeta0,
  const Real& zeta1,
  FlipOrClip negativeFix=CLIP_NEGATIVES )
{
    EL_DEBUG_CSE
    Real delta0MinusShift, delta1MinusShift;
    return
      TwoByTwo
      ( whichValue, delta0, delta1, rho, zeta0, zeta1,
        delta0MinusShift, delta1MinusShift, negativeFix );
}

} // namespace secular_evd
} // namespace El

#endif // ifndef EL_SECULAR_EVD_TWOBYTWO_HPP
