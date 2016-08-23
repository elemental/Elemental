/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SECULAR_SVD_TWOBYTWO_HPP
#define EL_SECULAR_SVD_TWOBYTWO_HPP

namespace El {
namespace secular_svd {

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real RelativeEigenvalueToRelativeSingularValue
( const Real& relEigval, const Real& shiftSqrt )
{
    // Since relEigval = singVal^2 - shiftSqrt^2,
    //
    //   relEigval / (shiftSqrt + Sqrt(shiftSqrt^2 + relEigval)) =
    //
    //   (singVal^2 - shiftSqrt^2) / (shiftSqrt + singVal) =
    //
    //   singVal - shiftSqrt.
    //
    return relEigval / (shiftSqrt + Sqrt(shiftSqrt*shiftSqrt + relEigval));
}

// Compute a single singular value corresponding to the square-root of the 
// eigenvalue of the (two-by-two) diagonal plus rank one matrix
//
//     D^2 + rho z z^T,
//
// where || z ||_2 = 1, with D = diag( delta_0, delta_1 ) such that 
//
//     0 <= delta_0 < delta_1
//
// and rho > 0.
//
// This routine loosely corresponds to LAPACK's {s,d}lasd5 [CITATION].
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
        Real& delta0PlusShift,
        Real& delta1PlusShift,
  FlipOrClip negativeFix=CLIP_NEGATIVES )
{
    DEBUG_CSE
    const Real zero(0), one(1), two(2), three(3), four(4);
    DEBUG_ONLY(
      if( whichValue < 0 || whichValue > 1 )
          LogicError("Invalid singular value request");
      if( delta0 < zero )
          LogicError("Assumption that delta0 >= 0 was broken");
      if( delta1 <= delta0 )
          LogicError("Assumption that delta0 < delta1 was broken");
      if( rho <= zero )
          LogicError("Assumption that rho > 0 was broken");
      // TODO(poulson): Check the assumption that || z ||_2 = 1
    )

    // Compute the difference of squares to high relative accuracy
    const Real diagDiff = delta1 - delta0;    
    const Real diagSum = delta1 + delta0;
    const Real diagSqDiff = diagDiff*diagSum; // delta_1^2 - delta_0^2

    if( whichValue == 0 )
    {
        // Find the singular value in (delta_0,delta_1)

        // Determine whether to shift the origin to delta0 or delta1 by 
        // testing whether the singular value occurs left or right of 
        // (delta0+delta1)/2 via testing the sign of the secular equation
        // at ((delta0+delta1)/2)^2.
        //
        // Since  
        //
        //   omega(x) = 1 + rho (zeta_0^2/(delta_0^2-x) +
        //                       zeta_1^2/(delta_1^2-x)),
        //
        // substituting x = ((delta_0+delta_1)/2)^2 yields
        //
        //   omega(x) = 1 + (4 rho /(delta_1-delta_0)) *
        //     (zeta_1^2/(3 delta_1 + delta_0) -
        //      zeta_0^2/(3 delta_0 + delta_1)).
        //
        const Real omega = one + four*rho*(
          zeta1*zeta1 / (delta0+three*delta1) -
          zeta0*zeta0 / (three*delta0+delta1)) / diagDiff; 

        if( omega > zero )
        {
            // The singular value is closer to delta_0 than delta_1, so shifting
            // the origin to delta_0 and setting
            //
            //   eta = sigma^2 - delta_0^2,
            //
            // the roots of the secular equation occur at
            //
            //   1 + rho (-zeta_0^2/eta +
            //             zeta_1^2/(delta_1^2-delta_0^2-eta)) = 0.
            //
            // Multiplying by eta*(delta_1^2-delta_0^2-eta) yield a quadratic
            // equation x^2 + b x + c = 0, with
            //
            //   b = -(delta_1^2 - delta_0^2) - rho (zeta_0^2 + zeta_1^2)
            //     = -(delta_1^2 - delta_0^2) - rho, and
            //
            //   c = rho zeta_0^2 (delta_1^2 - delta_0^2).
            //

            // Note that LAPACK's {s,d}lasd5 [CITATION] claims to assume
            // that || z ||_2 = 1, but it does not actively exploit this fact
            // when computing the analogue of our 'bNeg'.
            const Real bNeg = diagSqDiff + rho;
            const Real c = rho*zeta0*zeta0*diagSqDiff;

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
                eta = 2*c / (bNeg + Sqrt(discrim));
            }

            const Real sigmaRel =
              RelativeEigenvalueToRelativeSingularValue( eta, delta0 );

            delta0MinusShift = -sigmaRel;
            delta1MinusShift = diagDiff - sigmaRel;
            delta0PlusShift = two*delta0 + sigmaRel;
            delta1PlusShift = (delta0 + sigmaRel) + delta1;
            return sigmaRel + delta0;
        }
        else
        {
            // The singular value is at least as close to delta_1 as to delta_0,
            // so shifting the origin to delta_1 and setting
            //
            //   eta = sigma^2 - delta_1^2,
            //
            // the secular equation becomes a quadratic x^2 + b x + c = 0, with
            //
            //   b = (delta_1^2 - delta_0^2) - rho, and
            //   c = -rho zeta_1^2 (delta_1^2 - delta_0^2).
            //
            const Real bNeg = -diagSqDiff + rho;
            const Real c = -rho*zeta1*zeta1*diagSqDiff;

            const Real eta = SolveQuadraticMinus( bNeg, c, negativeFix );

            const Real sigmaRel =
              RelativeEigenvalueToRelativeSingularValue( eta, delta1 );

            delta0MinusShift = -(diagDiff + sigmaRel);
            delta1MinusShift = -sigmaRel;
            delta0PlusShift = delta0 + sigmaRel + delta1;
            delta1PlusShift = two*delta1 + sigmaRel;
            return sigmaRel + delta1;
        }
    }
    else
    {
        // Find the singular value above delta_1 by shifting the origin to 
        // delta_1 (similar to above, but with the '+' branch).
        const Real bNeg = -diagSqDiff + rho;
        const Real c = -rho*zeta1*zeta1*diagSqDiff;

        const Real eta = SolveQuadraticPlus( bNeg, c, negativeFix );

        const Real sigmaRel =
          RelativeEigenvalueToRelativeSingularValue( eta, delta1 );

        delta0MinusShift = -(diagDiff + sigmaRel);
        delta1MinusShift = -sigmaRel;
        delta0PlusShift = delta0 + sigmaRel + delta1;
        delta1PlusShift = two*delta1 + sigmaRel;
        return sigmaRel + delta1;
    }
}

// It is not worth the hassle to have a separate implementation whose only
// difference is not computing delta{0,1}{Minus,Plus}Shift.
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
    DEBUG_CSE
    Real delta0PlusShift, delta0MinusShift, delta1PlusShift, delta1MinusShift;
    return
      TwoByTwo
      ( whichValue, delta0, delta1, rho, zeta0, zeta1,
        delta0MinusShift, delta1MinusShift, delta0PlusShift, delta1PlusShift,
        negativeFix );
}

} // namespace secular_svd
} // namespace El

#endif // ifndef EL_SECULAR_SVD_TWOBYTWO_HPP
