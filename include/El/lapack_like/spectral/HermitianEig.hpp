/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SPECTRAL_HERMITIAN_EIG_HPP
#define EL_SPECTRAL_HERMITIAN_EIG_HPP

namespace El {
namespace herm_eig {

// It is guaranteed that |lambda0| >= |lambda1|.
//
// Given
//
//    | alpha00, alpha01 |
//    | alpha01, alpha11 |,
//
// the eigenvalues satisfy
//
//   lambda in [(alpha00+alpha11)+-sqrt((alpha00-alpha11)^2+(2 alpha01)^2)]/2,
//
// and it is useful to perform computation in terms of the three terms:
//
//   diagSum = alpha00 + alpha11,
//   diagDiff = alpha00 - alpha11, and
//   twiceAlpha01 = 2 alpha01,
//
// so that
//
//   lambda in [diagSum +- sqrt(diagDiff^2 + twiceAlpha01^2)] / 2.
//
// The eigenvalue corresponding to the "+" term is of maximal magnitude iff
// diagSum >= 0, so our computation branches based upon this condition. The
// smaller eigenvalue is *not* directly computed using the above formula, as
// this would potentially lead to catastrophic cancellation, and so we instead
// carefully divide the determinant by the dominant eigenvalue.
//
// After carefully computing the root of the discriminant, say
//
//   discrimRoot = sqrt(diagDiff^2 + twiceAlpha01^2),
//
// and supposing that diagSum >= 0, the dominant (unit) eigenvector, say [c; s],
// satisfies
//
//  |alpha00-(diagSum+discrimRoot)/2, alpha01                        | |c|=|0|,
//  |alpha01,                         alpha11-(diagSum+discrimRoot)/2| |s| |0|
//
// which simplifies to
//
//   |diagDiff-discrimRoot,  twiceAlpha01          | |c| = |0|.
//   |twiceAlpha01,         -(diagDiff+discrimRoot)| |s|   |0|
//
// In order to avoid catastrophic cancellation when choosing between solving
// the first or second equation, we pick the second equation if diagDiff >= 0
// and the first equation otherwise.
//
// [CITATION]
//
// The branching logic described above has been carefully condensed by LAPACK's
// {s,d}laev2, and we follow suit.
//
// TODO(poulson): Decide how/if to avoid redundancy in the following two
// routines without incurring unacceptable overheads.

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void TwoByTwo
( const Real& alpha00,
  const Real& alpha01,
  const Real& alpha11,
  Real& lambda0, Real& lambda1,
  bool fullAccuracy )
{
    const Real zero(0), one(1), two(2);

    // Compute our preferred representation of the 2x2 symmetric tridiagonal
    const Real diagSum = alpha00 + alpha11;
    const Real diagDiff = alpha00 - alpha11;
    const Real diagDiffAbs = Abs(diagDiff);
    const Real twiceAlpha01Abs = Abs(2*alpha01);

    Real discrimRoot;
    if( diagDiffAbs > twiceAlpha01Abs )
    {
        const Real ratio = twiceAlpha01Abs / diagDiffAbs;
        discrimRoot = diagDiffAbs*Sqrt(one + ratio*ratio);
    }
    else if( diagDiffAbs < twiceAlpha01Abs )
    {
        const Real ratio = diagDiffAbs / twiceAlpha01Abs;
        discrimRoot = twiceAlpha01Abs*Sqrt(one + ratio*ratio);
    }
    else
    {
        discrimRoot = twiceAlpha01Abs*Sqrt(two);
    }

    // Compute the two eigenvalues
    const bool nonNegDiagSum = (diagSum >= zero);
    if( diagSum != zero )
    {
        // Compute the dominant eigenvalue
        if( nonNegDiagSum )
        {
            lambda0 = (diagSum+discrimRoot) / two;
        }
        else
        {
            lambda0 = (diagSum-discrimRoot) / two;
        }

        // Avoid the naive formula for the smaller eigenvalue to prevent
        // the possibility of catastrophic cancellation. We instead carefully
        // divide the determinant by the dominant eigenvalue, as
        //
        //     lambda1 = det(A) / lambda1
        //             = (alpha00 alpha11 - alpha01^2) / lambda1
        //             = (alpha00 alpha11) / lambda1 - alpha01^2 / lambda1.
        //
        // We then evaluate the left term as
        //
        //    (max(alpha00,alpha11)/lambda1)*min(alpha00,alpha11),
        //
        // and the right term as
        //
        //    (alpha01/lambda1)*alpha01.
        //
        Real maxDiag, minDiag;
        if( Abs(alpha00) > Abs(alpha11) )
        {
            maxDiag = alpha00;
            minDiag = alpha11;
        }
        else
        {
            maxDiag = alpha11;
            minDiag = alpha00;
        }
        // [CITATION] LAPACK's {s,d}laev2
        //
        // Execute in higher precision for full accuracy
        if( fullAccuracy )
        {
            typedef Promote<Real> PReal;
            const PReal alpha01P(alpha01);
            const PReal lambda0P(lambda0);
            const PReal maxDiagP(maxDiag);
            const PReal minDiagP(minDiag);
            const PReal lambda1P = (maxDiagP/lambda0P)*minDiagP -
                                   (alpha01P/lambda0P)*alpha01P;
            lambda1 = Real(lambda1P);
        }
        else
        {
            lambda1 = (maxDiag/lambda0)*minDiag - (alpha01/lambda0)*alpha01;
        }
    }
    else
    {
        lambda0 = discrimRoot / two;
        lambda1 = -lambda0;
    }
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void TwoByTwo
( const Real& alpha00,
  const Real& alpha01,
  const Real& alpha11,
  Real& lambda0, Real& lambda1,
  Real& c, Real& s,
  bool fullAccuracy )
{
    const Real zero(0), one(1), two(2);

    // Compute our preferred representation of the 2x2 symmetric tridiagonal
    const Real diagSum = alpha00 + alpha11;
    const Real diagDiff = alpha00 - alpha11;
    const Real diagDiffAbs = Abs(diagDiff);
    const Real twiceAlpha01 = 2*alpha01;
    const Real twiceAlpha01Abs = Abs(2*alpha01);

    Real discrimRoot;
    if( diagDiffAbs > twiceAlpha01Abs )
    {
        const Real ratio = twiceAlpha01Abs / diagDiffAbs;
        discrimRoot = diagDiffAbs*Sqrt(one + ratio*ratio);
    }
    else if( diagDiffAbs < twiceAlpha01Abs )
    {
        const Real ratio = diagDiffAbs / twiceAlpha01Abs;
        discrimRoot = twiceAlpha01Abs*Sqrt(one + ratio*ratio);
    }
    else
    {
        discrimRoot = twiceAlpha01Abs*Sqrt(two);
    }

    // Compute the two eigenvalues
    const bool nonNegDiagSum = (diagSum >= zero);
    if( diagSum != zero )
    {
        // Compute the dominant eigenvalue
        if( nonNegDiagSum )
        {
            lambda0 = (diagSum+discrimRoot) / two;
        }
        else
        {
            lambda0 = (diagSum-discrimRoot) / two;
        }

        // Avoid the naive formula for the smaller eigenvalue to prevent
        // the possibility of catastrophic cancellation. We instead carefully
        // divide the determinant by the dominant eigenvalue, as
        //
        //     lambda1 = det(A) / lambda1
        //             = (alpha00 alpha11 - alpha01^2) / lambda1
        //             = (alpha00 alpha11) / lambda1 - alpha01^2 / lambda1.
        //
        // We then evaluate the left term as
        //
        //    (max(alpha00,alpha11)/lambda1)*min(alpha00,alpha11),
        //
        // and the right term as
        //
        //    (alpha01/lambda1)*alpha01.
        //
        Real maxDiag, minDiag;
        if( Abs(alpha00) > Abs(alpha11) )
        {
            maxDiag = alpha00;
            minDiag = alpha11;
        }
        else
        {
            maxDiag = alpha11;
            minDiag = alpha00;
        }
        // [CITATION] LAPACK's {s,d}laev2
        //
        // Execute in higher precision for full accuracy
        if( fullAccuracy )
        {
            typedef Promote<Real> PReal;
            const PReal alpha01P(alpha01);
            const PReal lambda0P(lambda0);
            const PReal maxDiagP(maxDiag);
            const PReal minDiagP(minDiag);
            const PReal lambda1P = (maxDiagP/lambda0P)*minDiagP -
                                   (alpha01P/lambda0P)*alpha01P;
            lambda1 = Real(lambda1P);
        }
        else
        {
            lambda1 = (maxDiag/lambda0)*minDiag - (alpha01/lambda0)*alpha01;
        }
    }
    else
    {
        lambda0 = discrimRoot / two;
        lambda1 = -lambda0;
    }

    // Choose between solving the first or second equation of
    //
    //   (A - lambda0 I) | c | = | 0 |
    //                   | s |   | 0 |
    //
    // based upon the sign of alpha00 - alpha11 in order to avoid catastrophic
    // cancellation.
    //
    const bool nonNegDiagDiff = (diagDiff >= zero);
    Real shiftedDiagTerm;
    if( nonNegDiagDiff )
    {
        shiftedDiagTerm = diagDiff + discrimRoot;
    }
    else
    {
        shiftedDiagTerm = diagDiff - discrimRoot;
    }
    // Carefully normalize based upon whether or not |shiftedDiagTerm| is
    // greater than |2 alpha01|.
    const Real shiftedDiagAbs = Abs(shiftedDiagTerm);
    if( shiftedDiagAbs > twiceAlpha01Abs )
    {
        const Real ratio = -twiceAlpha01/shiftedDiagTerm;
        s = one / Sqrt(one + ratio*ratio);
        c = ratio*s;
    }
    else
    {
        if( twiceAlpha01Abs == zero )
        {
            c = one;
            s = zero;
        }
        else
        {
            const Real ratio = -shiftedDiagTerm/twiceAlpha01;
            c = one / Sqrt(one + ratio*ratio);
            s = ratio*c;
        }
    }
    if( nonNegDiagSum == nonNegDiagDiff )
    {
        Real tmp = c;
        c = -s;
        s = tmp;
    }
}

} // namespace herm_eig
} // namespace El

#endif // ifndef EL_SPECTRAL_HERMITIAN_EIG_HPP
