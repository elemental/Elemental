/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SPECTRAL_SVD_HPP
#define EL_SPECTRAL_SVD_HPP

namespace El {
namespace svd {

// Compute the SVD of a two-by-two upper bidiagonal matrix
// =======================================================
// See the discussion in
//
//   J. Demmel and W. Kahan,
//   "Accurate singular values of bidiagonal matrices",
//   LAPACK Working Note 3, 1990, [CITATION]
//
// (which resulted in LAPACK's {s,d}las2 and {s,d}lasv2 [CITATION]),
// for details on why this routine is important for the convergence of the
// bidiagonal QR algorithm. A sketch of the routine for computing the 2x2 SVD
// is given in the appendix of
//
//   Z. Bai and J. Demmel,
//   "Computing the Generalized Singular Value Decomposition",
//   LAPACK Working Note 46, 1993 [CITATION].
//

// The following assumes that
//
//   alpha00, alpha01, alpha11 >= 0, and
//   alpha00 >= alpha11.
//
// Given that
//
//  A^T A = |    alpha00^2,       alpha00 alpha01    | = | tau00, tau01 | = T,
//          | alpha00 alpha01, alpha01^2 + alpha11^2 |   | tau01, tau11 |
//
// a gradeschool determinant calculation reveals that the eigenvalues of T are
//
//   (1/2) ((tau00+tau11) +- sqrt((tau00-tau11) + (2 tau01)^2)).
//
// Furthermore, it seems that the secret to decoding the strategy of the LAPACK
// routines {s,d}lasv2 is to interpret all expressions in terms of the variables
//
//   xi_+ = (alpha00+alpha11)/alpha01, and xi_- = (alpha00-alpha11)/alpha01,
//
// and to make use of the identity
//
//   inv(sqrt(1+xi^2)+xi) = sqrt(1+xi^2)-xi,
//
// which holds for any real value of xi, in order to avoid subtractions.
//
// To this end, one can verify that the singular values of A are given by
//
//   (alpha01/2) abs(sqrt(1+xi_+^2) +- sqrt(1+xi_-^2))
//
// by squaring and checking that it matches the gradeschool eigenvalue formula
// for T.
//
// Further, the decomposition
//
//   | alpha00, alpha01 | = | cU, -sU | | sigmaMax,    0     | |  cV, sV |,
//   |     0,   alpha11 |   | sU,  cU | |     0,    sigmaMin | | -sV, cV |
//
// or
//
//   A = U Sigma V^T,
//
// must satisfy
//
//   (A^T A) V = V Sigma^2,
//
// and so
//
//    | tau00 - sigmaMax^2,       tau01        | | cV | = | 0 |,
//    |       tau01,        tau11 - sigmaMax^2 | | sV |   | 0 |
//
// and the first row demands that
//
//    | cV | perpendicular | alpha00^2 - sigmaMax^2 |,
//    | sV |               |     alpha00 alpha01    |
//
// which is equivalent to the requirement that
//
//    | -sV | parallel | alpha00^2 - sigmaMax^2 |,
//    |  cV |          |     alpha00 alpha01    |
//
// or
//
//    | cV | parallel |    alpha00 alpha01     |.
//    | sV |          | sigmaMax^2 - alpha00^2 |
//
// After carefully computing [cV; sV] using a variant of the above formula
// which exploits the inv(sqrt(1+xi^2)+xi) = sqrt(1+xi^2)-xi identity,
// one can make use of a simple variant of
//
//    | cU | = (1/sigmaMax) | alpha00, alpha01 | | cV |
//    | sU |                |    0,    alpha11 | | sV |
//
// to complete the 2x2 SVD.
//
// On exit, this routine produces the SVD
//
//   | alpha00, alpha01 | = | cU, -sU | | sigmaMax,    0     | |  cV, sV |.
//   |     0,   alpha11 |   | sU,  cU | |     0,    sigmaMin | | -sV, cV |
//
template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void TwoByTwoUpperStandard
( const Real& alpha00,
  const Real& alpha01,
  const Real& alpha11,
        Real& sigmaMax,
        Real& sigmaMin,
        Real& cU,
        Real& sU,
        Real& cV,
        Real& sV )
{
    EL_DEBUG_CSE
    const Real zero(0), one(1), two(2), four(4);

    // Specially handle the diagonal matrix diag([alpha00;alpha11])
    if( alpha01 == zero )
    {
        sigmaMax = alpha00;
        sigmaMin = alpha11;
        cU = cV = one;
        sU = sV = zero;
        return;
    }

    // Specially handle very large values of alpha01
    const Real epsilon = limits::Epsilon<Real>();
    if( alpha01 > alpha00 && (alpha00/alpha01) < epsilon )
    {
        sigmaMax = alpha01;

        // Evaluate sigmaMin = det(A) / sigmaMax with moderate intermediates
        if( alpha11 > one)
            sigmaMin = alpha00/(alpha01/alpha11);
        else
            sigmaMin = (alpha00/alpha01)*alpha11;

        // alpha01 is very big, so these rotations are very close to the
        // identity; we follow LAPACK's {s,d}lasv2
        cU = one;
        sU = alpha11/alpha01;
        cV = alpha00/alpha01;
        sV = one;

        return;
    }

    const Real relOffdiag = alpha01/alpha00;
    const Real relOffdiagSq = relOffdiag*relOffdiag;

    const Real diagDiff = alpha00-alpha11;

    // relDiagDiff = (alpha00-alpha11)/alpha00
    //             = (alpha01/alpha00) xi_-
    //             = relOffdiag xi_-
    const Real relDiagDiff = ( diagDiff==alpha00 ? one : diagDiff/alpha00 );
    const Real relDiagDiffSq = relDiagDiff*relDiagDiff;

    // relDiagSum = (alpha00+alpha11)/alpha00
    //            = (alpha01/alpha00) xi_+
    //            = relOffdiag xi_+
    const Real relDiagSum = two - relDiagDiff;
    const Real relDiagSumSq = relDiagSum*relDiagSum;

    // The "left" term is
    //   (alpha01/alpha00) sqrt(1 + xi_+^2) = relOffdiag*sqrt(1 + xi_+^2)
    const Real leftTerm = Sqrt(relOffdiagSq+relDiagSumSq);

    // The "right" term is
    //   (alpha01/alpha00) sqrt(1 + xi_-^2) = relOffdiag*sqrt(1 + xi_-^2).
    const Real rightTerm =
      ( diagDiff==zero ? Abs(relOffdiag) : Sqrt(relOffdiagSq+relDiagDiffSq) );

    const Real relSigmaMax = (leftTerm+rightTerm)/two;

    sigmaMax = alpha00*relSigmaMax;
    sigmaMin = alpha11/relSigmaMax;

    // Following Demmel and Kahan, we go out of our way to avoid the subtraction
    // in the naive transformation of
    //
    //   | cV | parallel |    alpha00 alpha01     |
    //   | sV |          | sigmaMax^2 - alpha00^2 |
    //
    // into
    //
    //   cV = (alpha00 alpha01) /
    //     sqrt((sigmaMax^2 - alpha00^2)^2 + (alpha00 alpha01)^2).
    //
    // After multiplying the numerator and denominator b 2/(alpha00 alpha01),
    //
    //   cV = 2 / sqrt((2(sigmaMax^2 - alpha00^2)/(alpha00 alpha01))^2 + 4),
    //
    // and we define
    //
    //   phi = 2(sigmaMax^2 - alpha00^2)/(alpha00 alpha01)
    //
    // so that
    //
    //   cV = 2 / sqrt(phi^2 + 4).
    //
    // The expression for phi can be significantly simplified in terms of our
    // previously computed quantities by factoring out (sigmaMax+alpha00) as
    // alpha00*(relSigmaMax+1):
    //
    //   phi = (2/alpha01) (relSigmaMax + 1) (sigmaMax - alpha00).
    //
    // Lastly, we make use of the inv(sqrt(1+xi^2)+xi) = sqrt(1+xi^2) - xi
    // identity via the manipulation
    //
    //  (2/alpha01) (sigmaMax - alpha00) =
    //    (sqrt(1+xi_+^2) + sqrt(1+xi_-^2)) - (xi_+ - xi_-) =
    //    1/(sqrt(1+xi_+^2) + xi_+) + 1/(sqrt(1+xi_-^2) + xi_-) =
    //    relOffdiag/(leftTerm+relDiagSum) + relOffdiag/(rightTerm+relDiagDiff).
    //
    // Further, it can easily be seen that
    //
    //  (alpha00 alpha01) / (sigmaMax^2 - alpha00^2) = phi/2
    //
    // so that
    //
    //   sV = (phi/2) cV = phi / sqrt(phi^2 + 4).
    //
    // Putting
    //
    //   psi = sqrt(phi^2 + 4),
    //
    //   cV = 2 / psi, and sV = phi / psi.
    //
    // Note that this strategy is due to Demmel and Kahan, which led to
    // LAPACK's {s,d}lasv2 [CITATION], which was somewhat explained within
    // Bai and Demmel's "Computing the generalized SVD" [CITATION].
    //

    Real phi;
    if( relOffdiagSq == zero )
    {
        // relOffDiag is very small
        phi = ( diagDiff==zero ? two : alpha01/diagDiff+relOffdiag/relDiagSum );
    }
    else
    {
        phi = (relSigmaMax+1)*
              (relOffdiag/(leftTerm+relDiagSum)+
               relOffdiag/(rightTerm+relDiagDiff));
    }
    const Real psi = Sqrt(phi*phi + four);
    cV = two / psi;
    sV = phi / psi;

    // In terms of our existing variables,
    //
    //   | cU | = (1/sigmaMax) | alpha00, alpha01 | | cV |
    //   | sU |                |    0,    alpha11 | | sV |
    //
    // becomes
    //
    //   | cU | = | (alpha00/sigmaMax) cV + (alpha01/sigmaMax) sV |
    //   | sU |   | (alpha11/sigmaMax) sV                         |
    //
    //          = | (1/relSigmaMax) cV + (relOffdiag/relSigmaMax) sV |.
    //            | ((alpha11/alpha00)/relSigmaMax) sV               |
    cU = (cV + relOffdiag*sV) / relSigmaMax;
    sU = ((alpha11/alpha00)*sV) / relSigmaMax;
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void TwoByTwoUpper
( const Real& alpha00,
  const Real& alpha01,
  const Real& alpha11,
        Real& sigmaMax,
        Real& sgnMax,
        Real& sigmaMin,
        Real& sgnMin,
        Real& cU,
        Real& sU,
        Real& cV,
        Real& sV )
{
    EL_DEBUG_CSE
    const Real alpha00Abs = Abs(alpha00);
    const Real alpha01Abs = Abs(alpha01);
    const Real alpha11Abs = Abs(alpha11);

    enum LargestEntry { ALPHA00_LARGEST, ALPHA01_LARGEST, ALPHA11_LARGEST };

    LargestEntry largest;
    const bool alpha00LargestDiag = ( alpha00Abs >= alpha11Abs );
    if( alpha00LargestDiag )
    {
        largest =
          ( alpha01Abs > alpha00Abs ? ALPHA01_LARGEST : ALPHA00_LARGEST );
        TwoByTwoUpperStandard
        ( alpha00Abs, alpha01Abs, alpha11Abs,
          sigmaMax, sigmaMin, cU, sU, cV, sV );
        // sgn(sV) = sgn(phi) = sgn(alpha01/alpha00)
        sV = sV*Sgn(alpha01,false)*Sgn(alpha00,false);
        cU = cU*Sgn(alpha00*cV+alpha01*sV,false);
        sU = sU*Sgn(alpha11,false)*Sgn(sV,false);
    }
    else
    {
        largest =
          ( alpha01Abs > alpha11Abs ? ALPHA01_LARGEST : ALPHA11_LARGEST );
        // A decomposition of
        //
        //   | |alpha11|, |alpha01| |
        //   |     0,     |alpha00| |
        //
        // needs to be transposed and composed with a [0, 1; 1, 0] similarity
        // in order to provide a decomposition of
        //
        //   | |alpha00|, |alpha01| |
        //   |     0,     |alpha11| |.
        //
        // The transpose swaps (cU,sU) with (cV,sV) and the permutation swaps
        // cU with sU and cV with sV.
        //
        TwoByTwoUpperStandard
        ( alpha11Abs, alpha01Abs, alpha00Abs,
          sigmaMax, sigmaMin, sV, cV, sU, cU );
        // sgn(cU) = sgn(phi) = sgn(alpha01/alpha11)
        cU = cU*Sgn(alpha01,false)*Sgn(alpha11,false);
        sV = sV*Sgn(alpha11*sU+alpha01*cU,false);
        cV = cV*Sgn(alpha00,false)*Sgn(cU,false);
    }

    switch( largest )
    {
    case ALPHA00_LARGEST:
        // |alpha00| = sigmaMax cU cV + sigmaMin sU sV
        sgnMax = Sgn(cU,false)*Sgn(cV,false)*Sgn(alpha00,false);
        break;
    case ALPHA01_LARGEST:
        // |alpha01| = sigmaMax cU sV - sigmaMin sU cV
        sgnMax = Sgn(cU,false)*Sgn(sV,false)*Sgn(alpha01,false);
        break;
    case ALPHA11_LARGEST:
        // |alpha11| = sigmaMax sU sV + sigmaMin cU cV
        sgnMax = Sgn(sU,false)*Sgn(sV,false)*Sgn(alpha11,false);
        break;
    }
    // Make use of sgnMin*sgnMax = sgn(det(A))
    sgnMin = Sgn(alpha00,false)*Sgn(alpha11,false)*sgnMax;
}

// The following is a simplified version of the above that does not require
// singular vector computation (in the style of {s,d}las2 [CITATION]) and
// implicitly puts the upper bidiagonal matrix into the positive "standard"
// form described above.
template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void TwoByTwoUpper
( const Real& alpha00,
  const Real& alpha01,
  const Real& alpha11,
        Real& sigmaMax,
        Real& sigmaMin )
{
    EL_DEBUG_CSE
    const Real zero(0), one(1), two(2);
    const Real alpha00Abs = Abs(alpha00);
    const Real alpha01Abs = Abs(alpha01);
    const Real alpha11Abs = Abs(alpha11);
    const Real diagMaxAbs = Max( alpha00Abs, alpha11Abs );
    const Real diagMinAbs = Min( alpha00Abs, alpha11Abs );

    if( diagMinAbs == zero )
    {
        if( diagMaxAbs == zero )
            sigmaMax = alpha01Abs;
        else
            sigmaMax = SafeNormAbs( diagMaxAbs, alpha01Abs );
        sigmaMin = zero;
    }
    else
    {
        if( alpha01Abs < diagMaxAbs )
        {
            const Real relDiagSum = one + diagMinAbs/diagMaxAbs;
            const Real relDiagDiff = (diagMaxAbs-diagMinAbs)/diagMaxAbs;
            const Real relOffdiag = alpha01Abs/diagMaxAbs;
            const Real relOffdiagSq = relOffdiag*relOffdiag;
            const Real leftTerm = Sqrt(relOffdiagSq+relDiagSum*relDiagSum);
            const Real rightTerm = Sqrt(relOffdiagSq+relDiagDiff*relDiagDiff);
            // TODO(poulson): Understand the motivation (assuming any exists)
            // for LAPACK forming the inverse in {s,d}las2 relative to the
            // equivalent computation in {s,d}lasv2.
            const Real relSigmaMaxInv = two/(leftTerm+rightTerm);
            sigmaMax = diagMaxAbs/relSigmaMaxInv;
            sigmaMin = diagMinAbs*relSigmaMaxInv;
        }
        else
        {
            const Real relOffdiagInv = diagMaxAbs/alpha01Abs;
            if( relOffdiagInv == zero )
            {
                sigmaMax = alpha01Abs;
                // Exploit the equivalent determinants
                sigmaMin = (diagMaxAbs*diagMinAbs)/sigmaMax;
            }
            else
            {
                const Real relDiagSum = one + diagMinAbs/diagMaxAbs;
                const Real relDiagDiff = (diagMaxAbs-diagMinAbs)/diagMaxAbs;
                const Real leftTerm =
                  Sqrt(one+(relDiagSum*relOffdiagInv)*
                           (relDiagSum*relOffdiagInv));
                const Real rightTerm =
                  Sqrt(one+(relDiagDiff*relOffdiagInv)*
                           (relDiagDiff*relOffdiagInv));
                // TODO(poulson): Understand why LAPACK goes out of its way
                // to form one half of this quantity in {s,d}las2
                // NOTE: This is relative to alpha01 rather than alpha00
                // as in previous instances
                const Real relSigmaMaxInv = two/(leftTerm+rightTerm);
                sigmaMax = alpha01Abs/relSigmaMaxInv;
                // Exploit the equivalent determinants
                sigmaMin = (diagMinAbs*relSigmaMaxInv)*relOffdiagInv;
            }
        }
    }
}

} // namespace svd
} // namespace El

#endif // ifndef EL_SPECTRAL_SVD_HPP
