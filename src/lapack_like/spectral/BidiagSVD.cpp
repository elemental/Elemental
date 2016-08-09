/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./BidiagSVD/QR.hpp"
#include "./BidiagSVD/DivideAndConquer.hpp"

namespace El {

namespace bidiag_svd {

template<typename Real>
Real APosterioriThreshold
( Int m, Int n,
  const Real& twoNorm,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE
    Real thresh = Max(m,n)*twoNorm*limits::Epsilon<Real>();
    auto tol = ctrl.tol;
    auto tolType = ctrl.tolType;
    if( tolType == ABSOLUTE_SING_VAL_TOL )
    {
        if( tol == Real(0) )
        {
            // Nothing was specified...keep the default relative tolerance
        }
        else
        {
            thresh = tol;
        }
    }
    else if( tolType == RELATIVE_TO_MAX_SING_VAL_TOL )
    {
        if( tol == Real(0) )
        {
            // Nothing was specified...keep the default relative tolerance
        }
        else
        {
            thresh = tol*twoNorm;
        }
    }
    else
    {
        // Keeping all singular values relatively accurate means that none
        // can be truncated
        thresh = 0;
    }
    return thresh;
}

} // namespace bidiag_svd

// TODO(poulson): Generalize to complex bidiagonal via unitary diagonal
// rotations
template<typename Real,typename>
BidiagSVDInfo
BidiagSVD
( UpperOrLower uplo,
  Matrix<Real>& mainDiag,
  Matrix<Real>& offDiag,
  Matrix<Real>& s,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE
    if( mainDiag.Height() != offDiag.Height() &&
        mainDiag.Height() != offDiag.Height()+1 )
        LogicError("Invalid main and superdiagonal lengths");

    const Int m = ( uplo==UPPER ? mainDiag.Height() : offDiag.Height()+1 );
    const Int n = ( uplo==UPPER ? offDiag.Height()+1 : mainDiag.Height() );
    const Int minDim = Min(m,n);
    const bool square = ( m == n );
    BidiagSVDInfo info;

    UpperOrLower newUplo = uplo;
    if( !square )
    {
        if( uplo == UPPER )
        {
            // Rotate into lower bidiagonal form via Givens from the right to
            // expose [0,0,...,1] as a member of the null space. The reduction
            // occurs in two phases: 
            //
            //   | x x       | |-> | x         | ... | x         |,
            //   |   x x     |     | x x x     |     | x x       |
            //   |     x x   |     |     x x   |     |   x x     |
            //   |       x x |     |       x x |     |     x x x |
            //
            // followed by another single Givens rotation from the right to 
            // rotate the bottom-right entry into the entry to its left.
            if( ctrl.progress )
                Output
                ("Rotating non-square upper bidiagonal to lower bidiagonal");

            // Compute the Givens rotations for reducing to lower bidiagonal
            // with an extra entry in the bottom right.
            Real c, s, rho;
            for( Int i=0; i<n-2; ++i )
            {
                // | f g | | c -s | = | f c + g s,  -f s + g c | = | rho,  0  |,
                // | 0 h | | s  c |   |    h s,         h c    |   | h s, h c |
                //
                // where f = mainDiag(i), g = offDiag(i), and h = mainDiag(i+1).
                // After application, offDiag(i) will point to the bottom-left
                // corner of the result rather than the top-right (due to the
                // switch to lower bidiagonal).
                rho = Givens( mainDiag(i), offDiag(i), c, s );
                mainDiag(i) = rho;
                offDiag(i) = s*mainDiag(i+1);
                mainDiag(i+1) *= c;
            }
            // Rotate the bottom-right entry into its left neighbor
            mainDiag(n-1) = Givens( mainDiag(n-1), offDiag(n-1), c, s );
            newUplo = LOWER;
        }
        else
        {
            // Rotate into upper bidiagonal form via Givens from the left
            // TODO(poulson)
            LogicError("Non-square lower bidiagonal not yet supported");
            newUplo = UPPER;
        }
    }

    if( newUplo == LOWER )
    {
        // We are currently square and lower bidiagonal, so apply Givens
        // rotations from the left to become square and upper bidiagonal
        if( ctrl.progress )
            Output("Rotating square lower bidiagonal to upper bidiagonal");
        Real c, s, rho;
        for( Int i=0; i<m-1; ++i )
        {
            // |  c, s | | f, 0 | = |  c f + s g, s h | = | rho, s h |,
            // | -s, c | | g, h |   | -s f + c g, c h |   | 0,   c h |
            //
            // where f = mainDiag(i), g = offDiag(i), and h = mainDiag(i+1).
            // After application, offDiag(i) will point to the upper-right
            // corner since we are switching to upper-bidiagonal structure.
            rho = Givens( mainDiag(i), offDiag(i), c, s );
            mainDiag(i) = rho;
            offDiag(i) = s*mainDiag(i+1);
            mainDiag(i+1) *= c;
        }
    }
    // In all cases, we should now be deflated to square and upper bidiagonal

    if( square )
    {
        s = mainDiag;
        info.qrInfo = bidiag_svd::QRAlg( s, offDiag, ctrl );
        Sort( s, DESCENDING );
    }
    else if( uplo == LOWER )
    {
        // We were non-square and lower bidiagonal. 
        auto offDiag0 = offDiag( IR(0,n-1), ALL );
        s = mainDiag;
        info.qrInfo = bidiag_svd::QRAlg( s, offDiag0, ctrl );
        Sort( s, DESCENDING );
    }
    else
    {
        // We were non-square and upper bidiagonal. 
        auto offDiag0 = offDiag( IR(0,m-1), ALL );
        s = mainDiag;
        info.qrInfo = bidiag_svd::QRAlg( s, offDiag0, ctrl );
        Sort( s, DESCENDING );
    }

    if( ctrl.approach == THIN_SVD )
    {
        // This should be a no-op
    }
    else if( ctrl.approach == COMPACT_SVD )
    {
        // Determine the rank
        Int rank = minDim;
        for( Int i=0; i<minDim; ++i )
        {
            if( s(i) <= Real(0) )
            {
                rank = i; 
                break;
            } 
        }
        s.Resize( rank, 0 );
    }
    else if( ctrl.approach == FULL_SVD )
    {
        // This should be a no-op
    }
    else if( ctrl.approach == PRODUCT_SVD )
    {
        LogicError("Product SVD not yet supported for bidiagonal matrices");
    }

    return info;
}

template<typename F>
BidiagSVDInfo
BidiagSVD
( UpperOrLower uplo,
  Matrix<Base<F>>& mainDiag,
  Matrix<Base<F>>& offDiag,
  Matrix<F>& U,
  Matrix<Base<F>>& s,
  Matrix<F>& V,
  const BidiagSVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    typedef Base<F> Real;
    if( mainDiag.Height() != offDiag.Height() &&
        mainDiag.Height() != offDiag.Height()+1 )
        LogicError("Invalid main and superdiagonal lengths");

    const Int m = ( uplo==UPPER ? mainDiag.Height() : offDiag.Height()+1 );
    const Int n = ( uplo==UPPER ? offDiag.Height()+1 : mainDiag.Height() );
    const Int minDim = Min(m,n);
    const bool square = ( m == n );
    BidiagSVDInfo info;
    if( !ctrl.wantU && !ctrl.wantV )
    {
        U.Resize( m, 0 );
        V.Resize( n, 0 );
        return BidiagSVD( uplo, mainDiag, offDiag, s, ctrl );
    }

    // We may need to explicitly initialize U and V as identity matrices and
    // switch to accumulating if accumulateU and/or accumulateV were false.
    auto ctrlMod = ctrl;

    Matrix<Real> cList, sList;
    UpperOrLower newUplo = uplo;
    if( !square )
    {
        if( uplo == UPPER )
        {
            // Rotate into lower bidiagonal form via Givens from the right to
            // expose [0,0,...,1] as a member of the null space. The reduction
            // occurs in two phases: 
            //
            //   | x x       | |-> | x         | ... | x         |,
            //   |   x x     |     | x x x     |     | x x       |
            //   |     x x   |     |     x x   |     |   x x     |
            //   |       x x |     |       x x |     |     x x x |
            //
            // followed by another single Givens rotation from the right to 
            // rotate the bottom-right entry into the entry to its left.
            if( ctrlMod.wantV && !ctrlMod.accumulateV )
            {
                Identity( V, n, n );
                ctrlMod.accumulateV = true;
            }
            cList.Resize( n-1, 1 );
            sList.Resize( n-1, 1 );
            if( ctrl.progress )
                Output
                ("Rotating non-square upper bidiagonal to lower bidiagonal");

            // Compute the Givens rotations for reducing to lower bidiagonal
            // with an extra entry in the bottom right.
            Real c, s, rho;
            for( Int i=0; i<n-2; ++i )
            {
                // | f g | | c -s | = | f c + g s,  -f s + g c | = | rho,  0  |,
                // | 0 h | | s  c |   |    h s,         h c    |   | h s, h c |
                //
                // where f = mainDiag(i), g = offDiag(i), and h = mainDiag(i+1).
                // After application, offDiag(i) will point to the bottom-left
                // corner of the result rather than the top-right (due to the
                // switch to lower bidiagonal).
                rho = Givens( mainDiag(i), offDiag(i), c, s );
                mainDiag(i) = rho;
                offDiag(i) = s*mainDiag(i+1);
                mainDiag(i+1) *= c;

                cList(i) = c;
                sList(i) = s;
            }
            // Rotate the bottom-right entry into its left neighbor
            mainDiag(n-2) =
              Givens( mainDiag(n-2), offDiag(n-2), cList(n-2), sList(n-2) );
            
            if( ctrl.wantV ) 
            {
                ApplyGivensSequence
                ( RIGHT, VARIABLE_GIVENS_SEQUENCE, FORWARD, cList, sList, V );
            }
            newUplo = LOWER;
        }
        else
        {
            // Rotate into upper bidiagonal form via Givens from the left
            // TODO(poulson)
            LogicError("Non-square lower bidiagonal not yet supported");
            if( ctrlMod.wantU && !ctrlMod.accumulateU )
            {
                Identity( U, m, m );
                ctrlMod.accumulateU = true;
            }
            cList.Resize( m-1, 1 );
            sList.Resize( m-1, 1 );
            newUplo = UPPER;
        }
    }

    if( newUplo == LOWER )
    {
        // We are currently square and lower bidiagonal, so apply Givens
        // rotations from the left to become square and upper bidiagonal
        if( ctrlMod.wantU && !ctrlMod.accumulateU )
        {
            Identity( U, m, m );
            ctrlMod.accumulateU = true;
        }
        cList.Resize( m-1, 1 );
        sList.Resize( m-1, 1 );
        if( ctrl.progress )
            Output("Rotating square lower bidiagonal to upper bidiagonal");
        Real c, s, rho;
        for( Int i=0; i<m-1; ++i )
        {
            // |  c, s | | f, 0 | = |  c f + s g, s h | = | rho, s h |,
            // | -s, c | | g, h |   | -s f + c g, c h |   | 0,   c h |
            //
            // where f = mainDiag(i), g = offDiag(i), and h = mainDiag(i+1).
            // After application, offDiag(i) will point to the upper-right
            // corner since we are switching to upper-bidiagonal structure.
            rho = Givens( mainDiag(i), offDiag(i), c, s );
            mainDiag(i) = rho;
            offDiag(i) = s*mainDiag(i+1);
            mainDiag(i+1) *= c;
            if( ctrl.wantU )
            {
                cList(i) = c;
                sList(i) = s;
            }
        }
        if( ctrl.wantU )
        {
            ApplyGivensSequence
            ( RIGHT, VARIABLE_GIVENS_SEQUENCE, FORWARD, cList, sList, U );
        }
    }
    // In all cases, we should now be deflated to square and upper bidiagonal

    if( square )
    {
        s = mainDiag;
        info.qrInfo = bidiag_svd::QRAlg( s, offDiag, U, V, ctrlMod );
        auto sortPairs = TaggedSort( s, DESCENDING );
        ApplyTaggedSortToEachColumn( sortPairs, s );
        if( ctrl.wantU )
            ApplyTaggedSortToEachRow( sortPairs, U );
        if( ctrl.wantV )
            ApplyTaggedSortToEachRow( sortPairs, V );
    }
    else if( uplo == LOWER )
    {
        // We were non-square and lower bidiagonal. The last column of U has
        // been deflated.
        auto offDiag0 = offDiag( IR(0,n-1), ALL );
        auto U0 = U( ALL, IR(0,n) );

        s = mainDiag;
        info.qrInfo = bidiag_svd::QRAlg( s, offDiag0, U0, V, ctrlMod );
        auto sortPairs = TaggedSort( s, DESCENDING );
        ApplyTaggedSortToEachColumn( sortPairs, s );
        if( ctrl.wantU )
            ApplyTaggedSortToEachRow( sortPairs, U0 );
        if( ctrl.wantV )
            ApplyTaggedSortToEachRow( sortPairs, V );
    }
    else
    {
        // We were non-square and upper bidiagonal. The last column of V has
        // been deflated.
        auto offDiag0 = offDiag( IR(0,m-1), ALL );
        auto V0 = V( ALL, IR(0,m) );

        s = mainDiag;
        info.qrInfo = bidiag_svd::QRAlg( s, offDiag0, U, V0, ctrlMod );
        auto sortPairs = TaggedSort( s, DESCENDING );
        ApplyTaggedSortToEachColumn( sortPairs, s );
        if( ctrl.wantU )
            ApplyTaggedSortToEachRow( sortPairs, U );
        if( ctrl.wantV )
            ApplyTaggedSortToEachRow( sortPairs, V0 );
    }

    if( ctrlMod.approach == THIN_SVD )
    {
        if( ctrl.wantU )
            U.Resize( m, minDim );
        if( ctrl.wantV )
            V.Resize( n, minDim );
    }
    else if( ctrlMod.approach == COMPACT_SVD )
    {
        // Determine the rank
        Int rank = minDim;
        for( Int i=0; i<minDim; ++i )
        {
            if( s(i) <= Real(0) )
            {
                rank = i; 
                break;
            } 
        }
        s.Resize( rank, 0 );
        if( ctrl.wantU )
            U.Resize( m, rank );
        if( ctrl.wantV )
            V.Resize( n, rank );
    }
    else if( ctrlMod.approach == FULL_SVD )
    {
        // This should be a no-op
    }
    else if( ctrlMod.approach == PRODUCT_SVD )
    {
        LogicError("Product SVD not yet supported for bidiagonal matrices");
    }

    return info;
}

#define PROTO(F) \
  template BidiagSVDInfo BidiagSVD \
  ( UpperOrLower uplo, \
    Matrix<Base<F>>& mainDiag, \
    Matrix<Base<F>>& offDiag, \
    Matrix<F>& U, \
    Matrix<Base<F>>& s, \
    Matrix<F>& V, \
    const BidiagSVDCtrl<Base<F>>& ctrl ); \
  template bidiag_svd::QRInfo bidiag_svd::QRAlg \
  ( Matrix<Base<F>>& mainDiag, \
    Matrix<Base<F>>& superDiag, \
    Matrix<F>& U, \
    Matrix<F>& V, \
    const BidiagSVDCtrl<Base<F>>& ctrl );

#define PROTO_REAL(Real) \
  PROTO(Real) \
  template BidiagSVDInfo BidiagSVD \
  ( UpperOrLower uplo, \
    Matrix<Real>& mainDiag, \
    Matrix<Real>& offDiag, \
    Matrix<Real>& s, \
    const BidiagSVDCtrl<Real>& ctrl ); \
  template bidiag_svd::QRInfo bidiag_svd::QRAlg \
  ( Matrix<Real>& mainDiag, \
    Matrix<Real>& superDiag, \
    const BidiagSVDCtrl<Real>& ctrl ); \
  template Real bidiag_svd::APosterioriThreshold \
  ( Int m, Int n, \
    const Real& twoNorm, \
    const BidiagSVDCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
