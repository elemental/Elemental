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
  const Matrix<Real>& mainDiagOrig,
  const Matrix<Real>& offDiagOrig,
        Matrix<Real>& s,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE

    // TODO(poulson): Provide a means of making these a view rather than a copy?
    auto mainDiag( mainDiagOrig );
    auto offDiag( offDiagOrig );

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

namespace bidiag_svd {

template<typename Real,typename=EnableIf<IsReal<Real>>>
BidiagSVDInfo
Helper
( UpperOrLower uplo,
  const Matrix<Real>& mainDiagOrig,
  const Matrix<Real>& offDiagOrig,
        Matrix<Real>& U,
        Matrix<Real>& s,
        Matrix<Real>& V,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE
    if( !ctrl.useQR && (ctrl.accumulateU || ctrl.accumulateV) )
        LogicError
        ("Only the QR algorithm currently supports U or V accumulation");

    // TODO(poulson): Provide a means of making these a view rather than a copy?
    auto mainDiag( mainDiagOrig );
    auto offDiag( offDiagOrig );

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

    // TODO(poulson): Decide if accumulating into U/V is ever worthwhile for 
    // performance reasons given that there should be O(n^2 k) to subsequently
    // compose the original and output matrices, but O(n^3) with a significant 
    // coefficient to directly accumulate.

    UpperOrLower newUplo = uplo;
    Matrix<Real> cDeflateList, sDeflateList;
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
            if( ctrl.wantV )
            {
                cDeflateList.Resize( n-1, 1 );
                sDeflateList.Resize( n-1, 1 );
            }
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

                if( ctrl.wantV )
                {
                    cDeflateList(i) = c;
                    sDeflateList(i) = s;
                }
            }
            // Rotate the bottom-right entry into its left neighbor
            mainDiag(n-2) = Givens( mainDiag(n-2), offDiag(n-2), c, s );
            if( ctrl.wantV )
            {
                cDeflateList(n-2) = c;
                sDeflateList(n-2) = s;
            }

            // If we are not accumulating, we can delay the application until
            // the end of the routine.
            if( ctrl.wantV && ctrl.accumulateV )
            {
                ApplyGivensSequence
                ( RIGHT, VARIABLE_GIVENS_SEQUENCE, FORWARD,
                  cDeflateList, sDeflateList, V );
            }
            
            newUplo = LOWER;
        }
        else
        {
            // Rotate into upper bidiagonal form via Givens from the left
            // TODO(poulson)
            LogicError("Non-square lower bidiagonal not yet supported");
            if( ctrl.wantU )
            {
                cDeflateList.Resize( m-1, 1 );
                sDeflateList.Resize( m-1, 1 );
            }
            newUplo = UPPER;
        }
    }

    Matrix<Real> cFlipList, sFlipList;
    if( newUplo == LOWER )
    {
        // We are currently square and lower bidiagonal, so apply Givens
        // rotations from the left to become square and upper bidiagonal
        if( ctrl.wantU )
        {
            cFlipList.Resize( m-1, 1 );
            sFlipList.Resize( m-1, 1 );
        }
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
                cFlipList(i) = c;
                sFlipList(i) = s;
            }
        }

        // If not accumulating, we can delay the application until the end
        if( ctrl.wantU && ctrl.accumulateU )
        {
            ApplyGivensSequence
            ( RIGHT, VARIABLE_GIVENS_SEQUENCE, FORWARD, 
              cFlipList, sFlipList, U );
        }
    }
    // In all cases, we should now be deflated to square and upper bidiagonal

    if( square )
    {
        if( ctrl.useQR )
        {
            s = mainDiag;
            info.qrInfo = bidiag_svd::QRAlg( s, offDiag, U, V, ctrl );
        }
        else
        {
            info.dcInfo =
              bidiag_svd::DivideAndConquer( mainDiag, offDiag, U, s, V, ctrl );
        }

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
        // been (at least implicitly) deflated.
        auto offDiag0 = offDiag( IR(0,n-1), ALL );
        Matrix<Real> U0;
        if( ctrl.wantU )
        {
            if( !ctrl.accumulateU )
                Identity( U, m, m ); 
            View( U0, U, ALL, IR(0,n) );
        }

        if( ctrl.useQR )
        {
            s = mainDiag;
            info.qrInfo = bidiag_svd::QRAlg( s, offDiag0, U0, V, ctrl );
        }
        else
        {
            info.dcInfo =
              bidiag_svd::DivideAndConquer
              ( mainDiag, offDiag0, U0, s, V, ctrl );
        }

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
        // been (at least implicitly) deflated.
        auto offDiag0 = offDiag( IR(0,m-1), ALL );
        Matrix<Real> V0;
        if( ctrl.wantV )
        {
            if( !ctrl.accumulateV )
                Identity( V, n, n );
            View( V0, V, ALL, IR(0,m) );
        }

        if( ctrl.useQR )
        {
            s = mainDiag;
            info.qrInfo = bidiag_svd::QRAlg( s, offDiag0, U, V0, ctrl );
        }
        else
        {
            info.dcInfo =
              bidiag_svd::DivideAndConquer
              ( mainDiag, offDiag0, U, s, V0, ctrl );
        }

        auto sortPairs = TaggedSort( s, DESCENDING );
        ApplyTaggedSortToEachColumn( sortPairs, s );
        if( ctrl.wantU )
            ApplyTaggedSortToEachRow( sortPairs, U );
        if( ctrl.wantV )
            ApplyTaggedSortToEachRow( sortPairs, V0 );
    }

    if( ctrl.approach == THIN_SVD )
    {
        if( ctrl.wantU )
            U.Resize( m, minDim );
        if( ctrl.wantV )
            V.Resize( n, minDim );
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
        if( ctrl.wantU )
            U.Resize( m, rank );
        if( ctrl.wantV )
            V.Resize( n, rank );
    }
    else if( ctrl.approach == FULL_SVD )
    {
        // This should be a no-op
    }
    else if( ctrl.approach == PRODUCT_SVD )
    {
        LogicError("Product SVD not yet supported for bidiagonal matrices");
    }

    // Undo any bidiagonal deflations and flips now that U and V should have had
    // any unnecessary columns dropped.
    //
    // Recall that the application from the left involves the transpose
    // formulation as the application from the right. In particular, 
    //
    //     |  c,  s | | x |
    //     | -s,  c | | y |
    //
    // vs.
    //
    //     | x, y | | c, -s |
    //              | s,  c |.
    //
    // Thus, we need to negate s before the applications from the 
    // opposite side to cheaply effect the transpose.
    if( ctrl.wantU && !ctrl.accumulateU )
    {
        if( newUplo == LOWER )
        {
            // Undo the flip from lower to upper bidiagonal.
            sFlipList *= Real(-1);
            ApplyGivensSequence
            ( LEFT, VARIABLE_GIVENS_SEQUENCE, BACKWARD, 
              cFlipList, sFlipList, U );
        }
        if( uplo == LOWER && !square ) 
        {
            // TODO(poulson): Handle this after adding the original deflation
            LogicError("This case is not yet handled");
        }
    }
    if( ctrl.wantV && !ctrl.accumulateV )
    {
        if( uplo == UPPER && !square )
        {
            sDeflateList *= Real(-1);
            ApplyGivensSequence
            ( LEFT, VARIABLE_GIVENS_SEQUENCE, BACKWARD,
              cDeflateList, sDeflateList, V );
        }
    }

    return info;
}

template<typename F,typename=DisableIf<IsReal<F>>,typename=void>
BidiagSVDInfo
Helper
( UpperOrLower uplo,
  const Matrix<Base<F>>& mainDiagOrig,
  const Matrix<Base<F>>& offDiagOrig,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V,
  const BidiagSVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    typedef Base<F> Real;

    auto ctrlMod = ctrl;
    ctrlMod.accumulateU = false;
    ctrlMod.accumulateV = false;

    Matrix<Real> UReal, VReal;
    auto info =
      Helper( uplo, mainDiagOrig, offDiagOrig, UReal, s, VReal, ctrlMod );

    if( ctrl.wantU )
    {
        if( ctrl.accumulateU )
        {
            // TODO(poulson): Avoid performing in full complex and using so much
            // extra memory
            Matrix<F> UCpx;
            Copy( UReal, UCpx );
            auto UCopy( U );
            Gemm( NORMAL, NORMAL, F(1), UCopy, UCpx, U );
        }
        else
        {
            Copy( UReal, U );
        }
    }

    if( ctrl.wantV )
    {
        if( ctrl.accumulateV )
        {
            // TODO(poulson): Avoid performing in full complex and using so much
            // extra memory
            Matrix<F> VCpx;
            Copy( VReal, VCpx );
            auto VCopy( V );
            Gemm( NORMAL, NORMAL, F(1), VCopy, VCpx, V );
        }
        else
        {
            Copy( VReal, V );
        }
    }

    return info;
}

} // bidiag_svd

template<typename F>
BidiagSVDInfo
BidiagSVD
( UpperOrLower uplo,
  const Matrix<Base<F>>& mainDiagOrig,
  const Matrix<Base<F>>& offDiagOrig,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V,
  const BidiagSVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    return bidiag_svd::Helper( uplo, mainDiagOrig, offDiagOrig, U, s, V, ctrl );
}

#define PROTO(F) \
  template BidiagSVDInfo BidiagSVD \
  ( UpperOrLower uplo, \
    const Matrix<Base<F>>& mainDiag, \
    const Matrix<Base<F>>& offDiag, \
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
    const Matrix<Real>& mainDiag, \
    const Matrix<Real>& offDiag, \
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
