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

template<typename Real>
void PrepareBidiagonal
( UpperOrLower uplo,
  Int m,
  Int n,
  Matrix<Real>& mainDiag,
  Matrix<Real>& offDiag,
  Matrix<Real>& U,
  Matrix<Real>& V,
  Matrix<Real>& cDeflateList,
  Matrix<Real>& sDeflateList,
  Matrix<Real>& cFlipList,
  Matrix<Real>& sFlipList,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE
    if( m == 0 || n == 0 )
        return;

    const bool square = ( m == n );
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
}

template<typename Real>
void PrepareBidiagonal
( UpperOrLower uplo,
  Int m,
  Int n,
  Matrix<Real>& mainDiag,
  Matrix<Real>& offDiag,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE
    Matrix<Real> U, V, cDeflateList, sDeflateList, cFlipList, sFlipList;
    PrepareBidiagonal
    ( uplo, m, n, mainDiag, offDiag, U, V,
      cDeflateList, sDeflateList, cFlipList, sFlipList, ctrl );
}

template<typename Real>
BidiagSVDInfo
Helper
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
    if( minDim == 0 )
    {
        s.Resize( 0, 1 );
        return info;
    }

    PrepareBidiagonal( uplo, m, n, mainDiag, offDiag, ctrl );

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
        s.Resize( rank, 1 );
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

template<typename Real>
BidiagSVDInfo
Helper
( UpperOrLower uplo,
  DistMatrix<Real,STAR,STAR>& mainDiag,
  DistMatrix<Real,STAR,STAR>& offDiag,
  DistMatrix<Real,STAR,STAR>& s,
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
    if( minDim == 0 )
    {
        s.Resize( 0, 1 );
        return info;
    }

    PrepareBidiagonal
    ( uplo, m, n, mainDiag.Matrix(), offDiag.Matrix(), ctrl );

    if( square )
    {
        s = mainDiag;
        info.qrInfo = bidiag_svd::QRAlg( s.Matrix(), offDiag.Matrix(), ctrl );
        Sort( s, DESCENDING );
    }
    else if( uplo == LOWER )
    {
        // We were non-square and lower bidiagonal. 
        auto offDiag0 = offDiag( IR(0,n-1), ALL );
        s = mainDiag;
        info.qrInfo = bidiag_svd::QRAlg( s.Matrix(), offDiag0.Matrix(), ctrl );
        Sort( s, DESCENDING );
    }
    else
    {
        // We were non-square and upper bidiagonal. 
        auto offDiag0 = offDiag( IR(0,m-1), ALL );
        s = mainDiag;
        info.qrInfo = bidiag_svd::QRAlg( s.Matrix(), offDiag0.Matrix(), ctrl );
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
        auto& sLoc = s.Matrix();
        for( Int i=0; i<minDim; ++i )
        {
            if( sLoc(i) <= Real(0) )
            {
                rank = i; 
                break;
            } 
        }
        s.Resize( rank, 1 );
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

template<typename Real>
BidiagSVDInfo
Helper
( UpperOrLower uplo,
  const Matrix<Real>& mainDiagOrig,
  const Matrix<Real>& offDiagOrig,
        Matrix<Real>& s,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE
    auto mainDiag( mainDiagOrig );
    auto offDiag( offDiagOrig );
    return Helper( uplo, mainDiag, offDiag, s, ctrl );
}

template<typename Real>
BidiagSVDInfo
Helper
( UpperOrLower uplo,
  const AbstractDistMatrix<Real>& mainDiagPre,
  const AbstractDistMatrix<Real>& offDiagPre,
        AbstractDistMatrix<Real>& sPre,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE
    DistMatrix<Real,STAR,STAR> mainDiag(mainDiagPre), offDiag(offDiagPre);
    DistMatrixWriteProxy<Real,Real,STAR,STAR> sProx( sPre );
    auto& s = sProx.Get();
    return Helper( uplo, mainDiag, offDiag, s, ctrl );
}

template<typename Real>
BidiagSVDInfo
Helper
( UpperOrLower uplo,
  const Matrix<Complex<Real>>& mainDiag,
  const Matrix<Complex<Real>>& offDiag,
        Matrix<Real>& s,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE
    const Int mainDiagHeight = mainDiag.Height();
    const Int offDiagHeight = offDiag.Height();
    // We can implicitly rotate by the complex conjugates
    Matrix<Real> mainDiagReal(mainDiagHeight,1),
                 offDiagReal(offDiagHeight,1);
    for( Int i=0; i<mainDiagHeight; ++i )
        mainDiagReal(i) = Abs(mainDiag(i));
    for( Int i=0; i<offDiagHeight; ++i )
        offDiagReal(i) = Abs(offDiag(i));
    return Helper( uplo, mainDiagReal, offDiagReal, s, ctrl );
}

template<typename Real>
BidiagSVDInfo
Helper
( UpperOrLower uplo,
  const AbstractDistMatrix<Complex<Real>>& mainDiagPre,
  const AbstractDistMatrix<Complex<Real>>& offDiagPre,
        AbstractDistMatrix<Real>& s,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE
    DistMatrix<Complex<Real>,VC,STAR>
      mainDiag(mainDiagPre), offDiag(offDiagPre);
    const Grid& g = mainDiag.Grid();

    // We can implicitly rotate by the complex conjugates
    auto& mainDiagLoc = mainDiag.Matrix();
    auto& offDiagLoc = offDiag.Matrix();
    const Int mainDiagLocHeight = mainDiagLoc.Height();
    for( Int iLoc=0; iLoc<mainDiagLocHeight; ++iLoc )
        mainDiagLoc(iLoc) = Abs(mainDiagLoc(iLoc));
    const Int offDiagLocHeight = offDiagLoc.Height();
    for( Int iLoc=0; iLoc<offDiagLocHeight; ++iLoc ) 
        offDiagLoc(iLoc) = Abs(offDiagLoc(iLoc));

    DistMatrix<Real,STAR,STAR> mainDiagReal(g), offDiagReal(g);
    return Helper( uplo, mainDiagReal, offDiagReal, s, ctrl );
}

} // namespace bidiag_svd

template<typename F>
BidiagSVDInfo
BidiagSVD
( UpperOrLower uplo,
  const Matrix<F>& mainDiagOrig,
  const Matrix<F>& offDiagOrig,
        Matrix<Base<F>>& s,
  const BidiagSVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    return bidiag_svd::Helper( uplo, mainDiagOrig, offDiagOrig, s, ctrl );
}

template<typename F>
BidiagSVDInfo
BidiagSVD
( UpperOrLower uplo,
  const AbstractDistMatrix<F>& mainDiagOrig,
  const AbstractDistMatrix<F>& offDiagOrig,
        AbstractDistMatrix<Base<F>>& s,
  const BidiagSVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    return bidiag_svd::Helper( uplo, mainDiagOrig, offDiagOrig, s, ctrl );
}

namespace bidiag_svd {

template<typename Real>
BidiagSVDInfo
Helper
( UpperOrLower uplo,
  Matrix<Real>& mainDiag,
  Matrix<Real>& offDiag,
  Matrix<Real>& U,
  Matrix<Real>& s,
  Matrix<Real>& V,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE

    const Int m = ( uplo==UPPER ? mainDiag.Height() : offDiag.Height()+1 );
    const Int n = ( uplo==UPPER ? offDiag.Height()+1 : mainDiag.Height() );
    const Int minDim = Min(m,n);
    const bool square = ( m == n );
    if( mainDiag.Height() != offDiag.Height() &&
        mainDiag.Height() != offDiag.Height()+1 )
        LogicError("Invalid main and superdiagonal lengths");

    BidiagSVDInfo info;
    if( minDim == 0 )
    {
        // It would appear to be silly to form a "full" SVD of a 0 x n matrix
        // with an identity V matrix (or vice-versa for an m x 0 matrix with an
        // m x m identity U).
        U.Resize( 0, 0 );
        s.Resize( 0, 1 );
        V.Resize( 0, 0 );
        return info;
    }
    if( !ctrl.wantU && !ctrl.wantV )
    {
        U.Resize( 0, 0 );
        V.Resize( 0, 0 );
        return Helper( uplo, mainDiag, offDiag, s, ctrl );
    }

    // We currently only support QR and D&C
    bool useDC = !ctrl.useQR;
    auto ctrlMod( ctrl );
    if( useDC && m <= ctrl.dcCtrl.cutoff )
    {
        // The D&C algorithm will actually just be the QR algorithm, so take
        // advantage of the fact that the QR algorithm supports direct
        // accumulation.
        ctrlMod.useQR = true;
        useDC = false;
    }

    // TODO(poulson): Decide if accumulating into U/V is ever worthwhile for 
    // performance reasons given that there should be O(n^2 k) to subsequently
    // compose the original and output matrices, but O(n^3) with a significant 
    // coefficient to directly accumulate.

    Matrix<Real> cDeflateList, sDeflateList;
    Matrix<Real> cFlipList, sFlipList;
    PrepareBidiagonal
    ( uplo, m, n, mainDiag, offDiag, U, V,
      cDeflateList, sDeflateList, cFlipList, sFlipList, ctrlMod );

    // Since D&C does not support direct accumulation, we must make a copy.
    Matrix<Real> UCopy, VCopy;
    if( useDC )
    {
        if( ctrlMod.wantU && ctrlMod.accumulateU )
            UCopy = U;
        if( ctrlMod.wantV && ctrlMod.accumulateV )
            VCopy = V;
    }

    if( square )
    {
        if( ctrlMod.useQR )
        {
            s = mainDiag;
            info.qrInfo = bidiag_svd::QRAlg( s, offDiag, U, V, ctrlMod );
        }
        else
        {
            info.dcInfo =
              bidiag_svd::DivideAndConquer
              ( mainDiag, offDiag, U, s, V, ctrlMod );
        }

        auto sortPairs = TaggedSort( s, DESCENDING );
        ApplyTaggedSortToEachColumn( sortPairs, s );
        if( ctrlMod.wantU )
            ApplyTaggedSortToEachRow( sortPairs, U );
        if( ctrlMod.wantV )
            ApplyTaggedSortToEachRow( sortPairs, V );
    }
    else if( uplo == LOWER )
    {
        // We were non-square and lower bidiagonal. The last column of U has
        // been (at least implicitly) deflated.
        auto offDiag0 = offDiag( IR(0,n-1), ALL );
        Matrix<Real> U0;
        if( ctrlMod.wantU )
        {
            if( !ctrlMod.accumulateU )
                Identity( U, m, m ); 
            View( U0, U, ALL, IR(0,n) );
        }

        if( ctrlMod.useQR )
        {
            s = mainDiag;
            info.qrInfo = bidiag_svd::QRAlg( s, offDiag0, U0, V, ctrlMod );
        }
        else
        {
            info.dcInfo =
              bidiag_svd::DivideAndConquer
              ( mainDiag, offDiag0, U0, s, V, ctrlMod );
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
        if( ctrlMod.wantV )
        {
            if( !ctrlMod.accumulateV )
                Identity( V, n, n );
            View( V0, V, ALL, IR(0,m) );
        }

        if( ctrlMod.useQR )
        {
            s = mainDiag;
            info.qrInfo = bidiag_svd::QRAlg( s, offDiag0, U, V0, ctrlMod );
        }
        else
        {
            info.dcInfo =
              bidiag_svd::DivideAndConquer
              ( mainDiag, offDiag0, U, s, V0, ctrl );
        }

        auto sortPairs = TaggedSort( s, DESCENDING );
        ApplyTaggedSortToEachColumn( sortPairs, s );
        if( ctrlMod.wantU )
            ApplyTaggedSortToEachRow( sortPairs, U );
        if( ctrlMod.wantV )
            ApplyTaggedSortToEachRow( sortPairs, V0 );
    }

    if( ctrlMod.approach == THIN_SVD )
    {
        if( ctrlMod.wantU )
            U.Resize( m, minDim );
        if( ctrlMod.wantV )
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
                Output("rank=",rank);
                break;
            } 
        }
        s.Resize( rank, 1 );
        if( ctrlMod.wantU )
            U.Resize( m, rank );
        if( ctrlMod.wantV )
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
    if( ctrlMod.wantU && !ctrlMod.accumulateU )
    {
        if( uplo == UPPER && !square )
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
    if( ctrlMod.wantV && !ctrlMod.accumulateV )
    {
        if( uplo == UPPER && !square )
        {
            sDeflateList *= Real(-1);
            ApplyGivensSequence
            ( LEFT, VARIABLE_GIVENS_SEQUENCE, BACKWARD,
              cDeflateList, sDeflateList, V );
        }
    }

    // As noted above, D&C does not support direct accumulation, so directly
    // form the necessary product between the input matrices and the computed
    // singular vectors.
    if( useDC )
    {
        if( ctrlMod.wantU && ctrlMod.accumulateU )
        {
            // U := UCopy U
            Matrix<Real> temp;
            Gemm( NORMAL, NORMAL, Real(1), UCopy, U, temp ); 
            U = temp;
        }
        if( ctrlMod.wantV && ctrlMod.accumulateV )
        {
            // V := VCopy V
            Matrix<Real> temp;
            Gemm( NORMAL, NORMAL, Real(1), VCopy, V, temp );
            V = temp;
        }
    }

    return info;
}

template<typename Real>
BidiagSVDInfo
Helper
( UpperOrLower uplo,
  DistMatrix<Real,STAR,STAR>& mainDiag,
  DistMatrix<Real,STAR,STAR>& offDiag,
  DistMatrix<Real,VC,STAR>& U,
  DistMatrix<Real,STAR,STAR>& s,
  DistMatrix<Real,VC,STAR>& V,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE

    const Grid& g = mainDiag.Grid();
    const Int m = ( uplo==UPPER ? mainDiag.Height() : offDiag.Height()+1 );
    const Int n = ( uplo==UPPER ? offDiag.Height()+1 : mainDiag.Height() );
    const Int minDim = Min(m,n);
    const bool square = ( m == n );
    if( mainDiag.Height() != offDiag.Height() &&
        mainDiag.Height() != offDiag.Height()+1 )
        LogicError("Invalid main and superdiagonal lengths");

    BidiagSVDInfo info;
    if( minDim == 0 )
    {
        // It would appear to be silly to form a "full" SVD of a 0 x n matrix
        // with an identity V matrix (or vice-versa for an m x 0 matrix with an
        // m x m identity U).
        U.Resize( 0, 0 );
        s.Resize( 0, 1 );
        V.Resize( 0, 0 );
        return info;
    }
    if( !ctrl.wantU && !ctrl.wantV )
    {
        U.Resize( 0, 0 );
        V.Resize( 0, 0 );
        return Helper( uplo, mainDiag, offDiag, s, ctrl );
    }

    // We currently only support QR and D&C
    bool useDC = !ctrl.useQR;
    auto ctrlMod( ctrl );
    if( useDC && m <= ctrl.dcCtrl.cutoff )
    {
        // The D&C algorithm will actually just be the QR algorithm, so take
        // advantage of the fact that the QR algorithm supports direct
        // accumulation.
        ctrlMod.useQR = true;
        useDC = false;
    }

    // TODO(poulson): Decide if accumulating into U/V is ever worthwhile for 
    // performance reasons given that there should be O(n^2 k) to subsequently
    // compose the original and output matrices, but O(n^3) with a significant 
    // coefficient to directly accumulate.

    Matrix<Real> cDeflateList, sDeflateList;
    Matrix<Real> cFlipList, sFlipList;
    PrepareBidiagonal
    ( uplo, m, n, mainDiag.Matrix(), offDiag.Matrix(), U.Matrix(), V.Matrix(),
      cDeflateList, sDeflateList, cFlipList, sFlipList, ctrlMod );

    // Since D&C does not support direct accumulation, we must make a copy.
    DistMatrix<Real> UCopy(g), VCopy(g);
    if( useDC )
    {
        if( ctrlMod.wantU && ctrlMod.accumulateU )
            UCopy = U;
        if( ctrlMod.wantV && ctrlMod.accumulateV )
            VCopy = V;
    }

    if( square )
    {
        if( ctrlMod.useQR )
        {
            // TODO(poulson): Move this logic into a distributed version of 
            // bidiag_svd::QRAlg since a distributed version of the D&C is
            // needed anyway.
            auto ctrlModQR( ctrlMod );
            ctrlModQR.accumulateU = true;
            ctrlModQR.accumulateV = true;
            if( ctrlMod.wantU && !ctrlMod.accumulateU )
                Identity( U, m, m );
            if( ctrlMod.wantV && !ctrlMod.accumulateV )
                Identity( V, n, n );

            s = mainDiag;
            info.qrInfo =
              bidiag_svd::QRAlg
              ( s.Matrix(), offDiag.Matrix(), U.Matrix(), V.Matrix(),
                ctrlModQR );
        }
        else
        {
            DistMatrix<Real> U_MC_MR(g), V_MC_MR(g);
            if( ctrl.wantU )
                U_MC_MR = U;
            if( ctrl.wantV )
                V_MC_MR = V;
            const bool topLevel = true;
            info.dcInfo =
              bidiag_svd::DivideAndConquer
              ( mainDiag.Matrix(), offDiag.Matrix(),
                U_MC_MR, s, V_MC_MR, ctrlMod, topLevel );
            if( ctrl.wantU )
                U = U_MC_MR;
            if( ctrl.wantV )
                V = V_MC_MR;
        }

        auto sortPairs = TaggedSort( s, DESCENDING );
        ApplyTaggedSortToEachColumn( sortPairs, s );
        if( ctrlMod.wantU )
            ApplyTaggedSortToEachRow( sortPairs, U );
        if( ctrlMod.wantV )
            ApplyTaggedSortToEachRow( sortPairs, V );
    }
    else if( uplo == LOWER )
    {
        // We were non-square and lower bidiagonal. The last column of U has
        // been (at least implicitly) deflated.
        auto offDiag0 = offDiag( IR(0,n-1), ALL );
        DistMatrix<Real,VC,STAR> U0(g);
        if( ctrlMod.wantU )
        {
            if( !ctrlMod.accumulateU )
                Identity( U, m, m ); 
            View( U0, U, ALL, IR(0,n) );
        }

        if( ctrlMod.useQR )
        {
            auto ctrlModQR( ctrlMod );
            ctrlModQR.accumulateU = true;
            ctrlModQR.accumulateV = true;
            if( ctrlMod.wantV && !ctrlMod.accumulateV )
                Identity( V, n, n );

            s = mainDiag;
            info.qrInfo =
              bidiag_svd::QRAlg
              ( s.Matrix(), offDiag0.Matrix(), U0.Matrix(), V.Matrix(),
                ctrlModQR );
        }
        else
        {
            DistMatrix<Real> U0_MC_MR(g), V_MC_MR(g);
            if( ctrl.wantU )
                U0_MC_MR = U0;
            if( ctrl.wantV )
                V_MC_MR = V;
            const bool topLevel = true;
            info.dcInfo =
              bidiag_svd::DivideAndConquer
              ( mainDiag.Matrix(), offDiag0.Matrix(),
                U0_MC_MR, s, V_MC_MR, ctrlMod, topLevel );
            if( ctrl.wantU )
                U0 = U0_MC_MR;
            if( ctrl.wantV )
                V = V_MC_MR;
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
        DistMatrix<Real,VC,STAR> V0(g);
        if( ctrlMod.wantV )
        {
            if( !ctrlMod.accumulateV )
                Identity( V, n, n );
            View( V0, V, ALL, IR(0,m) );
        }

        if( ctrlMod.useQR )
        {
            auto ctrlModQR( ctrlMod );
            ctrlModQR.accumulateU = true;
            ctrlModQR.accumulateV = true;
            if( ctrlMod.wantU && !ctrlMod.accumulateU )
                Identity( U, m, m );

            s = mainDiag;
            info.qrInfo =
              bidiag_svd::QRAlg
              ( s.Matrix(), offDiag0.Matrix(), U.Matrix(), V0.Matrix(),
                ctrlMod );
        }
        else
        {
            DistMatrix<Real> U_MC_MR(g), V0_MC_MR(g);
            if( ctrl.wantU )
                U_MC_MR = U;
            if( ctrl.wantV )
                V0_MC_MR = V0;
            const bool topLevel = true;
            info.dcInfo =
              bidiag_svd::DivideAndConquer
              ( mainDiag.Matrix(), offDiag0.Matrix(),
                U_MC_MR, s, V0_MC_MR, ctrlMod, topLevel );
            if( ctrl.wantU )
                U = U_MC_MR;
            if( ctrl.wantV )
                V0 = V0_MC_MR;
        }

        auto sortPairs = TaggedSort( s, DESCENDING );
        ApplyTaggedSortToEachColumn( sortPairs, s );
        if( ctrlMod.wantU )
            ApplyTaggedSortToEachRow( sortPairs, U );
        if( ctrlMod.wantV )
            ApplyTaggedSortToEachRow( sortPairs, V0 );
    }

    if( ctrlMod.approach == THIN_SVD )
    {
        if( ctrlMod.wantU )
            U.Resize( m, minDim );
        if( ctrlMod.wantV )
            V.Resize( n, minDim );
    }
    else if( ctrlMod.approach == COMPACT_SVD )
    {
        // Determine the rank
        Int rank = minDim;
        auto& sLoc = s.Matrix();
        for( Int i=0; i<minDim; ++i )
        {
            if( sLoc(i) <= Real(0) )
            {
                rank = i; 
                break;
            } 
        }
        s.Resize( rank, 1 );
        if( ctrlMod.wantU )
            U.Resize( m, rank );
        if( ctrlMod.wantV )
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
    if( ctrlMod.wantU && !ctrlMod.accumulateU )
    {
        if( uplo == UPPER && !square )
        {
            // Undo the flip from lower to upper bidiagonal.
            sFlipList *= Real(-1);
            DistMatrix<Real,STAR,VR> U_STAR_VR( U );
            ApplyGivensSequence
            ( LEFT, VARIABLE_GIVENS_SEQUENCE, BACKWARD, 
              cFlipList, sFlipList, U_STAR_VR.Matrix() );
            U = U_STAR_VR;
        }
        if( uplo == LOWER && !square ) 
        {
            // TODO(poulson): Handle this after adding the original deflation
            LogicError("This case is not yet handled");
        }
    }
    if( ctrlMod.wantV && !ctrlMod.accumulateV )
    {
        if( uplo == UPPER && !square )
        {
            sDeflateList *= Real(-1);
            DistMatrix<Real,STAR,VR> V_STAR_VR( V );
            ApplyGivensSequence
            ( LEFT, VARIABLE_GIVENS_SEQUENCE, BACKWARD,
              cDeflateList, sDeflateList, V_STAR_VR.Matrix() );
            V = V_STAR_VR;
        }
    }

    // As noted above, D&C does not support direct accumulation, so directly
    // form the necessary product between the input matrices and the computed
    // singular vectors.
    if( useDC )
    {
        if( ctrlMod.wantU && ctrlMod.accumulateU )
        {
            // U := UCopy U
            DistMatrix<Real> temp(g);
            Gemm( NORMAL, NORMAL, Real(1), UCopy, U, temp ); 
            U = temp;
        }
        if( ctrlMod.wantV && ctrlMod.accumulateV )
        {
            // V := VCopy V
            DistMatrix<Real> temp(g);
            Gemm( NORMAL, NORMAL, Real(1), VCopy, V, temp );
            V = temp;
        }
    }

    return info;
}

template<typename Real>
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
    auto mainDiag( mainDiagOrig );
    auto offDiag( offDiagOrig );
    return Helper( uplo, mainDiag, offDiag, U, s, V, ctrl );
}

template<typename Real>
BidiagSVDInfo
Helper
( UpperOrLower uplo,
  const AbstractDistMatrix<Real>& mainDiagOrig,
  const AbstractDistMatrix<Real>& offDiagOrig,
        AbstractDistMatrix<Real>& UPre,
        AbstractDistMatrix<Real>& sPre,
        AbstractDistMatrix<Real>& VPre,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE
    const Grid& g = mainDiagOrig.Grid();
    DistMatrix<Real,STAR,STAR> mainDiag( mainDiagOrig );
    DistMatrix<Real,STAR,STAR> offDiag( offDiagOrig );

    DistMatrixWriteProxy<Real,Real,STAR,STAR> sProx( sPre );
    DistMatrixReadWriteProxy<Real,Real,VC,STAR> UProx( UPre ), VProx( VPre );
    auto& s = sProx.Get();
    auto& U = UProx.Get();
    auto& V = VProx.Get();

    return Helper( uplo, mainDiag, offDiag, U, s, V, ctrl );
}

template<typename Real>
BidiagSVDInfo
Helper
( UpperOrLower uplo,
  Matrix<Real>& mainDiag,
  Matrix<Real>& offDiag,
  Matrix<Complex<Real>>& U,
  Matrix<Real>& s,
  Matrix<Complex<Real>>& V,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE
    typedef Complex<Real> F;

    auto ctrlMod = ctrl;
    ctrlMod.accumulateU = false;
    ctrlMod.accumulateV = false;

    Matrix<Real> UReal, VReal;
    auto info = Helper( uplo, mainDiag, offDiag, UReal, s, VReal, ctrlMod );

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

template<typename Real>
BidiagSVDInfo
Helper
( UpperOrLower uplo,
  DistMatrix<Real,STAR,STAR>& mainDiag,
  DistMatrix<Real,STAR,STAR>& offDiag,
  DistMatrix<Complex<Real>,VC,STAR>& U,
  DistMatrix<Real,STAR,STAR>& s,
  DistMatrix<Complex<Real>,VC,STAR>& V,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE
    typedef Complex<Real> F;
    const Grid& g = mainDiag.Grid();

    auto ctrlMod = ctrl;
    ctrlMod.accumulateU = false;
    ctrlMod.accumulateV = false;

    DistMatrix<Real> UReal(g), VReal(g);
    auto info = Helper( uplo, mainDiag, offDiag, UReal, s, VReal, ctrlMod );

    if( ctrl.wantU )
    {
        if( ctrl.accumulateU )
        {
            // TODO(poulson): Avoid performing in full complex and using so much
            // extra memory
            DistMatrix<F> UCpx(g);
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
            DistMatrix<F> VCpx(g);
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

template<typename Real>
BidiagSVDInfo
Helper
( UpperOrLower uplo,
  const Matrix<Real>& mainDiagOrig,
  const Matrix<Real>& offDiagOrig,
        Matrix<Complex<Real>>& U,
        Matrix<Real>& s,
        Matrix<Complex<Real>>& V,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE
    auto mainDiag( mainDiagOrig );
    auto offDiag( offDiagOrig );
    return Helper( uplo, mainDiag, offDiag, U, s, V, ctrl );
}

template<typename Real>
BidiagSVDInfo
Helper
( UpperOrLower uplo,
  const AbstractDistMatrix<Real>& mainDiagOrig,
  const AbstractDistMatrix<Real>& offDiagOrig,
        AbstractDistMatrix<Complex<Real>>& UPre,
        AbstractDistMatrix<Real>& sPre,
        AbstractDistMatrix<Complex<Real>>& VPre,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE
    typedef Complex<Real> F;
    DistMatrix<Real,STAR,STAR> mainDiag( mainDiagOrig ),
      offDiag( offDiagOrig );

    DistMatrixWriteProxy<Real,Real,STAR,STAR> sProx( sPre );
    DistMatrixReadWriteProxy<F,F,VC,STAR> UProx( UPre ), VProx( VPre );
    auto& s = sProx.Get();
    auto& U = UProx.Get();
    auto& V = VProx.Get();

    return Helper( uplo, mainDiag, offDiag, U, s, V, ctrl );
}

template<typename Real>
BidiagSVDInfo
Helper
( UpperOrLower uplo,
  const Matrix<Complex<Real>>& mainDiag,
  const Matrix<Complex<Real>>& offDiag,
        Matrix<Complex<Real>>& U,
        Matrix<Real>& s,
        Matrix<Complex<Real>>& V,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE
    typedef Complex<Real> F;

    const Int mainDiagHeight = mainDiag.Height();
    const Int offDiagHeight = offDiag.Height();
    if( mainDiagHeight != offDiagHeight && mainDiagHeight != offDiagHeight+1 )
        LogicError("Invalid lengths of mainDiag and offDiag");
    const Int m = ( uplo==UPPER ? mainDiagHeight : offDiagHeight+1 );
    const Int n = ( uplo==UPPER ? offDiagHeight+1 : mainDiagHeight );

    Matrix<Complex<Real>> UPhase, VPhase;
    if( ctrl.wantU || ctrl.wantV )
    {
        Ones( UPhase, m, 1 );
        Ones( VPhase, n, 1 );
        F alphaNew, betaNew;
        if( uplo == UPPER )
        {
            for( Int i=0; i<mainDiagHeight; ++i ) 
            {
                alphaNew = mainDiag(i) / VPhase(i);
                UPhase(i) = Phase( alphaNew, false );
                if( offDiagHeight > i )
                {
                    betaNew = offDiag(i) / UPhase(i);
                    VPhase(i+1) = Phase( betaNew, false );
                }
            }
        }
        else
        {
            for( Int i=0; i<mainDiagHeight; ++i ) 
            {
                alphaNew = mainDiag(i) / UPhase(i);
                VPhase(i) = Phase( alphaNew, false );
                if( offDiagHeight > i )
                {
                    betaNew = offDiag(i) / VPhase(i);
                    UPhase(i+1) = Phase( betaNew, false ); 
                }
            }
        }
        Conjugate( VPhase );
    }
    Matrix<Real> mainDiagReal(mainDiagHeight,1), offDiagReal(offDiagHeight,1);
    for( Int i=0; i<mainDiagHeight; ++i )
    {
        mainDiagReal(i) = Abs(mainDiag(i));
        if( offDiagHeight > i )
            offDiagReal(i) = Abs(offDiag(i));
    }
    
    auto info = Helper( uplo, mainDiagReal, offDiagReal, U, s, V, ctrl );

    // Apply the phases as necessary
    if( ctrl.wantU )
        DiagonalScale( LEFT, NORMAL, UPhase, U );
    if( ctrl.wantV )
        DiagonalScale( LEFT, NORMAL, VPhase, V );

    return info;
}

template<typename Real>
BidiagSVDInfo
Helper
( UpperOrLower uplo,
  const AbstractDistMatrix<Complex<Real>>& mainDiagOrig,
  const AbstractDistMatrix<Complex<Real>>& offDiagOrig,
        AbstractDistMatrix<Complex<Real>>& UPre,
        AbstractDistMatrix<Real>& sPre,
        AbstractDistMatrix<Complex<Real>>& VPre,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE
    typedef Complex<Real> F;
    const Grid& g = mainDiagOrig.Grid();

    DistMatrix<F,STAR,STAR> mainDiag( mainDiagOrig ),
      offDiag( offDiagOrig );

    DistMatrixWriteProxy<Real,Real,STAR,STAR> sProx( sPre );
    DistMatrixReadWriteProxy<F,F,VC,STAR> UProx( UPre ), VProx( VPre );
    auto& s = sProx.Get();
    auto& U = UProx.Get();
    auto& V = VProx.Get();

    const Int mainDiagHeight = mainDiag.Height();
    const Int offDiagHeight = offDiag.Height();
    if( mainDiagHeight != offDiagHeight && mainDiagHeight != offDiagHeight+1 )
        LogicError("Invalid lengths of mainDiag and offDiag");
    const Int m = ( uplo==UPPER ? mainDiagHeight : offDiagHeight+1 );
    const Int n = ( uplo==UPPER ? offDiagHeight+1 : mainDiagHeight );

    auto& mainDiagLoc = mainDiag.Matrix();
    auto& offDiagLoc = offDiag.Matrix();
    DistMatrix<F,STAR,STAR> UPhase(g), VPhase(g);
    if( ctrl.wantU || ctrl.wantV )
    {
        Ones( UPhase, m, 1 );
        Ones( VPhase, n, 1 );
        auto& UPhaseLoc = UPhase.Matrix();
        auto& VPhaseLoc = VPhase.Matrix();

        F alphaNew, betaNew;
        if( uplo == UPPER )
        {
            for( Int i=0; i<mainDiagHeight; ++i ) 
            {
                alphaNew = mainDiagLoc(i) / VPhaseLoc(i);
                UPhaseLoc(i) = Phase( alphaNew, false );
                if( offDiagHeight > i )
                {
                    betaNew = offDiagLoc(i) / UPhaseLoc(i);
                    VPhaseLoc(i+1) = Phase( betaNew, false );
                }
            }
        }
        else
        {
            for( Int i=0; i<mainDiagHeight; ++i ) 
            {
                alphaNew = mainDiagLoc(i) / UPhaseLoc(i);
                VPhaseLoc(i) = Phase( alphaNew, false );
                if( offDiagHeight > i )
                {
                    betaNew = offDiagLoc(i) / VPhaseLoc(i);
                    UPhaseLoc(i+1) = Phase( betaNew, false ); 
                }
            }
        }
        Conjugate( VPhase );
    }
    DistMatrix<Real,STAR,STAR> mainDiagReal(mainDiagHeight,1,g),
                               offDiagReal(offDiagHeight,1,g);
    auto& mainDiagRealLoc = mainDiagReal.Matrix();
    auto& offDiagRealLoc = offDiagReal.Matrix();
    for( Int i=0; i<mainDiagHeight; ++i )
    {
        mainDiagRealLoc(i) = Abs(mainDiagLoc(i));
        if( offDiagHeight > i )
            offDiagRealLoc(i) = Abs(offDiagLoc(i));
    }
    
    auto info = Helper( uplo, mainDiagReal, offDiagReal, U, s, V, ctrl );

    // Apply the phases as necessary
    if( ctrl.wantU )
        DiagonalScale( LEFT, NORMAL, UPhase, U );
    if( ctrl.wantV )
        DiagonalScale( LEFT, NORMAL, VPhase, V );

    return info;
}

} // bidiag_svd

template<typename F>
BidiagSVDInfo
BidiagSVD
( UpperOrLower uplo,
  const Matrix<F>& mainDiag,
  const Matrix<F>& offDiag,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V,
  const BidiagSVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    return bidiag_svd::Helper( uplo, mainDiag, offDiag, U, s, V, ctrl );
}

template<typename F>
BidiagSVDInfo
BidiagSVD
( UpperOrLower uplo,
  const AbstractDistMatrix<F>& mainDiag,
  const AbstractDistMatrix<F>& offDiag,
        AbstractDistMatrix<F>& U,
        AbstractDistMatrix<Base<F>>& s,
        AbstractDistMatrix<F>& V,
  const BidiagSVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    return bidiag_svd::Helper( uplo, mainDiag, offDiag, U, s, V, ctrl );
}

template<typename F>
BidiagSVDInfo
BidiagSVD
( UpperOrLower uplo,
  const Matrix<Base<F>>& mainDiag,
  const Matrix<Base<F>>& offDiag,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V,
  const BidiagSVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    return bidiag_svd::Helper( uplo, mainDiag, offDiag, U, s, V, ctrl );
}

template<typename F>
BidiagSVDInfo
BidiagSVD
( UpperOrLower uplo,
  const AbstractDistMatrix<Base<F>>& mainDiag,
  const AbstractDistMatrix<Base<F>>& offDiag,
        AbstractDistMatrix<F>& U,
        AbstractDistMatrix<Base<F>>& s,
        AbstractDistMatrix<F>& V,
  const BidiagSVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    return bidiag_svd::Helper( uplo, mainDiag, offDiag, U, s, V, ctrl );
}

#define PROTO(F) \
  template BidiagSVDInfo BidiagSVD \
  ( UpperOrLower uplo, \
    const Matrix<F>& mainDiag, \
    const Matrix<F>& offDiag, \
          Matrix<Base<F>>& s, \
    const BidiagSVDCtrl<Base<F>>& ctrl ); \
  template BidiagSVDInfo BidiagSVD \
  ( UpperOrLower uplo, \
    const AbstractDistMatrix<F>& mainDiag, \
    const AbstractDistMatrix<F>& offDiag, \
          AbstractDistMatrix<Base<F>>& s, \
    const BidiagSVDCtrl<Base<F>>& ctrl ); \
  template BidiagSVDInfo BidiagSVD \
  ( UpperOrLower uplo, \
    const Matrix<F>& mainDiag, \
    const Matrix<F>& offDiag, \
          Matrix<F>& U, \
          Matrix<Base<F>>& s, \
          Matrix<F>& V, \
    const BidiagSVDCtrl<Base<F>>& ctrl ); \
  template BidiagSVDInfo BidiagSVD \
  ( UpperOrLower uplo, \
    const AbstractDistMatrix<F>& mainDiag, \
    const AbstractDistMatrix<F>& offDiag, \
          AbstractDistMatrix<F>& U, \
          AbstractDistMatrix<Base<F>>& s, \
          AbstractDistMatrix<F>& V, \
    const BidiagSVDCtrl<Base<F>>& ctrl );

#define PROTO_REAL(Real) \
  PROTO(Real) \
  template BidiagSVDInfo BidiagSVD \
  ( UpperOrLower uplo, \
    const Matrix<Real>& mainDiag, \
    const Matrix<Real>& offDiag, \
          Matrix<Complex<Real>>& U, \
          Matrix<Real>& s, \
          Matrix<Complex<Real>>& V, \
    const BidiagSVDCtrl<Real>& ctrl ); \
  template BidiagSVDInfo BidiagSVD \
  ( UpperOrLower uplo, \
    const AbstractDistMatrix<Real>& mainDiag, \
    const AbstractDistMatrix<Real>& offDiag, \
          AbstractDistMatrix<Complex<Real>>& U, \
          AbstractDistMatrix<Real>& s, \
          AbstractDistMatrix<Complex<Real>>& V, \
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
