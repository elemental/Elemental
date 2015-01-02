/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CHOLESKY_UVAR3PIVOTED_HPP
#define EL_CHOLESKY_UVAR3PIVOTED_HPP

#include "./LVar3Pivoted.hpp"

namespace El {
namespace cholesky {

template<typename F>
inline void
UUnblockedPivoted( Matrix<F>& A, Matrix<Int>& p )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::UUnblockedPivoted");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    const Int n = A.Height();

    // Initialize the permutation to the identity
    p.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
        p.Set( i, 0, i );
     
    for( Int k=0; k<n; ++k )
    {
        const Range<Int> ind1( k,   k+1 ),
                         ind2( k+1, n   ),
                         indB( k,   n   ),
                         indR( k,   n   );

        auto a12 = A( ind1, ind2 );
        auto A22 = A( ind2, ind2 );
        auto ABR = A( indB, indR );

        // Determine the pivot
        const LDLPivot pivot = pivot::Full( ABR );

        // Apply the pivot
        const Int from = k + pivot.from[0];
        HermitianSwap( UPPER, A, k, from );
        RowSwap( p, k, from );

        // a12 := a12 / sqrt(alpha11)
        const Base<F> delta11 = Sqrt(ABR.GetRealPart(0,0));
        const Base<F> delta11Inv = Base<F>(1)/delta11;
        ABR.Set(0,0,delta11);
        Scale( delta11Inv, a12 );

        // A22 := A22 - a12' a12
        // NOTE: This is silly, but currently necessary
        Conjugate( a12 );
        Her( UPPER, -F(1), a12, A22 );
        Conjugate( a12 );
    }
}

template<typename F>
inline void
UUnblockedPivoted( AbstractDistMatrix<F>& APre, AbstractDistMatrix<Int>& p )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::UUnblockedPivoted");
        if( APre.Height() != APre.Width() )
            LogicError("A must be square");
        AssertSameGrids( APre, p );
    )

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    // Initialize the permutation to the identity
    const Int n = A.Height();
    p.Resize( n, 1 );
    if( p.IsLocalCol(0) )
        for( Int iLoc=0; iLoc<p.LocalHeight(); ++iLoc )
            p.SetLocal( iLoc, 0, p.GlobalRow(iLoc) );

    for( Int k=0; k<n; ++k )
    {
        const Range<Int> ind1( k,   k+1 ),
                         ind2( k+1, n   ),
                         indB( k,   n   ),
                         indR( k,   n   );

        auto a12 = A( ind1, ind2 );
        auto A22 = A( ind2, ind2 );
        auto ABR = A( indB, indR );

        // Determine the pivot 
        const LDLPivot pivot = pivot::Full( ABR );

        // Apply the pivot
        const Int from = k + pivot.from[0];
        HermitianSwap( UPPER, A, k, from );
        RowSwap( p, k, from );

        // a12 := a12 / sqrt(alpha11)
        const Base<F> delta11 = Sqrt(ABR.GetRealPart(0,0));
        const Base<F> delta11Inv = Base<F>(1)/delta11;
        ABR.Set(0,0,delta11);
        Scale( delta11Inv, a12 );

        // A22 := A22 - a12' a12
        // NOTE: This is silly, but currently necessary
        Conjugate( a12 );
        Her( UPPER, -F(1), a12, A22 );
        Conjugate( a12 );
    }
}

// We must use a lazy algorithm so that the symmetric pivoting does not move
// data from a fully-updated to partially-updated region (and vice-versa)
template<typename F>
inline void
UPanelPivoted
( Matrix<F>& AFull, Matrix<Int>& p, 
  Matrix<F>& X, Matrix<F>& Y, Int bsize, Int off )
{
    DEBUG_ONLY(CallStackEntry cse("cholesky::UPanelPivoted"))
    const Int nFull = AFull.Height();
    auto A = AFull( IR(off,nFull), IR(off,nFull) );
    const Int n = A.Height();
    DEBUG_ONLY(
        if( A.Width() != n )
            LogicError("A must be square");
        if( p.Height() != n || p.Width() != 1 )
            LogicError("permutation vector is the wrong size");
    )
    Zeros( X, n, bsize );
    Zeros( Y, n, bsize );

    for( Int k=0; k<bsize; ++k )
    {
        const Range<Int> ind0( 0,   k   ),
                         ind1( k,   k+1 ),
                         ind2( k+1, n   ),
                         indB( k,   n   ),
                         indR( k,   n   );

        auto a12 = A( ind1, ind2 );
        auto a1R = A( ind1, indR );
        auto ABR = A( indB, indR );

        auto x10 = X( ind1, ind0 );
        auto x21 = X( ind2, ind1 );
        auto XB0 = X( indB, ind0 );

        auto y21 = Y( ind2, ind1 );
        auto YB0 = Y( indB, ind0 );

        // Determine the pivot 
        const auto pivot = pivot::PanelFull( ABR, XB0, YB0 );
        const Int from = k + pivot.from[0];

        // Apply the pivot
        HermitianSwap( UPPER, AFull, k+off, from+off );
        RowSwap( p, k, from );
        RowSwap( XB0, 0, pivot.from[0] );
        RowSwap( YB0, 0, pivot.from[0] );

        // Update A(k,k:end) -= X(k,0:k-1) Y(k:end,0:k-1)^T
        Gemv( NORMAL, F(-1), YB0, x10, F(1), a1R );
        a1R.MakeReal(0,0);

        // a12 := a12 / sqrt(alpha11)
        const Base<F> delta11 = Sqrt(A.GetRealPart(k,k));
        const Base<F> delta11Inv = Base<F>(1)/delta11;
        A.Set(k,k,delta11);
        Scale( delta11Inv, a12 );

        // Store x21 := a12' and y21 := a12^T
        Adjoint( a12, x21 );
        Transpose( a12, y21 );
    }
}

template<typename F>
inline void
UPanelPivoted
( DistMatrix<F>& AFull, AbstractDistMatrix<Int>& p, 
  DistMatrix<F,MC,STAR>& X, DistMatrix<F,MR,STAR>& Y, Int bsize, Int off )
{
    DEBUG_ONLY(CallStackEntry cse("cholesky::UPanelPivoted"))
    const Int nFull = AFull.Height();
    auto A = AFull( IR(off,nFull), IR(off,nFull) );
    const Int n = A.Height();
    DEBUG_ONLY(
        if( A.Width() != n )
            LogicError("A must be square");
        if( p.Height() != n || p.Width() != 1 )
            LogicError("pivot vector is the wrong size");
    )
    X.AlignWith( A );
    Y.AlignWith( A );
    Zeros( X, n, bsize );
    Zeros( Y, n, bsize );

    for( Int k=0; k<bsize; ++k )
    {
        const Range<Int> ind0( 0,   k   ),
                         ind1( k,   k+1 ),
                         ind2( k+1, n   ),
                         indB( k,   n   ),
                         indR( k,   n   );

        auto a12 = A( ind1, ind2 );
        auto a1R = A( ind1, indR );
        auto ABR = A( indB, indR );

        auto x10 = X( ind1, ind0 );
        auto x21 = X( ind2, ind1 );
        auto XB0 = X( indB, ind0 );

        auto y21 = Y( ind2, ind1 );
        auto YB0 = Y( indB, ind0 );

        // Determine the pivot
        const auto pivot = pivot::PanelFull( ABR, XB0, YB0 );
        const Int from = k + pivot.from[0];

        // Apply the pivot
        HermitianSwap( UPPER, AFull, k+off, from+off );
        RowSwap( p, k, from );
        RowSwap( XB0, 0, pivot.from[0] );
        RowSwap( YB0, 0, pivot.from[0] );

        // Update A(k,k:end) -= X(k,0:k-1) Y(k:n-off-1,0:k-1)^T
        if( a1R.ColAlign() == a1R.ColRank() )
            LocalGemv( NORMAL, F(-1), YB0, x10, F(1), a1R );
        a1R.MakeReal(0,0);

        // a12 := a12 / sqrt(alpha11)
        const Base<F> delta11 = Sqrt(A.GetRealPart(k,k));
        const Base<F> delta11Inv = Base<F>(1)/delta11;
        A.Set(k,k,delta11);
        Scale( delta11Inv, a12 );

        // Store x21 := a12' and y21 := a12^T
        Adjoint( a12, x21 );
        Transpose( a12, y21 );
    }
}

template<typename F>
inline void
UVar3( Matrix<F>& A, Matrix<Int>& p )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::UVar3");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    const Int n = A.Height();

    // Initialize the permutation to the identity
    p.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
        p.Set( i, 0, i );

    Matrix<F> XB1, YB1;
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> indB( k, n ),
                         indR( k, n );
        auto pB = p( indB, IR(0,1) );
        UPanelPivoted( A, pB, XB1, YB1, nb, k );

        // Update the bottom-right panel
        const Range<Int> ind2( k+nb, n ),
                         ind1Pan( 0,  nb  ),
                         ind2Pan( nb, n-k );
        auto A22 = A( ind2, ind2 );
        auto X21 = XB1( ind2Pan, ind1Pan );
        auto Y21 = YB1( ind2Pan, ind1Pan );
        Trrk( UPPER, NORMAL, TRANSPOSE, F(-1), X21, Y21, F(1), A22 );
    }
}

template<typename F>
inline void
UVar3( AbstractDistMatrix<F>& APre, AbstractDistMatrix<Int>& pPre )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::UVar3");
        if( APre.Height() != APre.Width() )
            LogicError("A must be square");
    )

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre ); auto& A = *APtr;
    auto pPtr = WriteProxy<Int,VC,STAR>( &pPre ); auto& p = *pPtr;

    // Initialize the permutation to the identity
    const Int n = A.Height();
    p.Resize( n, 1 );
    for( Int iLoc=0; iLoc<p.LocalHeight(); ++iLoc )
        p.SetLocal( iLoc, 0, p.GlobalRow(iLoc) );

    const Grid& g = A.Grid();
    DistMatrix<F,MC,STAR> XB1(g);
    DistMatrix<F,MR,STAR> YB1(g);
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> indB( k, n ),
                         indR( k, n );
        auto pB = p( indB, IR(0,1) );
        UPanelPivoted( A, pB, XB1, YB1, nb, k );

        // Update the bottom-right panel
        const Range<Int> ind2( k+nb, n ),
                         ind1Pan( 0,  nb  ),
                         ind2Pan( nb, n-k );
        auto A22 = A( ind2, ind2 );
        auto X21 = XB1( ind2Pan, ind1Pan );
        auto Y21 = YB1( ind2Pan, ind1Pan );
        LocalTrrk( UPPER, TRANSPOSE, F(-1), X21, Y21, F(1), A22 );
    }
}

} // namespace cholesky
} // namespace El

#endif // ifndef EL_CHOLESKY_UVAR3PIVOTED_HPP
