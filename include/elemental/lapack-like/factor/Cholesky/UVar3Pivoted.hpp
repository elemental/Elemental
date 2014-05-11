/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CHOLESKY_UVAR3PIVOTED_HPP
#define ELEM_CHOLESKY_UVAR3PIVOTED_HPP

#include "./LVar3Pivoted.hpp"

namespace elem {
namespace cholesky {

template<typename F>
inline void
UUnblockedPivoted( Matrix<F>& A, Matrix<Int>& pPerm )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::UUnblockedPivoted");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    const Int n = A.Height();

    // Initialize the permutation to the identity
    pPerm.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
        pPerm.Set( i, 0, i );
     
    for( Int k=0; k<n; ++k )
    {
        // Determine the pivot
        auto ABR = ViewRange( A, k, k, n, n );
        const LDLPivot pivot = pivot::Full( ABR );

        // Apply the pivot
        const Int from = k + pivot.from[0];
        HermitianSwap( UPPER, A, k, from );
        RowSwap( pPerm, k, from );

        // a12 := a12 / sqrt(alpha11)
        const Base<F> delta11 = Sqrt(ABR.GetRealPart(0,0));
        const Base<F> delta11Inv = Base<F>(1)/delta11;
        ABR.Set(0,0,delta11);
        auto a12 = ViewRange( ABR, 0, 1, 1, n-k );
        Scale( delta11Inv, a12 );

        // A22 := A22 - a12' a12
        auto A22 = ViewRange( ABR, 1, 1, n-k, n-k );
        // NOTE: This is silly, but currently necessary
        Conjugate( a12 );
        Her( UPPER, -F(1), a12, A22 );
        Conjugate( a12 );
    }
}

template<typename F,Dist UPerm>
inline void
UUnblockedPivoted( DistMatrix<F>& A, DistMatrix<Int,UPerm,STAR>& pPerm )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::UUnblockedPivoted");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Grid() != pPerm.Grid() )
            LogicError("A and pPerm must share the same grid");
    )
    const Int n = A.Height();

    // Initialize the permutation to the identity
    pPerm.Resize( n, 1 );
    for( Int iLoc=0; iLoc<pPerm.LocalHeight(); ++iLoc )
        pPerm.SetLocal( iLoc, 0, pPerm.GlobalRow(iLoc) );

    for( Int k=0; k<n; ++k )
    {
        // Determine the pivot 
        auto ABR = ViewRange( A, k, k, n, n );
        const LDLPivot pivot = pivot::Full( ABR );

        // Apply the pivot
        const Int from = k + pivot.from[0];
        HermitianSwap( UPPER, A, k, from );
        RowSwap( pPerm, k, from );

        // a12 := a12 / sqrt(alpha11)
        const Base<F> delta11 = Sqrt(ABR.GetRealPart(0,0));
        const Base<F> delta11Inv = Base<F>(1)/delta11;
        ABR.Set(0,0,delta11);
        auto a12 = ViewRange( ABR, 0, 1, 1,   n-k );
        Scale( delta11Inv, a12 );

        // A22 := A22 - a12' a12
        auto A22 = ViewRange( ABR, 1, 1, n-k, n-k );
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
( Matrix<F>& A, Matrix<Int>& pPerm, 
  Matrix<F>& X, Matrix<F>& Y, Int bsize, Int off=0 )
{
    DEBUG_ONLY(CallStackEntry cse("cholesky::UPanelPivoted"))
    const Int n = A.Height();
    DEBUG_ONLY(
        if( A.Width() != n )
            LogicError("A must be square");
        if( pPerm.Height() != n || pPerm.Width() != 1 )
            LogicError("permutation vector is the wrong size");
    )
    auto ABR = ViewRange( A, off, off, n, n );
    Zeros( X, n-off, bsize );
    Zeros( Y, n-off, bsize );

    for( Int k=0; k<bsize; ++k )
    {
        // Determine the pivot 
        auto XTrail = ViewRange( X, k, 0, n-off, k );
        auto YTrail = ViewRange( Y, k, 0, n-off, k );
        auto ATrail = ViewRange( A, off+k, off+k, n, n );
        const auto pivot = pivot::PanelFull( ATrail, XTrail, YTrail );
        const Int from = (off+k) + pivot.from[0];
        const Int to = off+k;

        // Apply the pivot
        HermitianSwap( UPPER, A, to, from );
        RowSwap( XTrail, 0, pivot.from[0] );
        RowSwap( YTrail, 0, pivot.from[0] );
        RowSwap( pPerm, off+k, from );

        // Update ABR(k,k:end) -= X(k,0:k-1) Y(k:n-off-1,0:k-1)^T
        auto x10 = LockedViewRange( X,   k, 0, k+1,   k     );
        auto YB0 = LockedViewRange( Y,   k, 0, n-off, k     );
        auto a1R =       ViewRange( ABR, k, k, k+1,   n-off );
        Gemv( NORMAL, F(-1), YB0, x10, F(1), a1R );
        a1R.MakeReal(0,0);

        // a12 := a12 / sqrt(alpha11)
        const Base<F> delta11 = Sqrt(ABR.GetRealPart(k,k));
        const Base<F> delta11Inv = Base<F>(1)/delta11;
        ABR.Set(k,k,delta11);
        auto a12 = ViewRange( ABR, k, k+1, k+1, n-off );
        Scale( delta11Inv, a12 );

        // Store x21 := a12' and y21 := a12^T
        auto x21 = ViewRange( X, k+1, k, n-off, k+1 );
        auto y21 = ViewRange( Y, k+1, k, n-off, k+1 );
        Adjoint( a12, x21 );
        Transpose( a12, y21 );
    }
}

template<typename F,Dist UPerm>
inline void
UPanelPivoted
( DistMatrix<F>& A, 
  DistMatrix<Int,UPerm,STAR>& pPerm, 
  DistMatrix<F,MC,STAR>& X, 
  DistMatrix<F,MR,STAR>& Y, Int bsize, Int off=0 )
{
    DEBUG_ONLY(CallStackEntry cse("cholesky::UPanelPivoted"))
    const Int n = A.Height();
    DEBUG_ONLY(
        if( A.Width() != n )
            LogicError("A must be square");
        if( pPerm.Height() != n || pPerm.Width() != 1 )
            LogicError("pivot vector is the wrong size");
    )
    auto ABR = ViewRange( A, off, off, n, n );
    X.AlignWith( ABR );
    Y.AlignWith( ABR );
    Zeros( X, n-off, bsize );
    Zeros( Y, n-off, bsize );

    for( Int k=0; k<bsize; ++k )
    {
        // Determine the pivot
        auto XTrail = ViewRange( X, k, 0, n-off, k );
        auto YTrail = ViewRange( Y, k, 0, n-off, k );
        auto ATrail = ViewRange( A, off+k, off+k, n, n );
        const auto pivot = pivot::PanelFull( ATrail, XTrail, YTrail );
        const Int from = (off+k) + pivot.from[0];
        const Int to = off+k;

        // Apply the pivot
        HermitianSwap( UPPER, A, to, from );
        RowSwap( XTrail, 0, pivot.from[0] );
        RowSwap( YTrail, 0, pivot.from[0] );
        RowSwap( pPerm, off+k, from );

        // Update ABR(k,k:end) -= X(k,0:k-1) Y(k:n-off-1,0:k-1)^T
        auto a1R = ViewRange( ABR, k, k, k+1, n-off );
        if( a1R.ColAlign() == a1R.ColRank() )
        {
            auto x10 = LockedViewRange( X, k, 0, k+1,   k );
            auto YB0 = LockedViewRange( Y, k, 0, n-off, k );
            LocalGemv( NORMAL, F(-1), YB0, x10, F(1), a1R );
        }
        a1R.MakeReal(0,0);

        // a12 := a12 / sqrt(alpha11)
        const Base<F> delta11 = Sqrt(ABR.GetRealPart(k,k));
        const Base<F> delta11Inv = Base<F>(1)/delta11;
        ABR.Set(k,k,delta11);
        auto a12 = ViewRange( ABR, k, k+1, k+1, n-off );
        Scale( delta11Inv, a12 );

        // Store x21 := a12' and y21 := a12^T
        auto x21 = ViewRange( X, k+1, k, n-off, k+1 );
        auto y21 = ViewRange( Y, k+1, k, n-off, k+1 );
        Adjoint( a12, x21 );
        Transpose( a12, y21 );
    }
}

template<typename F>
inline void
UVar3( Matrix<F>& A, Matrix<Int>& pPerm )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::UVar3");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    const Int n = A.Height();

    // Initialize the permutation to the identity
    pPerm.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
        pPerm.Set( i, 0, i );

    Matrix<F> X, Y;
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        UPanelPivoted( A, pPerm, X, Y, nb, k );

        // Update the bottom-right panel
        auto XTrail = ViewRange( X, nb,   0,    n-k, nb );
        auto YTrail = ViewRange( Y, nb,   0,    n-k, nb );
        auto ATrail = ViewRange( A, k+nb, k+nb, n,   n  );
        Trrk( UPPER, NORMAL, TRANSPOSE, F(-1), XTrail, YTrail, F(1), ATrail );
    }
}

template<typename F,Dist UPerm>
inline void
UVar3( DistMatrix<F>& A, DistMatrix<Int,UPerm,STAR>& pPerm )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::UVar3");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    const Int n = A.Height();

    // Initialize the permutation to the identity
    pPerm.Resize( n, 1 );
    for( Int iLoc=0; iLoc<pPerm.LocalHeight(); ++iLoc )
        pPerm.SetLocal( iLoc, 0, pPerm.GlobalRow(iLoc) );

    DistMatrix<F,MC,STAR> X( A.Grid() );
    DistMatrix<F,MR,STAR> Y( A.Grid() );
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        UPanelPivoted( A, pPerm, X, Y, nb, k );

        // Update the bottom-right panel
        auto XTrail = ViewRange( X, nb,   0,    n-k, nb );
        auto YTrail = ViewRange( Y, nb,   0,    n-k, nb );
        auto ATrail = ViewRange( A, k+nb, k+nb, n,   n  );
        LocalTrrk( UPPER, TRANSPOSE, F(-1), XTrail, YTrail, F(1), ATrail );
    }
}

} // namespace cholesky
} // namespace elem

#endif // ifndef ELEM_CHOLESKY_UVAR3PIVOTED_HPP
