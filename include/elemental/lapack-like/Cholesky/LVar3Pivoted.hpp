/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_CHOLESKY_LVAR3PIVOTED_HPP
#define ELEM_LAPACK_CHOLESKY_LVAR3PIVOTED_HPP

#include "elemental/blas-like/level1/Max.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level1/Swap.hpp"
#include "elemental/blas-like/level2/Her.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {
namespace cholesky {

namespace pivot {

template<typename F>
inline LDLPivot
Full( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("cholesky::pivot::Full"))
    const auto diagMax = DiagonalMax( A );
    LDLPivot pivot;
    pivot.nb = 1;
    pivot.from[0] = diagMax.index;
    return pivot;
}

template<typename F>
inline LDLPivot
Full( const DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("cholesky::pivot::Full"))
    const auto diagMax = DiagonalMax( A );
    LDLPivot pivot;
    pivot.nb = 1;
    pivot.from[0] = diagMax.index;
    return pivot;
}

template<typename F>
inline LDLPivot
PanelFull( const Matrix<F>& A, const Matrix<F>& X, const Matrix<F>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("cholesky::pivot::PanelFull"))
    // Form updated diagonal
    auto d = A.GetDiagonal();
    const Int height = d.Height();
    const Int k = X.Width();
    for( Int i=0; i<height; ++i )
        for( Int j=0; j<k; ++j )
            d.Update( i, 0, -X.Get(i,j)*Y.Get(i,j) );

    // Return maximum from it
    auto diagMax = VectorMax( d );
    LDLPivot pivot;
    pivot.nb = 1;
    pivot.from[0] = diagMax.index;
    return pivot;
}

template<typename F>
inline LDLPivot
PanelFull
( const DistMatrix<F>& A, 
  const DistMatrix<F,MC,STAR>& X, const DistMatrix<F,MR,STAR>& Y )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::pivot::PanelFull");
        if( A.ColAlign() != X.ColAlign() || A.RowAlign() != Y.ColAlign() )
            LogicError("A, X, and Y are not properly aligned");
    )
    // Form updated diagonal
    auto d = A.GetDiagonal();
    if( d.Participating() )
    {
        const Int dShift = d.ColShift();
        const Int XShift = X.ColShift();
        const Int YShift = Y.ColShift();
        const Int dStride = d.ColStride();
        const Int XStride = X.ColStride();
        const Int YStride = Y.ColStride();
        const Int dLocalHeight = d.LocalHeight();
        const Int k = X.Width();
        for( Int iLoc=0; iLoc<dLocalHeight; ++iLoc )
        {
            const Int i = dShift + iLoc*dStride;    
            const Int iLocX = (i-XShift) / XStride;
            const Int iLocY = (i-YShift) / YStride;
            for( Int j=0; j<k; ++j )
                d.UpdateLocal
                ( iLoc, 0, -X.GetLocal(iLocX,j)*Y.GetLocal(iLocY,j) );
            d.MakeReal( iLoc, 0 );
        }
    }

    // Return maximum from it
    auto diagMax = VectorMax( d );
    LDLPivot pivot;
    pivot.nb = 1;
    pivot.from[0] = diagMax.index;
    return pivot;
}

} // namespace pivot

template<typename F>
inline void
LUnblockedPivoted( Matrix<F>& A, Matrix<Int>& p )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::LUnblockedPivoted");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    const Int n = A.Height();
    p.Resize( n, 1 );
     
    for( Int k=0; k<n; ++k )
    {
        // Determine the pivot
        auto ABR = ViewRange( A, k, k, n, n );
        const LDLPivot pivot = pivot::Full( ABR );

        // Apply the pivot
        const Int from = k + pivot.from[0];
        HermitianSwap( LOWER, A, k, from );
        p.Set( k, 0, from );

        // a21 := a21 / sqrt(alpha11)
        const BASE(F) delta11 = Sqrt(ABR.GetRealPart(0,0));
        const BASE(F) delta11Inv = BASE(F)(1)/delta11;
        ABR.Set(0,0,delta11);
        auto a21 = ViewRange( ABR, 1, 0, n-k, 1 );
        Scale( delta11Inv, a21 );

        // A22 -= a21 a21'
        auto A22 = ViewRange( ABR, 1, 1, n-k, n-k );
        Her( LOWER, F(-1), a21, A22 );
    }
}

template<typename F>
inline void
LUnblockedPivoted( DistMatrix<F>& A, DistMatrix<Int,VC,STAR>& p )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::LUnblockedPivoted");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Grid() != p.Grid() )
            LogicError("A and p must share the same grid");
    )
    const Int n = A.Height();
    p.Resize( n, 1 );

    for( Int k=0; k<n; ++k )
    {
        // Determine the pivot
        auto ABR = ViewRange( A, k, k, n, n );
        const LDLPivot pivot = pivot::Full( ABR );

        // Apply the pivot
        const Int from = k + pivot.from[0];
        HermitianSwap( LOWER, A, k, from );
        p.Set( k, 0, from );

        // a21 := a21 / sqrt(alpha11)
        const BASE(F) delta11 = Sqrt(ABR.GetRealPart(0,0));
        const BASE(F) delta11Inv = BASE(F)(1)/delta11;
        ABR.Set(0,0,delta11);
        auto a21 = ViewRange( ABR, 1, 0, n-k, 1 );
        Scale( delta11Inv, a21 );

        // A22 -= a21 a21'
        auto A22 = ViewRange( ABR, 1, 1, n-k, n-k );
        Her( LOWER, F(-1), a21, A22 );
    }
}

// We must use a lazy algorithm so that the symmetric pivoting does not move
// data from a fully-updated to partially-updated region (and vice-versa)
template<typename F>
inline void
LPanelPivoted
( Matrix<F>& A, Matrix<Int>& p, 
  Matrix<F>& X, Matrix<F>& Y, Int bsize, Int off=0 )
{
    DEBUG_ONLY(CallStackEntry cse("cholesky::LPanelPivoted"))
    const Int n = A.Height();
    DEBUG_ONLY(
        if( A.Width() != n )
            LogicError("A must be square");
        if( p.Height() != n || p.Width() != 1 )
            LogicError("pivot vector is the wrong size");
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
        HermitianSwap( LOWER, A, to, from );
        RowSwap( XTrail, 0, pivot.from[0] );
        RowSwap( YTrail, 0, pivot.from[0] );
        p.Set( off+k, 0, from );

        // Update ABR(k:end,k) -= X(k:n-off-1,0:k-1) Y(k,0:k-1)^T
        auto XB0 = LockedViewRange( X,   k, 0, n-off, k   );
        auto y10 = LockedViewRange( Y,   k, 0, k+1,   k   );
        auto aB1 =       ViewRange( ABR, k, k, n-off, k+1 );
        Gemv( NORMAL, F(-1), XB0, y10, F(1), aB1 );
        aB1.MakeReal(0,0);

        // a21 := a21 / sqrt(alpha11)
        const BASE(F) delta11 = Sqrt(ABR.GetRealPart(k,k));
        const BASE(F) delta11Inv = BASE(F)(1)/delta11;
        ABR.SetRealPart(k,k,delta11);
        auto a21 = ViewRange( ABR, k+1, k, n-off, k+1 );
        Scale( delta11Inv, a21 );

        // Store x21 := a21 and y21 := conj(a21)
        auto x21 = ViewRange( X, k+1, k, n-off, k+1 );
        auto y21 = ViewRange( Y, k+1, k, n-off, k+1 );
        Conjugate( a21, y21 );
        x21 = a21;
    }
}

template<typename F>
inline void
LPanelPivoted
( DistMatrix<F>& A, DistMatrix<Int,VC,STAR>& p, 
  DistMatrix<F,MC,STAR>& X, DistMatrix<F,MR,STAR>& Y, Int bsize, Int off=0 )
{
    DEBUG_ONLY(CallStackEntry cse("cholesky::LPanelPivoted"))
    const Int n = A.Height();
    DEBUG_ONLY(
        if( A.Width() != n )
            LogicError("A must be square");
        if( p.Height() != n || p.Width() != 1 )
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
        HermitianSwap( LOWER, A, to, from );
        RowSwap( XTrail, 0, pivot.from[0] );
        RowSwap( YTrail, 0, pivot.from[0] );
        p.Set( off+k, 0, from );

        // ABR(k:end,k) -= X(k:n-off-1,0:k-1) Y(k,0:k-1)^T
        auto aB1 = ViewRange( ABR, k, k, n-off, k+1 );
        if( aB1.RowAlign() == aB1.RowRank() )
        {
            auto XB0 = LockedViewRange( X, k, 0, n-off, k );
            auto y10 = LockedViewRange( Y, k, 0, k+1,   k );
            LocalGemv( NORMAL, F(-1), XB0, y10, F(1), aB1 );
        }
        aB1.MakeReal(0,0);

        // a21 := a21 / sqrt(alpha11)
        const BASE(F) delta11 = Sqrt(ABR.GetRealPart(k,k));
        const BASE(F) delta11Inv = BASE(F)(1)/delta11;
        ABR.SetRealPart(k,k,delta11);
        auto a21 = ViewRange( ABR, k+1, k, n-off, k+1 );
        Scale( delta11Inv, a21 );

        // Store x21 := a21 and y21 := conj(a21)
        auto x21 = ViewRange( X, k+1, k, n-off, k+1 );
        auto y21 = ViewRange( Y, k+1, k, n-off, k+1 );
        Conjugate( a21, y21 );
        x21 = a21;
    }
}

template<typename F>
inline void
LVar3( Matrix<F>& A, Matrix<Int>& p )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::LVar3");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    const Int n = A.Height();
    p.Resize( n, 1 );

    Matrix<F> X, Y;
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        LPanelPivoted( A, p, X, Y, nb, k );

        // Update the bottom-right panel
        auto XTrail = ViewRange( X, nb,   0,    n-k, nb );
        auto YTrail = ViewRange( Y, nb,   0,    n-k, nb );
        auto ATrail = ViewRange( A, k+nb, k+nb, n,   n  );
        Trrk( LOWER, NORMAL, TRANSPOSE, F(-1), XTrail, YTrail, F(1), ATrail );
    }
}

template<typename F>
inline void
LVar3( DistMatrix<F>& A, DistMatrix<Int,VC,STAR>& p )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::LVar3");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    const Int n = A.Height();
    p.Resize( n, 1 );

    DistMatrix<F,MC,STAR> X( A.Grid() );
    DistMatrix<F,MR,STAR> Y( A.Grid() );
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        LPanelPivoted( A, p, X, Y, nb, k );

        // Update the bottom-right panel
        auto XTrail = ViewRange( X, nb,   0,    n-k, nb );
        auto YTrail = ViewRange( Y, nb,   0,    n-k, nb );
        auto ATrail = ViewRange( A, k+nb, k+nb, n,   n  );
        LocalTrrk( LOWER, TRANSPOSE, F(-1), XTrail, YTrail, F(1), ATrail );
    }
}

} // namespace cholesky
} // namespace elem

#endif // ifndef ELEM_LAPACK_CHOKESY_LVAR3PIVOTED_HPP
