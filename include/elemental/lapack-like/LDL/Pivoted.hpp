/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_LDL_PIVOTED_HPP
#define ELEM_LAPACK_LDL_PIVOTED_HPP

#include "elemental/blas-like/level1/Conjugate.hpp"
#include "elemental/blas-like/level1/Max.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level1/Swap.hpp"
#include "elemental/blas-like/level1/Symmetric2x2Solve.hpp"
#include "elemental/blas-like/level2/Syr.hpp"
#include "elemental/blas-like/level2/Trr.hpp"
#include "elemental/blas-like/level2/Trr2.hpp"
#include "elemental/matrices/Zeros.hpp"

// TODO: Reference LAPACK's dsytf2 and zhetf2

namespace elem {
namespace ldl {

// TODO: Add support for Algorithm D, Bunch-Parlett, etc.. 
// Currently using Algorithm A Bunch-Kaufman.
template<typename F>
inline Int
ChoosePivot( const Matrix<F>& A, LDLPivotType pivotType, Base<F> gamma )
{
#ifndef RELEASE
    CallStackEntry cse("ldl::ChoosePivot");
#endif
    typedef Base<F> Real;
    const Int n = A.Height();
    if( pivotType != BUNCH_KAUFMAN_A )
        LogicError("So far, only Bunch-Kaufman Algorithm A is supported");

    const Real alpha11Abs = Abs(A.Get(0,0));
    auto a21Pair = VectorMax( LockedViewRange(A,1,0,n,1) );
    const Int r = a21Pair.index + 1;
    const Real colMax = a21Pair.value;
    if( colMax == Real(0) && alpha11Abs == Real(0) )
        throw SingularMatrixException();

    if( alpha11Abs >= gamma*colMax )
        return 0;

    // Find maximum off-diag value in row r (exploit symmetry)
    auto leftPair   = VectorMax( LockedViewRange(A,r,  0,r+1,r  ) );
    auto bottomPair = VectorMax( LockedViewRange(A,r+1,r,n,  r+1) );
    const Real rowMax = Max(leftPair.value,bottomPair.value);

    if( alpha11Abs >= gamma*colMax*(colMax/rowMax) )
        return 0;

    if( Abs(A.Get(r,r)) >= gamma*rowMax )
        return r;

    // Default to a 2x2 pivot with 0 and r
    return -r;
}

template<typename F>
inline Int
ChoosePivot( const DistMatrix<F>& A, LDLPivotType pivotType, Base<F> gamma )
{
#ifndef RELEASE
    CallStackEntry cse("ldl::ChoosePivot");
#endif
    typedef Base<F> Real;
    const Int n = A.Height();
    if( pivotType != BUNCH_KAUFMAN_A )
        LogicError("So far, only Bunch-Kaufman Algorithm A is supported");

    const Real alpha11Abs = Abs(A.Get(0,0));
    auto a21Pair = VectorMax( LockedViewRange(A,1,0,n,1) );
    const Int r = a21Pair.index + 1;
    const Real colMax = a21Pair.value;
    if( colMax == Real(0) && alpha11Abs == Real(0) )
        throw SingularMatrixException();

    if( alpha11Abs >= gamma*colMax )
        return 0;

    // Find maximum off-diag value in row r (exploit symmetry)
    auto leftPair   = VectorMax( LockedViewRange(A,r,  0,r+1,r  ) );
    auto bottomPair = VectorMax( LockedViewRange(A,r+1,r,n,  r+1) );
    const Real rowMax = Max(leftPair.value,bottomPair.value);

    if( alpha11Abs >= gamma*colMax*(colMax/rowMax) )
        return 0;

    if( Abs(A.Get(r,r)) >= gamma*rowMax )
        return r;

    // Default to a 2x2 pivot with 0 and r
    return -r;
}

template<typename F>
inline Int
ChoosePivot
( const Matrix<F>& A, const Matrix<F>& X, const Matrix<F>& Y, 
  Int k, LDLPivotType pivotType, Base<F> gamma )
{
#ifndef RELEASE
    CallStackEntry cse("ldl::ChoosePivot");
#endif
    typedef Base<F> Real;
    const Int n = A.Height();
    if( pivotType != BUNCH_KAUFMAN_A )
        LogicError("So far, only Bunch-Kaufman Algorithm A is supported");

    auto aB1 = LockedViewRange( A, k, k, n, k+1 );
    auto zB1( aB1 );
    // A(k:n-1,k) -= X(k:n-1,0:k-1) Y(k,0:k-1)^T
    {
        auto XBL  = LockedViewRange( X, k, 0, n,   k );
        auto yRow = LockedViewRange( Y, k, 0, k+1, k );
        Gemv( NORMAL, F(-1), XBL, yRow, F(1), zB1 );
    } 

    const Real alpha11Abs = Abs(zB1.Get(0,0));
    auto a21Pair = VectorMax( LockedViewRange(zB1,1,0,n-k,1) );
    const Int r = a21Pair.index + (k+1);
    const Real colMax = a21Pair.value;
    if( colMax == Real(0) && alpha11Abs == Real(0) )
        throw SingularMatrixException();

    if( alpha11Abs >= gamma*colMax )
        return k;

    // Find maximum off-diag value in row r (exploit symmetry)
    auto aLeft   = LockedViewRange( A, r, k, r+1, r   );
    auto aBottom = LockedViewRange( A, r, r, n,   r+1 );
        
    auto zLeft( aLeft );
    auto zBottom( aBottom );
    auto zStrictBottom = ViewRange( zBottom, 1, 0, n-r, 1 );

    //
    // Update necessary components out-of-place
    //

    // A(r,k:r-1) -= X(r,0:k-1) Y(k:r-1,0:k-1)^T
    {
        auto xMid = LockedViewRange( X, r, 0, r+1, k );
        auto YBL = LockedViewRange( Y, k, 0, r, k );
        Gemv( NORMAL, F(-1), YBL, xMid, F(1), zLeft );
    }

    // A(r:n-1,r) -= X(r:n-1,0:k-1) Y(r,0:k-1)^T
    {
        auto XBL = LockedViewRange( X, r, 0, n, k );
        auto yRow = LockedViewRange( Y, r, 0, r+1, k );
        Gemv( NORMAL, F(-1), XBL, yRow, F(1), zBottom );
    } 

    auto leftPair   = VectorMax( zLeft );
    auto bottomPair = VectorMax( zStrictBottom );
    const Real rowMax = Max(leftPair.value,bottomPair.value);

    if( alpha11Abs >= gamma*colMax*(colMax/rowMax) )
        return k;

    if( Abs(zBottom.Get(0,0)) >= gamma*rowMax )
        return r;

    // Default to a 2x2 pivot with k and r
    return -r;
}

template<typename F>
inline Int
ChoosePivot
( const DistMatrix<F>& A, 
  const DistMatrix<F,MC,STAR>& X, const DistMatrix<F,MR,STAR>& Y, 
  Int k, LDLPivotType pivotType, Base<F> gamma )
{
#ifndef RELEASE
    CallStackEntry cse("ldl::ChoosePivot");
#endif
    typedef Base<F> Real;
    const Int n = A.Height();
    if( pivotType != BUNCH_KAUFMAN_A )
        LogicError("So far, only Bunch-Kaufman Algorithm A is supported");
    if( A.ColAlign() != X.ColAlign() || A.RowAlign() != Y.ColAlign() )
        LogicError("X and Y were not properly aligned with A");

    auto aB1 = LockedViewRange( A, k, k, n, k+1 );
    auto zB1( aB1 );
    // A(k:n-1,k) -= X(k:n-1,0:k-1) Y(k,0:k-1)^T
    if( aB1.RowAlign() == aB1.RowRank() )
    {
        auto XBL  = LockedViewRange( X, k, 0, n,   k );
        auto yRow = LockedViewRange( Y, k, 0, k+1, k );
        LocalGemv( NORMAL, F(-1), XBL, yRow, F(1), zB1 );
    } 

    const Real alpha11Abs = Abs(zB1.Get(0,0));
    auto a21Pair = VectorMax( LockedViewRange(zB1,1,0,n-k,1) );
    const Int r = a21Pair.index + (k+1);
    const Real colMax = a21Pair.value;
    if( colMax == Real(0) && alpha11Abs == Real(0) )
        throw SingularMatrixException();

    if( alpha11Abs >= gamma*colMax )
        return k;

    // Find maximum off-diag value in row r (exploit symmetry)
    auto aLeft   = LockedViewRange( A, r, k, r+1, r   );
    auto aBottom = LockedViewRange( A, r, r, n,   r+1 );
        
    auto zLeft( aLeft );
    auto zBottom( aBottom );
    auto zStrictBottom = ViewRange( zBottom, 1, 0, n-r, 1 );

    //
    // Update necessary components out-of-place
    //

    // A(r,k:r-1) -= X(r,0:k-1) Y(k:r-1,0:k-1)^T
    if( aLeft.ColAlign() == aLeft.ColRank() )
    {
        auto xMid = LockedViewRange( X, r, 0, r+1, k );
        auto YBL = LockedViewRange( Y, k, 0, r, k );
        LocalGemv( NORMAL, F(-1), YBL, xMid, F(1), zLeft );
    }

    // A(r:n-1,r) -= X(r:n-1,0:k-1) Y(r,0:k-1)^T
    if( aBottom.RowAlign() == aBottom.RowRank() )
    {
        auto XBL = LockedViewRange( X, r, 0, n, k );
        auto yRow = LockedViewRange( Y, r, 0, r+1, k );
        LocalGemv( NORMAL, F(-1), XBL, yRow, F(1), zBottom );
    } 

    auto leftPair   = VectorMax( zLeft );
    auto bottomPair = VectorMax( zStrictBottom );
    const Real rowMax = Max(leftPair.value,bottomPair.value);

    if( alpha11Abs >= gamma*colMax*(colMax/rowMax) )
        return k;

    if( Abs(zBottom.Get(0,0)) >= gamma*rowMax )
        return r;

    // Default to a 2x2 pivot with k and r
    return -r;
}

// Unblocked sequential pivoted LDL
template<typename F>
inline void
UnblockedPivoted
( Matrix<F>& A, Matrix<F>& dSub, Matrix<Int>& p, bool conjugate=false,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A,
  Base<F> gamma=(1+Sqrt(Base<F>(17)))/8 )
{
#ifndef RELEASE
    CallStackEntry entry("ldl::UnblockedPivoted");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
#endif
    const Int n = A.Height();
    if( n == 0 )
    {
        dSub.ResizeTo( 0, 1 );
        p.ResizeTo( 0, 1 );
        return;
    }
    Zeros( dSub, n-1, 1 );
    p.ResizeTo( n, 1 );
     
    Matrix<F> Y21;

    Int k=0;
    while( k < n )
    {
        // Determine the pivot (block)
        auto ABR = LockedViewRange( A, k, k, n, n );
        const Int pivot = ChoosePivot( ABR, pivotType, gamma );
        const Int nb   = ( pivot >= 0 ? 1       : 2      );
        const Int from = ( pivot >= 0 ? k+pivot : k-pivot );
        const Int to = k + (nb-1);

        // Apply the symmetric pivot
        SymmetricSwap( LOWER, A, to, from, conjugate );
        if( conjugate )
        {
            // Force the active diagonal entries to be real
            A.MakeReal( k,    k    );
            A.MakeReal( to,   to   );
            A.MakeReal( from, from );
        }

        // Update trailing submatrix and store pivots
        if( nb == 1 )
        {
            // Rank-one update: A22 -= a21 inv(delta11) a21'
            const F delta11Inv = F(1)/ABR.Get(0,0);
            auto a21 = ViewRange( ABR, 1, 0, n-k, 1   );
            auto A22 = ViewRange( ABR, 1, 1, n-k, n-k );
            Syr( LOWER, -delta11Inv, a21, A22, conjugate );
            Scale( delta11Inv, a21 );

            p.Set( k, 0, from );
        }
        else
        {
            // Rank-two update: A22 -= A21 inv(D11) A21'
            auto D11 = ViewRange( ABR, 0, 0, 2,   2   );
            auto A21 = ViewRange( ABR, 2, 0, n-k, 2   );
            auto A22 = ViewRange( ABR, 2, 2, n-k, n-k );
            Y21 = A21;
            Symmetric2x2Solve( RIGHT, LOWER, D11, A21, conjugate );
            Trr2( LOWER, F(-1), A21, Y21, A22, conjugate );

            // Only leave the main diagonal of D in A, so that routines like
            // Trsm can still be used. Thus, return the subdiagonal.
            dSub.Set( k, 0, D11.Get(1,0) );
            D11.Set( 1, 0, 0 );
            p.Set( k,   0, k    );
            p.Set( k+1, 0, from );
        }

        k += nb;
    }
}

template<typename F>
inline void
UnblockedPivoted
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& dSub, DistMatrix<Int,VC,STAR>& p, 
  bool conjugate=false, LDLPivotType pivotType=BUNCH_KAUFMAN_A,
  Base<F> gamma=(1+Sqrt(Base<F>(17)))/8 )
{
#ifndef RELEASE
    CallStackEntry entry("ldl::UnblockedPivoted");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( A.Grid() != dSub.Grid() || dSub.Grid() != p.Grid() )
        LogicError("A, dSub, and p must share the same grid");
#endif
    const Int n = A.Height();
    dSub.AlignWithDiagonal( A, -1 );
    Zeros( dSub, n-1, 1 );
    p.ResizeTo( n, 1 );

    DistMatrix<F> Y21( A.Grid() );
    DistMatrix<F,STAR,STAR> D11_STAR_STAR( A.Grid() );

    Int k=0;
    while( k < n )
    {
        // Determine the pivot (block)
        auto ABR = LockedViewRange( A, k, k, n, n );
        const Int pivot = ChoosePivot( ABR, pivotType, gamma );
        const Int nb   = ( pivot >= 0 ? 1       : 2      );
        const Int from = ( pivot >= 0 ? k+pivot : k-pivot );
        const Int to = k + (nb-1);

        // Apply the symmetric pivot
        SymmetricSwap( LOWER, A, to, from, conjugate );
        if( conjugate )
        {
            // Force the active diagonal entries to be real
            A.MakeReal( k,    k    );
            A.MakeReal( to,   to   );
            A.MakeReal( from, from );
        }

        // Update trailing submatrix and store pivots
        if( nb == 1 )
        {
            // Rank-one update: A22 -= a21 inv(delta11) a21'
            const F delta11Inv = F(1)/ABR.Get(0,0);
            auto a21 = ViewRange( ABR, 1, 0, n-k, 1   );
            auto A22 = ViewRange( ABR, 1, 1, n-k, n-k );
            Syr( LOWER, -delta11Inv, a21, A22, conjugate );
            Scale( delta11Inv, a21 );

            p.Set( k, 0, from );
        }
        else
        {
            // Rank-two update: A22 -= A21 inv(D11) A21'
            auto D11 = ViewRange( ABR, 0, 0, 2,   2   );
            auto A21 = ViewRange( ABR, 2, 0, n-k, 2   );
            auto A22 = ViewRange( ABR, 2, 2, n-k, n-k );
            Y21 = A21;
            D11_STAR_STAR = D11;
            Symmetric2x2Solve( RIGHT, LOWER, D11_STAR_STAR, A21, conjugate );
            Trr2( LOWER, F(-1), A21, Y21, A22, conjugate );

            // Only leave the main diagonal of D in A, so that routines like
            // Trsm can still be used. Thus, return the subdiagonal.
            dSub.Set( k, 0, D11_STAR_STAR.GetLocal(1,0) );
            D11.Set( 1, 0, 0 );
            p.Set( k,   0, k    );
            p.Set( k+1, 0, from );
        }

        k += nb;
    }
}

// We must use a lazy algorithm so that the symmetric pivoting does not move
// data from a fully-updated to partially-updated region (and vice-versa)
template<typename F>
inline void
PanelPivoted
( Matrix<F>& A, Matrix<F>& dSub, Matrix<Int>& p, 
  Matrix<F>& X, Matrix<F>& Y, Int bsize, Int off=0,
  bool conjugate=false, LDLPivotType pivotType=BUNCH_KAUFMAN_A,
  Base<F> gamma=(1+Sqrt(Base<F>(17)))/8 )
{
#ifndef RELEASE
    CallStackEntry entry("ldl::PanelPivoted");
#endif
    const Int n = A.Height();
    if( n == 0 )
        return;
#ifndef RELEASE
    if( A.Width() != n )
        LogicError("A must be square");
    if( dSub.Height() != n-1 || dSub.Width() != 1 )
        LogicError("dSub is the wrong size" );
    if( p.Height() != n || p.Width() != 1 )
        LogicError("pivot vector is the wrong size");
#endif
    auto ABR = ViewRange( A, off, off, n, n );
    Zeros( X, n-off, bsize );
    Zeros( Y, n-off, bsize );

    Int k=0;
    while( k < bsize )
    {
        // Determine the pivot (block)
        const Int pivot = ChoosePivot( ABR, X, Y, k, pivotType, gamma );
        const Int nb   = ( pivot >= 0 ? 1         : 2         );
        const Int from = ( pivot >= 0 ? off+pivot : off-pivot );
        const Int to = (off+k) + (nb-1);
        if( k+nb > bsize )
        {
            X.ResizeTo( n-off, bsize-1 );
            Y.ResizeTo( n-off, bsize-1 );
            break;
        }

        // Apply the symmetric pivot
        SymmetricSwap( LOWER, A, to, from, conjugate );
        RowSwap( X, to-off, from-off );
        RowSwap( Y, to-off, from-off );

        // Update the active columns and then store the new update factors
        // TODO: Reuse updates from pivot selection where possible
        if( nb == 1 ) 
        {
            // Update ABR(k:end,k) -= X(k:n-off-1,0:k-1) Y(k,0:k-1)^T
            auto XB0 = LockedViewRange( X,   k, 0, n-off, k   );
            auto y10 = LockedViewRange( Y,   k, 0, k+1,   k   );
            auto aB1 =       ViewRange( ABR, k, k, n-off, k+1 );
            Gemv( NORMAL, F(-1), XB0, y10, F(1), aB1 );
            if( conjugate )
                ABR.MakeReal(k,k);

            // Store x21 := a21/delta11 and y21 := a21
            const F delta11Inv = F(1)/ABR.Get(k,k);
            auto a21 = ViewRange( ABR, k+1, k, n-off, k+1 );
            auto x21 = ViewRange( X,   k+1, k, n-off, k+1 );
            auto y21 = ViewRange( Y,   k+1, k, n-off, k+1 );
            if( conjugate )
                Conjugate( a21, y21 );
            else
                y21 = a21;
            Scale( delta11Inv, a21 );
            x21 = a21;

            p.Set( off+k, 0, from );
        }
        else
        {
            // Update ABR(k:end,k:k+1) -= X(k:n-off-1,0:k-1) Y(k:k+1,0:k-1)^T
            // NOTE: top-right entry of AB1 is above-diagonal
            auto XB0 = LockedViewRange( X,   k, 0, n-off, k   );
            auto Y10 = LockedViewRange( Y,   k, 0, k+2,   k   );
            auto AB1 =       ViewRange( ABR, k, k, n-off, k+2 );
            const F psi = AB1.Get(0,1);
            Gemm( NORMAL, TRANSPOSE, F(-1), XB0, Y10, F(1), AB1 );
            AB1.Set(0,1,psi);
            if( conjugate )
            {
                ABR.MakeReal(k  ,k  );
                ABR.MakeReal(k+1,k+1);
            }

            // Store X21 := A21/D11 and Y21 := A21 or Y21 := Conj(A21)
            auto D11 = ViewRange( ABR, k,   k, k+2,   k+2 );
            auto A21 = ViewRange( ABR, k+2, k, n-off, k+2 );
            auto X21 = ViewRange( X,   k+2, k, n-off, k+2 );
            auto Y21 = ViewRange( Y,   k+2, k, n-off, k+2 );
            if( conjugate )
                Conjugate( A21, Y21 );
            else
                Y21 = A21;
            Symmetric2x2Solve( RIGHT, LOWER, D11, A21, conjugate );
            X21 = A21;

            // Only leave the main diagonal of D in A, so that routines like
            // Trsm can still be used. Thus, return the subdiagonal.
            dSub.Set( off+k, 0, D11.Get(1,0) );
            D11.Set( 1, 0, 0 );
            p.Set( off+k,   0, off+k );
            p.Set( off+k+1, 0, from  );
        }

        k += nb;
    }
}

template<typename F>
inline void
PanelPivoted
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& dSub, DistMatrix<Int,VC,STAR>& p, 
  DistMatrix<F,MC,STAR>& X, DistMatrix<F,MR,STAR>& Y, Int bsize, Int off=0,
  bool conjugate=false, LDLPivotType pivotType=BUNCH_KAUFMAN_A,
  Base<F> gamma=(1+Sqrt(Base<F>(17)))/8 )
{
#ifndef RELEASE
    CallStackEntry entry("ldl::PanelPivoted");
#endif
    const Int n = A.Height();
    if( n == 0 )
        return;
#ifndef RELEASE
    if( A.Width() != n )
        LogicError("A must be square");
    if( dSub.Height() != n-1 || dSub.Width() != 1 )
        LogicError("dSub is the wrong size" );
    if( p.Height() != n || p.Width() != 1 )
        LogicError("pivot vector is the wrong size");
#endif
    auto ABR = ViewRange( A, off, off, n, n );
    X.AlignWith( ABR );
    Y.AlignWith( ABR );
    Zeros( X, n-off, bsize );
    Zeros( Y, n-off, bsize );

    DistMatrix<F,STAR,STAR> D11_STAR_STAR( A.Grid() );

    Int k=0;
    while( k < bsize )
    {
        // Determine the pivot (block)
        const Int pivot = ChoosePivot( ABR, X, Y, k, pivotType, gamma );
        const Int nb   = ( pivot >= 0 ? 1         : 2         );
        const Int from = ( pivot >= 0 ? off+pivot : off-pivot );
        const Int to = (off+k) + (nb-1);
        if( k+nb > bsize )
        {
            X.ResizeTo( n-off, bsize-1 );
            Y.ResizeTo( n-off, bsize-1 );
            break;
        }

        // Apply the symmetric pivot
        SymmetricSwap( LOWER, A, to, from, conjugate );
        RowSwap( X, to-off, from-off );
        RowSwap( Y, to-off, from-off );

        // Update the active columns and then store the new update factors
        // TODO: Reuse updates from pivot selection where possible
        if( nb == 1 ) 
        {
            // Update ABR(k:end,k) -= X(k:n-off-1,0:k-1) Y(k,0:k-1)^T
            auto aB1 = ViewRange( ABR, k, k, n-off, k+1 );
            if( aB1.RowAlign() == aB1.RowRank() )
            {
                auto XB0 = LockedViewRange( X, k, 0, n-off, k );
                auto y10 = LockedViewRange( Y, k, 0, k+1,   k );
                LocalGemv( NORMAL, F(-1), XB0, y10, F(1), aB1 );
            }
            if( conjugate )
                ABR.MakeReal(k,k);

            // Store x21 := a21/delta11 and y21 := a21
            const F delta11Inv = F(1)/ABR.Get(k,k);
            auto a21 = ViewRange( ABR, k+1, k, n-off, k+1 );
            auto x21 = ViewRange( X,   k+1, k, n-off, k+1 );
            auto y21 = ViewRange( Y,   k+1, k, n-off, k+1 );
            if( conjugate )
                Conjugate( a21, y21 );
            else
                y21 = a21;
            Scale( delta11Inv, a21 );
            x21 = a21;

            p.Set( off+k, 0, from );
        }
        else
        {
            // Update ABR(k:end,k:k+1) -= X(k:n-off-1,0:k-1) Y(k:k+1,0:k-1)^T
            // NOTE: top-right entry of AB1 is above-diagonal
            auto XB0 = LockedViewRange( X,   k, 0, n-off, k   );
            auto Y10 = LockedViewRange( Y,   k, 0, k+2,   k   );
            auto AB1 =       ViewRange( ABR, k, k, n-off, k+2 );
            // TODO: Make Get and Set local
            const F psi = AB1.Get(0,1);
            LocalGemm( NORMAL, TRANSPOSE, F(-1), XB0, Y10, F(1), AB1 );
            AB1.Set(0,1,psi);
            if( conjugate )
            {
                ABR.MakeReal(k,  k  );
                ABR.MakeReal(k+1,k+1);
            }

            // Store X21 := A21/D11 and Y21 := A21 or Y21 := Conj(A21)
            auto D11 = ViewRange( ABR, k,   k, k+2,   k+2 );
            auto A21 = ViewRange( ABR, k+2, k, n-off, k+2 );
            auto X21 = ViewRange( X,   k+2, k, n-off, k+2 );
            auto Y21 = ViewRange( Y,   k+2, k, n-off, k+2 );
            if( conjugate )
                Conjugate( A21, Y21 );
            else
                Y21 = A21;
            D11_STAR_STAR = D11;
            Symmetric2x2Solve( RIGHT, LOWER, D11_STAR_STAR, A21, conjugate );
            X21 = A21;

            // Only leave the main diagonal of D in A, so that routines like
            // Trsm can still be used. Thus, return the subdiagonal.
            dSub.Set( off+k, 0, D11_STAR_STAR.GetLocal(1,0) );
            D11.Set( 1, 0, 0 );
            p.Set( off+k,   0, off+k );
            p.Set( off+k+1, 0, from  );
        }

        k += nb;
    }
}

template<typename F>
inline void
Pivoted
( Matrix<F>& A, Matrix<F>& dSub, Matrix<Int>& p, bool conjugate=false,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A,
  Base<F> gamma=(1+Sqrt(Base<F>(17)))/8 )
{
#ifndef RELEASE
    CallStackEntry entry("ldl::Pivoted");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
#endif
    const Int n = A.Height();
    if( n == 0 )
    {
        dSub.ResizeTo( 0, 1 );
        p.ResizeTo( 0, 1 );
        return;
    }
    Zeros( dSub, n-1, 1 );
    p.ResizeTo( n, 1 );

    Matrix<F> X, Y;
    const Int bsize = Blocksize();
    Int k=0;
    while( k < n )
    {
        const Int nbProp = Min(bsize,n-k);
        PanelPivoted
        ( A, dSub, p, X, Y, nbProp, k, conjugate, pivotType, gamma );
        const Int nb = X.Width();

        // Update the bottom-right panel
        auto X21B  = ViewRange( X, nb,   0,    n-k, nb );
        auto Y21B  = ViewRange( Y, nb,   0,    n-k, nb );
        auto A22BR = ViewRange( A, k+nb, k+nb, n,   n  );
        Trrk( LOWER, NORMAL, TRANSPOSE, F(-1), X21B, Y21B, F(1), A22BR );

        k += nb;
    }
}

template<typename F>
inline void
Pivoted
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& dSub, DistMatrix<Int,VC,STAR>& p, 
  bool conjugate=false, LDLPivotType pivotType=BUNCH_KAUFMAN_A,
  Base<F> gamma=(1+Sqrt(Base<F>(17)))/8 )
{
#ifndef RELEASE
    CallStackEntry entry("ldl::Pivoted");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
#endif
    const Grid& g = A.Grid();
    const Int n = A.Height();
    if( n == 0 )
    {
        dSub.ResizeTo( 0, 1 );
        p.ResizeTo( 0, 1 );
        return;
    }
    dSub.AlignWithDiagonal( A, -1 );
    Zeros( dSub, n-1, 1 );
    p.ResizeTo( n, 1 );

    DistMatrix<F,MC,STAR> X(g);
    DistMatrix<F,MR,STAR> Y(g);
    const Int bsize = Blocksize();
    Int k=0;
    while( k < n )
    {
        const Int nbProp = Min(bsize,n-k);
        PanelPivoted
        ( A, dSub, p, X, Y, nbProp, k, conjugate, pivotType, gamma );
        const Int nb = X.Width();

        // Update the bottom-right panel
        auto X21B  = ViewRange( X, nb,   0,    n-k, nb );
        auto Y21B  = ViewRange( Y, nb,   0,    n-k, nb );
        auto A22BR = ViewRange( A, k+nb, k+nb, n,   n  );
        LocalTrrk( LOWER, TRANSPOSE, F(-1), X21B, Y21B, F(1), A22BR );

        k += nb;
    }
}

} // namespace ldl
} // namespace elem

#endif // ifndef ELEM_LAPACK_LDL_PIVOTED_HPP
