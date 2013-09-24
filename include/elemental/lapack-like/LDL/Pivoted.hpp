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

#include "elemental/blas-like/level1/Max.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level1/Swap.hpp"
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
ChoosePivot( const Matrix<F>& A, Int k, LDLPivotType pivotType, BASE(F) gamma )
{
#ifndef RELEASE
    CallStackEntry cse("ldl::ChoosePivot");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();
    if( pivotType != BUNCH_KAUFMAN_A )
        LogicError("So far, only Bunch-Kaufman Algorithm A is supported");

    const Real alpha11Abs = Abs(A.Get(k,k));
    auto a21Pair = VectorMax( LockedViewRange(A,k+1,k,n,k+1) );
    const Int r = (k+1) + a21Pair.index;
    const Real colMax = a21Pair.value;
    if( colMax == Real(0) && alpha11Abs == Real(0) )
        throw SingularMatrixException();

    if( alpha11Abs >= gamma*colMax )
        return k;

    // Find maximum off-diag value in row r (exploit symmetry)
    auto leftPair   = VectorMax( LockedViewRange(A,r,  k,r+1,r  ) );
    auto bottomPair = VectorMax( LockedViewRange(A,r+1,r,n,  r+1) );
    const Real rowMax = Max(leftPair.value,bottomPair.value);

    if( alpha11Abs >= gamma*colMax*(colMax/rowMax) )
        return k;

    if( Abs(A.Get(r,r)) >= gamma*rowMax )
        return r;

    // Default to a 2x2 pivot with k and r
    return -r;
}

template<typename F>
inline Int
ChoosePivot
( const Matrix<F>& A, const Matrix<F>& X, const Matrix<F>& Y, 
  Int offset, Int k, LDLPivotType pivotType, BASE(F) gamma )
{
#ifndef RELEASE
    CallStackEntry cse("ldl::ChoosePivot");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();
    if( pivotType != BUNCH_KAUFMAN_A )
        LogicError("So far, only Bunch-Kaufman Algorithm A is supported");

    const Int kOff = offset + k;
    auto aB1 = LockedViewRange( A, kOff, kOff, n, kOff+1 );
    auto zB1( aB1 );
    // A(kOff:n-1,kOff) -= X(k:n-offset-1,0:k-1) Y(k,0:k-1)^T
    {
        auto XBL  = LockedViewRange( X, k, 0, n-offset, k );
        auto yRow = LockedViewRange( Y, k, 0, k+1,      k );
        Gemv( NORMAL, F(-1), XBL, yRow, F(1), zB1 );
    } 

    const Real alpha11Abs = Abs(zB1.Get(0,0));
    auto a21Pair = VectorMax( LockedViewRange(zB1,1,0,n-kOff,1) );
    const Int r = (kOff+1) + a21Pair.index;
    const Real colMax = a21Pair.value;
    if( colMax == Real(0) && alpha11Abs == Real(0) )
        throw SingularMatrixException();

    if( alpha11Abs >= gamma*colMax )
        return kOff;

    // Find maximum off-diag value in row r (exploit symmetry)
    auto aLeft   = LockedViewRange( A, r, kOff, r+1, r   );
    auto aBottom = LockedViewRange( A, r, r,    n,   r+1 );
        
    auto zLeft( aLeft );
    auto zBottom( aBottom );
    auto zStrictBottom = ViewRange( zBottom, 1, 0, n-r, 1 );

    //
    // Update necessary components out-of-place
    //

    // A(r,kOff:r-1) -= X(r-offset,0:k-1) Y(k:r-offset-1,0:k-1)^T
    {
        auto xMid = LockedViewRange( X, r-offset, 0, r-offset+1, k );
        auto YBL = LockedViewRange( Y, k, 0, r-offset, k );
        Gemv( NORMAL, F(-1), YBL, xMid, F(1), zLeft );
    }

    // A(r:n-1,r) -= X(r-offset:n-offset-1,0:k-1) Y(r-offset,0:k-1)^T
    {
        auto XBL = LockedViewRange( X, r-offset, 0, n-offset, k );
        auto yRow = LockedViewRange( Y, r-offset, 0, r-offset+1, k );
        Gemv( NORMAL, F(-1), XBL, yRow, F(1), zBottom );
    } 

    auto leftPair   = VectorMax( zLeft );
    auto bottomPair = VectorMax( zStrictBottom );
    const Real rowMax = Max(leftPair.value,bottomPair.value);

    if( alpha11Abs >= gamma*colMax*(colMax/rowMax) )
        return kOff;

    if( Abs(zBottom.Get(0,0)) >= gamma*rowMax )
        return r;

    // Default to a 2x2 pivot with k and r
    return -r;
}

// A := A inv(D)
template<typename F>
inline void
SolveAgainstSymmetric2x2
( UpperOrLower uplo, const Matrix<F>& D, Matrix<F>& A, bool conjugated=false )
{
#ifndef RELEASE    
    CallStackEntry cse("SolveAgainstSymmetric2x2");
    if( A.Width() != 2 )
        LogicError("A must have width 2");
#endif
    const Int m = A.Height();
    if( uplo == LOWER )
    {
        if( conjugated )     
        {
            const F delta11 = D.GetRealPart(0,0);
            const F delta21 = D.Get(1,0);
            const F delta22 = D.GetRealPart(1,1);
            const F delta21Abs = SafeAbs( delta21 );
            const F phi21To11 = delta22 / delta21Abs;
            const F phi21To22 = delta11 / delta21Abs;
            const F phi21 = delta21 / delta21Abs;
            const F xi = (F(1)/(phi21To11*phi21To22-F(1)))/delta21Abs;

            for( Int j=0; j<m; ++j )
            {
                const F eta0 = xi*(phi21To11*A.Get(j,0)-phi21      *A.Get(j,1));
                const F eta1 = xi*(phi21To22*A.Get(j,1)-Conj(phi21)*A.Get(j,0));
                A.Set( j, 0, eta0 );
                A.Set( j, 1, eta1 );
            }
        }
        else
        {
            const F delta11 = D.Get(0,0);
            const F delta21 = D.Get(1,0);
            const F delta22 = D.Get(1,1);
            const F chi21To11 = -delta22 / delta21;
            const F chi21To22 = -delta11 / delta21;
            const F chi21 = (F(1)/(F(1)-chi21To11*chi21To22))/delta21;

            for( Int j=0; j<m; ++j )
            {
                const F eta0 = chi21*(chi21To11*A.Get(j,0)+A.Get(j,1));
                const F eta1 = chi21*(chi21To22*A.Get(j,1)+A.Get(j,0));
                A.Set( j, 0, eta0 );
                A.Set( j, 1, eta1 );
            }
        }
    }
    else
        LogicError("This option not yet supported");
}

// Unblocked sequential pivoted LDL
template<typename F>
inline void
UnblockedPivoted
( Orientation orientation, Matrix<F>& A, Matrix<Int>& p, 
  LDLPivotType pivotType=BUNCH_KAUFMAN_A,
  BASE(F) gamma=(1+Sqrt(BASE(F)(17)))/8 )
{
#ifndef RELEASE
    CallStackEntry entry("ldl::UnblockedPivoted");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( orientation == NORMAL )
        LogicError("Can only perform LDL^T or LDL^H");
#endif
    const bool conjugate = ( orientation==ADJOINT );
    const Int n = A.Height();
    p.ResizeTo( n, 1 );

    Matrix<F> Y21;

    Int k=0;
    while( k < n )
    {
        // Determine the pivot (block)
        const Int pivot = ChoosePivot( A, k, pivotType, gamma );
        const Int nb   = ( pivot >= 0 ? 1     : 2      );
        const Int from = ( pivot >= 0 ? pivot : -pivot );
        const Int to = k + (nb-1);

        // Apply the symmetric pivot
        SymmetricSwap( LOWER, A, to, from, conjugate );
        if( conjugate )
        {
            // Force the active diagonal entries to be real
            A.Set( k,    k,    A.GetRealPart(k,   k   ) );
            A.Set( to,   to,   A.GetRealPart(to,  to  ) );
            A.Set( from, from, A.GetRealPart(from,from) );
        }

        // Update trailing submatrix and store pivots
        if( nb == 1 )
        {
            // Rank-one update: A22 -= a21 inv(delta11) a21'
            const F delta11Inv = F(1)/A.Get(k,k);
            auto a21 = ViewRange( A, k+1, k,   n, k+1 );
            auto A22 = ViewRange( A, k+1, k+1, n, n   );
            Syr( LOWER, -delta11Inv, a21, A22, conjugate );
            Scale( delta11Inv, a21 );

            p.Set( k, 0, pivot );
        }
        else
        {
            // Rank-two update: A22 -= A21 inv(D11) A21'
            auto D11 = ViewRange( A, k,   k,   k+2, k+2 );
            auto A21 = ViewRange( A, k+2, k,   n,   k+2 );
            auto A22 = ViewRange( A, k+2, k+2, n,   n   );
            Y21 = A21;
            SolveAgainstSymmetric2x2( LOWER, D11, A21, conjugate );
            Trr2( LOWER, F(-1), A21, Y21, A22, conjugate );

            p.Set( k,   0, pivot );
            p.Set( k+1, 0, pivot );
        }

        k += nb;
    }
}

// We must use a lazy algorithm so that the symmetric pivoting does not move
// data from a fully-updated to partially-updated region (and vice-versa)
template<typename F>
inline void
PanelPivoted
( Orientation orientation, Matrix<F>& A, Matrix<Int>& p, 
  Matrix<F>& X, Matrix<F>& Y, Int bsize, Int offset=0,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A,
  BASE(F) gamma=(1+Sqrt(BASE(F)(17)))/8 )
{
#ifndef RELEASE
    CallStackEntry entry("ldl::PanelPivoted");
#endif
    const Int n = A.Height();
#ifndef RELEASE
    if( A.Width() != n )
        LogicError("A must be square");
    if( p.Height() != n || p.Width() != 1 )
        LogicError("pivot vector is the wrong size");
    if( orientation == NORMAL )
        LogicError("Can only perform LDL^T or LDL^H");
#endif
    const bool conjugate = ( orientation==ADJOINT );
    Zeros( X, n-offset, bsize );
    Zeros( Y, n-offset, bsize );

    Int k=0;
    while( k < bsize )
    {
        // Determine the pivot (block)
        const Int pivot = ChoosePivot( A, X, Y, offset, k, pivotType, gamma );
        const Int kOff = k + offset;
        // const Int pivot = kOff;
        const Int nb   = ( pivot >= 0 ? 1     : 2      );
        const Int from = ( pivot >= 0 ? pivot : -pivot );
        const Int to = kOff + (nb-1);
        if( k+nb > bsize )
        {
            X.ResizeTo( n-offset, bsize-1 );
            Y.ResizeTo( n-offset, bsize-1 );
            break;
        }

        // Apply the symmetric pivot
        SymmetricSwap( LOWER, A, to, from, conjugate );
        RowSwap( X, to-offset, from-offset );
        RowSwap( Y, to-offset, from-offset );

        // Update the active columns and then store the new update factors
        // TODO: Reuse updates from pivot selection where possible
        if( nb == 1 ) 
        {
            // Update A(kOff:end,kOff) -= X(k:n-offset-1,0:k-1) Y(k,0:k-1)^T
            auto XB0 = LockedViewRange( X, k, 0, n-offset, k );
            auto y10 = LockedViewRange( Y, k, 0, k+1,      k );
            auto aB1 =  ViewRange( A, kOff, kOff, n, kOff+1 );
            Gemv( NORMAL, F(-1), XB0, y10, F(1), aB1 );

            // Store x21 := a21/delta11 and y21 := a21
            const F delta11Inv = F(1)/A.Get(kOff,kOff);
            auto a21 = ViewRange( A, kOff+1, kOff, n,        kOff+1 );
            auto x21 = ViewRange( X, k+1,    k,    n-offset, k+1    );
            auto y21 = ViewRange( Y, k+1,    k,    n-offset, k+1    );
            if( conjugate )
                Conjugate( a21, y21 );
            else
                y21 = a21;
            Scale( delta11Inv, a21 );
            x21 = a21;

            p.Set( kOff, 0, pivot );
        }
        else
        {
            // Update A(k:end,k:k+1) -= X(k:n-offset-1,0:k-1) Y(k:k+1,0:k-1)^T
            // NOTE: above-diagonal entry of A is modified
            auto XB0 = LockedViewRange( X, k, 0, n-offset, k );
            auto Y10 = LockedViewRange( Y, k, 0, k+2,      k );
            auto AB1 =       ViewRange( A, kOff, kOff,   n, kOff+2 );
            Gemm( NORMAL, TRANSPOSE, F(-1), XB0, Y10, F(1), AB1 );

            // Store X21 := A21/D11 and Y21 := A21 or Y21 := Conj(A21)
            auto D11 = ViewRange( A, kOff,   kOff, kOff+2, kOff+2 );
            auto A21 = ViewRange( A, kOff+2, kOff, n,      kOff+2 );
            auto X21 = ViewRange( X, k+2, k, n-offset, k+2 );
            auto Y21 = ViewRange( Y, k+2, k, n-offset, k+2 );
            if( conjugate )
                Conjugate( A21, Y21 );
            else
                Y21 = A21;
            SolveAgainstSymmetric2x2( LOWER, D11, A21, conjugate );
            X21 = A21;

            p.Set( kOff,   0, pivot );
            p.Set( kOff+1, 0, pivot );
        }

        if( conjugate )
        {
            // Force the active diagonal entries to be real
            A.Set( kOff, kOff, A.GetRealPart(kOff,kOff) );
            A.Set( to,   to,   A.GetRealPart(to,  to  ) );
            A.Set( from, from, A.GetRealPart(from,from) );
        }

        k += nb;
    }
}

template<typename F>
inline void
Pivoted
( Orientation orientation, Matrix<F>& A, Matrix<Int>& p, 
  LDLPivotType pivotType=BUNCH_KAUFMAN_A,
  BASE(F) gamma=(1+Sqrt(BASE(F)(17)))/8 )
{
#ifndef RELEASE
    CallStackEntry entry("ldl::Pivoted");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( orientation == NORMAL )
        LogicError("Can only perform LDL^T or LDL^H");
#endif
    const Int n = A.Height();
    p.ResizeTo( n, 1 );

    Matrix<F> X, Y;
    const Int bsize = Blocksize();
    Int k=0;
    while( k < n )
    {
        const Int nbProp = Min(bsize,n-k);
        PanelPivoted( orientation, A, p, X, Y, nbProp, k, pivotType, gamma );
        const Int nb = X.Width();

        // Update the bottom-right panel
        auto X21B  = ViewRange( X, nb,   0,    n-k, nb );
        auto Y21B  = ViewRange( Y, nb,   0,    n-k, nb );
        auto A22BR = ViewRange( A, k+nb, k+nb, n,   n  );
        Trrk( LOWER, NORMAL, TRANSPOSE, F(-1), X21B, Y21B, F(1), A22BR );

        k += nb;
    }
}

} // namespace ldl
} // namespace elem

#endif // ifndef ELEM_LAPACK_LDL_PIVOTED_HPP
