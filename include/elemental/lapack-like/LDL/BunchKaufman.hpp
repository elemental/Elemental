/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_LDL_BUNCHKAUFMAN_HPP
#define ELEM_LAPACK_LDL_BUNCHKAUFMAN_HPP

#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level2/Syr.hpp"

// TODO: Reference LAPACK's dsytf2 and zhetf2

namespace elem {
namespace ldl {

template<typename F>
inline ValueInt<BASE(F)>
FindMax( const Matrix<F>& x )
{
    typedef BASE(F) Real;
    const Int m = x.Height();
    const Int n = x.Width();

    ValueInt<Real> pivot;
    pivot.index = 0;
    pivot.value = 0;
    if( n == 1 )
    {
        for( Int i=0; i<m; ++i )
        {
            const Real abs = Abs(x.Get(i,0));
            if( abs > pivot.value )
            {
                pivot.index = i;
                pivot.value = abs;
            }
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            const Real abs = Abs(x.Get(0,j));
            if( abs > pivot.value )
            {
                pivot.index = j;
                pivot.value = abs;
            }
        }
    }
    return pivot;
}

template<typename F>
inline void Swap( Orientation orientation, Matrix<F>& X, Matrix<F>& Y )
{
#ifndef RELEASE
    CallStackEntry cse("ldl::Swap");
#endif
    const Int mX = X.Height();
    const Int nX = X.Width();

    if( orientation == NORMAL )
    {
#ifndef RELEASE
        if( Y.Height() != mX || Y.Width() != nX )
            LogicError("Invalid submatrix sizes");
#endif
        // TODO: Optimize memory access patterns
        for( Int j=0; j<nX; ++j )
        {
            for( Int i=0; i<mX; ++i )
            {
                const F alpha = X.Get(i,j);    
                X.Set( i, j, Y.Get(i,j) );
                Y.Set( i, j, alpha      );
            }
        }
    }
    else
    {
        const bool conjugate = ( orientation==ADJOINT );
#ifndef RELEASE
        if( Y.Width() != mX || Y.Height() != nX )
            LogicError("Invalid submatrix sizes");
#endif
        // TODO: Optimize memory access patterns
        for( Int j=0; j<nX; ++j )
        {
            for( Int i=0; i<mX; ++i )
            {
                const F alpha = X.Get(i,j);
                if( conjugate )
                {
                    X.Set( i, j, Conj(Y.Get(j,i)) ); 
                    Y.Set( j, i, Conj(alpha)      );
                }
                else
                {
                    X.Set( i, j, Y.Get(j,i) ); 
                    Y.Set( j, i, alpha      );
                }
            }
        }
    }
}

template<typename F>
inline void ApplyPivot
( Orientation orientation, Matrix<F>& A, int to, int from )
{
#ifndef RELEASE
    CallStackEntry cse("ldl::ApplyPivot");
    if( orientation == NORMAL )
        LogicError("Invalid orientation");
#endif
    const Int n = A.Height();
    if( to != from )
    {
        // Bottom swap
        auto aToBot   = ViewRange( A, from+1, to,   n, to+1   );
        auto aFromBot = ViewRange( A, from+1, from, n, from+1 );
        Swap( NORMAL, aToBot, aFromBot );
        // Inner swap
        auto aToInner   = ViewRange( A, to+1, to,   from,   to+1 );
        auto aFromInner = ViewRange( A, from, to+1, from+1, from );
        Swap( orientation, aToInner, aFromInner );
        // Corner swap
        if( orientation == ADJOINT )
            A.Set( from, to, Conj(A.Get(from,to)) );
        // Diagonal swap
        {
            const F value = A.Get(from,from);
            A.Set( from, from, A.Get(to,to) );
            A.Set( to,   to,   value        );
        }
        // Left swap
        // NOTE: LAPACK would only swap two entries here if nb=2
        //       (otherwise no left swap)
        auto aToLeft   = ViewRange( A, to,   0, to+1,   to );
        auto aFromLeft = ViewRange( A, from, 0, from+1, to );
        Swap( NORMAL, aToLeft, aFromLeft );
    }
}

// TODO: Add support for Algorithm D. Currently using Algorithm A
template<typename F>
inline Int
ChoosePivot( const Matrix<F>& A, Matrix<Int>& p, Int k, BASE(F) gamma )
{
#ifndef RELEASE
    CallStackEntry cse("ldl::ChoosePivot");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();

    const F alpha11 = A.Get( k, k );
    const Real alpha11Abs = Abs(alpha11);
    auto a21 = LockedViewRange( A, k+1, k, n, k+1 );
    auto a21Pair = FindMax( a21 );
    const Int r = (k+1) + a21Pair.index;
    const Real colMax = a21Pair.value;
    if( colMax == Real(0) && alpha11Abs == Real(0) )
        throw SingularMatrixException();

    if( alpha11Abs >= gamma*colMax )
    {
        p.Set( k, 0, k );
    }
    else
    {
        // Find maximum off-diag value in row r (exploit symmetry)
        auto aLeft   = LockedViewRange( A, r,   0, r+1, r   );
        auto aBottom = LockedViewRange( A, r+1, r, n,   r+1 );
        auto leftPair   = FindMax( aLeft );
        auto bottomPair = FindMax( aBottom );
        const Real rowMax = Max(leftPair.value,bottomPair.value);

        if( alpha11Abs >= gamma*colMax*(colMax/rowMax) )
        {
            p.Set( k, 0, k );
        }
        else if( Abs(A.Get(r,r)) >= gamma*rowMax )
        {
            p.Set( k, 0, r );
        }
        else
        {
            p.Set( k,   0, -r );
            p.Set( k+1, 0, -r );
        }
    }

    return p.Get( k, 0 );
}

// Unblocked sequential BunchKaufman
// TODO: Better documentation of LAPACK approach to diagonal-block inversion
template<typename F>
inline void
BunchKaufman
( Orientation orientation, Matrix<F>& A, Matrix<Int>& p, 
  BASE(F) gamma=(1+Sqrt(BASE(F)(17)))/8 )
{
#ifndef RELEASE
    CallStackEntry entry("ldl::BunchKaufman");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( orientation == NORMAL )
        LogicError("Can only perform LDL^T or LDL^H");
#endif
    const bool conjugate = ( orientation==ADJOINT );
    const Int n = A.Height();
    p.ResizeTo( n, 1 );

    Int k=0;
    while( k < n )
    {
        // Determine the pivot (block)
        const Int pivot = ChoosePivot( A, p, k, gamma );
        const Int nb   = ( pivot >= 0 ? 1     : 2      );
        const Int from = ( pivot >= 0 ? pivot : -pivot );
        const Int to = k + (nb-1);

        // Apply the symmetric pivot
        ApplyPivot( orientation, A, to, from );

        if( conjugate )
        {
            // Force the active diagonal entries to be real
            A.Set( k,    k,    A.GetRealPart(k,   k   ) );
            A.Set( to,   to,   A.GetRealPart(to,  to  ) );
            A.Set( from, from, A.GetRealPart(from,from) );
        }

        // Update trailing submatrix
        if( nb == 1 )
        {
            const F delta11 = F(1)/A.Get(k,k);
            auto a21 = ViewRange( A, k+1, k,   n, k+1 );
            auto A22 = ViewRange( A, k+1, k+1, n, n   );
            Syr( LOWER, -delta11, a21, A22, conjugate );
            Scale( delta11, a21 );
        }
        else
        {
            // Rank-two update
            // Alternative: compute W21=A21/D11 and then apply A22 -= W21 A21',
            // where A22=A(k+2:end,k+2:end)
            if( conjugate )
            {
                const F delta11 = A.GetRealPart(k,k);
                const F delta21 = A.Get(k+1,k);
                const F delta22 = A.GetRealPart(k+1,k+1);
                const F delta21Abs = SafeAbs( delta21 );
                const F phi21To11 = delta22 / delta21Abs;
                const F phi21To22 = delta11 / delta21Abs;
                const F phi21 = delta21 / delta21Abs;
                const F xi = (F(1)/(phi21To11*phi21To22-F(1)))/delta21Abs;

                for( Int j=k+2; j<n; ++j )
                {
                    const F etak = 
                        xi*(phi21To11*A.Get(j,k  )-phi21      *A.Get(j,k+1));
                    const F etakp1 = 
                        xi*(phi21To22*A.Get(j,k+1)-Conj(phi21)*A.Get(j,k  ));
                    for( Int i=j; i<n; ++i )
                        A.Update(i,j,-A.Get(i,k  )*Conj(etak  )
                                     -A.Get(i,k+1)*Conj(etakp1));
                    A.Set( j, k,   etak               );
                    A.Set( j, k+1, etakp1             );
                    A.Set( j, j,   A.GetRealPart(j,j) );
                }
            }
            else
            {
                const F delta11 = A.Get(k,k);
                const F delta21 = A.Get(k+1,k);
                const F delta22 = A.Get(k+1,k+1);
                const F chi21To11 = -delta22 / delta21;
                const F chi21To22 = -delta11 / delta21;
                const F chi21 = (F(1)/(F(1)-chi21To11*chi21To22))/delta21;

                for( Int j=k+2; j<n; ++j )
                {
                    const F etak = chi21*(chi21To11*A.Get(j,k)+A.Get(j,k+1));
                    const F etakp1 = chi21*(chi21To22*A.Get(j,k+1)+A.Get(j,k));
                    for( Int i=j; i<n; ++i )
                        A.Update(i,j,-A.Get(i,k)*etak-A.Get(i,k+1)*etakp1);
                    A.Set( j, k,   etak   );
                    A.Set( j, k+1, etakp1 );
                }
            }
        }

        k += nb;
    }
}

// TODO: Blocked algorithm

} // namespace ldl
} // namespace elem

#endif // ifndef ELEM_LAPACK_LDL_BUNCHKAUFMAN_HPP
