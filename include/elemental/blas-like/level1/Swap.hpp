/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_SWAP_HPP
#define ELEM_BLAS_SWAP_HPP

namespace elem {

template<typename F>
inline void Swap( Orientation orientation, Matrix<F>& X, Matrix<F>& Y )
{
#ifndef RELEASE
    CallStackEntry cse("Swap");
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
        if( mX > nX )
        {
            for( Int j=0; j<nX; ++j )
                blas::Swap( mX, X.Buffer(0,j), 1, Y.Buffer(0,j), 1 );
        }
        else
        {
            for( Int i=0; i<mX; ++i )
                blas::Swap
                ( nX, X.Buffer(i,0), X.LDim(), Y.Buffer(i,0), Y.LDim() );
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
            if( conjugate )
            {
                for( Int i=0; i<mX; ++i )
                {
                    const F alpha = X.Get(i,j);
                    X.Set( i, j, Conj(Y.Get(j,i)) );
                    Y.Set( j, i, Conj(alpha)      );
                }
            }
            else
            {
                blas::Swap( mX, X.Buffer(0,j), 1, Y.Buffer(j,0), Y.LDim() );
            }
        }
    }
}

template<typename F>
inline void RowSwap( Matrix<F>& A, Int to, Int from )
{
#ifndef RELEASE
    CallStackEntry cse("RowSwap");
#endif
    if( to == from )
        return;
    const Int n = A.Width();
    auto aToRow   = ViewRange( A, to,   0, to+1,   n );
    auto aFromRow = ViewRange( A, from, 0, from+1, n );
    Swap( NORMAL, aToRow, aFromRow );
}

template<typename F>
inline void ColumnSwap( Matrix<F>& A, Int to, Int from )
{
#ifndef RELEASE
    CallStackEntry cse("RowSwap");
#endif
    if( to == from )
        return;
    const Int m = A.Width();
    auto aToCol   = ViewRange( A, 0, to,   m, to+1   );
    auto aFromCol = ViewRange( A, 0, from, m, from+1 );
    Swap( NORMAL, aToCol, aFromCol );
}

template<typename F>
inline void SymmetricSwap
( UpperOrLower uplo, Matrix<F>& A, int to, int from, bool conjugate=false )
{
#ifndef RELEASE
    CallStackEntry cse("SymmetricSwap");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
#endif
    if( to == from )
        return;
    const Int n = A.Height();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    if( uplo == LOWER )
    { 
        // Bottom swap
        if( from+1 < n )
        {
            auto aToBot   = ViewRange( A, from+1, to,   n, to+1   );
            auto aFromBot = ViewRange( A, from+1, from, n, from+1 );
            Swap( NORMAL, aToBot, aFromBot );
        }
        // Inner swap
        if( to+1 < from )
        {
            auto aToInner   = ViewRange( A, to+1, to,   from,   to+1 );
            auto aFromInner = ViewRange( A, from, to+1, from+1, from );
            Swap( orientation, aToInner, aFromInner );
        }
        // Corner swap
        if( conjugate )
            A.Set( from, to, Conj(A.Get(from,to)) );
        // Diagonal swap
        {
            const F value = A.Get(from,from);
            A.Set( from, from, A.Get(to,to) );
            A.Set( to,   to,   value        );
        }
        // Left swap
        if( to > 0 )
        {
            auto aToLeft   = ViewRange( A, to,   0, to+1,   to );
            auto aFromLeft = ViewRange( A, from, 0, from+1, to );
            Swap( NORMAL, aToLeft, aFromLeft );
        }
    }
    else
    {
        LogicError("This option not yet supported");
    }
}

} // namespace elem

#endif // ifndef ELEM_BLAS_SWAP_HPP
