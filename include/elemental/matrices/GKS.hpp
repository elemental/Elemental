/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_GKS_HPP
#define ELEM_MATRICES_GKS_HPP

// The Golub Klema Stewart matrix is upper-triangular with 1/sqrt(j) on its 
// j'th diagonal entry and -1/sqrt(j) elsewhere in the upper triangle.
// 
// It was originally introduced as an example of where greedy RRQR fails.

namespace elem {

template<typename F> 
inline void
MakeGKS( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("MakeGKS");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square matrix GKS");

    MakeZeros( A );
    for( Int j=0; j<n; ++j )
    {
        const F jDiag = F(1)/Sqrt(F(j));
        for( Int i=0; i<j; ++i )
            A.Set( i, j, -jDiag );
        A.Set( j, j, jDiag );
    }
}

template<typename F,Distribution U,Distribution V>
inline void
MakeGKS( DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry cse("MakeGKS");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square matrix GKS");

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int colShift = A.ColShift();
    const Int rowShift = A.RowShift();
    const Int colStride = A.ColStride();
    const Int rowStride = A.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        const F jDiag = F(1)/Sqrt(F(j));
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
            if( i < j )
                A.SetLocal( iLoc, jLoc, -jDiag );
            else if( i == j )
                A.SetLocal( iLoc, jLoc, jDiag );
            else
                A.SetLocal( iLoc, jLoc, 0 );
        }
    }
}

template<typename F>
inline Matrix<F>
GKS( Int n )
{
    Matrix<F> A( n, n );
    MakeGKS( A );
    return A;
}

template<typename F,Distribution U,Distribution V>
inline void
GKS( DistMatrix<F,U,V>& A, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("GKS");
#endif
    A.ResizeTo( n, n );
    MakeGKS( A );
}

template<typename F,Distribution U=MC,Distribution V=MR>
inline DistMatrix<F,U,V>
GKS( const Grid& g, Int n )
{
    DistMatrix<F,U,V> A( n, n, g );
    MakeGKS( A );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_GKS_HPP
