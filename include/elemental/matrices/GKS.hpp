/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_GKS_HPP
#define MATRICES_GKS_HPP

// The Golub Klema Stewart matrix is upper-triangular with 1/sqrt(j) on its 
// j'th diagonal entry and -1/sqrt(j) elsewhere in the upper triangle.
// 
// It was originally introduced as an example of where greedy RRQR fails.

namespace elem {

template<typename F>
inline void
GKS( Matrix<F>& A, int n )
{
#ifndef RELEASE
    PushCallStack("GKS");
#endif
    A.ResizeTo( n, n );
    MakeGKS( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V>
inline void
GKS( DistMatrix<F,U,V>& A, int n )
{
#ifndef RELEASE
    PushCallStack("GKS");
#endif
    A.ResizeTo( n, n );
    MakeGKS( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
MakeGKS( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("MakeGKS");
#endif
    const int m = A.Height();
    const int n = A.Width();
    if( m != n )
        throw std::logic_error("Cannot make a non-square matrix GKS");

    MakeZeros( A );
    for( int j=0; j<n; ++j )
    {
        const F jDiag = F(1)/Sqrt(F(j));
        for( int i=0; i<j; ++i )
            A.Set( i, j, -jDiag );
        A.Set( j, j, jDiag );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V>
inline void
MakeGKS( DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("MakeGKS");
#endif
    const int m = A.Height();
    const int n = A.Width();
    if( m != n )
        throw std::logic_error("Cannot make a non-square matrix GKS");

    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();
    const int colStride = A.ColStride();
    const int rowStride = A.RowStride();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const int j = rowShift + jLocal*rowStride;
        const F jDiag = F(1)/Sqrt(F(j));
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = colShift + iLocal*colStride;
            if( i < j )
                A.SetLocal( iLocal, jLocal, -jDiag );
            else if( i == j )
                A.SetLocal( iLocal, jLocal, jDiag );
            else
                A.SetLocal( iLocal, jLocal, 0 );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef MATRICES_GKS_HPP
