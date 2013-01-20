/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_LEGENDRE_HPP
#define MATRICES_LEGENDRE_HPP

namespace elem {

template<typename F> 
inline void
Legendre( int n, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Legendre");
#endif
    A.ResizeTo( n, n );
    MakeLegendre( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V> 
inline void
Legendre( int n, DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("Legendre");
#endif
    A.ResizeTo( n, n );
    MakeLegendre( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
MakeLegendre( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("MakeLegendre");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make a non-square matrix Legendre");
    MakeZeros( A );

    const int n = A.Width();
    for( int j=0; j<n-1; ++j )
    {
        const F gamma = F(1) / Pow( F(2)*(j+1), F(2) );
        const F beta = F(1) / (2*Sqrt(F(1)-gamma));
        A.Set( j+1, j, beta );
        A.Set( j, j+1, beta );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V>
inline void
MakeLegendre( DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("MakeLegendre");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make a non-square matrix Legendre");
    MakeZeros( A );

    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();
    const int colStride = A.ColStride();
    const int rowStride = A.RowStride();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const int j = rowShift + jLocal*rowStride;
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = colShift + iLocal*colStride;
            if( j == i+1 || j == i-1 )
            {
                const int k = std::max( i, j );
                const F gamma = F(1) / Pow( F(2)*k, F(2) );
                const F beta = F(1) / (2*Sqrt(F(1)-gamma));
                A.SetLocal( iLocal, jLocal, beta );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef MATRICES_LEGENDRE_HPP
