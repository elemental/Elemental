/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_PARTER_HPP
#define MATRICES_PARTER_HPP

namespace elem {

template<typename F> 
inline void
Parter( Matrix<F>& P, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Parter");
#endif
    const F oneHalf = F(1)/F(2);
    P.ResizeTo( n, n );
    for( int j=0; j<n; ++j )
        for( int i=0; i<n; ++i )
            P.Set( i, j, F(1)/(i-j+oneHalf) );
}

template<typename F,Distribution U,Distribution V>
inline void
Parter( DistMatrix<F,U,V>& P, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Parter");
#endif
    const F oneHalf = F(1)/F(2);
    P.ResizeTo( n, n );
    const int localHeight = P.LocalHeight();
    const int localWidth = P.LocalWidth();
    const int colShift = P.ColShift();
    const int rowShift = P.RowShift();
    const int colStride = P.ColStride();
    const int rowStride = P.RowStride();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const int j = rowShift + jLocal*rowStride;
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = colShift + iLocal*colStride;
            P.SetLocal( iLocal, jLocal, F(1)/(i-j+oneHalf) );
        }
    }
}

} // namespace elem

#endif // ifndef MATRICES_PARTER_HPP
