/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_LEHMER_HPP
#define MATRICES_LEHMER_HPP

namespace elem {

template<typename F> 
inline void
Lehmer( Matrix<F>& L, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Lehmer");
#endif
    L.ResizeTo( n, n );
    for( int j=0; j<n; ++j )
    {
        for( int i=0; i<j; ++i )
            L.Set( i, j, F(i+1)/F(j+1) );
        for( int i=j; i<n; ++i )
            L.Set( i, j, F(j+1)/F(i+1) );
    }
}

template<typename F,Distribution U,Distribution V>
inline void
Lehmer( DistMatrix<F,U,V>& L, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Lehmer");
#endif
    L.ResizeTo( n, n );
    const int localHeight = L.LocalHeight();
    const int localWidth = L.LocalWidth();
    const int colShift = L.ColShift();
    const int rowShift = L.RowShift();
    const int colStride = L.ColStride();
    const int rowStride = L.RowStride();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const int i = colShift + iLoc*colStride;
            if( i < j )
                L.SetLocal( iLoc, jLoc, F(i+1)/F(j+1) );
            else
                L.SetLocal( iLoc, jLoc, F(j+1)/F(i+1) );
        }
    }
}

} // namespace elem

#endif // ifndef MATRICES_LEHMER_HPP
