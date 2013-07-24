/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_GCDMATRIX_HPP
#define ELEM_MATRICES_GCDMATRIX_HPP

namespace elem {

template<typename T>
inline void
MakeGCDMatrix( Matrix<T>& G )
{
#ifndef RELEASE
    CallStackEntry entry("MakeGCDMatrix");
#endif
    const int m = G.Height();
    const int n = G.Width();
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            G.Set( i, j, T(GCD(i+1,j+1)) );
}

template<typename T,Distribution U,Distribution V>
inline void
MakeGCDMatrix( DistMatrix<T,U,V>& G )
{
#ifndef RELEASE
    CallStackEntry entry("MakeGCDMatrix");
#endif
    const int localHeight = G.LocalHeight();
    const int localWidth = G.LocalWidth();
    const int colShift = G.ColShift();
    const int rowShift = G.RowShift();
    const int colStride = G.ColStride();
    const int rowStride = G.RowStride();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const int i = colShift + iLoc*colStride;
            G.SetLocal( iLoc, jLoc, T(GCD(i+1,j+1)) );
        }
    }
}

template<typename T>
inline void
GCDMatrix( Matrix<T>& G, int m, int n )
{
#ifndef RELEASE
    CallStackEntry entry("GCDMatrix");
#endif
    G.ResizeTo( m, n );
    MakeGCDMatrix( G );
}

template<typename T,Distribution U,Distribution V>
inline void
GCDMatrix( DistMatrix<T,U,V>& G, int m, int n )
{
#ifndef RELEASE
    CallStackEntry entry("GCDMatrix");
#endif
    G.ResizeTo( m, n );
    MakeGCDMatrix( G );
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_GCDMATRIX_HPP
