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
    CallStackEntry cse("MakeGCDMatrix");
#endif
    const Int m = G.Height();
    const Int n = G.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            G.Set( i, j, T(GCD(i+1,j+1)) );
}

template<typename T,Distribution U,Distribution V>
inline void
MakeGCDMatrix( DistMatrix<T,U,V>& G )
{
#ifndef RELEASE
    CallStackEntry cse("MakeGCDMatrix");
#endif
    const Int localHeight = G.LocalHeight();
    const Int localWidth = G.LocalWidth();
    const Int colShift = G.ColShift();
    const Int rowShift = G.RowShift();
    const Int colStride = G.ColStride();
    const Int rowStride = G.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
            G.SetLocal( iLoc, jLoc, T(GCD(i+1,j+1)) );
        }
    }
}

template<typename T>
inline void
GCDMatrix( Matrix<T>& G, Int m, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("GCDMatrix");
#endif
    G.ResizeTo( m, n );
    MakeGCDMatrix( G );
}

template<typename T>
inline Matrix<T>
GCDMatrix( Int m, Int n )
{
    Matrix<T> G( m, n );
    MakeGCDMatrix( G );
    return G;
}

template<typename T,Distribution U,Distribution V>
inline void
GCDMatrix( DistMatrix<T,U,V>& G, Int m, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("GCDMatrix");
#endif
    G.ResizeTo( m, n );
    MakeGCDMatrix( G );
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
GCDMatrix( const Grid& g, Int m, Int n )
{
    DistMatrix<T,U,V> G( m, n, g );
    MakeGCDMatrix( G );
    return G;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_GCDMATRIX_HPP
