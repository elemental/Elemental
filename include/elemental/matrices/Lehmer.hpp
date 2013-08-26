/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_LEHMER_HPP
#define ELEM_MATRICES_LEHMER_HPP

namespace elem {

template<typename F> 
inline void
Lehmer( Matrix<F>& L, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Lehmer");
#endif
    L.ResizeTo( n, n );
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<j; ++i )
            L.Set( i, j, F(i+1)/F(j+1) );
        for( Int i=j; i<n; ++i )
            L.Set( i, j, F(j+1)/F(i+1) );
    }
}

template<typename F> 
inline Matrix<F>
Lehmer( Int n )
{
    Matrix<F> L;
    Lehmer( L, n );
    return L;
}

template<typename F,Distribution U,Distribution V>
inline void
Lehmer( DistMatrix<F,U,V>& L, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Lehmer");
#endif
    L.ResizeTo( n, n );
    const Int localHeight = L.LocalHeight();
    const Int localWidth = L.LocalWidth();
    const Int colShift = L.ColShift();
    const Int rowShift = L.RowShift();
    const Int colStride = L.ColStride();
    const Int rowStride = L.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
            if( i < j )
                L.SetLocal( iLoc, jLoc, F(i+1)/F(j+1) );
            else
                L.SetLocal( iLoc, jLoc, F(j+1)/F(i+1) );
        }
    }
}

template<typename F,Distribution U=MC,Distribution V=MR>
inline DistMatrix<F,U,V>
Lehmer( const Grid& g, Int n )
{
    DistMatrix<F,U,V> L(g);
    Lehmer( L, n );
    return L;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_LEHMER_HPP
