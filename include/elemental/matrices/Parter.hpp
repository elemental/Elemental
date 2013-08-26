/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_PARTER_HPP
#define ELEM_MATRICES_PARTER_HPP

namespace elem {

template<typename F> 
inline void
Parter( Matrix<F>& P, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Parter");
#endif
    const F oneHalf = F(1)/F(2);
    P.ResizeTo( n, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i )
            P.Set( i, j, F(1)/(i-j+oneHalf) );
}

template<typename F> 
inline Matrix<F>
Parter( Int n )
{
    Matrix<F> P;
    Parter( P, n );
    return P;
}

template<typename F,Distribution U,Distribution V>
inline void
Parter( DistMatrix<F,U,V>& P, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Parter");
#endif
    const F oneHalf = F(1)/F(2);
    P.ResizeTo( n, n );
    const Int localHeight = P.LocalHeight();
    const Int localWidth = P.LocalWidth();
    const Int colShift = P.ColShift();
    const Int rowShift = P.RowShift();
    const Int colStride = P.ColStride();
    const Int rowStride = P.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
            P.SetLocal( iLoc, jLoc, F(1)/(i-j+oneHalf) );
        }
    }
}

template<typename F,Distribution U=MC,Distribution V=MR>
inline DistMatrix<F,U,V>
Parter( const Grid& g, Int n )
{
    DistMatrix<F,U,V> P(g);
    Parter( P, n );
    return P;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_PARTER_HPP
