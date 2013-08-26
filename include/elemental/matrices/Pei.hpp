/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_PEI_HPP
#define ELEM_MATRICES_PEI_HPP

namespace elem {

template<typename T> 
inline void
Pei( Matrix<T>& P, Int n, T alpha )
{
#ifndef RELEASE
    CallStackEntry cse("Pei");
#endif
    P.ResizeTo( n, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i )
            P.Set( i, j, T(1) );
    for( Int j=0; j<n; ++j )
        P.Update( j, j, alpha );
}

template<typename T> 
inline Matrix<T>
Pei( Int n, T alpha )
{
    Matrix<T> P;
    Pei( P, n, alpha );
    return P;
}

template<typename T,Distribution U,Distribution V>
inline void
Pei( DistMatrix<T,U,V>& P, Int n, T alpha )
{
#ifndef RELEASE
    CallStackEntry cse("Pei");
#endif
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
            P.SetLocal( iLoc, jLoc, T(1) );
            if( i == j )
                P.UpdateLocal( iLoc, jLoc, alpha );
        }
    }
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
Pei( const Grid& g, Int n, T alpha )
{
    DistMatrix<T,U,V> P(g);
    Pei( P, n, alpha );
    return P;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_PEI_HPP
