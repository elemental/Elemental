/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_IDENTITY_HPP
#define MATRICES_IDENTITY_HPP

#include "elemental/blas-like/level1/Zero.hpp"

namespace elem {

template<typename T>
inline void
Identity( Matrix<T>& I, int m, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Identity");
#endif
    I.ResizeTo( m, n );
    MakeIdentity( I );
}

template<typename T,Distribution U,Distribution V>
inline void
Identity( DistMatrix<T,U,V>& I, int m, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Identity");
#endif
    I.ResizeTo( m, n );
    MakeIdentity( I );
}

template<typename T> 
inline void
MakeIdentity( Matrix<T>& I )
{
#ifndef RELEASE
    CallStackEntry entry("MakeIdentity");
#endif
    Zero( I );
    const int m = I.Height();
    const int n = I.Width();
    for( int j=0; j<std::min(m,n); ++j )
        I.Set( j, j, T(1) );
}

template<typename T,Distribution U,Distribution V>
inline void
MakeIdentity( DistMatrix<T,U,V>& I )
{
#ifndef RELEASE
    CallStackEntry entry("MakeIdentity");
#endif
    Zero( I.Matrix() );

    const int localHeight = I.LocalHeight();
    const int localWidth = I.LocalWidth();
    const int colShift = I.ColShift();
    const int rowShift = I.RowShift();
    const int colStride = I.ColStride();
    const int rowStride = I.RowStride();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const int i = colShift + iLoc*colStride;
            if( i == j )
                I.SetLocal( iLoc, jLoc, T(1) );
        }
    }
}

} // namespace elem

#endif // ifndef MATRICES_IDENTITY_HPP
