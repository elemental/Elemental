/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_IDENTITY_HPP
#define ELEM_MATRICES_IDENTITY_HPP

#include "elemental/blas-like/level1/Zero.hpp"

namespace elem {

template<typename T> 
inline void
MakeIdentity( Matrix<T>& I )
{
#ifndef RELEASE
    CallStackEntry cse("MakeIdentity");
#endif
    Zero( I );
    const Int m = I.Height();
    const Int n = I.Width();
    for( Int j=0; j<std::min(m,n); ++j )
        I.Set( j, j, T(1) );
}

template<typename T,Distribution U,Distribution V>
inline void
MakeIdentity( DistMatrix<T,U,V>& I )
{
#ifndef RELEASE
    CallStackEntry cse("MakeIdentity");
#endif
    Zero( I.Matrix() );

    const Int localHeight = I.LocalHeight();
    const Int localWidth = I.LocalWidth();
    const Int colShift = I.ColShift();
    const Int rowShift = I.RowShift();
    const Int colStride = I.ColStride();
    const Int rowStride = I.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
            if( i == j )
                I.SetLocal( iLoc, jLoc, T(1) );
        }
    }
}

template<typename T>
inline void
Identity( Matrix<T>& I, Int m, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Identity");
#endif
    I.ResizeTo( m, n );
    MakeIdentity( I );
}

template<typename T>
inline Matrix<T>
Identity( Int m, Int n )
{
    Matrix<T> I( m, n );
    MakeIdentity( I );
    return I;
}

template<typename T,Distribution U,Distribution V>
inline void
Identity( DistMatrix<T,U,V>& I, Int m, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Identity");
#endif
    I.ResizeTo( m, n );
    MakeIdentity( I );
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
Identity( const Grid& g, Int m, Int n )
{
    DistMatrix<T,U,V> I( m, n, g );
    MakeIdentity( I );
    return I;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_IDENTITY_HPP
