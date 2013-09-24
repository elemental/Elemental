/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_RIEMANN_HPP
#define ELEM_MATRICES_RIEMANN_HPP

namespace elem {

template<typename T> 
inline void
Riemann( Matrix<T>& R, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Riemann");
#endif
    R.ResizeTo( n, n );
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<n; ++i )
        {
            if( ((j+2)%(i+2))==0 )
                R.Set( i, j, T(i+1) );
            else
                R.Set( i, j, T(-1) );
        }
    }
}

template<typename T> 
inline Matrix<T>
Riemann( Int n )
{
    Matrix<T> R;
    Riemann( R, n );
    return R;
}

template<typename T,Distribution U,Distribution V>
inline void
Riemann( DistMatrix<T,U,V>& R, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Riemann");
#endif
    R.ResizeTo( n, n );
    const Int localHeight = R.LocalHeight();
    const Int localWidth = R.LocalWidth();
    const Int colShift = R.ColShift();
    const Int rowShift = R.RowShift();
    const Int colStride = R.ColStride();
    const Int rowStride = R.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
            if( ((j+2)%(i+2))==0 )
                R.SetLocal( iLoc, jLoc, T(i+1) );
            else
                R.SetLocal( iLoc, jLoc, T(-1) );
        }
    }
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
Riemann( const Grid& g, Int n )
{
    DistMatrix<T,U,V> R(g);
    Riemann( R, n );
    return R;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_RIEMANN_HPP
