/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_REDHEFFER_HPP
#define ELEM_MATRICES_REDHEFFER_HPP

namespace elem {

template<typename T> 
inline void
Redheffer( Matrix<T>& R, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Redheffer");
#endif
    R.ResizeTo( n, n );
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<n; ++i )
        {
            if( j==0 || ((j+1)%(i+1))==0 )
                R.Set( i, j, T(1) );
            else
                R.Set( i, j, T(0) );
        }
    }
}

template<typename T> 
inline Matrix<T>
Redheffer( Int n )
{
    Matrix<T> R;
    Redheffer( R, n );
    return R;
}

template<typename T,Distribution U,Distribution V>
inline void
Redheffer( DistMatrix<T,U,V>& R, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Redheffer");
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
            if( j==0 || ((j+1)%(i+1))==0 )
                R.SetLocal( iLoc, jLoc, T(1) );
            else
                R.SetLocal( iLoc, jLoc, T(0) );
        }
    }
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
Redheffer( const Grid& g, Int n )
{
    DistMatrix<T,U,V> R(g);
    Redheffer( R, n );
    return R;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_REDHEFFER_HPP
