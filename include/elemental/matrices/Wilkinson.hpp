/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_WILKINSON_HPP
#define ELEM_MATRICES_WILKINSON_HPP

#include "elemental/matrices/Zeros.hpp"

namespace elem {

template<typename T> 
inline void
Wilkinson( Matrix<T>& A, Int k )
{
#ifndef RELEASE
    CallStackEntry cse("Wilkinson");
#endif
    const Int n = 2*k+1;
    A.ResizeTo( n, n );
    MakeZeros( A );

    for( Int j=0; j<n; ++j )
    {
        if( j <= k )
            A.Set( j, j, T(k-j) );
        else
            A.Set( j, j, T(j-k) );

        if( j > 0 )
            A.Set( j-1, j, T(1) );
        if( j < n-1 )
            A.Set( j+1, j, T(1) );
    }
}

template<typename T> 
inline Matrix<T>
Wilkinson( Int k )
{
    Matrix<T> A;
    Wilkinson( A, k );
    return A;
}

template<typename T,Distribution U,Distribution V>
inline void
Wilkinson( DistMatrix<T,U,V>& A, Int k )
{
#ifndef RELEASE
    CallStackEntry cse("Wilkinson");
#endif
    const Int n = 2*k+1;
    A.ResizeTo( n, n );
    MakeZeros( A );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int colShift = A.ColShift();
    const Int rowShift = A.RowShift();
    const Int colStride = A.ColStride();
    const Int rowStride = A.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
            if( i == j )
            {
                if( j <= k )
                    A.SetLocal( iLoc, jLoc, T(k-j) );
                else
                    A.SetLocal( iLoc, jLoc, T(j-k) );
            }
            else if( i == j-1 || i == j+1 )
                A.SetLocal( iLoc, jLoc, T(1) );
        }
    }
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
Wilkinson( const Grid& g, Int k )
{
    DistMatrix<T,U,V> A( g );
    Wilkinson( A, k );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_WILKINSON_HPP
