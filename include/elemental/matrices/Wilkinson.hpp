/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_WILKINSON_HPP
#define MATRICES_WILKINSON_HPP

#include "elemental/matrices/Zeros.hpp"

namespace elem {

template<typename T> 
inline void
Wilkinson( Matrix<T>& A, int k )
{
#ifndef RELEASE
    CallStackEntry entry("Wilkinson");
#endif
    const int n = 2*k+1;
    A.ResizeTo( n, n );
    MakeZeros( A );

    for( int j=0; j<n; ++j )
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

template<typename T,Distribution U,Distribution V>
inline void
Wilkinson( DistMatrix<T,U,V>& A, int k )
{
#ifndef RELEASE
    CallStackEntry entry("Wilkinson");
#endif
    const int n = 2*k+1;
    A.ResizeTo( n, n );
    MakeZeros( A );

    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();
    const int colStride = A.ColStride();
    const int rowStride = A.RowStride();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const int i = colShift + iLoc*colStride;
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

} // namespace elem

#endif // ifndef MATRICES_WILKINSON_HPP
