/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_JORDAN_HPP
#define MATRICES_JORDAN_HPP

#include "elemental/blas-like/level1/Zero.hpp"

namespace elem {

template<typename T>
inline void
Jordan( Matrix<T>& J, int n, T lambda )
{
#ifndef RELEASE
    CallStackEntry entry("Jordan");
#endif
    J.ResizeTo( n, n );
    MakeJordan( J, lambda );
}

template<typename T,Distribution U,Distribution V>
inline void
Jordan( DistMatrix<T,U,V>& J, int n, T lambda )
{
#ifndef RELEASE
    CallStackEntry entry("Jordan");
#endif
    J.ResizeTo( n, n );
    MakeJordan( J, lambda );
}

template<typename T> 
inline void
MakeJordan( Matrix<T>& J, T lambda )
{
#ifndef RELEASE
    CallStackEntry entry("MakeJordan");
#endif
    Zero( J );
    const int m = J.Height();
    const int n = J.Width();
    for( int j=0; j<std::min(m,n); ++j )
    {
        J.Set( j, j, lambda );
        if( j != 0 )
            J.Set( j-1, j, T(1) );
    }
}

template<typename T,Distribution U,Distribution V>
inline void
MakeJordan( DistMatrix<T,U,V>& J, T lambda )
{
#ifndef RELEASE
    CallStackEntry entry("MakeJordan");
#endif
    Zero( J.Matrix() );

    const int localHeight = J.LocalHeight();
    const int localWidth = J.LocalWidth();
    const int colShift = J.ColShift();
    const int rowShift = J.RowShift();
    const int colStride = J.ColStride();
    const int rowStride = J.RowStride();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const int i = colShift + iLoc*colStride;
            if( i == j )
                J.SetLocal( iLoc, jLoc, lambda );
            else if( i == j-1 )
                J.SetLocal( iLoc, jLoc, T(1) );
        }
    }
}

} // namespace elem

#endif // ifndef MATRICES_JORDAN_HPP
