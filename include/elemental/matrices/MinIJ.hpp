/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_MINIJ_HPP
#define ELEM_MATRICES_MINIJ_HPP

namespace elem {

template<typename T> 
inline void
MinIJ( Matrix<T>& M, int n )
{
#ifndef RELEASE
    CallStackEntry entry("MinIJ");
#endif
    M.ResizeTo( n, n );
    for( int j=0; j<n; ++j )
        for( int i=0; i<n; ++i )
            M.Set( i, j, std::min(i+1,j+1) );
}

template<typename T,Distribution U,Distribution V>
inline void
MinIJ( DistMatrix<T,U,V>& M, int n )
{
#ifndef RELEASE
    CallStackEntry entry("MinIJ");
#endif
    M.ResizeTo( n, n );
    const int localHeight = M.LocalHeight();
    const int localWidth = M.LocalWidth();
    const int colShift = M.ColShift();
    const int rowShift = M.RowShift();
    const int colStride = M.ColStride();
    const int rowStride = M.RowStride();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const int i = colShift + iLoc*colStride;
            M.SetLocal( iLoc, jLoc, std::min(i+1,j+1) );
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_MINIJ_HPP
