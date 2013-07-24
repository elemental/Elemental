/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_TOEPLITZ_HPP
#define ELEM_MATRICES_TOEPLITZ_HPP

namespace elem {

template<typename T> 
inline void
Toeplitz( Matrix<T>& A, int m, int n, const std::vector<T>& a )
{
#ifndef RELEASE
    CallStackEntry entry("Toeplitz");
#endif
    const int length = m+n-1;
    if( a.size() != (unsigned)length )
        throw std::logic_error("a was the wrong size");
    A.ResizeTo( m, n );

    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            A.Set( i, j, a[i-j+(n-1)] );
}

template<typename T,Distribution U,Distribution V>
inline void
Toeplitz( DistMatrix<T,U,V>& A, int m, int n, const std::vector<T>& a )
{
#ifndef RELEASE
    CallStackEntry entry("Toeplitz");
#endif
    const int length = m+n-1;
    if( a.size() != (unsigned)length )
        throw std::logic_error("a was the wrong size");
    A.ResizeTo( m, n );

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
            A.SetLocal( iLoc, jLoc, a[i-j+(n-1)] );
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_TOEPLITZ_HPP
