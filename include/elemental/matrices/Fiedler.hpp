/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_FIEDLER_HPP
#define MATRICES_FIEDLER_HPP

namespace elem {

template<typename F> 
inline void
Fiedler( Matrix<F>& A, const std::vector<F>& c )
{
#ifndef RELEASE
    CallStackEntry entry("Fiedler");
#endif
    const int n = c.size();
    A.ResizeTo( n, n );
    for( int j=0; j<n; ++j )
        for( int i=0; i<n; ++i )
            A.Set( i, j, Abs(c[i]-c[j]) );
}

template<typename F,Distribution U,Distribution V>
inline void
Fiedler( DistMatrix<F,U,V>& A, const std::vector<F>& c )
{
#ifndef RELEASE
    CallStackEntry entry("Fiedler");
#endif
    const int n = c.size();
    A.ResizeTo( n, n );
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
            A.SetLocal( iLoc, jLoc, Abs(c[i]-c[j]) );
        }
    }
}

} // namespace elem

#endif // ifndef MATRICES_FIEDLER_HPP
