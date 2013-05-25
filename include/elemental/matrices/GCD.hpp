/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_GCD_HPP
#define MATRICES_GCD_HPP

namespace elem {

template<typename T>
inline void
GCD( Matrix<T>& G, int m, int n )
{
#ifndef RELEASE
    CallStackEntry entry("GCD");
#endif
    G.ResizeTo( m, n );
    MakeGCD( G );
}

template<typename T,Distribution U,Distribution V>
inline void
GCD( DistMatrix<T,U,V>& G, int m, int n )
{
#ifndef RELEASE
    CallStackEntry entry("GCD");
#endif
    G.ResizeTo( m, n );
    MakeGCD( G );
}

template<typename T> 
inline void
MakeGCD( Matrix<T>& G )
{
#ifndef RELEASE
    CallStackEntry entry("MakeGCD");
#endif
    const int m = G.Height();
    const int n = G.Width();
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            G.Set( i, j, GCD(i+1,j+1) );
}

template<typename T,Distribution U,Distribution V>
inline void
MakeGCD( DistMatrix<T,U,V>& G )
{
#ifndef RELEASE
    CallStackEntry entry("MakeGCD");
#endif
    const int localHeight = G.LocalHeight();
    const int localWidth = G.LocalWidth();
    const int colShift = G.ColShift();
    const int rowShift = G.RowShift();
    const int colStride = G.ColStride();
    const int rowStride = G.RowStride();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const int i = colShift + iLoc*colStride;
            G.SetLocal( iLoc, jLoc, GCD(i+1,j+1) );
        }
    }
}

} // namespace elem

#endif // ifndef MATRICES_GCD_HPP
