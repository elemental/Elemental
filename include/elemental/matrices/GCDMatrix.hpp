/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_GCDMATRIX_HPP
#define ELEM_GCDMATRIX_HPP

namespace elem {

template<typename T>
inline void
MakeGCDMatrix( Matrix<T>& G )
{
    DEBUG_ONLY(CallStackEntry cse("MakeGCDMatrix"))
    const Int m = G.Height();
    const Int n = G.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            G.Set( i, j, T(GCD(i+1,j+1)) );
}

template<typename T,Dist U,Dist V>
inline void
MakeGCDMatrix( DistMatrix<T,U,V>& G )
{
    DEBUG_ONLY(CallStackEntry cse("MakeGCDMatrix"))
    const Int localHeight = G.LocalHeight();
    const Int localWidth = G.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = G.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = G.GlobalRow(iLoc);
            G.SetLocal( iLoc, jLoc, T(GCD(i+1,j+1)) );
        }
    }
}

template<typename T,Dist U,Dist V>
inline void
MakeGCDMatrix( BlockDistMatrix<T,U,V>& G )
{
    DEBUG_ONLY(CallStackEntry cse("MakeGCDMatrix"))
    const Int localHeight = G.LocalHeight();
    const Int localWidth = G.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = G.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = G.GlobalRow(iLoc);
            G.SetLocal( iLoc, jLoc, T(GCD(i+1,j+1)) );
        }
    }
}

template<typename T>
inline void
GCDMatrix( Matrix<T>& G, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("GCDMatrix"))
    G.Resize( m, n );
    MakeGCDMatrix( G );
}

template<typename T,Dist U,Dist V>
inline void
GCDMatrix( DistMatrix<T,U,V>& G, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("GCDMatrix"))
    G.Resize( m, n );
    MakeGCDMatrix( G );
}

template<typename T,Dist U,Dist V>
inline void
GCDMatrix( BlockDistMatrix<T,U,V>& G, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("GCDMatrix"))
    G.Resize( m, n );
    MakeGCDMatrix( G );
}

} // namespace elem

#endif // ifndef ELEM_GCDMATRIX_HPP
