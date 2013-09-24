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
MinIJ( Matrix<T>& M, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("MinIJ");
#endif
    M.ResizeTo( n, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i )
            M.Set( i, j, std::min(i+1,j+1) );
}

template<typename T> 
inline Matrix<T>
MinIJ( Int n )
{
    Matrix<T> M;
    MinIJ( M, n );
    return M;
}

template<typename T,Distribution U,Distribution V>
inline void
MinIJ( DistMatrix<T,U,V>& M, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("MinIJ");
#endif
    M.ResizeTo( n, n );
    const Int localHeight = M.LocalHeight();
    const Int localWidth = M.LocalWidth();
    const Int colShift = M.ColShift();
    const Int rowShift = M.RowShift();
    const Int colStride = M.ColStride();
    const Int rowStride = M.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
            M.SetLocal( iLoc, jLoc, std::min(i+1,j+1) );
        }
    }
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
MinIJ( const Grid& g, Int n )
{
    DistMatrix<T,U,V> M(g);
    MinIJ( M, n );
    return M;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_MINIJ_HPP
