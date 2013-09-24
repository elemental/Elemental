/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_DIAGONAL_HPP
#define ELEM_MATRICES_DIAGONAL_HPP

#include "elemental/matrices/Zeros.hpp"

namespace elem {

template<typename S,typename T> 
inline void
Diagonal( Matrix<S>& D, const std::vector<T>& d )
{
#ifndef RELEASE
    CallStackEntry cse("Diagonal");
#endif
    const Int n = d.size();
    Zeros( D, n, n );

    for( Int j=0; j<n; ++j )
        D.Set( j, j, d[j] );
}

template<typename T> 
inline Matrix<T>
Diagonal( const std::vector<T>& d )
{
    Matrix<T> D;
    Diagonal( D, d );
    return D;
}

template<typename S,typename T,Distribution U,Distribution V>
inline void
Diagonal( DistMatrix<S,U,V>& D, const std::vector<T>& d )
{
#ifndef RELEASE
    CallStackEntry cse("Diagonal");
#endif
    const Int n = d.size();
    Zeros( D, n, n );

    const Int localWidth = D.LocalWidth();
    const Int colShift = D.ColShift();
    const Int rowShift = D.RowShift();
    const Int colStride = D.ColStride();
    const Int rowStride = D.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        if( (j-colShift+colStride) % colStride == 0 )
        {
            const Int iLoc = (j-colShift) / colStride;
            D.SetLocal( iLoc, jLoc, d[j] );
        }
    }
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
Diagonal( const Grid& g, const std::vector<T>& d )
{
    DistMatrix<T,U,V> D(g);
    Diagonal( D, d );
    return D;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_DIAGONAL_HPP
