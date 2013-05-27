/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_DIAGONAL_HPP
#define MATRICES_DIAGONAL_HPP

#include "elemental/matrices/Zeros.hpp"

namespace elem {

template<typename T> 
inline void
Diagonal( Matrix<T>& D, const std::vector<T>& d )
{
#ifndef RELEASE
    CallStackEntry entry("Diagonal");
#endif
    const int n = d.size();
    Zeros( D, n, n );

    for( int j=0; j<n; ++j )
        D.Set( j, j, d[j] );
}

template<typename T,Distribution U,Distribution V>
inline void
Diagonal( DistMatrix<T,U,V>& D, const std::vector<T>& d )
{
#ifndef RELEASE
    CallStackEntry entry("Diagonal");
#endif
    const int n = d.size();
    Zeros( D, n, n );

    const int localWidth = D.LocalWidth();
    const int colShift = D.ColShift();
    const int rowShift = D.RowShift();
    const int colStride = D.ColStride();
    const int rowStride = D.RowStride();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        if( (j-colShift+colStride) % colStride == 0 )
        {
            const int iLoc = (j-colShift) / colStride;
            D.SetLocal( iLoc, jLoc, d[j] );
        }
    }
}

} // namespace elem

#endif // ifndef MATRICES_DIAGONAL_HPP
