/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_RIS_HPP
#define ELEM_MATRICES_RIS_HPP

namespace elem {

template<typename F> 
inline void
Ris( Matrix<F>& R, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Ris");
#endif
    const F oneHalf = F(1)/F(2);
    R.ResizeTo( n, n );
    for( int j=0; j<n; ++j )
        for( int i=0; i<n; ++i )
            R.Set( i, j, oneHalf/(n-i-j-oneHalf) );
}

template<typename F,Distribution U,Distribution V>
inline void
Ris( DistMatrix<F,U,V>& R, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Ris");
#endif
    const F oneHalf = F(1)/F(2);
    R.ResizeTo( n, n );
    const int localHeight = R.LocalHeight();
    const int localWidth = R.LocalWidth();
    const int colShift = R.ColShift();
    const int rowShift = R.RowShift();
    const int colStride = R.ColStride();
    const int rowStride = R.RowStride();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const int i = colShift + iLoc*colStride;
            R.SetLocal( iLoc, jLoc, oneHalf/(n-i-j-oneHalf) );
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_RIS_HPP
