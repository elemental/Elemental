/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_HILBERT_HPP
#define ELEM_MATRICES_HILBERT_HPP

namespace elem {

template<typename F>
inline void
Hilbert( Matrix<F>& A, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Hilbert");
#endif
    A.ResizeTo( n, n );
    MakeHilbert( A );
}

template<typename F,Distribution U,Distribution V>
inline void
Hilbert( DistMatrix<F,U,V>& A, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Hilbert");
#endif
    A.ResizeTo( n, n );
    MakeHilbert( A );
}

template<typename F> 
inline void
MakeHilbert( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("MakeHilbert");
#endif
    const int m = A.Height();
    const int n = A.Width();
    if( m != n )
        throw std::logic_error("Cannot make a non-square matrix Hilbert");

    const F one = F(1);
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            A.Set( i, j, one/(i+j+1) );
}

template<typename F,Distribution U,Distribution V>
inline void
MakeHilbert( DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry entry("MakeHilbert");
#endif
    const int m = A.Height();
    const int n = A.Width();
    if( m != n )
        throw std::logic_error("Cannot make a non-square matrix Hilbert");

    const F one = F(1);
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
            A.SetLocal( iLoc, jLoc, one/(i+j+1) );
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_HILBERT_HPP
