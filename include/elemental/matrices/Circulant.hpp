/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_CIRCULANT_HPP
#define ELEM_MATRICES_CIRCULANT_HPP

namespace elem {

template<typename T> 
inline void
Circulant( Matrix<T>& A, const std::vector<T>& a )
{
#ifndef RELEASE
    CallStackEntry cse("Circulant");
#endif
    const Int n = a.size();
    A.ResizeTo( n, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i )
            A.Set( i, j, a[(i-j+n)%n] );
}

template<typename T> 
inline Matrix<T>
Circulant( const std::vector<T>& a )
{
    Matrix<T> A;
    Circulant( A, a );
    return A;
}

template<typename T,Distribution U,Distribution V>
inline void
Circulant( DistMatrix<T,U,V>& A, const std::vector<T>& a )
{
#ifndef RELEASE
    CallStackEntry cse("Circulant");
#endif
    const Int n = a.size();
    A.ResizeTo( n, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int colShift = A.ColShift();
    const Int rowShift = A.RowShift();
    const Int colStride = A.ColStride();
    const Int rowStride = A.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
            A.SetLocal( iLoc, jLoc, a[(i-j+n)%n] );
        }
    }
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
Circulant( const Grid& g, const std::vector<T>& a )
{
    DistMatrix<T,U,V> A(g);
    Circulant( A, a );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_CIRCULANT_HPP
