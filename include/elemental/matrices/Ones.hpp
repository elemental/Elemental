/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_ONES_HPP
#define ELEM_MATRICES_ONES_HPP

namespace elem {

template<typename T> 
inline void
MakeOnes( Matrix<T>& A )
{
#ifndef RELEASE
    CallStackEntry cse("MakeOnes");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, T(1) );
}

template<typename T,Distribution U,Distribution V>
inline void
MakeOnes( DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry cse("MakeOnes");
#endif
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            A.SetLocal( iLoc, jLoc, T(1) );
}

template<typename T>
inline void
Ones( Matrix<T>& A, Int m, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Ones");
#endif
    A.ResizeTo( m, n );
    MakeOnes( A );
}

template<typename T>
inline Matrix<T>
Ones( Int m, Int n )
{
    Matrix<T> A( m, n );
    MakeOnes( A );
    return A;
}

template<typename T,Distribution U,Distribution V>
inline void
Ones( DistMatrix<T,U,V>& A, Int m, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Ones");
#endif
    A.ResizeTo( m, n );
    MakeOnes( A );
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
Ones( const Grid& g, Int m, Int n )
{
    DistMatrix<T,U,V> A( m, n, g );
    MakeOnes( A );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_ONES_HPP
