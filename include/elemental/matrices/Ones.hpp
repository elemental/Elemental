/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_ONES_HPP
#define MATRICES_ONES_HPP

namespace elem {

template<typename T>
inline void
Ones( Matrix<T>& A, int m, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Ones");
#endif
    A.ResizeTo( m, n );
    MakeOnes( A );
}

template<typename T,Distribution U,Distribution V>
inline void
Ones( DistMatrix<T,U,V>& A, int m, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Ones");
#endif
    A.ResizeTo( m, n );
    MakeOnes( A );
}

template<typename T> 
inline void
MakeOnes( Matrix<T>& A )
{
#ifndef RELEASE
    CallStackEntry entry("MakeOnes");
#endif
    const int m = A.Height();
    const int n = A.Width();
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            A.Set( i, j, T(1) );
}

template<typename T,Distribution U,Distribution V>
inline void
MakeOnes( DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry entry("MakeOnes");
#endif
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            A.SetLocal( iLocal, jLocal, T(1) );
}

} // namespace elem

#endif // ifndef MATRICES_ONES_HPP
