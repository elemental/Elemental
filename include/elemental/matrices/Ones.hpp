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
Ones( int m, int n, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Ones");
#endif
    A.ResizeTo( m, n );
    MakeOnes( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
Ones( int m, int n, DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("Ones");
#endif
    A.ResizeTo( m, n );
    MakeOnes( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T> 
inline void
MakeOnes( Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("MakeOnes");
#endif
    const int m = A.Height();
    const int n = A.Width();
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            A.Set( i, j, T(1) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
MakeOnes( DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("MakeOnes");
#endif
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            A.SetLocal( iLocal, jLocal, T(1) );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef MATRICES_ONES_HPP
