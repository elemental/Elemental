/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename T> 
inline void
OneTwoOne( int n, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("OneTwoOne");
#endif
    A.ResizeTo( n, n );
    MakeOneTwoOne( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V> 
inline void
OneTwoOne( int n, DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("OneTwoOne");
#endif
    A.ResizeTo( n, n );
    MakeOneTwoOne( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T> 
inline void
MakeOneTwoOne( Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("MakeOneTwoOne");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make a non-square matrix 1-2-1");
    MakeZeros( A );

    const int n = A.Width();
    for( int j=0; j<n; ++j )
    {
        A.Set( j, j, T(2) );
        if( j < n-1 )
        {
            A.Set( j+1, j, T(1) );
            A.Set( j, j+1, T(1) );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
MakeOneTwoOne( DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("MakeOneTwoOne");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make a non-square matrix 1-2-1");
    MakeZeros( A );

    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();
    const int colStride = A.ColStride();
    const int rowStride = A.RowStride();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const int j = rowShift + jLocal*rowStride;
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = colShift + iLocal*colStride;
            if( i == j )
                A.SetLocal( iLocal, jLocal, T(2) );
            else if( i == j-1 || i == j+1 )
                A.SetLocal( iLocal, jLocal, T(1) );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
