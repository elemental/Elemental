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
Circulant( const std::vector<T>& a, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Circulant");
#endif
    const int n = a.size();
    A.ResizeTo( n, n );
    for( int j=0; j<n; ++j )
        for( int i=0; i<n; ++i )
            A.Set( i, j, a[(i-j+n)%n] );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
Circulant( const std::vector<T>& a, DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("Circulant");
#endif
    const int n = a.size();
    A.ResizeTo( n, n );

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
            A.SetLocal( iLocal, jLocal, a[(i-j+n)%n] );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
