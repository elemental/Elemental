/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename T> 
inline void
Toeplitz( int m, int n, const std::vector<T>& a, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Toeplitz");
#endif
    const int length = m+n-1;
    if( a.size() != (unsigned)length )
        throw std::logic_error("a was the wrong size");
    A.ResizeTo( m, n );

    for( int j=0; j<n; ++j )
        for( int i=0; i<n; ++i )
            A.Set( i, j, a[i-j+(n-1)] );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
Toeplitz( int m, int n, const std::vector<T>& a, DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("Toeplitz");
#endif
    const int length = m+n-1;
    if( a.size() != (unsigned)length )
        throw std::logic_error("a was the wrong size");
    A.ResizeTo( m, n );

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
            A.SetLocal( iLocal, jLocal, a[i-j+(n-1)] );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
