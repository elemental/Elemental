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
Wilkinson( int k, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Wilkinson");
#endif
    const int n = 2*k+1;
    A.ResizeTo( n, n );
    MakeZeros( A );

    for( int j=0; j<n; ++j )
    {
        if( j <= k )
            A.Set( j, j, T(k-j) );
        else
            A.Set( j, j, T(j-k) );

        if( j > 0 )
            A.Set( j-1, j, T(1) );
        if( j < n-1 )
            A.Set( j+1, j, T(1) );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
Wilkinson( int k, DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("Wilkinson");
#endif
    const int n = 2*k+1;
    A.ResizeTo( n, n );
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
            {
                if( j <= k )
                    A.SetLocal( iLocal, jLocal, T(k-j) );
                else
                    A.SetLocal( iLocal, jLocal, T(j-k) );
            }
            else if( i == j-1 || i == j+1 )
                A.SetLocal( iLocal, jLocal, T(1) );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
