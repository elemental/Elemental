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
Identity( int m, int n, Matrix<T>& I )
{
#ifndef RELEASE
    PushCallStack("Identity");
#endif
    I.ResizeTo( m, n );
    MakeIdentity( I );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
Identity( int m, int n, DistMatrix<T,U,V>& I )
{
#ifndef RELEASE
    PushCallStack("Identity");
#endif
    I.ResizeTo( m, n );
    MakeIdentity( I );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T> 
inline void
MakeIdentity( Matrix<T>& I )
{
#ifndef RELEASE
    PushCallStack("MakeIdentity");
#endif
    Zero( I );
    const int m = I.Height();
    const int n = I.Width();
    for( int j=0; j<std::min(m,n); ++j )
        I.Set( j, j, T(1) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
MakeIdentity( DistMatrix<T,U,V>& I )
{
#ifndef RELEASE
    PushCallStack("MakeIdentity");
#endif
    Zero( I.LocalMatrix() );

    const int localHeight = I.LocalHeight();
    const int localWidth = I.LocalWidth();
    const int colShift = I.ColShift();
    const int rowShift = I.RowShift();
    const int colStride = I.ColStride();
    const int rowStride = I.RowStride();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const int j = rowShift + jLocal*rowStride;
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = colShift + iLocal*colStride;
            if( i == j )
                I.SetLocal( iLocal, jLocal, T(1) );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
