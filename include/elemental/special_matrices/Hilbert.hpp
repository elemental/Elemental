/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename F>
inline void
Hilbert( int n, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Hilbert");
#endif
    A.ResizeTo( n, n );
    MakeHilbert( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V>
inline void
Hilbert( int n, DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("Hilbert");
#endif
    A.ResizeTo( n, n );
    MakeHilbert( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
MakeHilbert( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("MakeHilbert");
#endif
    const int m = A.Height();
    const int n = A.Width();
    if( m != n )
        throw std::logic_error("Cannot make a non-square matrix Hilbert");

    const F one = F(1);
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            A.Set( i, j, one/(i+j+1) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V>
inline void
MakeHilbert( DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("MakeHilbert");
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
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const int j = rowShift + jLocal*rowStride;
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = colShift + iLocal*colStride;
            A.SetLocal( iLocal, jLocal, one/(i+j+1) );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
