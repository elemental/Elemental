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
Diagonal( const std::vector<T>& d, Matrix<T>& D )
{
#ifndef RELEASE
    PushCallStack("Diagonal");
#endif
    const int n = d.size();
    D.ResizeTo( n, n );
    MakeZeros( D );

    for( int j=0; j<n; ++j )
        D.Set( j, j, d[j] );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
Diagonal( const std::vector<T>& d, DistMatrix<T,U,V>& D )
{
#ifndef RELEASE
    PushCallStack("Diagonal");
#endif
    const int n = d.size();
    D.ResizeTo( n, n );
    MakeZeros( D );

    const int localWidth = D.LocalWidth();
    const int colShift = D.ColShift();
    const int rowShift = D.RowShift();
    const int colStride = D.ColStride();
    const int rowStride = D.RowStride();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const int j = rowShift + jLocal*rowStride;
        if( (j-colShift+colStride) % colStride == 0 )
        {
            const int iLocal = (j-colShift) / colStride;
            D.SetLocal( iLocal, jLocal, d[j] );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
