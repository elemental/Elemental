/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_KMS_HPP
#define MATRICES_KMS_HPP

namespace elem {

template<typename T> 
inline void
KMS( Matrix<T>& K, int n, T rho )
{
#ifndef RELEASE
    CallStackEntry entry("KMS");
#endif
    for( int j=0; j<n; ++j )
    {
        for( int i=0; i<j; ++i )
            K.Set( i, j, Pow(rho,j-i) );
        for( int i=j; i<n; ++i )
            K.Set( i, j, Conj(Pow(rho,i-j)) );
    }
}

template<typename T,Distribution U,Distribution V>
inline void
KMS( DistMatrix<T,U,V>& K, int n, T rho )
{
#ifndef RELEASE
    CallStackEntry entry("KMS");
#endif
    const int localHeight = K.LocalHeight();
    const int localWidth = K.LocalWidth();
    const int colShift = K.ColShift();
    const int rowShift = K.RowShift();
    const int colStride = K.ColStride();
    const int rowStride = K.RowStride();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const int j = rowShift + jLocal*rowStride;
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = colShift + iLocal*colStride;
            if( i < j )
                K.SetLocal( iLocal, jLocal, Pow(rho,j-i) );
            else
                K.SetLocal( iLocal, jLocal, Conj(Pow(rho,i-j)) );
        }
    }
}

} // namespace elem

#endif // ifndef MATRICES_KMS_HPP
